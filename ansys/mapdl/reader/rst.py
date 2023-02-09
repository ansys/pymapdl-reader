"""Read ANSYS binary result files (*.rst)

Used:
.../ansys/customize/include/fdresu.inc
"""
from collections.abc import Iterable, Sequence
from functools import wraps
import os
import pathlib
from threading import Thread
import time
from typing import Union
import warnings

import numpy as np
import pyvista as pv
from pyvista import _vtk as vtk
from pyvista.themes import DefaultTheme
from tqdm import tqdm

from ansys.mapdl.reader import _binary_reader, _reader, elements
from ansys.mapdl.reader._binary_reader import (
    AnsysFile,
    cells_with_all_nodes,
    cells_with_any_nodes,
    populate_surface_element_result,
)
from ansys.mapdl.reader._mp_keys import mp_keys
from ansys.mapdl.reader._rst_keys import (
    DOF_REF,
    STR_LIM_REF,
    boundary_condition_index_table,
    element_index_table_info,
    geometry_header_keys,
    result_header_keys,
    solution_data_header_keys,
    solution_header_keys_dp,
)
from ansys.mapdl.reader.common import (
    PRINCIPAL_STRESS_TYPES,
    STRAIN_TYPES,
    STRESS_TYPES,
    THERMAL_STRAIN_TYPES,
    AnsysBinary,
    parse_header,
    read_standard_header,
    read_table,
    rotate_to_global,
)
from ansys.mapdl.reader.mesh import Mesh
from ansys.mapdl.reader.misc import break_apart_surface, vtk_cell_info
from ansys.mapdl.reader.rst_avail import AvailableResults


def access_bit(data, num):
    """Access a single bit from a bitmask"""
    base = int(num // 8)
    shift = int(num % 8)
    return (data & (1 << shift)) >> shift


EMAIL_ME = """Please raise an issue at:
https://github.com/pyansys/pymapdl-reader/issues
Or email the PyAnsys support team at pyansys.support@ansys.com
"""
np.seterr(divide="ignore", invalid="ignore")


# Pointer information from ansys interface manual
# =============================================================================
# Individual element index table
ELEMENT_INDEX_TABLE_KEYS = [
    "EMS",
    "ENF",
    "ENS",
    "ENG",
    "EGR",
    "EEL",
    "EPL",
    "ECR",
    "ETH",
    "EUL",
    "EFX",
    "ELF",
    "EMN",
    "ECD",
    "ENL",
    "EHC",
    "EPT",
    "ESF",
    "EDI",
    "ETB",
    "ECT",
    "EXY",
    "EBA",
    "ESV",
    "MNL",
]

ELEMENT_RESULT_NCOMP = {
    "ENS": 6,
    "EEL": 7,
    "EPL": 7,
    "ECR": 7,
    "ETH": 8,
    "ENL": 10,
    "EDI": 7,
}


class Result(AnsysBinary):
    """Reads a binary ANSYS result file.

    Parameters
    ----------
    filename : str, pathlib.Path, optional
        Filename of the ANSYS binary result file.
    ignore_cyclic : bool, optional
        Ignores any cyclic properties.
    read_mesh : bool, optional
        Debug parameter.  Set to False to disable reading in the
        mesh from the result file.
    parse_vtk : bool, optional
       Set to ``False`` to skip the parsing the mesh as a VTK
       UnstructuredGrid, which might take a long time for large models.

    Examples
    --------
    >>> from ansys.mapdl import reader as pymapdl_reader
    >>> rst = pymapdl_reader.read_binary('file.rst')
    """

    ELEMENT_INDEX_TABLE_KEYS = ELEMENT_INDEX_TABLE_KEYS
    ELEMENT_RESULT_NCOMP = ELEMENT_RESULT_NCOMP

    def __init__(self, filename, read_mesh=True, parse_vtk=True, **kwargs):
        """Load basic result information from result file and init the rst object."""
        self._filename = pathlib.Path(filename)
        self._cfile = AnsysFile(str(filename))
        self._resultheader = self._read_result_header()
        self._animating = False
        self.__element_map = None
        self._ignore154 = True

        # Get the total number of results and log it
        self.nsets = len(self._resultheader["rpointers"])
        # LOG.debug('There are %d result(s) in this file', self.nsets)

        # Get indices to resort nodal and element results
        self._sidx = np.argsort(self._resultheader["neqv"])
        self._sidx_elem = np.argsort(self._resultheader["eeqv"])

        # Store node numbering in ANSYS
        self._neqv = self._resultheader["neqv"]
        self._eeqv = self._resultheader["eeqv"]  # unsorted

        # cache geometry header
        table = self.read_record(self._resultheader["ptrGEO"])
        self._geometry_header = parse_header(table, geometry_header_keys)
        self._element_table = self._load_element_table()

        self._materials = None
        self._section_data = None
        self._available_results = AvailableResults(
            self._resultheader["AvailData"], self._is_thermal
        )

        # store mesh for later retrieval
        self._mesh = None
        if read_mesh:
            self._store_mesh(parse_vtk)

    @property
    def filename(self) -> str:
        """String form of the filename. This property is read-only."""
        return str(self._filename)

    @property
    def pathlib_filename(self) -> pathlib.Path:
        """Return the ``pathlib.Path`` version of the filename. This property can not be set."""
        return self._filename

    @property
    def mesh(self):
        """Mesh from result file.

        Examples
        --------
        >>> rst.mesh
        ANSYS Mesh
          Number of Nodes:              1448
          Number of Elements:           226
          Number of Element Types:      1
          Number of Node Components:    0
          Number of Element Components: 0
        """
        if self._mesh is None:
            raise ValueError(
                "Pass ``read_mesh=True`` to store the mesh"
                " when initializing the result"
            )
        return self._mesh

    @property
    def n_sector(self):
        """Number of sectors"""
        return self._resultheader["nSector"]

    @property
    def n_results(self):
        """Number of results"""
        return len(self.time_values)

    def _read_result_header(self):
        """Returns pointers used to access results from an ANSYS result file.

        Parameters
        ----------
        filename : string
            Filename of result file.

        Returns
        -------
        resultheader : dictionary
            Result header
        """
        # consider moving this to the main class
        standard_header = read_standard_header(self.filename)

        # Read .RST FILE HEADER
        header = parse_header(self.read_record(103), result_header_keys)
        resultheader = merge_two_dicts(header, standard_header)

        # Read nodal equivalence table
        resultheader["neqv"] = self.read_record(resultheader["ptrNOD"])

        # Read nodal equivalence table
        resultheader["eeqv"] = self.read_record(resultheader["ptrELM"])

        # Read table of pointers to locations of results
        nsets = resultheader["nsets"]

        # Data sets index table. This record contains the record pointers
        # for the beginning of each data set. The first resmax records are
        # the first 32 bits of the index, the second resmax records are
        # the second 32 bits f.seek((ptrDSIl + 0) * 4)
        record = self.read_record(resultheader["ptrDSI"])
        raw0 = record[: resultheader["resmax"]].tobytes()
        raw1 = record[resultheader["resmax"] :].tobytes()

        # this combines the two ints, not that efficient
        subraw0 = [raw0[i * 4 : (i + 1) * 4] for i in range(nsets)]
        subraw1 = [raw1[i * 4 : (i + 1) * 4] for i in range(nsets)]
        longraw = [subraw0[i] + subraw1[i] for i in range(nsets)]
        longraw = b"".join(longraw)
        rpointers = np.frombuffer(longraw, "i8")

        assert (rpointers >= 0).all(), "Data set index table has negative pointers"
        resultheader["rpointers"] = rpointers

        # read in time values
        record = self.read_record(resultheader["ptrTIM"])
        resultheader["time_values"] = record[: resultheader["nsets"]]

        # load harmonic index of each result
        if resultheader["ptrCYC"]:
            record = self.read_record(resultheader["ptrCYC"])
            hindex = record[: resultheader["nsets"]]

            # ansys 15 doesn't track negative harmonic indices
            if not np.any(hindex < -1):
                # check if duplicate frequencies
                tvalues = resultheader["time_values"]
                for i in range(tvalues.size - 1):
                    # adjust tolarance(?)
                    if np.isclose(tvalues[i], tvalues[i + 1]):
                        hindex[i + 1] *= -1

            resultheader["hindex"] = hindex

        # load step table with columns:
        # [loadstep, substep, and cumulative]
        record = self.read_record(resultheader["ptrLSP"])
        resultheader["ls_table"] = record[: resultheader["nsets"] * 3].reshape(-1, 3)

        return resultheader

    def parse_coordinate_system(self):
        """Reads in coordinate system information from a binary result
        file.

        Returns
        -------
        c_systems : dict
            Dictionary containing one entry for each defined
            coordinate system.  If no non-standard coordinate systems
            have been defined, there will be only one None.  First
            coordinate system is assumed to be global cartesian.

        Notes
        -----
        euler angles : [THXY, THYZ, THZX]

        - First rotation about local Z (positive X toward Y).
        - Second rotation about local X (positive Y toward Z).
        - Third rotation about local Y (positive Z toward X).

        PAR1
        Used for elliptical, spheroidal, or toroidal systems. If KCS =
        1 or 2, PAR1 is the ratio of the ellipse Y-axis radius to
        X-axis radius (defaults to 1.0 (circle)). If KCS = 3, PAR1 is
        the major radius of the torus.

        PAR2
        Used for spheroidal systems. If KCS = 2, PAR2 = ratio of
        ellipse Z-axis radius to X-axis radius (defaults to 1.0
        (circle)).

        Coordinate system type:
            - 0: Cartesian
            - 1: Cylindrical (circular or elliptical)
            - 2: Spherical (or spheroidal)
            - 3: Toroidal
        """
        c_systems = [None]

        # load coordinate system index table
        ptr_csy = self._geometry_header["ptrCSY"]
        if ptr_csy:
            csy, sz = self.read_record(ptr_csy, return_bufsize=True)
        else:
            return c_systems

        # Assemble element record pointers relative to ptrETY
        if self._map_flag:
            csys_record_pointers = self.read_record(ptr_csy + sz)
        else:
            maxcsy = self._geometry_header["maxcsy"]
            csys_record_pointers = csy[:maxcsy]

        # parse each coordinate system
        # The items stored in each record:
        # * Items 1-9  are the transformation matrix.
        # * Items 10-12 are the coordinate system origin (XC,YC,ZC).
        # * Items 13-14 are the coordinate system parameters (PAR1, PAR2).
        # * Items 16-18 are the angles used to define the coordinate system.
        # * Items 19-20 are theta and phi singularity keys.
        # * Item 21 is the coordinate system type (0, 1, 2, or 3).
        # * Item 22 is the coordinate system reference number.
        for csys_record_pointer in csys_record_pointers:
            if not csys_record_pointer:
                c_system = None
            else:
                data = self.read_record(ptr_csy + csys_record_pointer)
                c_system = {
                    "transformation matrix": np.array(data[:9].reshape(-1, 3)),
                    "origin": np.array(data[9:12]),
                    "PAR1": data[12],
                    "PAR2": data[13],
                    "euler angles": data[15:18],  # may not be euler
                    "theta singularity": data[18],
                    "phi singularity": data[19],
                    "type": int(data[20]),
                    "reference num": int(
                        data[21],
                    ),
                }
            c_systems.append(c_system)

        return c_systems

    def _load_materials(self):
        """Reads the material properties from the result file.

        Materials property table is stored at ptrMAT, which contains
        pointers to each material property array for each material.
        Each array contains a material property (e.g. EX, or elastic
        modulus in the X direction) and there may be one or many
        values for each material (up to 101).
        """

        def read_mat_data(ptr):
            """reads double precision material data given an offset from ptrMAT"""
            arr, sz = self.read_record(self._geometry_header["ptrMAT"] + ptr, True)
            is_single_value = arr[arr != 0].size == 1

            if is_single_value:
                return arr[arr != 0][0]
            else:
                # We need to account for 0 at the temp and properties.
                temps_ = arr[:sz]
                values_ = arr[-1 : -1 * (sz + 1) : -1]

                prop_ = np.vstack((temps_, values_))
                for each in range(prop_.shape[1] - 1, -1, -1):
                    # checking that two records are empty starting from the end.
                    if (
                        temps_[each] == 0
                        and values_[each] == 0
                        and temps_[each - 1] == 0
                        and values_[each - 1] == 0
                    ):
                        continue  # they are zero, so we keep going.
                    else:
                        # found a non-zero report
                        break

                prop_ = prop_[:, :each]
                prop_[1, :] = prop_[1, ::-1]
                return prop_

        mat_table = self.read_record(self._geometry_header["ptrMAT"])
        if mat_table[0] != -101:  # pragma: no cover
            raise RuntimeError("Legacy record: Unable to read this material record.")
        self._materials = {}

        # size of record (from fdresu.inc)
        n_mat_prop = self._geometry_header.get("nMatProp")
        for i in range(self._geometry_header["nummat"]):
            # Material Record pointers for the Material Id at I th Position
            # is stored from location.
            # NOTE: somewhere between 182 and 194 this changed from:
            # 3+(I-1)*(158+1)+2 to 3+(I-1)*(158+1)+1+158
            # to
            # 3+(I-1)*(nMatProp+1)+2 to 3+(I-1)*(nMatProp+1)+1+nMatProp

            # pointers to the material data for each material
            legacy = not n_mat_prop  # < v194 (that we know)
            if legacy:
                mat_data_ptr = mat_table[3 + 158 * i : 3 + 159 * i + 159]
            else:
                idx_st = 3 + i * (n_mat_prop + 1)
                idx_en = 3 + i * (n_mat_prop + 1) + 1 + n_mat_prop
                mat_data_ptr = mat_table[idx_st:idx_en]

            self.read_record(self._geometry_header["ptrMAT"] + mat_table.size)

            material = {}
            for j, key in enumerate(mp_keys):
                ptr = mat_data_ptr[j + 1]
                if ptr:
                    if legacy:
                        # CPP reader does not work here since we are
                        # reading only a single value and not a table
                        # of data.
                        fs = (self._geometry_header["ptrMAT"] + ptr + 202) * 4
                        self._cfile._seekg(fs)
                        material[key] = np.fromstring(self._cfile._read(8))[0]
                    else:
                        material[key] = read_mat_data(ptr)

            # documentation for this is hard to follow, but it boils down to
            # the fact that "The record includes MP pointers followed by TB
            # pointers for a single Material ID"
            # So, after 80, these are all TB pointers, and fcCom.inc describes
            # maximum stresses:
            # co      XTEN,XCMP,YTEN,YCMP,ZTEN,ZCMP,XY,YZ,XZ,XYCP,YZCP,XZCP
            # co      XZIT,XZIC,YZIT,YZIC = Puck inclination parameters
            # only implemented for v194 and newer
            try:
                ptr = mat_data_ptr[163]
                if ptr:
                    arr = self.read_record(self._geometry_header["ptrMAT"] + ptr)
                    str_lim = {key: arr[index] for (index, key) in STR_LIM_REF.items()}
                    material["stress_failure_criteria"] = str_lim
            except:
                pass

            # store by material number
            self._materials[mat_data_ptr[0]] = material

    @property
    def materials(self):
        """Result file material properties.

        Returns
        -------
        dict
            Dictionary of Materials.  Keys are the material numbers,
            and each material is a dictionary of the material
            properrties of that material with only the valid entries filled.

        Notes
        -----
        Material properties:

        - EX : Elastic modulus, element x direction (Force/Area)
        - EY : Elastic modulus, element y direction (Force/Area)
        - EZ : Elastic modulus, element z direction (Force/Area)
        - ALPX : Coefficient of thermal expansion, element x direction (Strain/Temp)
        - ALPY : Coefficient of thermal expansion, element y direction (Strain/Temp)
        - ALPZ : Coefficient of thermal expansion, element z direction (Strain/Temp)
        - REFT : Reference temperature (as a property) [TREF]
        - PRXY : Major Poisson's ratio, x-y plane
        - PRYZ : Major Poisson's ratio, y-z plane
        - PRX  Z : Major Poisson's ratio, x-z plane
        - NUXY : Minor Poisson's ratio, x-y plane
        - NUYZ : Minor Poisson's ratio, y-z plane
        - NUXZ : Minor Poisson's ratio, x-z plane
        - GXY : Shear modulus, x-y plane (Force/Area)
        - GYZ : Shear modulus, y-z plane (Force/Area)
        - GXZ : Shear modulus, x-z plane (Force/Area)
        - DAMP : K matrix multiplier for damping [BETAD] (Time)
        - MU : Coefficient of friction (or, for FLUID29 and FLUID30
               elements, boundary admittance)
        - DENS : Mass density (Mass/Vol)
        - C : Specific heat (Heat/Mass*Temp)
        - ENTH : Enthalpy (e DENS*C d(Temp)) (Heat/Vol)
        - KXX : Thermal conductivity, element x direction
                (Heat*Length / (Time*Area*Temp))
        - KYY : Thermal conductivity, element y direction
                (Heat*Length / (Time*Area*Temp))
        - KZZ : Thermal conductivity, element z direction
                (Heat*Length / (Time*Area*Temp))
        - HF : Convection (or film) coefficient (Heat / (Time*Area*Temp))
        - EMIS : Emissivity
        - QRATE : Heat generation rate (MASS71 element only) (Heat/Time)
        - VISC : Viscosity (Force*Time / Length2)
        - SONC : Sonic velocity (FLUID29 and FLUID30 elements only) (Length/Time)
        - RSVX : Electrical resistivity, element x direction (Resistance*Area / Length)
        - RSVY : Electrical resistivity, element y direction (Resistance*Area / Length)
        - RSVZ : Electrical resistivity, element z direction (Resistance*Area / Length)
        - PERX : Electric permittivity, element x direction (Charge2 / (Force*Length))
        - PERY : Electric permittivity, element y direction (Charge2 / (Force*Length))
        - PERZ : Electric permittivity, element z direction (Charge2 / (Force*Length))
        - MURX : Magnetic relative permeability, element x direction
        - MURY : Magnetic relative permeability, element y direction
        - MURZ : Magnetic relative permeability, element z direction
        - MGXX : Magnetic coercive force, element x direction (Charge / (Length*Time))
        - MGYY : Magnetic coercive force, element y direction (Charge / (Length*Time))
        - MGZZ : Magnetic coercive force, element z direction (Charge / (Length*Time))

        Materials may contain the key ``"stress_failure_criteria"``, which
        contains failure criteria information for temperature-dependent stress
        limits. This includes the following keys:

        - XTEN : Allowable tensile stress or strain in the x-direction. (Must
          be positive.)

        - XCMP : Allowable compressive stress or strain in the
          x-direction. (Defaults to negative of XTEN.)

        - YTEN : Allowable tensile stress or strain in the y-direction. (Must
          be positive.)

        - YCMP : Allowable compressive stress or strain in the
          y-direction. (Defaults to negative of YTEN.)

        - ZTEN : Allowable tensile stress or strain in the z-direction. (Must
          be positive.)

        - ZCMP : Allowable compressive stress or strain in the
          z-direction. (Defaults to negative of ZTEN.)

        - XY : Allowable XY stress or shear strain. (Must be positive.)

        - YZ : Allowable YZ stress or shear strain. (Must be positive.)

        - XZ : Allowable XZ stress or shear strain. (Must be positive.)

        - XYCP : XY coupling coefficient (Used only if Lab1 = S). Defaults to -1.0. [1]

        - YZCP : YZ coupling coefficient (Used only if Lab1 = S). Defaults to -1.0. [1]

        - XZCP : XZ coupling coefficient (Used only if Lab1 = S). Defaults to -1.0. [1]

        - XZIT : XZ tensile inclination parameter for Puck failure index (default =
          0.0)

        - XZIC : XZ compressive inclination parameter for Puck failure index
          (default = 0.0)

        - YZIT : YZ tensile inclination parameter for Puck failure index
          (default = 0.0)

        - YZIC : YZ compressive inclination parameter for Puck failure index
          (default = 0.0)

        Examples
        --------
        Return the material properties from the example result
        file. Note that the keys of ``rst.materials`` is the material
        type.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> rst = pymapdl_reader.read_binary(examples.rstfile)
        >>> rst.materials
        {1: {'EX': 16900000.0, 'NUXY': 0.31, 'DENS': 0.00041408}}

        """
        if self._materials is None:
            self._load_materials()
        return self._materials

    def _load_section_data(self):
        """Loads the section data from the result file"""
        ptr_sec = self._geometry_header["ptrSEC"]
        sec_table, sz = self.read_record(ptr_sec, True)

        # Assemble section pointers relative to ptrSEC
        if self._map_flag:  # required as of 2021R1
            sec_pointers = self.read_record(ptr_sec + sz)
        else:
            sec_pointers = sec_table

        self._section_data = {}
        for offset in sec_pointers:
            if offset:
                table = self.read_record(ptr_sec + offset)
                self._section_data[int(table[0])] = table[1:]

        # it might be possible to interpret the section data...
        # sec[3] # total thickness
        # sec[4] # number of layers (?)
        # sec[5] # total integration points (?)

        # sec[29] # start of a layer (thickness)
        # - MatID
        # - Ori. Angle
        # - Num Intg.
        # sec[33] # start of another layer
        # sec[37] # start of another layer

    @property
    def section_data(self):
        """The section data from the result file

        Returns
        -------
        section_data : dict
            Dictionary of the section data with the section numbers as
            keys.

        Notes
        -----
        There is limited documentation on how ANSYS stores the
        sections within a result file, and as such it may be difficult
        to interpret the section data for a given model.
        """
        if self._section_data is None:
            self._load_section_data()
        return self._section_data

    def plot(
        self, node_components=None, element_components=None, sel_type_all=True, **kwargs
    ):
        """Plot result geometry

        Parameters
        ----------
        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``.

        Returns
        -------
        cpos : list
            List of camera position, focal point, and view up.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> rst.plot()

        Plot just the element component 'ROTOR_SHAFT'

        >>> rst.plot(element_components='ROTOR_SHAFT')

        Plot two node components
        >>> rst.plot(node_components=['MY_COMPONENT', 'MY_OTHER_COMPONENT'])
        """
        kwargs.setdefault("show_edges", True)

        if node_components:
            grid, _ = self._extract_node_components(node_components, sel_type_all)
        elif element_components:
            grid, _ = self._extract_element_components(element_components)
        else:
            grid = self.grid

        return self._plot_point_scalars(
            None,
            grid=grid,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            **kwargs,
        )

    def plot_nodal_solution(
        self,
        rnum,
        comp=None,
        show_displacement=False,
        displacement_factor=1.0,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plots the nodal solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Display component to display.  Options are ``'X'``,
            ``'Y'``, ``'Z'``, ``'NORM'``, or an available degree of
            freedom.  Result may also include other degrees of
            freedom, check ``rst.result_dof`` for available degrees of
            freedoms for a given result.  Defaults to ``"NORM"`` for a
            structural displacement result, and ``"TEMP"`` for a
            thermal result.

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``.

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot the nodal solution result 0 of verification manual
        example

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> result = examples.download_verification_result(33)
        >>> result.plot_nodal_solution(0)

        Plot with a white background and showing edges

        >>> result.plot_nodal_solution(0, background='w', show_edges=True)

        """
        rnum = self.parse_step_substep(rnum)
        nnum, result = self.nodal_solution(rnum)
        dof = self.result_dof(rnum)

        if self._is_thermal:
            label = "Temperature"
        else:
            label = "Displacement"

        if comp is None:
            if self._is_thermal:
                comp = "TEMP"
            else:
                comp = "NORM"
        if not isinstance(comp, str):
            raise TypeError("``comp`` must be a string")
        comp = comp.upper()

        # check component is in available degrees of freedom
        if comp in ("X", "Y", "Z"):
            comp = "U" + comp

        if comp == "NORM":
            if dof[:3] != ["UX", "UY", "UZ"] and dof[:2] != ["UX", "UY"]:
                raise AttributeError(
                    "Unable to compute norm given the DOF(s) %s" % str(dof)
                )
            # Normalize displacement
            scalars = np.linalg.norm(result[:, :3], axis=1)
            comp = "Normalized"
        else:
            if comp not in dof:
                raise ValueError(
                    "Component %s not in the available " % comp
                    + "degrees of freedom:\n``%s``" % dof
                )
            scalars = result[:, dof.index(comp)]

        # sometimes there are less nodes in the result than in the geometry
        npoints = self.grid.number_of_points
        if nnum.size != npoints:
            new_scalars = np.zeros(npoints)
            nnum_grid = self.grid.point_data["ansys_node_num"]
            new_scalars[np.in1d(nnum_grid, nnum, assume_unique=True)] = scalars
            scalars = new_scalars

        ind = None
        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            scalars = scalars[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            scalars = scalars[ind]
        else:
            grid = self.grid

        kwargs.setdefault("scalar_bar_args", {"title": f"{comp}\n{label}\n"})
        return self._plot_point_scalars(
            scalars,
            rnum=rnum,
            grid=grid,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            treat_nan_as_zero=treat_nan_as_zero,
            **kwargs,
        )

    @wraps(plot_nodal_solution)
    def plot_nodal_displacement(self, *args, **kwargs):
        """wraps plot_nodal_solution"""
        if self._is_thermal:
            raise AttributeError(
                "Thermal solution does not contain nodal " "displacement results"
            )
        return self.plot_nodal_solution(*args, **kwargs)

    @property
    def node_components(self):
        """Dictionary of ansys node components from the result file.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> rst = pymapdl_reader.read_binary(examples.rstfile)
        >>> rst.node_components.keys()
        dict_keys(['ECOMP1', 'ECOMP2', 'ELEM_COMP'])
        >>> rst.node_components['NODE_COMP']
        array([ 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
              20], dtype=int32)
        """
        return self._mesh.node_components

    @property
    def element_components(self):
        """Dictionary of ansys element components from the result file.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> rst = pymapdl_reader.read_binary(examples.rstfile)
        >>> rst.element_components
        {'ECOMP1': array([17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], dtype=int32),
        'ECOMP2': array([ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15, 16, 17, 18, 19, 20, 23, 24], dtype=int32),
        'ELEM_COMP': array([ 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20], dtype=int32)}
        """
        return self._mesh.element_components

    def _extract_node_components(self, node_components, sel_type_all=True, grid=None):
        """Return the part of the grid matching node components"""
        if grid is None:
            grid = self.grid

        if not self._mesh.node_components:  # pragma: no cover
            raise Exception(
                "Missing component information.\n"
                + "Either no components have been stored, or "
                + "the version of this result file is <18.2"
            )

        if isinstance(node_components, str):
            node_components = [node_components]

        mask = np.zeros(grid.n_points, np.bool_)
        for component in node_components:
            component = component.upper()
            if component not in grid.point_data:
                raise KeyError(
                    "Result file does not contain node " + 'component "%s"' % component
                )

            mask += grid.point_data[component].view(np.bool_)

        # need to extract the mesh
        cells, offset = vtk_cell_info(grid)
        if sel_type_all:
            cell_mask = cells_with_all_nodes(offset, cells, mask.view(np.uint8))
        else:
            cell_mask = cells_with_any_nodes(offset, cells, mask.view(np.uint8))

        if not cell_mask.any():
            raise RuntimeError("Empty component")
        reduced_grid = grid.extract_cells(cell_mask)

        if not reduced_grid.n_cells:
            raise RuntimeError(
                "Empty mesh due to component selection\n" + 'Try "sel_type_all=False"'
            )

        ind = reduced_grid.point_data["vtkOriginalPointIds"]
        if not ind.any():
            raise RuntimeError("Invalid selection.\n" + "Try ``sel_type_all=False``")

        return reduced_grid, ind

    def _extract_element_components(self, element_components, grid=None):
        """Return the part of the grid matching element components"""
        if grid is None:
            grid = self.grid

        if not self._mesh.element_components:  # pragma: no cover
            raise Exception(
                "Missing component information.\n"
                + "Either no components have been stored, or "
                + "the version of this result file is <18.2"
            )

        if isinstance(element_components, str):
            element_components = [element_components]

        cell_mask = np.zeros(grid.n_cells, np.bool_)
        for component in element_components:
            component = component.upper()
            if component not in grid.cell_data:
                raise KeyError(
                    "Result file does not contain element "
                    + 'component "%s"' % component
                )
            cell_mask += grid.cell_data[component].view(np.bool_)

        if not cell_mask.any():
            raise RuntimeError("Empty component")
        reduced_grid = grid.extract_cells(cell_mask)

        ind = reduced_grid.point_data["vtkOriginalPointIds"]
        if not ind.any():
            raise RuntimeError("Invalid component")

        return reduced_grid, ind

    @property
    def time_values(self):
        return self._resultheader["time_values"]

    def animate_nodal_solution(
        self,
        rnum,
        comp="norm",
        node_components=None,
        element_components=None,
        sel_type_all=True,
        add_text=True,
        displacement_factor=0.1,
        n_frames=100,
        loop=True,
        movie_filename=None,
        progress_bar=True,
        **kwargs,
    ):
        """Animate nodal solution.

        Assumes nodal solution is a displacement array from a modal or static
        solution.

        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, default: "norm"
            Scalar component to display.  Options are ``'x'``,
            ``'y'``, ``'z'``, and ``'norm'``, and ``None``.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        add_text : bool, optional
            Adds information about the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        n_frames : int, optional
            Number of "frames" between each full cycle.

        loop : bool, optional
            Loop the animation.  Default ``True``.  Disable this to
            animate once and close.  Automatically disabled when
            ``off_screen=True`` and ``movie_filename`` is set.

        movie_filename : str, pathlib.Path, optional
            Filename of the movie to open.  Filename should end in
            ``'mp4'``, but other filetypes may be supported like
            ``"gif"``.  See ``imagio.get_writer``.  A single loop of
            the mode will be recorded.

        progress_bar : bool, default: True
            Displays a progress bar when generating a movie while
            ``off_screen=True``.

        kwargs : optional keyword arguments
            See :func:`pyvista.plot` for additional keyword arguments.

        Examples
        --------
        Animate the first result interactively.

        >>> rst.animate_nodal_solution(0)

        Animate second result while displaying the x scalars
        without looping

        >>> rst.animate_nodal_solution(1, comp='x', loop=False)

        Animate the second result and save as a movie.

        >>> rst.animate_nodal_solution(0, movie_filename='disp.mp4')

        Animate the second result and save as a movie in the background.

        >>> rst.animate_nodal_solution(0, movie_filename='disp.mp4', off_screen=True)

        Disable plotting within the notebook.

        >>> rst.animate_nodal_solution(0, notebook=False)

        """
        if "nangles" in kwargs:  # pragma: no cover
            n_frames = kwargs.pop("nangles")
            warnings.warn(
                "The ``nangles`` kwarg is depreciated and ``n_frames`` "
                "should be used instead."
            )

        scalars = None
        if comp:
            _, disp = self.nodal_solution(rnum)
            disp = disp[:, :3]

            if comp == "x":
                axis = 0
            elif comp == "y":
                axis = 1
            elif comp == "z":
                axis = 2
            else:
                axis = None

            if axis is not None:
                scalars = disp[:, axis]
            else:
                scalars = (disp * disp).sum(1) ** 0.5

        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            if comp:
                scalars = scalars[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            if comp:
                scalars = scalars[ind]
        else:
            grid = self.grid

        return self._plot_point_scalars(
            scalars,
            rnum=rnum,
            grid=grid,
            add_text=add_text,
            animate=True,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            n_frames=n_frames,
            displacement_factor=displacement_factor,
            movie_filename=movie_filename,
            loop=loop,
            progress_bar=progress_bar,
            **kwargs,
        )

    @wraps(animate_nodal_solution)
    def animate_nodal_displacement(self, *args, **kwargs):
        """wraps animate_nodal_solution"""
        return self.animate_nodal_solution(*args, **kwargs)

    def animate_nodal_solution_set(
        self,
        rnums=None,
        comp="norm",
        node_components=None,
        element_components=None,
        sel_type_all=True,
        loop=True,
        movie_filename=None,
        add_text=True,
        fps=20,
        **kwargs,
    ):
        """Animate a set of nodal solutions.

        Animates the scalars of all the result sets.  Best when used
        with a series of static analyses.

        rnums : collection.Iterable
            Range or list containing the zero based indexed cumulative
            result numbers to animate.

        comp : str, optional
            Scalar component to display.  Options are ``'x'``,
            ``'y'``, ``'z'``, and ``'norm'``, and ``None``.  Not
            applicable for a thermal analysis.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        loop : bool, optional
            Loop the animation.  Default ``True``.  Disable this to
            animate once and close.

        movie_filename : str, optional
            Filename of the movie to open.  Filename should end in
            ``'mp4'``, but other filetypes may be supported.  See
            ``imagio.get_writer``.  A single loop of the mode will be
            recorded.

        add_text : bool, optional
            Adds information about the result to the animation.

        fps : int, optional
            Frames per second.  Defaults to 20 and limited to hardware
            capabilities and model density. Carries over to movies
            created by providing the ``movie_filename`` argument,
            but *not* to gifs.

        kwargs : optional keyword arguments, optional
            See help(pyvista.Plot) for additional keyword arguments.

        Examples
        --------
        Animate all results

        >>> rst.animate_nodal_solution_set()

        Animate every 50th result in a set of results and save to a
        gif.  Use the "zx" camera position to view the ZX plane from
        the top down.

        >>> rsets = range(0, rst.nsets, 50)
        >>> rst.animate_nodal_solution_set(rsets,
        ...                                scalar_bar_args={'title': 'My Animation'},
        ...                                lighting=False, cpos='zx',
        ...                                movie_filename='example.gif')
        """
        if self.nsets == 1:
            raise RuntimeError(
                "Unable to animate.  Solution contains only one " "result set"
            )

        if rnums is None:
            rnums = range(self.nsets)
        elif not isinstance(rnums, Iterable):
            raise TypeError("``rnums`` must be a range or list of result numbers")

        def get_scalars(rnum):
            _, data = self.nodal_solution(rnum)
            if not self._is_thermal:
                data = data[:, :3]

                axis = None
                if comp == "x":
                    axis = 0
                elif comp == "y":
                    axis = 1
                elif comp == "z":
                    axis = 2

                if axis is not None:
                    data = data[:, axis]
                else:
                    data = (data * data).sum(1) ** 0.5
            else:
                data = data.ravel()

            return data

        ind = None
        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
        else:
            grid = self.grid

        scalars = []
        for rnum in tqdm(rnums, desc="Caching scalars"):
            if ind is None:
                scalars.append(get_scalars(rnum))
            else:
                scalars.append(get_scalars(rnum)[ind])

        text = None
        if add_text:
            time_values = self.time_values
            text = ["Time Value: %12.4f" % time_values[i] for i in rnums]

        return self._animate_point_scalars(
            scalars,
            grid=grid,
            text=text,
            movie_filename=movie_filename,
            loop=loop,
            fps=fps,
            **kwargs,
        )

    def nodal_time_history(self, solution_type="NSL", in_nodal_coord_sys=False):
        """The DOF solution for each node for all result sets.

        The nodal results are returned returned in the global
        cartesian coordinate system or nodal coordinate system.

        Parameters
        ----------
        solution_type: str, optional
            The solution type.  Must be either nodal displacements
            (``'NSL'``), nodal velocities (``'VEL'``) or nodal
            accelerations (``'ACC'``).

        in_nodal_coord_sys : bool, optional
            When ``True``, returns results in the nodal coordinate system.
            Default ``False``.

        Returns
        -------
        nnum : int np.ndarray
            Node numbers associated with the results.

        result : float np.ndarray
            Nodal solution for all result sets.  Array is sized
            ``rst.nsets x nnod x Sumdof``, which is the number of
            time steps by number of nodes by degrees of freedom.
        """
        if not isinstance(solution_type, str):
            raise TypeError("Solution type must be a string")

        if solution_type == "NSL":
            func = self.nodal_solution
        elif solution_type == "VEL":
            func = self.nodal_velocity
        elif solution_type == "ACC":
            func = self.nodal_acceleration
        else:
            raise ValueError(
                "Argument 'solution type' must be either 'NSL', " "'VEL', or 'ACC'"
            )

        # size based on the first result
        nnum, sol = func(0)
        data = np.empty((self.nsets, sol.shape[0], sol.shape[1]), np.float64)
        for i in range(self.nsets):
            data[i] = func(i)[1]

        return nnum, data

    def nodal_solution(self, rnum, in_nodal_coord_sys=False, nodes=None):
        """Returns the DOF solution for each node in the global
        cartesian coordinate system or nodal coordinate system.

        Solution may be nodal temperatures or nodal displacements
        depending on the type of the solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        in_nodal_coord_sys : bool, optional
            When ``True``, returns results in the nodal coordinate
            system.  Default ``False``.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : int np.ndarray
            Node numbers associated with the results.

        result : float np.ndarray
            Array of nodal displacements or nodal temperatures.  Array
            is (``nnod`` x ``sumdof``), the number of nodes by the
            number of degrees of freedom which includes ``numdof`` and
            ``nfldof``

        Examples
        --------
        Return the nodal solution (in this case, displacement) for the
        first result of ``"file.rst"``.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, data = rst.nodal_solution(0)

        Return the nodal solution just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, data = rst.nodal_solution(0, nodes='MY_COMPONENT')

        Return the nodal solution just for the nodes from 20 through 50.

        >>> nnum, data = rst.nodal_solution(0, nodes=range(20, 51))

        Notes
        -----
        Some solution results may not include results for each node.
        These results are removed by and the node numbers of the
        solution results are reflected in ``nnum``.
        """
        return self._nodal_solution_result(rnum, "NSL", in_nodal_coord_sys, nodes)

    def nodal_velocity(self, rnum, in_nodal_coord_sys=False):
        """Nodal velocities for a given result set.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        in_nodal_coord_sys : bool, optional
            When ``True``, returns results in the nodal coordinate
            system.  Default False.

        Returns
        -------
        nnum : int np.ndarray
            Node numbers associated with the results.

        result : float np.ndarray
            Array of nodal velocities.  Array is (``nnod`` x
            ``sumdof``), the number of nodes by the number of degrees
            of freedom which includes ``numdof`` and ``nfldof``

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, data = rst.nodal_velocity(0)

        Notes
        -----
        Some solution results may not include results for each node.
        These results are removed by and the node numbers of the
        solution results are reflected in ``nnum``.
        """
        return self._nodal_solution_result(rnum, "VSL", in_nodal_coord_sys)

    def nodal_acceleration(self, rnum, in_nodal_coord_sys=False):
        """Nodal velocities for a given result set.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        in_nodal_coord_sys : bool, optional
            When ``True``, returns results in the nodal coordinate
            system.  Default False.

        Returns
        -------
        nnum : int np.ndarray
            Node numbers associated with the results.

        result : float np.ndarray
            Array of nodal accelerations.  Array is (``nnod`` x
            ``sumdof``), the number of nodes by the number of degrees
            of freedom which includes ``numdof`` and ``nfldof``

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, data = rst.nodal_acceleration(0)

        Notes
        -----
        Some solution results may not include results for each node.
        These results are removed by and the node numbers of the
        solution results are reflected in ``nnum``.
        """
        return self._nodal_solution_result(rnum, "ASL", in_nodal_coord_sys)

    def _nodal_solution_result(
        self, rnum, solution_type, in_nodal_coord_sys=False, nodes=None
    ):
        """Returns the DOF solution for each node in the global
        cartesian coordinate system or nodal coordinate system.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        solution_type: str, optional
            Specify, whether nodal displacements (``'NSL'``), nodal
            velocities (``'VEL'``) or nodal accelerations (``'ACC'``)
            will be read.

        in_nodal_coord_sys : bool, optional
            When True, returns results in the nodal coordinate system.
            Default False.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : int np.ndarray
            Node numbers associated with the results.

        result : float np.ndarray
            Result is (``nnod`` x ``sumdof``), or number of nodes by degrees
            of freedom which includes ``numdof`` and ``nfldof``.

        Notes
        -----
        Some solution results may not include results for each node.
        These results are removed by and the node numbers of the
        solution results are reflected in ``nnum``.
        """
        if solution_type not in ("NSL", "VSL", "ASL"):  # pragma: no cover
            raise ValueError(
                "Argument ``solution type`` must be either " "'NSL', 'VSL', or 'ASL'"
            )

        # check if nodal solution exists
        if not self.available_results[solution_type]:
            raise AttributeError(
                'Result file is missing "%s"'
                % self.available_results.description[solution_type]
            )

        # result pointer
        rnum = self.parse_step_substep(rnum)  # convert to cumulative index
        ptr_rst = self._resultheader["rpointers"][rnum]
        result_solution_header = self._result_solution_header(rnum)

        nnod = result_solution_header["nnod"]
        numdof = result_solution_header["numdof"]
        nfldof = result_solution_header["nfldof"]
        sumdof = numdof + nfldof

        # # Size of the results is dependent on result type
        if solution_type == "NSL":
            dof = sumdof  # NSL : nnod*Sumdof
        elif solution_type == "VSL":  # for transient
            # nnod*numvdof
            dof = self._result_solution_header_ext(rnum)[126]
        elif solution_type == "ASL":
            # nnod*numadof
            dof = self._result_solution_header_ext(rnum)[127]

        ptr = result_solution_header["ptr" + solution_type]
        if ptr == 0:  # shouldn't get here...
            raise AttributeError(
                'Result file is missing "%s"'
                % self.available_results.description[solution_type]
            )

        # Read the nodal results
        result, bufsz = self.read_record(ptr + ptr_rst, True)
        result = result.reshape(-1, dof)

        # additional entries are sometimes written for no discernible
        # reason
        if result.shape[0] > nnod:
            result = result[:nnod]

        # it's possible that not all results are written
        if result.shape[0] < nnod:
            # read second buffer containing the node indices of the
            # results and convert from fortran to zero indexing
            sidx = self.read_record(ptr + ptr_rst + bufsz) - 1
            unsort_nnum = self._resultheader["neqv"][sidx]

            # now, sort using the new sorted node numbers indices
            new_sidx = np.argsort(unsort_nnum)
            nnum = unsort_nnum[new_sidx]
            result = result[new_sidx]
        else:
            nnum = self._neqv[self._sidx]
            result = result.take(self._sidx, 0)

        # Convert result to the global coordinate system
        if not in_nodal_coord_sys:
            euler_angles = self._mesh.node_angles[self._insolution].T
            rotate_to_global(result, euler_angles)

        # check for invalid values (mapdl writes invalid values as 2*100)
        result[result == 2**100] = 0

        # subselect from component or range
        # TODO: inefficient and should be able to only read from this subselection.
        if nodes is not None:
            if isinstance(nodes, str):
                if nodes not in self.grid.point_data:
                    raise ValueError(f"Invalid node component {nodes}")
                mask = self.grid.point_data[nodes].view(bool)
            elif isinstance(nodes, Sequence):
                if isinstance(nodes[0], int):
                    mask = np.in1d(self.mesh.nnum, nodes)
                elif isinstance(nodes[0], str):
                    mask = np.zeros(self.grid.n_points, bool)
                    for node_component in nodes:
                        if node_component not in self.grid.point_data:
                            raise ValueError(f"Invalid node component {nodes}")
                        mask = np.logical_or(
                            mask, self.grid.point_data[node_component].view(bool)
                        )
                else:
                    raise TypeError(
                        f"Invalid type for {nodes}.  Expected Sequence of " "str or int"
                    )
            else:
                raise TypeError(
                    f"Invalid type for {nodes}.  Expected Sequence of " "str or int"
                )

            # mask is for global nodes, need for potentially subselectd nodes
            submask = np.in1d(nnum, self.grid.point_data["ansys_node_num"][mask])
            nnum, result = nnum[submask], result[submask]

        return nnum, result

    @wraps(nodal_solution)
    def nodal_displacement(self, *args, **kwargs):
        """wraps plot_nodal_solution"""
        if self._is_thermal:
            raise AttributeError(
                "Thermal solution does not contain nodal " "displacement results"
            )
        return self.nodal_solution(*args, **kwargs)

    def _read_components(self):
        """Read components from an ANSYS result file

        Returns
        -------
        components : dict
            Dictionary of components

        Notes
        -----
        Only available as of ~v194 and newer.

        From fdresu.inc:

        Component/assembly data.  The first word will describe the
        type of component or assembly (1=node component, 2=elem
        component, 11->15 assembly) For components, values 2->9 are
        the component name converted to integers.  The remaining
        values are the list of nodes/elems in the component.  For
        assemblies, the remaining values (in groups of 8) are the
        component names converted to integers.
        """
        node_comp = {}
        elem_comp = {}
        ncomp = self._geometry_header["maxcomp"]

        # Read through components
        file_ptr = self._geometry_header["ptrCOMP"]
        for _ in range(ncomp):
            table, sz = self.read_record(file_ptr, True)
            file_ptr += sz  # increment file_pointer
            name = table[1:9].tobytes().split(b"\x00")[0].decode("latin")
            name = (
                name[:4][::-1]
                + name[4:8][::-1]
                + name[8:12][::-1]
                + name[12:16][::-1]
                + name[16:20][::-1]
                + name[20:24][::-1]
                + name[24:28][::-1]
                + name[28:32][::-1]
            )
            name = name.strip()
            data = table[9:]
            if data.any():
                if table[0] == 1:  # node component
                    node_comp[name] = _reader.component_interperter(data)
                elif table[0] == 2:
                    elem_comp[name] = _reader.component_interperter(data)

        return node_comp, elem_comp

    @property
    def _map_flag(self):
        """Flag to indicate format of mapping index vectors for
        element types, real constants, coordinate systems, and sections.

        When ``True``, uses two vectors, and ``False`` using the old
        format (prior to 2021R1).

        """
        return bool(self._geometry_header["mapFlag"])

    def _load_element_table(self):
        """Store element type information"""
        # pointer to the element type index table
        ptrety = self._geometry_header["ptrETY"]
        e_type_table, sz = self.read_record(ptrety, True)

        # store information for each element type
        nodelm = {}  # n nodes for this element type
        nodfor = {}  # n nodes per element having nodal forces
        nodstr = {}  # n nodes per element having nodal stresses
        ekey = []
        keyopts = {}  # key options

        # Assemble element record pointers relative to ptrETY
        if self._map_flag:
            elem_record_pointers = self.read_record(ptrety + sz)
        else:
            elem_record_pointers = []
            for i in range(self._geometry_header["maxety"]):
                if e_type_table[i]:
                    elem_record_pointers.append(e_type_table[i])

        for elem_record_pointer in elem_record_pointers:
            ptr = ptrety + elem_record_pointer
            einfo = self.read_record(ptr)

            etype_ref = einfo[0]
            ekey.append(einfo[:2])

            # Items 3-14 - element type option keys (keyopts)
            keyopts[etype_ref] = einfo[2:13]

            # Item 61 - number of nodes for this element type (nodelm)
            nodelm[etype_ref] = einfo[60]

            # Item 63 - number of nodes per element having nodal
            # forces, etc. (nodfor)
            nodfor[etype_ref] = einfo[62]

            # Item 94 - number of nodes per element having nodal
            # stresses, etc. (nodstr).  This number is the number
            # of corner nodes for higher-ordered elements.
            nodstr[etype_ref] = einfo[93]

            # with KEYOPT(8)=0, the record contains stresses at
            # each corner node (first at the bottom shell surface,
            # then the top surface)
            #
            # Only valid for SHELL181 or SHELL281 elements.
            if einfo[1] == 181 or einfo[1] == 281:
                if keyopts[etype_ref][7] == 0:
                    nodstr[etype_ref] *= 2

        def dict_to_arr(index_dict):
            """Convert an index dictionary to an array

            For example:
            {1: 20, 2: 30} --> [UNDEF, 20, 30]

            """
            arr = np.empty(max(index_dict.keys()) + 1, np.int32)
            for key, value in index_dict.items():
                arr[key] = value
            return arr

        # rest of pymapdl-reader expects the following as arrays.
        # store element table data
        return {
            "nodelm": dict_to_arr(nodelm),
            "nodfor": dict_to_arr(nodfor),
            "nodstr": dict_to_arr(nodstr),
            "keyopts": keyopts,
            "ekey": np.array(ekey),
        }

    def _store_mesh(self, parse_vtk=True):
        """Store the mesh from the result file.

        Parameters
        ----------
        parse_vtk : bool, optional
            Set to ``False`` to skip the parsing the mesh as a VTK
            UnstructuredGrid, which might take a long time for large models.
        """
        # Node information
        nnod = self._geometry_header["nnod"]
        nnum = np.empty(nnod, np.int32)
        nodes = np.empty((nnod, 6), np.float64)
        _binary_reader.load_nodes(
            self.filename, self._geometry_header["ptrLOC"], nnod, nodes, nnum
        )

        # the element description table
        # must view this record as int64, even though ansys reads
        # it in as a 32 bit int
        ptr_eid = self._geometry_header["ptrEID"]
        e_disp_table = self.read_record(ptr_eid).view(np.int64)

        # get pointer to start of element table and adjust element pointers
        ptr_elem = self._geometry_header["ptrEID"] + e_disp_table[0]
        e_disp_table -= e_disp_table[0]

        # read in coordinate systems, material properties, and sections
        self._c_systems = self.parse_coordinate_system()

        # load elements
        nelm = self._geometry_header["nelm"]
        elem, elem_off = _binary_reader.load_elements(
            self.filename, ptr_elem, nelm, e_disp_table
        )

        # Store geometry and parse to VTK quadradic and null unallowed
        ncomp, ecomp = self._read_components()
        self._mesh = Mesh(
            nnum,
            nodes,
            elem,
            elem_off,
            self._element_table["ekey"],
            node_comps=ncomp,
            elem_comps=ecomp,
        )
        if parse_vtk:
            self.quadgrid = self._mesh._parse_vtk(
                null_unallowed=True, fix_midside=False
            )
            self.grid = self.quadgrid.linear_copy()
        else:
            self.quadgrid = None
            self.grid = None

        # identify nodes that are actually in the solution
        self._insolution = np.in1d(
            self._mesh.nnum, self._resultheader["neqv"], assume_unique=True
        )

    def solution_info(self, rnum):
        """Return an informative dictionary of solution data for a
        result.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        header : dict
            Double precision solution header data.

        Examples
        --------
        Extract the solution info from a sample example result file.

        >>> from ansys.mapdl.reader import examples
        >>> rst = examples.download_pontoon()
        >>> rst.solution_info(0)
        {'cgcent': [],
         'fatjack': [],
         'timfrq': 44.85185724963714,
         'lfacto': 1.0,
         'lfactn': 1.0,
         'cptime': 3586.4873046875,
         'tref': 71.6,
         'tunif': 71.6,
         'tbulk': 293.0,
         'volbase': 0.0,
         'tstep': 0.0,
         '__unused': 0.0,
         'accel_x': 0.0,
         'accel_y': 0.0,
         'accel_z': 0.0,
         'omega_v_x': 0.0,
         'omega_v_y': 0.0,
         'omega_v_z': 0.0,
         'omega_a_x': 0.0,
         'omega_a_y': 0.0,
         'omega_a_z': 0.0,
         'omegacg_v_x': 0.0,
         'omegacg_v_y': 0.0,
         'omegacg_v_z': 0.0,
         'omegacg_a_x': 0.0,
         'omegacg_a_y': 0.0,
         'omegacg_a_z': 0.0,
         'dval1': 0.0,
         'pCnvVal': 0.0}


        Notes
        -----
        The keys of the solution header are described below:

        - timfrq : Time value (or frequency value, for a modal or
                   harmonic analysis)
        - lfacto : the "old" load factor (used in ramping a load
                    between old and new values)
        - lfactn  : The "new" load factor
        - cptime  : Elapsed CPU time (in seconds)
        - tref    : The reference temperature
        - tunif   : The uniform temperature
        - tbulk   : Bulk temp for FLOTRAN film coefs.
        - VolBase : Initial total volume for VOF
        - tstep   : Time Step size for FLOTRAN analysis
        - 0.0     : Position not used
        - accel   : Linear acceleration terms
        - omega   : Angular velocity (first 3 terms) and angular acceleration
                    (second 3 terms)
        - omegacg : Angular velocity (first 3 terms) and angular
                    acceleration (second 3 terms) these
                    velocity/acceleration terms are computed about the
                    center of gravity
        - cgcent  : (X,y,z) location of center of gravity
        - fatjack : Fatjack ocean wave data (wave height and period)
        - dval1   : If pmeth=0: FATJACK ocean wave direction
                    if pmeth=1: p-method convergence values
        - pCnvVal : P-method convergence values
        """
        # convert to cumulative index
        rnum = self.parse_step_substep(rnum)

        # skip solution data header
        ptr = self._resultheader["rpointers"][rnum]
        _, sz = self.read_record(ptr, True)

        table = self.read_record(ptr + sz)
        return parse_header(table, solution_header_keys_dp)

    def _element_solution_header(self, rnum):
        """Return the element solution header information.

        Parameters
        ----------
        rnum : int
            Result number.

        Returns
        -------
        ele_ind_table :
            A table of pointers relative to the element result table
            offset containing element result pointers.

        etype : np.ndarray
            Element type array

        element_rst_ptr : int
            Offset to the start of the element result table.
        """
        # Get the header information from the header dictionary
        rpointers = self._resultheader["rpointers"]

        # read result solution header
        solution_header = self._result_solution_header(rnum)

        # key to extrapolate integration
        # = 0 - move
        # = 1 - extrapolate unless active
        # non-linear
        # = 2 - extrapolate always
        if solution_header["rxtrap"] == 0:
            warnings.warn(
                "Strains and stresses are being evaluated at "
                "gauss points and not extrapolated"
            )

        # 64-bit pointer to element solution
        if not solution_header["ptrESL"]:
            raise ValueError(
                "No element solution in result set %d\n" % (rnum + 1)
                + 'Try running with "MXPAND,,,,YES"'
            )

        # Seek to element result header
        element_rst_ptr = rpointers[rnum] + solution_header["ptrESL"]
        ele_ind_table = self.read_record(element_rst_ptr).view(np.int64)
        # ele_ind_table += element_rst_ptr

        return ele_ind_table, self._mesh._ans_etype, element_rst_ptr

    def result_dof(self, rnum):
        """Return a list of degrees of freedom for a given result number.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        dof : list
            List of degrees of freedom.

        Examples
        --------
        >>> rst.result_dof(0)
        ['UX', 'UY', 'UZ']
        """
        rnum = self.parse_step_substep(rnum)
        return [DOF_REF[dof_idx] for dof_idx in self._solution_header(rnum)["DOFS"]]

    def nodal_boundary_conditions(self, rnum):
        """Nodal boundary conditions for a given result number.

        These nodal boundary conditions are generally set with the
        APDL command ``D``.  For example, ``D, 25, UX, 0.001``

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        nnum : np.ndarray
            Node numbers of the nodes with boundary conditions.

        dof : np.ndarray
            Array of indices of the degrees of freedom of the nodes
            with boundary conditions.  See ``rst.result_dof`` for the
            degrees of freedom associated with each index.

        bc : np.ndarray
            Boundary conditions.

        Examples
        --------
        Print the boundary conditions where:
        - Node 3 is fixed
        - Node 25 has UX=0.001
        - Node 26 has UY=0.0011
        - Node 27 has UZ=0.0012

        >>> rst.nodal_boundary_conditions(0)
        (array([ 3,  3,  3, 25, 26, 27], dtype=int32),
        array([1, 2, 3, 1, 2, 3], dtype=int32),
        array([0.    , 0.    , 0.    , 0.001 , 0.0011, 0.0012]))
        """
        bc_ptr, bc_header = self._bc_header(rnum)

        if not bc_header["ptrDIX"]:
            raise IndexError("This result contains no nodal input forces")

        # nodes with boundary conditions
        ptr_dix = bc_ptr + bc_header["ptrDIX"]
        bc_nnum, sz = self.read_record(ptr_dix, return_bufsize=True)
        if bc_header["format"]:
            # there is a record immediately following containing the DOF data
            dof = self.read_record(ptr_dix + sz)
        else:
            raise NotImplementedError('Format "0" is unhandled.\n%s' % EMAIL_ME)

        # RECORD : DIS     dp       1    4*numdis
        # Nodal constraints.  This record contains present and
        # previous values (real and imaginary) of the nodal
        # constraints at each DOF.
        #
        # Repeated 4 times for real/imag for present/prior
        # numdis : number of nodal constraints
        bc_record = self.read_record(bc_ptr + bc_header["ptrDIS"])
        bc = bc_record.reshape(bc_nnum.size, -1)[:, 0]  # use only present values

        return bc_nnum, dof, bc

    def nodal_input_force(self, rnum):
        """Nodal input force for a given result number.

        Nodal input force is generally set with the APDL command
        ``F``.  For example, ``F, 25, FX, 0.001``

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        nnum : np.ndarray
            Node numbers of the nodes with nodal forces.

        dof : np.ndarray
            Array of indices of the degrees of freedom of the nodes
            with input force.  See ``rst.result_dof`` for the degrees
            of freedom associated with each index.

        force : np.ndarray
            Nodal input force.

        Examples
        --------
        Print the nodal input force where:
        - Node 25 has FX=20
        - Node 26 has FY=30
        - Node 27 has FZ=40

        >>> rst.nodal_input_force(0)
        (array([ 71,  52, 127], dtype=int32),
         array([2, 1, 3], dtype=int32),
         array([30., 20., 40.]))
        """
        bc_ptr, bc_header = self._bc_header(rnum)

        if not bc_header["ptrFIX"]:
            raise IndexError("This result contains no nodal input forces")

        # nodes with input force
        ptr_fix = bc_ptr + bc_header["ptrFIX"]
        f_nnum, sz = self.read_record(ptr_fix, return_bufsize=True)
        if bc_header["format"]:
            # there is a record immediately following containing the DOF data
            dof = self.read_record(ptr_fix + sz)
        else:
            raise NotImplementedError('Format "0" is unhandled.\n%s' % EMAIL_ME)

        # Repeated 4 times for real/imag for present/prior
        # numdis : number of nodal constraints
        f_record = self.read_record(bc_ptr + bc_header["ptrFOR"])
        force = f_record.reshape(f_nnum.size, -1)[:, 0]  # use only present values
        return f_nnum, dof, force

    def _bc_header(self, rnum):
        """Return the boundary conditions ptr and header for a given result number"""
        rnum = self.parse_step_substep(rnum)
        rptr = self._result_pointers[rnum]

        # read bc_header
        bc_ptr = rptr + self._solution_header(rnum)["ptrBC"]
        bc_header = parse_header(
            self.read_record(bc_ptr), boundary_condition_index_table
        )

        return bc_ptr, bc_header

    def _solution_header(self, rnum):
        """The solution header for a given cumulative result"""
        record = self.read_record(self._result_pointers[rnum])
        return parse_header(record, solution_data_header_keys)

    @property
    def _result_pointers(self):
        """Array of result pointers (position within file)"""
        return self._resultheader["rpointers"]

    @property
    def version(self):
        """The version of MAPDL used to generate this result file.

        Examples
        --------
        >>> rst.version
        20.1
        """
        return float(self._resultheader["verstring"])

    def element_stress(
        self, rnum, principal=False, in_element_coord_sys=False, **kwargs
    ):
        """Retrieves the element component stresses.

        Equivalent ANSYS command: PRESOL, S

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        principal : bool, optional
            Returns principal stresses instead of component stresses.
            Default False.

        in_element_coord_sys : bool, optional
            Returns the results in the element coordinate system.
            Default False and will return the results in the global
            coordinate system.

        **kwargs : optional keyword arguments
            Hidden options for distributed result files.

        Returns
        -------
        enum : np.ndarray
            ANSYS element numbers corresponding to each element.

        element_stress : list
            Stresses at each element for each node for Sx Sy Sz Sxy
            Syz Sxz or SIGMA1, SIGMA2, SIGMA3, SINT, SEQV when
            principal is True.

        enode : list
            Node numbers corresponding to each element's stress
            results.  One list entry for each element.

        Examples
        --------
        Element component stress for the first result set.

        >>> rst.element_stress(0)

        Element principal stress for the first result set.

        >>> enum, element_stress, enode = result.element_stress(0, principal=True)

        Notes
        -----
        Shell stresses for element 181 are returned for top and bottom
        layers.  Results are ordered such that the top layer and then
        the bottom layer is reported.
        """
        rnum = self.parse_step_substep(rnum)
        ele_ind_table, etype, ptr_off = self._element_solution_header(rnum)

        # certain element types do not output stress
        elemtype = self._mesh.etype

        # load in raw results
        nnode = self._nodstr[etype]
        nelemnode = nnode.sum()

        # bitmask (might use this at some point)
        # bitmask = bin(int(hex(self._resultheader['rstsprs']), base=16)).lstrip('0b')
        # description maybe in resucm.inc

        if self.version >= 14.5:
            if self._resultheader["rstsprs"] != 0:
                nitem = 6
            else:
                nitem = 11

            # add extra elements to data array.  Sometimes there are
            # more items than listed in the result header (or there's
            # a mistake here)
            ele_data_arr = np.empty((nelemnode + 50, nitem), np.float64)
            ele_data_arr[:] = np.nan  # necessary?  should do this in read stress

            _binary_reader.read_element_stress(
                self.filename,
                ele_ind_table,
                self._nodstr.astype(np.int64),
                etype,
                ele_data_arr,
                nitem,
                elemtype,
                ptr_off,
                as_global=not in_element_coord_sys,
            )

            if nitem != 6:
                ele_data_arr = ele_data_arr[:, :6]

        else:
            raise NotImplementedError("Not implemented for ANSYS older than v14.5")

        # trim off extra data
        ele_data_arr = ele_data_arr[:nelemnode]

        if principal:
            ele_data_arr, isnan = _binary_reader.compute_principal_stress(ele_data_arr)
            ele_data_arr[isnan] = np.nan

        splitind = np.cumsum(nnode)
        element_stress = np.split(ele_data_arr, splitind[:-1])

        # just return the distributed result if requested
        if kwargs.get("is_dist_rst", False):
            return element_stress

        # reorder list using sorted indices
        enum = self._eeqv
        sidx = np.argsort(enum)
        element_stress = [element_stress[i] for i in sidx]

        enode = []
        for i in sidx:
            enode.append(self._mesh.elem[i][10 : 10 + nnode[i]])

        # Get element numbers
        elemnum = self._eeqv[self._sidx_elem]
        return elemnum, element_stress, enode

    @property
    def _element_map(self):
        if self.__element_map is None:
            self.__element_map = dict(zip(self.mesh.enum, np.arange(self.mesh.n_elem)))
        return self.__element_map

    def element_lookup(self, element_id):
        """Index of the element the element within the result mesh"""
        try:
            return self._element_map[element_id]
        except KeyError:
            raise KeyError(
                f"Element ID {element_id} not in this result file."
            ) from None

    def overwrite_element_solution_record(self, data, rnum, solution_type, element_id):
        """Overwrite element solution record.

        This method replaces solution data for of an element at a
        result index for a given solution type.  The number of items
        in ``data`` must match the number of items in the record.

        If you are not sure how many records are in a given record,
        use ``element_solution_data`` to retrieve all the records for
        a given ``solution_type`` and check the number of items in the
        record.

        Note: The record being replaced cannot be a compressed record.
        If the result file uses compression (default sparse
        compression as of 2019R1), you can disable this within MAPDL
        with:
        ``/FCOMP, RST, 0``

        Parameters
        ----------
        data : list or np.ndarray
            Data that will replace the existing records.

        rnum : int
            Zero based result number.

        solution_type : str
            Element data type to overwrite.

            - EMS: misc. data
            - ENF: nodal forces
            - ENS: nodal stresses
            - ENG: volume and energies
            - EGR: nodal gradients
            - EEL: elastic strains
            - EPL: plastic strains
            - ECR: creep strains
            - ETH: thermal strains
            - EUL: euler angles
            - EFX: nodal fluxes
            - ELF: local forces
            - EMN: misc. non-sum values
            - ECD: element current densities
            - ENL: nodal nonlinear data
            - EHC: calculated heat generations
            - EPT: element temperatures
            - ESF: element surface stresses
            - EDI: diffusion strains
            - ETB: ETABLE items
            - ECT: contact data
            - EXY: integration point locations
            - EBA: back stresses
            - ESV: state variables
            - MNL: material nonlinear record

        element_id : int
            Ansys element number (e.g. ``1``)

        Examples
        --------
        Overwrite the elastic strain record for element 1 for the
        first result with random data.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> data = np.random.random(56)
        >>> rst.overwrite_element_solution_data(data, 0, 'EEL', 1)
        """
        if not isinstance(data, np.ndarray):
            data = np.asarray(data)

        table_ptr = solution_type.upper()
        if table_ptr not in ELEMENT_INDEX_TABLE_KEYS:
            err_str = "Data type %s is invalid\n" % str(solution_type)
            err_str += "\nAvailable types:\n"
            for key in ELEMENT_INDEX_TABLE_KEYS:
                err_str += "\t%s: %s\n" % (key, element_index_table_info[key])
            raise ValueError(err_str)

        if table_ptr in self.available_results._parsed_bits:
            if table_ptr not in self.available_results:
                raise ValueError(
                    "Result %s is not available in this result file" % table_ptr
                )

        # location of data pointer within each element result table
        table_index = ELEMENT_INDEX_TABLE_KEYS.index(table_ptr)

        rnum = self.parse_step_substep(rnum)
        ele_ind_table, etype, ptr_off = self._element_solution_header(rnum)

        # index of the element in the element table
        idx = self.element_lookup(element_id)
        elem_ptr = ele_ind_table[idx]
        if elem_ptr == 0:
            raise ValueError(
                "This element does not have any data associated with the "
                f" solution_type {solution_type}"
            )

        # new c stuff here
        if data.dtype == np.float32:
            self._cfile.overwrite_element_data_float(
                elem_ptr + ptr_off, table_index, data
            )
        elif data.dtype == np.double:
            self._cfile.overwrite_element_data_double(
                elem_ptr + ptr_off, table_index, data
            )
        else:
            self._cfile.overwrite_element_data_double(
                elem_ptr + ptr_off, table_index, data.astype(np.double)
            )

        ### LEGACY ###
        # py_overwrite = False
        # if py_overwrite:
        #     ptr = self.read_record(elem_ptr + ptr_off)[table_index]
        #     if ptr <= 0:
        #         raise ValueError('This element does not have any data associated with the '
        #                          f' solution_type {solution_type}')

        #     with open(self.filename, 'r+b') as f:
        #         # record size
        #         f.seek((ptr_off + elem_ptr + ptr)*4)
        #         sz = int.from_bytes(f.read(4), 'little')

        #         # verify data being written is the same size
        #         if sz != data.size:
        #             raise ValueError(f'Size of the data being overwritten is {sz}, while '
        #                              f'size of the data input is {data.size}.')

        #         # verify no compression (8th byte from start of record)
        #         f.seek(3, 1)
        #         info_byte = f.read(1)
        #         compression = access_bit(info_byte, 3) or \
        #             access_bit(info_byte, 4) or \
        #             access_bit(info_byte, 5)
        #         # Note: file is now at the start of the data

        #         if compression:
        #             raise RuntimeError('Unable to overwrite this compressed record.  '
        #                                'Please disable result file compression and try '
        #                                'again.')

        #         # make input data match data type of the record
        #         prec_flag = access_bit(info_byte, 6)
        #         type_flag = access_bit(info_byte, 7)

        #         if type_flag:
        #             if prec_flag:
        #                 record_dtype = np.int16
        #             else:
        #                 record_dtype = np.int32
        #         else:
        #             if prec_flag:
        #                 record_dtype = np.float32
        #             else:
        #                 record_dtype = np.float64

        #         if data.dtype != record_dtype:
        #             data = data.astype(record_dtype)

        #         # write the record
        #         f.write(data.tobytes())

    def overwrite_element_solution_records(self, element_data, rnum, solution_type):
        """Overwrite element solution record.

        This method replaces solution data for a set of elements at a
        result index for a given solution type.  The number of items
        in ``data`` must match the number of items in the record.

        If you are not sure how many records are in a given record,
        use ``element_solution_data`` to retrieve all the records for
        a given ``solution_type`` and check the number of items in the
        record.

        Note: The record being replaced cannot be a compressed record.
        If the result file uses compression (default sparse
        compression as of 2019R1), you can disable this within MAPDL
        with:
        ``/FCOMP, RST, 0``

        Parameters
        ----------
        element_data : dict
            Dictionary of results that will replace the existing records.

        rnum : int
            Zero based result number.

        solution_type : str
            Element data type to overwrite.

            - EMS: misc. data
            - ENF: nodal forces
            - ENS: nodal stresses
            - ENG: volume and energies
            - EGR: nodal gradients
            - EEL: elastic strains
            - EPL: plastic strains
            - ECR: creep strains
            - ETH: thermal strains
            - EUL: euler angles
            - EFX: nodal fluxes
            - ELF: local forces
            - EMN: misc. non-sum values
            - ECD: element current densities
            - ENL: nodal nonlinear data
            - EHC: calculated heat generations
            - EPT: element temperatures
            - ESF: element surface stresses
            - EDI: diffusion strains
            - ETB: ETABLE items
            - ECT: contact data
            - EXY: integration point locations
            - EBA: back stresses
            - ESV: state variables
            - MNL: material nonlinear record

        Examples
        --------
        Overwrite the elastic strain record for elements 1 and 2 with
        for the first result with random data.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> data = {1: np.random.random(56),
                    2: np.random.random(56)}
        >>> rst.overwrite_element_solution_data(data, 0, 'EEL')
        """
        if not isinstance(element_data, dict):
            raise TypeError("``element_data`` must be a dictionary")

        table_ptr = solution_type.upper()
        if table_ptr not in ELEMENT_INDEX_TABLE_KEYS:
            err_str = "Data type %s is invalid\n" % str(solution_type)
            err_str += "\nAvailable types:\n"
            for key in ELEMENT_INDEX_TABLE_KEYS:
                err_str += "\t%s: %s\n" % (key, element_index_table_info[key])
            raise ValueError(err_str)

        if table_ptr in self.available_results._parsed_bits:
            if table_ptr not in self.available_results:
                raise ValueError(
                    "Result %s is not available in this result file" % table_ptr
                )

        # location of data pointer within each element result table
        table_index = ELEMENT_INDEX_TABLE_KEYS.index(table_ptr)

        rnum = self.parse_step_substep(rnum)
        ele_ind_table, etype, ptr_off = self._element_solution_header(rnum)

        # index of the element in the element table
        for element_id, data in element_data.items():
            if not isinstance(data, np.ndarray):
                data = np.asarray(data)

            elem_ptr = ele_ind_table[self.element_lookup(element_id)]
            if elem_ptr == 0:
                raise ValueError(
                    "This element does not have any data associated "
                    f"with the solution_type {solution_type}."
                )

            # new c stuff here
            if data.dtype == np.float32:
                self._cfile.overwrite_element_data_float(
                    elem_ptr + ptr_off, table_index, data
                )
            elif data.dtype == np.double:
                self._cfile.overwrite_element_data_double(
                    elem_ptr + ptr_off, table_index, data
                )
            else:  # catch-all
                self._cfile.overwrite_element_data_double(
                    elem_ptr + ptr_off, table_index, data.astype(np.double)
                )

    def element_solution_data(self, rnum, datatype, sort=True, **kwargs):
        """Retrieves element solution data.  Similar to ETABLE.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        datatype : str
            Element data type to retrieve.

            - EMS: misc. data
            - ENF: nodal forces
            - ENS: nodal stresses
            - ENG: volume and energies
            - EGR: nodal gradients
            - EEL: elastic strains
            - EPL: plastic strains
            - ECR: creep strains
            - ETH: thermal strains
            - EUL: euler angles
            - EFX: nodal fluxes
            - ELF: local forces
            - EMN: misc. non-sum values
            - ECD: element current densities
            - ENL: nodal nonlinear data
            - EHC: calculated heat generations
            - EPT: element temperatures
            - ESF: element surface stresses
            - EDI: diffusion strains
            - ETB: ETABLE items
            - ECT: contact data
            - EXY: integration point locations
            - EBA: back stresses
            - ESV: state variables
            - MNL: material nonlinear record

        sort : bool
            Sort results by element number.  Default ``True``.

        **kwargs : optional keyword arguments
            Hidden options for distributed result files.

        Returns
        -------
        enum : np.ndarray
            Element numbers.

        element_data : list
            List with one data item for each element.

        enode : list
            Node numbers corresponding to each element.
            results.  One list entry for each element.

        Notes
        -----
        See ANSYS element documentation for available items for each
        element type.  See:

        https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_elem/

        Examples
        --------
        Retrieve "LS" solution results from an PIPE59 element for result set 1

        >>> enum, edata, enode = result.element_solution_data(0, datatype='ENS')
        >>> enum[0]  # first element number
        >>> enode[0]  # nodes belonging to element 1
        >>> edata[0]  # data belonging to element 1
        array([ -4266.19   ,   -376.18857,  -8161.785  , -64706.766  ,
                -4266.19   ,   -376.18857,  -8161.785  , -45754.594  ,
                -4266.19   ,   -376.18857,  -8161.785  ,      0.     ,
                -4266.19   ,   -376.18857,  -8161.785  ,  45754.594  ,
                -4266.19   ,   -376.18857,  -8161.785  ,  64706.766  ,
                -4266.19   ,   -376.18857,  -8161.785  ,  45754.594  ,
                -4266.19   ,   -376.18857,  -8161.785  ,      0.     ,
                -4266.19   ,   -376.18857,  -8161.785  , -45754.594  ,
                -4274.038  ,   -376.62527,  -8171.2603 ,   2202.7085 ,
               -29566.24   ,   -376.62527,  -8171.2603 ,   1557.55   ,
               -40042.613  ,   -376.62527,  -8171.2603 ,      0.     ,
               -29566.24   ,   -376.62527,  -8171.2603 ,  -1557.55   ,
                -4274.038  ,   -376.62527,  -8171.2603 ,  -2202.7085 ,
                21018.164  ,   -376.62527,  -8171.2603 ,  -1557.55   ,
                31494.537  ,   -376.62527,  -8171.2603 ,      0.     ,
                21018.164  ,   -376.62527,  -8171.2603 ,   1557.55   ],
              dtype=float32)

        This data corresponds to the results you would obtain directly
        from MAPDL with ESOL commands:

        >>> ansys.esol(nvar='2', elem=enum[0], node=enode[0][0], item='LS', comp=1)
        >>> ansys.vget(par='SD_LOC1', ir='2', tstrt='1') # store in a variable
        >>> ansys.read_float_parameter('SD_LOC1(1)')
        -4266.19
        """
        table_ptr = datatype.upper()
        if table_ptr not in ELEMENT_INDEX_TABLE_KEYS:
            err_str = "Data type %s is invalid\n" % str(datatype)
            err_str += "\nAvailable types:\n"
            for key in ELEMENT_INDEX_TABLE_KEYS:
                err_str += "\t%s: %s\n" % (key, element_index_table_info[key])
            raise ValueError(err_str)

        if table_ptr in self.available_results._parsed_bits:
            if table_ptr not in self.available_results:
                raise ValueError(
                    "Result %s is not available in this result file" % table_ptr
                )

        # location of data pointer within each element result table
        table_index = ELEMENT_INDEX_TABLE_KEYS.index(table_ptr)

        rnum = self.parse_step_substep(rnum)
        ele_ind_table, etype, ptr_off = self._element_solution_header(rnum)

        # read element data
        if self._cfile is not None:
            element_data = self._cfile.read_element_data(
                ele_ind_table, table_index, ptr_off
            )
        else:
            element_data = []
            for ind in ele_ind_table:
                if ind == 0:
                    element_data.append(None)
                else:
                    # read element table index pointer to data
                    ptr = self.read_record(ind + ptr_off)[table_index]
                    if ptr > 0:
                        record = self.read_record(ptr_off + ind + ptr)
                        element_data.append(record)
                    else:
                        element_data.append(None)

        # just return the distributed result if requested
        if kwargs.get("is_dist_rst", False):
            return element_data

        enum = self._eeqv
        if sort:
            sidx = np.argsort(enum)
            enum = enum[sidx]
            element_data = [element_data[i] for i in sidx]

        enode = []
        nnode = self._nodstr[etype]
        if sort:
            enode = [self._mesh.elem[i][10 : 10 + nnode[i]] for i in sidx]
        else:
            enode = [self._mesh.elem[i][10 : 10 + nnode[i]] for i in range(enum.size)]

        return enum, element_data, enode

    def principal_nodal_stress(self, rnum, nodes=None):
        """Computes the principal component stresses for each node in
        the solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        pstress : numpy.ndarray
            Principal stresses, stress intensity, and equivalent stress.
            [sigma1, sigma2, sigma3, sint, seqv]

        Examples
        --------
        Load the principal nodal stress for the first solution.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, stress = rst.principal_nodal_stress(0)

        Notes
        -----
        ANSYS equivalent of:
        PRNSOL, S, PRIN

        which returns:
        S1, S2, S3 principal stresses, SINT stress intensity, and SEQV
        equivalent stress.

        Internal averaging algorithm averages the component values
        from the elements at a common node and then calculates the
        principal using the averaged value.

        See the MAPDL ``AVPRIN`` command for more details.
        ``ansys-mapdl-reader`` uses the default ``AVPRIN, 0`` option.

        """
        # get component stress
        nodenum, stress = self.nodal_stress(rnum, nodes)
        pstress, isnan = _binary_reader.compute_principal_stress(stress)
        pstress[isnan] = np.nan
        return nodenum, pstress

    def plot_principal_nodal_stress(
        self,
        rnum,
        comp=None,
        show_displacement=False,
        displacement_factor=1.0,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot the principal stress.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : string
            Stress component to plot.  S1, S2, S3 principal stresses, SINT
            stress intensity, and SEQV equivalent stress.

            Stress type must be a string from the following list:
            ``['S1', 'S2', 'S3', 'SINT', 'SEQV']``

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Returns
        -------
        cpos : list
            VTK camera position.

        Examples
        --------
        Plot the equivalent von mises stress.

        >>> rst.plot_principal_nodal_stress(0, comp='SEQV')

        """
        # get the correct component of the principal stress
        comp = comp.upper()
        idx = check_comp(PRINCIPAL_STRESS_TYPES, comp)
        stress = self.principal_nodal_stress(rnum)[1][:, idx]

        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            stress = stress[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            stress = stress[ind]
        else:
            grid = self.grid

        stype_stitle_map = {
            "S1": "Principal Stress 1",
            "S2": "Principal Stress 2",
            "S3": "Principal Stress 3",
            "SINT": "Stress Intensity",
            "SEQV": "von Mises Stress",
        }
        kwargs.setdefault("scalar_bar_args", {"title": stype_stitle_map[comp]})
        return self._plot_point_scalars(
            stress,
            grid=grid,
            rnum=rnum,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            treat_nan_as_zero=treat_nan_as_zero,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            **kwargs,
        )

    def cs_4x4(self, cs_cord, as_vtk_matrix=False):
        """return a 4x4 transformation array for a given coordinate system"""
        # assemble 4 x 4 matrix
        csys = self._c_systems[cs_cord]
        trans = np.hstack(
            (csys["transformation matrix"], csys["origin"].reshape(-1, 1))
        )
        matrix = trans_to_matrix(trans)

        if as_vtk_matrix:
            return matrix

        return pv.array_from_vtkmatrix(matrix)

    def _plot_point_scalars(
        self,
        scalars,
        rnum=None,
        grid=None,
        show_displacement=False,
        displacement_factor=1,
        add_text=True,
        animate=False,
        n_frames=100,
        overlay_wireframe=False,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        movie_filename=None,
        treat_nan_as_zero=True,
        progress_bar=True,
        **kwargs,
    ):
        """Plot point scalars on active mesh.

        Parameters
        ----------
        scalars : np.ndarray
            Node scalars to plot.

        rnum : int, optional
            Cumulative result number.  Used for adding informative
            text.

        grid : pyvista.PolyData or pyvista.UnstructuredGrid, optional
            Uses ``self.grid`` by default.  When specified, uses this
            grid instead.

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        add_text : bool, optional
            Adds information about the result when rnum is given.

        overlay_wireframe : bool, optional
            Overlay a wireframe of the original undeformed mesh.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        progress_bar : bool, optional
            Displays a progress bar when generating a movie while
            ``off_screen=True``.  Default is ``True``.

        kwargs : keyword arguments
            Additional keyword arguments.

        Returns
        -------
        cpos : list
            Camera position.
        """
        if rnum:
            rnum = self.parse_step_substep(rnum)
        loop = kwargs.pop("loop", False)

        if grid is None:
            grid = self.grid

        if treat_nan_as_zero and scalars is not None:
            scalars[np.isnan(scalars)] = 0

        disp = None
        if show_displacement and not animate:
            disp = self.nodal_solution(rnum)[1][:, :3] * displacement_factor
            if node_components:
                _, ind = self._extract_node_components(node_components, sel_type_all)
                disp = disp[ind]
            elif element_components:
                _, ind = self._extract_element_components(element_components)
                disp = disp[ind]

            if disp.shape[1] == 2:  # ignore Z
                new_points = grid.points.copy()
                new_points[:, :-1] += disp
            else:
                new_points = disp + grid.points
            grid = grid.copy()
            grid.points = new_points

        elif animate:
            disp = self.nodal_solution(rnum)[1][:, :3] * displacement_factor

        # extract mesh surface
        mapped_indices = None
        if "vtkOriginalPointIds" in grid.point_data:
            mapped_indices = grid.point_data["vtkOriginalPointIds"]

        # weird bug in some edge cases, have to roll our own indices here
        if "_temp_id" not in grid.point_data:
            grid.point_data["_temp_id"] = np.arange(grid.n_points)
        mesh = grid.extract_surface()
        ind = mesh.point_data["_temp_id"]

        if disp is not None:
            if mapped_indices is not None:
                disp = disp[mapped_indices][ind]
            else:
                disp = disp[ind]

        if scalars is not None:
            if scalars.ndim == 2:
                scalars = scalars[:, ind]
            else:
                scalars = scalars[ind]

            if animate:
                smax = np.abs(scalars).max()
                rng = kwargs.pop("rng", [-smax, smax])
            else:
                rng = kwargs.pop("rng", [scalars.min(), scalars.max()])
        else:
            rng = kwargs.pop("rng", None)

        return_cpos = kwargs.pop("return_cpos", False)
        cmap = kwargs.pop("cmap", "jet")
        window_size = kwargs.pop("window_size", None)
        full_screen = kwargs.pop("full_screen", None)
        notebook = kwargs.pop("notebook", None)
        off_screen = kwargs.pop("off_screen", None)
        cpos = kwargs.pop("cpos", None)
        screenshot = kwargs.pop("screenshot", None)
        interactive = kwargs.pop("interactive", True)
        text_color = kwargs.pop("text_color", None)

        kwargs.setdefault("smooth_shading", None)
        kwargs.setdefault("color", "w")
        kwargs.setdefault("interpolate_before_map", True)

        plotter = pv.Plotter(off_screen=off_screen, notebook=notebook)

        # set axes
        if kwargs.pop("show_axes", True):
            plotter.add_axes()

        # set background
        theme = DefaultTheme()
        theme.background = kwargs.pop("background", None)

        # remove extra keyword args
        kwargs.pop("node_components", None)
        kwargs.pop("sel_type_all", None)

        if overlay_wireframe:
            plotter.add_mesh(self.grid, style="wireframe", color="w", opacity=0.5)

        copied_mesh = mesh.copy()
        plotter.add_mesh(copied_mesh, scalars=scalars, rng=rng, cmap=cmap, **kwargs)

        # set scalar bar text colors
        if text_color:
            theme = DefaultTheme()
            theme.color = text_color
            pv.global_theme.load_theme(theme)

        # NAN/missing data are white
        # plotter.renderers[0].SetUseDepthPeeling(1)  # <-- for transparency issues
        theme.nan_color = [1, 1, 1, 1]

        if cpos:
            plotter.camera_position = cpos

        if movie_filename:
            movie_filename = str(movie_filename)
            if movie_filename.strip()[-3:] == "gif":
                plotter.open_gif(movie_filename)
            else:
                plotter.open_movie(movie_filename)

        # add table
        if add_text and rnum is not None:
            result_text = self.text_result_table(rnum)
            plotter.add_text(result_text, font_size=20, color=text_color)

        # camera position added in 0.32.0
        show_kwargs = {}
        if pv._version.version_info >= (0, 32):
            show_kwargs["return_cpos"] = return_cpos

        if animate:
            if off_screen:  # otherwise this never exits
                loop = False
            pbar = None
            if progress_bar and off_screen and movie_filename is not None:
                pbar = tqdm(total=n_frames, desc="Rendering animation")

            orig_pts = copied_mesh.points.copy()

            plotter.show(
                interactive=False,
                auto_close=False,
                window_size=window_size,
                full_screen=full_screen,
                interactive_update=not off_screen,
                **show_kwargs,
            )

            self._animating = True

            def q_callback():
                """exit when user wants to leave"""
                self._animating = False

            def exit_callback(plotter, RenderWindowInteractor, event):
                """exit when user wants to leave"""
                self._animating = False
                plotter.close()

            plotter.add_key_event("q", q_callback)
            if os.name == "nt":
                # Adding closing window callback
                plotter.iren.add_observer(
                    vtk.vtkCommand.ExitEvent,
                    lambda render, event: exit_callback(plotter, render, event),
                )

            first_loop = True
            cached_normals = [None for _ in range(n_frames)]
            while self._animating:
                for j, angle in enumerate(np.linspace(0, np.pi * 2, n_frames + 1)[:-1]):
                    mag_adj = np.sin(angle)
                    if scalars is not None:
                        copied_mesh.active_scalars[:] = scalars * mag_adj
                    copied_mesh.points[:] = orig_pts + disp * mag_adj

                    # normals have to be updated on the fly
                    if kwargs["smooth_shading"]:
                        if cached_normals[j] is None:
                            copied_mesh.compute_normals(
                                cell_normals=False, inplace=True
                            )
                            cached_normals[j] = copied_mesh.point_data["Normals"]
                        else:
                            copied_mesh.point_data["Normals"][:] = cached_normals[j]

                    if add_text:
                        phase = angle * 180 / np.pi
                        plotter.add_text(f"{result_text} \nPhase {phase} Degrees")

                    # at max supported framerate
                    plotter.update(1, force_redraw=True)
                    if not self._animating:
                        break

                    if movie_filename and first_loop:
                        if pbar is not None:
                            pbar.update()
                        plotter.write_frame()

                first_loop = False
                if not loop:
                    break

            plotter.close()
            cpos = plotter.camera_position

        elif screenshot is True:
            cpos, img = plotter.show(
                interactive=interactive,
                window_size=window_size,
                full_screen=full_screen,
                screenshot=screenshot,
                **show_kwargs,
            )

        else:
            cpos = plotter.show(
                interactive=interactive,
                window_size=window_size,
                full_screen=full_screen,
                screenshot=screenshot,
                **show_kwargs,
            )

        if screenshot is True:
            return cpos, img

        return cpos

    def _animate_point_scalars(
        self,
        scalars,
        grid=None,
        text=None,
        movie_filename=None,
        loop=True,
        fps=20,
        **kwargs,
    ):
        """Animate a set of point scalars on active mesh.

        Parameters
        ----------
        scalars : List
            List of node scalars to plot.  One set of values for each frame.

        grid : pyvista.PolyData or pyvista.UnstructuredGrid, optional
            Uses self.grid by default.  When specified, uses this grid
            instead.

        text : List, optional
            One string for each value in ``scalars``.

        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Returns
        -------
        cpos : list
            Camera position.
        """
        # if rnum:
        #     rnum = self.parse_step_substep(rnum)

        if grid is None:
            grid = self.grid

        # disp = None
        # if show_displacement and not animate:
        #     disp = self.nodal_solution(rnum)[1][:, :3]*displacement_factor
        #     if node_components:
        #         _, ind = self._extract_node_components(node_components, sel_type_all)
        #         disp = disp[ind]
        #     elif element_components:
        #         _, ind = self._extract_element_components(element_components)
        #         disp = disp[ind]
        #     new_points = disp + grid.points
        #     grid = grid.copy()
        #     grid.points = new_points
        # elif animate:
        #     disp = self.nodal_solution(rnum)[1][:, :3]*displacement_factor

        # extract mesh surface
        # mapped_indices = None
        # if 'vtkOriginalPointIds' in grid.point_data:
        #     mapped_indices = grid.point_data['vtkOriginalPointIds']

        # weird bug in some edge cases, have to roll our own indices here
        if "_temp_id" not in grid.point_data:
            grid["_temp_id"] = np.arange(grid.n_points)
        mesh = grid.extract_surface()
        ind = mesh.point_data["_temp_id"]

        # if disp is not None:
        #     if mapped_indices is not None:
        #         disp = disp[mapped_indices][ind]
        #     else:
        #         disp = disp[ind]

        # remap scalars
        for i in range(len(scalars)):
            scalars[i] = scalars[i][ind]

        rng = kwargs.pop("rng", [np.min(scalars), np.max(scalars)])
        cmap = kwargs.pop("cmap", "jet")
        window_size = kwargs.pop("window_size", [1024, 768])
        full_screen = kwargs.pop("full_screen", False)
        notebook = kwargs.pop("notebook", False)
        off_screen = kwargs.pop("off_screen", None)
        cpos = kwargs.pop("cpos", None)
        # screenshot = kwargs.pop('screenshot', None)
        # interactive = kwargs.pop('interactive', True)
        text_color = kwargs.pop("text_color", None)

        kwargs.setdefault("smooth_shading", True)
        kwargs.setdefault("color", "w")
        kwargs.setdefault("interpolate_before_map", True)

        plotter = pv.Plotter(off_screen=off_screen, notebook=notebook)

        # set axes
        if kwargs.pop("show_axes", True):
            plotter.add_axes()

        # set background
        theme = DefaultTheme()
        theme.background = kwargs.pop("background", None)

        # remove extra keyword args
        kwargs.pop("node_components", None)
        kwargs.pop("sel_type_all", None)

        # if overlay_wireframe:
        #     plotter.add_mesh(self.grid, style='wireframe', color='w',
        #                      opacity=0.5)

        copied_mesh = mesh.copy()
        plotter.add_mesh(copied_mesh, scalars=scalars[0], rng=rng, cmap=cmap, **kwargs)

        # set scalar bar text colors
        if text_color:
            theme.color = text_color
            pv.global_theme.load_theme(theme)

        # NAN/missing data are white
        # plotter.renderers[0].SetUseDepthPeeling(1)  # <-- for transparency issues
        theme.nan_color = [1, 1, 1, 1]

        if cpos:
            plotter.camera_position = cpos

        if movie_filename:
            movie_filename = str(movie_filename)
            if movie_filename.strip()[-3:] == "gif":
                plotter.open_gif(movie_filename)
            else:
                plotter.open_movie(movie_filename, framerate=fps)

        # add table
        if text is not None:
            if len(text) != len(scalars):
                raise ValueError(
                    "Length of ``text`` must be the same as " "``scalars``"
                )
            plotter.add_text(text[0], font_size=20, color=text_color)

        # orig_pts = copied_mesh.points.copy()
        plotter.show(
            interactive=False,
            auto_close=False,
            window_size=window_size,
            full_screen=full_screen,
            interactive_update=not off_screen,
        )

        self._animating = True

        def q_callback():
            """exit when user wants to leave"""
            self._animating = False

        plotter.add_key_event("q", q_callback)

        if off_screen:
            twait = 0
        else:
            twait = 1 / fps

        first_loop = True
        while self._animating:
            for i, data in enumerate(scalars):
                tstart_render = time.time()
                copied_mesh.active_scalars[:] = data

                if text is not None:
                    plotter.add_text(text[i])

                # at max supported framerate
                plotter.update(1, force_redraw=True)
                if not self._animating:
                    break

                if movie_filename and first_loop:
                    plotter.write_frame()

                # wait to maintain desired frame rate
                telap = time.time() - tstart_render
                if telap < twait:
                    time.sleep(twait - telap)

            first_loop = False
            if not loop:
                break

        plotter.close()
        return plotter.camera_position

    def text_result_table(self, rnum):
        """Returns a text result table for plotting"""
        ls_table = self._resultheader["ls_table"]
        timevalue = self.time_values[rnum]
        text = "Cumulative Index: {:3d}\n".format(ls_table[rnum, 2])
        if self._resultheader["nSector"] > 1:
            hindex = self._resultheader["hindex"][rnum]
            text += "Harmonic Index    {:3d}\n".format(hindex)
        text += "Loadstep:         {:3d}\n".format(ls_table[rnum, 0])
        text += "Substep:          {:3d}\n".format(ls_table[rnum, 1])
        text += "Time Value:     {:10.4f}".format(timevalue)

        return text

    def plot_nodal_stress(
        self,
        rnum,
        comp=None,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plots the stresses at each node in the solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Stress component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Returns
        -------
        cpos : list
            3 x 3 vtk camera position.

        Examples
        --------
        Plot the X component nodal stress while showing displacement.

        >>> rst.plot_nodal_stress(0, comp='x', show_displacement=True)
        """
        kwargs.setdefault(
            "scalar_bar_args", {"title": f"{comp} Component Nodal Stress"}
        )
        return self._plot_nodal_result(
            rnum,
            "ENS",
            comp,
            STRESS_TYPES,
            show_displacement,
            displacement_factor,
            node_components,
            element_components,
            sel_type_all,
            treat_nan_as_zero=treat_nan_as_zero,
            **kwargs,
        )

    def save_as_vtk(
        self, filename, rsets=None, result_types=["ENS"], progress_bar=True
    ):
        """Writes results to a vtk readable file.

        Nodal results will always be written.

        The file extension will select the type of writer to use.
        ``'.vtk'`` will use the legacy writer, while ``'.vtu'`` will
        select the VTK XML writer.

        Parameters
        ----------
        filename : str, pathlib.Path
            Filename of grid to be written.  The file extension will
            select the type of writer to use.  ``'.vtk'`` will use the
            legacy writer, while ``'.vtu'`` will select the VTK XML
            writer.

        rsets : collections.Iterable
            List of result sets to write.  For example ``range(3)`` or
            [0].

        result_types : list
            Result type to write.  For example ``['ENF', 'ENS']``
            List of some or all of the following:

            - EMS: misc. data
            - ENF: nodal forces
            - ENS: nodal stresses
            - ENG: volume and energies
            - EGR: nodal gradients
            - EEL: elastic strains
            - EPL: plastic strains
            - ECR: creep strains
            - ETH: thermal strains
            - EUL: euler angles
            - EFX: nodal fluxes
            - ELF: local forces
            - EMN: misc. non-sum values
            - ECD: element current densities
            - ENL: nodal nonlinear data
            - EHC: calculated heat generations
            - EPT: element temperatures
            - ESF: element surface stresses
            - EDI: diffusion strains
            - ETB: ETABLE items
            - ECT: contact data
            - EXY: integration point locations
            - EBA: back stresses
            - ESV: state variables
            - MNL: material nonlinear record

        progress_bar : bool, optional
            Display a progress bar using ``tqdm``.

        Notes
        -----
        Binary files write much faster than ASCII, but binary files
        written on one system may not be readable on other systems.
        Binary can only be selected for the legacy writer.

        Examples
        --------
        Write nodal results as a binary vtk file.

        >>> rst.save_as_vtk('results.vtk')

        Write using the xml writer

        >>> rst.save_as_vtk('results.vtu')

        Write only nodal and elastic strain for the first result

        >>> rst.save_as_vtk('results.vtk', [0], ['EEL', 'EPL'])

        Write only nodal results (i.e. displacements) for the first result.

        >>> rst.save_as_vtk('results.vtk', [0], [])

        """
        # Copy grid as to not write results to original object
        grid = self.quadgrid.copy()

        if rsets is None:
            rsets = range(self.nsets)
        elif isinstance(rsets, int):
            rsets = [rsets]
        elif not isinstance(rsets, Iterable):
            raise TypeError("rsets must be an iterable like [0, 1, 2] or range(3)")

        if result_types is None:
            result_types = ELEMENT_INDEX_TABLE_KEYS
        elif not isinstance(result_types, list):
            raise TypeError("result_types must be a list of solution types")
        else:
            for item in result_types:
                if item not in ELEMENT_INDEX_TABLE_KEYS:
                    raise ValueError(f'Invalid result type "{item}"')

        pbar = None
        if progress_bar:
            pbar = tqdm(total=len(rsets), desc="Saving to file")

        for i in rsets:
            # Nodal results
            _, val = self.nodal_solution(i)
            grid.point_data["Nodal Solution {:d}".format(i)] = val

            # Nodal results
            for rtype in self.available_results:
                if rtype in result_types:
                    _, values = self._nodal_result(i, rtype)
                    desc = element_index_table_info[rtype]
                    grid.point_data["{:s} {:d}".format(desc, i)] = values

            if pbar is not None:
                pbar.update(1)

        grid.save(str(filename))
        if pbar is not None:
            pbar.close()

    def write_tables(self, filename: Union[str, pathlib.Path]):
        """Write binary tables to ASCII.  Assumes int32

        Parameters
        ----------
        filename : str, pathlib.Path
            Filename to write the tables to.

        Examples
        --------
        >>> rst.write_tables('tables.txt')
        """
        rawresult = open(self.pathlib_filename, "rb")
        with open(filename, "w") as f:
            while True:
                try:
                    table = read_table(rawresult)
                    f.write("*** %d ***\n" % len(table))
                    for item in table:
                        f.write(str(item) + "\n")
                    f.write("\n\n")
                except:
                    break
        rawresult.close()

    def parse_step_substep(self, user_input):
        """Converts (step, substep) to a cumulative index"""
        if is_int(user_input):
            # check if result exists
            if user_input > self.nsets - 1:
                raise ValueError(
                    "There are only %d result(s) in the result file." % self.nsets
                )
            return user_input

        elif isinstance(user_input, (list, tuple)):
            if len(user_input) != 2:
                raise ValueError(
                    "Input must contain (step, loadstep) using  "
                    + "1 based indexing (e.g. (1, 1))."
                )
            ls_table = self._resultheader["ls_table"]
            step, substep = user_input
            mask = np.logical_and(ls_table[:, 0] == step, ls_table[:, 1] == substep)

            if not np.any(mask):
                raise ValueError(
                    "Load step table does not contain "
                    + "step %d and substep %d" % tuple(user_input)
                )

            index = mask.nonzero()[0]
            assert index.size == 1, "Multiple cumulative index matches"
            return index[0]
        else:
            raise TypeError("Input must be either an int or a list")

    def __repr__(self):
        if self._is_distributed:
            rst_info = ["PyMAPDL Reader Distributed Result"]
        else:
            rst_info = ["PyMAPDL Result"]
        keys = ["title", "subtitle", "units"]
        for key in keys:
            value = self._resultheader[key]
            if value:
                rst_info.append("{:<12s}: {:s}".format(key.capitalize(), value))

        value = self._resultheader["verstring"]
        rst_info.append("{:<12s}: {:s}".format("Version", value))

        value = str(self._resultheader["nSector"] > 1)
        rst_info.append("{:<12s}: {:s}".format("Cyclic", value))

        value = self._resultheader["nsets"]
        rst_info.append("{:<12s}: {:d}".format("Result Sets", value))

        value = self._resultheader["nnod"]
        rst_info.append("{:<12s}: {:d}".format("Nodes", value))

        value = self._resultheader["nelm"]
        rst_info.append("{:<12s}: {:d}".format("Elements", value))

        rst_info.append("\n")
        rst_info.append(str(self.available_results))
        return "\n".join(rst_info)

    def _nodal_result(self, rnum, result_type, nnum_of_interest=None):
        """Generic load nodal result

        Parameters
        ----------
        rnum : int
            Result number.

        result_type : int
            EMS: misc. data
            ENF: nodal forces
            ENS: nodal stresses
            ENG: volume and energies
            EGR: nodal gradients
            EEL: elastic strains
            EPL: plastic strains
            ECR: creep strains
            ETH: thermal strains
            EUL: euler angles
            EFX: nodal fluxes
            ELF: local forces
            EMN: misc. non-sum values
            ECD: element current densities
            ENL: nodal nonlinear data
            EHC: calculated heat
            EPT: element temperatures
            ESF: element surface stresses
            EDI: diffusion strains
            ETB: ETABLE items (post1 only)
            ECT: contact data
            EXY: integration point locations
            EBA: back stresses
            ESV: state variables
            MNL: material nonlinear record

        nnum_of_interest : list, np.ndarray, optional
            Array of nodes numbers to extract results for.  All
            adjacent elements are still considered, so averaging will
            still account for unselected nodes.

        Returns
        -------
        np.ndarray
            Ansys node numbers

        np.ndarray
            Array of result data
        """
        if not self.available_results[result_type]:
            rtype_str = element_index_table_info[result_type]
            raise ValueError(
                f"Result {rtype_str} ({result_type}) is not available in "
                "this result file."
            )

        # element header
        rnum = self.parse_step_substep(rnum)
        ele_ind_table, ans_etype, ptr_off = self._element_solution_header(rnum)

        result_type = result_type.upper()
        nitem = self._result_nitem(rnum, result_type)
        result_index = ELEMENT_INDEX_TABLE_KEYS.index(result_type)

        if nnum_of_interest is not None:
            nnum_sel = np.unique(nnum_of_interest)
            mask = np.in1d(self.mesh.nnum, nnum_sel)

            # extract any elements containing these values
            # note that we need this for accurate averaging of results
            grid = self.grid.extract_points(
                np.nonzero(mask)[0], adjacent_cells=True, include_cells=True
            )

            # mask relative to global eeqv array
            ele_mask = np.in1d(self._eeqv, grid["ansys_elem_num"], assume_unique=True)

            # verify that all nodes of interest actually occur within
            node_mask = np.in1d(nnum_sel, grid["ansys_node_num"], assume_unique=True)
            if not node_mask.all():
                warnings.warn(
                    "Not all nodes IDs in the ``nnum_of_interest`` "
                    "are within the result"
                )
            etype = self._mesh.etype[ele_mask]
            ele_ind_table = ele_ind_table[ele_mask]

        else:
            grid = self.grid
            etype = self._mesh.etype

        # Ignore contributions from SURF154 for temperature results
        ignore_154 = result_type in ["EPT"]

        # NOTE: For nodal temperatures:
        # For solid elements and SHELL41, the record contains nodal
        # temperatures at each node and the number of items in this
        # record is nodfor.
        # nodfor
        nodstr = self._nodstr.copy()  # copy here since we may change it inplace
        if result_type == "EPT":
            # replace values in nodstr
            is_solid = elements.SOLID_ELEMENTS[self.mesh.ekey[:, 1]]
            is_solid += self.mesh.ekey[:, 1] == 41  # must include SHELL41
            nodstr[1:][is_solid] = self._nodfor[1:][is_solid]

        # call cython to extract and assemble the data
        cells, offset = vtk_cell_info(grid)
        data, ncount = _binary_reader.read_nodal_values(
            self.filename,
            grid.celltypes,
            ele_ind_table,
            offset,
            cells,
            nitem,
            grid.n_points,
            nodstr,
            ans_etype,
            etype,
            result_index,
            ptr_off,
            ignore_154,
        )

        # sanity check
        if not np.any(ncount):
            raise ValueError(
                "Result file contains no %s records for result %d"
                % (element_index_table_info[result_type.upper()], rnum)
            )

        if result_type == "ENS" and nitem != 6:
            data = data[:, :6]

        if nnum_of_interest is not None:
            if grid.n_points != nnum_sel.size:
                mask = np.in1d(grid["ansys_node_num"], nnum_sel)
                nnum = grid["ansys_node_num"][mask]
                ncount = ncount[mask]
                data = data[mask]
            else:
                nnum = self.grid.point_data["ansys_node_num"]
        else:
            nnum = self.grid.point_data["ansys_node_num"]

        # average across nodes
        ncount = ncount.reshape(-1, 1)
        return nnum, data / ncount

    @property
    def _nodstr(self):
        """Number of nodes per element having stresses or strains.

        For higher-order elements, nodstr equals to the number of
        corner nodes (e.g., for 20-noded SOLID186, nodstr = 8). NOTE-
        the /config NST1 option may suppress all but one node output
        of nodal ENS EEL EPL ECR ETH and ENL for any given element.
        """
        return self._element_table["nodstr"]

    @property
    def _nodfor(self):
        """Number of nodes per element having nodal forces"""
        return self._element_table["nodfor"]

    def _result_nitem(self, rnum, result_type):
        """Return the number of items for a given result type"""
        if self._resultheader["rstsprs"] == 0 and result_type == "ENS":
            nitem = 11
        elif result_type == "ENF":
            # number of items is nodfor*NDOF*M where:
            # NDOF is the number of DOFs/node for this element, nodfor
            # is the number of nodes per element having nodal forces
            # (defined in element type description record), and M may
            # be 1, 2, or 3.  For a static analysis, M=1 only.  For a
            # transient analysis, M can be 1, 2, or 3.

            solution_header = self._result_solution_header(rnum)
            numdof = solution_header["numdof"]
            nitem = numdof * 1  # m = 1  for static, 1, 2 or 3 for transient

        else:
            nitem = ELEMENT_RESULT_NCOMP.get(result_type, 1)
        return nitem

    def nodal_reaction_forces(self, rnum):
        """Nodal reaction forces.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        Returns
        -------
        rforces : np.ndarray
            Nodal reaction forces for each degree of freedom.

        nnum : np.ndarray
            Node numbers corresponding to the reaction forces.  Node
            numbers may be repeated if there is more than one degree
            of freedom for each node.

        dof : np.ndarray
            Degree of freedom corresponding to each node using the
            MAPDL degree of freedom reference table.  See
            ``rst.result_dof`` for the corresponding degrees of
            freedom for a given solution.

        Examples
        --------
        Get the nodal reaction forces for the first result and print
        the reaction forces of a single node.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> rforces, nnum, dof = rst.nodal_reaction_forces(0)
        >>> dof_ref = rst.result_dof(0)
        >>> rforces[:3], nnum[:3], dof[:3], dof_ref
        (array([  24102.21376091, -109357.01854005,   22899.5303263 ]),
         array([4142, 4142, 4142]),
         array([1, 2, 3], dtype=int32),
         ['UX', 'UY', 'UZ'])

        """
        # PROGRAMMERS NOTE:
        # RF     LONG      1       nrf     Reaction force DOFs.  This
        #                                  index is calculated as
        #                                  (N-1)*numdof+DOF, where N
        #                                  is the position number of
        #                                  the node in the nodal
        #                                  equivalence table, and DOF
        #                                  is the DOF reference
        #                                  number.
        #
        #  ---     dp       1       nrf    Reaction forces.  The force
        #                                  values are ordered
        #                                  according to the DOF order
        #                                  shown above in the DOF
        #                                  number reference table.

        # pointer to result reaction forces
        rnum = self.parse_step_substep(rnum)
        rpointers = self._resultheader["rpointers"]
        solution_header = self._solution_header(rnum)
        ptr = rpointers[rnum] + solution_header["ptrRF"]
        n_rf = solution_header["nrf"]
        numdof = solution_header["numdof"]

        # Reaction force DOFs
        # table is always ANSYS LONG (INT64)
        table, bufsz = self.read_record(ptr, True)
        table = table.view(np.int64)

        # reaction forces
        rforces, bufsz = self.read_record(ptr + bufsz, return_bufsize=True)
        rforces = rforces[:n_rf]

        solution_header = self._result_solution_header(rnum)

        # following ``(N - 1)*numdof + DOF``
        shifted_table = (table - 1) / numdof
        index = np.array(shifted_table, np.int64)
        dof = np.round(1 + numdof * (shifted_table - index)).astype(np.int32)

        # index appears to be sorted
        idx = (table - dof) // numdof
        nnum = self._resultheader["neqv"][idx]
        return rforces, nnum, dof

    def plot_element_result(
        self, rnum, result_type, item_index, in_element_coord_sys=False, **kwargs
    ):
        """Plot an element result.

        Parameters
        ----------
        rnum : int
            Result number.

        result_type : str
            Element data type to retrieve.

            - EMS: misc. data
            - ENF: nodal forces
            - ENS: nodal stresses
            - ENG: volume and energies
            - EGR: nodal gradients
            - EEL: elastic strains
            - EPL: plastic strains
            - ECR: creep strains
            - ETH: thermal strains
            - EUL: euler angles
            - EFX: nodal fluxes
            - ELF: local forces
            - EMN: misc. non-sum values
            - ECD: element current densities
            - ENL: nodal nonlinear data
            - EHC: calculated heat generations
            - EPT: element temperatures
            - ESF: element surface stresses
            - EDI: diffusion strains
            - ETB: ETABLE items
            - ECT: contact data
            - EXY: integration point locations
            - EBA: back stresses
            - ESV: state variables
            - MNL: material nonlinear record

        item_index : int
            Index of the data item for each node within the element.

        in_element_coord_sys : bool, optional
            Returns the results in the element coordinate system.
            Default False and will return the results in the global
            coordinate system.

        Returns
        -------
        nnum : np.ndarray
            ANSYS node numbers

        result : np.ndarray
            Array of result data
        """
        # check result exists
        result_type = result_type.upper()
        if not self.available_results[result_type]:
            raise ValueError(
                f"Result {result_type} is not available in this result file"
            )

        if result_type not in ELEMENT_INDEX_TABLE_KEYS:
            raise ValueError(f'Result "{result_type}" is not an element result')

        bsurf = self._extract_surface_element_result(
            rnum, result_type, item_index, in_element_coord_sys
        )

        desc = self.available_results.description[result_type].capitalize()
        kwargs.setdefault("scalar_bar_args", {"title": desc})
        return bsurf.plot(scalars="_scalars", **kwargs)

    def _extract_surface_element_result(
        self, rnum, result_type, item_index, in_element_coord_sys
    ):
        """Return the surface of the grid with the active scalars
        being the element result"""
        # element header
        rnum = self.parse_step_substep(rnum)
        ele_ind_table, etype, ptr_off = self._element_solution_header(rnum)

        # the number of items per node
        nitem = self._result_nitem(rnum, result_type)
        if item_index > nitem - 1:
            raise ValueError(
                "Item index greater than the number of items in "
                "this result type %s" % result_type
            )

        # extract the surface and separate the surface faces
        # TODO: add element/node components
        surf = self.grid.extract_surface()
        bsurf = break_apart_surface(surf, force_linear=True)
        nnum_surf = surf.point_data["ansys_node_num"][bsurf["orig_ind"]]
        faces = bsurf.faces
        if faces.dtype != np.int64:
            faces = faces.astype(np.int64)

        elem_ind = surf.cell_data["vtkOriginalCellIds"]
        # index within the element table pointing to the data of interest
        result_index = ELEMENT_INDEX_TABLE_KEYS.index(result_type)

        data = populate_surface_element_result(
            self.filename,
            ele_ind_table,
            self._nodstr,
            etype,
            nitem,
            ptr_off,  # start of result data
            result_index,
            bsurf.n_points,
            faces,
            bsurf.n_faces,
            nnum_surf,
            elem_ind,
            self._mesh._elem,
            self._mesh._elem_off,
            item_index,
            as_global=not in_element_coord_sys,
        )
        bsurf["_scalars"] = data
        return bsurf

    def _result_solution_header(self, rnum):
        """Return the solution header for a given cumulative result index."""
        ptr = self._resultheader["rpointers"][rnum]
        return parse_header(self.read_record(ptr), solution_data_header_keys)

    def _result_solution_header_ext(self, rnum):
        """Return the extended solution header for a given cumulative result index.

        This header remains unparsed:

        Header extension
        positions   1-32  - current DOF for this
                            result set
        positions  33-64  - current DOF labels for
                            this result set
        positions  65-84  - The third title, in
                            integer form
        positions  85-104 - The fourth title, in
                            integer form
        positions 105-124 - The fifth title, in
                            integer form
        position 125 - unused
        position 126 - unused
        position 127 - numvdof, number of velocity
                       items per node (ANSYS
                       transient)
        position 128 - numadof, number of
                       acceleration items per
                       node (ANSYS transient)
        position 131-133 - position of velocity
                           in DOF record
                           (ANSYS transient)
        position 134-136 - position of acceleration
                           in DOF record
                           (ANSYS transient)
        position 137-142 - velocity and
                           acceleration labels
                           (ANSYS transient)
        position 143 - number of stress items
                       (6 or 11); a -11 indicates
                       to use principals directly
                       and not recompute (for PSD)
        position 144-146 - position of rotational
                           velocity in DOF record
                           (ANSYS transient)
        position 147-149 - position of rotational
                           accel. in DOF record
                           (ANSYS transient)
        position 150-155 - rotational velocity and
                           acceleration labels
                           (ANSYS transient)
        position 160 - pointer to CINT results
        position 161 - size of CINT record
        position 162 - Pointer to  XFEM/SMART-crack surface information
        position 163 - size of crk surface records
        position 164-200 - unused

        """
        # adding 203 here as the header is 200 words + 3 for padding and the
        # dp solution header is 203 + 3
        ptr = self._resultheader["rpointers"][rnum] + (203 + 203)
        return self.read_record(ptr)

    def nodal_stress(self, rnum, nodes=None):
        """Retrieves the component stresses for each node in the
        solution.

        The order of the results corresponds to the sorted node
        numbering.

        Computes the nodal stress by averaging the stress for each
        element at each node.  Due to the discontinuities across
        elements, stresses will vary based on the element they are
        evaluated from.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : numpy.ndarray
            Node numbers of the result.

        stress : numpy.ndarray
            Stresses at ``X, Y, Z, XY, YZ, XZ`` averaged at each corner
            node.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, stress = rst.nodal_stress(0)

        Return the nodal stress just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, stress = rst.nodal_stress(0, nodes='MY_COMPONENT')

        Return the nodal stress just for the nodes from 20 through 50.

        >>> nnum, stress = rst.nodal_solution(0, nodes=range(20, 51))

        Notes
        -----
        Nodes without a stress value will be NAN.
        Equivalent ANSYS command: PRNSOL, S
        """
        return self._nodal_result(rnum, "ENS", nnum_of_interest=nodes)

    def cylindrical_nodal_stress(self, rnum, nodes=None):
        """Retrieves the stresses for each node in the solution in the
        cylindrical coordinate system as the following values:

        ``R``, ``THETA``, ``Z``, ``RTHETA``, ``THETAZ``, and ``RZ``

        The order of the results corresponds to the sorted node
        numbering.

        Computes the nodal stress by averaging the stress for each
        element at each node.  Due to the discontinuities across
        elements, stresses will vary based on the element they are
        evaluated from.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : numpy.ndarray
            Node numbers of the result.

        stress : numpy.ndarray
            Stresses at ``R, THETA, Z, RTHETA, THETAZ, RZ`` averaged
            at each corner node where ``R`` is radial.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, stress = rst.cylindrical_nodal_stress(0)

        Return the cylindrical nodal stress just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, stress = rst.cylindrical_nodal_stress(0, nodes='MY_COMPONENT')

        Return the nodal stress just for the nodes from 20 through 50.

        >>> nnum, stress = rst.cylindrical_nodal_stress(0, nodes=range(20, 51))

        Notes
        -----
        Nodes without a stress value will be NAN.
        Equivalent ANSYS commands:
        RSYS, 1
        PRNSOL, S
        """
        nnum, stress = self._nodal_result(rnum, "ENS", nnum_of_interest=nodes)

        # angles relative to the XZ plane
        if nnum.size != self._mesh.nodes.shape[0]:
            mask = np.in1d(nnum, self._mesh.nnum)
            angle = np.arctan2(self._mesh.nodes[mask, 1], self._mesh.nodes[mask, 0])
        else:
            angle = np.arctan2(self._mesh.nodes[:, 1], self._mesh.nodes[:, 0])

        _binary_reader.euler_cart_to_cyl(stress, angle)  # mod stress inplace
        return nnum, stress

    def nodal_temperature(self, rnum, nodes=None, **kwargs):
        """Retrieves the temperature for each node in the
        solution.

        The order of the results corresponds to the sorted node
        numbering.

        Equivalent MAPDL command: PRNSOL, TEMP

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : numpy.ndarray
            Node numbers of the result.

        temperature : numpy.ndarray
            Temperature at each node.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, temp = rst.nodal_temperature(0)

        Return the temperature just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, temp = rst.nodal_stress(0, nodes='MY_COMPONENT')

        Return the temperature just for the nodes from 20 through 50.

        >>> nnum, temp = rst.nodal_solution(0, nodes=range(20, 51))

        """
        if self._is_thermal:
            nnum, temp = self.nodal_solution(rnum, nodes=nodes)
        else:
            nnum, temp = self._nodal_result(rnum, "EPT", nnum_of_interest=nodes)
        return nnum, temp.ravel()

    def plot_cylindrical_nodal_stress(
        self,
        rnum,
        comp=None,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal_stress in the cylindrical coordinate system.

        Parameters
        ----------
        rnum : int
            Result number

        comp : str, optional
            Stress component to display.  Available options:
            - ``"R"``
            - ``"THETA"``
            - ``"Z"``
            - ``"RTHETA"``
            - ``"THETAZ"``
            - ``"RZ"``

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``

        Examples
        --------
        Plot nodal stress in the radial direction.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> result = pymapdl_reader.read_binary('file.rst')
        >>> result.plot_cylindrical_nodal_stress(0, 'R')
        """
        available_comps = ["R", "THETA", "Z", "RTHETA", "THETAZ", "RZ"]
        idx = check_comp(available_comps, comp)
        _, scalars = self.cylindrical_nodal_stress(rnum)
        scalars = scalars[:, idx]
        grid = self.grid

        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            scalars = scalars[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            scalars = scalars[ind]

        kwargs.setdefault(
            "scalar_bar_args", {"title": f"{comp} Cylindrical\nNodal Stress"}
        )
        return self._plot_point_scalars(
            scalars,
            grid=grid,
            rnum=rnum,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            treat_nan_as_zero=treat_nan_as_zero,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            **kwargs,
        )

    def plot_nodal_temperature(
        self,
        rnum,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal temperature

        Parameters
        ----------
        rnum : int
            Result number

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``

        Examples
        --------
        Plot temperature of a result.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> result = pymapdl_reader.read_binary('file.rst')
        >>> result.plot_nodal_temperature(0)

        Plot while showing edges and disabling lighting

        >>> result.plot_nodal_temperature(0, show_edges=True, lighting=False)
        """
        _, scalars = self.nodal_temperature(rnum)

        grid = self.grid
        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            scalars = scalars[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            scalars = scalars[ind]

        return self._plot_point_scalars(
            scalars,
            grid=grid,
            rnum=rnum,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            scalar_bar_args={"title": "Nodal Temperature"},
            treat_nan_as_zero=treat_nan_as_zero,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            **kwargs,
        )

    def nodal_thermal_strain(self, rnum, nodes=None):
        """Nodal component thermal strain.

        This record contains strains in the order X, Y, Z, XY, YZ, XZ,
        EQV, and eswell (element swelling strain).  Thermal strains
        are always values at the integration points moved to the
        nodes.

        Equivalent MAPDL command: PRNSOL, EPTH, COMP

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : np.ndarray
            MAPDL node numbers.

        thermal_strain : np.ndarray
            Nodal component plastic strains.  Array is in the order
            ``X, Y, Z, XY, YZ, XZ, EQV, ESWELL``

        Examples
        --------
        Load the nodal thermal strain for the first solution.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, thermal_strain = rst.nodal_thermal_strain(0)

        Return the nodal thermal strain just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, thermal_strain = rst.nodal_thermal_strain(0, nodes='MY_COMPONENT')

        Return the nodal thermal strain just for the nodes from 20 through 50.

        >>> nnum, thermal_strain = rst.nodal_thermal_strain(0, nodes=range(20, 51))
        """
        return self._nodal_result(rnum, "ETH", nnum_of_interest=nodes)

    def plot_nodal_thermal_strain(
        self,
        rnum,
        comp=None,
        scalar_bar_args={"title": "Nodal Thermal Strain"},
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal component thermal strains.

        Equivalent MAPDL command: PLNSOL, EPTH, COMP

        Parameters
        ----------
        rnum : int
            Result number

        comp : str, optional
            Thermal strain component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``
            - ``"EQV"``
            - ``"ESWELL"``

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``

        Examples
        --------
        Plot thermal strain for result 0 of verification manual example 33.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> result = examples.download_verification_result(33)
        >>> result.plot_nodal_thermal_strain(0)
        """
        return self._plot_nodal_result(
            rnum,
            "ETH",
            comp,
            THERMAL_STRAIN_TYPES,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            scalar_bar_args=scalar_bar_args,
            treat_nan_as_zero=treat_nan_as_zero,
            **kwargs,
        )

    def nodal_elastic_strain(self, rnum, nodes=None):
        """Nodal component elastic strains.  This record contains
        strains in the order ``X, Y, Z, XY, YZ, XZ, EQV``.

        Elastic strains can be can be nodal values extrapolated from
        the integration points or values at the integration points
        moved to the nodes.

        Equivalent MAPDL command: ``PRNSOL, EPEL``

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : np.ndarray
            MAPDL node numbers.

        elastic_strain : np.ndarray
            Nodal component elastic strains.  Array is in the order
            ``X, Y, Z, XY, YZ, XZ, EQV``.

        Examples
        --------
        Load the nodal elastic strain for the first result.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, elastic_strain = rst.nodal_elastic_strain(0)

        Return the nodal elastic strain just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, elastic_strain = rst.nodal_elastic_strain(0, nodes='MY_COMPONENT')

        Return the nodal elastic strain just for the nodes from 20 through 50.

        >>> nnum, elastic_strain = rst.nodal_elastic_strain(0, nodes=range(20, 51))

        Notes
        -----
        Nodes without a strain will be NAN.
        """
        return self._nodal_result(rnum, "EEL", nnum_of_interest=nodes)

    def plot_nodal_elastic_strain(
        self,
        rnum,
        comp,
        scalar_bar_args={"title": "Nodal Elastic Strain"},
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal elastic strain.

        Parameters
        ----------
        rnum : int
            Result number

        comp : str, optional
            Elastic strain component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``
            - ``"EQV"``

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT']``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT']``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``

        Examples
        --------
        Plot nodal elastic strain for a static pontoon model

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> result = examples.download_pontoon()
        >>> result.plot_nodal_elastic_strain(0)
        """
        scalar_bar_args["title"] = " ".join([comp.upper(), scalar_bar_args["title"]])
        return self._plot_nodal_result(
            rnum,
            "EEL",
            comp,
            STRAIN_TYPES,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            scalar_bar_args=scalar_bar_args,
            treat_nan_as_zero=treat_nan_as_zero,
            **kwargs,
        )

    def nodal_plastic_strain(self, rnum, nodes=None):
        """Nodal component plastic strains.

        This record contains strains in the order:
        ``X, Y, Z, XY, YZ, XZ, EQV``.

        Plastic strains are always values at the integration points
        moved to the nodes.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : np.ndarray
            MAPDL node numbers.

        plastic_strain : np.ndarray
            Nodal component plastic strains.  Array is in the order
            ``X, Y, Z, XY, YZ, XZ, EQV``.

        Examples
        --------
        Load the nodal plastic strain for the first solution.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> rst = pymapdl_reader.read_binary('file.rst')
        >>> nnum, plastic_strain = rst.nodal_plastic_strain(0)

        Return the nodal plastic strain just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, plastic_strain = rst.nodal_plastic_strain(0, nodes='MY_COMPONENT')

        Return the nodal plastic strain just for the nodes from 20
        through 50.

        >>> nnum, plastic_strain = rst.nodal_plastic_strain(0, nodes=range(20, 51))

        """
        return self._nodal_result(rnum, "EPL", nnum_of_interest=nodes)

    def plot_nodal_plastic_strain(
        self,
        rnum,
        comp,
        scalar_bar_args={"title": "Nodal Plastic Strain"},
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal component plastic strain.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Plastic strain component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``
            - ``"EQV"``

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        element_components : list, optional
            Accepts either a string or a list strings of element
            components to plot.  For example:
            ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``

        sel_type_all : bool, optional
            If node_components is specified, plots those elements
            containing all nodes of the component.  Default ``True``.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        **kwargs : keyword arguments
            Optional keyword arguments.  See ``help(pyvista.plot)``.

        Examples
        --------
        Plot plastic strain for a static pontoon model

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> result = examples.download_pontoon()
        >>> result.plot_nodal_plastic_strain(0)
        """
        scalar_bar_args["title"] = " ".join([comp.upper(), scalar_bar_args["title"]])
        return self._plot_nodal_result(
            rnum,
            "EPL",
            comp,
            STRAIN_TYPES,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            node_components=node_components,
            element_components=element_components,
            sel_type_all=sel_type_all,
            scalar_bar_args=scalar_bar_args,
            treat_nan_as_zero=treat_nan_as_zero,
            **kwargs,
        )

    def _plot_nodal_result(
        self,
        rnum,
        result_type,
        comp,
        available_comps,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=True,
        **kwargs,
    ):
        """Plot nodal result"""
        component_index = check_comp(available_comps, comp)
        _, result = self._nodal_result(rnum, result_type)
        scalars = result[:, component_index]

        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            scalars = scalars[ind]
        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            scalars = scalars[ind]
        else:
            grid = self.grid

        return self._plot_point_scalars(
            scalars,
            grid=grid,
            rnum=rnum,
            show_displacement=show_displacement,
            displacement_factor=displacement_factor,
            treat_nan_as_zero=treat_nan_as_zero,
            sel_type_all=sel_type_all,
            **kwargs,
        )

    def _animate_time_solution(
        self,
        result_type,
        index=0,
        frame_rate=10,
        show_displacement=True,
        displacement_factor=1,
        off_screen=None,
    ):
        """Animate time solution result"""
        # load all results
        results = []
        for i in range(self.nsets):
            results.append(self._nodal_result(i, result_type)[1][:, index])

        if show_displacement:
            disp = []
            for i in range(self.nsets):
                disp.append(self.nodal_solution(i)[1][:, :3] * displacement_factor)

        mesh = self.grid.copy()
        results = np.array(results)
        if np.all(np.isnan(results)):
            raise ValueError(
                "Result file contains no %s records"
                % element_index_table_info[result_type.upper()]
            )

        # prepopulate mesh with data
        mesh["data"] = results[0]

        # set default range
        rng = [results.min(), results.max()]
        t_wait = 1 / frame_rate

        def q_callback():
            """exit when user wants to leave"""
            self._animating = False

        self._animating = True

        def plot_thread():
            plotter = pv.Plotter(off_screen=off_screen)
            plotter.add_key_event("q", q_callback)
            plotter.add_mesh(mesh, scalars="data", rng=rng)
            plotter.show(auto_close=False, interactive_update=True, interactive=False)
            text_actor = plotter.add_text("Result 1")
            while self._animating:
                for i in range(self.nsets):
                    mesh["data"] = results[i]

                    if show_displacement:
                        mesh.points = self.grid.points + disp[i]

                    # if interactive:
                    plotter.update(30, force_redraw=True)
                    if hasattr(text_actor, "SetInput"):
                        text_actor.SetInput("Result %d" % (i + 1))
                    else:
                        text_actor.SetText(0, "Result %d" % (i + 1))

                    time.sleep(t_wait)

                if off_screen:
                    break

            plotter.close()

        thread = Thread(target=plot_thread)
        thread.start()

    @property
    def available_results(self):
        """Available result types.

        Examples
        --------
        >>> rst.available_results
        Available Results:
        ENS : Nodal stresses
        ENG : Element energies and volume
        EEL : Nodal elastic strains
        EPL : Nodal plastic strains
        ETH : Nodal thermal strains (includes swelling strains)
        EUL : Element euler angles
        ENL : Nodal nonlinear items, e.g. equivalent plastic strains
        EPT : Nodal temperatures
        NSL : Nodal displacements
        RF  : Nodal reaction forces
        """
        return self._available_results

    def nodal_static_forces(self, rnum, nodes=None):
        """Return the nodal forces averaged at the nodes.

        Nodal forces are computed on an element by element basis, and
        this method averages the nodal forces for each element for
        each node.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        nodes : str, sequence of int or str, optional
            Select a limited subset of nodes.  Can be a nodal
            component or array of node numbers.  For example

            * ``"MY_COMPONENT"``
            * ``['MY_COMPONENT', 'MY_OTHER_COMPONENT]``
            * ``np.arange(1000, 2001)``

        Returns
        -------
        nnum : np.ndarray
            MAPDL node numbers.

        forces : np.ndarray
           Averaged nodal forces.  Array is sized ``[nnod x numdof]``
           where ``nnod`` is the number of nodes and ``numdof`` is the
           number of degrees of freedom for this solution.

        Examples
        --------
        Load the nodal static forces for the first result using the
        example hexahedral result file.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> rst = pymapdl_reader.read_binary(examples.rstfile)
        >>> nnum, forces = rst.nodal_static_forces(0)

        Return the nodal static forces just for the nodal component
        ``'MY_COMPONENT'``.

        >>> nnum, forces = rst.nodal_static_forces(0, nodes='MY_COMPONENT')

        Return the nodal static forces just for the nodes from 20 through 50.

        >>> nnum, forces = rst.nodal_static_forces(0, nodes=range(20, 51))

        Notes
        -----
        Nodes without a a nodal will be NAN.  These are generally
        midside (quadratic) nodes.
        """
        return self._nodal_result(rnum, "ENF", nnum_of_interest=nodes)

    @property
    def _is_distributed(self):
        """True when this result file is part of a distributed result

        Only True when Global number of nodes does not equal the
        number of nodes in this file.

        Notes
        -----
        Not a reliabile indicator if a cyclic result.
        """
        return self._resultheader["Glbnnod"] != self._resultheader["nnod"]

    def _is_main(self):
        """Result file written from the main process"""
        return bool(self._resultheader["ptrGNOD"])

    @property
    def _is_thermal(self):
        """True when result file is a rth file"""
        return self.pathlib_filename.suffix == ".rth"

    @property
    def _is_cyclic(self):
        """True when the result file is from a cyclic analysis"""
        return self.n_sector > 1


def pol2cart(rho, phi):
    """Convert cylindrical to cartesian"""
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def is_int(value):
    """Return true if can be parsed as an int"""
    try:
        int(value)
        return True
    except:
        return False


def trans_to_matrix(trans):
    """Convert a numpy.ndarray to a vtk.vtkMatrix4x4"""
    matrix = vtk.vtkMatrix4x4()
    for i in range(trans.shape[0]):
        for j in range(trans.shape[1]):
            matrix.SetElement(i, j, trans[i, j])
    return matrix


def transform(points, trans):
    """In-place 3d transformation of a points array given a 4x4
    transformation matrix.

    Parameters
    ----------
    points : np.ndarray or vtk.vtkTransform
        Points to transform.

    transform : np.ndarray or vtk.vtkTransform
        4x4 transformation matrix.

    """
    if isinstance(trans, vtk.vtkMatrix4x4):
        trans = pv.array_from_vtkmatrix(trans)

    _binary_reader.affline_transform(points, trans)


def merge_two_dicts(x, y):
    merged = x.copy()  # start with x's keys and values
    merged.update(y)  # modifies z with y's keys and values & returns None
    return merged


def check_comp(available_comps, comp):
    """Check if a component is in available components and return a
    helpful error message

    Parameters
    ----------
    available_comps : list
        List of available components.  For example:
        ``['R', 'THETA', 'Z', 'RTHETA', 'THETAZ', 'RZ']``
        Case should be all caps

    comp : str
        Component to check.  Any case.

    Returns
    -------
    idx : int
        Index in ``available_comps``.

    """
    if comp is None:
        raise ValueError(
            'Missing "comp" parameter.  Please select'
            " from the following:\n%s" % available_comps
        )
    comp = comp.upper()
    if comp not in available_comps:
        raise ValueError(
            'Invalid "comp" parameter %s.  Please select' % comp
            + " from the following:\n%s" % available_comps
        )
    return available_comps.index(comp)
