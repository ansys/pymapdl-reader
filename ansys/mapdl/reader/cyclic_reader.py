"""Supports reading cyclic structural result files from ANSYS"""

from collections.abc import Iterable
from functools import wraps
import warnings

import numpy as np
import pyvista as pv
from tqdm import tqdm
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersCore import vtkAppendFilter

from ansys.mapdl.reader import _binary_reader
from ansys.mapdl.reader._rst_keys import element_index_table_info
from ansys.mapdl.reader.common import (
    PRINCIPAL_STRESS_TYPES,
    STRAIN_TYPES,
    STRESS_TYPES,
    THERMAL_STRAIN_TYPES,
    axis_rotation,
)
from ansys.mapdl.reader.rst import ELEMENT_INDEX_TABLE_KEYS, Result, check_comp

np.seterr(divide="ignore", invalid="ignore")

# MAPDL results that are tensors
RESULT_TENSORS_TYPES = [
    "ENS",  # nodal stresses
    "EEL",  # elastic strains
    "EPL",  # plastic strains
    "ECR",  # creep strains
    "ETH",  # thermal strains
    "EDI",  # diffusion strains
    "EBA",  # back stresses
]

# MAPDL results that are stresses
RESULT_STRESS_TYPES = [
    "ENS",  # nodal stresses
    "ESF",  # element surface stresses
    "EBA",  # back stresses
]


class CyclicResult(Result):
    """Adds cyclic functionality to the result class"""

    def __init__(self, filename, read_mesh: bool = True) -> None:
        """Initializes cyclic result"""
        super().__init__(filename, read_mesh=read_mesh)

        # sanity check
        if not self._is_cyclic:
            raise TypeError("Result is not a cyclic model")

        self._rotor_cache = None
        self._has_duplicate_sector = None
        self._is_repeated_mode = np.empty(0)
        self._repeated_index = np.empty(0)

        if read_mesh:
            self._add_cyclic_properties()

    def plot_sectors(self, **kwargs):
        """Plot the full rotor and individually color the sectors.

        Parameters
        ----------
        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Examples
        --------
        >>> from ansys.mapdl.reader import examples
        >>> rst = examples.download_academic_rotor()
        >>> rst.plot_sectors()

        Save a screenshot of the sectors

        >>> rst.plot_sectors(screenshot='sectors.png')

        """
        scalars = np.empty((self.n_sector, self._mas_grid.n_points))
        scalars[:] = np.arange(self.n_sector).reshape(-1, 1)
        kwargs.setdefault("show_edges", True)
        kwargs.setdefault("n_colors", self.n_sector)
        return self._plot_cyclic_point_scalars(scalars, None, add_text=False, **kwargs)

    def plot(self, **kwargs):
        """Plot the full rotor geometry.

        Parameters
        ----------
        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Returns
        -------
        cpos : list
            List of camera position, focal point, and view up.

        Examples
        --------
        >>> from ansys.mapdl.reader import examples
        >>> rst = examples.download_academic_rotor()
        >>> rst.plot()

        Save a screenshot of the rotor

        >>> rst.plot(screenshot='rotor.png')
        """
        kwargs.setdefault("color", "w")
        kwargs.setdefault("show_edges", True)
        return self._plot_cyclic_point_scalars(None, None, add_text=False, **kwargs)

    def _add_cyclic_properties(self):
        """Add cyclic properties"""

        # identity the sector based on number of elements in master sector
        cs_els = self._resultheader["csEls"]
        mask = self.quadgrid.cell_data["ansys_elem_num"] <= cs_els

        self.master_cell_mask = mask
        self._mas_grid = self.grid.extract_cells(mask)

        # NOTE: number of nodes in sector may not match number of
        # nodes in geometry
        node_mask = self._neqv[self._sidx] <= self._resultheader["csNds"]
        self._mas_ind = np.nonzero(node_mask)[0]
        self._dup_ind = np.nonzero(~node_mask)[0]
        self._has_duplicate_sector = np.any(self._dup_ind)

        # determine repeated modes
        mask_a = np.isclose(self.time_values, np.roll(self.time_values, 1))
        mask_b = np.isclose(self.time_values, np.roll(self.time_values, -1))
        self._is_repeated_mode = np.logical_or(mask_a, mask_b)

        # edge case for single pair of repeated modes
        if self._is_repeated_mode.size == 2 and self._is_repeated_mode.all():
            self._repeated_index = np.array([1, 0])
            return
        elif self._is_repeated_mode.size == 1:  # edge case single result
            self._is_repeated_mode = np.array([False])
            return

        self._repeated_index = np.empty(self._is_repeated_mode.size, np.int_)
        self._repeated_index[:] = -1
        if np.any(self._is_repeated_mode):
            self._repeated_index[mask_a] = np.nonzero(mask_b)[0]
            self._repeated_index[mask_b] = np.nonzero(mask_a)[0]

    def nodal_solution(self, rnum, phase=0, full_rotor=False, as_complex=False):
        """Return the DOF solution for each node in the global cartesian coordinate system.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float, optional
            Phase to rotate sector result in radians.

        full_rotor : bool, optional
            Expands the single sector solution for the full rotor.
            Sectors are rotated counter-clockwise about the axis of
            rotation.  Default False.

        as_complex : bool, optional
            Returns result as a complex number, otherwise as the real
            part rotated by phase.  Default False.

        Returns
        -------
        nnum : numpy.ndarray
            Node numbers of master sector.

        result : numpy.ndarray
            Result is (nnod x numdof), nnod is the number of nodes in
            a sector and numdof is the number of degrees of freedom.
            When full_rotor is True the array will be (nSector x nnod
            x numdof).

        Examples
        --------
        Visualize the 1st nodal diameter mode.

        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> result = examples.download_academic_rotor()
        >>> result.nodal_solution((2, 1))

        Same result but uses Python (zero based) cumulative indexing

        >>> result.nodal_solution(2)

        Notes
        -----
        Somewhere between v15.0 and v18.2 ANSYS stopped writing the
        duplicate sector to the result file and instead records results in
        pairs (i.e. harmonic index 1, -1).
        """

        func = super().nodal_solution
        return self._get_full_result(
            rnum, func, phase, full_rotor, as_complex, tensor=False
        )

    @wraps(nodal_solution)
    def nodal_displacement(self, *args, **kwargs):
        """wraps nodal_solution"""
        return self.nodal_solution(*args, **kwargs)

    def _expand_cyclic_static(self, result, tensor=False, stress=True):
        """Expand cyclic static result for a full rotor"""
        cs_cord = self._resultheader["csCord"]
        if cs_cord > 1:
            matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix.Invert()
        else:
            matrix = vtkMatrix4x4()
            i_matrix = vtkMatrix4x4()

        shp = (self.n_sector, result.shape[0], result.shape[1])
        full_result = np.empty(shp)
        full_result[:] = result

        rang = 360.0 / self.n_sector
        for i in range(1, self.n_sector):
            # transform to standard position, rotate about Z axis,
            # transform back
            transform = vtkTransform()
            transform.RotateZ(rang * i)
            transform.Update()
            rot_matrix = transform.GetMatrix()

            if cs_cord > 1:
                temp_matrix = vtkMatrix4x4()
                rot_matrix.Multiply4x4(i_matrix, rot_matrix, temp_matrix)
                rot_matrix.Multiply4x4(temp_matrix, matrix, rot_matrix)

            trans = pv.array_from_vtkmatrix(rot_matrix)
            if tensor:
                if stress:
                    _binary_reader.tensor_arbitrary(full_result[i], trans)
                else:
                    _binary_reader.tensor_strain_arbitrary(full_result[i], trans)
            else:
                _binary_reader.affline_transform(full_result[i], trans)

        return full_result

    def _expand_cyclic_modal(
        self, result, result_r, hindex, phase, as_complex, full_rotor
    ):
        """Combines repeated results from ANSYS"""
        if as_complex or full_rotor:
            # Matches ansys direction
            # if self._positive_cyclic_dir:
            result_combined = result + result_r * 1j
            # else:
            # result_combined = result - result_r*1j
            if phase:
                result_combined *= 1 * np.cos(phase) - 1j * np.sin(phase)
        else:  # convert to real
            result_combined = result * np.cos(phase) - result_r * np.sin(phase)

        # just return single sector
        if not full_rotor:
            return result_combined

        # Generate full rotor solution
        result_expanded = []
        angles = np.linspace(0, 2 * np.pi, self.n_sector + 1)[:-1] + phase
        for angle in angles:
            # need to rotate solution and rotate direction
            result_expanded.append(
                axis_rotation(result_combined, angle, deg=False, axis="z")
            )

        result_expanded = np.asarray(result_expanded, dtype=np.complex128)

        # ANSYS scales the result
        if hindex == 0 or hindex == self.n_sector / 2:
            result_expanded /= self.n_sector**0.5
        else:
            result_expanded /= (self.n_sector / 2) ** 0.5

        # adjust phase of the full result based on the harmonic index
        f_arr = np.zeros(self.n_sector)
        f_arr[hindex] = 1
        jang = np.fft.ifft(f_arr)[: self.n_sector] * self.n_sector
        cjang = jang * (np.cos(phase) - np.sin(phase) * 1j)  # 14-233

        result_expanded *= cjang.reshape(-1, 1, 1)
        if as_complex:
            return result_expanded
        return np.real(result_expanded)

    def _expand_cyclic_modal_tensor(
        self, result, result_r, hindex, phase, as_complex, full_rotor, stress=True
    ):
        """Combines repeated results from ANSYS and optionally
        duplicates/rotates it around the axis of rotation"""

        # must scale value

        if as_complex or full_rotor:
            # if self._positive_cyclic_dir:
            result_combined = result + result_r * 1j
            # else:
            # result_combined = result - result_r*1j

            if phase:
                result_combined *= 1 * np.cos(phase) - 1j * np.sin(phase)

        else:  # convert to real
            result_combined = result * np.cos(phase) - result_r * np.sin(phase)

        # just return single sector
        if not full_rotor:
            return result_combined

        # Generate full rotor solution
        shp = (self.n_sector, result.shape[0], result.shape[1])
        result_expanded = np.empty(shp, np.complex128)
        result_expanded[:] = result_combined

        # convert hindex to nodal content
        f_arr = np.zeros(self.n_sector)
        f_arr[hindex] = self.n_sector

        jang = np.fft.ifft(f_arr)[: self.n_sector]

        # must be adjusted of to expand
        if hindex == 0 or hindex == self.n_sector / 2:
            jang *= self.n_sector ** (-1 / 2)
        else:
            jang *= 2 * (2 * self.n_sector) ** (-1 / 2)

        cjang = jang * (np.cos(phase) - np.sin(phase) * 1j)

        full_result = np.empty(shp)
        full_result[:] = np.real(result_expanded * cjang.reshape(-1, 1, 1))

        cs_cord = self._resultheader["csCord"]
        if cs_cord > 1:
            matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix.Invert()
        else:
            matrix = vtkMatrix4x4()
            i_matrix = vtkMatrix4x4()

        rang = 360.0 / self.n_sector
        for i in range(1, self.n_sector):
            # transform to standard position, rotate about Z axis,
            # transform back
            transform = vtkTransform()
            transform.RotateZ(rang * i)
            transform.Update()
            rot_matrix = transform.GetMatrix()

            if cs_cord > 1:
                temp_matrix = vtkMatrix4x4()
                rot_matrix.Multiply4x4(i_matrix, rot_matrix, temp_matrix)
                rot_matrix.Multiply4x4(temp_matrix, matrix, rot_matrix)

            trans = pv.array_from_vtkmatrix(rot_matrix)
            if stress:
                _binary_reader.tensor_arbitrary(full_result[i], trans)
            else:
                _binary_reader.tensor_strain_arbitrary(full_result[i], trans)

        return full_result

    def harmonic_index_to_cumulative(self, hindex, mode):
        """Converts a harmonic index and a 0 index mode number to a
        cumulative result index.

        Harmonic indices are stored as positive and negative pairs for
        modes other than 0 and N/nsectors.

        Parameters
        ----------
        hindex : int
            Harmonic index.  Must be less than or equal to nsectors/2.
            May be positive or negative

        mode : int
            Mode number.  0 based indexing.  Access mode pairs by with
            a negative/positive harmonic index.

        Returns
        -------
        rnum : int
            Cumulative index number.  Zero based indexing.

        """
        hindex_table = self._resultheader["hindex"]
        if not np.any(abs(hindex) == np.abs(hindex_table)):
            raise Exception(
                "Invalid harmonic index.\n"
                + "Available indices: %s" % np.unique(hindex_table)
            )

        mask = np.logical_and(hindex == hindex_table, mode == self.mode_table)

        if not mask.any():
            mode_mask = abs(hindex) == np.abs(hindex_table)
            avail_modes = np.unique(self.mode_table[mode_mask])
            raise Exception(
                "Invalid mode for harmonic index %d\n" % hindex
                + "Available modes: %s" % avail_modes
            )

        index = mask.nonzero()[0]
        assert index.size == 1, "Multiple cumulative index matches"
        return index[0]

    @property
    def mode_table(self):
        """Unique modes for cyclic results"""
        hindex_table = self._resultheader["hindex"]
        diff = np.diff(np.abs(hindex_table))
        freqs = self.time_values
        mode_table = [0]
        c = 0
        for i in range(1, freqs.size):
            if diff[i - 1]:
                c = 0
                mode_table.append(c)
            elif np.isclose(freqs[i], freqs[i - 1]):
                mode_table.append(c)
            else:
                c += 1
                mode_table.append(c)
        return np.asarray(mode_table)

    @property
    def harmonic_indices(self):
        """Harmonic indices of the result file.

        Harmonic index is simply the Nodal Diameter of the mode.  This
        is defined as the number of complete sine waves that pass
        through the circumference.

        Examples
        --------
        >>> rst.harmonic_indices
        array([ 0,  0,  0,  0,  0,  0, -1,  1, -1,  1,  1, -1,
               -2,  2, -2,  2, -2,  2,  3,  3,  3,  3,  3,  3], dtype=int32)
        """
        return self._resultheader["hindex"]

    def nodal_stress(self, rnum, phase=0, as_complex=False, full_rotor=False):
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

        phase : float
            Phase adjustment of the stress in degrees.

        as_complex : bool, optional
            Reports stress as a complex result.  Real and imaginary
            stresses correspond to the stress of the main and repeated
            sector.  Stress can be "rotated" using the phase
            parameter.

        full_rotor : bool, optional
            Expands the results to the full rotor when True.  Default
            False.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        stress : numpy.ndarray
            Stresses at Sx Sy Sz Sxy Syz Sxz averaged at each corner
            node.  For the corresponding node numbers, see where
            result is the result object.

        Examples
        --------
        >>> nnum, stress = rst.nodal_stress(0)

        Notes
        -----
        Nodes without a stress value will be NAN.

        """
        func = super().nodal_stress
        return self._get_full_result(
            rnum, func, phase, full_rotor, as_complex, tensor=True, stress=True
        )

    def _get_full_result(
        self, rnum, func, phase, full_rotor, as_complex, tensor=True, stress=False
    ):
        """Return the full rotor result or the complex result for a cyclic model.

        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float
            Phase adjustment of the stress in radians.

        tensor : bool
            True when the result is a stress/strain tensor.  False
            when a scalar or displacement value.

        stress : bool
            True when tensor is a stress.  False when tensor is a
            strain.  Ignored when not a tensor.

        """
        rnum = self.parse_step_substep(rnum)
        nnum, full_result = func(rnum)

        # full result may or may not contain the duplicate sector
        if self._has_duplicate_sector:
            if nnum.size < self._mas_ind.size:
                # node numbers of the master sector
                nnum_mas = self._neqv[self._mas_ind]
                mask = np.isin(nnum, nnum_mas)
                nnum = nnum[mask]
                result = full_result[mask]
            else:
                result = full_result[self._mas_ind]
                nnum = nnum[self._mas_ind]
        else:
            result = full_result

        if self._resultheader["kan"] == 0:  # static result
            if full_rotor:
                expanded_result = self._expand_cyclic_static(
                    result, tensor=tensor, stress=stress
                )
            else:
                return nnum, result
        elif self._resultheader["kan"] == 2:  # modal analysis
            # combine modal solution results
            hindex_table = self._resultheader["hindex"]
            hindex = hindex_table[rnum]  # move this to expand_modal_tensor
            result_r = self._get_complex_result(func, rnum, result)
            if tensor:
                expanded_result = self._expand_cyclic_modal_tensor(
                    result,
                    result_r,
                    hindex,
                    phase,
                    as_complex,
                    full_rotor,
                    stress=stress,
                )
            else:
                assert result.shape == result_r.shape
                expanded_result = self._expand_cyclic_modal(
                    result, result_r, hindex, phase, as_complex, full_rotor
                )
        else:
            raise RuntimeError("Unsupported analysis type")

        return nnum, expanded_result

    def _get_complex_result(self, func, rnum, full_result):
        """Acquire the duplicate sector or repeated result.

        Depending on the version of MAPDL, this may mean using the
        result from the duplicate sector or the mode pair when there
        are duplicate modes.

        When there is no repeated mode or duplicate sector, returns an
        all zero array.

        Parameters
        ----------
        func : function
            Function to acquire the sector only result.

        rnum : int
            Cumulative result number.

        full_result : np.ndarray
            Full result (may include duplicate sector).

        Returns
        --------
        result_r : np.ndarray
            Repeated result
        """
        has_dup_result = False
        if self._has_duplicate_sector:
            has_dup_result = self._dup_ind[-1] <= full_result.shape[0] - 1

        if self._is_repeated_mode[rnum]:  # get mode pair result
            _, result_r = func(self._repeated_index[rnum])
            if result_r.shape[0] != self._mas_ind.size:
                result_r = result_r[self._mas_ind]
        elif has_dup_result:  # use the duplicate sector
            result_r = full_result[self._dup_ind]
        else:  # otherwise, a standing wave (no complex component)
            result_r = np.zeros(
                (self._mas_ind.size, full_result.shape[1]), dtype=full_result.dtype
            )

        return result_r

    def nodal_thermal_strain(self, rnum, phase=0, as_complex=False, full_rotor=False):
        """Nodal component thermal strains.  This record contains
        strains in the order X, Y, Z, XY, YZ, XZ, EQV, and eswell
        (element swelling strain).  Thermal strains are always values
        at the integration points moved to the nodes.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float
            Phase adjustment of the stress in degrees.

        as_complex : bool, optional
            Reports stress as a complex result.  Real and imaginary
            stresses correspond to the stress of the main and repeated
            sector.  Stress can be "rotated" using the phase
            parameter.

        full_rotor : bool, optional
            Expands the results to the full rotor when True.  Default
            False.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        thermal_strain : np.ndarray
            Nodal component plastic strains.  Array is in the order
            X, Y, Z, XY, YZ, XZ, EQV, ESWELL

        Examples
        --------
        Load the nodal thermal strain for the first result.

        >>> nnum, thermal_strain = rst.nodal_thermal_strain(0)

        Notes
        -----
        Nodes without a strain will be NAN.

        """
        func = super().nodal_thermal_strain
        return self._get_full_result(
            rnum, func, phase, full_rotor, as_complex, tensor=True, stress=False
        )

    def plot_nodal_thermal_strain(
        self,
        rnum,
        comp=None,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        add_text=True,
        overlay_wireframe=False,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot nodal thermal strain.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Thermal strain component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``
            - ``"EQV"``
            - ``"ESWELL"`` (element swelling strain)

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0.

        full_rotor : bool, optional
            Expand the sector solution to the full rotor.

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot nodal thermal strain for an academic rotor

        >>> rst.plot_nodal_thermal_strain(0)

        """
        if not full_rotor:
            return super().plot_nodal_thermal_strain(
                rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
            )

        idx = check_comp(THERMAL_STRAIN_TYPES, comp)
        _, strain = self.nodal_thermal_strain(rnum, phase, False, True)
        scalars = strain[:, :, idx]

        kwargs.setdefault("scalar_bar_args", {"title": f"{comp} Nodal Thermal Strain"})
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def nodal_elastic_strain(self, rnum, phase=0, as_complex=False, full_rotor=False):
        """Nodal component elastic strains.  This record contains
        strains in the order X, Y, Z, XY, YZ, XZ, EQV.

        Elastic strains can be can be nodal values extrapolated from
        the integration points or values at the integration points
        moved to the nodes.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float
            Phase adjustment of the stress in radians.

        as_complex : bool, optional
            Reports stress as a complex result.  Real and imaginary
            stresses correspond to the stress of the main and repeated
            sector.  Stress can be "rotated" using the phase
            parameter.

        full_rotor : bool, optional
            Expands the results to the full rotor when True.  Default
            False.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        elastic_strain : numpy.ndarray
            Nodal component elastic strains.  Array is in the order
            X, Y, Z, XY, YZ, XZ, EQV.

        Examples
        --------
        Load the nodal elastic strain for the first result.

        >>> nnum, elastic_strain = rst.nodal_stress(0)

        Notes
        -----
        Nodes without a strain will be NAN.

        """
        func = super().nodal_elastic_strain
        return self._get_full_result(
            rnum, func, phase, full_rotor, as_complex, tensor=True, stress=False
        )

    def plot_nodal_elastic_strain(
        self,
        rnum,
        comp=None,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        add_text=True,
        overlay_wireframe=False,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot nodal elastic strain.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Elastic strain component to display.  Available options:
            - ``"X"``
            - ``"Y"``
            - ``"Z"``
            - ``"XY"``
            - ``"YZ"``
            - ``"XZ"``

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

        full_rotor : bool, optional
            Expand the sector solution to the full rotor.

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot nodal elastic strain for an academic rotor.

        >>> result.plot_nodal_elastic_strain(0, 'X')

        """
        idx = check_comp(STRAIN_TYPES[:-1], comp)
        _, strain = self.nodal_elastic_strain(rnum, phase, False, full_rotor)
        scalars = strain[:, :, idx]

        kwargs.setdefault("scalar_bar_args", {"title": f"{comp} Nodal Elastic Strain"})
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def nodal_plastic_strain(self, rnum, phase=0, as_complex=False, full_rotor=False):
        """Nodal component plastic strains.  This record contains
        strains in the order X, Y, Z, XY, YZ, XZ, EQV.

        Plastic strains can be can be nodal values extrapolated from
        the integration points or values at the integration points
        moved to the nodes.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float
            Phase adjustment of the stress in degrees.

        as_complex : bool, optional
            Reports stress as a complex result.  Real and imaginary
            stresses correspond to the stress of the main and repeated
            sector.  Stress can be "rotated" using the phase
            parameter.

        full_rotor : bool, optional
            Expands the results to the full rotor when True.  Default
            False.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        plastic_strain : numpy.ndarray
            Nodal component plastic strains.  Array is in the order
            X, Y, Z, XY, YZ, XZ, EQV.

        Examples
        --------
        Load the nodal plastic strain for the first result.

        >>> nnum, plastic_strain = rst.nodal_stress(0)

        Notes
        -----
        Nodes without a strain will be NAN.

        """
        func = super().nodal_plastic_strain
        return self._get_full_result(
            rnum, func, phase, full_rotor, as_complex, tensor=True, stress=False
        )

    def plot_nodal_plastic_strain(
        self,
        rnum,
        comp=None,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        add_text=True,
        overlay_wireframe=False,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot nodal plastic strain.

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

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

        full_rotor : bool, optional
            Expand the sector solution to the full rotor.

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot nodal plastic strain for an academic rotor

        >>> result.plot_nodal_plastic_strain(0)

        """
        idx = check_comp(STRAIN_TYPES[:-1], comp)
        _, strain = self.nodal_plastic_strain(rnum, phase, False, full_rotor)
        scalars = strain[:, :, idx]

        kwargs.setdefault("scalar_bar_args", {"title": f"{comp} Nodal Plastic Strain"})
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def principal_nodal_stress(self, rnum, phase=0, as_complex=False, full_rotor=False):
        """Computes the principal component stresses for each node in
        the solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        phase : float
            Phase adjustment of the stress in degrees.

        as_complex : bool, optional
            Returns result as a complex number, otherwise as the real
            part rotated by phase.  Default False.

        full_rotor : bool, optional
            Expand sector solution to full rotor.

        Returns
        -------
        nodenum : numpy.ndarray
            Node numbers of the result.

        pstress : numpy.ndarray
            Principal stresses, stress intensity, and equivalent stress.
            [sigma1, sigma2, sigma3, sint, seqv]

        Notes
        -----
        ANSYS equivalent of:
        PRNSOL, S, PRIN

        which returns:
        S1, S2, S3 principal stresses, SINT stress intensity, and SEQV
        equivalent stress.
        """
        if as_complex and full_rotor:
            raise ValueError("`complex` and `full_rotor` cannot both be True")

        # get component stress
        nnum, stress = self.nodal_stress(rnum, phase, as_complex, full_rotor)

        # compute principle stress
        if as_complex:
            stress_r = np.imag(stress)
            stress = np.real(stress)

            if not stress.flags["C_CONTIGUOUS"]:
                stress = np.ascontiguousarray(stress)
            if not stress_r.flags["C_CONTIGUOUS"]:
                stress_r = np.ascontiguousarray(stress_r)

            pstress, isnan = _binary_reader.compute_principal_stress(stress)
            pstress[isnan] = np.nan
            pstress_r, isnan = _binary_reader.compute_principal_stress(stress_r)
            pstress_r[isnan] = np.nan

            return nnum, pstress + 1j * pstress_r

        if full_rotor:
            # compute principle stress for each sector
            pstress = np.empty((self.n_sector, stress.shape[1], 5), np.float64)
            for i in range(stress.shape[0]):
                pstress[i], isnan = _binary_reader.compute_principal_stress(stress[i])
                pstress[i, isnan] = np.nan
            return nnum, pstress

        pstress, isnan = _binary_reader.compute_principal_stress(stress)
        pstress[isnan] = np.nan
        return nnum, pstress

    def plot_nodal_solution(
        self,
        rnum,
        comp="norm",
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1.0,
        node_components=None,
        element_components=None,
        overlay_wireframe=False,
        add_text=True,
        sel_type_all=True,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot the nodal solution (generally displacement).

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Display component to display.  Options are 'x', 'y', 'z',
            and 'norm', corresponding to the x direction, y direction,
            z direction, and the normalized direction:
            ``(x**2 + y**2 + z**2)**0.5``

        full_rotor : bool, optional
            Expand sector solution to full rotor.

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot the displacement of the first cyclic result.

        >>> result.plot_nodal_solution(0)

        """
        # Load result from file
        if not full_rotor:
            return super().plot_nodal_solution(
                rnum,
                comp,
                show_displacement=show_displacement,
                displacement_factor=displacement_factor,
                node_components=node_components,
                element_components=element_components,
                sel_type_all=sel_type_all,
                treat_nan_as_zero=treat_nan_as_zero,
                **kwargs,
            )

        rnum = self.parse_step_substep(rnum)
        _, result = self.nodal_solution(rnum, phase, full_rotor, as_complex=False)

        # Process result
        label = "Cyclic Rotor\nDisplacement"
        if comp == "x":
            scalars = result[:, :, 0]
            title = f"X {label}\n"
        elif comp == "y":
            scalars = result[:, :, 1]
            title = f"Y {label}\n"
        elif comp == "z":
            scalars = result[:, :, 2]
            title = f"Z {label}\n"
        else:
            # Normalize displacement
            scalars = np.linalg.norm(result, axis=-1)
            title = f"Normalized\n{label}\n"

        kwargs.setdefault("scalar_bar_args", {"title": title})
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def plot_nodal_stress(
        self,
        rnum,
        comp=None,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        overlay_wireframe=False,
        add_text=True,
        sel_type_all=True,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot nodal stress of a given component

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

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

        full_rotor : bool, optional
            Expand the sector solution to the full rotor.

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot the ``"Z"`` nodal stress of the first cyclic result.

        >>> result.plot_nodal_stress(0, comp="Z")

        """
        if not full_rotor:
            super().plot_nodal_stress(
                rnum,
                comp,
                show_displacement,
                displacement_factor,
                node_components,
                element_components,
                sel_type_all,
                treat_nan_as_zero=treat_nan_as_zero,
                **kwargs,
            )

        idx = check_comp(STRESS_TYPES, comp)
        _, stress = self.nodal_stress(rnum, phase, False, full_rotor=True)
        scalars = stress[:, :, idx]

        kwargs.setdefault(
            "scalar_bar_args", {"title": f"Cyclic Rotor\nNodal Stress\n{comp}\n"}
        )
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def plot_principal_nodal_stress(
        self,
        rnum,
        comp=None,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        add_text=True,
        overlay_wireframe=False,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot the nodal principal stress.

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

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

        full_rotor : bool, optional
            Expand sector solution to full rotor.

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        node_components : list, optional
            Accepts either a string or a list strings of node
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
        Plot the von Mises stress of the first cyclic result.

        >>> result.plot_principal_nodal_stress(0, comp='SEQV')

        """
        if not full_rotor:
            return super().plot_principal_nodal_stress(
                rnum,
                comp,
                show_displacement,
                displacement_factor,
                node_components,
                element_components,
                sel_type_all,
                treat_nan_as_zero=treat_nan_as_zero,
                **kwargs,
            )

        # get the correct component of the principal stress for the rotor
        idx = check_comp(PRINCIPAL_STRESS_TYPES, comp)
        _, pstress = self.principal_nodal_stress(rnum, phase, full_rotor=True)
        scalars = pstress[:, :, idx]

        kwargs.setdefault(
            "scalar_bar_args",
            {"title": f"Cyclic Rotor\nPrincipal Nodal Stress\n{comp}\n"},
        )
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        self._plot_cyclic_point_scalars(
            scalars, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def nodal_temperature(self, rnum, full_rotor=False):
        """Retrieves the temperature for each node in the solution.

        The order of the results corresponds to the sorted node
        numbering.

        Equivalent MAPDL commands:
        PRNSOL, TEMP
        PRNSOL, BFE

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        full_rotor : bool, optional
            Expand sector solution to full rotor.

        Returns
        -------
        nnum : numpy.ndarray
            Node numbers of the result.

        temperature : numpy.ndarray
            Temperature at each node.

        Examples
        --------
        >>> nnum, stress = rst.nodal_temperature(0)

        """
        nnum, temp = super()._nodal_result(rnum, "EPT")
        nnum = nnum[self._mas_ind]
        temp = temp[self._mas_ind]
        if not full_rotor:  # return only the master sector
            return nnum, temp.ravel()

        # otherwise, duplicate and repeat as temperature is constant across sectors
        return nnum, temp.T.repeat(self.n_sector, axis=0)

    def plot_nodal_temperature(
        self,
        rnum,
        phase=0,
        full_rotor=True,
        show_displacement=False,
        displacement_factor=1.0,
        node_components=None,
        overlay_wireframe=False,
        add_text=True,
        element_components=None,
        sel_type_all=True,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot the nodal temperature.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        full_rotor : bool, optional
            Expand the sector solution and plot the full rotor.

        phase : float, optional
            Phase angle of the modal result in radians.  Only valid
            when full_rotor is True.  Default 0

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

        Returns
        -------
        cpos : list
            Camera position from vtk render window.

        Examples
        --------
        Plot the nodal temperature of a rotor for the first result.

        >>> result.plot_nodal_temperature(0)

        """
        # Load result from file
        if not full_rotor:
            return super().plot_nodal_temperature(
                rnum,
                show_displacement=show_displacement,
                displacement_factor=displacement_factor,
                node_components=node_components,
                element_components=element_components,
                sel_type_all=sel_type_all,
                treat_nan_as_zero=treat_nan_as_zero,
                **kwargs,
            )

        _, temp = self.nodal_temperature(rnum, True)

        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["show_displacement"] = show_displacement
        kwargs["displacement_factor"] = displacement_factor
        kwargs["overlay_wireframe"] = overlay_wireframe
        kwargs["add_text"] = add_text
        kwargs["node_components"] = node_components
        kwargs["element_components"] = element_components
        kwargs["sel_type_all"] = sel_type_all
        kwargs["phase"] = phase
        return self._plot_cyclic_point_scalars(
            temp, rnum, treat_nan_as_zero=treat_nan_as_zero, **kwargs
        )

    def animate_nodal_solution(
        self,
        rnum,
        comp="norm",
        displacement_factor=0.1,
        n_frames=180,
        add_text=True,
        loop=True,
        movie_filename=None,
        **kwargs,
    ):
        """Animate nodal solution.

        Assumes nodal solution is a displacement array from a modal solution.

        Parameters
        ----------
        rnum : int or list
            Cumulative result number with zero based indexing, or a
            list containing (step, substep) of the requested result.

        comp : str, optional
            Component of displacement to display.  Options are 'x',
            'y', 'z', or 'norm', which correspond to the x , y, z, or
            the normalized direction ``(x**2 + y**2 + z**2)**0.5``

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        n_frames : int, optional
            Number of "frames" between each full cycle.

        show_phase : bool, optional
            Shows the phase at each frame.

        add_text : bool, optional
            Includes result information at the top left-hand corner of the
            plot. Set font size with the ``font_size`` parameter.

        interpolate_before_map : bool, optional
            Leaving this at default generally results in a better plot.

        movie_filename : str, optional
            Filename of the movie to open.  Filename should end in mp4,
            but other filetypes may be supported.  See ``imagio.get_writer``.
            A single loop of the mode will be recorded.

        kwargs : optional keyword arguments, optional
            See :func:`pyvista.plot` for additional keyword arguments.

        Examples
        --------
        Generate a movie of a mode shape while plotting off-screen.

        >>> from ansys.mapdl.reader import read_binary
        >>> rst = read_binary("academic_rotor.rst")
        >>> rst.animate_nodal_displacement(
        ...     (3, 2),
        ...     displacement_factor=0.02,
        ...     movie_filename="movie.mp4",
        ...     off_screen=True
        ... )

        """
        # Avoid infinite while loop by ensure looping is disabled if off screen
        # and writing a movie
        if movie_filename and kwargs.get("off_screen", False):
            loop = False

        if "nangles" in kwargs:  # pragma: no cover
            n_frames = kwargs.pop("nangles")
            warnings.warn(
                "The ``nangles`` kwarg is depreciated and ``n_frames`` "
                "should be used instead."
            )

        rnum = self.parse_step_substep(rnum)  # need cumulative
        if "full_rotor" in kwargs:
            raise NotImplementedError("``full_rotor`` keyword argument not supported")

        # normalize nodal solution
        _, complex_disp = self.nodal_solution(rnum, as_complex=True, full_rotor=True)
        complex_disp *= displacement_factor
        complex_disp = complex_disp.reshape(-1, 3)

        if comp == "x":
            axis = 0
        elif comp == "y":
            axis = 1
        elif comp == "z":
            axis = 2
        else:
            axis = None

        result_info = ""
        if add_text:
            result_info = self.text_result_table(rnum)

        # need only the surface of the full rotor
        plot_mesh = self.full_rotor.extract_surface()
        orig_pt = plot_mesh.points.copy()

        # reduce the complex displacement to just the surface points
        ind = plot_mesh.point_data["vtkOriginalPointIds"]
        complex_disp = np.take(complex_disp, ind, axis=0)

        if axis is not None:
            scalars = complex_disp[:, axis]
        else:
            scalars = (complex_disp * complex_disp).sum(1) ** 0.5

        # initialize plotter
        font_size = kwargs.pop("font_size", 16)
        text_color = kwargs.pop("text_color", None)
        cpos = kwargs.pop("cpos", None)
        off_screen = kwargs.pop("off_screen", None)
        plotter = pv.Plotter(off_screen=off_screen)
        if kwargs.pop("show_axes", True):
            plotter.add_axes()

        if "rng" not in kwargs:
            smax = np.abs(scalars).max()
            if comp == "norm":
                kwargs["rng"] = [0, smax]
            else:
                kwargs["rng"] = [-smax, smax]

        background = kwargs.pop("background", None)
        if background:
            plotter.set_background(background)

        plotter.add_mesh(plot_mesh, scalars=np.real(scalars), **kwargs)

        # setup text
        if add_text:
            # results in a corner annotation actor
            text_actor = plotter.add_text(" ", font_size=font_size, color=text_color)

        if cpos:
            plotter.camera_position = cpos

        if movie_filename:
            movie_filename = str(movie_filename)
            if movie_filename.strip()[-3:] == "gif":
                plotter.open_gif(movie_filename)
            else:
                plotter.open_movie(movie_filename)

        self._animating = True

        def q_callback():
            """exit when user wants to leave"""
            self._animating = False

        plotter.add_key_event("q", q_callback)

        # run until q is pressed
        plotter.show(
            interactive=False, auto_close=False, interactive_update=not off_screen
        )

        first_loop = True
        while self._animating:
            for angle in np.linspace(0, np.pi * 2, n_frames):
                padj = 1 * np.cos(angle) - 1j * np.sin(angle)
                complex_disp_adj = np.real(complex_disp * padj)

                if axis is not None:
                    scalars = complex_disp_adj[:, axis]
                else:
                    scalars = (complex_disp_adj * complex_disp_adj).sum(1) ** 0.5

                # pyvista defaults to 'Data' scalars name
                plot_mesh["Data"][:] = scalars
                plot_mesh.points[:] = orig_pt + complex_disp_adj

                if add_text:
                    text_actor.set_text(
                        2,  # place in the upper left
                        f"{result_info}\nPhase {np.rad2deg(angle):.1f} Degrees",
                    )

                plotter.update(1, force_redraw=True)
                if not self._animating:
                    break

                if movie_filename and first_loop:
                    plotter.write_frame()

            first_loop = False
            if not loop:
                break

        cpos = plotter.camera_position
        plotter.close()
        return cpos

    @wraps(animate_nodal_solution)
    def animate_nodal_displacement(self, *args, **kwargs):
        """wraps animate_nodal_solution"""
        return self.animate_nodal_solution(*args, **kwargs)

    def _gen_full_rotor(self):
        """Create full rotor vtk unstructured grid"""
        grid = self._mas_grid.copy()
        # transform to standard coordinate system
        cs_cord = self._resultheader["csCord"]
        if cs_cord > 1:
            matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            grid.transform(matrix, inplace=True)

        # consider forcing low and high to be exact
        # self._mas_grid.point_data['CYCLIC_M01H'] --> rotate and match

        vtkappend = vtkAppendFilter()
        # vtkappend.MergePointsOn()
        # vtkappend.SetTolerance(1E-3)  # not available until vtk 9+
        rang = 360.0 / self.n_sector
        for i in range(self.n_sector):
            # Transform mesh
            sector = grid.copy()
            sector_id = np.empty(grid.n_points)
            sector_id[:] = i
            sector.point_data["sector_id"] = sector_id
            sector.rotate_z(rang * i, inplace=True)
            vtkappend.AddInputData(sector)

        vtkappend.Update()
        full_rotor = pv.wrap(vtkappend.GetOutput())

        if cs_cord > 1:
            matrix.Invert()
            full_rotor.transform(matrix, inplace=True)

        return full_rotor

    @property
    def full_rotor(self):
        """UnstructuredGrid of the full replicated rotor"""
        if self._rotor_cache is None:
            self._rotor_cache = self._gen_full_rotor()
        return self._rotor_cache

    def _plot_cyclic_point_scalars(
        self,
        scalars,
        rnum,
        show_displacement=False,
        displacement_factor=1,
        overlay_wireframe=False,
        add_text=True,
        node_components=None,
        element_components=None,
        sel_type_all=True,
        phase=None,
        treat_nan_as_zero=False,
        **kwargs,
    ):
        """Plot point scalars on active mesh.

        Parameters
        ----------
        scalars : numpy.ndarray
            Node scalars to plot.

        rnum : int, optional
            Cumulative result number.  Used for adding informative
            text.

        grid : pyvista.PolyData or pyvista.UnstructuredGrid, optional
            Uses self.grid by default.  When specified, uses this grid
            instead.

        show_displacement : bool, optional
            Deforms mesh according to the result.

        displacement_factor : float, optional
            Increases or decreases displacement by a factor.

        overlay_wireframe : bool, optional
            Overlay a wireframe of the original undeformed mesh.

        add_text : bool, optional
            Adds information about the result when rnum is given.

        treat_nan_as_zero : bool, optional
            Treat NAN values (i.e. stresses at midside nodes) as zero
            when plotting.

        kwargs : keyword arguments
            Additional keyword arguments.  See ``help(pyvista.plot)``

        Returns
        -------
        cpos : list
            Camera position.

        """
        # extract a sub-component of the grid if requested
        grid = self._mas_grid
        if node_components:
            grid, ind = self._extract_node_components(node_components, sel_type_all)
            if scalars is not None:
                scalars = scalars[:, ind]

        elif element_components:
            grid, ind = self._extract_element_components(element_components)
            if scalars is not None:
                scalars = scalars[:, ind]

        if treat_nan_as_zero and scalars is not None:
            scalars[np.isnan(scalars)] = 0

        # must be removed before add_mesh **kwargs
        window_size = kwargs.pop("window_size", None)
        full_screen = kwargs.pop("full_screen", False)
        screenshot = kwargs.pop("screenshot", None)
        text_color = kwargs.pop("text_color", None)
        kwargs.setdefault("cmap", "jet")
        if scalars is not None:
            kwargs.setdefault("rng", [np.nanmin(scalars), np.nanmax(scalars)])

        # Plot off screen when not interactive
        off_screen = kwargs.pop("off_screen", None)
        plotter = pv.Plotter(off_screen=off_screen)

        # various plotting properties that must be removed before add_mesh
        if kwargs.pop("show_axes", True):
            plotter.add_axes()
        plotter.background_color = kwargs.pop("background", None)
        cpos = kwargs.pop("cpos", None)

        cs_cord = self._resultheader["csCord"]
        if cs_cord > 1:
            matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix = self.cs_4x4(cs_cord, as_vtk_matrix=True)
            i_matrix.Invert()
        else:
            matrix = vtkMatrix4x4()
            i_matrix = vtkMatrix4x4()

        if overlay_wireframe:
            rang = 360.0 / self.n_sector
            for i in range(self.n_sector):
                actor = plotter.add_mesh(
                    grid.copy(False),
                    color="w",
                    style="wireframe",
                    opacity=0.5,
                    **kwargs,
                )

                # transform to standard position, rotate about Z axis,
                # transform back
                transform = vtkTransform()
                transform.RotateZ(rang * i)
                transform.Update()
                rot_matrix = transform.GetMatrix()

                if cs_cord > 1:
                    temp_matrix = vtkMatrix4x4()
                    rot_matrix.Multiply4x4(i_matrix, rot_matrix, temp_matrix)
                    rot_matrix.Multiply4x4(temp_matrix, matrix, rot_matrix)
                    transform.SetMatrix(rot_matrix)

                actor.SetUserTransform(transform)

        # add main mesh
        if show_displacement:
            _, disp = self.nodal_solution(rnum, phase, full_rotor=True)
            disp *= displacement_factor
            if node_components:
                _, ind = self._extract_node_components(node_components, sel_type_all)
                disp = disp[:, ind]
            elif element_components:
                _, ind = self._extract_element_components(element_components)
                disp = disp[:, ind]

            disp = disp.reshape(-1, 3)

            rotor = self.full_rotor.copy()
            rotor.points += disp
            actor = plotter.add_mesh(rotor, scalars=scalars.reshape(-1, 3), **kwargs)

        else:
            surf_sector = grid.extract_surface()
            ind = surf_sector.point_data["vtkOriginalPointIds"]
            rang = 360.0 / self.n_sector
            for i in range(self.n_sector):
                if scalars is not None:
                    sector_scalars = scalars[i, ind]
                else:
                    sector_scalars = None

                actor = plotter.add_mesh(
                    surf_sector.copy(False), scalars=sector_scalars, **kwargs
                )

                # NAN/missing data are white
                plotter.mapper.GetLookupTable().SetNanColor(1, 1, 1, 1)

                # transform to standard position, rotate about Z axis,
                # transform back
                transform = vtkTransform()
                transform.RotateZ(rang * i)
                transform.Update()
                rot_matrix = transform.GetMatrix()

                if cs_cord > 1:
                    temp_matrix = vtkMatrix4x4()
                    rot_matrix.Multiply4x4(i_matrix, rot_matrix, temp_matrix)
                    rot_matrix.Multiply4x4(temp_matrix, matrix, rot_matrix)
                    transform.SetMatrix(rot_matrix)

                actor.SetUserTransform(transform)

        # add table
        if isinstance(add_text, str):
            plotter.add_text(add_text, font_size=20, position=[0, 0], color=text_color)
        elif add_text:
            rnum = self.parse_step_substep(rnum)
            plotter.add_text(
                self.text_result_table(rnum), font_size=20, color=text_color
            )

        # must set camera position at the ended
        if cpos is not None:
            plotter.camera_position = cpos

        if screenshot:
            cpos = plotter.show(
                auto_close=False, window_size=window_size, full_screen=full_screen
            )
            plotter.screenshot(screenshot)
            plotter.close()
        else:
            cpos = plotter.show(window_size=window_size, full_screen=full_screen)

        return cpos

    def save_as_vtk(
        self,
        filename,
        rsets=None,
        result_types=["ENS"],
        progress_bar=True,
        expand_cyclic=True,
        merge_sectors=True,
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

        expand_cyclic : bool, default: True.
            When ``True``, expands cyclic results by writing out the result as
            a full cyclic result rather than as a single cyclic sector.

        merge_sectors : bool, default: False
            When ``expand_cyclic`` is True and this parameter is ``True``,
            sectors will be merged to create one unified grid. Set this to
            ``False`` to not merge nodes between sectors.

        Notes
        -----
        Nodal solutions are stored within the ``point_data`` attribute of the
        unstructured grid and can be accessed after reading in the result with
        pyvista with:

        .. code::

            >>> grid.point_data
            pyvista DataSetAttributes
            Association     : POINT
            Active Scalars  : Nodal stresses (0, -2)-2
            Active Vectors  : None
            Active Texture  : None
            Active Normals  : None
            Contains arrays :
                Nodal solution (0, 0)   float64    (18864, 3)
                Nodal stresses (0, 0)   float64    (18864, 6)
                Nodal solution (1, 0)   float64    (18864, 3)
                Nodal stresses (1, 0)   float64    (18864, 6)
                Nodal solution (0, -1)  float64    (18864, 3)
                Nodal stresses (0, -1)  float64    (18864, 6)
                Nodal solution (0, 1)   float64    (18864, 3)
                Nodal stresses (0, 1)   float64    (18864, 6)

        See the examples section for more details.

        Examples
        --------
        Write nodal results as a binary vtk file. Larger file size, loads quickly.

        >>> rst.save_as_vtk('results.vtk')

        Write using the xml writer. This file is more compressed compressed but
        will load slower.

        >>> rst.save_as_vtk('results.vtu')

        Write only nodal and elastic strain for the first result:

        >>> rst.save_as_vtk('results.vtk', [0], ['EEL', 'EPL'])

        Write only nodal results (i.e. displacements) for the first result:

        >>> rst.save_as_vtk('results.vtk', [0], [])

        Read in the results using ``pyvista.read()``. Plot the 'Z' component of
        the first mode's -2 nodal diameter nodal displacement.

        >>> import pyvista as pv
        >>> grid = pv.read('results.vtk')
        >>> grid.plot(scalars="Nodal solution (0, -2)", component=2)

        Do not merge sectors when saving the results and separate sectors into
        multiple blocks within pyvista.

        >>> rst.save_as_vtk('results.vtk', merge_sectors=False)
        >>> grid = pv.read('results.vtk')
        >>> mblock = grid.split_bodies()

        """

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

        # Copy grid as to not write results to original object
        if not expand_cyclic:
            super().save_as_vtk(filename, rsets, result_types, progress_bar)

        # generate the full rotor with separate sectors
        grid = self._gen_full_rotor()
        grid.cell_data.pop("vtkOriginalCellIds", None)

        ansys_node_num = None
        if "ansys_node_num" in grid.point_data:
            ansys_node_num = grid.point_data.pop("ansys_node_num")

        grid.point_data.clear()
        if ansys_node_num is not None:
            grid.point_data["ansys_node_num"] = ansys_node_num

        for i in rsets:
            # convert the result number to a harmonic index and mode number
            mode, hindex = self.mode_table[i], self.harmonic_indices[i]
            sol_name = f"({mode}, {hindex})"

            # Nodal results
            # NOTE: val is shaped (n_blades, n_nodes_sector, 3)
            _, val = self.nodal_solution(i, phase=0, full_rotor=True, as_complex=False)
            grid.point_data[f"Nodal solution {sol_name}"] = np.vstack(val)

            # Nodal results
            for rtype in self.available_results:
                if rtype in result_types:

                    def sector_func(rnum):
                        return self._nodal_result(rnum, rtype)

                    _, values = self._get_full_result(
                        i,
                        sector_func,
                        0,
                        True,
                        False,
                        tensor=rtype in RESULT_TENSORS_TYPES,
                        stress=rtype in RESULT_STRESS_TYPES,
                    )

                    desc = element_index_table_info[rtype]
                    grid.point_data[f"{desc} {sol_name}"] = np.vstack(values)

            if pbar is not None:
                pbar.update(1)

        if merge_sectors:
            grid = grid.clean(tolerance=1e-6)

        grid.save(str(filename))
        if pbar is not None:
            pbar.close()
