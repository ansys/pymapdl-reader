# cython: language_level=3
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: embedsignature=True

from libc.stdio cimport fopen, fclose, FILE

ctypedef unsigned char uint8_t
from libc.stdint cimport int64_t


cdef extern from 'archive.h' nogil:
    int write_nblock(FILE*, const int, const int, const int*, const double*,
                     const double*, int)
    int write_eblock(FILE*, const int, const int*, const int*, const int*,
                     const int*, const int*, const uint8_t*, const int64_t*,
                     const int64_t*, const int*, const int*, const int);


def py_write_nblock(filename, const int [::1] node_id, int max_node_id,
                    const double [:, ::1] pos, const double [:, ::1] angles,
                    mode='w'):
    """Write a node block to a file.

    Parameters
    ----------
    fid : _io.TextIOWrapper
        Opened Python file object.

    node_id : np.ndarray
        Array of node ids.

    pos : np.ndarray
        Double array of node coordinates

    angles : np.ndarray, optional
        

    """
    # attach the stream to the python file
    cdef FILE* cfile = fopen(filename.encode(), mode.encode())

    cdef int n_nodes = pos.shape[0]
    cdef double [::1] dummy_arr

    cdef int has_angles = 0
    if angles.size == pos.size:
        has_angles = 1
    write_nblock(cfile, n_nodes, max_node_id, &node_id[0], &pos[0, 0],
                 &angles[0, 0], has_angles);
    fclose(cfile)


def py_write_eblock(filename,
                    int [::1] elem_id,
                    int [::1] etype,
                    int [::1] mtype,
                    int [::1] rcon,
                    int [::1] elem_nnodes,
                    int64_t [::1] cells,
                    int64_t [::1] offset,
                    uint8_t [::1] celltypes,
                    int [::1] typenum,
                    int [::1] nodenum,
                    vtk9,
                    mode='w'):
    cdef FILE* cfile = fopen(filename.encode(), mode.encode())
    write_eblock(cfile,
                 celltypes.size,
                 &elem_id[0],
                 &etype[0],
                 &mtype[0],
                 &rcon[0],
                 &elem_nnodes[0],
                 &celltypes[0],
                 &offset[0],
                 &cells[0],
                 &typenum[0],
                 &nodenum[0],
                 int(vtk9))
    fclose(cfile)
