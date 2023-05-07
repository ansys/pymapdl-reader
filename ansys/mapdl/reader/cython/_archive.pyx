# cython: language_level=3
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: embedsignature=True

from libc.math cimport abs
from libc.stdio cimport FILE, fclose, fdopen, fopen, fwrite

import numpy as np

cimport numpy as np

ctypedef unsigned char uint8_t
from libc.stdint cimport int64_t


cdef extern from 'archive.h' nogil:
    int write_nblock(FILE*, const int, const int, const int*, const double*,
                     const double*, int)
    int write_nblock_float(FILE*, const int, const int, const int*, const float*,
                           const float*, int)
    int write_eblock(FILE*, const int, const int*, const int*, const int*,
                     const int*, const int*, const uint8_t*, const int64_t*,
                     const int64_t*, const int*, const int*);


cdef extern from "stdio.h":
    FILE *fdopen(int, const char *)



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

    cdef int has_angles = 0
    if angles.size == pos.size:
        has_angles = 1
    else:
        angles = np.zeros((1, 1), np.double)
    write_nblock(cfile, n_nodes, max_node_id, &node_id[0], &pos[0, 0],
                 &angles[0, 0], has_angles);
    fclose(cfile)


def py_write_nblock_float(filename, const int [::1] node_id, int max_node_id,
                          const float [:, ::1] pos, const float [:, ::1] angles,
                          mode='w'):
    """Write a float32 node block to a file.

    Parameters
    ----------
    fid : _io.TextIOWrapper
        Opened Python file object.

    node_id : np.ndarray
        Array of node ids.

    pos : np.float32 np.ndarray
        Double array of node coordinates

    angles : np.ndarray, optional
    
    """
    # attach the stream to the python file
    cdef FILE* cfile = fopen(filename.encode(), mode.encode())

    cdef int n_nodes = pos.shape[0]
    cdef int has_angles = 0
    if angles.size == pos.size:
        has_angles = 1
    else:
        angles = np.zeros((1, 1), np.float32)
    write_nblock_float(cfile, n_nodes, max_node_id, &node_id[0], &pos[0, 0],
                       &angles[0, 0], has_angles);
    fclose(cfile)


def py_write_eblock(filename,
                    const int [::1] elem_id,
                    const int [::1] etype,
                    const int [::1] mtype,
                    const int [::1] rcon,
                    const int [::1] elem_nnodes,
                    const int64_t [::1] cells,
                    const int64_t [::1] offset,
                    const uint8_t [::1] celltypes,
                    const int [::1] typenum,
                    const int [::1] nodenum,
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
                 &nodenum[0])
    fclose(cfile)


def cmblock_items_from_array(int [::1] array):
    """Given a list of items, convert to a ANSYS formatted CMBLOCK.

    For example ``1, 2, 3, 4, 8``

    will be converted to

    ``1, -4, 8``

    Where the -4 indicates all the items between 1 and -4.
    """

    # first, verify items in array are sorted
    cdef int i
    cdef is_sorted = 1
    for i in range(array.size - 1):
        if array[i] > array[i + 1]:
            is_sorted = 0
            break

    if not is_sorted:
        array = np.unique(array)

    cdef int [::1] items = np.empty_like(array)
    items[0] = array[0]

    cdef int c = 1
    cdef int in_list = 0
    for i in range(array.size - 1):
        # check if part of a range
        if array[i + 1] - array[i] == 1:
            in_list = 1
        elif array[i + 1] - array[i] > 1:
            if in_list:
                items[c] = -array[i]; c += 1
                items[c] = array[i + 1]; c += 1
            else:
                items[c] = array[i + 1]; c += 1
            in_list = 0

    # catch if last item is part of a list
    if items[c - 1] != abs(array[array.size - 1]):
        items[c] = -array[array.size - 1]; c += 1

    return np.array(items[:c])
