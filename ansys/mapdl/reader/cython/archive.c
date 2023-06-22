#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

// necessary for ubuntu build on azure
#ifdef __linux__
#include <stdint.h>
#endif

// VTK cell types
#define VTK_TRIANGLE 5
#define VTK_QUAD 9
#define VTK_TETRA 10
#define VTK_VOXEL 11
#define VTK_HEXAHEDRON 12
#define VTK_WEDGE 13
#define VTK_PYRAMID 14
#define VTK_QUADRATIC_TETRA 24
#define VTK_QUADRATIC_HEXAHEDRON 25
#define VTK_QUADRATIC_WEDGE 26
#define VTK_QUADRATIC_PYRAMID 27

// Write node IDs, node coordinates, and angles to file as a NBLOCK
int write_nblock(FILE *file, const int n_nodes, const int max_node_id,
                 const int *node_id, const double *nodes, const double *angles,
                 int has_angles) {
  // Header
  // Tell ANSYS to start reading the node block with 6 fields,
  // associated with a solid, the maximum node number and the number
  // of lines in the node block
  fprintf(file, "/PREP7\n");
  fprintf(file, "NBLOCK,6,SOLID,%10d,%10d\n", max_node_id, n_nodes);
  fprintf(file, "(3i8,6e20.13)\n");

  int i;
  if (has_angles) {
    for (i = 0; i < n_nodes; i++) {
      fprintf(file,
              "%8d       0       0%20.12E%20.12E%20.12E%20.12E%20.12E%20.12E\n",
              node_id[i], nodes[i * 3 + 0], nodes[i * 3 + 1], nodes[i * 3 + 2],
              angles[i * 3 + 0], angles[i * 3 + 1], angles[i * 3 + 2]);
    }
  } else {
    for (i = 0; i < n_nodes; i++) {
      fprintf(file, "%8d       0       0%20.12E%20.12E%20.12E\n", node_id[i],
              nodes[i * 3 + 0], nodes[i * 3 + 1], nodes[i * 3 + 2]);
    }
  }

  fprintf(file, "N,R5.3,LOC,       -1,\n");
  return 0;
}

int write_nblock_float(FILE *file, const int n_nodes, const int max_node_id,
                       const int *node_id, const float *nodes,
                       const float *angles, int has_angles) {
  // Header
  // Tell ANSYS to start reading the node block with 6 fields,
  // associated with a solid, the maximum node number and the number
  // of lines in the node block
  fprintf(file, "/PREP7\n");
  fprintf(file, "NBLOCK,6,SOLID,%10d,%10d\n", max_node_id, n_nodes);
  fprintf(file, "(3i8,6e20.13)\n");

  int i;
  if (has_angles) {
    for (i = 0; i < n_nodes; i++) {
      fprintf(file,
              "%8d       0       0%20.12E%20.12E%20.12E%20.12E%20.12E%20.12E\n",
              node_id[i], nodes[i * 3 + 0], nodes[i * 3 + 1], nodes[i * 3 + 2],
              angles[i * 3 + 0], angles[i * 3 + 1], angles[i * 3 + 2]);
    }
  } else {
    for (i = 0; i < n_nodes; i++) {
      fprintf(file, "%8d       0       0%20.12E%20.12E%20.12E\n", node_id[i],
              nodes[i * 3 + 0], nodes[i * 3 + 1], nodes[i * 3 + 2]);
    }
  }

  fprintf(file, "N,R5.3,LOC,       -1,\n");
  return 0;
}

// Write an ANSYS EBLOCK to file
int write_eblock(
    FILE *file,
    const int n_elem,         // number of elements
    const int *elem_id,       // element ID array
    const int *etype,         // element type ID array
    const int *mtype,         // material ID array
    const int *rcon,          // real constant ID array
    const int *elem_nnodes,   // number of nodes per element
    const uint8_t *celltypes, // VTK celltypes array
    const int *offset,        // VTK offset array
    const int *cells,         // VTK cell connectivity array
    const int *typenum,       // ANSYS type number (e.g. 187 for SOLID187)
    const int *nodenum) {     // ANSYS node numbering

  // Write header
  fprintf(file, "EBLOCK,19,SOLID,%10d,%10d\n", elem_id[n_elem - 1], n_elem);
  fprintf(file, "(19i8)\n");

  int c; // position within offset array
  for (int i = 0; i < n_elem; i++) {
    c = offset[i];

    // Write cell info
    fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d",
            mtype[i],       // Field 1: material reference number
            etype[i],       // Field 2: element type number
            rcon[i],        // Field 3: real constant reference number
            1,              // Field 4: section number
            0,              // Field 5: element coordinate system
            0,              // Field 6: Birth/death flag
            0,              // Field 7:
            0,              // Field 8:
            elem_nnodes[i], // Field 9: Number of nodes
            0,              // Field 10: Not Used
            elem_id[i]);    // Field 11: Element number

    // Write element nodes
    switch (celltypes[i]) {
    case VTK_QUADRATIC_TETRA:
      if (typenum[i] == 187) {
        fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n%8d%8d\n",
                nodenum[cells[c + 0]],  // 0, I
                nodenum[cells[c + 1]],  // 1, J
                nodenum[cells[c + 2]],  // 2, K
                nodenum[cells[c + 3]],  // 3, L
                nodenum[cells[c + 4]],  // 4, M
                nodenum[cells[c + 5]],  // 5, N
                nodenum[cells[c + 6]],  // 6, O
                nodenum[cells[c + 7]],  // 7, P
                nodenum[cells[c + 8]],  // 8, Q
                nodenum[cells[c + 9]]); // 9, R
      } else {
        // Using SOLID186-like format
        fprintf(
            file,
            "%8d%8d%8d%8d%8d%8d%8d%8d\n%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
            nodenum[cells[c + 0]],  // 0,  I
            nodenum[cells[c + 1]],  // 1,  J
            nodenum[cells[c + 2]],  // 2,  K
            nodenum[cells[c + 2]],  // 3,  L (duplicate of K)
            nodenum[cells[c + 3]],  // 4,  M
            nodenum[cells[c + 3]],  // 5,  N (duplicate of M)
            nodenum[cells[c + 3]],  // 6,  O (duplicate of M)
            nodenum[cells[c + 3]],  // 7,  P (duplicate of M)
            nodenum[cells[c + 4]],  // 8,  Q
            nodenum[cells[c + 5]],  // 9,  R
            nodenum[cells[c + 3]],  // 10, S (duplicate of K)
            nodenum[cells[c + 6]],  // 11, T
            nodenum[cells[c + 3]],  // 12, U (duplicate of M)
            nodenum[cells[c + 3]],  // 13, V (duplicate of M)
            nodenum[cells[c + 3]],  // 14, W (duplicate of M)
            nodenum[cells[c + 3]],  // 15, X (duplicate of M)
            nodenum[cells[c + 7]],  // 16, Y
            nodenum[cells[c + 8]],  // 17, Z
            nodenum[cells[c + 9]],  // 18, A
            nodenum[cells[c + 9]]); // 19, B (duplicate of A)
      }
      break;
    case VTK_TETRA: // point
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",
              nodenum[cells[c + 0]],  // 0,  I
              nodenum[cells[c + 1]],  // 1,  J
              nodenum[cells[c + 2]],  // 2,  K
              nodenum[cells[c + 2]],  // 3,  L (duplicate of K)
              nodenum[cells[c + 3]],  // 4,  M
              nodenum[cells[c + 3]],  // 5,  N (duplicate of M)
              nodenum[cells[c + 3]],  // 6,  O (duplicate of M)
              nodenum[cells[c + 3]]); // 7,  P (duplicate of M)
      break;
    case VTK_WEDGE:
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",
              nodenum[cells[c + 2]],  // 0,  I
              nodenum[cells[c + 1]],  // 1,  J
              nodenum[cells[c + 0]],  // 2,  K
              nodenum[cells[c + 0]],  // 3,  L (duplicate of K)
              nodenum[cells[c + 5]],  // 4,  M
              nodenum[cells[c + 4]],  // 5,  N
              nodenum[cells[c + 3]],  // 6,  O
              nodenum[cells[c + 3]]); // 7,  P (duplicate of O)
      break;
    case VTK_QUADRATIC_WEDGE:
      fprintf(
          file,
          "%8d%8d%8d%8d%8d%8d%8d%8d\n%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          nodenum[cells[c + 2]],   // 0,  I
          nodenum[cells[c + 1]],   // 1,  J
          nodenum[cells[c + 0]],   // 2,  K
          nodenum[cells[c + 0]],   // 3,  L (duplicate of K)
          nodenum[cells[c + 5]],   // 4,  M
          nodenum[cells[c + 4]],   // 5,  N
          nodenum[cells[c + 3]],   // 6,  O
          nodenum[cells[c + 3]],   // 7,  P (duplicate of O)
          nodenum[cells[c + 7]],   // 8,  Q
          nodenum[cells[c + 6]],   // 9,  R
          nodenum[cells[c + 0]],   // 10, S   (duplicate of K)
          nodenum[cells[c + 8]],   // 11, T
          nodenum[cells[c + 10]],  // 12, U
          nodenum[cells[c + 9]],   // 13, V
          nodenum[cells[c + 3]],   // 14, W (duplicate of O)
          nodenum[cells[c + 11]],  // 15, X
          nodenum[cells[c + 14]],  // 16, Y
          nodenum[cells[c + 13]],  // 17, Z
          nodenum[cells[c + 12]],  // 18, A
          nodenum[cells[c + 12]]); // 19, B (duplicate of A)
      break;
    case VTK_QUADRATIC_PYRAMID:
      fprintf(
          file,
          "%8d%8d%8d%8d%8d%8d%8d%8d\n%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          nodenum[cells[c + 0]],   // 0,  I
          nodenum[cells[c + 1]],   // 1,  J
          nodenum[cells[c + 2]],   // 2,  K
          nodenum[cells[c + 3]],   // 3,  L
          nodenum[cells[c + 4]],   // 4,  M
          nodenum[cells[c + 4]],   // 5,  N (duplicate of M)
          nodenum[cells[c + 4]],   // 6,  O (duplicate of M)
          nodenum[cells[c + 4]],   // 7,  P (duplicate of M)
          nodenum[cells[c + 5]],   // 8,  Q
          nodenum[cells[c + 6]],   // 9,  R
          nodenum[cells[c + 7]],   // 10, S
          nodenum[cells[c + 8]],   // 11, T
          nodenum[cells[c + 4]],   // 12, U (duplicate of M)
          nodenum[cells[c + 4]],   // 13, V (duplicate of M)
          nodenum[cells[c + 4]],   // 14, W (duplicate of M)
          nodenum[cells[c + 4]],   // 15, X (duplicate of M)
          nodenum[cells[c + 9]],   // 16, Y
          nodenum[cells[c + 10]],  // 17, Z
          nodenum[cells[c + 11]],  // 18, A
          nodenum[cells[c + 12]]); // 19, B (duplicate of A)
      break;
    case VTK_PYRAMID:
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",
              nodenum[cells[c + 0]],  // 0,  I
              nodenum[cells[c + 1]],  // 1,  J
              nodenum[cells[c + 2]],  // 2,  K
              nodenum[cells[c + 3]],  // 3,  L
              nodenum[cells[c + 4]],  // 4,  M
              nodenum[cells[c + 4]],  // 5,  N (duplicate of M)
              nodenum[cells[c + 4]],  // 6,  O (duplicate of M)
              nodenum[cells[c + 4]]); // 7,  P (duplicate of M)
      break;
    case VTK_VOXEL:
      // note the flipped order for nodes (K, L) and (O, P)
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",
              nodenum[cells[c + 0]],  // 0, I
              nodenum[cells[c + 1]],  // 1, J
              nodenum[cells[c + 3]],  // 2, K
              nodenum[cells[c + 2]],  // 3, L
              nodenum[cells[c + 4]],  // 4, M
              nodenum[cells[c + 5]],  // 5, N
              nodenum[cells[c + 7]],  // 6, O
              nodenum[cells[c + 6]]); // 7, P
      break;
    case VTK_HEXAHEDRON:
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n", nodenum[cells[c + 0]],
              nodenum[cells[c + 1]], nodenum[cells[c + 2]],
              nodenum[cells[c + 3]], nodenum[cells[c + 4]],
              nodenum[cells[c + 5]], nodenum[cells[c + 6]],
              nodenum[cells[c + 7]]);
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      fprintf(
          file,
          "%8d%8d%8d%8d%8d%8d%8d%8d\n%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          nodenum[cells[c + 0]], nodenum[cells[c + 1]], nodenum[cells[c + 2]],
          nodenum[cells[c + 3]], nodenum[cells[c + 4]], nodenum[cells[c + 5]],
          nodenum[cells[c + 6]], nodenum[cells[c + 7]], nodenum[cells[c + 8]],
          nodenum[cells[c + 9]], nodenum[cells[c + 10]], nodenum[cells[c + 11]],
          nodenum[cells[c + 12]], nodenum[cells[c + 13]],
          nodenum[cells[c + 14]], nodenum[cells[c + 15]],
          nodenum[cells[c + 16]], nodenum[cells[c + 17]],
          nodenum[cells[c + 18]], nodenum[cells[c + 19]]);
      break;
    case VTK_TRIANGLE:
      fprintf(file, "%8d%8d%8d%8d\n",
              nodenum[cells[c + 0]],  // 0,  I
              nodenum[cells[c + 1]],  // 1,  J
              nodenum[cells[c + 2]],  // 2,  K
              nodenum[cells[c + 2]]); // 3,  L (duplicate of K)
      break;
    case VTK_QUAD:
      fprintf(file, "%8d%8d%8d%8d\n",
              nodenum[cells[c + 0]],  // 0,  I
              nodenum[cells[c + 1]],  // 1,  J
              nodenum[cells[c + 2]],  // 2,  K
              nodenum[cells[c + 3]]); // 3,  L
      break;
    }
  }
  fprintf(file, "      -1\n");
  return 0;
}
