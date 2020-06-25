#include <assert.h>
#include <ctype.h>
#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define _POSIX_C_SOURCE 2 // To use GNU's getopt.
#define PRINT_ERROR(fmt, ...) fprintf(stderr, "\033[0;31mERROR:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...) fprintf(stderr, "\033[0;35mWARNING:\033[0m " fmt "\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...) fprintf(stderr, "\033[0;32mINFO:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_STATUS(status) fprintf(stderr, "Status: %s\n", dpu_api_status_to_string(status))
#define PRINT_DEBUG(fmt, ...) printf("\033[0;34mDEBUG:\033[0m    " fmt "\n", ##__VA_ARGS__)

#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

enum Algorithm {
  SrcVtx = 0,
  DstVtx = 1,
  Edge = 2,
};

enum Partition {
  None = 0,
  _1D = 1,
  _2D = 2,
};

struct COOMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowIdxs;
  uint32_t *colIdxs;
};

struct CSRMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowPtrs;
  uint32_t *colIdxs;
};

struct CSCMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *colPtrs;
  uint32_t *rowIdxs;
};

// Parse CLI args and options.
void parse_args(int argc, char **argv, int *num_DPUs, enum Algorithm *alg, enum Partition *prt, char **file) {
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "n:a:p:")) != -1)
    switch (c) {
    case 'n':
      *num_DPUs = atoi(optarg);
      if (num_DPUs == 0 || *num_DPUs % 8 != 0) {
        fprintf(stderr, "Number of DPUs must be a multiple of 8.\n");
        exit(1);
      }
      break;
    case 'a':
      if (strcmp(optarg, "src_vtx") == 0)
        *alg = SrcVtx;
      else if (strcmp(optarg, "dst_vtx") == 0)
        *alg = DstVtx;
      else if (strcmp(optarg, "edge") == 0)
        *alg = Edge;
      else {
        fprintf(stderr, "Incorrect -a argument. Supported algorithms: src_vtx | dst_vtx | edge\n");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "none") == 0)
        *prt = None;
      else if (strcmp(optarg, "1d") == 0)
        *prt = _1D;
      else if (strcmp(optarg, "2d") == 0)
        *prt = _2D;
      else {
        fprintf(stderr, "Incorrect -p argument. Supported partitioning: none | 1d | 2d\n");
        exit(1);
      }
      break;
    case '?':
    default:
      fprintf(stderr, "Bad args. Usage: -n <num_DPUs> -a <src_vtx|dst_vtx|edge> -p <none|1d|2d>\n");
    }

  int numargs = argc - optind;
  if (numargs != 1) {
    if (numargs > 1)
      fprintf(stderr, "Too many arguments!\n");
    else
      fprintf(stderr, "Too few arguments! Please provide data file name (0-indexed COO-formatted matrix).\n");
    exit(1);
  }
  *file = argv[optind];
}

// Reads a coo-formated file into memory.
struct COOMatrix read_coo_matrix(char *file) {

  PRINT_INFO("Loading COO-formated graph from %s.", file);

  struct COOMatrix coo;

  // Initialize fields.
  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &coo.numRows);
  fscanf(fp, "%u", &coo.numCols);
  fscanf(fp, "%u", &coo.numNonzeros);

  if (coo.numRows % 2 == 1) {
    PRINT_WARNING("Number of rows must be even. Padding with an extra row.");
    coo.numRows++;
  }
  if (coo.numCols % 2 == 1) {
    PRINT_WARNING("Number of columns must be even. Padding with an extra column.");
    coo.numCols++;
  }

  coo.rowIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));
  coo.colIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));

  PRINT_INFO("Reading COO-formated matrix - %u rows, %u columns, %u nonzeros.", coo.numRows, coo.numCols, coo.numNonzeros);

  // Read nonzeros.
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx;
    uint32_t colIdx;
    fscanf(fp, "%u", &rowIdx);
    fscanf(fp, "%u", &colIdx);
    coo.rowIdxs[i] = rowIdx;
    coo.colIdxs[i] = colIdx;
  }

  return coo;
}

// Converts COO matrix to CSR format.
struct CSRMatrix coo_to_csr(struct COOMatrix coo) {

  struct CSRMatrix csr;

  // Initialize fields
  csr.numRows = coo.numRows;
  csr.numCols = coo.numCols;
  csr.numNonzeros = coo.numNonzeros;
  csr.rowPtrs = malloc((csr.numRows + 1) * sizeof(uint32_t));
  csr.colIdxs = malloc(csr.numNonzeros * sizeof(uint32_t));

  // Histogram rowIdxs
  memset(csr.rowPtrs, 0, (csr.numRows + 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    csr.rowPtrs[rowIdx]++;
  }

  // Prefix sum rowPtrs
  uint32_t sumBeforeNextRow = 0;
  for (uint32_t rowIdx = 0; rowIdx < csr.numRows; ++rowIdx) {
    uint32_t sumBeforeRow = sumBeforeNextRow;
    sumBeforeNextRow += csr.rowPtrs[rowIdx];
    csr.rowPtrs[rowIdx] = sumBeforeRow;
  }
  csr.rowPtrs[csr.numRows] = sumBeforeNextRow;

  // Bin the nonzeros
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    uint32_t nnzIdx = csr.rowPtrs[rowIdx]++;
    csr.colIdxs[nnzIdx] = coo.colIdxs[i];
  }

  // Restore rowPtrs
  for (uint32_t rowIdx = csr.numRows - 1; rowIdx > 0; --rowIdx) {
    csr.rowPtrs[rowIdx] = csr.rowPtrs[rowIdx - 1];
  }
  csr.rowPtrs[0] = 0;

  return csr;
}

// Converts COO matrix to CSC format.
struct CSCMatrix coo_to_csc(struct COOMatrix coo) {

  // Transpose COO matrix.
  struct COOMatrix cooTrs = coo;
  cooTrs.colIdxs = coo.rowIdxs;
  cooTrs.rowIdxs = coo.colIdxs;
  cooTrs.numCols = coo.numRows;
  cooTrs.numRows = coo.numCols;

  struct CSRMatrix csrMatrix = coo_to_csr(cooTrs);
  struct CSCMatrix csc = {
      .colPtrs = csrMatrix.rowPtrs,
      .rowIdxs = csrMatrix.colIdxs,
      .numCols = csrMatrix.numCols,
      .numRows = csrMatrix.numRows};

  return csc;
}

void free_coo_matrix(struct COOMatrix coo) {
  free(coo.rowIdxs);
  free(coo.colIdxs);
}

void free_csr_matrix(struct CSRMatrix csr) {
  free(csr.rowPtrs);
  free(csr.colIdxs);
}

void free_csc_matrix(struct CSCMatrix csc) {
  free(csc.colPtrs);
  free(csc.rowIdxs);
}

/**
 * @fn dpu_insert_to_mram_uint32_array
 * @brief Inserts data into the MRAM of a DPU at the last used MRAM address.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol where to copy the pointer of the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_to_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint32_t);
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src, size, 0));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));
}

/**
 * @fn dpu_copy_to_mram_uint32_array
 * @brief Copy data to the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_copy_to_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                      // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_from_mram_uint32_array
 * @brief Copy data from the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param dst the host buffer where the data is copied.
 * @param length the number of elements in the array.
 */
void dpu_copy_from_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                  // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_to_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param src the host buffer containing the data to copy
 */
void dpu_copy_to_uint32(struct dpu_set_t dpu, const char *symbol_name, uint32_t src) {
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &src, sizeof(uint32_t)));
}

/**
 * @fn dpu_copy_from_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param dst the host buffer where the data is copied
 */
void dpu_copy_from_uint32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst) {
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, dst, sizeof(uint32_t)));
}

// Each DPU gets chunks of 32 nodes, distributed as evenly as possible.
uint32_t chunkize_nodes(uint32_t numRows, uint32_t num_DPUs, uint32_t **dpu_node_chunks) {
  uint32_t totalChunks = (numRows + 31) / 32; // ceil.
  uint32_t minChunksPerDPU = totalChunks / num_DPUs;
  uint32_t chunksRemainder = totalChunks % num_DPUs;
  uint32_t *chunksPerDPU = malloc(num_DPUs * sizeof(uint32_t));
  for (uint32_t i = 0; i < num_DPUs; ++i)
    if (i < chunksRemainder)
      chunksPerDPU[i] = minChunksPerDPU + 1;
    else
      chunksPerDPU[i] = minChunksPerDPU;

  // Node idxs for each DPU.
  uint32_t *chunks = malloc((num_DPUs + 1) * sizeof(uint32_t));
  chunks[0] = 0;
  for (uint32_t i = 1; i < num_DPUs; ++i)
    chunks[i] = chunks[i - 1] + chunksPerDPU[i - 1] * 32;
  chunks[num_DPUs] = numRows;

  *dpu_node_chunks = chunks;

  free(chunksPerDPU);
  return totalChunks;
}

void populate_mram_csr(struct dpu_set_t *set, struct dpu_set_t *dpu, struct CSRMatrix csr, uint32_t *dpu_node_chunks, uint32_t totalChunks) {

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(*set, "totalChunks", totalChunks);

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    uint32_t from = dpu_node_chunks[i];   // index of rowPtrs for this dpu.
    uint32_t to = dpu_node_chunks[i + 1]; // index of rowPtrs for the next dpu.
    uint32_t numNodes = to - from;
    uint32_t numNeighbors = csr.rowPtrs[to] - csr.rowPtrs[from];

    uint32_t nodeChunkFrom = from / 32;
    uint32_t nodeChunkTo = (to + 31) / 32;            // ceil.
    uint32_t numChunks = nodeChunkTo - nodeChunkFrom; // == chunksPerDPU[i]
    uint32_t *nodePtrs = &csr.rowPtrs[from];
    uint32_t origin = nodePtrs[0]; // index of first neighbor.
    uint32_t *neighbors = &csr.colIdxs[origin];

    PRINT_INFO("DPU %d: nodeChunks = [%u, %u[", i, nodeChunkFrom, nodeChunkTo);

    dpu_copy_to_uint32(*dpu, "numNodes", numNodes);
    dpu_copy_to_uint32(*dpu, "numNeighbors", numNeighbors);
    dpu_copy_to_uint32(*dpu, "numChunks", numChunks);
    dpu_copy_to_uint32(*dpu, "nodeChunkFrom", nodeChunkFrom);
    dpu_copy_to_uint32(*dpu, "nodeChunkTo", nodeChunkTo);
    dpu_copy_to_uint32(*dpu, "origin", origin);

    dpu_insert_to_mram_uint32_array(*dpu, "nodePtrs", nodePtrs, numNodes);
    dpu_insert_to_mram_uint32_array(*dpu, "neighbors", neighbors, numNeighbors);
    dpu_insert_to_mram_uint32_array(*dpu, "nodeLevels", 0, numNodes);

    dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, totalChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", nextFrontier, totalChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", 0, numChunks);
  }

  free(nextFrontier);
}

void populate_mram_csc(struct dpu_set_t *set, struct dpu_set_t *dpu, struct CSCMatrix csc, uint32_t *dpu_node_chunks, uint32_t totalChunks) {

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(*set, "totalChunks", totalChunks);

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    uint32_t from = dpu_node_chunks[i];   // index of rowPtrs for this dpu.
    uint32_t to = dpu_node_chunks[i + 1]; // index of rowPtrs for the next dpu.
    uint32_t numNodes = to - from;
    uint32_t numNeighbors = csc.colPtrs[to] - csc.colPtrs[from];

    uint32_t nodeChunkFrom = from / 32;
    uint32_t nodeChunkTo = (to + 31) / 32;            // ceil.
    uint32_t numChunks = nodeChunkTo - nodeChunkFrom; // == chunksPerDPU[i]
    uint32_t *nodePtrs = &csc.colPtrs[from];
    uint32_t origin = nodePtrs[0]; // index of first neighbor.
    uint32_t *neighbors = &csc.rowIdxs[origin];

    PRINT_INFO("DPU %d: nodeChunks = [%u, %u[", i, nodeChunkFrom, nodeChunkTo);

    dpu_copy_to_uint32(*dpu, "numNodes", numNodes);
    dpu_copy_to_uint32(*dpu, "numNeighbors", numNeighbors);
    dpu_copy_to_uint32(*dpu, "numChunks", numChunks);
    dpu_copy_to_uint32(*dpu, "nodeChunkFrom", nodeChunkFrom);
    dpu_copy_to_uint32(*dpu, "nodeChunkTo", nodeChunkTo);
    dpu_copy_to_uint32(*dpu, "origin", origin);

    dpu_insert_to_mram_uint32_array(*dpu, "nodePtrs", nodePtrs, numNodes);
    dpu_insert_to_mram_uint32_array(*dpu, "neighbors", neighbors, numNeighbors);
    dpu_insert_to_mram_uint32_array(*dpu, "nodeLevels", 0, numNodes);

    dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, numChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", 0, numChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", nextFrontier, totalChunks);
  }

  free(nextFrontier);
}

void bfs_src_vtx(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  uint32_t *nf_dpu = calloc(totalChunks, sizeof(uint32_t));

  uint32_t currentLevel = 0;
  uint32_t done = true;
  int i = 0;

  while (true) {
    PRINT_INFO("Level %u", currentLevel);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Get union of nextFrontiers.
    _DPU_FOREACH_I(set, dpu, i) {

      dpu_copy_from_mram_uint32_array(dpu, "nextFrontier", nf_dpu, totalChunks);
      for (uint32_t c = 0; c < totalChunks; ++c) {
        nextFrontier[c] |= nf_dpu[c];
        if (nextFrontier[c] != 0)
          done = false;
      }

      // Get DPU logs.
      // PRINT_INFO("DPU %d:", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    if (done)
      break;

    done = true;
    ++currentLevel;

    // Update currentLevel and nextFrontier of DPUs.
    dpu_copy_to_uint32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_copy_to_mram_uint32_array(dpu, "nextFrontier", nextFrontier, totalChunks);
    }

    // Clear nextFrontier.
    for (uint32_t c = 0; c < totalChunks; ++c)
      nextFrontier[c] = 0;
  }

  free(nextFrontier);
  free(nf_dpu);
}

void bfs_dst_vtx(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {

  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  uint32_t done = true;
  uint32_t currentLevel = 0;

  while (true) {
    PRINT_INFO("Level %u", currentLevel);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Concatenate all nextFrontiers.
    uint32_t writeIdx = 0;
    int i = 0;
    _DPU_FOREACH_I(set, dpu, i) {
      uint32_t numChunks;
      dpu_copy_from_uint32(dpu, "numChunks", &numChunks);
      dpu_copy_from_mram_uint32_array(dpu, "nextFrontier",
                                      &nextFrontier[writeIdx], numChunks);
      writeIdx += numChunks;

      // Get DPU logs.
      // PRINT_INFO("DPU %d", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    // Check if done.
    for (uint32_t c = 0; c < totalChunks; ++c)
      if (nextFrontier[c] != 0) {
        done = false;
        break;
      }
    if (done)
      break;

    done = true;
    ++currentLevel;

    // Update currentLevel and currentFrontier of DPUs.
    dpu_copy_to_uint32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_copy_to_mram_uint32_array(dpu, "currFrontier", nextFrontier, totalChunks);
    }
  }

  // Free resources.
  free(nextFrontier);
}

void bfs_edge(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {
}

void print_node_levels(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t numNodes) {
  // Get nodeLevels from each DPU.
  uint32_t *nodeLevels = malloc(numNodes * sizeof(uint32_t));
  uint32_t writeIdx = 0;

  int i = 0;
  _DPU_FOREACH_I(set, dpu, i) {
    uint32_t numNodes;
    dpu_copy_from_uint32(dpu, "numNodes", &numNodes);
    dpu_copy_from_mram_uint32_array(dpu, "nodeLevels", &nodeLevels[writeIdx],
                                    numNodes);
    writeIdx += numNodes;
  }

  // Output.
  PRINT_INFO("Output:");
  for (uint32_t node = 0; node < numNodes; ++node)
    printf("nodeLevels[%u]=%u\n", node, nodeLevels[node]);

  free(nodeLevels);
}

int main(int argc, char **argv) {

  int num_DPUs = 8;
  enum Algorithm alg = SrcVtx;
  enum Partition prt = None;
  char *file = NULL;

  parse_args(argc, argv, &num_DPUs, &alg, &prt, &file);

  char *bin_path;
  switch (alg) {
  case SrcVtx:
    bin_path = "bin/src-vtx";
    break;
  case DstVtx:
    bin_path = "bin/dst-vtx";
    break;
  case Edge:
    bin_path = "bin/edge";
    break;
  }

  struct dpu_set_t set, dpu;
  PRINT_INFO("Allocating %d DPUs.", num_DPUs);
  DPU_ASSERT(dpu_alloc(num_DPUs, NULL, &set));
  DPU_ASSERT(dpu_load(set, bin_path, NULL));

  struct COOMatrix coo = read_coo_matrix(file); // Load coo-matrix from file.
  uint32_t *dpu_node_chunks;                    // Node-chunk indexes for each DPU.

  switch (alg) {
  case SrcVtx: {
    switch (prt) {
    case None: {
      struct CSRMatrix csr = coo_to_csr(coo);
      uint32_t totalChunks = chunkize_nodes(csr.numRows, num_DPUs, &dpu_node_chunks);
      populate_mram_csr(&set, &dpu, csr, dpu_node_chunks, totalChunks);
      bfs_src_vtx(set, dpu, totalChunks);
      free_csr_matrix(csr);
    } break;
    case _1D: {
    } break;

    case _2D: {

    } break;
    }

  } break;
  case DstVtx: {
    // struct CSCMatrix csc = coo2csc(coo);
    // uint32_t totalChunks = chunkize_nodes(csc.numCols, num_DPUs, &dpu_node_chunks);
    // populate_mram_csc(&set, &dpu, csc, dpu_node_chunks, totalChunks);
    // bfs_dst_vtx(set, dpu, totalChunks);
    // free_csc_matrix(csc);
  } break;
  case Edge:
    break;
  }

  // Distribute nodes to DPUs as fairly as possible.
  print_node_levels(set, dpu, coo.numRows);

  // Free resources.
  free(dpu_node_chunks);
  free_coo_matrix(coo);
  DPU_ASSERT(dpu_free(set));
}
