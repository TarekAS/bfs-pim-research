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

enum algorithm {
  src_vtx = 0,
  dst_vtx = 1,
  edge = 2,
};

enum partition {
  none = 0,
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

struct COOMatrix readCOOMatrix(char *fileName) {

  PRINT_INFO("Loading graph from file %s", fileName);

  struct COOMatrix cooMatrix;

  // Initialize fields
  FILE *fp = fopen(fileName, "r");
  fscanf(fp, "%u", &cooMatrix.numRows);

  if (cooMatrix.numRows % 2 == 1) {
    PRINT_WARNING("Reading matrix %s: number of rows must be even. Padding with an extra row.", fileName);
    cooMatrix.numRows++;
  }

  fscanf(fp, "%u", &cooMatrix.numCols);
  fscanf(fp, "%u", &cooMatrix.numNonzeros);
  cooMatrix.rowIdxs = malloc(cooMatrix.numNonzeros * sizeof(uint32_t));
  cooMatrix.colIdxs = malloc(cooMatrix.numNonzeros * sizeof(uint32_t));

  PRINT_INFO("Reading matrix %s: %u rows, %u columns, %u nonzeros", fileName, cooMatrix.numRows, cooMatrix.numCols, cooMatrix.numNonzeros);

  // Read the nonzeros
  for (uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
    uint32_t rowIdx;
    fscanf(fp, "%u", &rowIdx);
    cooMatrix.rowIdxs[i] = rowIdx;
    uint32_t colIdx;
    fscanf(fp, "%u", &colIdx);
    cooMatrix.colIdxs[i] = colIdx;
  }

  return cooMatrix;
}

struct CSRMatrix coo2csr(struct COOMatrix cooMatrix) {

  struct CSRMatrix csrMatrix;

  // Initialize fields
  csrMatrix.numRows = cooMatrix.numRows;
  csrMatrix.numCols = cooMatrix.numCols;
  csrMatrix.numNonzeros = cooMatrix.numNonzeros;
  csrMatrix.rowPtrs = malloc((csrMatrix.numRows + 1) * sizeof(uint32_t));
  csrMatrix.colIdxs = malloc(csrMatrix.numNonzeros * sizeof(uint32_t));

  // Histogram rowIdxs
  memset(csrMatrix.rowPtrs, 0, (csrMatrix.numRows + 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
    uint32_t rowIdx = cooMatrix.rowIdxs[i];
    csrMatrix.rowPtrs[rowIdx]++;
  }

  // Prefix sum rowPtrs
  uint32_t sumBeforeNextRow = 0;
  for (uint32_t rowIdx = 0; rowIdx < csrMatrix.numRows; ++rowIdx) {
    uint32_t sumBeforeRow = sumBeforeNextRow;
    sumBeforeNextRow += csrMatrix.rowPtrs[rowIdx];
    csrMatrix.rowPtrs[rowIdx] = sumBeforeRow;
  }
  csrMatrix.rowPtrs[csrMatrix.numRows] = sumBeforeNextRow;

  // Bin the nonzeros
  for (uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
    uint32_t rowIdx = cooMatrix.rowIdxs[i];
    uint32_t nnzIdx = csrMatrix.rowPtrs[rowIdx]++;
    csrMatrix.colIdxs[nnzIdx] = cooMatrix.colIdxs[i];
  }

  // Restore rowPtrs
  for (uint32_t rowIdx = csrMatrix.numRows - 1; rowIdx > 0; --rowIdx) {
    csrMatrix.rowPtrs[rowIdx] = csrMatrix.rowPtrs[rowIdx - 1];
  }
  csrMatrix.rowPtrs[0] = 0;

  return csrMatrix;
}

void freeCOOMatrix(struct COOMatrix cooMatrix) {
  free(cooMatrix.rowIdxs);
  free(cooMatrix.colIdxs);
}

void freeCSRMatrix(struct CSRMatrix csrMatrix) {
  free(csrMatrix.rowPtrs);
  free(csrMatrix.colIdxs);
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

// Parse CLI args and options.
void parse_args(int argc, char **argv, int *num_DPUs, enum algorithm *alg, enum partition *prt, char **filename) {
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
        *alg = src_vtx;
      else if (strcmp(optarg, "dst_vtx") == 0)
        *alg = dst_vtx;
      else if (strcmp(optarg, "edge") == 0)
        *alg = edge;
      else {
        fprintf(stderr, "Incorrect -a argument. Supported algorithms: src_vtx | dst_vtx | edge\n");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "none") == 0)
        *prt = none;
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
  *filename = argv[optind];
}

void allocate_dpus(struct dpu_set_t *set, struct dpu_set_t *dpu, int num_DPUs, enum algorithm alg) {
  PRINT_INFO("Allocating %d DPUs.", num_DPUs);
  DPU_ASSERT(dpu_alloc(num_DPUs, NULL, set));
  char *binPath;
  switch (alg) {
  case src_vtx:
    binPath = "bin/src-vtx";
    break;
  case dst_vtx:
    binPath = "bin/dst-vtx";
    break;
  case edge:
    binPath = "bin/edge";
    break;
  }
  DPU_ASSERT(dpu_load(*set, binPath, NULL));
}

void populate_mram(struct dpu_set_t *set, struct dpu_set_t *dpu, struct CSRMatrix csr, uint32_t *dpu_node_chunks, uint32_t totalChunks, enum algorithm alg) {

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(*set, "totalChunks", totalChunks);

  uint32_t i;
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

    switch (alg) {
    case src_vtx:
      dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, totalChunks);
      dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", nextFrontier, totalChunks);
      dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", 0, numChunks);
      break;
    case dst_vtx:
      dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, numChunks);
      dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", 0, numChunks);
      dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", nextFrontier, totalChunks);
      break;
    case edge:
      break;
    }
  }

  free(nextFrontier);
}

// Each DPU gets chunks of 32 nodes, distributed as evenly as possible.
void chunkize_nodes(struct CSRMatrix csr, uint32_t num_DPUs, uint32_t **dpu_node_chunks, uint32_t *totalChunks) {
  *totalChunks = (csr.numRows + 31) / 32; // ceil.
  uint32_t minChunksPerDPU = *totalChunks / num_DPUs;
  uint32_t chunksRemainder = *totalChunks % num_DPUs;
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
  chunks[num_DPUs] = csr.numRows;

  *dpu_node_chunks = chunks;

  free(chunksPerDPU);
}

void print_node_levels(struct dpu_set_t set, struct dpu_set_t dpu, struct CSRMatrix csr) {
  // Get nodeLevels from each DPU.
  uint32_t *nodeLevels = malloc(csr.numRows * sizeof(uint32_t));
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
  for (uint32_t node = 0; node < csr.numRows; ++node)
    printf("nodeLevels[%u]=%u\n", node, nodeLevels[node]);

  free(nodeLevels);
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

int main(int argc, char **argv) {

  int num_DPUs = 8;
  enum algorithm alg = src_vtx;
  enum partition prt = none;
  char *filename = NULL;

  parse_args(argc, argv, &num_DPUs, &alg, &prt, &filename);

  struct COOMatrix coo = readCOOMatrix(filename); // Load coo-matrix from file.
  struct CSRMatrix csr = coo2csr(coo);

  struct dpu_set_t set, dpu;
  allocate_dpus(&set, &dpu, num_DPUs, alg);

  // Distribute nodes to DPUs as fairly as possible.
  uint32_t *dpu_node_chunks; // Node-chunk indexes for each DPU.
  uint32_t totalChunks;
  chunkize_nodes(csr, num_DPUs, &dpu_node_chunks, &totalChunks);
  populate_mram(&set, &dpu, csr, dpu_node_chunks, totalChunks, alg);

  switch (alg) {
  case src_vtx:
    bfs_src_vtx(set, dpu, totalChunks);
    break;
  case dst_vtx:
    bfs_dst_vtx(set, dpu, totalChunks);
    break;
  case edge:
    bfs_edge(set, dpu, totalChunks);
    break;
  }

  print_node_levels(set, dpu, csr);

  // Free resources.
  PRINT_INFO("Freeing resources");
  free(dpu_node_chunks);
  freeCOOMatrix(coo);
  freeCSRMatrix(csr);
  DPU_ASSERT(dpu_free(set));
}
