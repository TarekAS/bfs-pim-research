#include <assert.h>
#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef DPU_BINARY
#define DPU_BINARY "./bin/task"
#endif

#define PRINT_ERROR(fmt, ...)                                                  \
  fprintf(stderr, "\033[0;31mERROR:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...)                                                \
  fprintf(stderr, "\033[0;35mWARNING:\033[0m " fmt "\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...)                                                   \
  fprintf(stderr, "\033[0;32mINFO:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_STATUS(status)                                                   \
  fprintf(stderr, "Status: %s\n", dpu_api_status_to_string(status))
#define PRINT_DEBUG(fmt, ...)                                                  \
  printf("\033[0;34mDEBUG:\033[0m    " fmt "\n", ##__VA_ARGS__)

#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

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

struct COOMatrix readCOOMatrix(const char *fileName) {

  struct COOMatrix cooMatrix;

  // Initialize fields
  FILE *fp = fopen(fileName, "r");
  fscanf(fp, "%u", &cooMatrix.numRows);
  if (cooMatrix.numRows % 2 == 1) {
    PRINT_WARNING("Reading matrix %s: number of rows must be even. Padding "
                  "with an extra row.",
                  fileName);
    cooMatrix.numRows++;
  }
  fscanf(fp, "%u", &cooMatrix.numCols);
  fscanf(fp, "%u", &cooMatrix.numNonzeros);
  cooMatrix.rowIdxs = malloc(cooMatrix.numNonzeros * sizeof(uint32_t));
  cooMatrix.colIdxs = malloc(cooMatrix.numNonzeros * sizeof(uint32_t));

  PRINT_INFO("Reading matrix %s: %u rows, %u columns, %u nonzeros", fileName,
             cooMatrix.numRows, cooMatrix.numCols, cooMatrix.numNonzeros);

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
 * @param symbol_name the name of the DPU symbol where to copy the pointer of
 * the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_to_mram_uint32_array(struct dpu_set_t dpu,
                                     const char *symbol_name, uint32_t *src,
                                     uint32_t length) {
  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end,
                           sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(
      dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint32_t);
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src,
                              size, 0));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end,
                         sizeof(mram_addr_t)));
}

/**
 * @fn dpu_copy_to_mram_uint32_array
 * @brief Copy data to the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM
 * destination.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_copy_to_mram_uint32_array(struct dpu_set_t dpu,
                                   const char *symbol_name, uint32_t *src,
                                   uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(
      dpu_copy_from(dpu, symbol_name, 0, &p_array,
                    sizeof(mram_addr_t))); // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src,
                              length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_from_mram_uint32_array
 * @brief Copy data from the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM
 * destination.
 * @param dst the host buffer where the data is copied.
 * @param length the number of elements in the array.
 */
void dpu_copy_from_mram_uint32_array(struct dpu_set_t dpu,
                                     const char *symbol_name, uint32_t *dst,
                                     uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(
      dpu_copy_from(dpu, symbol_name, 0, &p_array,
                    sizeof(mram_addr_t))); // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array,
                                length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_to_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param src the host buffer containing the data to copy
 */
void dpu_copy_to_uint32(struct dpu_set_t dpu, const char *symbol_name,
                        uint32_t src) {
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &src, sizeof(uint32_t)));
}

/**
 * @fn dpu_copy_from_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param dst the host buffer where the data is copied
 */
void dpu_copy_from_uint32(struct dpu_set_t dpu, const char *symbol_name,
                          uint32_t *dst) {
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, dst, sizeof(uint32_t)));
}

void bfs() {}

int main() {

  // Load coo-matrix from file and convert to csr.
  struct COOMatrix coo = readCOOMatrix("data/loc-gowalla_edges.txt");
  struct CSRMatrix csr = coo2csr(coo);

  // Each DPU gets chunks of 32 nodes, distributed as evenly as possible.
  uint32_t numDPUs = 8;
  uint32_t totalChunks = (csr.numRows + 31) / 32; // ceil.
  uint32_t minChunksPerDPU = totalChunks / numDPUs;
  uint32_t chunksRemainder = totalChunks % numDPUs;
  uint32_t *chunksPerDPU = malloc(numDPUs * sizeof(uint32_t));
  for (uint32_t i = 0; i < numDPUs; ++i)
    if (i < chunksRemainder)
      chunksPerDPU[i] = minChunksPerDPU + 1;
    else
      chunksPerDPU[i] = minChunksPerDPU;

  // Partition nodes for each DPU.
  uint32_t *node_partitions = malloc((numDPUs + 1) * sizeof(uint32_t));
  node_partitions[0] = 0;
  for (uint32_t i = 1; i < numDPUs; ++i)
    node_partitions[i] = node_partitions[i - 1] + chunksPerDPU[i - 1] * 32;
  node_partitions[numDPUs] = csr.numRows;

  // Allocate DPUs.
  struct dpu_set_t set, dpu;
  PRINT_INFO("Allocating %d DPUs.", numDPUs);
  DPU_ASSERT(dpu_alloc(numDPUs, NULL, &set));
  DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(set, "totalChunks", totalChunks);

  uint32_t i;
  _DPU_FOREACH_I(set, dpu, i) {
    uint32_t from = node_partitions[i];   // index of rowPtrs for this dpu.
    uint32_t to = node_partitions[i + 1]; // index of rowPtrs for the next dpu.
    uint32_t numNodes = to - from;
    uint32_t numNeighbors = csr.rowPtrs[to] - csr.rowPtrs[from];
    uint32_t nodeChunkFrom = from / 32;
    uint32_t nodeChunkTo = (to + 31) / 32;            // ceil.
    uint32_t numChunks = nodeChunkTo - nodeChunkFrom; // == chunksPerDPU[i]
    uint32_t *nodePtrs = &csr.rowPtrs[from];
    uint32_t origin = nodePtrs[0]; // index of first neighbor.
    uint32_t *neighbors = &csr.colIdxs[origin];

    PRINT_INFO("DPU %d: nodeChunks = [%u, %u[", i, nodeChunkFrom, nodeChunkTo);

    dpu_copy_to_uint32(dpu, "numNodes", numNodes);
    dpu_copy_to_uint32(dpu, "numNeighbors", numNeighbors);
    dpu_copy_to_uint32(dpu, "numChunks", numChunks);
    dpu_copy_to_uint32(dpu, "nodeChunkFrom", nodeChunkFrom);
    dpu_copy_to_uint32(dpu, "nodeChunkTo", nodeChunkTo);
    dpu_copy_to_uint32(dpu, "origin", origin);

    dpu_insert_to_mram_uint32_array(dpu, "nodePtrs", nodePtrs, numNodes);
    dpu_insert_to_mram_uint32_array(dpu, "neighbors", neighbors, numNeighbors);
    dpu_insert_to_mram_uint32_array(dpu, "nodeLevels", 0, numNodes);
    dpu_insert_to_mram_uint32_array(dpu, "visited", 0, numChunks);
    dpu_insert_to_mram_uint32_array(dpu, "nextFrontier", nextFrontier,
                                    numChunks);
    dpu_insert_to_mram_uint32_array(dpu, "currFrontier", 0, totalChunks);
  }

  // BFS.
  uint32_t *nf = calloc(totalChunks, sizeof(uint32_t));
  uint32_t currentLevel = 0;
  uint32_t done = true;
  nextFrontier[0] = 0;

  while (true) {
    PRINT_INFO("Level %u", currentLevel);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Concatenate all nextFrontiers.
    uint32_t writeIdx;
    _DPU_FOREACH_I(set, dpu, i) {
      uint32_t nfNumChunks;
      // dpu_copy_from_uint32(dpu, "nfNumChunks", &nfNumChunks); // TODO
      dpu_copy_from_mram_uint32_array(dpu, "nextFrontier", &nf[writeIdx],
                                      nfNumChunks);
      writeIdx += nfNumChunks;

      // Get DPU logs.
      // PRINT_INFO("DPU %d:", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    // Check if done.
    for (uint32_t c = 0; c < totalChunks; ++c)
      if (nf[c] != 0) {
        done = false;
        break;
      }

    if (done)
      break;

    done = true;
    ++currentLevel;

    // Update currentLevel and nextFrontier of DPUs.
    dpu_copy_to_uint32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_copy_to_mram_uint32_array(dpu, "nextFrontier", nextFrontier,
                                    totalChunks);
    }

    // Clear nextFrontier. TODO: Do i need this?
    // for (uint32_t c = 0; c < totalChunks; ++c)
    //   nextFrontier[c] = 0;
  }

  // Get nodeLevels from each DPU.
  uint32_t *nodeLevels = malloc(csr.numRows * sizeof(uint32_t));
  uint32_t writeIdx;
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
    printf("%u\n", nodeLevels[node]);

  // Free resources.
  PRINT_INFO("Freeing resources");
  freeCOOMatrix(coo);
  freeCSRMatrix(csr);
  DPU_ASSERT(dpu_free(set));
  free(chunksPerDPU);
  free(node_partitions);
  free(nextFrontier);
  free(nf_dpu);
  free(nodeLevels);
}
