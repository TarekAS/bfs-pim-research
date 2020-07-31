#include <assert.h>
#include <ctype.h>
#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <math.h>
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
#define PRINT_DEBUG(fmt, ...) fprintf(stderr, "\033[0;34mDEBUG:\033[0m   " fmt "\n", ##__VA_ARGS__)

#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

// Note: this is overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 16
#endif

enum Algorithm {
  SrcVtx = 0,
  DstVtx = 1,
  Edge = 2,
};

enum Partition {
  Row = 0,
  Col = 1,
  _2D = 2,
};

struct COO {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_edges;
  uint32_t *row_idxs;
  uint32_t *col_idxs;
};

struct CSR {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_edges;
  uint32_t *row_ptrs;
  uint32_t *col_idxs;
};

struct CSC {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_edges;
  uint32_t *col_ptrs;
  uint32_t *row_idxs;
};

/**
 * @fn dpu_insert_mram_array_u32
 * @brief Inserts data into the MRAM of a DPU at the last used MRAM address.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol where to copy the pointer of the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
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
 * @fn dpu_set_mram_array_u32
 * @brief Copy data to the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_set_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                      // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_get_mram_array_u32
 * @brief Copy data from the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param dst the host buffer where the data is copied.
 * @param length the number of elements in the array.
 */
void dpu_get_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                  // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_set_u32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param src the host buffer containing the data to copy
 */
void dpu_set_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t src) {
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &src, sizeof(uint32_t)));
}

/**
 * @fn dpu_get_u32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param dst the host buffer where the data is copied
 */
void dpu_get_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst) {
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, dst, sizeof(uint32_t)));
}

// Parse CLI args and options.
void parse_args(int argc, char **argv, int *num_dpu, enum Algorithm *alg, enum Partition *prt, char **bin_path, char **file) {
  bool is_prt_set = false;
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "n:a:p:")) != -1)
    switch (c) {
    case 'n':
      *num_dpu = atoi(optarg);
      if (num_dpu == 0 || *num_dpu % 8 != 0) {
        PRINT_ERROR("Number of DPUs must be a multiple of 8.");
        exit(1);
      }
      break;
    case 'a':
      if (strcmp(optarg, "src") == 0) {
        PRINT_INFO("Algorithm: Source-Vertex-based BFS.");
        *alg = SrcVtx;
        if (!is_prt_set)
          *prt = Row;
      } else if (strcmp(optarg, "dst") == 0) {
        PRINT_INFO("Algorithm: Destination-Vertex-based BFS.");
        *alg = DstVtx;
        if (!is_prt_set)
          *prt = Col;
      } else if (strcmp(optarg, "edge") == 0) {
        PRINT_INFO("Algorithm: Edge-based BFS.");
        *alg = Edge;
        if (!is_prt_set)
          *prt = _2D;
      } else {
        PRINT_ERROR("Incorrect -a argument. Supported algorithms: src | dst | edge");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "row") == 0) {
        PRINT_INFO("1D Row partitioning (source-nodes).");
        *prt = Row;
      } else if (strcmp(optarg, "col") == 0) {
        PRINT_INFO("1D Column partitioning (destination-nodes/neighbors).");
        *prt = Col;
      } else if (strcmp(optarg, "2d") == 0) {
        PRINT_INFO("2D partitioning (both source-nodes and destination-nodes).");
        *prt = _2D;
      } else {
        PRINT_ERROR("Incorrect -p argument. Supported partitioning: row | col | 2d");
        exit(1);
      }
      is_prt_set = true;
      break;
    case '?':
    default:
      PRINT_ERROR("Bad args. Usage: -n <num_dpu> -a <src|dst|edge> -p <row|col|2d>");
      exit(1);
    }

  int numargs = argc - optind;
  if (numargs != 1) {
    if (numargs > 1)
      PRINT_ERROR("Too many arguments!");
    else
      PRINT_ERROR("Too few arguments! Please provide data file name (COO-formatted matrix).");
    exit(1);
  }
  *file = argv[optind];

  if (*alg == SrcVtx && *prt == Row)
    *bin_path = "bin/src-vtx-row";
  else if (*alg == SrcVtx && *prt == Col)
    *bin_path = "bin/src-vtx-col";
  else if (*alg == SrcVtx && *prt == _2D)
    *bin_path = "bin/src-vtx-2d";
  else if (*alg == DstVtx && *prt == Row)
    *bin_path = "bin/dst-vtx-row";
  else if (*alg == DstVtx && *prt == Col)
    *bin_path = "bin/dst-vtx-col";
  else if (*alg == DstVtx && *prt == _2D)
    *bin_path = "bin/dst-vtx-2d";
  else if (*alg == Edge && *prt == Row)
    *bin_path = "bin/edge-row";
  else if (*alg == Edge && *prt == Col)
    *bin_path = "bin/edge-col";
  else if (*alg == Edge && *prt == _2D)
    *bin_path = "bin/edge-2d";
}

// Load coo-formated file into memory.
// Pads nodes to multiple of lcm of 32 and n, and to a minimum of n * 32.
struct COO load_coo(char *file, uint32_t n) {

  if (access(file, F_OK) == -1) {
    PRINT_ERROR("Could not find file %s.", file);
    exit(1);
  }

  PRINT_INFO("Loading COO-formated graph from %s.", file);
  struct COO coo;

  // Initialize COO from file.
  uint32_t num_nodes = 0;
  uint32_t num_edges = 0;

  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &num_nodes);
  fscanf(fp, "%u", &num_edges);
  coo.num_edges = num_edges;
  coo.row_idxs = malloc(num_edges * sizeof(uint32_t));
  coo.col_idxs = malloc(num_edges * sizeof(uint32_t));

  // Pad number of nodes to a minimum of n * 32.
  uint32_t min = n * 32;
  if (num_nodes < min) {
    PRINT_WARNING("Number of nodes too low. Setting number of nodes to %u.", min);
    num_nodes = min;
  }

  // Find Least Common Multiple of n and 32.
  uint32_t lcm = 0;
  uint32_t lar = (n > 32) ? n : 32;
  uint32_t small = (n > 32) ? 32 : n;

  for (int i = lar;; i += lar)
    if (i % small == 0) {
      lcm = i;
      break;
    }

  if (num_nodes % lcm != 0) {

    uint32_t padding = lcm - num_nodes % lcm;
    if (padding != 0) {
      PRINT_WARNING("Number of nodes must be multiple of %u. Padding with %u extra nodes.", lcm, padding);
      num_nodes += padding;
    }
  }

  coo.num_rows = num_nodes;
  coo.num_cols = num_nodes;

  // Read nonzeros.
  PRINT_INFO("Reading COO-formated graph - %u nodes, %u edges.", num_nodes, num_edges);

  uint32_t row_idx, col_idx;
  fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx); // Read first 2 integers of each line.

  uint32_t row_offset = row_idx; // Guarantee 0-indexed COO.
  coo.row_idxs[0] = row_idx - row_offset;
  coo.col_idxs[0] = col_idx - row_offset;

  for (uint32_t i = 1; i < num_edges; ++i) {
    fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx);
    coo.row_idxs[i] = row_idx - row_offset;
    coo.col_idxs[i] = col_idx - row_offset;
  }
  fclose(fp);

  return coo;
}

// Partition COO matrix into n COO matrices by col, or by row, or both (2D). Assumes n is even.
struct COO *partition_coo(struct COO coo, int n, enum Partition prt) {

  PRINT_INFO("Partitioning COO-matrix into %u parts.", n);

  struct COO *prts = malloc(n * sizeof(struct COO));

  // Initialize num_edges.
  for (int i = 0; i < n; ++i)
    prts[i].num_edges = 0;

  uint32_t num_rows = coo.num_rows;
  uint32_t num_cols = coo.num_cols;
  uint32_t row_div = 1;
  uint32_t col_div = 1;
  bool offset_row = false;
  bool offset_col = false;

  // Determine num_rows, num_cols, and num_edges per partition.
  switch (prt) {
  case Row:
    offset_row = true;
    row_div = n;
    num_rows /= row_div;
    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t row_idx = coo.row_idxs[i];
      prts[row_idx / num_rows].num_edges++;
    }
    break;

  case Col:
    offset_col = true;
    col_div = n;
    num_cols /= col_div;
    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t col_idx = coo.col_idxs[i];
      prts[col_idx / num_cols].num_edges++;
    }
    break;

  case _2D:
    offset_row = true;
    offset_col = true;

    // Find the two nearest factors of n.
    row_div = (uint32_t)sqrt(n);
    while (n % row_div != 0)
      row_div--;
    col_div = n / row_div;

    num_rows /= row_div;
    num_cols /= col_div;

    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t p_row = coo.row_idxs[i] / num_rows; // Partition row index.
      uint32_t p_col = coo.col_idxs[i] / num_cols; // Partition col index.
      uint32_t p = p_row * col_div + p_col;        // col-major index of coo.
      prts[p].num_edges++;
    }
    break;
  }

  // Initialize COO partitions.
  for (uint32_t i = 0; i < n; ++i) {
    prts[i].num_rows = num_rows;
    prts[i].num_cols = num_cols;
    prts[i].row_idxs = malloc(prts[i].num_edges * sizeof(uint32_t));
    prts[i].col_idxs = malloc(prts[i].num_edges * sizeof(uint32_t));
    prts[i].num_edges = 0; // We'll re-increment as we append data.
  }

  // Bin row and col pairs.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    uint32_t col_idx = coo.col_idxs[i];

    uint32_t p;

    if (prt == Row)
      p = row_idx / num_rows;
    else if (prt == Col)
      p = col_idx / num_cols;
    else if (prt == _2D) {
      uint32_t p_row = row_idx / num_rows;
      uint32_t p_col = col_idx / num_cols;
      p = p_row * col_div + p_col;
    }

    uint32_t idx = prts[p].num_edges;
    prts[p].row_idxs[idx] = row_idx;
    prts[p].col_idxs[idx] = col_idx;
    prts[p].num_edges++;
  }

  // Offset nodes.
  for (uint32_t p = 0; p < n; ++p) {
    uint32_t row_offset = offset_row ? p / col_div * num_rows : 0;
    uint32_t col_offset = offset_col ? p % col_div * num_cols : 0;

    for (uint32_t i = 0; i < prts[p].num_edges; ++i) {
      prts[p].row_idxs[i] -= row_offset;
      prts[p].col_idxs[i] -= col_offset;
    }
  }

  return prts;
}

// Converts COO matrix to CSR format.
struct CSR coo_to_csr(struct COO coo) {

  struct CSR csr;

  // Initialize fields.
  csr.num_rows = coo.num_rows;
  csr.num_cols = coo.num_cols;
  csr.num_edges = coo.num_edges;
  csr.row_ptrs = calloc((csr.num_rows + 1), sizeof(uint32_t));
  csr.col_idxs = malloc(csr.num_edges * sizeof(uint32_t));

  // Histogram row_idxs.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    csr.row_ptrs[row_idx]++;
  }

  // Prefix sum row_ptrs.
  uint32_t sum_before_next_row = 0;
  for (uint32_t row_idx = 0; row_idx < csr.num_rows; ++row_idx) {
    uint32_t sum_before_row = sum_before_next_row;
    sum_before_next_row += csr.row_ptrs[row_idx];
    csr.row_ptrs[row_idx] = sum_before_row;
  }
  csr.row_ptrs[csr.num_rows] = sum_before_next_row;

  // Bin the nonzeros.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    uint32_t nnzIdx = csr.row_ptrs[row_idx]++;
    csr.col_idxs[nnzIdx] = coo.col_idxs[i];
  }

  // Restore rowPtrs.
  for (uint32_t row_idx = csr.num_rows - 1; row_idx > 0; --row_idx)
    csr.row_ptrs[row_idx] = csr.row_ptrs[row_idx - 1];
  csr.row_ptrs[0] = 0;

  return csr;
}

// Converts COO matrix to CSC format.
struct CSC coo_to_csc(struct COO coo) {

  // Transpose COO matrix.
  struct COO coo_trs = {
      .col_idxs = coo.row_idxs,
      .row_idxs = coo.col_idxs,
      .num_cols = coo.num_rows,
      .num_rows = coo.num_cols};

  // Convert to CSR, then CSC.
  struct CSR csr = coo_to_csr(coo_trs);
  struct CSC csc = {
      .col_ptrs = csr.row_ptrs,
      .row_idxs = csr.col_idxs,
      .num_cols = csr.num_cols,
      .num_rows = csr.num_rows};

  return csc;
}

// Frees COO matrix.
void free_coo(struct COO coo) {
  free(coo.row_idxs);
  free(coo.col_idxs);
}

// Frees CSR matrix.
void free_csr(struct CSR csr) {
  free(csr.row_ptrs);
  free(csr.col_idxs);
}

// Frees CSC matrix.
void free_csc(struct CSC csc) {
  free(csc.col_ptrs);
  free(csc.row_idxs);
}

void bfs_dst_vtx(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t total_chunks) {

  uint32_t *next_frontier = calloc(total_chunks, sizeof(uint32_t));
  uint32_t done = true;
  uint32_t level = 0;

  while (true) {
    PRINT_INFO("Level %u", level);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Concatenate all next_frontiers.
    uint32_t write_idx = 0;
    int i = 0;
    _DPU_FOREACH_I(set, dpu, i) {
      uint32_t num_chunks;
      dpu_get_u32(dpu, "num_chunks", &num_chunks);
      dpu_get_mram_array_u32(dpu, "next_frontier",
                             &next_frontier[write_idx], num_chunks);
      write_idx += num_chunks;

      // Get DPU logs.
      // PRINT_INFO("DPU %d", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    // Check if done.
    for (uint32_t c = 0; c < total_chunks; ++c)
      if (next_frontier[c] != 0) {
        done = false;
        break;
      }
    if (done)
      break;

    done = true;
    ++level;

    // Update level and currentFrontier of DPUs.
    dpu_set_u32(set, "level", level);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_set_mram_array_u32(dpu, "currFrontier", next_frontier, total_chunks);
    }
  }

  // Free resources.
  free(next_frontier);
}

void bfs_src_vtx_row(struct dpu_set_t *set, struct dpu_set_t *dpu, int num_dpu, struct CSR *csr, enum Partition prt) {

  // Create BFS metadata.
  uint32_t num_nodes = csr[0].num_rows;       // Equal for all partitions.
  uint32_t num_chunks = num_nodes / 32;       // Number of chunks of 32 nodes per DPU.
  uint32_t total_nodes = num_nodes * num_dpu; // Total number of nodes.
  uint32_t total_chunks = total_nodes / 32;   // Total chunks of 32 nodes.

  if (num_chunks % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: num_chunks=%u not divisible by NR_TASKLETS=%u.", num_chunks, NR_TASKLETS);
    exit(1);
  }
  if (total_chunks % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: total_chunks=%u not divisible by NR_TASKLETS=%u.", total_chunks, NR_TASKLETS);
    exit(1);
  }

  uint32_t *next_frontier = calloc(total_chunks, sizeof(uint32_t));
  next_frontier[0] = 1; // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {

    // Copy required CSR data.
    dpu_set_u32(*dpu, "num_nodes", num_nodes);
    dpu_insert_mram_array_u32(*dpu, "node_ptrs", csr[i].row_ptrs, num_nodes + 1);
    dpu_insert_mram_array_u32(*dpu, "edges", csr[i].col_idxs, csr[i].num_edges);

    // Chunks data.
    dpu_set_u32(*dpu, "len_nf", total_chunks);
    dpu_set_u32(*dpu, "len_cf", num_chunks);
    dpu_set_u32(*dpu, "cf_from", num_chunks * i);
    dpu_set_u32(*dpu, "cf_to", num_chunks * (i + 1));

    // Initialize BFS data.
    dpu_insert_mram_array_u32(*dpu, "visited", 0, total_chunks);
    dpu_insert_mram_array_u32(*dpu, "curr_frontier", 0, num_chunks);
    dpu_insert_mram_array_u32(*dpu, "next_frontier", next_frontier, total_chunks);
    dpu_insert_mram_array_u32(*dpu, "node_levels", 0, num_nodes);
  }

  // Start BFS.
  uint32_t *nf_tmp = calloc(total_chunks, sizeof(uint32_t));
  next_frontier[0] = 0;
  uint32_t level = 0;
  int done = true;

  while (true) {

    // Launch DPUs.
    PRINT_INFO("Level %u", level);
    DPU_ASSERT(dpu_launch(*set, DPU_SYNCHRONOUS));

    // Union next_frontiers.
    uint32_t i = 0;
    _DPU_FOREACH_I(*set, *dpu, i) {

      dpu_get_mram_array_u32(*dpu, "next_frontier", nf_tmp, total_chunks);
      for (uint32_t c = 0; c < total_chunks; ++c) {
        next_frontier[c] |= nf_tmp[c];
        if (next_frontier[c] != 0)
          done = false;
      }

      // DPU_ASSERT(dpu_log_read(*dpu, stdout));
    }

    if (done)
      break;
    done = true;
    ++level;

    // Update level and next_frontier of DPUs.
    dpu_set_u32(*set, "level", level);
    _DPU_FOREACH_I(*set, *dpu, i) {
      dpu_set_mram_array_u32(*dpu, "next_frontier", next_frontier, total_chunks);
    }

    // Clear next_frontier.
    for (uint32_t c = 0; c < total_chunks; ++c)
      next_frontier[c] = 0;
  }

  // Get node_levels from each DPU.
  uint32_t *node_levels = calloc(total_nodes, sizeof(uint32_t));
  uint32_t write_idx = 0;

  i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    uint32_t num_nodes;
    dpu_get_u32(*dpu, "num_nodes", &num_nodes);
    dpu_get_mram_array_u32(*dpu, "node_levels", &node_levels[write_idx], num_nodes);
    write_idx += num_nodes;
  }

  // Print node levels.
  PRINT_INFO("Output:");
  for (uint32_t node = 0; node < total_nodes; ++node) {
    uint32_t level = node_levels[node];
    if (node != 0 && level == 0) // Filters out "padded" rows.
      continue;
    printf("node_levels[%u]=%u\n", node, node_levels[node]);
  }

  free(nf_tmp);
  free(next_frontier);
}

void bfs_src_vtx_col(struct dpu_set_t *set, struct dpu_set_t *dpu, int num_dpu, struct CSR *csr, enum Partition prt) {

  // Create BFS metadata.
  uint32_t num_nodes = csr[0].num_rows;     // Num destination nodes (not to be confused with num edges).
  uint32_t num_neighbors = csr[0].num_cols; // Num destination nodes (not to be confused with num edges).
  uint32_t total_neighbors = num_neighbors * num_dpu;
  uint32_t num_chunks = num_neighbors / 32;
  uint32_t total_chunks = num_neighbors * num_dpu / 32;

  if (num_chunks % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: num_chunks=%u not divisible by NR_TASKLETS=%u.", num_chunks, NR_TASKLETS);
    exit(1);
  }
  if (total_chunks % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: total_chunks=%u not divisible by NR_TASKLETS=%u.", total_chunks, NR_TASKLETS);
    exit(1);
  }

  uint32_t *next_frontier = calloc(total_chunks, sizeof(uint32_t));
  next_frontier[0] = 1; // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {

    dpu_set_u32(*dpu, "dpu_id", i);

    // Copy required CSR data.
    dpu_set_u32(*dpu, "num_neighbors", num_neighbors);
    dpu_insert_mram_array_u32(*dpu, "node_ptrs", csr[i].row_ptrs, num_nodes + 1);
    dpu_insert_mram_array_u32(*dpu, "edges", csr[i].col_idxs, csr[i].num_edges);

    // Chunks data.
    dpu_set_u32(*dpu, "len_nf", num_chunks);
    dpu_set_u32(*dpu, "len_cf", total_chunks);
    dpu_set_u32(*dpu, "cf_from", num_chunks * i);
    // dpu_set_u32(*dpu, "cf_to", len_nf * (i + 1));

    // Initialize BFS data.
    dpu_insert_mram_array_u32(*dpu, "visited", 0, num_chunks);
    dpu_insert_mram_array_u32(*dpu, "curr_frontier", next_frontier, total_chunks);
    dpu_insert_mram_array_u32(*dpu, "next_frontier", 0, num_chunks);
    dpu_insert_mram_array_u32(*dpu, "node_levels", 0, num_neighbors);
  }

  // Start BFS.
  uint32_t level = 0;
  int done = true;

  while (true) {

    // Launch DPUs.
    PRINT_INFO("Level %u", level);
    DPU_ASSERT(dpu_launch(*set, DPU_SYNCHRONOUS));

    // Concatenate all next_frontiers.
    uint32_t i = 0;
    uint32_t write_idx = 0;
    _DPU_FOREACH_I(*set, *dpu, i) {

      dpu_get_mram_array_u32(*dpu, "next_frontier", &next_frontier[write_idx], num_chunks);
      write_idx += num_chunks;

      // DPU_ASSERT(dpu_log_read(*dpu, stdout));
    }

    // Check if done.
    for (uint32_t c = 0; c < total_chunks; ++c)
      if (next_frontier[c] != 0) {
        done = false;
        break;
      }

    if (done)
      break;
    done = true;
    ++level;

    // Update level and curr_frontier of DPUs.
    dpu_set_u32(*set, "level", level);
    _DPU_FOREACH_I(*set, *dpu, i) {
      dpu_set_mram_array_u32(*dpu, "curr_frontier", next_frontier, total_chunks);
    }
  }

  // Print node levels.
  uint32_t *node_levels = calloc(total_neighbors, sizeof(uint32_t));

  i = 0;
  uint32_t write_idx = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    dpu_get_mram_array_u32(*dpu, "node_levels", &node_levels[write_idx], num_neighbors);
    write_idx += num_neighbors;
  }

  PRINT_INFO("Output:");
  for (uint32_t node = 0; node < total_neighbors; ++node) {
    uint32_t level = node_levels[node];
    if (node != 0 && level == 0) // Filters out "padded" rows.
      continue;
    printf("node_levels[%u]=%u\n", node, level);
  }

  // Free resources.
  free(next_frontier);
  free(node_levels);
}

void bfs_src_vtx_2d(struct dpu_set_t *set, struct dpu_set_t *dpu, int num_dpu, struct CSR *csr, enum Partition prt) {

  exit(1);

  // Create BFS metadata.
  uint32_t num_nodes = csr[0].num_rows;
  uint32_t num_neighbors = csr[0].num_cols;
  uint32_t total_nodes = num_nodes * num_dpu;
  uint32_t total_neighbors = num_neighbors * num_dpu;
  uint32_t len_frontier = total_nodes / 32;
  uint32_t len_nf = num_neighbors / 32;
  uint32_t len_cf = num_nodes / 32;

  if (len_nf % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: len_nf=%u not divisible by NR_TASKLETS=%u.", len_nf, NR_TASKLETS);
    exit(1);
  }
  if (len_cf % NR_TASKLETS != 0) {
    PRINT_ERROR("Error: len_cf=%u not divisible by NR_TASKLETS=%u.", len_cf, NR_TASKLETS);
    exit(1);
  }

  uint32_t *next_frontier = calloc(len_frontier, sizeof(uint32_t));
  next_frontier[0] = 1; // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {

    dpu_set_u32(*dpu, "dpu_id", i);

    // Copy required CSR data.
    dpu_insert_mram_array_u32(*dpu, "node_ptrs", csr[i].row_ptrs, num_nodes + 1);
    dpu_insert_mram_array_u32(*dpu, "edges", csr[i].col_idxs, csr[i].num_edges);

    // Chunks data.
    dpu_set_u32(*dpu, "len_nf", len_nf);
    dpu_set_u32(*dpu, "len_cf", len_cf);

    // Initialize BFS data.
    dpu_insert_mram_array_u32(*dpu, "visited", 0, len_nf);
    dpu_insert_mram_array_u32(*dpu, "curr_frontier", 0, len_cf);
    dpu_insert_mram_array_u32(*dpu, "next_frontier", next_frontier, len_nf);
    dpu_insert_mram_array_u32(*dpu, "node_levels", 0, num_neighbors);
  }
}

void start_src_vtx(struct dpu_set_t *set, struct dpu_set_t *dpu, struct COO *coo, int num_dpu, enum Partition prt) {

  // Convert COO partitions to CSR.
  struct CSR *csr = malloc(num_dpu * sizeof(struct CSR));
  for (int i = 0; i < num_dpu; ++i) {
    csr[i] = coo_to_csr(coo[i]);
    free_coo(coo[i]);
  }

  if (prt == Row)
    bfs_src_vtx_row(set, dpu, num_dpu, csr, prt);
  else if (prt == Col)
    bfs_src_vtx_col(set, dpu, num_dpu, csr, prt);
  else
    bfs_src_vtx_2d(set, dpu, num_dpu, csr, prt);
}

int main(int argc, char **argv) {

  int num_dpu = 8;
  enum Algorithm alg = SrcVtx;
  enum Partition prt = Row;
  char *bin_path;
  char *file = NULL;
  parse_args(argc, argv, &num_dpu, &alg, &prt, &bin_path, &file);

  PRINT_INFO("Allocating %u DPUs, %u tasklets each.", num_dpu, NR_TASKLETS);
  struct dpu_set_t set, dpu;
  DPU_ASSERT(dpu_alloc(num_dpu, NULL, &set));
  DPU_ASSERT(dpu_load(set, bin_path, NULL));

  struct COO coo = load_coo(file, num_dpu);
  struct COO *coo_prts = partition_coo(coo, num_dpu, prt);
  free_coo(coo);

  if (alg == SrcVtx)
    start_src_vtx(&set, &dpu, coo_prts, num_dpu, prt);
  else if (alg == DstVtx) {
    // TODO
  } else if (alg == Edge) {
    // TODO
  }

  DPU_ASSERT(dpu_free(set));
  return 0;
}
