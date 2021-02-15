
# What is this?

Project to benchmark DPU-accelerated Breadth-First Search.

**Note:** Tested with `UPMEM DPU SDK 2020.3.x`.



```
  $ make
```
Optional environment variables for make:
- `BENCHMARK_TIME=true` benchmarks the BFS duration in seconds.
- `BENCHMARK_CYCLES=true` counts the number of DPU cycles per BFS iteration.
- `NR_TASKLETS=<integer>` sets the number of tasklets per DPU (max 24, recommended 11).
- `BLOCK_SIZE=<multiple_of_8>` sets the MRAM DMA block size (multiple of 8, max 512 bytes).

```
  $ ./bin/bfs -n <num_dpu> -a <base_algorithm> -p <partitioning> -o <output_result_path> <datafile>
```
Notes:
- `num_dpu` must be a multiple of 8. Emulator has a limit of 64.
- `datafile` COO-formated graph (adjacency list) that is tab separated, and sorted by the first column then the second column. The first line contains the number of nodes followed by the number of edges. See example below.
- `base_algorithm` is the base BFS algorithm to use, with options:
  - `top` for vertex-centric top-down BFS.
  - `bot` for vertex-centric bottom-up BFS.
  - `edge` for edge-centric BFS.
- `partitioning` the way the adjacency matrix is partitioned over the DPUs, with options:
  - `row` partition the source nodes (i.e. nodes).
  - `col` partition the destination nodes (i.e. neighbors).
  - `2d` partition both source nodes and destination nodes in tiles.

Example datafile:
```
<NUM_NODES> <NUM_EDGES>
0 1
0 2
0 11
1 3
1 4
...
```

# Directory Structure

```
bfs-dpu/
  host/      # Host code (CPU side)
  dpu/       # Task code (DPU side). Contains optimized DMA versions and non-optimized reader-friendly versions.
```
