
# What is this?

Project to benchmark DPU-accelerated Breadth-First Search, compared to classic CPU-only parallel BFS implementations. Part of my MSc Thesis on Near-Memory Acceleration of Graph Algorithms.

**Note:** Requires `UPMEM DPU SDK 2020.x.x` to be installed.

```
  $ make all
  $ ./bin/app <datafile>  # bench DPU-accelerated version.
  $ ./bin/cpu <datafile>  # bench CPU-only version.
  $ make clean
```
Where `datafile` is a COO-formated graph datafile, such as the first line contains the number of nodes and edges, and the subsequent lines contain the row/col index pairs separated by a space.

For example:
```
<NUM_NODES> <NUM_EDGES>
0 1
0 2
0 11
1 3
1 4
...
```


You can pass these optional arguments to `./bin/app`
```
-n <int>        Number of DPUs to use. Must be a multiple of 8. Max 64 if using emulator.
-a <string>     BFS algorithm variation to use: src | dst | edge.
-p <partition>  The manner in which to partition the nodes per DPU: row | col | 2d.
```

- `src` is for source-vertex BFS (i.e. by nodes). Uses `row` partitioning by default. This is the standard way.
- `dst` is for destination-vertex BFS (i.e. by neighbors). Uses `col` partitioning by default.
- `edge` is for edge-based BFS. Uses `2d` partitioning by default.

# Directory Structure

```
bfs-cpu/  # classic (cpu-only) BFS implementation.
bfs-dpu/  # dpu-accelerated BFS implementation.
  dpu/       # task code (runs on DPU)
  host/      # main app code (runs on CPU)
```

