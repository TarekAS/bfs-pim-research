
# What is this?

Project to benchmark DPU-accelerated Breadth-First Search, compared to classic CPU-only parallel BFS implementations. Part of my MSc Thesis on Near-Memory Acceleration of Graph Algorithms.

**Note:** Requires `UPMEM DPU SDK 2020.x.x` to be installed.

```
  $ make all
  $ ./bin/app <datafile>  # bench DPU-accelerated version.
  $ ./bin/cpu <datafile>  # bench CPU-only version.
  $ make clean
```
Where datafile is a 0-indexed COO-formated graph.

You can pass these optional arguments to `./bin/app`
```
-n <int>        Number of DPUs to use. Must be a multiple of 8.
-a <string>     BFS algorithm variation to use: src_vtx (default) or dst_vtx or edge.
-p <partition>  Whether to partition the nodes: none (default) or 1d or 2d.
```

# Directory Structure

```
bfs-cpu/  # classic (cpu-only) BFS implementation.
bfs-dpu/  # dpu-accelerated BFS implementation.
  dpu/       # task code (runs on DPU)
  host/      # main app code (runs on CPU)
```

