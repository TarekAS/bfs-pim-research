NR_TASKLETS ?= 11
BLOCK_SIZE ?= 32
BENCHMARK_CYCLES ?= false
BENCHMARK_TIME ?= false

all:
	gcc --std=c99 bfs-cpu/bfs-cpu.c -o bin/bfs-cpu
	gcc --std=c11 bfs-dpu/host/bfs.c -Wall -Wextra -g -O3 -D "_POSIX_C_SOURCE=2" -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -o bin/bfs -lm `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/src-vtx bfs-dpu/dpu/src-vtx.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/dst-vtx bfs-dpu/dpu/dst-vtx.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/edge bfs-dpu/dpu/edge.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/src-vtx-dma bfs-dpu/dpu/src-vtx-dma.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/dst-vtx-dma bfs-dpu/dpu/dst-vtx-dma.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/edge-dma bfs-dpu/dpu/edge-dma.c

clean:
	rm -f bin/bfs-cpu
	rm -f bin/bfs
	rm -f bin/src-vtx
	rm -f bin/dst-vtx
	rm -f bin/edge
	rm -f bin/src-vtx-dma
	rm -f bin/dst-vtx-dma
	rm -f bin/edge-dma
