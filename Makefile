NR_TASKLETS ?= 11
BLOCK_SIZE ?= 32
BENCHMARK_CYCLES ?= false
BENCHMARK_TIME ?= false

all:
	gcc --std=c11 bfs-dpu/host/bfs.c -Wall -Wextra -g -O3 -D "_POSIX_C_SOURCE=2" -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -o bin/bfs -lm `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/top-down bfs-dpu/dpu/top-down.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/bottom-up bfs-dpu/dpu/bottom-up.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/edge bfs-dpu/dpu/edge.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/top-down-dma bfs-dpu/dpu/top-down-dma.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/bottom-up-dma bfs-dpu/dpu/bottom-up-dma.c
	dpu-upmem-dpurte-clang -Wall -Wextra -g -O2 -DNR_TASKLETS=$(NR_TASKLETS) -DBLOCK_SIZE=$(BLOCK_SIZE) -DBENCHMARK_CYCLES=$(BENCHMARK_CYCLES)  -DBENCHMARK_TIME=$(BENCHMARK_TIME) -O2 -o bin/edge-dma bfs-dpu/dpu/edge-dma.c

clean:
	rm -f bin/bfs
	rm -f bin/top-down
	rm -f bin/top-down-dma
	rm -f bin/bottom-up
	rm -f bin/bottom-up-dma
	rm -f bin/edge
	rm -f bin/edge-dma
