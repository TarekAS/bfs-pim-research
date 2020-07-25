all:
	gcc --std=c99 bfs-cpu/cpu.c -o bin/cpu
	gcc --std=c99 bfs-dpu/host/app.c -D "_POSIX_C_SOURCE=2" -o bin/app -lm `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/src-vtx-row bfs-dpu/dpu/src-vtx-row.c
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/src-vtx-col bfs-dpu/dpu/src-vtx-col.c
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/dst-vtx bfs-dpu/dpu/dst-vtx.c
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/edge bfs-dpu/dpu/edge.c

test:
	gcc --std=c99 bfs-dpu/host/test.c -D "_POSIX_C_SOURCE=2" -o bin/test -lm `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -DNR_TASKLETS=1 -O2 -o bin/testdpu bfs-dpu/dpu/test.c


clean:
	rm -f bin/cpu
	rm -f bin/app
	rm -f bin/src-vtx-row
	rm -f bin/src-vtx-col
	rm -f bin/dst-vtx
	rm -f bin/test
	rm -f bin/testdpu
