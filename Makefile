all:
	gcc --std=c99 bfs-cpu/cpu.c -o bin/cpu
	gcc --std=c99 bfs-dpu/host/app.c -D "_POSIX_C_SOURCE=2" -o bin/app `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/src-vtx bfs-dpu/dpu/src-vtx.c
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/dst-vtx bfs-dpu/dpu/dst-vtx.c

test:
	gcc --std=c99 bfs-dpu/host/test.c -D "_POSIX_C_SOURCE=2" -o bin/test `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -DNR_TASKLETS=1 -O2 -o bin/testdpu bfs-dpu/dpu/test.c


clean:
	rm bin/cpu
	rm bin/app
	rm bin/src-vtx
	rm bin/dst-vtx
