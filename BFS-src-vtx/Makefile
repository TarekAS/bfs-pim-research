all:
	gcc --std=c99 host/src/app.c -o bin/app `dpu-pkg-config --cflags --libs dpu`
	dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o bin/task dpu/src/task.c

test:
	./bin/app

target_cpu:
	gcc --std=c99 cpu/cpu.c -o bin/cpu

test_cpu:
	./bin/cpu

clean:
	rm bin/task
	rm bin/app
	rm bin/cpu
