
all:
	gcc -I$(UPMEM_HOME)/usr/share/upmem/include/host host/src/app.c -o bin/app -L$(UPMEM_HOME)/usr/lib/ -ldpu -ldpuelf -lpthread -DDPU_TARGET_TYPE=FUNCTIONAL_SIMULATOR -DDPU_TARGET_PROFILE=\"\"
	clang --target=dpu-upmem-dpurte -o bin/task.bin dpu/src/task.c

test:
	./bin/app

clean:
	rm bin/task.bin
	rm bin/app

