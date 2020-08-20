#!/bin/bash

NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p row data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p col data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p 2d data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p row data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p col data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p 2d data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p row data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p col data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p 2d data/loc-gowalla_edges.txt > res && diff -q res res_true_gowalla && rm res


NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p row data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p col data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a src -p 2d data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p row data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p col data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a dst -p 2d data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p row data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p col data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res
# NR_TASKLETS=11 BLOCK_SIZE=512 make all && ./bin/app -n 64 -a edge -p 2d data/loc-brightkite_edges.txt > res && diff -q res res_true_brightkite && rm res

