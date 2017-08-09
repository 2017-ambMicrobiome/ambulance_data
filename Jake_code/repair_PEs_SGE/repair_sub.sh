#!/bin/bash
# qsub script for repair_PE

for f in `seq 1 8`; do 
	echo -e "qsubing repair_PE${f}"
	qsub repair_PE${f}.sh
	sleep 20s
done
