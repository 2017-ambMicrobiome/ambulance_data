#!/bin/bash
# qsub_clark.sh
# Script for submitting other clark scripts

for f in {1..4}; do
	echo -e "\nSubmitting script: clark0${f}.sh"
	qsub clark0${f}_orphan.sh
	sleep 10s
done
