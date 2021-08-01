#!/bin/bash

for window in {21..84..7}; do 
	for releffic in {0..14..1}; do	
		sbatch sensibility_serial.sh ${window} ${releffic} &
	done
done


