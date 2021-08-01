#!/bin/bash

for window in 21 49 84; do 
	for reff in {9..14..1}; do
		for releffic in {0..100..10}; do	
			sbatch sensibility_serial.sh ${reff} ${window} ${releffic} &
		done
	done
done

