#!/bin/bash
#SBATCH -n 1                #Neste caso estão sendo pedidos 64 cores, configure como desejar. 
#SBATCH --ntasks-per-node=1 #Para jobs de até 64 cores manter o mesmo número que o numero de cores pedidos, para jobs com mais de 64 cores, manter sempre em 64.
#SBATCH -p serial            #Pedido para fila média, se desejar outra fila, modifique.
#SBATCH --job-name=Vax_serial   #Mudar o JOBNAME para o nome desejado para o job.

module purge
ml restore R

#cd $SLURM_SUBMIT_DIR

#mkdir /scratch-local/$SLURM_JOB_ID
#cp -rf * /scratch-local/$SLURM_JOB_ID
#cd /scratch-local/$SLURM_JOB_ID


#i=$1   ####reff
#j=$2   ####vax.window.days
#k=$3   ####rel.effic
#Rscript sensibility.R $i $j $k

i=$1 ###vax.window.days
j=$2 ###rel.effic

Rscript sensibility_grid.R $i $j

#cp -rf * $SLURM_SUBMIT_DIR
#cd $SLURM_SUBMIT_DIR
#rm -rf /scratch-local/$SLURM_JOB_ID
