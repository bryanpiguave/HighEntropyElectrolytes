#!/bin/csh

#$ -pe smp 4                 # Request 4 cores
#$ -N dmc_only               # Job name
#$ -cwd                      # Run from current directory
#$ -V                        # Export environment variables

# --- Load modules ---
module load intel
module load intelmpi
module load lammps

# --- Define output directory with timestamp ---
set start_time = `date +%Y%m%d_%H%M%S`
set outdir = "${JOB_NAME}_${JOB_ID}_${start_time}"
mkdir -p $outdir

# --- Run LAMMPS ---
lmp < in.dmc_only > ${outdir}/out.dmc_only_${JOB_ID}

# --- Move additional files into output directory ---
mv log.lammps $outdir/ >& /dev/null
mv dump* $outdir/ >& /dev/null
mv traj* $outdir/ >& /dev/null
mv *.txt $outdir/ >& /dev/null
mv *.dat $outdir/ >& /dev/null

# --- Completion Message ---
echo "Job finished at `date`" | tee -a $outdir/finished.log
echo "Job ID: $JOB_ID"        | tee -a $outdir/finished.log
echo "Job Name: $JOB_NAME"    | tee -a $outdir/finished.log
