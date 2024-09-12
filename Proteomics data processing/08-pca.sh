#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N PCA

cd /rds/general/user/rw1317/projects/hda_23-24/live/TDS/Project3/Proteomics_data_processing/Scripts

module load anaconda3/personal
source activate Renv

Rscript 08-pca.R

