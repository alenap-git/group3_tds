#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N QC

cd /rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Scripts

module load anaconda3/personal
source activate r413

Rscript 02-quality_control.R

