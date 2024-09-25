#!/bin/bash
#SBATCH --account=def-jrober27
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=502G

module load StdEnv/2020 gcc python/3.9.6
module load arrow rust
source /home/shar1105/projects/def-jrober27/shar1105/scenicplus/venv/bin/activate

#mkdir -p scenicplus_obj && cd scenicplus_obj
#rankings_database_url='https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather'
#scores_database_url='https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather'
#motif_database_url = 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
#wget "$rankings_database_url"
#wget "$scores_database_url"
#wget "$motif_database_url"
#cd ..
#python scenic.py
python scenic.py
