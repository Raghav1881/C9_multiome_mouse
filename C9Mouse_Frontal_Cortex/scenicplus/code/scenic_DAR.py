import os
import pickle
import sys
from pycisTopic.diff_features import *

work_dir = '/home/shar1105/projects/def-jrober27/shar1105/scenicplus/scenic_ftctx/'
tmp_dir = '/home/shar1105/scratch/'

cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
imputed_acc_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/imputed_acc_obj.pkl'), 'rb'))
normalized_imputed_acc_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/normalized_imputed_acc_obj.pkl'), 'rb'))
variable_regions = pickle.load(open(os.path.join(work_dir, 'scATAC/variable_regions.pkl'), 'rb'))
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))

print('Calculating DARs for each region...')
markers_dict_region = find_diff_features(cistopic_obj, 
                                         imputed_acc_obj, 
                                         variable='subclass', 
                                         var_features=variable_regions, 
                                         split_pattern = '-',
                                         log2fc_thr=0.3)
print('Calculating DARs for each genotype...')
markers_dict_genotype = find_diff_features(cistopic_obj,
                                           imputed_acc_obj,
                                           variable='genotype',
                                           var_features=variable_regions,
                                           split_pattern = '-',
                                           log2fc_thr=0.3,
                                           contrasts = [[["C9HET"], ["C9WT"]], 
                                                        [["C9KO"], ["C9WT"]]])
print('Calculating DARs for the contrast C9WT vs C9HET and C9KO')
contrasts = [
    [["Astro_C9HET"], ["Astro_C9WT"]],
    [["Astro_C9KO"], ["Astro_C9WT"]],
    [["L2/3 IT CTX_C9HET"], ["L2/3 IT CTX_C9WT"]],
    [["L2/3 IT CTX_C9KO"], ["L2/3 IT CTX_C9WT"]],
    [["L4/5 IT CTX_C9HET"], ["L4/5 IT CTX_C9WT"]],
    [["L4/5 IT CTX_C9KO"], ["L4/5 IT CTX_C9WT"]],
    [["L5 IT CTX_C9HET"], ["L5 IT CTX_C9WT"]],
    [["L5 IT CTX_C9KO"], ["L5 IT CTX_C9WT"]],
    [["L5 PT CTX_C9HET"], ["L5 PT CTX_C9WT"]],
    [["L5 PT CTX_C9KO"], ["L5 PT CTX_C9WT"]],
    [["L5/6 NP CTX_C9HET"], ["L5/6 NP CTX_C9WT"]],
    [["L5/6 NP CTX_C9KO"], ["L5/6 NP CTX_C9WT"]],
    [["L6 CT CTX_C9HET"], ["L6 CT CTX_C9WT"]],
    [["L6 CT CTX_C9KO"], ["L6 CT CTX_C9WT"]],
    [["L6 IT CTX_C9HET"], ["L6 IT CTX_C9WT"]],
    [["L6 IT CTX_C9KO"], ["L6 IT CTX_C9WT"]],
    [["L6b CTX_C9HET"], ["L6b CTX_C9WT"]],
    [["L6b CTX_C9KO"], ["L6b CTX_C9WT"]],
    [["Lamp5_C9HET"], ["Lamp5_C9WT"]],
    [["Lamp5_C9KO"], ["Lamp5_C9WT"]],
    [["Micro_C9HET"], ["Micro_C9WT"]],
    [["Micro_C9KO"], ["Micro_C9WT"]],
    [["Oligo_C9HET"], ["Oligo_C9WT"]],
    [["Oligo_C9KO"], ["Oligo_C9WT"]],
    [["OPC_C9HET"], ["OPC_C9WT"]],
    [["OPC_C9KO"], ["OPC_C9WT"]],
    [["Peri_C9HET"], ["Peri_C9WT"]],
    [["Peri_C9KO"], ["Peri_C9WT"]],
    [["Pvalb_C9HET"], ["Pvalb_C9WT"]],
    [["Pvalb_C9KO"], ["Pvalb_C9WT"]],
    [["Pvalb Vipr2_C9HET"], ["Pvalb Vipr2_C9WT"]],
    [["Pvalb Vipr2_C9KO"], ["Pvalb Vipr2_C9WT"]],
    [["Sncg_C9HET"], ["Sncg_C9WT"]],
    [["Sncg_C9KO"], ["Sncg_C9WT"]],
    [["Sst_C9HET"], ["Sst_C9WT"]],
    [["Sst_C9KO"], ["Sst_C9WT"]],
    [["Sst Chodl_C9HET"], ["Sst Chodl_C9WT"]],
    [["Sst Chodl_C9KO"], ["Sst Chodl_C9WT"]],
    [["Vip_C9HET"], ["Vip_C9WT"]],
    [["Vip_C9KO"], ["Vip_C9WT"]],
    [["VLMC_C9HET"], ["VLMC_C9WT"]],
    [["VLMC_C9KO"], ["VLMC_C9WT"]]
]

markers_dict_C9_genotypes = find_diff_features(cistopic_obj,
                                               imputed_acc_obj,
                                               variable='subclass_genotype',
                                               var_features=variable_regions,
                                               split_pattern = '-', 
                                               log2fc_thr=0.3,
                                               contrasts = contrasts)
pickle.dump(markers_dict_C9_genotypes, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict_C9_genotypes_raw.pkl'), 'wb'))
markers_dict_C9_genotypes = {key: df for key, df in markers_dict_C9_genotypes.items() if not df.empty}
pickle.dump(markers_dict_C9_genotypes, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict_C9_genotypes.pkl'), 'wb'))

# Save
if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(markers_dict_region, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict_region.pkl'), 'wb'))
pickle.dump(markers_dict_genotype, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict_genotype.pkl'), 'wb'))
pickle.dump(markers_dict_C9_genotypes, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict_C9_genotypes.pkl'), 'wb'))

# Perform motif enrichment
# Regions of interest: Astro, L2/3 IT CTX, L4/5 IT CTX, L5 IT CTX, L6 CT CTX, Micro, Oligo, OPC
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['DARs_region'] = {}
region_sets['DARs_genotype'] = {}
region_sets['DARs_C9_genotypes'] = {}

for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict_region.keys():
    regions = markers_dict_region[DAR].index[markers_dict_region[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_region'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict_genotype.keys():
    regions = markers_dict_genotype[DAR].index[markers_dict_genotype[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_genotype'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict_C9_genotypes.keys():
    regions = markers_dict_C9_genotypes[DAR].index[markers_dict_C9_genotypes[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_C9_genotypes'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

pickle.dump(region_sets, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_sets.pkl'), 'wb'))
rankings_db = os.path.join(work_dir, 'scenicplus_obj/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db = os.path.join(work_dir, 'scenicplus_obj/mm10_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(work_dir, 'scenicplus_obj/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl')

# Make directory for motifs
if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget
sys.stderr = open(os.devnull, "w")  # silence stderr
run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = False,
    n_cpu = 30,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust')
sys.stderr = sys.__stderr__  # unsilence stderr
