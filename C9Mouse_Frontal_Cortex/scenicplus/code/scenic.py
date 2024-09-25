import os
import pickle

work_dir = '/home/shar1105/projects/def-jrober27/shar1105/scenicplus/scenic_ftctx/'

cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

#from pycisTopic.clust_vis import *
#run_umap(cistopic_obj, target  = 'cell', rna_components=pd.DataFrame(adata_RNA.obsm['pca'], index = cistopic_obj.cell_names), rna_weight=0.8)
#plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['region'])

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu', plot=False)
binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=False)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj,
                                       selected_cells=None,
                                       selected_regions=None,
                                       scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj)
                                                 

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))
pickle.dump(imputed_acc_obj,
            open(os.path.join(work_dir, 'scATAC/imputed_acc_obj.pkl'), 'wb'))
pickle.dump(normalized_imputed_acc_obj,
            open(os.path.join(work_dir, 'scATAC/normalized_imputed_acc_obj.pkl'), 'wb'))

pickle.dump(variable_regions,
            open(os.path.join(work_dir, 'scATAC/variable_regions.pkl'), 'wb'))
