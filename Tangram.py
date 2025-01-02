import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg
import diopy
import os
os.environ['NUMEXPR_MAX_THREADS'] = '24'

#Read the input data
ad_sc=diopy.input.read_h5("snrna.merged.h5")
ad_sp=diopy.input.read_h5("scStereo.h5")

#Read the top marker genes
marker=pd.read_csv("snrna.top.marker.txt")
markers=np.reshape(marker.values,(-1,))
tg.pp_adatas(ad_sc,ad_sp,genes=markers,gene_to_lowercase=False)
assert ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']

#Mapping at the cluster level
ad_map=tg.map_cells_to_space(ad_sc,ad_sp,mode='clusters',cluster_label='seurat_clusters',device='cuda:0',num_epochs=500)
tg.project_cell_annotations(ad_map,ad_sp,annotation='seurat_clusters',threshold=0.5)
annotation_list=list(pd.unique(ad_sc.obs['seurat_clusters']))
ad_sp.obsm['tangram_ct_pred'].to_csv("scStereo2snrna.celltype.txt",sep="\t",index=True,header=True)

#Mapping at the gene level
ad_map=tg.map_cells_to_space(ad_sc,ad_sp,mode='cell',device='cuda:0',num_epochs=500)
ad_ge=tg.project_genes(adata_map=ad_map,adata_sc=ad_sc)
sc.write('snrna2scStereo.gene.expression.h5ad',ad_ge)