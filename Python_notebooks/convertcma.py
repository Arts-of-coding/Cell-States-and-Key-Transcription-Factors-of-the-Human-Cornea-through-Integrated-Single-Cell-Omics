# Move to new environment with newer version of scanpy (latestscanpy python 3.11)
import scanpy as sc
import numpy as np

adata = sc.read_h5ad('cornea_v1_pnas_nexus.h5ad')
influence_TFs = adata.obs.filter(like='_influence').columns.tolist() # get a list that ends_with "influence" collist=adata.obs.filter(like='_influence').columns.tolist()

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

sc.pp.log1p(adata)

# Rename columns and keep interesting ones
adata.obs["integrated_cell_states"] = adata.obs["latestjoined"]
adata.obs["integrated_clusters"] = adata.obs["latest"]
adata.obs["studies"] = adata.obs["batch"]


all_genes=adata.var.index.to_list()
relevant_states = ["integrated_cell_states", "integrated_clusters","studies","leiden","phase"] #,"predicted_doublets" gives error
list_qc = ["n_genes_by_counts", "total_counts","log10GenesPerUMI","pct_counts_mt","pct_counts_ribo","pct_counts_hb","G2M_score","S_score"]
genes_cell_cycle = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"]

# ['MLF1IP', 'FAM64A', 'HN1'] not in columns and thus removed from the full list

print(adata.obs)

adata.obs = adata.obs.round(1)

from scipy import sparse
adata.X = sparse.csr_matrix(adata.X)
adata.X

adata.obs = adata.obs[[*relevant_states, *list_qc, *influence_TFs]]

import pandas as pd
expression_treshold = 324 # 2x 162; Two times number of smallest cell state    
cluster_id = "integrated_cell_states"
res = pd.DataFrame(columns=adata.var_names.tolist(), index=adata.obs[cluster_id].astype("category").unique())

for clust in adata.obs[cluster_id].astype("category").unique():
	if adata.raw is not None:
		res.loc[clust] = adata[adata.obs[cluster_id].isin([clust]),:].raw.X.sum(0)
	else:
		res.loc[clust] = adata[adata.obs[cluster_id].isin([clust]),:].X.sum(0)
res.loc["sum"]=np.sum(res,axis=0).tolist()
res=res.transpose()
res=res.loc[res["sum"] > expression_treshold]

genes_expressed = res.index.tolist()
genes_expressed = list(set(genes_expressed)-set(genes_cell_cycle))

print("Amount of genes that remain: " + str(len(genes_expressed)))

#print(genes_expressed)

del res



plotdf = sc.get.obs_df(
        adata,
        keys=[*relevant_states, *list_qc, *influence_TFs, *genes_cell_cycle], #*all_genes, #, *s_genes,*g2m_genes ,"doublet_scores"
        obsm_keys=[("X_umap", 0), ("X_umap", 1)]
    )

# Convert adata.X to dataframe
import scipy.sparse
import pandas as pd

adata = adata[:, genes_expressed] # genes_expressed only use if needed otherwise use all genes!
df_genes = pd.DataFrame.sparse.from_spmatrix(adata.X)

del adata.X
del adata

df_genes.index = plotdf.index
df_genes.columns = genes_expressed
df_genes = df_genes.sparse.to_dense().round(1)

plotdf = plotdf.join(df_genes)

del df_genes

plotdf.to_parquet("cornea_v1_umap_clusres_scVI.parquet")

import polars as pl
import polars.selectors as cs

plotdf = pl.read_parquet("cornea_v1_umap_clusres_scVI.parquet")
plotdf = plotdf.with_columns(pl.col(genes_expressed).cast(pl.Float32).round(1))
plotdf = plotdf.with_columns(pl.col(genes_cell_cycle).cast(pl.Float32).round(1))
plotdf = plotdf.with_columns(pl.col(relevant_states).cast(pl.String))
plotdf = plotdf.with_columns(pl.col(influence_TFs).cast(pl.Float32).round(1))
plotdf = plotdf.with_columns(pl.col(list_qc).cast(pl.Float32).round(1))
plotdf = plotdf.with_columns(pl.col("total_counts").round(0))
plotdf = plotdf.with_columns(pl.col("n_genes_by_counts").round(0))
plotdf = plotdf.select([*relevant_states, *list_qc, *genes_expressed, *genes_cell_cycle, *influence_TFs, "X_umap-0", "X_umap-1"])

# resave with float32 for everything and check difference in size
plotdf.write_parquet("cornea_v1_umap_clusres_scVI_polars.parquet")

#print(plotdf.collect_schema())