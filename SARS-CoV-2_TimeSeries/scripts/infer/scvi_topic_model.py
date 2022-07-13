
import pandas as pd
import pyarrow.parquet as pq
import scvi
import scanpy as sc
import anndata as ad

cell_feature_columns_TS = pd.read_csv(
    "product/cell_feature_columns_TS_202008.tsv",
    sep = "\t")
cell_metadata_columns_TS = pd.read_csv(
    "product/cell_metadata_columns_TS_202008.tsv",
    sep = "\t")
cell_features_path = "product/TS2_2M_Cell_MasterDataTable.parquet"

cell_features = pq.read_table(source=cell_features_path).to_pandas()
cell_features = ad.AnnData(
    cell_features[celL_feature_columns_TS],
    obs = cell_meta_data)

#
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=10e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`:

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer="counts", n_top_genes=1000, subset=True)

n_topics = 10

scvi.model.AmortizedLDA.setup_anndata(adata, layer = "counts")
model = scvi.model.AmortizedLDA(adata, n_topics = n_topics)

model.train()


topic_prop = model.get_latent_representation()
topic_prop.head()
