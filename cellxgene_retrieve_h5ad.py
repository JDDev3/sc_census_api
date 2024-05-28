# Information on how to create conda environment for 
# create new conda environment
# conda create -n cellxgene python=3.11
# conda install -n cellxgene pip
# conda activate cellxgene
# pip install -U cellxgene-census
# python
import cellxgene_census
import scanpy as sc
with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census = census,
        organism = "Mus Musculus",
        obs_value_filter = "tissue == 'spleen' and cell_type in ['T cell']"
    )
adata.write_h5ad("test.h5ad")
