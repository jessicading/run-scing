# Run Single Cell INtegrative Gene regulatory network inference (SCING) from Seurat object/h5ad file and Leiden community detection

This repository provides a python script to run SCING from an h5ad file (pre-subset to one cell type) and Leiden community detection from the generated network. Below is a vignette to run these scripts and also convert a Seurat object to an AnnData object preceding running SCING (SCING uses scanpy). Leiden community detection outputs communities at resolutions from 0.0005 to 0.006 at 0.0005 increments. Review modules at different resolutions to choose a desired resolution. Usually, choose a resolution where no modules are larger than ~200 genes and consider modules no less than 10 genes (there will be several modules with just 1 gene, and I normally don't consider those as being part of a community, but they could be merged with the community of a neighbor to include all network genes in downstream analysis). The final output will be output to merged_networks/ (unless otherwise specified) and network communities (modules) will be output to modules/.

### Install SCING
First, clone SCING repo and set up conda environment following instructions here:
https://github.com/XiaYangLabOrg/SCING

### Convert Seurat object to AnnData object

```R
library(Seurat)
library(SeuratDisk)

seuratObject <- readRDS("seuratObject.rds")
# so that counts data goes properly into X of anndata
out <- CreateSeuratObject(counts = seuratObject@assays$RNA@counts)

# save to data directory of downloaded github repo^
# or move h5ad file to SCING/data/ after saving
SaveH5Seurat(out, filename = "SCING/data/sc_counts.h5Seurat")
Convert("SCING/data/sc_counts.h5Seurat", dest = "h5ad")
```

### Download RunSCING.py script

```bash
conda activate scing
cd SCING
mkdir modules # to write Leiden communites into
mkdir scripts
cd scripts
wget https://github.com/jessicading/run-scing/blob/master/RunSCING.py
```

### Set parameters (see parameters list below) if desired and submit script

Make sure h5ad file is in the SCING/data/ directory!

Run with default parameters (if necessary, adjust ```n_core``` and ```mem_per_core``` to match your computing environment):<br>
(First argument is the name of the counts file (without .h5ad extension)
```bash
python3 RunSCING.py sc_counts
```
For ~30000 cells with default parameters, the process will run for ~3.5 hours.

If desired, modify parameters (a full list of parameters is below).
```bash
python3 RunSCING.py sc_counts --n_supercells 500 --n_networks 100 --n_genes_net 4000 --knn 100 --n_pcs_net 10 --subsample_pct 0.7 --merge_consensus_pct 0.2 --n_genes_sc 2000 --n_pcs_sc 20 --n_core 12 --mem_per_core 2000000000
```

Alternatively, submit job using example script, RunSCING.sh (made for Sun Grid Engine queuing system).<br>
Parameters and requested computational resources can be adjusted in the RunSCING.sh file.
```bash
qsub RunSCING.sh sc_counts
```

Parameters:
- ```n_supercells```: number of supercells (meta-cells/aggregated cells) from which to build the network to handle sparsity issue
- ```n_networks```: number of intermediate networks to build (these networks are merged to form the consensus network)
- ```n_genes_net```: number of genes to utilize in creating the network, set to -1 to use all genes (recommended 2000-5000)
- ```knn```: k nearest neighbors - maximum potential gene neighbors (number of genes tested to be a neighbor to another gene in an edge)
- ```n_pcs_net```: number of PCs to be factored into gene network inference
- ```subsample_pct```: proportion of supercells to consider in the gene network inference (a different set of cells are used for each intermediate network)
- ```merge_consensus_pct```: minimum percentage of edges shared across networks to be included in the final merged network
- ```n_genes_sc```: number of genes to include in building supercells (top highly variable genes are chosen)
- ```n_pcs_sc```: number of PCs to include in building supercells
- ```n_core```: number of cores to use - if using HPC, this should be the same as the number of cores requested
- ```mem_per_core```: memory in bytes to use - if using HPC, this should be the same as the memory requested per core
- ```grn_prefix```: characters to append to beginning of individual network files
- ```grn_out_dir```: output directory for individual network files (default: ../networks/)
- ```merge_prefix```: characters to append to beginning of merged network file
- ```merge_out_dir```: output directory for merged network file (default: ../merged_networks/)

