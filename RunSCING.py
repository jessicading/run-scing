# Example run:
# python3 RunSCING.py mic --n_supercells 500 --n_networks 100 --n_genes_net 4000 --knn 100 \
#   --n_pcs_net 10 --subsample_pct 0.7 --merge_consensus_pct 0.2 --n_genes_sc 2000 --n_pcs_sc 20 --n_core 12 \
#   --mem_per_core 2000000000

import os
import warnings
import sys
import leidenalg
import igraph as ig
import argparse
import pickle

sys.path.insert(1, '../src/')

from supercellHelpers import *
from buildGRNHelpers import *
from MergeNetworksHelpers import *

warnings.filterwarnings("ignore")


def get_network(ct_name, n_scs, n_nets, n_genes_network, n_knn, n_pcs_network, sub_pct, merge_con_pct,
                n_genes_supercells, n_pcs_supercells, grn_prefix_name, grn_out_directory, merge_prefix_name,
                merge_out_directory, cores, core_mem):
    n_threads = cores
    os.environ["MKL_NUM_THREADS"] = str(n_threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(n_threads)
    os.environ["OMP_NUM_THREADS"] = str(n_threads)

    adata = sc.read('../data/' + ct_name + '.h5ad')

    adata_merged = supercell_pipeline(adata,
                                      ngenes=n_genes_supercells,
                                      npcs=n_pcs_supercells,
                                      ncell=n_scs,
                                      verbose=True)

    all_edges = []
    for i in range(n_nets):
        print(i)
        scing = grnBuilder(adata_merged,
                           ngenes=n_genes_network,
                           nneighbors=n_knn,
                           npcs=n_pcs_network,
                           subsample_perc=sub_pct,
                           prefix=grn_prefix_name,
                           outdir=grn_out_directory,
                           ncore=cores,
                           mem_per_core=core_mem,
                           verbose=True)

        scing.subsample_cells()

        scing.filter_genes()
        scing.filter_gene_connectivities()
        scing.build_grn()

        all_edges.append(scing.edges)

    with open('../intermediate_data/' + ct_name + '_edges.pickle', 'wb') as handle:
        pickle.dump(all_edges, handle, protocol=pickle.HIGHEST_PROTOCOL)

    merger = NetworkMerger(adata_merged,
                           networks=all_edges,
                           minimum_edge_appearance_threshold=merge_con_pct,
                           prefix=merge_prefix_name,
                           outdir=merge_out_directory,
                           ncore=cores,
                           mem_per_core=core_mem,
                           verbose=True)

    merger.preprocess_network_files()
    merger.remove_reversed_edges()
    merger.remove_cycles()
    merger.get_triads()
    merger.remove_redundant_edges()

    merger.edge_df.sort_values(by='importance',
                               ascending=False)

    merger.edge_df.to_csv(merge_out_directory + "/" + ct_name + "_edges.csv.gz")

    network = merger.edge_df.iloc[:, [0, 1]]

    g = ig.Graph.TupleList([tuple(x) for x in network.values],
                           directed=True)

    for res in np.linspace(start=0.0005, stop=0.006, num=12):
        partition = leidenalg.find_partition(g,
                                             leidenalg.CPMVertexPartition,
                                             resolution_parameter=res)
        groups = np.array(partition.membership)
        genes = g.vs['name']
        gene_membership = pd.DataFrame(np.array([genes, groups]).T,
                                       columns=['genes', 'cluster_membership'])
        gene_membership = gene_membership.sort_values('genes')
        gene_membership.cluster_membership = gene_membership.cluster_membership.astype(str)
        gene_membership.to_csv('../modules/' + ct_name + '.res' + str(res) + '.csv.gz',
                               index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build GRN with SCING and find network modules with leiden clustering"
    )

    parser.add_argument('cell_type_name',
                        metavar='cell_type_name',
                        type=str, nargs=1, action='store',
                        default=None)

    parser.add_argument('--n_supercells',
                        metavar='n_supercells',
                        type=int, nargs=1, action='store',
                        default=500)

    parser.add_argument('--n_networks',
                        metavar='n_networks',
                        type=int, nargs=1, action='store',
                        default=100)

    parser.add_argument('--n_genes_net',
                        metavar='n_genes_net',
                        type=int, nargs=1, action='store',
                        default=4000)

    parser.add_argument('--knn',
                        metavar='knn',
                        type=int, nargs=1, action='store',
                        default=100)

    parser.add_argument('--n_pcs_net',
                        metavar='n_pcs_net',
                        type=int, nargs=1, action='store',
                        default=10)

    parser.add_argument('--subsample_pct',
                        metavar='subsample_pct',
                        type=float, nargs=1, action='store',
                        default=0.7)

    parser.add_argument('--merge_consensus_pct',
                        metavar='merge_consensus_pct',
                        type=float, nargs=1, action='store',
                        default=0.2)

    parser.add_argument('--n_genes_sc',
                        metavar='n_genes_sc',
                        type=int, nargs=1, action='store',
                        default=2000)

    parser.add_argument('--n_pcs_sc',
                        metavar='n_pcs_sc',
                        type=int, nargs=1, action='store',
                        default=20)

    parser.add_argument('--grn_prefix',
                        metavar='grn_prefix',
                        type=str, nargs=1, action='store',
                        default=None)

    parser.add_argument('--grn_out_dir',
                        metavar='grn_out_dir',
                        type=str, nargs=1, action='store',
                        default=None)

    parser.add_argument('--merge_prefix',
                        metavar='merge_prefix',
                        type=str, nargs=1, action='store',
                        default=None)

    parser.add_argument('--merge_out_dir',
                        metavar='merge_out_dir',
                        type=str, nargs=1, action='store',
                        default=None)

    parser.add_argument('--n_core',
                        metavar='n_core',
                        type=int, nargs=1, action='store',
                        default=12)

    parser.add_argument('--mem_per_core',
                        metavar='mem_per_core',
                        type=int, nargs=1, action='store',
                        default=int(2e9))

    args = parser.parse_args()

    cell_type_name = args.cell_type_name[0]
    n_supercells = args.n_supercells
    n_networks = args.n_networks
    n_genes_net = args.n_genes_net
    knn = args.knn
    n_pcs_net = args.n_pcs_net
    subsample_pct = args.subsample_pct
    merge_consensus_pct = args.merge_consensus_pct
    n_genes_sc = args.n_genes_sc
    n_pcs_sc = args.n_pcs_sc
    grn_prefix = args.grn_prefix
    if grn_prefix is None:
        grn_prefix = cell_type_name
    grn_out_dir = args.grn_out_dir
    if grn_out_dir is None:
        grn_out_dir = '../networks/' + cell_type_name
    merge_prefix = args.merge_prefix
    if merge_prefix is None:
        merge_prefix = cell_type_name
    merge_out_dir = args.merge_out_dir
    if merge_out_dir is None:
        merge_out_dir = '../merged_networks'
    n_core = args.n_core
    mem_per_core = args.mem_per_core

    get_network(cell_type_name, n_supercells, n_networks, n_genes_net, knn, n_pcs_net, subsample_pct,
                merge_consensus_pct, n_genes_sc, n_pcs_sc, grn_prefix, grn_out_dir, merge_prefix, merge_out_dir,
                n_core, mem_per_core)
