import os
import numpy as np
import pandas as pd
import prody
from scipy.sparse import csr_matrix
from geometricus import MomentInvariants, SplitType


"""
spektral.data.Graph
spektral.data.graph.Graph(x=None, a=None, e=None, y=None)
x: the node features
a: the adjacency matrix
e: the edge features
y: the labels

Edge attributes can be stored in a dense format as arrays of shape
(n_nodes, n_nodes, n_edge_features)

x: np.array, the node features (shape (n_nodes, n_node_features));
a: np.array or scipy.sparse matrix, the adjacency matrix (shape (n_nodes, n_nodes));
e: np.array, the edge features (shape (n_nodes, n_nodes, n_edge_features) or (n_edges, n_edge_features));
y: np.array, the node or graph labels (shape (n_nodes, n_labels) or (n_labels, ));
"""
def add_residue_col(ppi):
    residue_df = ppi[[c for c in ppi.columns if "residue_aa" in c]]
    int_labels  = np.where(residue_df.values == 1)[1]
    int_2_res = dict( zip(range(max(int_labels)+1), [c.split('_')[-1] for c in residue_df.columns]))
    ppi['res_label'] = int_labels
    ppi['res_label'] = ppi['res_label'].replace(int_2_res)
    return ppi

def parse_df(ppi_df):
    atom_groups = []
    invariants_kmer = []
    invariants_radius = []
    # iterates through all rows of dataframe and adds geometricus features
    for dom_str in ppi_df.domain.unique():
        pdb_id = dom_str[:-3]
        chain = dom_str[-3].upper()
        atom_grp = prody.parsePDB(pdb_id, chain=chain)
        atom_groups.append(atom_grp)
        invariants_kmer.append(
            MomentInvariants.from_prody_atomgroup(dom_str, atom_grp, split_type=SplitType.KMER, split_size=16))
        invariants_radius.append(
            MomentInvariants.from_prody_atomgroup(dom_str, atom_grp, split_type=SplitType.RADIUS, split_size=10))
    return atom_groups, invariants_kmer, invariants_radius

def construct_distance_matrix(moment_invariant):
    # this function constructs a distance matrix for all features that
    n_res = moment_invariant.length
    distances = np.zeros(n_res, n_res)
    for i in range(n_res):
        for j in range(i+1, n_res):
            res1, res2 = moment_invariant.coordinates[i], moment_invariant.coordinates[j]
            distances[i,j] = distances[j,i] = np.linalg.norm(res1 - res2)

def make_sparse_matrix(moment_invariants):
    n_domains = len(moment_invariants)
    n_res = sum([invariant.length for invariant in moment_invariants])
    shape = (n_res, n_res) # size of the matrix is the number of amino acids across all proteins
    row = []
    col = []
    data = []
    row_counter = 0
    for dom_idx in range(n_domains):
        domain = moment_invariants[dom_idx]
        for i in range(domain.length):
            for j in range(i+1, i+domain.length):
                res1, res2 = domain.coordinates[i], domain.coordinates[j]
                dist = np.linalg.norm(res1 - res2)
                row.append(i)
                col.append(j)
                data.append(dist)
                if i!=j:
                    row.append(row_counter)
                    col.append(j)
                    data.append(dist)
            row_counter +=1
    sparse_matrix = csr_matrix((data, (row, col)), shape=shape)
    return sparse_matrix
                    





def make_sparse_matrix_example():
    row = np.array([0, 0, 1, 2, 2, 2])
    col = np.array([0, 2, 2, 0, 1, 2])
    data = np.array([1, 2, 3, 4, 5, 6])
    csr_matrix((data, (row, col)), shape=(3, 3))

def make_amino_dicts():
    residue_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    res_to_idx = dict(zip(residue_letters, range(len(residue_letters))))
    idx_to_res = dict(zip(range(len(residue_letters)),  residue_letters))
    return res_to_idx, idx_to_res

def make_one_hot_res(one_protein, res_to_idx):
    sequence = [res_to_idx[res] for res in list(one_protein.sequence)]
    shape = (one_protein.length, len(res_to_idx))
    one_hot = np.zeros(shape)
    one_hot[range(one_protein.length), sequence] = 1
    return one_hot

def construct_invariant_feature_set(invariants_kmer, invariants_radius, res_to_idx):
    # one_hot_res will act as labels

    for i in range(len(invariants_kmer)):
        one_hot_res = make_one_hot_res(invariants_kmer[i], res_to_idx)
        one_protein_features = np.append(invariants_kmer[i].moments, invariants_radius[i].moments, axis=1)
        one_protein_features = np.append()




if __name__ == '__main__':
    features_ff = ['scons', 'avg_scons', 'sc5_gs', 'sc5_scons', 'conserved_hotspot_struc_neighbourhood',
                   'conserved_surface_hotspot_struc_neighbourhood', 'highly_conserved_struc_neighbourhood',
                   'highly_conserved_surface_struc_neighbourhood', 'pocket_conserved_struc_neighbourhood',
                   'pocket_surface_conserved_struc_neighbourhood',
                   'avg_charged', 'avg_cx', 'avg_dpx', 'avg_electric_effect', 'avg_flexibility', 'avg_hydropathicity',
                   'avg_hydrophobicity', 'avg_polarity', 'avg_surface_residues', 'avg_surrounding_hydrophobicity',
                   'dist_to_hotspot', 'dist_to_surface', 'hotspot_struc_neighbourhood', 'pocket_struc_neighbourhood',
                   'surface_residues_struc_neighbourhood', 'min_dist_to_cleft123', 'min_dist_to_cleft_1',
                   'min_dist_to_cleft_2', 'min_dist_to_cleft_3', 'surface_residue_rsa', 'cleft_residue',
                   'hydrophobic_aa', 'polar_aa',
                   'alpha', 'betweenness', 'bulkiness', 'charge', 'cleft_depth', 'cleft_num', 'closeness', 'degree',
                   'foldx_alascan', 'free_energy_solution', 'hydration_potential', 'hydropathicity', 'hydropathy_index',
                   'hydrophobicity', 'hydrophobicity_psaia', 'kappa', 'localised_electrical_effect', 'max_cx',
                   'max_dpx', 'min_cx', 'min_dpx', 'nhBonds_ptr', 'oBonds_ptr', 'phi', 'polarity', 'psi', 'resTco',
                   'res_bfactor_n', 'rsa_allatoms', 'rsa_mainchain', 'rsa_nonpolar', 'rsa_polar', 'rsa_totside',
                   'van_der_waals_vol_normalised',
                   'dssp_type_B', 'dssp_type_H', 'dssp_type_NO_PRED', 'dssp_type_T',
                   'A_pssm_ff', 'A_wop_ff', 'C_pssm_ff', 'C_wop_ff', 'D_pssm_ff', 'D_wop_ff', 'E_pssm_ff', 'E_wop_ff',
                   'F_pssm_ff', 'F_wop_ff', 'G_pssm_ff', 'G_wop_ff', 'H_pssm_ff', 'H_wop_ff', 'I_pssm_ff', 'I_wop_ff',
                   'K_pssm_ff', 'K_wop_ff', 'L_pssm_ff', 'L_wop_ff', 'M_pssm_ff', 'M_wop_ff', 'N_pssm_ff', 'N_wop_ff',
                   'P_pssm_ff', 'P_wop_ff', 'Q_pssm_ff', 'Q_wop_ff', 'R_pssm_ff', 'R_wop_ff', 'S_pssm_ff', 'S_wop_ff',
                   'T_pssm_ff', 'T_wop_ff', 'V_pssm_ff', 'V_wop_ff', 'W_pssm_ff', 'W_wop_ff', 'Y_pssm_ff', 'Y_wop_ff',
                   'gapless_match_to_pseudocounts_ff', 'info_per_pos_ff',
                   ]

    basedir = '../cath-funsite-predictor/datasets/PPI/'
    df = pd.read_csv(basedir + 'PPI_validation_dataset.csv')
    df = add_residue_col(df)
    x = df[features_ff].values
    os.chdir('../protein')
    atom_groups, invariants_kmer, invariants_radius = parse_df(df)
    res_to_idx, idx_to_res = make_amino_dicts(invariants_kmer)
    breakpoint_var = True



