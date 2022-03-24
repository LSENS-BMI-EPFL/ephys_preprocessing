#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: merge_ks_phy_output.py
@time: 18/03/2022 15:20
"""

# Imports
import os
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
PATH_ANALYSIS = 'M:/analysis/Axel_Bisi/mice_data/'

def merge_ks_phy_output(m_name):
    """
    Merge KS cluster info with phy updated cluster info into one dataframe matching KS output dimensions/indices.

    :param m_name: Mouse name.
    :return:
    """

    path_ks_output = os.path.join(PATH_ANALYSIS,'{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/ks25'.format(m_name))

    # Load cluster labels KS output
    cluster_label_ks = pd.read_csv(os.path.join(path_ks_output, 'cluster_KSlabel.tsv'), sep='\t')

    # Load cluster labels curated post-phy
    cluster_label_afterphy = pd.read_csv(os.path.join(path_ks_output, 'cluster_group.tsv'), sep='\t')

    # Load cluster total info phy
    cluster_info_phy = pd.read_csv(os.path.join(path_ks_output, 'cluster_info.tsv'), sep='\t')

    n_labels_ks = len(cluster_label_ks)
    n_labels_phy = len(cluster_label_afterphy)

    if n_labels_phy < n_labels_ks:
        print('Phy curation INCOMPLETE: some {} phy labels missing'.format(n_labels_ks-n_labels_phy))

        # Perform merge
        cluster_info_phy = cluster_info_phy.set_index(keys='cluster_id', drop=False, inplace=True)
        cluster_info_merge = cluster_info_phy.reindex(range(np.max(cluster_info_phy.cluster_id) + 1), fill_value=0, copy=True)
        cluster_info_df = cluster_info_merge

    else:
        print('Phy curation COMPLETE: just use phy output')
        cluster_info_df = cluster_info_phy

    # Save
    f_name = 'cluster_info_postphy.pkl' #keep same extension
    cluster_info_df.to_pickle(path=os.path.join(path_ks_output, f_name))

    return cluster_info_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--m_name', type=str, nargs='?', default='ABXXX', required=True)
    args = parser.parse_args()

    merge_ks_phy_output(args.m_name)
