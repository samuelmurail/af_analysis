#!/usr/bin/env python3
# coding: utf-8

import os
import logging
from tqdm.auto import tqdm
import pandas as pd
import numpy as np
import json

logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)


def read_dir(directory):
    """Extract pdb list from a directory.

    Parameters
    ----------
    directory : str
        Path to the directory containing the pdb files.

    Returns
    -------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the directory.

    """

    logger.info(f"Reading {directory}")

    log_dict_list = []

    out_folders = [
        os.path.join(directory, folder)
        for folder in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, folder))
    ]

    trunk_folders = [
        os.path.join(out_folder, f)
        for out_folder in out_folders
        for f in os.listdir(out_folder)
        if os.path.isdir(os.path.join(out_folder, f)) and "trunk" in f
    ]

    model_path_list = [
        os.path.join(trunk_folder, f)
        for trunk_folder in trunk_folders
        for f in os.listdir(trunk_folder)
        if f.endswith(".cif")
    ]

    for model_path in model_path_list:

        query = model_path.split("/")[-3]
        trunk_id = int(model_path.split("/")[-2][6:])
        model_id = int(model_path.split("/")[-1][15:-4])
        data_file = model_path.replace("/pred.", "/scores.").replace(".cif", ".npz")

        info_dict = {
            "pdb": model_path,
            "query": query,
            "trunk_id": trunk_id,
            "model_id": model_id,
            "data_file": data_file,
        }

        npz_dict = np.load(data_file)
        for key in npz_dict.keys():
            info_dict[key] = npz_dict[key][0]


        # for key in npz_dict.keys():
        #    print(key, npz_dict[key])
        # print('ratio', ( 0.2* npz_dict['ptm'] + 0.8*npz_dict['iptm'] )/npz_dict['aggregate_score'])

        """
            "pTM": npz_dict['ptm'],
            "ipTM": npz_dict['iptm'],
            "ranking_confidence": npz_dict['aggregate_score'],
            'per_chain_ptm',
            'per_chain_pair_iptm', 'has_inter_chain_clashes', 'chain_chain_clashes'
            
        }
        """
        log_dict_list.append(info_dict)

    log_pd = pd.DataFrame(log_dict_list)
    # log_pd.loc[:, "query"] = os.path.basename(os.path.normpath(pred_dir))

    # Update column names
    log_pd = log_pd.rename(
        columns={
            # https://github.com/jwohlwend/boltz/issues/73
            # "confidence_score": "ranking_confidence",
            "ptm": "pTM",
            "iptm": "ipTM",
            "aggregate_score": "ranking_confidence",
        }
    )

    # To ensure that tests are consistent across different systems
    # we sort the dataframe by pdb
    log_pd = log_pd.sort_values(by=["pdb"]).reset_index(drop=True)
    return log_pd
