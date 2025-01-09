#!/usr/bin/env python3
# coding: utf-8

import os
import logging
from tqdm.auto import tqdm
import pandas as pd
import numpy as np

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

    pred_dir = os.path.join(directory)

    for file in os.listdir(pred_dir):
        if file.endswith(".cif"):
            token = file[4:-4].split("_")[-1]
            model = int(token[-1])

            npz_dict = np.load(os.path.join(pred_dir, f"scores.model_idx_{model}.npz"))

            # for key in npz_dict.keys():
            #    print(key, npz_dict[key])
            # print('ratio', ( 0.2* npz_dict['ptm'] + 0.8*npz_dict['iptm'] )/npz_dict['aggregate_score'])

            info_dict = {
                "pdb": os.path.join(pred_dir, file),
                "model": model,
            }

            for key in npz_dict.keys():
                info_dict[key] = npz_dict[key][0]

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
    log_pd.loc[:, "query"] = os.path.basename(os.path.normpath(pred_dir))

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
