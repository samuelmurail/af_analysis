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

    pred_dir = directory

    for file in os.listdir(pred_dir):
        if file.endswith(".pdb"):
            tokens = file[: -4].split("_")
            model = int(tokens[4])
            weight = tokens[5]+'_'+tokens[6]
            seed = tokens[-1]

            pkl_score = os.path.join(
                pred_dir, f"result_model_{model}_{weight}_pred_{seed}.pkl",
            )
            np_score = np.load(pkl_score, allow_pickle=True)

            info_dict = {
                "pdb": os.path.join(pred_dir, file),
                "query": "unknown",
                "seed": seed,
                "model": model,
                "weight": weight,
                "recycle": "unknown",
                "pLDDT": np_score['plddt'].mean(),
                "pTM": np_score['ptm'],
                "ipTM": np_score['iptm'],
                "ranking_confidence": np_score['ranking_confidence'],
                "data_file": os.path.join(
                    pred_dir, f"result_model_{model}_{weight}_pred_{seed}.pkl",
                ),
            }

            log_dict_list.append(info_dict)

    log_pd = pd.DataFrame(log_dict_list)

    # To ensure that tests are consistent across different systems
    # we sort the dataframe by pdb
    log_pd = log_pd.sort_values(by=["pdb"]).reset_index(drop=True)
    return log_pd
