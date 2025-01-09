#!/usr/bin/env python3
# coding: utf-8

import os
import logging
import json
from tqdm.auto import tqdm
import pandas as pd

logger = logging.getLogger()
logger.setLevel(logging.INFO)


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

    for file in os.listdir(directory):
        if file.endswith(".cif"):
            token = file[:-4].split("_")

            model = int(token[-1])
            query = "_".join(token[1:-2])
            json_score = os.path.join(
                directory, f"fold_{query}_summary_confidences_{model}.json"
            )
            with open(json_score, "r") as f_in:
                json_dict = json.load(f_in)

            info_dict = {
                "pdb": os.path.join(directory, file),
                "query": query,
                "model": model,
                "data_file": os.path.join(
                    directory, f"fold_{query}_full_data_{model}.json"
                ),
            }
            info_dict.update(json_dict)
            log_dict_list.append(info_dict)

    log_pd = pd.DataFrame(log_dict_list)

    # Update column names
    log_pd = log_pd.rename(
        columns={
            "ranking_score": "ranking_confidence",
            "ptm": "pTM",
            "iptm": "ipTM",
        }
    )

    # To ensure that tests are consistent across different systems
    # we sort the dataframe by pdb
    log_pd = log_pd.sort_values(by=["pdb"]).reset_index(drop=True)
    return log_pd
