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

    out_folders = [
        os.path.join(directory, folder)
        for folder in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, folder))
    ]
    model_path_list = [
        os.path.join(out_folder, f) + "/model.cif"
        for out_folder in out_folders
        for f in os.listdir(out_folder)
        if os.path.isdir(os.path.join(out_folder, f)) and "seed" in f and "sample" in f
    ]

    for model_path in model_path_list:

        seed_sample = model_path.split("/")[-2]  # for instance seed-1_sample-3
        query = model_path.split("/")[-3]
        json_score = model_path.replace(
            "/model.cif", "/summary_confidences.json"
        )  # summary_confidences.json file
        with open(json_score, "r") as f_in:
            json_dict = json.load(f_in)

        info_dict = {
            "pdb": model_path,
            "query": query,
            "seed_sample": seed_sample,
            "data_file": model_path.replace("/model.cif", "/confidences.json"),
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
