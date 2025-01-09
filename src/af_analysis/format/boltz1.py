#!/usr/bin/env python3
# coding: utf-8

import os
import logging
from tqdm.auto import tqdm
import pandas as pd
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

    pred_dir = os.path.join(directory, "predictions")

    for query in os.listdir(pred_dir):
        for file in os.listdir(os.path.join(os.path.join(pred_dir, query))):
            if file.endswith(".cif"):
                token = file[len(query) + 1 : -4].split("_")[-1]
                model = int(token[-1])

                json_score = os.path.join(
                    pred_dir,
                    os.path.join(query, f"confidence_{query}_model_{model}.json"),
                )
                with open(json_score, "r") as f_in:
                    json_dict = json.load(f_in)

                info_dict = {
                    "pdb": os.path.join(pred_dir, os.path.join(query, file)),
                    "query": query,
                    "model": model,
                }
                plddt_file = os.path.join(
                    pred_dir, os.path.join(query, f"plddt_{query}_model_{model}.npz")
                )
                if os.path.isfile(plddt_file):
                    info_dict["plddt"] = plddt_file
                pae_file = os.path.join(
                    pred_dir, os.path.join(query, f"pae_{query}_model_{model}.npz")
                )
                if os.path.isfile(pae_file):
                    info_dict["data_file"] = pae_file
                pde_file = os.path.join(
                    pred_dir, os.path.join(query, f"pae_{query}_model_{model}.npz")
                )
                if os.path.isfile(pde_file):
                    info_dict["pde"] = pde_file

                info_dict.update(json_dict)
                log_dict_list.append(info_dict)

    log_pd = pd.DataFrame(log_dict_list)

    # Update column names
    log_pd = log_pd.rename(
        columns={
            # https://github.com/jwohlwend/boltz/issues/73
            # "confidence_score": "ranking_confidence",
            "ptm": "pTM",
            "iptm": "ipTM",
        }
    )
    log_pd["ranking_confidence"] = 0.2 * log_pd["pTM"] + 0.8 * log_pd["ipTM"]

    # To ensure that tests are consistent across different systems
    # we sort the dataframe by pdb
    log_pd = log_pd.sort_values(by=["pdb"]).reset_index(drop=True)
    return log_pd
