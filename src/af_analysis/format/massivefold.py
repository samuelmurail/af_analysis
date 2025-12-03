#!/usr/bin/env python3
# coding: utf-8

import os
import logging
from tqdm.auto import tqdm
import pandas as pd
import numpy as np
import json

from af_analysis import data

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
    query = pred_dir.split("/")[-1]

    for file in os.listdir(pred_dir):
        logger.info(f"Processing file: {file}")
        if file.endswith(".pdb") or file.endswith(".cif"):

            confidence_score = None

            if not file.find("af3") != -1:
                tokens = file[:-4].split("_")
                model = int(tokens[4])
                weight = tokens[5] + "_" + tokens[6]
                seed = int(tokens[-1])

                pkl_score = os.path.join(
                    pred_dir,
                    f"result_model_{model}_{weight}_pred_{seed}.pkl",
                )
                if not os.path.exists(pkl_score):
                    pkl_score = os.path.join(
                        pred_dir,
                        f"light_pkl/result_model_{model}_{weight}_pred_{seed}.pkl",
                    )
            else: # Special case for af3
                tokens = file[:-4].split("_")
                model = int(tokens[4]) # In reality it is a the seed
                weight = "af3"
                seed = int(tokens[-3])
                pred_num = int(tokens[-1])
                rank = int(tokens[1])

                pkl_score = os.path.join(
                    pred_dir,
                    f"result_model_{model}_{weight}_pred_{seed}.pkl",
                )
                if not os.path.exists(pkl_score):
                    pkl_score = os.path.join(
                        pred_dir,
                        f"light_pkl/result_af3_seed_{model}_sample_{seed}_pred_{pred_num}.pkl",
                    )
                    confidence_score = os.path.join(
                        pred_dir,
                        f"confidences/ranked_{rank}_af3_seed_{model}_sample_{seed}_pred_{pred_num}.json",
                    )
                    with open(confidence_score, "r") as json_file:
                        json_data = json.load(json_file)
                        # print(json_data.keys())

            if not os.path.exists(pkl_score):
                logger.warning(f"Score file {pkl_score} does not exist.")
                continue

            np_score = np.load(pkl_score, allow_pickle=True)

            if "plddt" in np_score:
                plddt = np_score["plddt"].mean()
            else:
                plddt = None
            
            if "ptm" in np_score:
                ptm = float(np_score["ptm"])
            elif confidence_score is not None:
                ptm = json_data["ptm"]
            else:
                ptm = None
            
            if "iptm" in np_score:
                iptm = float(np_score["iptm"])
            elif confidence_score is not None:
                iptm = json_data["iptm"]
            else:
                iptm = None
                
            if "ranking_confidence" in np_score:
                ranking_confidence = float(np_score["ranking_confidence"])
            elif confidence_score is not None:
                ranking_confidence = json_data["ranking_score"]
            else:
                ranking_confidence = None
            
            if "num_recycles" not in np_score:
                np_score["num_recycles"] = -1

            info_dict = {
                "pdb": os.path.join(pred_dir, file),
                "query": query,
                "seed": seed,
                "model": model,
                "weight": weight,
                "recycle": int(np_score["num_recycles"]),
                "pLDDT": plddt,
                "pTM": ptm,
                "ipTM": iptm,
                "ranking_confidence": ranking_confidence,
                "data_file": pkl_score,
            }

            log_dict_list.append(info_dict)

    log_pd = pd.DataFrame(log_dict_list)

    # To ensure that tests are consistent across different systems
    # we sort the dataframe by pdb
    log_pd = log_pd.sort_values(by=["pdb"]).reset_index(drop=True)
    return log_pd


def read_full_directory(directory):
    """Extract pdb list from a directory and return as a dictionary.

    Parameters
    ----------
    directory : str
        Path to the directory containing the pdb files.

    Returns
    -------
    log_dict : dict
        Dictionary containing the information extracted from the directory.

    """

    logger.info(f"Reading full MassiveFold {directory}")

    subfolders = [f.path for f in os.scandir(directory) if f.is_dir()]
    log_pd_list = []

    for folder in subfolders:
        # print(f"Reading {folder}")
        if not os.path.basename(folder).startswith("msa"):
            log_pd_list.append(read_dir(folder))

    log_dict = pd.concat(log_pd_list, ignore_index=True)
    if os.path.basename(directory) != "":
        log_dict["query"] = os.path.basename(directory)
    else:
        log_dict["query"] = os.path.basename(os.path.dirname(directory))

    return log_dict
