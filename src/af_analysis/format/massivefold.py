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
    ranking_files = [ 
        os.path.join(pred_dir, i) 
        for i in os.listdir(pred_dir) 
        if i.startswith('ranking_') and i.endswith('.json') 
    ]

    # get metrics preferably from massivefold ranking files instead of pickles
    global_df = pd.DataFrame()
    for ranking in ranking_files:
        json_ranking = json.load(open(ranking, 'r'))
        ranking_metric = list(json_ranking.keys())[0]
        models_metric = json_ranking[ranking_metric]
        metric_df = pd.DataFrame(models_metric.items(), columns=['prediction', ranking_metric])
        if global_df.empty:
            global_df = metric_df
        else:
            global_df = global_df.merge(metric_df, on="prediction")

    # compute ranking_confidence for AF3 which has a different 'ranking_score'
    if ({'iptm', 'ptm'}.issubset(set(global_df.columns))) and ('ranking_confidence' not in global_df.columns):
        global_df["ranking_confidence"] = 0.8*global_df["iptm"] + 0.2*global_df["ptm"]

    pkl_not_found = []
    for file in os.listdir(pred_dir):
        logger.info(f"Processing file: {file}")
        if file.endswith(".pdb") or file.endswith(".cif"):
            is_pkl = True
            # print(f"Processing file: {file}")

            confidence_score = None

            if not file.find("af3") != -1:
                tokens = file[:-4].split("_")
                model = int(tokens[4])
                pred_index = tokens.index("pred") if "pred" in tokens else len(tokens)
                weight = "_".join(tokens[5:pred_index])
                seed = int(tokens[-1])

                # print(tokens)
                # print(f"model: {model}, weight: {weight}, seed: {seed}")

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
                prediction = f"af3_seed_{model}_sample_{seed}_pred_{pred_num}"

                pkl_score = os.path.join(
                    pred_dir,
                    f"result_{prediction}.pkl",
                )
                if not os.path.exists(pkl_score):
                    pkl_score = os.path.join(
                        pred_dir,
                        f"light_pkl/result_{prediction}.pkl",
                    )
                    confidence_score = os.path.join(
                        pred_dir,
                        f"confidences/ranked_{rank}_{prediction}.json",
                    )
                    with open(confidence_score, "r") as json_file:
                        json_data = json.load(json_file)
                        # print(json_data.keys())

            ranking_scores = (
                global_df
                .set_index("prediction")
                .to_dict(orient='index')
            )
            if not os.path.exists(pkl_score):
                pkl_not_found.append(pkl_score)
                is_pkl = False
            else:
                np_score = np.load(pkl_score, allow_pickle=True)
                print(f"{np_score=}")

            if "mean_plddt" in global_df.columns:
                plddt = ranking_scores[prediction]["mean_plddt"]
            elif is_pkl and "plddt" in np_score:
                plddt = np_score["plddt"].mean()
            else:
                plddt = None

            if "ptm" in global_df.columns:
                ptm = ranking_scores[prediction]["ptm"]
            elif is_pkl and "ptm" in np_score:
                ptm = float(np_score["ptm"])
            elif is_pkl and confidence_score is not None:
                ptm = json_data["ptm"]
            else:
                ptm = None

            if "iptm" in global_df.columns:
                iptm = ranking_scores[prediction]["iptm"]
            elif is_pkl and "iptm" in np_score:
                iptm = float(np_score["iptm"])
            elif is_pkl and confidence_score is not None:
                iptm = json_data["iptm"]
            else:
                iptm = None

            if "ranking_confidence" in global_df.columns:
                ranking_confidence = ranking_scores[prediction]["ranking_confidence"]
            elif is_pkl and "ranking_confidence" in np_score:
                ranking_confidence = float(np_score["ranking_confidence"])
            elif is_pkl and confidence_score is not None:
                ranking_confidence = json_data["ranking_score"]
            else:
                ranking_confidence = None

            if is_pkl and "num_recycles" in np_score:
                num_recycles = np_score["num_recycles"]
            else:
                num_recycles = -1

            info_dict = {
                "pdb": os.path.join(pred_dir, file),
                "query": query,
                "seed": seed,
                "model": model,
                "weight": weight,
                "recycle": int(num_recycles),
                "pLDDT": plddt,
                "pTM": ptm,
                "ipTM": iptm,
                "ranking_confidence": ranking_confidence,
                "data_file": pkl_score,
            }

            log_dict_list.append(info_dict)

    if len(pkl_not_found) < 30:
        logger.warning(f"Non-existing score files (.pkl) ({len(pkl_not_found)}):\n {', '.join(pkl_not_found)}.")
    elif pkl_not_found:
        logger.warning(f"Detected {len(pkl_not_found)} non-existing score files (.pkl).")

    log_pd = pd.DataFrame(log_dict_list)

    # Dirty fix of "ranking_confidence" between 0 and 1.
    # Seems to be ok because ranking confidence is rarely below 1.
    log_pd["ranking_confidence"] = log_pd["ranking_confidence"].apply(
        lambda x: x * 100 if x < 1 else x
    )

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
