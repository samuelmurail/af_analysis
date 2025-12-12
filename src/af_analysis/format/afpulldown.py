#!/usr/bin/env python3
# coding: utf-8

import os
import logging
import json
from tqdm.auto import tqdm
import pandas as pd
import numpy as np

logger = logging.getLogger()
logger.setLevel(logging.INFO)

weigths = [
    "multimer_v1",
    "multimer_v2",
    "multimer_v3",
]


def read_dir(directory, query=None):
    """Extract pdb list from a directory.

    Parameters
    ----------
    directory : str
        Path to the directory containing the pdb files.
    query : str
        Name of the query, default is None.

    Returns
    -------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the directory.
    """

    logger.info(f"Reading {directory}")

    log_dict_list = []

    json_score = os.path.join(directory, f"ranking_debug.json")
    with open(json_score, "r") as f_in:
        json_dict = json.load(f_in)

    model_list = list(json_dict["iptm+ptm"].keys())

    for model in model_list:
        pdb_name = os.path.join(directory, "unrelaxed_" + model + ".pdb")
        assert os.path.isfile(pdb_name)

        token = model.split("_")

        if query is None:
            query_local = os.path.basename(os.path.normpath(directory))
        else:
            query_local = query
        # print(f"query_local: {query_local} {directory}")
        weight = "_".join(token[2:4])
        assert weight in weigths

        json_file = os.path.join(directory, "pae_" + model + ".json")
        assert os.path.isfile(json_file)

        iptm = float(json_dict["iptm"][model])
        conf = float(json_dict["iptm+ptm"][model])
        ptm = (conf - 0.8 * iptm) / 0.2

        local_json_score = os.path.join(directory, "confidence_" + model + ".json")
        with open(local_json_score, "r") as f_in:
            local_json_dict = json.load(f_in)
        pLDDT = np.mean(local_json_dict["confidenceScore"])

        log_dict_list.append(
            {
                "pdb": pdb_name,
                "query": query_local,
                "rank": int(json_dict["order"].index(model) + 1),
                "state": "unrelaxed_",
                "model": int(token[1]),
                "weight": weight,
                "data_file": json_file,
                "ipTM": iptm,
                "ranking_confidence": conf,
                "pTM": ptm,
                "pLDDT": pLDDT,
            }
        )

    log_pd = pd.DataFrame(log_dict_list)
    return log_pd


def read_full_dir(directory):
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

    logger.info(f"Reading full AlphaPulldown {directory}")

    subfolders = [f.path for f in os.scandir(directory) if f.is_dir()]
    log_pd_list = []

    for folder in subfolders:

        # print(f"Reading {folder}")
        log_pd = read_dir(folder)
        # print(log_pd)
        log_pd_list.append(log_pd)

    log_dict = pd.concat(log_pd_list, ignore_index=True)

    return log_dict
