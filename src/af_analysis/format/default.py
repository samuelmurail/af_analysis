#!/usr/bin/env python3
# coding: utf-8

import os
import re
import logging
from tqdm.auto import tqdm
import pandas as pd

logger = logging.getLogger()
logger.setLevel(logging.INFO)

weigths = [
    "alphafold2",
    "alphafold2_ptm",
    "alphafold2_multimer_v1",
    "alphafold2_multimer_v2",
    "alphafold2_multimer_v3",
]


def read_dir(directory):
    """Extract pdb list from a directory."""

    logger.info(f"Reading {directory}")

    log_dict_list = []

    for file in os.listdir(directory):
        if file.endswith(".pdb"):
            token = file.split("_")

            for i, t in enumerate(token):
                if t in ["relaxed", "unrelaxed"]:
                    start_index = i
                    # print(start_index)
                    break

            name_token = "_".join(token[:start_index])
            state = token[start_index]
            rank = int(token[start_index + 2])
            weight = "_".join(token[start_index + 3 : -4])
            assert weight in weigths
            model = int(token[-3])
            seed = int(token[-1][:-4])

            log_dict_list.append(
                {
                    "pdb": os.path.join(directory, file),
                    "query": name_token,
                    "rank": rank,
                    "state": state,
                    "seed": seed,
                    "model": model,
                    "weight": weight,
                }
            )

    log_pd = pd.DataFrame(log_dict_list)
    return log_pd


def add_json(log_pd, directory):
    """Find json files in the directory.

    Parameters
    ----------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    directory : str
        Path to the directory containing the json files.

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    logger.info(f"Extracting json files location")

    raw_list = os.listdir(directory)
    file_list = []
    for file in raw_list:
        if file.endswith(".json"):
            file_list.append(file)

    json_list = []

    # Get the last recycle:
    if "recycle" not in log_pd.columns:
        last_recycle = log_pd.groupby(["query", "seed", "model", "weight", "state"])
    else:
        last_recycle = (
            log_pd.groupby(["query", "seed", "model", "weight"])["recycle"].transform(
                "max"
            )
            == log_pd["recycle"]
        )

    for i, last in tqdm(enumerate(last_recycle), total=len(last_recycle)):
        row = log_pd.iloc[i]

        reg = rf"{row['query']}_scores_.*_{row['weight']}_model_{row['model']}_seed_{row['seed']:03d}\.json"
        r = re.compile(reg)

        res = list(filter(r.match, file_list))

        if len(res) == 1:
            json_list.append(os.path.join(directory, res[0]))
        elif len(res) == 0:
            logger.warning(f"Not founded : {reg}")
            json_list.append(None)
        else:
            logger.warning(f"Multiple json file for {reg}: {res}")
            json_list.append(res[0])

    log_pd.loc[:, "data_file"] = json_list
