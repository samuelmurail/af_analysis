#!/usr/bin/env python3
# coding: utf-8

import os
import re
import pathlib
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
    "multimer_v1",
    "multimer_v2",
    "multimer_v3",
]


def read_dir(directory):
    """Extract pdb list from a directory."""

    logger.info(f"Reading {directory}")

    log_dict_list = []

    for file in os.listdir(directory):
        if file.endswith(".pdb"):
            token = file.split("_")
            if token[0] == "ranked":
                continue

            for i, t in enumerate(token):
                if t in ["relaxed", "unrelaxed"]:
                    start_index = i
                    # print(start_index)
                    break

            if start_index == 0:
                path_dir = pathlib.Path(directory)
                name_token = (
                    path_dir.parent.name if not path_dir.is_dir() else path_dir.name
                )
                state = token[start_index]
                rank = None
                weight = "_".join(token[start_index + 3 : -2])
                assert weight in weigths
                model = int(token[start_index + 2])
                seed = int(token[-1][:-4])
            else:
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


def add_json(log_pd, directory, verbose=True):
    """Find json files in the directory.

    Parameters
    ----------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    directory : str
        Path to the directory containing the json files.
    verbose : bool
        If True, show a progress bar.
        If False, no progress bar is shown.

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    logger.info(f"Extracting json files location")

    raw_list = os.listdir(directory)
    file_list = []
    for file in raw_list:
        if file.endswith(".json") or file.endswith(".pkl"):
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

    disable = False if verbose else True

    for i, last in tqdm(
        enumerate(last_recycle), total=len(last_recycle), disable=disable
    ):
        row = log_pd.iloc[i]

        reg = rf"{row['query']}_scores_.*_{row['weight']}_model_{row['model']}_seed_{row['seed']:03d}\.json"
        r = re.compile(reg)
        res = list(filter(r.match, file_list))

        reg_pkl = (
            rf"result_model_{row['model']}_{row['weight']}_pred_{row['seed']}\.pkl"
        )
        r_pkl = re.compile(reg_pkl)
        res_pkl = list(filter(r_pkl.match, file_list))

        if len(res) == 1:
            json_list.append(os.path.join(directory, res[0]))
        elif len(res) == 0:
            if len(res_pkl) == 1:
                json_list.append(os.path.join(directory, res_pkl[0]))
            else:
                logger.warning(f"Not founded : {reg} or {reg_pkl}")
                json_list.append(None)
        else:
            logger.warning(f"Multiple json file for {reg}: {res}")
            json_list.append(res[0])

    log_pd.loc[:, "data_file"] = json_list
