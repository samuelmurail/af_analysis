#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np

import af_analysis
from .data_files import TEST_FILE_PATH


def test_cf_1_5_5_relax():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "colabfold_1.5"
    assert len(my_data.df) == 40
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(
            [
                "query",
                "seed",
                "model",
                "weight",
                "recycle",
                "pLDDT",
                "pTM",
                "ipTM",
                "ranking_confidence",
                "pdb",
                "relaxed_pdb",
                "json",
            ]
        )
    ).all()

    query = my_data.df.iloc[0]["query"]

    assert my_data.chain_length[query] == [42, 42]
    assert my_data.chains[query] == ["A", "B"]

    # There should be only 5 relaxed structures

    relaxed_num = sum(my_data.df["relaxed_pdb"].notna())

    assert relaxed_num == 5
    assert list(my_data.df["recycle"]) == [
        9,
        16,
        3,
        5,
        5,
        14,
        15,
        15,
        48,
        12,
        7,
        14,
        7,
        8,
        9,
        4,
        18,
        6,
        10,
        8,
        11,
        14,
        7,
        11,
        9,
        11,
        7,
        7,
        19,
        9,
        30,
        7,
        5,
        9,
        7,
        12,
        7,
        6,
        6,
        6,
    ]

    assert list(my_data.df["ipTM"]) == [
        0.0812,
        0.0685,
        0.161,
        0.158,
        0.541,
        0.117,
        0.0698,
        0.239,
        0.0648,
        0.331,
        0.0789,
        0.0815,
        0.145,
        0.306,
        0.604,
        0.0997,
        0.0662,
        0.143,
        0.219,
        0.589,
        0.0794,
        0.0684,
        0.15,
        0.299,
        0.559,
        0.0797,
        0.0662,
        0.147,
        0.0609,
        0.318,
        0.0964,
        0.0683,
        0.151,
        0.274,
        0.584,
        0.0776,
        0.0693,
        0.14,
        0.199,
        0.598,
    ]


def test_af3_webserver():
    data_path = os.path.join(TEST_FILE_PATH, "fold_2024_07_01_12_14_prot_dna_zn")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "AF3_webserver"
    assert len(my_data.df) == 5
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(
            [
                "pdb",
                "query",
                "model",
                "json",
                "chain_iptm",
                "chain_pair_iptm",
                "chain_pair_pae_min",
                "chain_ptm",
                "fraction_disordered",
                "has_clash",
                "iptm",
                "num_recycles",
                "ptm",
                "ranking_score",
            ]
        )
    ).all()

    query = my_data.df.iloc[0]["query"]

    assert my_data.chain_length[query] == [90, 1, 1, 1, 11, 11]
    assert my_data.chains[query] == ["A", "B", "C", "D", "E", "F"]

    # There should be 0 relaxed structures

    assert "relaxed_pdb" not in my_data.df.columns
    print(my_data.df.iloc[:, :])
    assert list(my_data.df["num_recycles"]) == [10] * 5

    assert list(my_data.df["iptm"]) == [0.93, 0.94, 0.93, 0.93, 0.93]
