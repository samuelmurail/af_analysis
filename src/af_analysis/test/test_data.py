#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pandas as pd
import pytest

import af_analysis
from af_analysis import analysis
from .data_files import TEST_FILE_PATH


def test_cf_1_5_5_write_read_csv(tmp_path):
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)
    my_data.export_csv(os.path.join(tmp_path, "test.csv"))

    my_data2 = af_analysis.Data(csv=os.path.join(tmp_path, "test.csv"))

    assert my_data2.format == "csv"
    assert len(my_data2.df) == 40
    print(my_data2.df.columns)
    assert (
        my_data2.df.columns
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
                "format",
                "pdb",
                "relaxed_pdb",
                "data_file",
            ]
        )
    ).all()

    query = my_data2.df.iloc[0]["query"]

    assert my_data2.chain_length[query] == [42, 42]
    assert my_data2.chains[query] == ["A", "B"]

    # There should be only 5 relaxed structures

    relaxed_num = sum(my_data2.df["relaxed_pdb"].notna())

    assert relaxed_num == 5
    assert list(my_data2.df["recycle"]) == [
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

    assert list(my_data2.df["ipTM"]) == [
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


def test_get_plddt():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)

    plddt_array = my_data.get_plddt(0)

    print(plddt_array)

    assert len(plddt_array) == 84

    expected_plddt = [
        36.12,
        39.41,
        41.03,
        42.34,
        41.69,
        43.5,
        42.78,
        44.25,
        38.16,
        39.84,
        37.75,
        36.56,
        36.06,
        39.0,
        38.66,
        38.44,
        39.34,
        39.06,
        41.03,
        37.22,
        39.31,
        38.75,
        40.72,
        43.12,
        39.62,
        41.19,
        39.22,
        39.41,
        39.31,
        37.56,
        39.06,
        41.12,
        40.12,
        37.81,
        36.56,
        36.75,
        37.97,
        38.56,
        34.19,
        32.66,
        33.19,
        33.59,
        36.28,
        38.75,
        40.5,
        42.44,
        41.28,
        42.72,
        42.72,
        43.91,
        37.81,
        40.03,
        37.84,
        36.28,
        35.72,
        38.91,
        38.38,
        38.62,
        39.19,
        39.03,
        40.72,
        37.31,
        39.41,
        38.56,
        40.31,
        42.84,
        39.75,
        41.06,
        39.25,
        39.16,
        39.34,
        37.66,
        38.91,
        40.66,
        40.47,
        37.78,
        36.41,
        36.5,
        37.69,
        38.25,
        34.09,
        32.03,
        33.06,
        33.16,
    ]
    precision = 0.001
    assert plddt_array == pytest.approx(expected_plddt, precision)


def test_get_plddt_dna_ions():
    data_path = os.path.join(TEST_FILE_PATH, "fold_2024_07_01_12_14_prot_dna_zn")
    my_data = af_analysis.Data(data_path)

    plddt_array = my_data.get_plddt(0)

    print(plddt_array)

    assert len(plddt_array) == 115

    expected_plddt = [
        83.82,
        89.41,
        97.62,
        98.41,
        98.81,
        98.88,
        98.9,
        98.82,
        98.7,
        97.65,
        98.17,
        98.85,
        98.46,
        98.89,
        98.86,
        98.94,
        98.85,
        98.84,
        98.84,
        98.9,
        98.91,
        98.92,
        98.85,
        98.96,
        98.98,
        98.99,
        98.88,
        98.71,
        98.8,
        98.8,
        98.98,
        98.94,
        98.99,
        98.83,
        98.98,
        98.99,
        98.91,
        98.97,
        98.9,
        98.95,
        98.99,
        98.91,
        98.98,
        98.92,
        98.99,
        98.93,
        98.94,
        98.94,
        98.95,
        98.92,
        98.9,
        98.98,
        98.92,
        98.97,
        98.92,
        98.98,
        98.98,
        98.91,
        98.96,
        98.93,
        98.96,
        98.93,
        98.94,
        98.93,
        98.94,
        98.82,
        98.78,
        98.97,
        98.98,
        98.92,
        98.99,
        98.94,
        98.9,
        98.96,
        98.92,
        98.9,
        98.98,
        98.95,
        98.9,
        98.94,
        98.84,
        98.78,
        98.41,
        98.15,
        98.4,
        96.75,
        86.69,
        82.44,
        81.47,
        78.44,
        98.9,
        98.94,
        98.98,
        93.57,
        98.57,
        98.89,
        98.94,
        98.91,
        98.87,
        98.86,
        98.83,
        98.75,
        98.71,
        98.02,
        81.79,
        94.31,
        98.18,
        98.67,
        98.84,
        98.81,
        98.73,
        98.57,
        98.55,
        98.43,
        98.27,
    ]
    precision = 0.001
    assert plddt_array == pytest.approx(expected_plddt, precision)


def test_concat():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)

    data_path_2 = os.path.join(TEST_FILE_PATH, "fold_2024_07_01_12_14_prot_dna_zn")
    my_data_2 = af_analysis.Data(data_path_2)

    my_data_all = af_analysis.data.concat_data([my_data, my_data_2])
    assert len(my_data_all.df) == 45
    assert len(my_data_all.chain_length) == 2
    assert len(my_data_all.chains) == 2


def test_extract_data():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)
    print(my_data.df.columns)
    new_column = ["plddt", "max_pae", "pae", "ptm", "iptm"]
    assert np.all([not (col in my_data.df.columns) for col in new_column])
    my_data.extract_data()
    print(my_data.df.columns)
    assert np.all([col in my_data.df.columns for col in new_column])

    old_column = ["pLDDT", "pTM", "ipTM"]

    rtol = 1e-02
    assert np.all(np.abs(my_data.df["pTM"] - my_data.df["ptm"]) < rtol)
    assert np.all(np.abs(my_data.df["ipTM"] - my_data.df["iptm"]) < rtol)
    rtol = 1e-01
    assert np.all(
        np.abs(
            np.array([np.mean(plddts) for plddts in my_data.df["plddt"]])
            - my_data.df["pLDDT"]
        )
        < rtol
    )


def test_read_data_dict():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    pdb_list = [
        "beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002.pdb",
        "beta_amyloid_dimer_d2fa3_0_relaxed_rank_002_alphafold2_multimer_v3_model_5_seed_007.pdb",
        "beta_amyloid_dimer_d2fa3_0_relaxed_rank_003_alphafold2_multimer_v3_model_5_seed_003.pdb",
        "beta_amyloid_dimer_d2fa3_0_relaxed_rank_004_alphafold2_multimer_v3_model_5_seed_006.pdb",
        "beta_amyloid_dimer_d2fa3_0_relaxed_rank_005_alphafold2_multimer_v3_model_5_seed_004.pdb",
    ]
    pdb_list = [os.path.join(data_path, pdb) for pdb in pdb_list]
    json_list = [
        "beta_amyloid_dimer_d2fa3_0_scores_rank_001_alphafold2_multimer_v3_model_5_seed_002.json",
        "beta_amyloid_dimer_d2fa3_0_scores_rank_002_alphafold2_multimer_v3_model_5_seed_007.json",
        "beta_amyloid_dimer_d2fa3_0_scores_rank_003_alphafold2_multimer_v3_model_5_seed_003.json",
        "beta_amyloid_dimer_d2fa3_0_scores_rank_004_alphafold2_multimer_v3_model_5_seed_006.json",
        "beta_amyloid_dimer_d2fa3_0_scores_rank_005_alphafold2_multimer_v3_model_5_seed_004.json",
    ]
    json_list = [os.path.join(data_path, json) for json in json_list]

    data_dict = {
        "data_file": json_list,
        "pdb": pdb_list,
        "query": "beta_amyloid_dimer_d2fa3_0",
    }

    my_data = af_analysis.Data(data_dict=data_dict)
    my_data.extract_fields(["plddt", "ptm", "iptm"])

    assert len(my_data.df) == 5
    assert my_data.chain_length == {"beta_amyloid_dimer_d2fa3_0": [42, 42]}
    assert my_data.chains == {"beta_amyloid_dimer_d2fa3_0": ["A", "B"]}
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(["data_file", "pdb", "query", "format", "plddt", "ptm", "iptm"])
    ).all()
    print(my_data.df["iptm"])
    assert list(my_data.df["iptm"]) == [0.6, 0.6, 0.59, 0.58, 0.56]
