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
                "format",
                "pdb",
                "relaxed_pdb",
                "data_file",
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
    assert (
        my_data.df.columns
        == np.array(
            [
                "pdb",
                "query",
                "model",
                "data_file",
                "chain_iptm",
                "chain_pair_iptm",
                "chain_pair_pae_min",
                "chain_ptm",
                "fraction_disordered",
                "has_clash",
                "ipTM",
                "num_recycles",
                "pTM",
                "ranking_confidence",
                "format",
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

    assert list(my_data.df["ipTM"]) == [0.93, 0.94, 0.93, 0.93, 0.93]


def test_boltz1():
    data_path = os.path.join(TEST_FILE_PATH, "boltz_results_prot_dna_ligand")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "boltz1"
    assert len(my_data.df) == 2
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(
            [
                "pdb",
                "query",
                "model",
                "plddt",
                "data_file",
                "pde",
                "confidence_score",
                "pTM",
                "ipTM",
                "ligand_iptm",
                "protein_iptm",
                "complex_plddt",
                "complex_iplddt",
                "complex_pde",
                "complex_ipde",
                "chains_ptm",
                "pair_chains_iptm",
                "ranking_confidence",
                "format",
            ]
        )
    ).all()

    query = my_data.df.iloc[0]["query"]

    assert my_data.chain_length[query] == [586, 26, 13, 13, 1]
    assert my_data.chains[query] == ["A", "B", "C", "D", "E"]

    # There should be 0 relaxed structures

    assert "relaxed_pdb" not in my_data.df.columns
    print(my_data.df.iloc[:, :])
    assert list(my_data.df["model"]) == list(range(2))
    print(my_data.df["ipTM"])
    expected_iptm = [0.967552, 0.968439]

    precision = 0.01

    np.testing.assert_allclose(
        np.array(list(my_data.df["ipTM"])), np.array(expected_iptm), atol=precision
    )


def test_chai1():
    data_path = os.path.join(TEST_FILE_PATH, "chai1_prot_dna_ligand")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "chai1"
    assert len(my_data.df) == 5
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(
            [
                "pdb",
                "model",
                "ranking_confidence",
                "pTM",
                "ipTM",
                "per_chain_ptm",
                "per_chain_pair_iptm",
                "has_inter_chain_clashes",
                "chain_chain_clashes",
                "query",
                "format",
            ]
        )
    ).all()

    query = my_data.df.iloc[0]["query"]

    assert my_data.chain_length[query] == [586, 26, 13, 13, 1]
    assert my_data.chains[query] == ["A", "B", "C", "D", "E"]

    # There should be 0 relaxed structures

    assert "relaxed_pdb" not in my_data.df.columns
    print(my_data.df.iloc[:, :])
    assert list(my_data.df["model"]) == list(range(5))
    print(my_data.df["ipTM"])
    expected_iptm = [0.954016, 0.953741, 0.953701, 0.953696, 0.953734]

    precision = 0.01

    np.testing.assert_allclose(
        np.array(list(my_data.df["ipTM"])), np.array(expected_iptm), atol=precision
    )
