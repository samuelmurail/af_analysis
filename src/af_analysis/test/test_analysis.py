#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest

import af_analysis
from af_analysis import analysis
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

    analysis.pdockq(my_data)
    print([round(i, 4) for i in my_data.df["pdockq"]])
    expected_pdockq = [
        0.0332,
        0.0205,
        0.03,
        0.0297,
        0.1209,
        0.0285,
        0.0225,
        0.0485,
        0.0198,
        0.0715,
        0.0332,
        0.0238,
        0.0276,
        0.0558,
        0.1383,
        0.0242,
        0.0211,
        0.0252,
        0.0419,
        0.1415,
        0.0295,
        0.0212,
        0.0285,
        0.0524,
        0.137,
        0.0291,
        0.0204,
        0.0284,
        0.0207,
        0.0823,
        0.0231,
        0.0203,
        0.0282,
        0.0509,
        0.1392,
        0.0294,
        0.0206,
        0.0254,
        0.0362,
        0.1426,
    ]

    precision = 0.01
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq"] == pytest.approx(expected_pdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )


def test_af3_webserver():
    data_path = os.path.join(TEST_FILE_PATH, "fold_2024_07_01_12_14_prot_dna_zn")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "AF3_webserver"

    analysis.pdockq(my_data)

    expected_pdockq = [0.2756, 0.2621, 0.2755, 0.2754, 0.2758]

    # print([round(i, 4) for i in my_data.df["pdockq"]])
    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq"] == pytest.approx(expected_pdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.mpdockq(my_data)
    expected_mpdockq = [0.262, 0.262, 0.262, 0.262, 0.262]

    # print([round(i, 6) for i in my_data.df["mpdockq"]])
    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["mpdockq"]
            == pytest.approx(expected_mpdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.pdockq2(my_data)
    # print([round(i, 4) for i in my_data.df["pdockq2_A"]])
    expected_pdockq2 = [0.9148, 0.9187, 0.9151, 0.913, 0.9154]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_A"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )
    # print([round(i, 4) for i in my_data.df["pdockq2_D"]])
    expected_pdockq2 = [0.8972, 0.8925, 0.8884, 0.889, 0.8785]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_D"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.LIS_matrix(my_data)
    expected_LIS_0 = [
        [0.83139, 0.8075, 0.85381, 0.85251, 0.85559, 0.85551],
        [0.82717, 0.93666, 0.7975, 0.83166, 0.84568, 0.83962],
        [0.80911, 0.745, 0.93666, 0.82333, 0.83409, 0.83060],
        [0.84268, 0.84166, 0.84, 0.93666, 0.865, 0.85886],
        [0.83427, 0.84507, 0.83712, 0.84522, 0.87633, 0.87050],
        [0.81519, 0.79833, 0.83053, 0.82, 0.85831, 0.86331],
    ]

    np.testing.assert_allclose(
        np.array(my_data.df["LIS"][0]), np.array(expected_LIS_0), atol=precision
    )

    analysis.inter_chain_pae(my_data)

    expected_PAE_A_B = [2.8373, 2.6611, 2.8013, 2.8286, 2.7292]
    # print([round(i, 4) for i in my_data.df["PAE_A_B"]])
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_B"]
            == pytest.approx(expected_PAE_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    expected_PAE_A_E = [2.7772, 2.6177, 2.8398, 2.8672, 2.7849]
    # print([round(i, 4) for i in my_data.df["PAE_A_E"]])
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_E"]
            == pytest.approx(expected_PAE_A_E[i], precision)
            for i in range(len(my_data.df))
        ]
    )


def test_af3_boltz1():
    data_path = os.path.join(TEST_FILE_PATH, "boltz_results_prot_dna_ligand")

    my_data = af_analysis.Data(data_path)

    assert my_data.format == "boltz1"

    analysis.pdockq(my_data)

    expected_pdockq = [0.018281, 0.018281]

    # print([round(i, 6) for i in my_data.df["pdockq"]])
    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq"] == pytest.approx(expected_pdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.mpdockq(my_data)
    expected_mpdockq = [0.262, 0.262]

    # print([round(i, 6) for i in my_data.df["mpdockq"]])
    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["mpdockq"]
            == pytest.approx(expected_mpdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.pdockq2(my_data)
    # print([round(i, 6) for i in my_data.df["pdockq2_A"]])
    expected_pdockq2 = [0.007527, 0.007527]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_A"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    # print([round(i, 6) for i in my_data.df["pdockq2_D"]])
    expected_pdockq2 = [0.007526, 0.007526]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_D"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.LIS_matrix(my_data)
    expected_LIS_0 = [
        [0.781716, 0.77767, 0.78187, 0.774672, 0.833509],
        [0.72203, 0.769744, 0.767006, 0.767054, 0.753537],
        [0.711757, 0.757241, 0.857382, 0.654726, 0.741389],
        [0.707814, 0.750932, 0.662734, 0.854446, 0.770577],
        [0.280363, 0.323815, 0.455396, 0.289744, 0.978317],
    ]
    for j in range(len(my_data.df["LIS"][0])):
        print([round(i, 6) for i in my_data.df["LIS"][0][j]])

    np.testing.assert_allclose(
        np.array(my_data.df["LIS"][0]), np.array(expected_LIS_0), atol=precision
    )

    analysis.inter_chain_pae(my_data)

    expected_PAE_A_B = [2.893401, 2.936897]
    # print([round(i, 6) for i in my_data.df["PAE_A_B"]])
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_B"]
            == pytest.approx(expected_PAE_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    # print([round(i, 6) for i in my_data.df["PAE_A_E"]])
    expected_PAE_A_E = [2.181038, 2.197674]
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_E"]
            == pytest.approx(expected_PAE_A_E[i], precision)
            for i in range(len(my_data.df))
        ]
    )
