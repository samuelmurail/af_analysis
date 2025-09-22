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
    expected_pdockq2 = [0.9599, 0.9614, 0.9606, 0.959, 0.9607]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_A"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )
    
    # print([round(i, 4) for i in my_data.df["pdockq2_D"]])
    expected_pdockq2 = [0.9632, 0.9646, 0.9642, 0.9638, 0.964]
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


def test_cf_1_5_5_ftdmp():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af_analysis.Data(data_path)

    ftdmp_path = os.path.join(TEST_FILE_PATH, "ftdmp_beta_amyloid_dimer")

    ftdmp_df_list = analysis.extract_ftdmp(ftdmp_path)
    assert len(ftdmp_df_list) == 1

    my_data.df["ID"] = [
        os.path.basename(file_path) if file_path is not None else None
        for file_path in my_data.df["pdb"]]
    my_data.df = my_data.df.merge(ftdmp_df_list [0], on="ID", how="inner")


    print([round(i, 4) for i in my_data.df["raw_FGV_full_light_score"]])
    expected_FGV_full_light_score = [
        0.0992,
        0.0838,
        0.2851,
        0.258,
        0.3597,
        0.1517,
        0.0859,
        0.3844,
        0.1133,
        0.2971,
        0.1068,
        0.0879,
        0.2687,
        0.376,
        0.3755,
        0.1128,
        0.0894,
        0.2679,
        0.3602,
        0.3529,
        0.0923,
        0.0853,
        0.2877,
        0.41,
        0.3564,
        0.098,
        0.0842,
        0.2678,
        0.1171,
        0.3017,
        0.1077,
        0.0874,
        0.2723,
        0.3824,
        0.3651,
        0.1024,
        0.0857,
        0.2716,
        0.348,
        0.3926,
    ]

    precision = 0.01
    assert np.all(
        [
            my_data.df.iloc[i]["raw_FGV_full_light_score"]
            == pytest.approx(expected_FGV_full_light_score[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    # raw_FIGNN_average_gnn_score
    print([round(i, 4) for i in my_data.df["raw_FIGNN_average_gnn_score"]])
    expected_FIGNN_average_gnn_score = [
        -1.2378,
        -0.3711,
        -0.9664,
        -1.4222,
        1.0206,
        1.9895,
        -0.8248,
        1.5024,
        -0.2968,
        -0.0947,
        -1.5941,
        -0.1801,
        0.0933,
        2.1569,
        0.7643,
        0.0444,
        -1.3547,
        -0.8445,
        1.0281,
        0.7819,
        -1.6051,
        -0.7928,
        -0.8217,
        1.9333,
        0.3474,
        -1.5848,
        -0.6718,
        -0.7082,
        -1.2727,
        -0.6135,
        0.1668,
        -1.1248,
        -0.6884,
        1.6645,
        0.4051,
        -1.736,
        -1.002,
        -0.809,
        1.3703,
        0.5616,
    ]

    precision = 0.01
    assert np.all(
        [
            my_data.df.iloc[i]["raw_FIGNN_average_gnn_score"]
            == pytest.approx(expected_FIGNN_average_gnn_score[i], precision)
            for i in range(len(my_data.df))
        ]
    )


def test_iptm_d0():
    """Test for ipTM and D0 calculation from colabfold 1.5.5 data.

    ``` python
    python ~/Documents/Code/IPSAE/ipsae.py src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_scores_rank_001_alphafold2_multimer_v3_model_5_seed_002.json src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002.pdb 10 10
    ```

    The output should be similar to the following:
    ```
    Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model
    A    B     10   10   asym  0.309192    0.529505    0.492353    0.600    0.518804      0.1443     0.1843     0.4630      41     84     73    1.87    3.29    3.00     31      42      27      27   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    B    A     10   10   asym  0.309226    0.529211    0.503041    0.600    0.518532      0.1443     0.1839     0.4598      41     84     76    1.87    3.29    3.08     34      42      27      26   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    A    B     10   10   max   0.309226    0.529505    0.503041    0.600    0.518804      0.1443     0.1843     0.4614      41     84     76    1.87    3.29    3.08     42      42      27      27   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    ```

    """

    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)

    analysis.ipTM_d0(my_data)

    print([round(i, 4) for i in my_data.df["ipTM_d0_A_B"]])
    expected_ipTM_d0_A_B = [
        0.0427,
        0.039,
        0.0555,
        0.0611,
        0.442,
        0.0625,
        0.0396,
        0.0891,
        0.0331,
        0.1421,
        0.0418,
        0.0443,
        0.0518,
        0.1504,
        0.5188,  # The one tested previously
        0.0539,
        0.037,
        0.0528,
        0.0883,
        0.5003,
        0.0434,
        0.0391,
        0.0536,
        0.1346,
        0.4595,
        0.0435,
        0.0375,
        0.0537,
        0.0332,
        0.1323,
        0.0515,
        0.0384,
        0.0542,
        0.1222,
        0.4938,
        0.0419,
        0.0388,
        0.0538,
        0.0774,
        0.5212,
    ]

    precision = 0.01
    assert np.all(
        [
            my_data.df.iloc[i]["ipTM_d0_A_B"]
            == pytest.approx(expected_ipTM_d0_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )


def test_ipSAE():
    """Test for ipSAE calculation from colabfold 1.5.5 data.

    ``` python
    python ~/Documents/Code/IPSAE/ipsae.py src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_scores_rank_001_alphafold2_multimer_v3_model_5_seed_002.json src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002.pdb 10 10
    ```

    The output should be similar to the following:
    ```
    Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model
    A    B     10   10   asym  0.309192    0.529505    0.492353    0.600    0.518804      0.1443     0.1843     0.4630      41     84     73    1.87    3.29    3.00     31      42      27      27   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    B    A     10   10   asym  0.309226    0.529211    0.503041    0.600    0.518532      0.1443     0.1839     0.4598      41     84     76    1.87    3.29    3.08     34      42      27      26   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    A    B     10   10   max   0.309226    0.529505    0.503041    0.600    0.518804      0.1443     0.1843     0.4614      41     84     76    1.87    3.29    3.08     42      42      27      27   src/af_analysis/test/inputs/beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_relaxed_rank_001_alphafold2_multimer_v3_model_5_seed_002
    ```

    """

    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af_analysis.Data(data_path)

    analysis.ipSAE(my_data)

    print([round(i, 4) for i in my_data.df["ipSAE_A_B"]])
    expected_ipSAE_A_B = [
        0.0156,
        0.0133,
        0.0158,
        0.0161,
        0.2411,
        0.0192,
        0.0141,
        0.0183,
        0.0,
        0.0314,
        0.0147,
        0.0138,
        0.0133,
        0.0289,
        0.3092, # The one tested previously
        0.0162,
        0.0125,
        0.0137,
        0.0169,
        0.2949,
        0.0153,
        0.0137,
        0.0141,
        0.0224,
        0.2559,
        0.0154,
        0.0127,
        0.015,
        0.0,
        0.0236,
        0.0153,
        0.0132,
        0.0144,
        0.0231,
        0.2902,
        0.0145,
        0.0135,
        0.0141,
        0.0183,
        0.3137,
    ]

    precision = 0.01
    assert np.all(
        [
            my_data.df.iloc[i]["ipSAE_A_B"]
            == pytest.approx(expected_ipSAE_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )
