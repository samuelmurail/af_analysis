#!/usr/bin/env python3
# coding: utf-8

import os
import pytest
import numpy as np
import af_analysis
from af_analysis import clustering

from .data_files import TEST_FILE_PATH
from unittest.mock import patch


@patch("matplotlib.pyplot.show")
def test_cf_1_5_5_relax_new_clust(mock_show):
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af_analysis.Data(data_path)

    assert "cluster" not in my_data.df.columns
    clustering.hierarchical(my_data.df, threshold=0.3, MDS_coors=False, rmsd_scale=True)
    assert "cluster" in my_data.df.columns

    clust_num = len(my_data.df["cluster"].unique())
    assert clust_num == 6
    assert list(my_data.df["cluster"]) == [
        3,
        2,
        1,
        1,
        1,
        4,
        2,
        1,
        5,
        6,
        3,
        2,
        1,
        1,
        1,
        2,
        2,
        1,
        1,
        1,
        3,
        2,
        1,
        1,
        1,
        3,
        2,
        1,
        5,
        6,
        4,
        2,
        1,
        1,
        1,
        3,
        2,
        1,
        1,
        1,
    ]

    clustering.hierarchical(
        my_data.df,
        MDS_coors=False,
    )
    clust_num = len(my_data.df["cluster"].unique())
    assert clust_num == 18
    print([i for i in my_data.df["cluster"]])
    assert list(my_data.df["cluster"]) == [
        6,
        7,
        1,
        9,
        2,
        10,
        4,
        11,
        12,
        8,
        6,
        13,
        1,
        3,
        2,
        14,
        15,
        1,
        3,
        2,
        5,
        4,
        1,
        3,
        2,
        5,
        7,
        1,
        16,
        8,
        17,
        18,
        1,
        3,
        2,
        5,
        4,
        1,
        3,
        2,
    ]

    assert "MDS 1" not in my_data.df.columns
    assert "MDS 2" not in my_data.df.columns
    assert "MDS 3" not in my_data.df.columns
    clustering.hierarchical(
        my_data.df,
    )
    assert "MDS 1" in my_data.df.columns
    assert "MDS 2" in my_data.df.columns

    # Can't check results of MDS 1 and MDS 2 as they seem to be random

    clustering.hierarchical(
        my_data.df,
        MDS_coors=False,
        align_selection="backbone and chainID A",
        distance_selection="backbone and chainID B",
    )
    clust_num = len(my_data.df["cluster"].unique())
    assert clust_num == 23
