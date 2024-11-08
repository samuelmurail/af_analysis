#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest
import random

import af_analysis
from af_analysis import analysis
from .data_files import TEST_FILE_PATH
from af_analysis.sequence import parse_a3m


def test_parse_a3m_with_lines():
    a3m_lines = [
        ">seq1",
        "ACDEFGHIKLMNPQRSTVWY",
        ">seq2",
        "ACDEFGHIKLMNPQRSTVWY",
        ">seq3",
        "ACDEFGHIKLMNPQRSTVWY",
    ]
    seqs, mtx, nams = parse_a3m(a3m_lines=a3m_lines)

    assert len(seqs) == 3
    assert len(mtx) == 3
    assert len(nams) == 3
    assert seqs[0] == "ACDEFGHIKLMNPQRSTVWY"
    assert nams[0] == "X_seq1"


def test_parse_a3m_with_file(tmp_path):
    a3m_content = ">seq1\nACDEFGHIKLMNPQRSTVWY\n>seq2\nACDEFGHIKLMNPQRSTVWY\n>seq3\nACDEFGHIKLMNPQRSTVWY\n"
    a3m_file = tmp_path / "test.a3m"
    a3m_file.write_text(a3m_content)

    seqs, mtx, nams = parse_a3m(a3m_file=str(a3m_file))

    assert len(seqs) == 3
    assert len(mtx) == 3
    assert len(nams) == 3
    assert seqs[0] == "ACDEFGHIKLMNPQRSTVWY"
    assert nams[0] == "X_seq1"


def test_parse_a3m_filtering():
    a3m_lines = [
        ">seq1",
        "ACDEFGHIKLMNPQRSTVWY",
        ">seq2",
        "ACDEFGHIKLMNPQRSTVWY",
        ">seq3",
        "ACDEFGHIKLMNPQRSTVWY",
        ">seq4",
        "acdefghiklmnpqrstvwy",  # Lowercase should be filtered out
    ]
    seqs, mtx, nams = parse_a3m(a3m_lines=a3m_lines)

    assert len(seqs) == 3
    assert len(mtx) == 3
    assert len(nams) == 3
    assert "X_seq4" not in nams


def test_parse_a3m_file():
    a3m_file = os.path.join(
        TEST_FILE_PATH,
        "beta_amyloid_dimer_cf_1.5.5/beta_amyloid_dimer_d2fa3_0_env/uniref.a3m",
    )

    seqs, mtx, nams = parse_a3m(a3m_file=a3m_file, filter_qid=0)

    assert len(seqs) == 271
    assert len(mtx) == 271
    assert len(nams) == 271

    seqs, mtx, nams = parse_a3m(a3m_file=a3m_file, filter_qid=0, N=100)

    assert len(seqs) == 101
    assert len(mtx) == 101
    assert len(nams) == 101

    seqs, mtx, nams = parse_a3m(a3m_file=a3m_file, filter_qid=0.8)

    assert len(seqs) == 239
    assert len(mtx) == 239
    assert len(nams) == 239

    seqs, mtx, nams = parse_a3m(a3m_file=a3m_file, filter_qid=0.8, filter_cov=0.8)

    assert len(seqs) == 198
    assert len(mtx) == 198
    assert len(nams) == 198
