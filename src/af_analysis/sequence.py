from string import ascii_lowercase
import numpy as np
import logging

# Logging
logger = logging.getLogger(__name__)


def convert_aa_msa(seqs):
    """
    Convert amino acid sequences to numbers.
    """

    convert_dict = {
        "A": 0,
        "R": 1,
        "N": 2,
        "D": 3,
        "C": 4,
        "Q": 5,
        "E": 6,
        "G": 7,
        "H": 8,
        "I": 9,
        "L": 10,
        "K": 11,
        "M": 12,
        "F": 13,
        "P": 14,
        "S": 15,
        "T": 16,
        "W": 17,
        "V": 19,
        "Y": 18,
        "X": 20,
        "-": 21,
    }

    seqnums = []

    for seq in seqs:
        seqnums.append([convert_dict[letter] for letter in seq])

    return np.array(seqnums)


def parse_a3m(a3m_lines=None, a3m_file=None, filter_qid=0.15, filter_cov=0.5, N=100000):
    """
    Parses an A3M file or list of A3M lines and filters sequences based on sequence identity and coverage.

    Parameters
    ----------
    a3m_lines: list of str, optional
        List of lines from an A3M file. Default is None.
    a3m_file: str, optional
        Path to an A3M file. Default is None.
    filter_qid: float, optional
        Minimum sequence identity threshold for filtering. Default is 0.15.
    filter_cov:  float, optional
        Minimum coverage threshold for filtering. Default is 0.5.
    N: int, optional
        Maximum number of sequences to return. Default is 100000.

    Returns
    -------
        tuple: A tuple containing:
            - seqs (list of str): List of filtered sequences.
            - mtx (list of list of int): List of deletion matrices corresponding to the sequences.
            - nams (list of str): List of sequence names.
    """

    def seqid(a, b):
        return sum(c1 == c2 for c1, c2 in zip(a, b))

    def nongaps(a):
        return sum(c != "-" for c in a)

    def chk(seq, ref_seq):
        rL = len(ref_seq)
        L = nongaps(seq)
        return not (L > filter_cov * rL and seqid(seq, ref_seq) > filter_qid * L)

    rm_lower = str.maketrans("", "", ascii_lowercase)

    if a3m_file is None and a3m_lines is None:
        raise ValueError("Either a3m_file or a3m_lines must be provided.")

    # prep inputs
    if a3m_lines is None:
        a3m_lines = open(a3m_file, "r")
    # else: a3m_lines = a3m_lines.splitlines()

    # parse inputs
    n, nams, seqs, mtx = 0, [], [], []

    def do_filter():
        seq = seqs[-1].translate(rm_lower)
        if "_UPI" in nams[-1] or chk(seq, ref_seq):
            nams.pop()
            seqs.pop()
        else:
            # deletion matrix
            deletion_vec = []
            deletion_count = 0
            for j in seqs[-1]:
                if j.islower():
                    deletion_count += 1
                else:
                    deletion_vec.append(deletion_count)
                    deletion_count = 0
            mtx.append(deletion_vec)
            seqs[-1] = seq

    for line in a3m_lines:
        line = line.rstrip()
        if line.startswith(">"):
            if n == 1:
                ref_seq = seqs[0].translate(rm_lower)
            if n >= 1:
                # filter previous entry
                do_filter()
            # start new sequence entry
            nam = line.split()[0][1:]
            if "_" not in nam:
                nam = f"X_{nam}"
            nams.append(nam)
            seqs.append("")
            n += 1
        else:
            seqs[-1] += line

    logger.info(f"- Reading {n:6} sequences.")

    # filter last entry
    do_filter()

    if len(seqs) > N + 1:
        logger.info(
            f"found too many sequences ({len(seqs)}), taking the top{N} (sorted by qid)"
        )
        sid = np.argsort([seqid(seq, ref_seq) for seq in seqs])[::-1][: N + 1]
        seqs = [seqs[i] for i in sid]
        mtx = [mtx[i] for i in sid]
        nams = [nams[i] for i in sid]
    # return seqs[1:], mtx[1:], nams[1:]
    return seqs, mtx, nams
