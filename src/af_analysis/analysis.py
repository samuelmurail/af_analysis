import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import json

import pdb_numpy

# Autorship information
__author__ = "Alaa Reguei"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License version 2"
__version__ = "0.1.3"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"


def get_pae(data_file):
    """Get the PAE matrix from a json/npz file.

    Parameters
    ----------
    data_file : str
        Path to the json/npz file.

    Returns
    -------
    np.array
        PAE matrix.
    """

    if data_file is None:
        return None

    if data_file.endswith(".json"):
        return extract_pae_json(data_file)
    elif data_file.endswith(".npz"):
        return extract_pae_npz(data_file)
    else:
        raise ValueError("Unknown file format.")


def extract_pae_json(json_file):
    """Get the PAE matrix from a json file.

    Parameters
    ----------
    json_file : str
        Path to the json file.

    Returns
    -------
    np.array
        PAE matrix.
    """

    with open(json_file) as f:
        local_json = json.load(f)

    if "pae" in local_json:
        pae_array = np.array(local_json["pae"])
    elif "predicted_aligned_error" in local_json[0]:
        pae_array = np.array(local_json[0]["predicted_aligned_error"])
    else:
        raise ValueError("No PAE found in the json file.")

    return pae_array


def extract_pae_npz(npz_file):
    """Get the PAE matrix from a json file.

    Parameters
    ----------
    npz_file : str
        Path to the npz file.

    Returns
    -------
    np.array
        PAE matrix.
    """

    data_npz = np.load(npz_file)
    pae_array = data_npz["pae"]

    return pae_array


def extract_fields_file(data_file, fields):
    """Get the PAE matrix from a json/pickle file.

    Parameters
    ----------
    file : str
        Path to the json file.
    fields : list
        List of fields to extract.

    Returns
    -------
    value
    """

    if data_file is None:
        return None

    if data_file.endswith(".json"):
        with open(data_file) as f:
            local_data = json.load(f)
    elif data_file.endswith(".npz"):
        local_data = np.load(data_file)

    values = []
    for field in fields:
        if field in local_data:
            values.append(local_data[field])
        else:
            raise ValueError(f"No field {field} found in the json/npz file.")

    return values


def pdockq(data, verbose=True):
    r"""Compute the pDockq [1]_ from the pdb file.

    .. math::
        pDockQ = \frac{L}{1 + e^{-k (x-x_{0})}} + b

    where:

    .. math::
        x = \overline{plDDT_{interface}} \cdot log(number \: of \: interface \: contacts)

    :math:`L = 0.724` is the maximum value of the sigmoid,
    :math:`k = 0.052` is the slope of the sigmoid, :math:`x_{0} = 152.611`
    is the midpoint of the sigmoid, and :math:`b = 0.018` is the y-intercept
    of the sigmoid.

    Implementation was inspired from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py

    Parameters
    ----------
    data : AFData
        object containing the data
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.


    References
    ----------

    .. [1] Bryant P, Pozzati G and Elofsson A. Improved prediction of protein-protein
        interactions using AlphaFold2. *Nature Communications*. vol. 13, 1265 (2022)
        https://www.nature.com/articles/s41467-022-28865-w


    """

    from pdb_numpy.analysis import compute_pdockQ

    pdockq_list = []

    disable = False if verbose else True

    for pdb in tqdm(data.df["pdb"], total=len(data.df["pdb"]), disable=disable):
        if pdb is None or pdb is np.nan:
            pdockq_list.append(None)
            continue

        model = pdb_numpy.Coor(pdb)
        pdockq_list += compute_pdockQ(model)

    data.df["pdockq"] = pdockq_list


def mpdockq(data, verbose=True):
    r"""Compute the mpDockq [2]_ from the pdb file.

    .. math::
        pDockQ = \frac{L}{1 + e^{-k (x-x_{0})}} + b

    where:

    .. math::
        x = \overline{plDDT_{interface}} \cdot log(number \: of \: interface \: contacts)

    :math:`L = 0.728`, :math:`x0 = 309.375`, :math:`k = 0.098` and :math:`b = 0.262`.

    Implementation was inspired from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py

    Parameters
    ----------
    data : AFData
        object containing the data
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.


    References
    ----------

    .. [2] Bryant P, Pozzati G, Zhu W, Shenoy A, Kundrotas P & Elofsson A.
        Predicting the structure of large protein complexes using AlphaFold and Monte
        Carlo tree search. *Nature Communications*. vol. 13, 6028 (2022)
        https://www.nature.com/articles/s41467-022-33729-4
    """

    from pdb_numpy.analysis import compute_pdockQ

    pdockq_list = []
    disable = False if verbose else True

    for pdb in tqdm(data.df["pdb"], total=len(data.df["pdb"]), disable=disable):
        if pdb is None or pdb is np.nan:
            pdockq_list.append(None)
            continue

        model = pdb_numpy.Coor(pdb)
        pdockq_list += compute_pdockQ(
            model, cutoff=8.0, L=0.728, x0=309.375, k=0.098, b=0.262
        )

    data.df.loc[:, "mpdockq"] = pdockq_list


def pdockq2(data, verbose=True):
    r"""
    Compute pdockq2 from the pdb file [3]_.

    .. math::
        pDockQ_2 = \frac{L}{1 + exp [-k*(X_i-X_0)]} + b

    with

    .. math::
        X_i = \langle \frac{1}{1+(\frac{PAE_{int}}{d_0})^2} \rangle * \langle pLDDT \rangle_{int}

    References:

    .. [3] : https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714

    """

    from pdb_numpy.analysis import compute_pdockQ2

    pdockq_list = []

    max_chain_num = 0
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num

    for i in range(max_chain_num):
        pdockq_list.append([])

    disable = False if verbose else True

    if "data_file" not in data.df.columns:
        raise ValueError(
            "No \`data_file\` column found in the dataframe. pae scores are required to compute pdockq2."
        )

    for pdb, data_path in tqdm(
        zip(data.df["pdb"], data.df["data_file"]),
        total=len(data.df["pdb"]),
        disable=disable,
    ):
        if (
            pdb is not None
            and pdb is not np.nan
            and data_path is not None
            and data_path is not np.nan
        ):
            model = pdb_numpy.Coor(pdb)
            # with open(json_path) as f:
            #     local_json = json.load(f)
            # pae_array = np.array(local_json["pae"])
            pae_array = get_pae(data_path)

            pdockq2 = compute_pdockQ2(model, pae_array)

            for i in range(max_chain_num):
                if i < len(pdockq2):
                    pdockq_list[i].append(pdockq2[i][0])
                else:
                    pdockq_list[i].append(None)

        else:
            for i in range(max_chain_num):
                pdockq_list[i].append(None)

    # print(pdockq_list)
    for i in range(max_chain_num):
        data.df.loc[:, f"pdockq2_{chr(65+i)}"] = pdockq_list[i]


def inter_chain_pae(data, fun=np.mean, verbose=True):
    """Read the PAE matrix and extract the average inter chain PAE.

    Parameters
    ----------
    data : AFData
        object containing the data
    fun : function
        function to apply to the PAE scores
    verbose : bool
        print progress bar

    Returns
    -------
    None
    """
    pae_list = []

    disable = False if verbose else True

    if "data_file" not in data.df.columns:
        raise ValueError(
            "No 'data_file' column found in the dataframe. pae scores are required to compute pdockq2."
        )

    for query, data_path in tqdm(
        zip(data.df["query"], data.df["data_file"]),
        total=len(data.df["data_file"]),
        disable=disable,
    ):
        if data_path is not None and data_path is not np.nan:
            pae_array = get_pae(data_path)

            chain_lens = data.chain_length[query]
            chain_len_sums = np.cumsum([0] + chain_lens)
            chain_ids = data.chains[query]

            pae_dict = {}

            for i in range(len(chain_lens)):
                for j in range(len(chain_lens)):
                    pae_val = fun(
                        pae_array[
                            chain_len_sums[i] : chain_len_sums[i + 1],
                            chain_len_sums[j] : chain_len_sums[j + 1],
                        ]
                    )
                    pae_dict[f"PAE_{chain_ids[i]}_{chain_ids[j]}"] = pae_val

            pae_list.append(pae_dict)
        else:
            pae_list.append({})

    pae_df = pd.DataFrame(pae_list)

    for col in pae_df.columns:
        data.df.loc[:, col] = pae_df.loc[:, col].to_numpy()


def compute_LIS_matrix(
    pae_array,
    chain_length,
    pae_cutoff=12.0,
):
    r"""Compute the LIS score as define in [1]_.

    Implementation was inspired from implementation in https://github.com/flyark/AFM-LIS

    Parameters
    ----------
    pae_array : np.array
        array of predicted PAE
    chain_length : list
        list of chain lengths
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A

    Returns
    -------
    list
        LIS scores

    References
    ----------

    .. [1] Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. "Enhanced
        Protein-Protein Interaction Discovery via AlphaFold-Multimer" bioRxiv (2024).
        https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1

    """

    if pae_array is None:
        return None

    chain_len_sums = np.cumsum([0] + chain_length)

    # Use list instead of array, because
    # df[column].iloc[:] = LIS_list does not work with numpy array
    LIS_list = []

    trans_matrix = np.zeros_like(pae_array)
    mask = pae_array < pae_cutoff
    trans_matrix[mask] = 1 - pae_array[mask] / pae_cutoff

    for i in range(len(chain_length)):
        i_start = chain_len_sums[i]
        i_end = chain_len_sums[i + 1]
        local_LIS_list = []
        for j in range(len(chain_length)):
            j_start = chain_len_sums[j]
            j_end = chain_len_sums[j + 1]

            submatrix = trans_matrix[i_start:i_end, j_start:j_end]

            if np.any(submatrix > 0):
                local_LIS_list.append(submatrix[submatrix > 0].mean())
            else:
                local_LIS_list.append(0)
        LIS_list.append(local_LIS_list)

    return LIS_list


def LIS_matrix(data, pae_cutoff=12.0, verbose=True):
    """
    Compute the LIS score as define in [2]_.

    Implementation was inspired from implementation in:

    .. [2] https://github.com/flyark/AFM-LIS

    Parameters
    ----------
    data : AFData
        object containing the data
    pae_cutoff : float
        cutoff for PAE matrix values, default is 12.0 A
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The dataframe is modified in place.
    """
    LIS_matrix_list = []

    max_chain_num = 0
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num

    disable = False if verbose else True

    for query, data_path in tqdm(
        zip(data.df["query"], data.df["data_file"]),
        total=len(data.df["query"]),
        disable=disable,
    ):
        if data.chain_length[query] is None:
            LIS_matrix_list.append(None)
            continue

        pae_array = get_pae(data_path)
        LIS_matrix = compute_LIS_matrix(pae_array, data.chain_length[query], pae_cutoff)
        LIS_matrix_list.append(LIS_matrix)

    assert len(LIS_matrix_list) == len(data.df["query"])
    data.df.loc[:, "LIS"] = LIS_matrix_list
