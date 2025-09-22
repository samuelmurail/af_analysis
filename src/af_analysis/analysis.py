import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import json
import os
import logging
import itertools
import pickle

import pdb_numpy
from pdb_numpy.geom import distance_matrix

# Logging
logger = logging.getLogger(__name__)

# Autorship information
__author__ = "Alaa Reguei"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License version 2"
__version__ = "0.1.4"
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

    if data_file is None or data_file == "" or pd.isna(data_file):
        return None

    if data_file.endswith(".json"):
        return extract_pae_json(data_file)
    elif data_file.endswith(".npz"):
        return extract_pae_npz(data_file)
    elif data_file.endswith(".npy"):
        return extract_pae_npy(data_file)
    elif data_file.endswith(".pkl"):
        return extract_pae_pkl(data_file)
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
    """Get the PAE matrix from a npz file.

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


def extract_pae_npy(npy_file):
    """Get the PAE matrix from a npy file.

    Parameters
    ----------
    npy_file : str
        Path to the npy file.

    Returns
    -------
    np.array
        PAE matrix.
    """

    pae_array = np.load(npy_file)

    return pae_array


def extract_pae_pkl(pkl_file):
    """Get the PAE matrix from a pkl file.

    Parameters
    ----------
    pkl_file : str
        Path to the pkl file.

    Returns
    -------
    np.array
        PAE matrix.
    """


    try:
        data_pkl = np.load(pkl_file, allow_pickle=True)
    except pickle.UnpicklingError as e:
        logger.error(f"Error loading PAE from {pkl_file}: {e}")
        return None
    
    pae_array = data_pkl["predicted_aligned_error"]
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
    max_chain_val = []
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num
            max_chain_val = data.chains[query]

    for i in range(max_chain_num):
        pdockq_list.append([])

    disable = False if verbose else True

    if "data_file" not in data.df.columns:
        raise ValueError(
            "No \"data_file\" column found in the dataframe. pae scores are required to compute pdockq2."
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
        data.df.loc[:, f"pdockq2_{max_chain_val[i]}"] = pdockq_list[i]




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


def PAE_matrix(data, verbose=True, fun=np.average):
    """
    Compute the average (or something else) PAE matrix.

    Parameters
    ----------
    data : AFData
        object containing the data
    verbose : bool
        print progress bar
    fun : function
        function to apply to the PAE scores

    Returns
    -------
    None
        The dataframe is modified in place.
    """

    PAE_avg_list = []

    disable = False if verbose else True

    for query, data_path in tqdm(
        zip(data.df["query"], data.df["data_file"]),
        total=len(data.df["query"]),
        disable=disable,
    ):
        if data.chain_length[query] is None:
            PAE_avg_list.append(None)
            continue

        pae_array = get_pae(data_path)
        chain_len_cum = np.cumsum([data.chain_length[query]])
        chain_len_cum = np.insert(chain_len_cum, 0, 0)

        avg_matrix = np.zeros((len(chain_len_cum) - 1, len(chain_len_cum) - 1))

        for i in range(len(chain_len_cum) - 1):
            for j in range(len(chain_len_cum) - 1):
                avg_matrix[i, j] = fun(
                    pae_array[
                        chain_len_cum[i] : chain_len_cum[i + 1],
                        chain_len_cum[j] : chain_len_cum[j + 1],
                    ]
                )

        PAE_avg_list.append(avg_matrix)

    assert len(PAE_avg_list) == len(data.df["query"])
    data.df.loc[:, "PAE_fun"] = PAE_avg_list


def read_ftdmp_raw_score(raw_path):
    """Read raw ftdmp score files

    Parameters
    ----------
    raw_path : str
        Path to the raw score file

    Returns
    -------
    raw_score : pandas.DataFrame
        Dataframe containing the raw score data
    """

    with open(raw_path, "r") as filin:
        # extract the header
        header = filin.readline().strip().split()

        # extract the data
        data = filin.readlines()
        # extract the data
        raw_data = [line.strip().split() for line in data]
        data = []
        for line in raw_data:
            data_line = []
            for i, val in enumerate(line):
                if i == 0:
                    data_line.append(val)
                else:
                    data_line.append(float(val))
            data.append(data_line)

        # convert to pandas dataframe
        score_df = pd.DataFrame(data, columns=header)

    return score_df


def extract_ftdmp(
    ftdmp_result_path, score_list=["raw_scoring_results_without_ranks.txt"]
):
    """Read ftdmp output files

    Parameters
    ----------
    ftdmp_result_path : str
        Path to the ftdmp output directory

    Returns
    -------
    my_data : AFData
        object containing the data
    """

    jobs_dir = os.listdir(os.path.join(ftdmp_result_path, "jobs"))
    if len(jobs_dir) > 1:
        logger.warning(f"Two outputs directory founded in {jobs_dir}")

    df_list = []

    for job_dir in jobs_dir:
        local_job_dir = os.path.join(ftdmp_result_path, "jobs", job_dir)

        file_list = os.listdir(local_job_dir)

        for file in file_list:
            if file in score_list:
                logger.info(f"Reading ftdmp score file : {file}")
                local_score = os.path.join(ftdmp_result_path, "jobs", job_dir, file)
                df_list.append(read_ftdmp_raw_score(local_score))

    # my_data.df["ID"] = [
    #     os.path.basename(file_path) if file_path is not None else None
    #     for file_path in my_data.df["pdb"]
    # ]

    # return df_list

    # for df in df_list:
    #     if len(df) == 0:
    #         continue

    #     my_data.df = my_data.df.merge(df, on="ID", how="inner")
    return df_list

def compute_ftdmp(
    my_data,
    ftdmp_path=None,
    out_path="tmp_ftdmp",
    score_list=["raw_scoring_results_without_ranks.txt"],
    env=None,
    keep_tmp=False,
):
    """Compute ftdmp scores

    Parameters
    ----------
    ftdmp_path : str
        Path to the ftdmp output directory

    Returns
    -------
    my_data : AFData
        object containing the data
    """

    import shutil
    import subprocess

    if ftdmp_path is None:
        ftdmp_exe_path = shutil.which("ftdmp-qa-all")
        if ftdmp_exe_path is None:
            logger.warning("Software ftdmp-qa-all not found in PATH")
            return
    else:
        ftdmp_exe_path = os.path.expanduser(os.path.join(ftdmp_path, "ftdmp-qa-all"))

    # Test if pytorch is installed
    try:
        import torch
    except ImportError:
        logger.warning("Pytorch not found, ftdmp will not work")
        return

    if not os.path.exists(out_path):
        os.makedirs(out_path)


    # Check that all pdb files are in the same directory
    if not all(
        os.path.dirname(pdb) == os.path.dirname(my_data.df["pdb"].iloc[0])
        for pdb in my_data.df["pdb"]
    ):

        pdb_run_list = []
        pdb_dir_list = []
        for pdb in my_data.df["pdb"].tolist():
            pdb_dirname = os.path.dirname(pdb)
            if pdb_dirname not in pdb_dir_list:
                pdb_dir_list.append(pdb_dirname)
                pdb_run_list.append([])

            # Add the pdb file to the corresponding directory
            index = pdb_dir_list.index(pdb_dirname)
            pdb_run_list[index].append(pdb)

        logger.info(f"For Ftdmp, all PDB files must be in the same directory. Ftdmp will be launched {len(pdb_run_list)} times.")
    else:
        # If all pdb files are in the same directory, we can run ftdmp on all of them at once
        pdb_run_list = [my_data.df["pdb"].tolist()]
        logger.info("All PDB files are in the same directory. Ftdmp will be launched once.")

    # If the pdb files are not in the same directory, we need to copy them to the out_path
    
    #pdb_list = [pdb for pdb in my_data.df["pdb"].tolist() if pdb is not None and not pd.isna(pdb)]

    df_list = []
    # print(pdb_run_list)

    for i, pdb_list in enumerate(pdb_run_list):
        logger.info(f"Running ftdmp on {len(pdb_list)} PDB files, step {i+1}/{len(pdb_run_list)}")

        # ls MY_AF_DIRECTORY/*.pdb | ~/Documents/Code/ftdmp/ftdmp-qa-all --workdir ftdmp_beta_amyloid_dimer

        cmd = [ftdmp_exe_path, "--workdir", out_path]

        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )

        #pdb_list = [pdb for pdb in my_data.df["pdb"].tolist() if pdb is not None and not pd.isna(pdb)]
        com_input = "\n".join(pdb_list)
        com_input += "\n"

        (stdout_data, stderr_data) = proc.communicate(com_input.encode())

        # print(stdout_data)
        # print(stderr_data)

        df_list += extract_ftdmp(ftdmp_result_path=out_path, score_list=score_list)

        if not keep_tmp:
            shutil.rmtree(out_path)

    logger.info("Ftdmp scores computed. Merging data frames.")

    # print(df_list)
    ftdmp_df = pd.DataFrame()
    
    for df in df_list:
        if len(df) == 0:
            continue

        # Add the ID column to the dataframe
        df["ID"] = df["ID"].astype(str)
        ftdmp_df = pd.concat([ftdmp_df, df], ignore_index=True)

    my_data.df["ID"] = [
        os.path.basename(file_path) if file_path is not None else None
        for file_path in my_data.df["pdb"]]  
    my_data.df = my_data.df.merge(ftdmp_df, on="ID", how="inner")


    return


def compute_dockq(data, ref_dict, verbose=True, fun=np.average, dockq_thresold = 0.3):
    """Compute the DockQ score from the PAE matrix.

    Parameters
    ----------
    data : AFData
        object containing the data
    ref_dict : dict
        dictionary containing the reference PAE matrix for each query
    verbose : bool
        print progress bar
    fun : function
        function to apply to the PAE scores
    dockq_thresold : float
        threshold with multiple chain to recompute DockQ score, default is 0.3

    Returns
    -------
    None
        The dataframe is modified in place.
    """
    from pdb_numpy import analysis

    dockq_list = []
    lrmsd_list = []
    fnat_list = []
    old_query = ""

    disable = False if verbose else True

    for query, pdb in tqdm(
        zip(data.df["query"], data.df["pdb"]),
        total=len(data.df["query"]),
        disable=disable,
    ):
        if (query not in ref_dict
            or data.chain_length[query] is None
            or pdb is None
            or pdb is np.nan
        ):
            if query not in ref_dict:
                logger.warning(f"No reference pdb structure found for query {query}.")
            if pdb is None or pdb is np.nan:
                logger.warning(f"No PDB file found for query {query}.")

            dockq_list.append(None)
            lrmsd_list.append(None)
            fnat_list.append(None)
            continue

        if query != old_query:
            ref_coor = pdb_numpy.Coor(ref_dict[query])
            native_seq = ref_coor.get_aa_seq()
            old_query = query

        model = pdb_numpy.Coor(pdb)
        model_seq = model.get_aa_seq()
        lig_chains = [
            min(model_seq.items(), key=lambda x: len(x[1].replace("-", "")))[0]
        ]
        rec_chains = [chain for chain in model_seq if chain not in lig_chains]
        native_lig_chains = [
            min(native_seq.items(), key=lambda x: len(x[1].replace("-", "")))[0]
        ]
        native_rec_chains = [
            chain for chain in native_seq if chain not in native_lig_chains
        ]

        dockq_score = analysis.dockQ(
            model,
            ref_coor,
            rec_chains=rec_chains,
            lig_chains=lig_chains,
            native_rec_chains=native_rec_chains,
            native_lig_chains=native_lig_chains,
        )

        if len(rec_chains) > 1 and dockq_score["DockQ"][0] < dockq_thresold:
            # If the model has multiple chains and the DockQ score is below the threshold,
            # recompute the DockQ score using different alignement modes
            new_results = [dockq_score]
            for chain_perm in list(itertools.permutations(rec_chains, len(rec_chains)))[1:]:
                new_results.append(analysis.dockQ(
                    model,
                    ref_coor,
                    rec_chains=list(chain_perm),
                    lig_chains=lig_chains,
                    native_rec_chains=native_rec_chains,
                    native_lig_chains=native_lig_chains,
                ))
            # Select the best result
            dockq_score = max(new_results, key=lambda x: x["DockQ"][0])


        #print(dockq_score)
        # print(f"dockq: {dockq_score['DockQ'][0]:.3f} Lrms:  {dockq_score['LRMS'][0]:.2f}")
        dockq_list.append(dockq_score["DockQ"][0])
        lrmsd_list.append(dockq_score["LRMS"][0])
        fnat_list.append(dockq_score["Fnat"][0])

    assert len(dockq_list) == len(data.df["query"])
    data.df.loc[:, "dockq"] = dockq_list
    data.df.loc[:, "lrmsd"] = lrmsd_list
    data.df.loc[:, "fnat"] = fnat_list


def ipTM_d0(data, verbose=True):
    """Compute the ipTM_d0 score from the PAE matrix.

    Implementation is based on the ipTM_d0 function from the IPSAE package
    https://github.com/DunbrackLab/IPSAE/blob/main/ipsae.py

    Cite:
    .. [1] Dunbrack RL Jr. "Rēs ipSAE loquunt: What’s wrong with AlphaFold’s
    ipTM score and how to fix it" bioRxiv (2025).

    Parameters
    ----------
    data : AFData
        object containing the data
    ref_dict : dict
        dictionary containing the reference PAE matrix for each query
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The dataframe is modified in place.
    """

    iptm_d0_list = []

    disable = False if verbose else True

    for query, data_file in tqdm(
        zip(data.df["query"], data.df["data_file"]),
        total=len(data.df["query"]),
        disable=disable,
    ):

        PAE_matrix = get_pae(data_file)
        if PAE_matrix is None:
            logger.warning(f"No PAE matrix found for query {query}.")
            iptm_d0_list.append({f"ipTM_d0_{data.chains[query][0]}_{data.chains[query][1]}":None})
            continue
    
        # Check if the PAE matrix is square
        iptm_d0_values = compute_iptm_d0_values(
            pae_array=PAE_matrix,
            chain_ids=data.chains[query],
            chain_length=data.chain_length[query],
            chain_type=data.chain_type[query],
        )

        iptm_d0_list.append(iptm_d0_values)

    assert len(iptm_d0_list) == len(data.df["query"])

    iptm_d0_df = pd.DataFrame(iptm_d0_list)

    for col in iptm_d0_df.columns:
        data.df.loc[:, col] = iptm_d0_df.loc[:, col].to_numpy()


def compute_iptm_d0_values(pae_array, chain_ids, chain_length, chain_type):
    """Compute the ipTM_d0 score from the PAE matrix.

    Parameters
    ----------
    pae_array : np.array
        array of predicted PAE
    chain_ids : list
        list of chain IDs
    chain_length : list
        list of chain lengths
    chain_type : list
        list of chain types (e.g. "protein", "nucleic_acid")

    Returns
    -------
    list
        ipTM_d0 score
    """

    # Define the ptm and d0 functions
    def ptm_func(x, d0):
        return 1.0 / (1 + (x / d0) ** 2.0)

    ptm_func_vec = np.vectorize(ptm_func)  # vector version

    # Define the d0 functions for numbers and arrays; minimum value = 1.0; from Yang and Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
    def calc_d0(L, pair_type):
        L = float(L)
        if L < 27:
            L = 27
        min_value = 1.0
        if pair_type == "nucleic_acid":
            min_value = 2.0
        d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
        return max(min_value, d0)

    def fun(matrix):
        """Function to apply to the ipTM_d0 matrix.
        This function computes the mean of each row and returns the maximum value.
        """

        matrix_res = np.zeros(matrix.shape[0])

        for i in range(matrix.shape[0]):
            matrix_res[i] = np.mean(matrix[i, :])

        max_index = np.argmax(matrix_res)
        return matrix_res[max_index]

    chain_len_sums = np.cumsum([0] + chain_length)

    iptm_d0_dict = {}
    iptm_do_sum = 0.0
    iptm_size = 0.0
    for i in range(len(chain_length)):
        for j in range(len(chain_length)):
            if i != j:
                do_chain = chain_length[i] + chain_length[j]
                type = (
                    "nucleic_acid"
                    if (chain_type[i] != "protein" or chain_type[j] != "protein")
                    else "protein"
                )

                d0 = calc_d0(do_chain, type)
                iptm_d0_matrix = ptm_func_vec(
                    pae_array[
                        chain_len_sums[i] : chain_len_sums[i + 1],
                        chain_len_sums[j] : chain_len_sums[j + 1],
                    ],
                    d0,
                )
                iptm_d0_mean = fun(iptm_d0_matrix)
                iptm_d0_dict[f"ipTM_d0_{chain_ids[i]}_{chain_ids[j]}"] = iptm_d0_mean
                iptm_do_sum += iptm_d0_mean * chain_length[i] * chain_length[j]
                iptm_size += chain_length[i] * chain_length[j]

    iptm_d0_dict[f"ipTM_d0"] = (
        iptm_do_sum / iptm_size if iptm_size > 0 else None
    )
    return iptm_d0_dict


def ipSAE(data, verbose=True, pae_cutoff=10.0, dist_cutoff=10.0):
    """Compute the ipSAE score from the PAE matrix.

        Implementation is based on the ipTM_d0 function from the IPSAE package
        https://github.com/DunbrackLab/IPSAE/blob/main/ipsae.py
    # apply cutoff
        Cite:
        .. [1] Dunbrack RL Jr. "Rēs ipSAE loquunt: What’s wrong with AlphaFold’s
        ipTM score and how to fix it" bioRxiv (2025).

        Parameters
        ----------
        data : AFData
            object containing the dipSAE(ata
        ref_dict : dict
            dictionary containing the reference PAE matrix for each query
        verbose : bool
            print progress bar

        Returns
        -------
        None
            The dataframe is modified in place.
    """

    ipSAE_list = []

    disable = False if verbose else True

    for query, data_file in tqdm(
        zip(data.df["query"], data.df["data_file"]),
        total=len(data.df["query"]),
        disable=disable,
    ):

        PAE_matrix = get_pae(data_file)
        if PAE_matrix is None:
            logger.warning(f"No PAE matrix found for query {query}.")
            ipSAE_list.append({f"ipSAE_{data.chains[query][0]}_{data.chains[query][1]}": None})
            continue
        ipSAE_matrix = compute_ipSAE_matrix(
            pae_array=PAE_matrix,
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff,
            chain_ids=data.chains[query],
            chain_length=data.chain_length[query],
            chain_type=data.chain_type[query],
        )

        ipSAE_list.append(ipSAE_matrix)

    assert len(ipSAE_list) == len(data.df["query"])

    ipSAE_df = pd.DataFrame(ipSAE_list)

    for col in ipSAE_df.columns:
        data.df.loc[:, col] = ipSAE_df.loc[:, col].to_numpy()


def compute_ipSAE_matrix(
    pae_array, pae_cutoff, dist_cutoff, chain_ids, chain_length, chain_type
):
    """Compute the ipSAE score from the PAE matrix.

    Parameters
    ----------
    pdb : str
        path to the pdb file
    pae_array : np.array
        array of predicted PAE
    pae_cutoff : float
        cutoff for PAE matrix values, default is 10.0 A
    dist_cutoff : float
        cutoff for distance between atoms, default is 10.0 A
    chain_ids : list
        list of chain IDs
    chain_length : list
        list of chain lengths
    chain_type : list
        list of chain types (e.g. "protein", "nucleic_acid")

    Returns
    -------
    list
        ipSAE score matrix
    """

    # Define the ptm and d0 functions
    def ptm_func(x, d0):
        return 1.0 / (1 + (x / d0) ** 2.0)

    ptm_func_vec = np.vectorize(ptm_func)  # vector version

    def calc_d0_array(L, pair_type):
        # Convert L to a NumPy array if it isn't already one (enables flexibility in input types)
        L = np.array(L, dtype=float)
        L = np.maximum(27, L)
        min_value = 1.0

        if pair_type == "nucleic_acid":
            min_value = 2.0

        # Calculate d0 using the vectorized operation
        return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)

    def fun(matrix):
        """Function to apply to the ipTM_d0 matrix."""
        # return np.mean(x)

        matrix_res = np.zeros(matrix.shape[0])

        for i in range(matrix.shape[0]):
            matrix_res[i] = np.mean(matrix[i, :])
            # print("matrix_res[i]", matrix_res[i])

        max_index = np.argmax(matrix_res)
        return matrix_res[max_index]

    chain_len_sums = np.cumsum([0] + chain_length)
    ipSAE_dict = {}
    for i in range(len(chain_length)):
        for j in range(len(chain_length)):
            if i != j:
                sub_pae_array = pae_array[
                    chain_len_sums[i] : chain_len_sums[i + 1],
                    chain_len_sums[j] : chain_len_sums[j + 1],
                ]

                type = (
                    "nucleic_acid"
                    if (chain_type[i] != "protein" or chain_type[j] != "protein")
                    else "protein"
                )
                valid_pairs_matrix = sub_pae_array <= pae_cutoff
                n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
                d0_res = calc_d0_array(n0res_byres_all, type)

                ipsae_d0res_byres = np.zeros(chain_length[i])

                # print(
                #     "sub_pae_array:",
                #     sub_pae_array.shape,
                #     "d0_res.shape:",
                #     d0_res.shape,
                #     "chain_length[i]:",
                #     chain_length[i],
                # )

                for k in range(chain_length[i]):
                    ptm_row_d0res = ptm_func_vec(sub_pae_array[k], d0_res[k])
                    if valid_pairs_matrix[k].any():
                        ipsae_d0res_byres[k] = ptm_row_d0res[
                            valid_pairs_matrix[k]
                        ].mean()
                    else:
                        ipsae_d0res_byres[k] = 0.0

                max_index = np.argmax(ipsae_d0res_byres)
                ipSAE_dict[f"ipSAE_{chain_ids[i]}_{chain_ids[j]}"] = ipsae_d0res_byres[
                    max_index
                ]
                # print("ipsae_d0res_byres[max_index]", ipsae_d0res_byres[max_index])

    return ipSAE_dict


def ipTM_d0_interface(data, verbose=True):
    """Compute the ipTM_d0 score from the PAE matrix.

    Implementation is based on the ipTM_d0 function from the IPSAE package
    https://github.com/DunbrackLab/IPSAE/blob/main/ipsae.py

    Cite:
    .. [1] Dunbrack RL Jr. "Rēs ipSAE loquunt: What’s wrong with AlphaFold’s
    ipTM score and how to fix it" bioRxiv (2025).

    Parameters
    ----------
    data : AFData
        object containing the data
    ref_dict : dict
        dictionary containing the reference PAE matrix for each query
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The dataframe is modified in place.
    """

    iptm_d0_list = []

    disable = False if verbose else True

    for query, pdb, data_file in tqdm(
        zip(data.df["query"], data.df["pdb"], data.df["data_file"]),
        total=len(data.df["query"]),
        disable=disable,
    ):

        PAE_matrix = get_pae(data_file)
        if PAE_matrix is None:
            logger.warning(f"No PAE matrix found for query {query}.")
            iptm_d0_list.append({f"ipTM_interface_{data.chains[query][0]}_{data.chains[query][1]}": None})
            continue
        iptm_d0_values = compute_iptm_d0_interface_values(
            pdb=pdb,
            pae_array=PAE_matrix,
            chain_ids=data.chains[query],
            chain_length=data.chain_length[query],
            chain_type=data.chain_type[query],
        )

        iptm_d0_list.append(iptm_d0_values)

    assert len(iptm_d0_list) == len(data.df["query"])

    iptm_d0_df = pd.DataFrame(iptm_d0_list)

    for col in iptm_d0_df.columns:
        data.df.loc[:, col] = iptm_d0_df.loc[:, col].to_numpy()


def compute_iptm_d0_interface_values(pdb, pae_array, chain_ids, chain_length, chain_type):
    """Compute the ipTM_d0 score from the PAE matrix.

    Parameters
    ----------
    pdb : str
        path to the pdb file
    pae_array : np.array
        array of predicted PAE
    chain_ids : list
        list of chain IDs
    chain_length : list
        list of chain lengths
    chain_type : list
        list of chain types (e.g. "protein", "nucleic_acid")

    Returns
    -------
    list
        ipTM_d0 score
    """

    # Define the ptm and d0 functions
    def ptm_func(x, d0):
        return 1.0 / (1 + (x / d0) ** 2.0)

    ptm_func_vec = np.vectorize(ptm_func)  # vector version

    # Define the d0 functions for numbers and arrays; minimum value = 1.0; from Yang and Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
    def calc_d0(L, pair_type):
        L = float(L)
        if L < 27:
            L = 27
        min_value = 1.0
        if pair_type == "nucleic_acid":
            min_value = 2.0
        d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
        return max(min_value, d0)

    def fun(matrix):
        """Function to apply to the ipTM_d0 matrix.
        This function computes the mean of each row and returns the maximum value.
        """
        if np.isnan(matrix).all():
            return 0.0
        
        matrix_res = np.zeros(matrix.shape[0])

        for i in range(matrix.shape[0]):
            matrix_res[i] = np.nanmean(matrix[i, :])

        max_index = np.nanargmax(matrix_res)
        return matrix_res[max_index]

    chain_len_sums = np.cumsum([0] + chain_length)

    model = pdb_numpy.Coor(pdb)
    model_cb = model.select_atoms("name CB C3 or (resname GLY and name CA)")
    assert model_cb.len == sum(chain_length), "Number of CB atoms does not match the sum of chain lengths."
    assert model_cb.len == pae_array.shape[0], "Number of CB atoms does not match the number of rows in the PAE matrix."
    distance = distance_matrix(model_cb.xyz, model_cb.xyz)

    iptm_d0_dict = {}
    for i in range(len(chain_length)):
        for j in range(len(chain_length)):
            if i != j:
                do_chain = chain_length[i] + chain_length[j]
                type = (
                    "nucleic_acid"
                    if (chain_type[i] != "protein" or chain_type[j] != "protein")
                    else "protein"
                )

                d0 = calc_d0(do_chain, type)
                iptm_d0_matrix = ptm_func_vec(
                    pae_array[
                        chain_len_sums[i] : chain_len_sums[i + 1],
                        chain_len_sums[j] : chain_len_sums[j + 1],
                    ],
                    d0,
                )
                sub_distance = distance[
                    chain_len_sums[i] : chain_len_sums[i + 1],
                    chain_len_sums[j] : chain_len_sums[j + 1],
                ]

                # Apply distance cutoff
                iptm_d0_matrix[sub_distance > 10.0] = np.nan
                iptm_d0_mean = fun(iptm_d0_matrix)
                iptm_d0_dict[f"ipTM_interface_{chain_ids[i]}_{chain_ids[j]}"] = iptm_d0_mean
                # print(f"ipTM_interface_{chain_ids[i]}_{chain_ids[j]}", iptm_d0_mean)

    return iptm_d0_dict
