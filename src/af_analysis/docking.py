import numpy as np
import pdb_numpy

from tqdm.auto import tqdm
from . import data, analysis


# Autorship information
__author__ = "Alaa Reguei, Samuel Murail"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License version 2"
__version__ = "0.1.4"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"

"""
The module contains functions to extract and compute docking scores.

.. warning:
    The ligand chain is assumed to be the last chain in the list of chains.
"""


def pae_pep(my_data, fun=np.mean, verbose=True):
    """Extract the PAE score for the receptor(s)-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the PAE scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    pep_rec_pae_list = []
    rec_pep_pae_list = []

    disable = False if verbose else True

    for i, (query, data_file) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["data_file"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        pae = data.get_pae(data_file)

        # print(f"0:{cum_sum_chain[-2]} , {cum_sum_chain[-2]}:{cum_sum_chain[-1]}")
        # print(pae.shape)

        if pae is None:
            pep_rec_pae_list.append(None)
            rec_pep_pae_list.append(None)
            continue

        rec_pep_pae = fun(
            pae[0 : cum_sum_chain[-2], cum_sum_chain[-2] : cum_sum_chain[-1]]
        )
        pep_rec_pae = fun(
            pae[cum_sum_chain[-2] : cum_sum_chain[-1], 0 : cum_sum_chain[-2]]
        )

        pep_rec_pae_list.append(pep_rec_pae)
        rec_pep_pae_list.append(rec_pep_pae)

    my_data.df.loc[:, "PAE_pep_rec"] = pep_rec_pae_list
    my_data.df.loc[:, "PAE_rec_pep"] = rec_pep_pae_list


def pae_contact_pep(my_data, fun=np.mean, cutoff=8.0, verbose=True, max_pae=30.98):
    """Extract the PAE score for the receptor(s)-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the PAE scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    pep_rec_pae_list = []
    rec_pep_pae_list = []

    disable = False if verbose else True

    for i, (query, data_file, pdb) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["data_file"], my_data.df["pdb"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chains = my_data.chains[query]

        if pdb is None or data_file is None:
            pep_rec_pae_list.append(None)
            rec_pep_pae_list.append(None)
            continue

        pae = data.get_pae(data_file)

        model = pdb_numpy.Coor(pdb)
        model_CA = model.select_atoms("name CA")
        contact_lig = model_CA.select_atoms(
            f"chain {chains[-1]} and within {cutoff} of chain {' '.join(chains[:-1])}"
        )
        contact_rec = model_CA.select_atoms(
            f"chain {' '.join(chains[:-1])} and within {cutoff} of chain {chains[-1]}"
        )

        if contact_lig.len == 0 or contact_rec.len == 0:
            pep_rec_pae_list.append(max_pae)
            rec_pep_pae_list.append(max_pae)
            continue

        rec_mask = np.zeros(pae.shape)
        lig_mask = np.zeros(pae.shape)

        rec_mask[contact_rec.residue, :] = 1
        lig_mask[:, contact_lig.residue] = 1
        pair_mask = np.logical_and(rec_mask, lig_mask)

        rec_pep_pae = fun(pae[pair_mask])
        pep_rec_pae = fun(pae[pair_mask.T])

        pep_rec_pae_list.append(pep_rec_pae)
        rec_pep_pae_list.append(rec_pep_pae)

    my_data.df.loc[:, "PAE_contact_pep_rec"] = pep_rec_pae_list
    my_data.df.loc[:, "PAE_contact_rec_pep"] = rec_pep_pae_list


def plddt_pep(my_data, fun=np.mean, verbose=True):
    """Extract the pLDDT score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the pLDDT scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """
    pep_plddt_list = []

    disable = False if verbose else True

    for i, (query, pdb) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["pdb"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        plddt = my_data.get_plddt(i)
        if plddt is None:
            pep_plddt_list.append(None)
            continue
        pep_plddt_list.append(fun(plddt[cum_sum_chain[-2] : cum_sum_chain[-1]]))

    my_data.df.loc[:, "plddt_pep"] = pep_plddt_list


def plddt_contact_pep(my_data, fun=np.mean, cutoff=8.0, verbose=True):
    """Extract the pLDDT score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the pLDDT scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """
    lig_plddt_list = []
    rec_plddt_list = []

    disable = False if verbose else True

    for i, (query, pdb) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["pdb"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        chains = my_data.chains[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        if pdb is None:
            lig_plddt_list.append(None)
            rec_plddt_list.append(None)
            continue

        model = pdb_numpy.Coor(pdb)
        model_CA = model.select_atoms("name CA")
        contact_lig = model_CA.select_atoms(
            f"chain {chains[-1]} and within {cutoff} of chain {' '.join(chains[:-1])}"
        )
        contact_rec = model_CA.select_atoms(
            f"chain {' '.join(chains[:-1])} and within {cutoff} of chain {chains[-1]}"
        )

        if contact_lig.len > 0:
            lig_plddt_list.append(fun(contact_lig.beta))
        else:
            lig_plddt_list.append(0)

        if contact_rec.len > 0:
            rec_plddt_list.append(fun(contact_rec.beta))
        else:
            rec_plddt_list.append(0)

    my_data.df.loc[:, "plddt_contact_lig"] = lig_plddt_list
    my_data.df.loc[:, "plddt_contact_rec"] = rec_plddt_list


def LIS_pep(my_data, pae_cutoff=12.0, fun=np.max, verbose=True):
    """Compute the LIS score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A
    fun : function
        function to apply to the LIS matrix

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.

    """

    analysis.LIS_matrix(my_data, pae_cutoff=pae_cutoff, verbose=verbose)

    pep_LIS_list = []
    pep_LIS2_list = []

    for query, LIS in zip(my_data.df["query"], my_data.df["LIS"]):
        if LIS is None:
            pep_LIS_list.append(None)
            pep_LIS2_list.append(None)
            continue

        chain_num = len(my_data.chains[query])
        LIS_array = np.array(LIS)
        pep_LIS_list.append(fun(LIS_array[0 : chain_num - 1, chain_num - 1]))
        pep_LIS2_list.append(fun(LIS_array[chain_num - 1, 0 : chain_num - 1]))

    my_data.df.loc[:, "LIS_rec_pep"] = pep_LIS2_list
    my_data.df.loc[:, "LIS_pep_rec"] = pep_LIS_list


def pdockq2_lig(my_data, verbose=True):
    """Compute the LIS score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A
    fun : function
        function to apply to the LIS matrix

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.

    """

    analysis.pdockq2(my_data, verbose=verbose)

    old_query = ""
    pdockq2_list = []

    for index, row in my_data.df.iterrows():
        if row["query"] != old_query:
            old_query = row["query"]
            chains = my_data.chains[old_query]
            lig_chain = chains[-1]
            rec_chains = chains[:-1]

        pdockq2_list.append(row[f"pdockq2_{lig_chain}"])

    my_data.df.loc[:, "pdockq2_lig"] = pdockq2_list
