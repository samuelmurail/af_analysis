#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pandas as pd
import json
import pdb_numpy
import seaborn as sns
import matplotlib.pyplot as plt
from cmcrameri import cm
from tqdm.auto import tqdm
import json
import logging
import ipywidgets as widgets

from .format import colabfold_1_5, af3_webserver, afpulldown, boltz1, chai1, default
from . import sequence, plot
from .analysis import get_pae, extract_fields_file


# Autorship information
__author__ = "Alaa Reguei, Samuel Murail"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License version 2"
__version__ = "0.1.3"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"

# Logging
logger = logging.getLogger(__name__)


plddt_main_atom_list = [
    "CA",
    "P",
    "ZN",
    "MG",
    "CL",
    "CA",
    "NA",
    "MN",
    "K",
    "FE",
    "CU",
    "CO",
]


class Data:
    """Data class

    Parameters
    ----------
    dir : str
        Path to the directory containing the `log.txt` file.
    format : str
        Format of the data.
    df : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    chains : dict
        Dictionary containing the chains of each query.
    chain_length : dict
        Dictionary containing the length of each chain of each query.

    Methods
    -------
    read_directory(directory, keep_recycles=False)
        Read a directory.
    export_csv(path)
        Export the dataframe to a csv file.
    import_csv(path)
        Import a csv file to the dataframe.
    add_json()
        Add json files to the dataframe.
    extract_data()
        Extract json/npz files to the dataframe.
    add_pdb()
        Add pdb files to the dataframe.
    add_fasta(csv)
        Add fasta sequence to the dataframe.
    keep_last_recycle()
        Keep only the last recycle for each query.
    plot_maxscore_as_col(score, col, hue='query')
        Plot the maxscore as a function of a column.
    plot_pae(index, cmap=cm.vik)
        Plot the PAE matrix.
    plot_plddt(index_list)
        Plot the pLDDT.
    show_3d(index)
        Show the 3D structure.
    plot_msa(filter_qid=0.15, filter_cov=0.4)
        Plot the msa from the a3m file.
    show_plot_info()
        Show the plot info.

    """

    def __init__(
        self, directory=None, data_dict=None, csv=None, verbose=True, format=None
    ):
        """ """

        if directory is not None:
            self.read_directory(directory, verbose=verbose, format=format)
        elif csv is not None:
            self.format = "csv"
            self.import_csv(csv)
        elif data_dict is not None:
            assert "pdb" in data_dict.keys()
            assert "query" in data_dict.keys()
            assert "data_file" in data_dict.keys()

            self.df = pd.DataFrame(data_dict)
            self.dir = None
            self.df["format"] = "custom"
            self.set_chain_length()

    def read_directory(self, directory, keep_recycles=False, verbose=True, format=None):
        """Read a directory.

        If the directory contains a `log.txt` file, the format is set to `colabfold_1.5`.

        Parameters
        ----------
        directory : str
            Path to the directory containing the `log.txt` file.
        keep_recycles : bool
            Keep only the last recycle for each query.
        verbose : bool
            Print information about the directory.

        Returns
        -------
        None
        """
        self.dir = directory

        if format == "colabfold_1.5" or os.path.isfile(
            os.path.join(directory, "log.txt")
        ):
            self.format = "colabfold_1.5"
            self.df = colabfold_1_5.read_log(directory, keep_recycles)
            self.df["format"] = "colabfold_1.5"
            self.add_pdb(verbose=verbose)
            self.add_json(verbose=verbose)
        elif format == "AF3_webserver" or os.path.isfile(
            os.path.join(directory, "terms_of_use.md")
        ):
            self.format = "AF3_webserver"
            self.df = af3_webserver.read_dir(directory)
            self.df["format"] = "AF3_webserver"
        elif format == "AlphaPulldown" or os.path.isfile(
            os.path.join(directory, "ranking_debug.json")
        ):
            self.format = "AlphaPulldown"
            self.df = afpulldown.read_dir(directory)
            self.df["format"] = "AlphaPulldown"
        elif format == "boltz1" or (
            os.path.isdir(os.path.join(directory, "predictions"))
        ):
            self.format = "boltz1"
            self.df = boltz1.read_dir(directory)
            self.df["format"] = "boltz1"
        elif format == "chai1" or os.path.isfile(
            os.path.join(directory, "msa_depth.pdf")
        ):
            self.format = "chai1"
            self.df = chai1.read_dir(directory)
            self.df["format"] = "chai1"
        else:
            self.format = "default"
            self.df = default.read_dir(directory)
            self.df["format"] = "default"
            self.add_json()

        self.set_chain_length()

    def set_chain_length(self):
        """Find chain information from the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.chains = {}
        self.chain_length = {}
        for querie in self.df["query"].unique():
            # print(querie, self.df[self.df['query'] == querie])
            first_model = pdb_numpy.Coor(
                self.df[self.df["query"] == querie].iloc[0]["pdb"]
            )
            self.chains[querie] = list(np.unique(first_model.models[0].chain))
            self.chain_length[querie] = [
                len(
                    np.unique(
                        first_model.models[0].uniq_resid[
                            first_model.models[0].chain == chain
                        ]
                    )
                )
                for chain in self.chains[querie]
            ]

    def export_csv(self, path):
        """Export the dataframe to a csv file.

        Parameters
        ----------
        path : str
            Path to the csv file.

        Returns
        -------
        None
        """

        self.df.to_csv(path, index=False)

    def import_csv(self, path):
        """Import a csv file to the dataframe.

        Parameters
        ----------
        path : str
            Path to the csv file.

        Returns
        -------
        None
        """

        self.df = pd.read_csv(path)
        self.dir = os.path.dirname(self.df["pdb"][0])

        self.chains = {}
        self.chain_length = {}
        for querie in self.df["query"].unique():
            # print(querie, self.df[self.df['query'] == querie])
            first_model = pdb_numpy.Coor(
                self.df[self.df["query"] == querie].iloc[0]["pdb"]
            )
            self.chains[querie] = list(np.unique(first_model.models[0].chain))
            self.chain_length[querie] = [
                len(
                    np.unique(
                        first_model.models[0].uniq_resid[
                            first_model.models[0].chain == chain
                        ]
                    )
                )
                for chain in self.chains[querie]
            ]

    def add_json(self, verbose=True):
        """Add json files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_json(self.df, self.dir, verbose=verbose)
        if self.format == "default":
            default.add_json(self.df, self.dir, verbose=verbose)

    def extract_data(self):
        """Extract json/npz files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        index_list = []
        data_list = []
        for index, data_path in zip(self.df.index, self.df["data_file"]):
            if data_path is not None:
                if data_path.endswith(".json"):
                    with open(data_path, "r") as f:
                        data = json.load(f)
                    data_list.append(data)
                    index_list.append(index)
                elif data_path.endswith(".npz"):
                    data_npz = np.load(data_path)
                    data = {}
                    for key in data_npz.keys():
                        data[key] = data_npz[key]
                    data_list.append(data)
                    index_list.append(index)

        new_column = {}
        for key in data_list[0].keys():
            new_column[key] = []
        for data in data_list:
            for key in data.keys():
                new_column[key].append(data[key])

        for key in new_column.keys():
            self.df.loc[:, key] = None
            self.df.loc[index_list, key] = pd.Series(new_column[key], index=index_list)

    def extract_fields(self, fields, disable=False):
        """Extract fields from data files to the dataframe.

        Parameters
        ----------
        fields : list
            List of fields to extract.
        disable : bool
            Disable the progress bar.

        Returns
        -------
        None
        """

        values_list = []
        for field in fields:
            values_list.append([])
        for data_path in tqdm(
            self.df["data_file"], total=len(self.df["data_file"]), disable=disable
        ):
            if data_path is not None:
                local_values = extract_fields_file(data_path, fields)

                for i in range(len(fields)):
                    values_list[i].append(local_values[i])

            else:
                for i in range(len(fields)):
                    values_list[i].append(None)

        for i, field in enumerate(fields):
            self.df[field] = None
            new_col = pd.Series(values_list[i])
            self.df.loc[:, field] = new_col

    def add_pdb(self, verbose=True):
        """Add pdb files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_pdb(self.df, self.dir, verbose=verbose)

    def add_fasta(self, csv):
        """Add fasta sequence to the dataframe.

        Parameters
        ----------
        csv : str
            Path to the csv file containing the fasta sequence.

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_fasta(self.df, csv)

    def keep_last_recycle(self):
        """Keep only the last recycle for each query."""

        idx = (
            self.df.groupby(["query", "seed", "model", "weight"])["recycle"].transform(
                "max"
            )
            == self.df["recycle"]
        )
        self.df = self.df[idx]

    def plot_maxscore_as_col(self, score, col, hue="query"):
        col_list = self.df[col].unique()
        query_list = self.df[hue].unique()
        # print(col_list)
        # print(query_list)

        out_list = []

        for query in query_list:
            # print(query)
            query_pd = self.df[self.df[hue] == query]

            for column in col_list:
                # print(column)
                # ~print()

                col_pd = query_pd[query_pd[col] <= column]
                # print(col_pd[score])
                # print(column, len(col_pd))
                # print(col, col_pd.columns)

                if len(col_pd) > 0:
                    out_list.append(
                        {hue: query, score: col_pd[score].max(), col: column}
                    )
                    # print(column, len(col_pd), col_pd[score].max())

        max_pd = pd.DataFrame(out_list)

        fig, ax = plt.subplots()
        sns.lineplot(max_pd, x=col, y=score, hue=hue)

        return (fig, ax)

    def plot_pae(self, index, cmap=cm.vik):
        row = self.df.iloc[index]

        if row["json"] is None:
            return None
        pae_array = get_pae(row["json"])

        fig, ax = plt.subplots()
        res_max = sum(self.chain_length[row["query"]])

        img = ax.imshow(
            pae_array,
            cmap=cmap,
            vmin=0.0,
            vmax=30.0,
        )

        plt.hlines(
            np.cumsum(self.chain_length[row["query"]][:-1]) - 0.5,
            xmin=-0.5,
            xmax=res_max,
            colors="black",
        )

        plt.vlines(
            np.cumsum(self.chain_length[row["query"]][:-1]) - 0.5,
            ymin=-0.5,
            ymax=res_max,
            colors="black",
        )

        plt.xlim(-0.5, res_max - 0.5)
        plt.ylim(res_max - 0.5, -0.5)
        chain_pos = []
        len_sum = 0
        for longueur in self.chain_length[row["query"]]:
            chain_pos.append(len_sum + longueur / 2)
            len_sum += longueur

        ax.set_yticks(chain_pos)
        ax.set_yticklabels(self.chains[row["query"]])
        cbar = plt.colorbar(img)
        cbar.set_label("Predicted Aligned Error (Å)", rotation=270)
        cbar.ax.get_yaxis().labelpad = 15

        return (fig, ax)

    def get_plddt(self, index):
        """Extract the pLDDT array either from the pdb file or form the
        json/plddt files.

        Parameters
        ----------
        index : int
            Index of the dataframe.

        Returns
        -------
        np.array
            pLDDT array.
        """

        row = self.df.iloc[index]

        if row["format"] in ["AF3_webserver", "csv", "AlphaPulldown"]:
            model = pdb_numpy.Coor(row["pdb"])
            plddt_array = model.models[0].beta[
                np.isin(model.models[0].name, plddt_main_atom_list)
            ]
            return plddt_array

        if row["format"] in ["boltz1"]:
            data_npz = np.load(row["plddt"])
            plddt_array = data_npz["plddt"]
            return plddt_array * 100

        if row["data_file"] is None:
            return None
        elif row["data_file"].endswith(".json"):
            with open(row["data_file"]) as f:
                local_json = json.load(f)

            if "plddt" in local_json:
                plddt_array = np.array(local_json["plddt"])
            else:
                return None
        elif row["data_file"].endswith(".npz"):
            data_npz = np.load(row["data_file"])
            if "plddt" in data_npz:
                plddt_array = data_npz["plddt"]
            else:
                return None

        return plddt_array

    def plot_plddt(self, index_list=None):
        if index_list is None:
            index_list = range(len(self.df))

        fig, ax = plt.subplots()

        for index in index_list:
            plddt_array = self.get_plddt(index)
            plt.plot(np.arange(1, len(plddt_array) + 1), plddt_array)

        plt.vlines(
            np.cumsum(self.chain_length[self.df.iloc[index_list[0]]["query"]][:-1]),
            ymin=0,
            ymax=100.0,
            colors="black",
        )
        plt.ylim(0, 100)
        plt.xlim(0, sum(self.chain_length[self.df.iloc[index_list[0]]["query"]]))
        plt.xlabel("Residue")
        plt.ylabel("predicted LDDT")

        return (fig, ax)

    def show_3d(self, index):
        row = self.df.iloc[index]

        if row["pdb"] is None:
            return (None, None)

        import nglview as nv

        # Bug with show_file
        # view = nv.show_file(row['pdb'])
        view = nv.show_structure_file(row["pdb"])
        # view.add_component(ref_coor[0])
        # view.clear_representations(1)
        # view[1].add_cartoon(selection="protein", color='blue')
        # view[1].add_licorice(selection=":A", color='blue')
        # view[0].add_licorice(selection=":A")
        return view

    def plot_msa(self, filter_qid=0.15, filter_cov=0.4):
        """
        Plot the msa from the a3m file.

        Parameters
        ----------
        filter_qid : float
            Minimal sequence identity to keep a sequence.
        filter_cov : float
            Minimal coverage to keep a sequence.

        Returns
        -------
        None

        ..Warning only tested with colabfold 1.5
        """

        raw_list = os.listdir(self.dir)
        file_list = []
        for file in raw_list:
            if file.endswith(".a3m"):
                file_list.append(file)

        for a3m_file in file_list:
            logger.info(f"Reading MSA file:{a3m_file}")
            querie = a3m_file.split("/")[-1].split(".")[0]

            a3m_lines = open(os.path.join(self.dir, a3m_file), "r").readlines()[1:]
            seqs, mtx, nams = sequence.parse_a3m(
                a3m_lines=a3m_lines, filter_qid=filter_qid, filter_cov=filter_cov
            )
            logger.info(f"- Keeping {len(seqs):6} sequences for plotting.")
            feature_dict = {}
            feature_dict["msa"] = sequence.convert_aa_msa(seqs)
            feature_dict["num_alignments"] = len(seqs)

            if len(seqs) == sum(self.chain_length[querie]):
                feature_dict["asym_id"] = []
                for i, chain_len in enumerate(self.chain_length[querie]):
                    feature_dict["asym_id"] += [i + 1.0] * chain_len
                feature_dict["asym_id"] = np.array(feature_dict["asym_id"])

            fig = plot.plot_msa_v2(feature_dict)
            plt.show()

    def count_msa_seq(self):
        """
        Count for each chain the number of sequences in the MSA.

        Parameters
        ----------
        None

        Returns
        -------
        None

        ..Warning only tested with colabfold 1.5
        """

        raw_list = os.listdir(self.dir)
        file_list = []
        for file in raw_list:
            if file.endswith(".a3m"):
                file_list.append(file)

        alignement_len = {}

        for a3m_file in file_list:
            logger.info(f"Reading MSA file:{a3m_file}")
            querie = a3m_file.split("/")[-1].split(".")[0]

            a3m_lines = open(os.path.join(self.dir, a3m_file), "r").readlines()[1:]
            seqs, mtx, nams = sequence.parse_a3m(
                a3m_lines=a3m_lines, filter_qid=0, filter_cov=0
            )
            feature_dict = {}
            feature_dict["msa"] = sequence.convert_aa_msa(seqs)
            feature_dict["num_alignments"] = len(seqs)

            seq_dict = {}
            for chain in self.chains[querie]:
                seq_dict[chain] = 0

            chain_len_list = self.chain_length[querie]
            chain_list = self.chains[querie]
            seq_len = sum(chain_len_list)

            # Treat the cases of homomers
            # I compare the length of each sequence with the other ones
            # It is wrong and should be FIXED
            # The original sequence should be retrieved from eg. the pdb file
            if len(seqs[0]) != seq_len:
                new_chain_len = []
                new_chain_list = []
                for i, seq_len in enumerate(chain_len_list):
                    if seq_len not in chain_len_list[:i]:
                        new_chain_len.append(seq_len)
                        new_chain_list.append(chain_list[i])

                chain_len_list = new_chain_len
                chain_list = new_chain_list
                seq_len = sum(chain_len_list)

            assert (
                len(seqs[0]) == seq_len
            ), f"len(seqs[0])={len(seqs[0])} != seq_len={seq_len}"

            for seq in seqs:
                start = 0
                for i, num in enumerate(chain_len_list):
                    gap_num = seq[start : start + num].count("-")
                    if gap_num < num:
                        seq_dict[chain_list[i]] += 1
                    start += num

            alignement_len[
                querie
            ] = seq_dict  # [seq_dict[chain] for chain in self.chains[querie]]
        return alignement_len

    def show_plot_info(self, cmap=cm.vik):
        """
        Need to solve the issue with:

        ```
        %matplotlib ipympl
        ```

        plots don´t update when changing the model number.

        """

        model_widget = widgets.IntSlider(
            value=1,
            min=1,
            max=len(self.df),
            step=1,
            description="model:",
            disabled=False,
        )
        display(model_widget)

        def show_model(rank_num):
            fig, (ax_plddt, ax_pae) = plt.subplots(1, 2, figsize=(10, 4))
            plddt_array = self.get_plddt(rank_num - 1)
            (plddt_plot,) = ax_plddt.plot(plddt_array)
            query = self.df.iloc[model_widget.value - 1]["query"]
            json_file = self.df.iloc[model_widget.value - 1]["json"]
            ax_plddt.vlines(
                np.cumsum(self.chain_length[query][:-1]),
                ymin=0,
                ymax=100.0,
                colors="black",
            )
            ax_plddt.set_ylim(0, 100)
            res_max = sum(self.chain_length[query])
            ax_plddt.set_xlim(0, res_max)
            ax_plddt.set_xlabel("Residue")
            ax_plddt.set_ylabel("predicted LDDT")

            pae_array = get_pae(json_file)
            ax_pae.imshow(
                pae_array,
                cmap=cmap,
                vmin=0.0,
                vmax=30.0,
            )
            ax_pae.vlines(
                np.cumsum(self.chain_length[query][:-1]),
                ymin=-0.5,
                ymax=res_max,
                colors="yellow",
            )
            ax_pae.hlines(
                np.cumsum(self.chain_length[query][:-1]),
                xmin=-0.5,
                xmax=res_max,
                colors="yellow",
            )
            ax_pae.set_xlim(-0.5, res_max - 0.5)
            ax_pae.set_ylim(res_max - 0.5, -0.5)
            chain_pos = []
            len_sum = 0
            for longueur in self.chain_length[query]:
                chain_pos.append(len_sum + longueur / 2)
                len_sum += longueur
            ax_pae.set_yticks(chain_pos)
            ax_pae.set_yticklabels(self.chains[query])
            plt.show(fig)

        output = widgets.Output(layout={"width": "95%"})
        display(output)

        with output:
            show_model(model_widget.value)
            # logger.info(results['metric'][0][rank_num - 1]['print_line'])

        def on_value_change(change):
            output.clear_output()
            with output:
                show_model(model_widget.value)

        model_widget.observe(on_value_change, names="value")


def concat_data(data_list):
    """Concatenate data from a list of Data objects.

    Parameters
    ----------
    data_list : list
        List of Data objects.

    Returns
    -------
    Data
        Concatenated Data object.
    """

    concat = Data(directory=None, csv=None)

    concat.df = pd.concat([data.df for data in data_list], ignore_index=True)
    concat.chains = data_list[0].chains
    concat.chain_length = data_list[0].chain_length
    concat.format = data_list[0].format
    for i in range(1, len(data_list)):
        concat.chains.update(data_list[i].chains)
        concat.chain_length.update(data_list[i].chain_length)

    return concat


def read_multiple_alphapulldown(directory):
    """Read multiple directories containing AlphaPulldown data.

    Parameters
    ----------
    directory : str
        Path to the directory containing the directories.

    Returns
    -------
    Data
        Concatenated Data object.
    """

    dir_list = [
        name
        for name in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, name))
    ]
    data_list = []

    for dir in dir_list:
        if "ranking_debug.json" in os.listdir(os.path.join(directory, dir)):
            data_list.append(Data(os.path.join(directory, dir)))

    if len(data_list) == 0:
        raise ValueError("No AlphaPulldown data found in the directory.")
    return concat_data(data_list)
