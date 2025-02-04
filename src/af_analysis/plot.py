import matplotlib.pyplot as plt
import numpy as np
import ipywidgets as widgets
from cmcrameri import cm

from .analysis import get_pae


def plot_msa_v2(feature_dict, sort_lines=True, dpi=100):
    """
    Taken from:
    https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
    """

    seq = feature_dict["msa"][0]
    # print("len(seq), seq", len(seq), seq)

    if "asym_id" in feature_dict:
        Ls = [0]
        k = feature_dict["asym_id"][0]
        for i in feature_dict["asym_id"]:
            if i == k:
                Ls[-1] += 1
            else:
                Ls.append(1)
            k = i
    else:
        Ls = [len(seq)]
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"]

    # print("asym_id:", feature_dict["asym_id"])
    # print(len(feature_dict["asym_id"]))
    # print(f"Ln: {Ln}  Ls:{Ls}")

    msa = feature_dict["msa"][:N]
    gap = msa != 21
    qid = msa == seq

    # print("gap.shape:", gap.shape)
    # for i in range(len(Ls)):
    #    print(f"i: {i}  Ls[i]: {Ls[i]}  Ln[i+1]: {Ln[i+1]}")
    #    print(f"gap[:, Ln[i]: Ln[i+1]]: {gap[:, Ln[i]: Ln[i+1]]}")
    #    print(f"gap[:, Ln[i]: Ln[i+1]].max(-1): {gap[:, Ln[i]: Ln[i+1]].max(-1)}")

    gapid = np.stack([gap[:, Ln[i] : Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack(
            [qid_[:, Ln[i] : Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1
        ).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)
    fig = plt.figure(figsize=(8, 5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(
        lines,
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
        extent=(0, lines.shape[1], 0, lines.shape[0]),
    )
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color="black")
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    return fig


def show_info(
    data_af2, cmap=cm.vik, score_list=["pLDDT", "pTM", "ipTM", "ranking_confidence"]
):
    """
    Use with
    ```
    %matplotlib widget
    ```
    """

    model_widget = widgets.IntSlider(
        value=1,
        min=1,
        max=len(data_af2.df),
        step=1,
        description="model:",
        disabled=False,
    )
    display(model_widget)

    rank_num = 1

    fig, (ax_plddt, ax_pae) = plt.subplots(1, 2, figsize=(10, 4))
    plddt_array = data_af2.get_plddt(rank_num - 1)
    (plddt_plot,) = ax_plddt.plot(plddt_array)
    query = data_af2.df.iloc[model_widget.value - 1]["query"]
    data_file = data_af2.df.iloc[model_widget.value - 1]["data_file"]
    vline_plot = ax_plddt.vlines(
        np.cumsum(data_af2.chain_length[query][:-1]),
        ymin=0,
        ymax=100.0,
        colors="black",
    )
    ax_plddt.set_ylim(0, 100)
    res_max = sum(data_af2.chain_length[query])
    ax_plddt.set_xlim(0, res_max)
    ax_plddt.set_xlabel("Residue")
    ax_plddt.set_ylabel("predicted LDDT")

    pae_array = get_pae(data_file)
    pae_plot = ax_pae.imshow(
        pae_array,
        cmap=cmap,
        vmin=0.0,
        vmax=30.0,
    )
    vline_pae = ax_pae.vlines(
        np.cumsum(data_af2.chain_length[query][:-1]),
        ymin=-0.5,
        ymax=res_max,
        colors="yellow",
    )
    hline_pae = ax_pae.hlines(
        np.cumsum(data_af2.chain_length[query][:-1]),
        xmin=-0.5,
        xmax=res_max,
        colors="yellow",
    )
    ax_pae.set_xlim(-0.5, res_max - 0.5)
    ax_pae.set_ylim(res_max - 0.5, -0.5)
    chain_pos = []
    len_sum = 0
    for longueur in data_af2.chain_length[query]:
        chain_pos.append(len_sum + longueur / 2)
        len_sum += longueur
    ax_pae.set_yticks(chain_pos)
    ax_pae.set_yticklabels(data_af2.chains[query])
    plt.show(fig)

    # out_score = widgets.Output(layout={'border': '1px solid black'})
    out_score = widgets.HTML()
    display(out_score)

    pattern = "<p style='display: inline-block; width:100px'> <strong>{score_name:15} : </strong> {score_value:7.2f} </p>"

    for score in score_list:
        if score in data_af2.df.columns:
            out_score.value += pattern.format(
                score_name=score,
                score_value=data_af2.df.iloc[model_widget.value - 1][score],
            )
            # (f"<div> <strong>{score:15} : </strong> {data_af2.df.iloc[model_widget.value - 1][score]:7.2f} </div>")

    def update_model(change):
        rank_num = model_widget.value
        # print("Update")
        plddt_array = data_af2.get_plddt(rank_num - 1)
        res_num = len(plddt_array)
        plddt_plot.set_data(range(res_num), plddt_array)
        ax_plddt.set_xlim(0, len(plddt_array))

        query = data_af2.df.iloc[rank_num - 1]["query"]
        vline_plot.set_segments(
            [
                np.array([[x, 0], [x, 100]])
                for x in np.cumsum(data_af2.chain_length[query][:-1])
            ]
        )
        # ax_plddt.set_title(self.chain_length[query][:-1])

        data_file = data_af2.df.iloc[rank_num - 1]["data_file"]
        pae_array = get_pae(data_file)
        pae_plot.set_extent((0, res_num, 0, res_num))
        pae_plot.set_data(pae_array)
        ax_pae.set_xlim(0, res_num)
        ax_pae.set_ylim(0, res_num)

        vline_pae.set_segments(
            [
                np.array([[x, -0.5], [x, res_num]])
                for x in np.cumsum(data_af2.chain_length[query][:-1])
            ]
        )
        hline_pae.set_segments(
            [
                np.array([[-0.5, res_num - x], [res_num, res_num - x]])
                for x in np.cumsum(data_af2.chain_length[query][:-1])
            ]
        )
        chain_pos = []
        len_sum = 0
        for longueur in data_af2.chain_length[query]:
            chain_pos.append(res_num - (len_sum + longueur / 2))
            len_sum += longueur
        ax_pae.set_yticks(chain_pos)
        ax_pae.set_yticklabels(data_af2.chains[query])
        fig.canvas.draw()

        new_out_score = ""
        for score in score_list:
            if score in data_af2.df.columns:
                # new_out_score += (f"<div> <strong>{score:15} : </strong> {data_af2.df.iloc[model_widget.value - 1][score]:7.2f} </div>")
                new_out_score += pattern.format(
                    score_name=score,
                    score_value=data_af2.df.iloc[model_widget.value - 1][score],
                )

        out_score.value = new_out_score

        # out_score.clear_output()
        # with out_score:
        #    for score in score_list:
        #        if score in data_af2.df.columns:
        #            print(f"{score:15} : {data_af2.df.iloc[model_widget.value - 1][score]:7.2f}")

    model_widget.observe(update_model, names="value")
