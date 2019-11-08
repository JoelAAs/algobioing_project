import seaborn as sns
import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt
import matplotlib

matplotlib.use('Agg')


class AssociationGroupAdjustment:
    def __init__(self, assoc_file, grouping_file, output_path):
        self.grouping = pd.read_csv(grouping_file, sep="\t")
        self.assoc_df = pd.read_csv(assoc_file, sep="\t")
        self.assoc_df = self.assoc_df.sort_values("P", ascending=False)
        self.assoc_df = self.assoc_df.drop_duplicates(subset=["external_gene_name", "refsnp_id"])
        self.current_selection = None
        self.plot_path = output_path

    def set_current_selection(self, group):
        current_genes = self.grouping[self.grouping[group]].index.values
        self.current_selection = self.assoc_df[
            self.assoc_df["external_gene_name"].isin(current_genes)].copy()

    def get_pi_zero(self):
        """
        Cubic spline from p-values to estimate null-dist
        """

        def _calc_pi_zero(lam, p_values):
            over = sum(p_values >= lam)
            under = len(p_values) * (1 - lam)
            return over / under

        # lambda_range = np.arange(0.001, 0.95, 0.05)
        # pi_0_range = [_calc_pi_zero(x, self.current_selection["P"]) for x in lambda_range]
        # f = sp.interpolate.CubicSpline(x=lambda_range, y=pi_0_range)
        # return f([1])[0]

        p_values = self.current_selection["P"]
        selected = p_values[p_values > 0.75]
        selected = selected[selected != 1]
        return np.mean([_calc_pi_zero(l, p_values) for l in selected])

    def set_fdr(self):
        """
        Calculate local min FDR q-values
        """

        def _fdr(t, pi0, p_values):
            return pi0 * len(p_values) * t / sum(p_values <= t)

        def _redef(q_raw):
            prev = q_raw[-1] + 1
            for i, value in enumerate(q_raw):
                if value > prev:
                    q_raw[i] = prev
                else:
                    prev = value
            return q_raw

        pi0 = self.get_pi_zero()
        self.current_selection["q_values"] = self.current_selection["P"].apply(
            lambda x: _fdr(x, pi0, self.current_selection["P"]))
        self.current_selection["q_values"] = _redef(self.current_selection["q_values"].to_list())

    def get_dataframe(self):
        return self.current_selection.copy()

    def get_fdr_per_group(self, plot=False):
        """ q-value adjustment depending on group"""
        all_groups = self.grouping.columns.values
        all_data = [None] * len(all_groups)
        for i, group in enumerate(all_groups):
            print("calculating q-values for " + group)
            self.set_current_selection(group)
            self.set_fdr()
            df_out = self.get_dataframe()
            all_data[i] = df_out
            if plot:
                self.plot_part(df=df_out, group_name=group)

        return all_data

    def plot_part(self, df, group_name):

        # Plot config
        output_name = self.plot_path + str(group_name).strip().replace(" ", "_") + ".png"
        fig, axs = plt.subplots(
            nrows=4,
            figsize=(15, 15))

        # Manhattan plot
        manhattan_df = df.copy()
        manhattan_df = manhattan_df.sort_values(["CHR", "BP"])
        manhattan_df["ind"] = range(len(manhattan_df))
        manhattan_df["-log(p)"] = -np.log10(manhattan_df["P"])
        chr_group = manhattan_df.groupby(("CHR"))
        colours = ["blue", "red"]
        x_lab = []
        x_lab_pos = []
        for num, (name, group) in enumerate(chr_group):
            group.plot(
                kind="scatter",
                x="ind", y="-log(p)",
                color=colours[num % len(colours)],
                ax=axs[0])
            x_lab.append(name)
            x_lab_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

        axs[0].set_xticks(x_lab_pos)
        axs[0].set_xticklabels(x_lab)
        axs[0].set_xlim([0, len(manhattan_df)])
        axs[0].set_ylim([0, 10])
        axs[0].set_xlabel("Chromosome", fontsize=15)
        axs[0].set_ylabel("-log(p)", fontsize=15)

        # P-values distribution
        axs[1].set_xlabel("Freq", fontsize=15)
        axs[1].set_ylabel("p-value", fontsize=15)
        sns.distplot(df["P"], ax=axs[1])

        # P-values vs Q-values
        axs[2].set_xlabel("q-value", fontsize=15)
        axs[2].set_ylabel("p-value", fontsize=15)
        sns.scatterplot(x="P", y="q_values", data=df, ax=axs[2])

        # Q-values manhattan
        manhattan_df = df.copy()
        manhattan_df = manhattan_df.sort_values(["CHR", "BP"])
        manhattan_df["ind"] = range(len(manhattan_df))
        manhattan_df["-log(q)"] = -np.log10(manhattan_df["q_values"])
        chr_group = manhattan_df.groupby(("CHR"))
        colours = ["blue", "red"]
        x_lab = []
        x_lab_pos = []
        for num, (name, group) in enumerate(chr_group):
            group.plot(
                kind="scatter",
                x="ind", y="-log(q)",
                color=colours[num % len(colours)],
                ax=axs[3])
            x_lab.append(name)
            x_lab_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

        axs[3].set_xticks(x_lab_pos)
        axs[3].set_xticklabels(x_lab)
        axs[3].set_xlim([0, len(manhattan_df)])
        axs[3].set_xlabel("Chromosome", fontsize=15)
        axs[3].set_ylabel("-log(q)", fontsize=15)

        fig.savefig(re.sub("[\[\]()]", "", output_name))
        plt.close()
