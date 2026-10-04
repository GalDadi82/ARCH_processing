"""Mutation QC plots, collected into a single PDF.

Port of the notebook's plotting cells: mutations per sample, recurrence vs. average VAF
(labelled by gene), and the number of samples carrying a mutation in each gene.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    valid_samples = pd.read_csv(snakemake.input.samples)
    df = pd.read_table(snakemake.input.listed, sep="\t")

    emit(f"listed table: {df.shape}")
    emit(f"Found {df['Sample_Name'].nunique()}/{valid_samples.shape[0]} samples with mutations")

    with PdfPages(snakemake.output.pdf) as pdf:
        # 1. Mutations per sample.
        sample_mut_counts = df["Sample_Name"].value_counts().sort_index().reset_index()
        sample_mut_counts.columns = ["Sample_Name", "Count"]

        plt.figure(figsize=(10, 6))
        sns.histplot(
            data=sample_mut_counts,
            x="Count",
            color="steelblue",
            bins=max(int(sample_mut_counts["Count"].max()), 1),
        )
        plt.xlabel("Number of Mutations")
        plt.ylabel("Number of Samples")
        plt.title("Number of Mutations per Sample")
        plt.xticks(rotation=90)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # 2. Recurrence vs. average VAF, labelled by gene.
        df_grouped = (
            df.groupby(["CHR", "POS", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene"])
            .agg(
                number_of_samples=("Sample_Name", "nunique"),
                avg_vaf=("avg_VAF", "mean"),
                max_vaf=("avg_VAF", "max"),
            )
            .reset_index()
        )

        plt.figure(figsize=(12, 7))
        ax = sns.scatterplot(
            data=df_grouped,
            x="number_of_samples",
            y="avg_vaf",
            s=80,
            color="skyblue",
            edgecolor="gray",
        )

        # Only label the interesting points, else the plot is unreadable.
        to_show = df_grouped[
            (df_grouped["number_of_samples"] >= 3) | (df_grouped["avg_vaf"] > 0.5)
        ]
        for _, row in to_show.iterrows():
            ax.text(
                row["number_of_samples"],
                row["avg_vaf"],
                str(row["Gene.refGene"]),
                fontsize=10,
                ha="center",
                va="top",
            )

        plt.xlabel("Number of Samples")
        plt.ylabel("Average VAF")
        plt.title("Number of Samples vs. Average VAF (labeled by gene)")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # 3. Samples carrying a mutation, per gene.
        df_u = df.drop_duplicates(subset=["Sample_Name", "Gene.refGene"])
        gene_counts = df_u["Gene.refGene"].value_counts().sort_index().reset_index()
        gene_counts.columns = ["Gene", "Count"]
        gene_counts = gene_counts.sort_values("Count", ascending=False)

        plt.figure(figsize=(10, 6))
        sns.barplot(data=gene_counts, x="Gene", y="Count", color="darkred")
        plt.xlabel("Gene")
        plt.ylabel("Number of samples")
        plt.title("Number of samples with mutations per Gene")
        plt.xticks(rotation=90)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # The notebook printed this table to eyeball recurrent / high-VAF calls; keep it in
    # the log rather than losing it.
    emit("Recurrent or high-VAF variants:")
    emit(to_show.sort_values("number_of_samples", ascending=False).to_string(index=False))
    emit(f"Saved plots to: {snakemake.output.pdf}")
