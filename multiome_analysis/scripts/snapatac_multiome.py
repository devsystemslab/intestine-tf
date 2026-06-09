import sys
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

import json
import snapatac2 as snap
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import marsilea as ma
import marsilea.plotter as mp
from pathlib import Path
from snapatac2.genome import Genome


def main():
    # Load the genome and annotation
    genomes_dir = Path("/projects/site/pred/ihb-intestine-evo/lukas_area/genomes")
    genome = Genome(fasta=genomes_dir / "hg38.fa", annotation=genomes_dir / "hg38.sorted.gtf")

    # Remove any chromosomes with a dot in their name
    for key in list(genome.chrom_sizes.keys()):
        if "." in key:
            del genome.chrom_sizes[key]

    # Set the base project directory
    base_path = Path("/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/processed")

    # List of fragment files to process
    files = [
        "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/multiome_data/atac/ITBOGE001_Fujii.tsv.gz",
        "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/parse_data/atac/ITBOGE013_YBP2_14_Multiome_1.tsv",
        "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/parse_data/atac/ITBOGE013_YBP2_14_Multiome_2.tsv",
    ]
    files = {Path(f).stem.replace(".tsv.gz", ""): f for f in files}

    # Import the fragments and create AnnData objects
    adatas = snap.pp.import_fragments(
        [files[fl] for fl in files.keys()],
        file=[Path(base_path) / (name + ".atac.raw.h5ad") for name in files.keys()],
        chrom_sizes=genome,
        min_num_fragments=1000,
        sorted_by_barcode=False,
        n_jobs=18,
    )


if __name__ == "__main__":
    main()
