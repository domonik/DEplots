import os

import numpy as np
import pandas as pd
import pysam
import pickle


def get_coverage_data(design: pd.DataFrame, output_file, callback=None):
    hash = pd.util.hash_pandas_object(design)

    if not os.path.exists(output_file):
        coverage = precompute_from_design(design, callback)

        with open(output_file, "wb") as handle:
            pickle.dump((coverage, hash), handle)

    with open(output_file, "rb") as handle:
        coverage, d_hash = pickle.load(handle)
    e = ValueError(f"Hash values do not match please delete {os.path.abspath(output_file)} and rerun")
    try:
        if not all(hash == d_hash):
            raise e
    except ValueError:
        raise e


def precompute_from_design(design, callback=None):
    coverage = {

    }
    for idx, row in design.iterrows():
        file = row["File"]
        index = row["index"]
        pysam.index(file)
        samfile = pysam.AlignmentFile(file, "rb")
        for contig in samfile.references:
            yp = samfile.count_coverage(contig=contig, read_callback=filter_strand("+"))
            ym = samfile.count_coverage(contig=contig, read_callback=filter_strand("-"))
            yp = np.asarray(yp).sum(axis=0)
            ym = np.asarray(ym).sum(axis=0)
            if "Size factors" in design.columns:
                yp *= row["Size factors"]
                ym *= row["Size factors"]
            if contig not in coverage:
                coverage[contig] = {"+": np.empty((m_idx, yp.shape[0])), "-": np.empty((m_idx, ym.shape[0]))}
            coverage[contig]["+"][index] = yp
            coverage[contig]["-"][index] = ym
        if callback:
            progress = str(int(((idx+1) / m_idx) * 100))
            print(f"progress: {progress}")
            callback(progress)
    if callback:
        callback("100")
    return coverage


def read_contigs(design):
    contigs = {}
    for idx, row in design.iterrows():
        file = row["File"]
        #pysam.index(file)
        samfile = pysam.AlignmentFile(file, "rb")
        for (ref, size) in zip(samfile.references, samfile.lengths):
            if ref in contigs:
                assert contigs[ref] == size
            else:
                contigs[ref] = size

    return contigs


def filter_strand(strand):
    if strand == "+":
        def filter_read(read):
         return not read.is_reverse
    elif strand == "-":
        def filter_read(read):
            return read.is_reverse
    else:
        return "all"
    return filter_read
