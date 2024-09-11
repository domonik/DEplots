import pandas as pd
import gffutils
from tempfile import TemporaryDirectory
import os

GFF_COLNAMES = ['seqid', 'source', 'featuretype', 'start', 'end', 'score', 'strand', 'frame', 'attributes']


def parse_attributes(attributes):
    """
    Parse the attributes field of a GFF3 file and return a dictionary of key-value pairs.

    :param attributes: The attributes string from a GFF3 file.
    :return: A dictionary of key-value pairs.
    """
    if pd.isna(attributes):
        return {}
    return dict(item.split('=') for item in attributes.split(';') if '=' in item)


def read_gff3(filepath):
    """
    Read a GFF3 file into a pandas DataFrame with appropriate column names.

    :param filepath: Path to the GFF3 file.
    :return: A pandas DataFrame containing the GFF3 data.
    """
    column_names = GFF_COLNAMES
    df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=column_names)
    attributes_df = df['attributes'].apply(parse_attributes).apply(pd.Series)

    # Combine the original DataFrame with the expanded attributes DataFrame
    df = pd.concat([df.drop(columns=['attributes']), attributes_df], axis=1)
    return df


def filter_gff_by_interval(df, seqid, start, end):
    """
    Filter the GFF3 DataFrame to include only rows where the coordinates overlap with the specified interval.

    :param df: The GFF3 DataFrame.
    :param start: The start of the interval.
    :param end: The end of the interval.
    :return: A filtered pandas DataFrame.
    """
    # Filter based on the overlap condition

    condition = (df['start'] <= end) & (df['end'] >= start)
    df = df[df["seqid"] == seqid]
    filtered_df = df[condition]

    return filtered_df


def find_gene_coordinates(gff, gene, anno_type):
    df = gff[gff[anno_type].str.contains(gene, case=False) == True]
    return df


def coords_from_gff_row(row):
    return row["seqid"], row["start"], row.end, row.strand


def read_gff_via_gffutils(file, dbname=None):
    if dbname is None:
        dir = TemporaryDirectory()
        os.path.join(dir.name, "GFF.db")
    if os.path.exists(dbname):
        db = gffutils.FeatureDB(dbfn=dbname)
    else:
        db = gffutils.create_db(file, dbfn=dbname, force=True, keep_order=True, merge_strategy="merge", id_spec="ID")
    return db


