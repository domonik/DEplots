import gffutils
import yaml
import os
import pandas as pd
from plotly import express as px
import pickle
from DEplots.readCount import precompute_from_design
from DEplots.gff_helper import read_gff3, read_gff_via_gffutils

DIRPATH = os.path.dirname(os.path.abspath(__file__))


def get_data(config_file: str = None, run_dir: str = None,):
    if config_file is None:
        config_file = os.path.join(DIRPATH, "default_config.yaml")
    with open(config_file, "r") as handle:
        config = yaml.safe_load(handle)

    if run_dir:
        config["run_dir"] = run_dir
    return read_files(config)


def get_coverage_data(config_file):
    if config_file is None:
        config_file = os.path.join(DIRPATH, "default_config.yaml")
    with open(config_file, "r") as handle:
        config = yaml.safe_load(handle)

    if not config["coverage"]["design"]:
        return None, None, None, None
    else:
        design = pd.read_csv(config["coverage"]["design"], sep="\t")
        design["index"] = list(range(len(design)))
        file = config["coverage"]["precomputed_file"]

        if not os.path.exists(file):
            coverage = precompute_from_design(design)

            with open(file, "wb") as handle:
                pickle.dump(coverage, handle)
        else:
            with open(file, "rb") as handle:
                coverage = pickle.load(handle)
        if config["coverage"]["use_gffutils"]:
            dbname = config["coverage"]["dbname"] if config["coverage"]["dbname"] else None
            gff = read_gff_via_gffutils(config["coverage"]["gff"], dbname=dbname)
            mapping = _compute_gff_lines(gff)
            gff = dbname
        else:
            gff = read_gff3(config["coverage"]["gff"])

        return design, coverage, gff, mapping


def _compute_gff_lines(gff: gffutils.FeatureDB):
    genes = gff.execute(
        """
        SELECT f.*
        FROM features f
        LEFT JOIN relations r ON f.id = r.child
        WHERE f.featuretype = 'gene' OR r.parent IS NULL
        ORDER BY f.seqid, f.start, f.strand;
        """
    )
    lines = []
    mapping = {}
    old_seqid = None
    for row in genes:
        seq_id = row["seqid"]
        if seq_id != old_seqid:
            lines = []
            old_seqid = seq_id
        idx = row["id"]
        start, end = row['start'], row['end']
        placed = False

        # Try to place the span in an existing line
        for i, line_end in enumerate(lines):
            if start > line_end:
                mapping[idx] = i
                lines[i] = end
                placed = True
                break
        # If no suitable line was found, create a new line
        if not placed:
            mapping[idx] = len(lines)
            lines.append(end)
    return mapping



def read_files(config):
    if config["run_dir"]:
        runs = os.listdir(config["run_dir"])
        config_files = {run: os.path.join(config["run_dir"], run, "config.yml") for run in runs}
    else:
        assert config["config_files"], "No files specified for display"
        config_files = config["config_files"]
    name_mapping = config["name_mapping"]
    if config["add_data"]:
        add_data = config["add_data"]
        add_data = pd.read_csv(add_data, sep="\t", index_col=0)
        add_data.columns = add_data.columns.str.replace('_', ' ')

    else:
        add_data = None

    dash_data = {}
    multiindex_data = {}
    multiindex_count = {}
    for name, file in config_files.items():
        with open(file, "r") as handle:
            d = yaml.safe_load(handle)
            dirname = os.path.dirname(file)
            dash_data[name] = {
                "comparisons": {},
                "config": d
            }
            for (condition, baseline) in zip(d["conditions"], d["baselines"]):
                comp_file = os.path.join(dirname, f"PipelineData/DESeqResults/DESeqResult_c{condition}_vs_b{baseline}.tsv")

                assert os.path.isfile(comp_file), f"File {comp_file} not found"
                df = pd.read_csv(comp_file, sep="\t")

                if name_mapping:
                    cname = name_mapping[condition] if condition in name_mapping else condition
                    bname = name_mapping[baseline] if baseline in name_mapping else baseline
                else:
                    cname = condition
                    bname = baseline
                comp_str = f"{cname} vs {bname}"
                dash_data[name]["comparisons"][comp_str] = {}
                if comp_str not in multiindex_data:
                    multiindex_data[comp_str] = [[], []]
                    multiindex_count[comp_str] = 0
                multiindex_data[comp_str][0] = multiindex_data[comp_str][0] + [(name, col) for col in df.columns]
                multiindex_data[comp_str][1].append(df)
                old = multiindex_count[comp_str]
                multiindex_count[comp_str] = old + df.shape[1]
                dash_data[name]["comparisons"][comp_str]["deseq"] = (old, multiindex_count[comp_str])
                dash_data[name]["comparisons"][comp_str]["enrich"] = {}
                dash_data[name]["comparisons"][comp_str]["condition"] = cname
                dash_data[name]["comparisons"][comp_str]["baseline"] = bname

                for enrich in ["GO", "KEGG"]:
                    enrich_file_up = os.path.join(dirname,
                                                  f"PipelineData/Enrichment/{enrich}Enrichment_up_c{condition}_vs_b{baseline}.tsv")
                    enrich_file_down = os.path.join(dirname,
                                                    f"PipelineData/Enrichment/{enrich}Enrichment_down_c{condition}_vs_b{baseline}.tsv")
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich] = {}
                    if os.path.isfile(enrich_file_up):
                        enrich_up = pd.read_csv(enrich_file_up, sep="\t")
                    else:
                        enrich_up = None
                        print(f"File {enrich_file_up} not found")
                    if os.path.isfile(enrich_file_down):
                        enrich_down = pd.read_csv(enrich_file_down, sep="\t")
                    else:
                        enrich_down = None
                        print(f"File {enrich_file_down} not found")
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich][cname] = enrich_up
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich][bname] = enrich_down

                gsea_file = os.path.join(
                    dirname,
                    f"PipelineData/Enrichment/GSEAGO_c{condition}_vs_b{baseline}.tsv"
                )
                gsea_plot_data = os.path.join(
                    dirname,
                    f"PipelineData/Enrichment/GSEAGO_plot_data_c{condition}_vs_b{baseline}.tsv"
                )
                if os.path.isfile(gsea_file) and os.path.isfile(gsea_plot_data):
                    gsea_plot_data = pd.read_csv(gsea_plot_data, sep="\t")
                    gsea_file = pd.read_csv(gsea_file, sep="\t")
                    gsea_file = gsea_file.drop("core_enrichment", axis=1)
                    dash_data[name]["comparisons"][comp_str]["gsea"] = {}
                    dash_data[name]["comparisons"][comp_str]["gsea"]["plot_data"] = gsea_plot_data
                    dash_data[name]["comparisons"][comp_str]["gsea"]["df"] = gsea_file

    for key, value in multiindex_data.items():
        names, dfs = value
        if add_data is not None:
            dfs.append(add_data)
            names += [("Additional Data", col) for col in add_data.columns]
        multi_index = pd.MultiIndex.from_tuples(names)

        df = pd.concat(dfs, axis=1, join="outer")
        df.columns = multi_index

        df[("Name", "Name")] = df.index
        cols = [("Name", "Name")] + [col for col in df.columns if col != ("Name", "Name")]
        df = df[cols]
        multiindex_data[key] = df
    return dash_data, multiindex_data


DEFAULT_PLOTLY_COLORS = {
    "not-enriched": "#344A9A",
    "enriched": "#00a082",
    "Selected": "#8f6b30",
    "placeholder": "#ffe863",
    "placeholder2": "#f5c2ed"

}
DEFAULT_PLOTLY_COLORS_LIST = list(DEFAULT_PLOTLY_COLORS.values()) + px.colors.DEFAULT_PLOTLY_COLORS
LAYOUT = {
    "template": "plotly_white",
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': 'rgb(219, 219, 219)',
    "font": {"color": "black", "size": 16},
    "xaxis": {"showline": True, "mirror": True, "linecolor": "black"},
    "yaxis": {"showline": True, "mirror": True, "linecolor": "black"},
    "margin": {"b": 10, "t": 10},

}
DARK_LAYOUT = {
    "template": "plotly_white",
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': 'rgba(0,0,0,0)',
    "font": {"color": "white", "size": 16},
    "xaxis": {"showline": True, "mirror": True, "linecolor": "white"},
    "yaxis": {"showline": True, "mirror": True, "linecolor": "white"},
    "margin": {"b": 10, "t": 10}
}
UP_COLOR_LIGHT = "#00a082"
UP_COLOR_DARK = "#00a082"
DOWN_COLOR_LIGHT = "#344A9A"
DOWN_COLOR_DARK = "#f5c2ed"
DASH_DATA = None
COVERAGE_DATA = None
COVERAGE_DESIGN = None
GFF = None
LINE_MAPPING = None
