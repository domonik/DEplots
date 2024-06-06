import yaml
import os
import pandas as pd


def read_files(config_file):
    with open(config_file, "r") as handle:
        config = yaml.safe_load(handle)

    config_files = config["config_files"]
    name_mapping = config["name_mapping"]
    add_data = config["add_data"]
    add_data = pd.read_csv(add_data, sep="\t", index_col=0)

    dash_data = {}
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
                df["Name"] = df.index
                df = df[["Name"] + [col for col in df.columns if col != "Name"]]

                df = df.merge(add_data, left_index=True, right_index=True)
                if name_mapping:
                    cname = name_mapping[condition]
                    bname = name_mapping[baseline]
                else:
                    cname = condition
                    bname = baseline
                comp_str = f"{cname} vs {bname}"
                dash_data[name]["comparisons"][comp_str] = {}
                dash_data[name]["comparisons"][comp_str]["deseq"] = df
                dash_data[name]["comparisons"][comp_str]["enrich"] = {}
                dash_data[name]["comparisons"][comp_str]["condition"] = cname
                dash_data[name]["comparisons"][comp_str]["baseline"] = bname

                for enrich in ["GO", "KEGG"]:
                    enrich_file_up = os.path.join(dirname,
                                                  f"PipelineData/Enrichment/{enrich}Enrichment_up_c{condition}_vs_b{baseline}.tsv")
                    enrich_file_down = os.path.join(dirname,
                                                    f"PipelineData/Enrichment/{enrich}Enrichment_down_c{condition}_vs_b{baseline}.tsv")

                    assert os.path.isfile(enrich_file_up), f"File {enrich_file_up} not found"
                    assert os.path.isfile(enrich_file_down), f"File {enrich_file_down} not found"

                    enrich_up = pd.read_csv(enrich_file_up, sep="\t")
                    enrich_down = pd.read_csv(enrich_file_down, sep="\t")
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich] = {}
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich][cname] = enrich_up
                    dash_data[name]["comparisons"][comp_str]["enrich"][enrich][bname] = enrich_down
    return dash_data


