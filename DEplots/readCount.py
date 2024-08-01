import pysam
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS
import plotly.express as px
import pandas as pd
import re

from DEplots.gff_helper import read_gff3, filter_gff_by_interval, find_gene_coordinates, coords_from_gff_row

# Regex patterns to match different CSS color formats
hex_pattern = re.compile(r'^#([0-9a-fA-F]{3,8})$')
rgb_pattern = re.compile(r'^rgb\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})\)$')
rgba_pattern = re.compile(r'^rgba\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3}),\s*(\d+(\.\d+)?)\)$')
hsl_pattern = re.compile(r'^hsl\((\d{1,3}),\s*(\d{1,3})%,\s*(\d{1,3})%\)$')
hsla_pattern = re.compile(r'^hsla\((\d{1,3}),\s*(\d{1,3})%,\s*(\d{1,3})%,\s*(\d+(\.\d+)?)\)$')


GFF3_TYPE_TO_COLOR = {
    "gene": px.colors.qualitative.Light24[0],
    "mRNA": px.colors.qualitative.Light24[1],
    "exon": px.colors.qualitative.Light24[2],
    "CDS": px.colors.qualitative.Light24[3],
    "start_codon": px.colors.qualitative.Light24[4],
    "stop_codon": px.colors.qualitative.Light24[5],
    "five_prime_UTR": px.colors.qualitative.Light24[6],
    "5UTR": px.colors.qualitative.Light24[6],
    "5'UTR": px.colors.qualitative.Light24[6],
    "three_prime_UTR": px.colors.qualitative.Light24[7],
    "3UTR": px.colors.qualitative.Light24[7],
    "3'UTR": px.colors.qualitative.Light24[7],
    "ncRNA": px.colors.qualitative.Light24[8],
    "rRNA": px.colors.qualitative.Light24[9],
    "tRNA": px.colors.qualitative.Light24[10],
    "repeat_region": px.colors.qualitative.Light24[11],
    "default": "#7CB9E8"
}


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    length = len(hex_color)
    if length == 3:
        return tuple(int(hex_color[i] * 2, 16) for i in range(3))
    elif length == 4:
        return tuple(int(hex_color[i] * 2, 16) for i in range(3)) + (int(hex_color[3] * 2, 16) / 255,)
    elif length == 6:
        return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    elif length == 8:
        return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4)) + (int(hex_color[6:8], 16) / 255,)


def hsl_to_rgb(h, s, l):
    s /= 100
    l /= 100
    c = (1 - abs(2 * l - 1)) * s
    x = c * (1 - abs((h / 60) % 2 - 1))
    m = l - c / 2
    if 0 <= h < 60:
        r, g, b = c, x, 0
    elif 60 <= h < 120:
        r, g, b = x, c, 0
    elif 120 <= h < 180:
        r, g, b = 0, c, x
    elif 180 <= h < 240:
        r, g, b = 0, x, c
    elif 240 <= h < 300:
        r, g, b = x, 0, c
    elif 300 <= h < 360:
        r, g, b = c, 0, x
    else:
        raise ValueError
    return (int((r + m) * 255), int((g + m) * 255), int((b + m) * 255))


def css_color_to_rgba(color, opacity):
    """
    Convert any acceptable CSS color to the same color with the specified opacity.

    :param color: A string representing the CSS color.
    :param opacity: A float representing the opacity value (between 0 and 1).
    :return: A string representing the color in RGBA format.
    """
    # Ensure opacity is a valid float between 0 and 1
    if not (0 <= opacity <= 1):
        raise ValueError("Opacity must be a float between 0 and 1.")


    match = hex_pattern.match(color)
    if match:
        rgba = hex_to_rgb(match.group(1))
        if len(rgba) == 4:
            return f'rgba({rgba[0]}, {rgba[1]}, {rgba[2]}, {opacity})'
        else:
            return f'rgba({rgba[0]}, {rgba[1]}, {rgba[2]}, {opacity})'

    match = rgb_pattern.match(color)
    if match:
        return f'rgba({match.group(1)}, {match.group(2)}, {match.group(3)}, {opacity})'

    match = rgba_pattern.match(color)
    if match:
        return f'rgba({match.group(1)}, {match.group(2)}, {match.group(3)}, {opacity})'

    match = hsl_pattern.match(color)
    if match:
        rgb = hsl_to_rgb(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        return f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})'

    match = hsla_pattern.match(color)
    if match:
        rgb = hsl_to_rgb(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        return f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})'

    raise ValueError("Invalid CSS color format.")

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


def get_trace(file, name, color, strand: str = None, size_factor: float = None,  **kwargs):
    samfile = pysam.AlignmentFile(file, "rb")

    y = samfile.count_coverage(read_callback=filter_strand(strand), **kwargs)
    y = np.asarray(y).sum(axis=0)
    if size_factor:
        y = y * size_factor
    s = kwargs["start"] if "start" in kwargs else 0
    trace = go.Scatter(
        x=np.arange(s, s+len(y)),
        y=y,
        name=name,
        line=dict(color=color),
        mode="lines"
    )
    return trace


def get_summary_trace(files, name, color, strand, size_factors = None, **kwargs):
    t = []
    for idx, file in enumerate(files):
        samfile = pysam.AlignmentFile(file, "rb")
        y = samfile.count_coverage(read_callback=filter_strand(strand), **kwargs)
        y = np.asarray(y).sum(axis=0)
        if size_factors:
            y = y * size_factors[idx]

        t.append(y)
    y = np.asarray(t)
    mean_y = np.mean(y, axis=0)
    median_y = np.median(y, axis=0)
    q25 = np.quantile(y, q=0.25, axis=0)
    q75 = np.quantile(y, q=0.75, axis=0)
    quantiles = np.concatenate((q75, np.flip(q25)), axis=0)
    s = kwargs["start"] if "start" in kwargs else 0
    x = np.arange(s, s + len(mean_y))
    qx = np.concatenate((x, np.flip(x)), axis=0)

    traces = [
        go.Scatter(
            x=x,
            y=median_y,
            name=f"{name} - Median",
            line=dict(dash="dot", color=color),
            mode="lines",

        ),
        go.Scatter(
            x=x,
            y=mean_y,
            name=f"{name} - Mean",
            line=dict(color=color),
            mode="lines",

        )
    ]
    a_color = css_color_to_rgba(color, 0.4)

    errors = [
        go.Scatter(
            x=qx,
            y=quantiles,
            marker=dict(color=color),
            name=f"{name} Q.25-Q.75",
            fill="toself",
            fillcolor=a_color,
            line=dict(color='rgba(255,255,255,0)')
        )
    ]
    return traces, errors


def plot_coverage(design, contig, start, stop, extend: int = 50, strand: str = None, colors = None, mode: str = "replicate", gff: pd.DataFrame = None, gff_types = None, gff_name: str = "locus_tag", type_colors = None, **kwargs):
    if colors is None:
        colors = DEFAULT_PLOTLY_COLORS
    if type_colors is None:
        type_colors = GFF3_TYPE_TO_COLOR
    else:
        type_colors = GFF3_TYPE_TO_COLOR | type_colors
    sf = True if "Size factors" in design else False

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, **kwargs)
    wstart = start - extend
    wend = stop + extend
    gff = filter_gff_by_interval(gff, contig, wstart, wend)
    if gff_types:
        gff = gff[gff["type"].isin(gff_types)]
    if strand:
        gff = gff[gff["strand"] == strand]

    color_idx = {treat: idx for idx, treat in enumerate(design["Treatment"].unique())}
    if mode == "replicate":
        for idx, row in design.iterrows():
            name = f"{row['Treatment']} - {row['Replicate']}"
            color = colors[color_idx[row["Treatment"]]]
            size_factor = row["Size factor"] if sf else None
            trace = get_trace(row["File"], name=name, contig=contig, start=wstart, stop=wend, strand=strand, color=color, size_factor=size_factor)
            fig.add_trace(trace)
    else:
        if sf:
            tmp = design.groupby(["Treatment"], as_index=False, observed=False).agg({"File": list, "Size factors": list}
                ).dropna().reset_index()
        else:
            tmp = design.groupby(["Treatment"], as_index=False, observed=False)["File"].agg(
                list).dropna().reset_index()

        t = []
        e = []
        for idx, row in tmp.iterrows():
            name = row["Treatment"]
            color = colors[color_idx[row["Treatment"]]]
            sfs = row["Size factors"] if sf else None

            traces, errors = get_summary_trace(row["File"], name=name, contig=contig,  start=wstart, stop=wend, strand=strand, color=color, size_factors=sfs)
            t += traces
            e += errors
        fig.add_traces(
            e+t,
            rows=1,
            cols=1
        )
    gff = gff[~gff[gff_name].isna()]
    indices = {name: idx for idx, name in enumerate(gff[gff_name].unique())}
    for (_, row) in gff.iterrows():
        idx = indices[row[gff_name]]
        v = 0.4
        arrow = 0.01
        if row["strand"] == "+":
            y = [idx - v, idx-v, idx, idx+v, idx+v, idx - v]
            pe = row["end"] - arrow * (wend - wstart)
            x = [row["start"], pe, row["end"], pe, row["start"], row["start"]]
        elif row["strand"] == "-":
            y = [idx, idx - v, idx - v, idx + v, idx + v, idx]
            ps = row["start"] + arrow * (wend - wstart)
            x = [row["start"], ps, row["end"], row["end"], ps, row["start"]]
        else:
            y = [idx - v, idx - v, idx + v, idx + v, idx - v]
            x = [row["start"], row["end"], row["end"], row["start"], row["start"]]
        try:
            color = type_colors[row["type"]]
        except KeyError:
            color = type_colors["default"]
        fig.add_trace(
            go.Scatter(
                y=y,
                x=x,
                mode="lines",
                fillcolor=color,
                line=dict(color="black"),
                fill="toself",
                name=f'{row["type"]}-{row[gff_name]}',

            ),
            row=2,
            col=1

        )
        fig.add_annotation(
            text=f'{row["type"]}',
            x=(row["end"] - row["start"]) / 2 + row["start"],
            y=idx,
            row=2,
            col=1,
            showarrow=False,

        )



    fig.update_xaxes(range=[wstart, wend])
    fig.update_yaxes(
        tickvals=list(range(len(indices))),
        ticktext=list(indices.keys()),
        tickmode="array",
        row=2
    )
    return fig


if __name__ == '__main__':

    chr = "BA000022.2"
    start = 7179
    stop = 8311

    file = "../testData/19_-P_M_R3_S16_UMI_included.fastq.gz.bam"
    file2 = "../testData/28_-P_M_R4_S19_UMI_included.fastq.gz.bam"
    file3 = "../testData/1_-P_M_R1_S13_UMI_included.fastq.gz.bam"
    file4 = "../testData/5_-P_C_R1_S14_UMI_included.fastq.gz.bam"
    file5 = "../testData/32_-P_C_R4_S20_UMI_included.fastq.gz.bam"
    file6 = "../testData/23_-P_C_R3_S17_UMI_included.fastq.gz.bam"
    files = [
        file,
        file2,
        file3,
        file4,
        file5,
        file6,
    ]
    design = {
        "File": [],
        "Replicate": [],
        "Treatment": [],

    }
    for file in files:
        s = file.split("_")
        rep = s[3]
        treat = s[2]
        design["Replicate"].append(rep)
        design["Treatment"].append(treat)
        design["File"].append(file)

    design = pd.DataFrame(design)
    #design["Size factors"] = list(range(len(design)))
    gff = "../testData/20210217_syne_onlyUnique_withFeat.gff3"
    lt2name = "../testData/lt2gname.tsv"
    lt2name = pd.read_csv(lt2name, sep="\t")
    gff = read_gff3(gff)
    gff = gff.merge(lt2name, on="locus_tag")
    sl = "psaF"
    d = find_gene_coordinates(gff, gene=sl, anno_type="gene_name")
    #d = d.iloc[1]
    d = d[d["type"] == "CDS"]
    for _, row in d.iterrows():
        _, start, stop, strand = coords_from_gff_row(row)
        fig = plot_coverage(design, chr, start, stop, extend=1000, strand=strand, mode="f", gff=gff, gff_types=None,
                            gff_name="gene_name",
                            colors=None,  # ("rgb(0, 142, 151)", "rgb(252, 76, 2)"),
                            type_colors={"CDS": "rgb(252, 76, 2)"},
                            vertical_spacing=0)
        fig.show()
        break
