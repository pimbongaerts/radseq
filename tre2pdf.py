#!/usr/bin/env python
"""
tre2pdf.py: Generates a single-page PDF visualization of a phylogenetic tree.
Tip labels are color-coded by a specified underscore-separated field in the
sample name (default: last field).
"""
import sys
import os
import glob
import argparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import Bio.Phylo

__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2016 Pim Bongaerts"
__license__ = "GPL"


def get_group(name, field_index):
    """Extract group from underscore-separated sample name"""
    parts = name.split("_")
    try:
        return parts[field_index]
    except IndexError:
        return name


def assign_colors(tip_groups):
    """Assign distinct colors to each unique group"""
    unique_groups = sorted(set(tip_groups.values()))
    n = len(unique_groups)
    if n <= 10:
        cmap = plt.cm.tab10
    elif n <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.gist_rainbow
    colors = {g: cmap(i / max(n - 1, 1)) for i, g in enumerate(unique_groups)}
    return colors


def main(tree_filename, output_filename=None, color_field=-1, tree_format=None):
    if tree_format is None:
        ext = os.path.splitext(tree_filename)[1].lower()
        fmt_map = {
            ".tre": "nexus",
            ".nex": "nexus",
            ".nexus": "nexus",
            ".nwk": "newick",
            ".newick": "newick",
            ".tree": "newick",
        }
        tree_format = fmt_map.get(ext, "nexus")

    if output_filename is None:
        output_filename = os.path.splitext(tree_filename)[0] + ".pdf"

    tree = Bio.Phylo.read(tree_filename, tree_format)
    tree.root_at_midpoint()
    tree.ladderize()

    # Hide internal node labels
    for clade in tree.get_nonterminals():
        clade.name = None

    tips = [t for t in tree.get_terminals() if t.name]
    tip_groups = {t.name: get_group(t.name, color_field) for t in tips}
    group_colors = assign_colors(tip_groups)
    label_colors = {name: group_colors[g] for name, g in tip_groups.items()}

    n_tips = len(tips)
    fig_height = max(10, n_tips * 0.35)
    fig_width = max(12, fig_height * 0.7)
    font_size = max(8, min(14, 800 / n_tips))
    marker_size = max(10, min(20, 1200 / n_tips))

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

    Bio.Phylo.draw(
        tree,
        axes=ax,
        label_colors=label_colors,
        label_func=lambda c: '  ' + c.name if c.name else '',
        do_show=False,
    )

    # Add colored circles at tip positions
    for text in ax.texts:
        label = text.get_text().strip()
        if label in label_colors:
            x, y = text.get_position()
            ax.plot(x, y, 'o', color=label_colors[label],
                    markersize=marker_size, markeredgewidth=0.5,
                    markeredgecolor='black', zorder=5,
                    clip_on=False)
            text.set_fontsize(font_size)

    ax.set_ylabel('')
    ax.set_xlabel('Distance')

    unique_groups = sorted(set(tip_groups.values()))
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', label=g,
                   markerfacecolor=group_colors[g], markersize=8,
                   markeredgewidth=0.5, markeredgecolor='black')
        for g in unique_groups
    ]
    ax.legend(handles=legend_handles, loc='upper left',
              fontsize=font_size, title='Group', framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.close(fig)

    sys.stderr.write("Tree PDF written to {}\n".format(output_filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "tree_files",
        metavar="tree_file",
        nargs="+",
        help="input tree file(s); supports glob patterns " '(e.g. "*.tre")',
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="output PDF filename (only valid for a single "
        "input file; default: <tree_file>.pdf)",
    )
    parser.add_argument(
        "--color_field",
        type=int,
        default=-1,
        help="which underscore-separated field in the sample "
        "name to use for color coding; supports negative "
        "indexing (default: -1, i.e. last field)",
    )
    parser.add_argument(
        "--format",
        choices=["nexus", "newick"],
        default=None,
        help="tree file format (default: auto-detect from " "extension)",
    )
    args = parser.parse_args()

    files = []
    for pattern in args.tree_files:
        expanded = glob.glob(pattern)
        files.extend(sorted(expanded) if expanded else [pattern])

    if args.output and len(files) > 1:
        sys.exit("Error: --output cannot be used with multiple input files")

    for tree_file in files:
        main(
            tree_file,
            output_filename=args.output,
            color_field=args.color_field,
            tree_format=args.format,
        )
