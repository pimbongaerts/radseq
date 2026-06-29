#!/usr/bin/env python
"""
Explores the deeper genetic clusters (lineages) in a `.vcf` dataset in a way
that is NOT biased by the presence of clones. Where `vcf_clone_detect.py` asks
"which individuals are near-identical (clones)", this script asks "what are the
genetic groups, how many of them are there, and how differentiated are they".

Clones are non-independent samples that distort per-group allele frequencies and
therefore bias essentially every population-genetic measure (Fst, private and
fixed-private alleles, heterozygosity, ordination, and even the clustering
itself). By default the script does NOT remove clones - it runs on all
individuals as given (use this when the input is already clone-corrected, or to
inspect the raw structure). Optional CLONE-CORRECTION reduces each clonal genet
to a single representative ramet (the one with the least missing data) and runs
the entire analysis - similarity, tree, K-evaluation and all differentiation
statistics - on that clone-corrected (genet) set. Enable it with `--auto-clone`
(clonal genets detected internally by reusing `vcf_clone_detect.py` at an auto-
inferred threshold), `--clone-threshold PCT` (manual threshold), or
`--clone-list FILE` (external list of samples to drop, e.g. the "individuals to
remove" output of `vcf_clone_detect.py`).

The script (1) computes pairwise genetic similarities (`--method`, default
`ibs`; `dosage` is the documented alternative), (2) builds a UPGMA / average-
linkage tree, (3) sweeps K = 2, 3, ... determining for each K the genetic-
similarity cut-off that splits the tree into K groups and a separation gap,
stopping once a clean split can no longer be made (or `--max-k` is reached),
(4) assigns every genet to a cluster at each K, and (5) produces a multi-panel
PDF: on the left a tree with per-K cluster-assignment columns and a % genotyped
bar aligned to the tips, and on the right a set of differentiation views
(cluster-combination similarity histogram, private / fixed-private alleles both
per-cluster and per-pair, a pairwise Fst heatmap, a fixed-difference matrix, a
PCA/PCoA scatter and a per-cluster diversity bar).

The UPGMA linkage is the single source of truth: it is drawn as the tree AND
cut to give every K-assignment, so the tree and the columns are always coherent.
With `--tree nj` a neighbour-joining tree (as in `vcf_clone_detect.py`) is drawn
for display instead, but the cluster assignments still come from UPGMA.

Example:
  python3 vcf_cluster_explore.py --vcf vcf_file.vcf --pop pop_file.txt \
      --output clusters.csv
"""
import sys
import os
import argparse
import numpy as np

import vcf_clone_detect

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2024 Pim Bongaerts'
__license__ = 'GPL'


DEFAULT_METHOD = 'ibs'
DEF_MAX_K = 10
DEF_MIN_CLUSTER_SIZE = 2
DEF_MIN_GAP = 0.05
CLONE_FLOOR = vcf_clone_detect.DEF_THRESHOLD   # only pairs >= this can be clones

CSV_SUFFIX = '_clusters.csv'
PDF_SUFFIX = '_clusters.pdf'


# --------------------------------------------------------------------------- #
#  Clone-correction (the defining step)
# --------------------------------------------------------------------------- #
def read_clone_list(filename):
    """ Read a list of sample names to drop (one per line, '#' comments). """
    remove = set()
    with open(filename, 'r') as handle:
        for line in handle:
            name = line.strip()
            if name and not name.startswith('#'):
                remove.add(name.split()[0])
    return remove


def infer_clone_threshold(comparisons):
    """ Quietly infer a clonal similarity threshold from the largest gap among
    the highest matches (mirrors `vcf_clone_detect.output_highest_matches`).
    Returns None when there is no distinct high-similarity (clonal) signal. """
    last = highest_diff = max_p = min_p = 0.0
    for row in np.nditer(comparisons):
        match_perc = float(row[vcf_clone_detect.C_MATCH_PERC])
        if match_perc < CLONE_FLOOR:
            break
        if last != 0.0:
            diff = round(last - match_perc, 2)
            if diff > highest_diff:
                highest_diff = diff
                min_p = match_perc
                max_p = last
        last = match_perc
    if max_p == 0.0:
        return None
    if int(max_p) > min_p:
        return float(int(max_p))
    return (max_p + min_p) / 2.0


def subset_data(data, keep_idx):
    """ Return a copy of the load_vcf_arrays dict restricted to sample columns
    in `keep_idx` (sample-indexed arrays sliced; site-level fields unchanged). """
    new = dict(data)
    new['names'] = [data['names'][i] for i in keep_idx]
    new['state_code'] = data['state_code'][:, keep_idx]
    new['dosage'] = data['dosage'][:, keep_idx]
    if data.get('ad_p') is not None:
        new['ad_p'] = data['ad_p'][:, keep_idx]
    if data.get('ad_present') is not None:
        new['ad_present'] = data['ad_present'][:, keep_idx]
    return new


def clone_correct(data, method, clone_list_file, clone_threshold, auto_clone):
    """ Optionally reduce each clonal genet to one representative ramet. Clone-
    correction only happens with --clone-list, --clone-threshold or --auto-clone;
    otherwise the data is returned unchanged. Returns (data, removed_set, info). """
    names = data['names']
    if clone_list_file:
        remove = read_clone_list(clone_list_file) & set(names)
        info = ('supplied --clone-list `{0}` ({1} of its samples present)'
                .format(clone_list_file, len(remove)))
    elif clone_threshold is not None or auto_clone:
        comparisons = vcf_clone_detect.build_comparisons(data, {}, method)
        if clone_threshold is not None:
            threshold = float(clone_threshold)
            source = 'manual --clone-threshold'
        else:
            threshold = infer_clone_threshold(comparisons)
            source = 'auto-inferred threshold'
        if threshold is None:
            return data, set(), 'no distinct clonal signal detected (--auto-clone)'
        groups = vcf_clone_detect.cluster_clones(comparisons, threshold)
        remove = set()
        for group in groups:
            remove.update(group.get_samples_to_remove())
        info = ('{0} = {1}%, {2} clonal group(s), {3} ramet(s) removed'
                .format(source, round(threshold, 2), len(groups), len(remove)))
    else:
        return data, set(), ('clone-correction disabled (default; enable with '
                             '--auto-clone, --clone-threshold or --clone-list)')

    if not remove:
        return data, set(), info
    keep_idx = [i for i, nm in enumerate(names) if nm not in remove]
    return subset_data(data, keep_idx), remove, info


# --------------------------------------------------------------------------- #
#  Similarity -> distance -> UPGMA linkage
# --------------------------------------------------------------------------- #
def compute_distance_matrix(data, method):
    """ Build the NxN similarity (%) and distance (1 - sim/100) matrices.
    Returns (sim, dist, n_no_overlap). """
    match_mat, denom_mat = vcf_clone_detect._match_denom_matrices(data, method)
    with np.errstate(divide='ignore', invalid='ignore'):
        sim = np.where(denom_mat > 0, 100.0 * match_mat / denom_mat, np.nan)
    np.fill_diagonal(sim, 100.0)
    dist = 1.0 - sim / 100.0
    # Pairs with no shared genotyped sites are NaN; scipy rejects those, so set
    # them to the maximum distance (1.0) and count how many were affected.
    nan_mask = ~np.isfinite(dist)
    np.fill_diagonal(nan_mask, False)
    n_no_overlap = int(nan_mask.sum() // 2)
    dist[~np.isfinite(dist)] = 1.0
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0          # enforce exact symmetry for squareform
    np.fill_diagonal(dist, 0.0)
    return sim, dist, n_no_overlap


def build_linkage(dist):
    """ UPGMA / average-linkage condensed-distance clustering (scipy). """
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage
    return linkage(squareform(dist, checks=False), method='average')


# --------------------------------------------------------------------------- #
#  K evaluation (cut-off + separation gap per K, with stopping rule)
# --------------------------------------------------------------------------- #
def evaluate_k(linkage_matrix, n, max_k, min_cluster_size, min_gap):
    """ For each K from 2 upward determine the similarity cut-off that yields K
    clusters and the separation gap supporting it, stopping once a clean split
    can no longer be made. Returns (records, k_stop, best_k).

    Indexing: `Z[:,2]` are the n-1 merge heights ascending. K clusters occur
    for a cut height in (Z[n-1-K, 2], Z[n-K, 2]); the gap supporting K is the
    width of that interval. """
    from scipy.cluster.hierarchy import fcluster
    heights = linkage_matrix[:, 2]
    root_height = heights[n - 2] if n >= 2 else 1.0
    kmax = min(max_k, n - 1)
    records = []
    for k in range(2, kmax + 1):
        h_low = heights[n - 1 - k]
        h_high = heights[n - k]
        cutoff_sim = 100.0 * (1.0 - 0.5 * (h_low + h_high))
        abs_gap = h_high - h_low
        rel_gap = abs_gap / root_height if root_height > 0 else 0.0
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        sizes = np.bincount(labels)[1:]
        records.append({'K': k, 'labels': labels, 'cutoff': cutoff_sim,
                        'abs_gap': abs_gap, 'rel_gap': rel_gap,
                        'min_size': int(sizes.min()),
                        'n_clusters': int((sizes > 0).sum())})
        # Stop scanning once this split is no longer clean (record it, then halt)
        if rel_gap < min_gap or int(sizes.min()) < min_cluster_size:
            break

    clean = [r for r in records
             if r['rel_gap'] >= min_gap and r['min_size'] >= min_cluster_size]
    k_stop = max((r['K'] for r in clean), default=2)
    candidates = [r for r in records if r['K'] <= k_stop]
    best_k = max(candidates, key=lambda r: r['abs_gap'])['K'] if candidates else 2
    return records, k_stop, best_k


def select_record(records, k):
    """ Return the record for a given K (or the closest available). """
    for record in records:
        if record['K'] == k:
            return record
    return records[0]


# --------------------------------------------------------------------------- #
#  Per-cluster / pairwise differentiation statistics (clone-corrected data)
# --------------------------------------------------------------------------- #
def cluster_allele_stats(dosage, labels):
    """ Per-cluster x per-site allele summaries from the alt-dosage matrix. """
    clusters = sorted(set(int(x) for x in labels))
    n_clusters = len(clusters)
    n_sites = dosage.shape[0]
    freq = np.full((n_clusters, n_sites), np.nan)
    ncall = np.zeros((n_clusters, n_sites))
    altcnt = np.zeros((n_clusters, n_sites))
    for ci, cluster in enumerate(clusters):
        sub = dosage[:, labels == cluster]
        valid = sub >= 0
        nseg = valid.sum(axis=1)
        alt = np.where(valid, sub, 0).sum(axis=1)
        ncall[ci] = nseg
        altcnt[ci] = alt
        seg = nseg > 0
        freq[ci, seg] = alt[seg] / (2.0 * nseg[seg])
    return {'clusters': clusters, 'n_clusters': n_clusters, 'freq': freq,
            'ncall': ncall, 'altcnt': altcnt,
            'alt_present': altcnt > 0,
            'ref_present': (2 * ncall - altcnt) > 0,
            'has_data': ncall > 0}


def private_alleles(stats):
    """ Per-cluster private alleles (present in this cluster, absent in all
    others). Returns list of (total_private, fixed_private) per cluster. """
    alt_present, ref_present = stats['alt_present'], stats['ref_present']
    freq = stats['freq']
    total_alt = alt_present.sum(axis=0)
    total_ref = ref_present.sum(axis=0)
    out = []
    for ci in range(stats['n_clusters']):
        priv_alt = alt_present[ci] & ((total_alt - alt_present[ci]) == 0)
        priv_ref = ref_present[ci] & ((total_ref - ref_present[ci]) == 0)
        total = int(priv_alt.sum() + priv_ref.sum())
        fixed = int((priv_alt & (freq[ci] == 1.0)).sum() +
                    (priv_ref & (freq[ci] == 0.0)).sum())
        out.append((total, fixed))
    return out


def pairwise_matrices(stats):
    """ Pairwise fixed-difference counts and directional private-allele counts.
    Returns (fixed_diff[KxK], priv_pair[KxK]) where priv_pair[a, b] = alleles
    present in cluster a but absent in cluster b. """
    n_clusters = stats['n_clusters']
    freq, has_data = stats['freq'], stats['has_data']
    alt_present, ref_present = stats['alt_present'], stats['ref_present']
    fixed_diff = np.zeros((n_clusters, n_clusters), dtype=int)
    priv_pair = np.zeros((n_clusters, n_clusters), dtype=int)
    for a in range(n_clusters):
        for b in range(n_clusters):
            if a == b:
                continue
            valid = has_data[a] & has_data[b]
            if a < b:
                fdiff = int((valid & (np.abs(freq[a] - freq[b]) == 1.0)).sum())
                fixed_diff[a, b] = fixed_diff[b, a] = fdiff
            priv_pair[a, b] = int(
                (alt_present[a] & ~alt_present[b] & has_data[b]).sum() +
                (ref_present[a] & ~ref_present[b] & has_data[b]).sum())
    return fixed_diff, priv_pair


def pairwise_fst(stats):
    """ Hudson's Fst between every cluster pair (ratio of averages, allele-copy
    counts so a single diploid still gives n-1 = 1). Returns KxK (NaN diag). """
    n_clusters = stats['n_clusters']
    freq, ncall = stats['freq'], stats['ncall']
    fst = np.full((n_clusters, n_clusters), np.nan)
    for a in range(n_clusters):
        for b in range(a + 1, n_clusters):
            pa, pb = freq[a], freq[b]
            na, nb = 2.0 * ncall[a], 2.0 * ncall[b]
            valid = (na >= 2) & (nb >= 2) & np.isfinite(pa) & np.isfinite(pb)
            with np.errstate(invalid='ignore', divide='ignore'):
                num = ((pa - pb) ** 2
                       - pa * (1 - pa) / np.where(na > 1, na - 1, np.nan)
                       - pb * (1 - pb) / np.where(nb > 1, nb - 1, np.nan))
                den = pa * (1 - pb) + pb * (1 - pa)
            use = valid & np.isfinite(num) & np.isfinite(den)
            den_sum = den[use].sum()
            value = num[use].sum() / den_sum if den_sum > 0 else np.nan
            fst[a, b] = fst[b, a] = value
    return fst


def cluster_diversity(state_code, state_is_hom, labels, stats):
    """ Per-cluster observed heterozygosity and % polymorphic sites. """
    present = state_code >= 0
    safe = np.where(present, state_code, 0)
    is_het = present & ~state_is_hom[safe]
    het_obs, poly = [], []
    for ci, cluster in enumerate(stats['clusters']):
        member = labels == cluster
        denom = present[:, member].sum()
        het_obs.append(is_het[:, member].sum() / denom if denom > 0 else np.nan)
        freq = stats['freq'][ci]
        segregating = (freq > 0) & (freq < 1)
        has_data = stats['has_data'][ci].sum()
        poly.append(100.0 * segregating.sum() / has_data if has_data else np.nan)
    return het_obs, poly


def ordinate(dosage, dist, kind):
    """ 2D ordination of individuals. Returns (coords[n, 2], pct_var[2]). """
    if kind == 'pca':
        geno = dosage.T.astype(float)            # samples x sites
        missing = geno < 0
        valid = ~missing
        col_n = valid.sum(axis=0)
        col_mean = np.where(col_n > 0,
                            np.where(valid, geno, 0).sum(axis=0)
                            / np.maximum(col_n, 1), 0.0)
        geno = np.where(missing, col_mean, geno)
        geno = geno - geno.mean(axis=0)
        u_mat, sing, _ = np.linalg.svd(geno, full_matrices=False)
        coords = u_mat[:, :2] * sing[:2]
        total = (sing ** 2).sum()
        pct = 100.0 * sing[:2] ** 2 / total if total > 0 else np.zeros(2)
        return coords, pct
    # PCoA on the distance matrix (coherent with the tree, but D is non-Euclidean)
    n = dist.shape[0]
    centering = np.eye(n) - 1.0 / n
    gram = -0.5 * centering @ (dist ** 2) @ centering
    eigval, eigvec = np.linalg.eigh(gram)
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    coords = eigvec[:, :2] * np.sqrt(np.clip(eigval[:2], 0, None))
    positive = np.clip(eigval, 0, None).sum()
    pct = (100.0 * np.clip(eigval[:2], 0, None) / positive
           if positive > 0 else np.zeros(2))
    return coords, pct


# --------------------------------------------------------------------------- #
#  Output: assignment CSV, filenames
# --------------------------------------------------------------------------- #
def write_assignment_csv(names, records, filename):
    """ Write per-sample cluster assignment across all evaluated K (popfile-like:
    any single `K*` column is a usable 2-column popfile). """
    headers = ['sample'] + ['K{0}'.format(r['K']) for r in records]
    with open(filename, 'w') as handle:
        handle.write(','.join(headers) + '\n')
        for i, name in enumerate(names):
            row = [name] + [str(int(r['labels'][i])) for r in records]
            handle.write(','.join(row) + '\n')
    print('Cluster assignments outputted to file: `{0}`'.format(filename))


def derive_outputs(output_filename, vcf_filename, pdf_output):
    """ Decide the CSV and PDF output filenames. With `-o foo.csv` the PDF
    mirrors it (`foo.pdf`); otherwise both derive from the vcf basename. """
    if output_filename:
        csv_filename = output_filename
        pdf_base = os.path.splitext(output_filename)[0]
    else:
        base = os.path.splitext(vcf_filename)[0] if vcf_filename \
            else 'cluster_explore'
        csv_filename = base + CSV_SUFFIX
        pdf_base = base + '_clusters'
    pdf_filename = pdf_output if pdf_output else pdf_base + '.pdf'
    return csv_filename, pdf_filename


# --------------------------------------------------------------------------- #
#  PDF report
# --------------------------------------------------------------------------- #
def cluster_palette(n):
    """ A categorical colour per cluster id (1..n). """
    import matplotlib.pyplot as plt
    if n <= 10:
        cmap = plt.cm.tab10
        cols = [cmap(i) for i in range(n)]
    elif n <= 20:
        cmap = plt.cm.tab20
        cols = [cmap(i) for i in range(n)]
    else:
        cmap = plt.cm.gist_rainbow
        cols = [cmap(i / max(n - 1, 1)) for i in range(n)]
    return cols


def build_upgma_tree(linkage_matrix, names):
    """ Convert the scipy UPGMA linkage into a Bio.Phylo tree (so the drawn tree
    is exactly the object that is cut for the K-assignments). """
    import Bio.Phylo.BaseTree as BaseTree
    from scipy.cluster.hierarchy import to_tree
    root = to_tree(linkage_matrix, rd=False)

    def build(node):
        if node.is_leaf():
            return BaseTree.Clade(name=names[node.id], branch_length=0.0)
        left, right = build(node.get_left()), build(node.get_right())
        left.branch_length = node.dist - node.get_left().dist
        right.branch_length = node.dist - node.get_right().dist
        return BaseTree.Clade(clades=[left, right])

    tree = BaseTree.Tree(root=build(root))
    tree.ladderize()
    return tree


def build_nj_tree(dist, names):
    """ Neighbour-joining tree for display (as in vcf_clone_detect). """
    import Bio.Phylo.TreeConstruction as TC
    lower = [[dist[i][j] for j in range(i + 1)] for i in range(len(names))]
    dmatrix = TC._DistanceMatrix(names=list(names), matrix=lower)
    tree = TC.DistanceTreeConstructor().nj(dmatrix)
    tree.ladderize()
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree


def write_pdf_report(sim, dist, linkage_matrix, names, records, col_max,
                     selected, perc_genotyped, data, method, tree_mode,
                     ordination, pdf_filename):
    """ Multi-panel PDF: left = tree + per-K cluster columns + %genotyped bar;
    right = cluster-combination histogram, private/fixed alleles (per-cluster &
    per-pair), Fst heatmap, fixed-difference matrix, ordination, diversity. """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    import Bio.Phylo

    name_idx = {nm: i for i, nm in enumerate(names)}
    n_tips = len(names)
    display = [r for r in records if r['K'] <= col_max]
    n_cols = len(display)

    # Selected-K clustering drives the right-hand panels and the tip colours
    labels = selected['labels']
    stats = cluster_allele_stats(data['dosage'], labels)
    clusters = stats['clusters']
    sel_palette = cluster_palette(stats['n_clusters'])
    cl_index = {c: i for i, c in enumerate(clusters)}
    sel_color = {nm: sel_palette[cl_index[int(labels[name_idx[nm]])]]
                 for nm in names}

    # --- canvas geometry (inches) -> figure fractions via ax_in() ------------
    left_margin, right_margin = 0.6, 0.6
    top_margin, bottom_margin = 0.5, 0.5
    tree_w = 4.2
    kcol_w = max(0.9, 0.22 * n_cols)
    bar_w = 1.0
    right_w = 6.6
    gap = 0.5
    per_tip = 0.13

    left_band_w = (left_margin + tree_w + 0.1 + kcol_w + 0.15 + bar_w)
    fig_w = left_band_w + gap + right_w + right_margin
    x_right = left_band_w + gap

    panel_h = [1.8, 1.9, 1.7, 1.7, 2.0, 1.6]      # the 6 right-hand panels
    right_needed = top_margin + sum(panel_h) + 0.5 * len(panel_h) + bottom_margin
    tree_h = max(3.0, n_tips * per_tip)
    left_needed = top_margin + tree_h + bottom_margin
    fig_h = max(left_needed, right_needed)

    fig = plt.figure(figsize=(fig_w, fig_h))

    def ax_in(x_in, ytop_in, w_in, h_in):
        return fig.add_axes([x_in / fig_w,
                             (fig_h - ytop_in - h_in) / fig_h,
                             w_in / fig_w, h_in / fig_h])

    # --- tree ----------------------------------------------------------------
    ax_tree = ax_in(left_margin, top_margin, tree_w, tree_h)
    if tree_mode == 'nj':
        tree = build_nj_tree(dist, names)
    else:
        tree = build_upgma_tree(linkage_matrix, names)
    sys.stderr.write('Drawing {0} tree ({1} tips)...\n'.format(tree_mode, n_tips))
    font_size = max(1.0, min(8, 500.0 / n_tips))
    Bio.Phylo.draw(tree, axes=ax_tree, do_show=False,
                   label_func=lambda c: '  ' + c.name if c.name else '')
    ax_tree.set_ylabel('')
    ax_tree.set_xlabel('Genetic distance (1 - {0} similarity)'.format(method),
                       fontsize=7)
    ax_tree.tick_params(labelsize=6)

    tip_y = {}
    for text in ax_tree.texts:
        label = text.get_text().strip()
        if label in name_idx:
            tip_y[label] = text.get_position()[1]
            text.set_fontsize(font_size)
            text.set_color(sel_color[label])
    if not tip_y:
        sys.stderr.write('Warning: could not recover tip positions.\n')
        plt.close(fig)
        return
    ys = sorted(tip_y.values())
    row_h = (ys[1] - ys[0]) if len(ys) > 1 else 1.0
    if tree_mode == 'nj':
        title = '{0} tree (display); clusters = UPGMA'.format(tree_mode.upper())
    else:
        title = 'UPGMA tree (clusters = tree cuts)'
    ax_tree.set_title(title, fontsize=8)

    # --- per-K cluster-assignment columns ------------------------------------
    ax_k = ax_in(left_margin + tree_w + 0.1, top_margin, kcol_w, tree_h)
    rects, colors = [], []
    for col, record in enumerate(display):
        palette = (sel_palette if record['K'] == selected['K']
                   else cluster_palette(record['n_clusters']))
        idx = {c: i for i, c in enumerate(sorted(set(int(x)
               for x in record['labels'])))}
        for nm, y in tip_y.items():
            cid = int(record['labels'][name_idx[nm]])
            rects.append(Rectangle((col, y - row_h * 0.45), 1, row_h * 0.9))
            colors.append(palette[idx[cid]])
    ax_k.add_collection(PatchCollection(rects, facecolors=colors,
                                        edgecolors='none'))
    ax_k.set_xlim(0, max(n_cols, 1))
    ax_k.set_ylim(ax_tree.get_ylim())
    ax_k.set_yticks([])
    ax_k.set_xticks([c + 0.5 for c in range(n_cols)])
    ax_k.set_xticklabels(
        ['K={0}{1}'.format(r['K'], '*' if r['K'] == selected['K'] else '')
         for r in display], fontsize=6, rotation=90)
    ax_k.tick_params(length=0)
    ax_k.set_title('cluster @ K', fontsize=7)

    # --- % genotyped bar -----------------------------------------------------
    ax_bar = ax_in(left_margin + tree_w + 0.1 + kcol_w + 0.15, top_margin,
                   bar_w, tree_h)
    bar_names = [nm for nm in names if nm in tip_y]
    bar_y = [tip_y[nm] for nm in bar_names]
    bar_v = [perc_genotyped[nm] for nm in bar_names]
    ax_bar.barh(bar_y, bar_v, height=row_h * 0.8,
                color=[sel_color[nm] for nm in bar_names], edgecolor='none')
    ax_bar.set_ylim(ax_tree.get_ylim())
    ax_bar.set_yticks([])
    ax_bar.set_xlim(max(0, min(bar_v) - 5), 100)
    ax_bar.set_title('% genotyped', fontsize=7)
    ax_bar.tick_params(labelsize=6)
    for spine in ('top', 'right'):
        ax_bar.spines[spine].set_visible(False)

    cluster_labels = ['C{0}'.format(c) for c in clusters]

    # --- (a) cluster-combination similarity histogram ------------------------
    y_cur = top_margin
    ax_hist = ax_in(x_right, y_cur, right_w, panel_h[0])
    triu = np.triu_indices(n_tips, 1)
    sim_vals = sim[triu]
    li, lj = labels[triu[0]], labels[triu[1]]
    finite = np.isfinite(sim_vals)
    within = li == lj
    lo = np.floor(np.nanmin(sim_vals[finite])) if finite.any() else 0
    bins = np.linspace(lo, 100, 40)
    combos = stats['n_clusters'] * (stats['n_clusters'] + 1) // 2
    if combos <= 8:
        for ci, c in enumerate(clusters):
            mask = finite & within & (li == c)
            if mask.any():
                ax_hist.hist(sim_vals[mask], bins=bins, alpha=0.55,
                             color=sel_palette[ci],
                             label='within C{0}'.format(c))
        between_cmap = plt.cm.Greys
        pairs = [(a, b) for ai, a in enumerate(clusters)
                 for b in clusters[ai + 1:]]
        for pi, (a, b) in enumerate(pairs):
            mask = finite & (((li == a) & (lj == b)) | ((li == b) & (lj == a)))
            if mask.any():
                ax_hist.hist(sim_vals[mask], bins=bins, histtype='step',
                             linewidth=1.2,
                             color=between_cmap(0.4 + 0.5 * pi / max(len(pairs)
                                                - 1, 1)),
                             label='C{0}-C{1}'.format(a, b))
    else:
        ax_hist.hist(sim_vals[finite & ~within], bins=bins, color='0.6',
                     label='between clusters')
        ax_hist.hist(sim_vals[finite & within], bins=bins, color='#1f6f6f',
                     alpha=0.8, label='within clusters')
    ax_hist.axvline(selected['cutoff'], color='red', ls='--', lw=1)
    ax_hist.set_title('Similarity by cluster combination (K={0})'
                      .format(selected['K']), fontsize=8)
    ax_hist.set_xlabel('Genetic similarity (%)', fontsize=7)
    ax_hist.set_ylabel('Pairs', fontsize=7)
    ax_hist.legend(fontsize=5, ncol=2, framealpha=0.9)
    ax_hist.tick_params(labelsize=6)
    for spine in ('top', 'right'):
        ax_hist.spines[spine].set_visible(False)

    # --- (b) private + fixed-private alleles (per-cluster & per-pair) ---------
    y_cur += panel_h[0] + 0.5
    priv = private_alleles(stats)
    fixed_diff, priv_pair = pairwise_matrices(stats)
    ax_pc = ax_in(x_right, y_cur, right_w * 0.46, panel_h[1])
    seg = [t - f for t, f in priv]
    fix = [f for _, f in priv]
    xpos = np.arange(stats['n_clusters'])
    ax_pc.bar(xpos, seg, color=sel_palette, label='private (segregating)')
    ax_pc.bar(xpos, fix, bottom=seg, color=sel_palette, hatch='//',
              edgecolor='black', linewidth=0.3, label='private (fixed)')
    ax_pc.set_xticks(xpos)
    ax_pc.set_xticklabels(cluster_labels, fontsize=6, rotation=90)
    ax_pc.set_title('Private alleles per cluster', fontsize=8)
    ax_pc.set_ylabel('# alleles', fontsize=7)
    ax_pc.legend(fontsize=5)
    ax_pc.tick_params(labelsize=6)

    ax_pp = ax_in(x_right + right_w * 0.54, y_cur, right_w * 0.46, panel_h[1])
    pairs = [(a, b) for ai, a in enumerate(range(stats['n_clusters']))
             for b in range(ai + 1, stats['n_clusters'])]
    if pairs:
        xpp = np.arange(len(pairs))
        priv_a = [priv_pair[a, b] for a, b in pairs]
        priv_b = [priv_pair[b, a] for a, b in pairs]
        ax_pp.bar(xpp, priv_a, color=[sel_palette[a] for a, _ in pairs],
                  label='private to first')
        ax_pp.bar(xpp, priv_b, bottom=priv_a,
                  color=[sel_palette[b] for _, b in pairs], alpha=0.6,
                  label='private to second')
        ax_pp.set_xticks(xpp)
        ax_pp.set_xticklabels(['C{0}-C{1}'.format(clusters[a], clusters[b])
                               for a, b in pairs], fontsize=6, rotation=90)
    ax_pp.set_title('Private alleles per cluster pair', fontsize=8)
    ax_pp.tick_params(labelsize=6)

    # --- (c) pairwise Fst heatmap --------------------------------------------
    y_cur += panel_h[1] + 0.5
    ax_fst = ax_in(x_right, y_cur, right_w * 0.46, panel_h[2])
    fst = pairwise_fst(stats)
    disp = np.where(np.isfinite(fst), np.clip(fst, 0, None), np.nan)
    image = ax_fst.imshow(disp, cmap='viridis')
    ax_fst.set_xticks(range(stats['n_clusters']))
    ax_fst.set_yticks(range(stats['n_clusters']))
    ax_fst.set_xticklabels(cluster_labels, fontsize=6, rotation=90)
    ax_fst.set_yticklabels(cluster_labels, fontsize=6)
    for a in range(stats['n_clusters']):
        for b in range(stats['n_clusters']):
            if np.isfinite(fst[a, b]):
                ax_fst.text(b, a, '{0:.2f}'.format(fst[a, b]), ha='center',
                            va='center', fontsize=5,
                            color='white' if disp[a, b] > np.nanmax(disp) / 2
                            else 'black')
    ax_fst.set_title('Pairwise Fst (Hudson)', fontsize=8)
    fig.colorbar(image, ax=ax_fst, fraction=0.046, pad=0.04)

    # --- (d) fixed-difference matrix -----------------------------------------
    ax_fd = ax_in(x_right + right_w * 0.54, y_cur, right_w * 0.46, panel_h[2])
    fd_disp = fixed_diff.astype(float)
    np.fill_diagonal(fd_disp, np.nan)
    image2 = ax_fd.imshow(fd_disp, cmap='magma')
    ax_fd.set_xticks(range(stats['n_clusters']))
    ax_fd.set_yticks(range(stats['n_clusters']))
    ax_fd.set_xticklabels(cluster_labels, fontsize=6, rotation=90)
    ax_fd.set_yticklabels(cluster_labels, fontsize=6)
    for a in range(stats['n_clusters']):
        for b in range(stats['n_clusters']):
            if a != b:
                ax_fd.text(b, a, str(fixed_diff[a, b]), ha='center',
                           va='center', fontsize=5, color='0.9')
    ax_fd.set_title('Fixed differences', fontsize=8)
    fig.colorbar(image2, ax=ax_fd, fraction=0.046, pad=0.04)

    # --- (e) ordination scatter ----------------------------------------------
    y_cur += panel_h[2] + 0.5
    ax_ord = ax_in(x_right, y_cur, right_w, panel_h[4])
    coords, pct = ordinate(data['dosage'], dist, ordination)
    for ci, c in enumerate(clusters):
        member = labels == c
        ax_ord.scatter(coords[member, 0], coords[member, 1], s=14,
                       color=sel_palette[ci], edgecolor='none',
                       label='C{0}'.format(c))
    ax_ord.set_xlabel('{0}1 ({1:.1f}%)'.format(ordination.upper(), pct[0]),
                      fontsize=7)
    ax_ord.set_ylabel('{0}2 ({1:.1f}%)'.format(ordination.upper(), pct[1]),
                      fontsize=7)
    ax_ord.set_title('Ordination ({0})'.format(ordination.upper()), fontsize=8)
    ax_ord.legend(fontsize=5, ncol=2)
    ax_ord.tick_params(labelsize=6)
    for spine in ('top', 'right'):
        ax_ord.spines[spine].set_visible(False)

    # --- (f) per-cluster diversity -------------------------------------------
    y_cur += panel_h[4] + 0.5
    ax_div = ax_in(x_right, y_cur, right_w, panel_h[5])
    het_obs, poly = cluster_diversity(data['state_code'], data['state_is_hom'],
                                      labels, stats)
    xpos = np.arange(stats['n_clusters'])
    ax_div.bar(xpos - 0.2, het_obs, 0.4, color=sel_palette, label='Ho')
    ax_div.set_ylabel('Observed heterozygosity', fontsize=7)
    ax_div.set_xticks(xpos)
    ax_div.set_xticklabels(cluster_labels, fontsize=6, rotation=90)
    ax_div.tick_params(labelsize=6)
    ax_poly = ax_div.twinx()
    ax_poly.bar(xpos + 0.2, poly, 0.4, color='0.6', label='% polymorphic')
    ax_poly.set_ylabel('% polymorphic sites', fontsize=7)
    ax_poly.tick_params(labelsize=6)
    ax_div.set_title('Per-cluster diversity', fontsize=8)
    lines = (ax_div.get_legend_handles_labels()[0]
             + ax_poly.get_legend_handles_labels()[0])
    labs = (ax_div.get_legend_handles_labels()[1]
            + ax_poly.get_legend_handles_labels()[1])
    ax_div.legend(lines, labs, fontsize=5, loc='upper center', ncol=2,
                  framealpha=0.9)

    fig.savefig(pdf_filename, format='pdf')
    plt.close(fig)
    sys.stderr.write('PDF report written to {0}\n'.format(pdf_filename))


# --------------------------------------------------------------------------- #
#  Text-output helpers
# --------------------------------------------------------------------------- #
def print_k_table(records, k_stop, best_k, selected_k):
    """ ###3 - K evaluation table. """
    print('{0:>3}  {1:>10}  {2:>8}  {3:>8}  {4:>9}  {5:>10}  {6}'.format(
        'K', 'cutoff(%)', 'abs_gap', 'rel_gap', 'min_clus', 'n_clusters',
        'flag'))
    for record in records:
        flags = []
        if record['K'] == best_k:
            flags.append('best')
        if record['K'] == selected_k:
            flags.append('selected')
        if record['K'] > k_stop:
            flags.append('unclean')
        print('{0:>3}  {1:>10.2f}  {2:>8.4f}  {3:>8.3f}  {4:>9}  {5:>10}  {6}'
              .format(record['K'], record['cutoff'], record['abs_gap'],
                      record['rel_gap'], record['min_size'],
                      record['n_clusters'], ', '.join(flags)))


def cluster_size_map(labels, clusters):
    return {c: int((labels == c).sum()) for c in clusters}


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main(vcf_filename, pop_filename, output_filename, method, max_k,
         min_cluster_size, min_gap, force_k, tree_mode, ordination,
         clone_list, clone_threshold, auto_clone, make_pdf, pdf_output):

    print('###1 - Loading VCF (method: {0})'.format(method))
    if not vcf_filename:
        sys.exit('Error: please provide a vcf_file (--vcf).')
    need_ad = (method == 'single-read')
    data = vcf_clone_detect.load_vcf_arrays(vcf_filename, need_ad=need_ad)
    if data['multiallelic']:
        sys.stderr.write('Warning: multiallelic sites present; allele counts '
                         'use alt-dosage only.\n')
    n_total = len(data['names'])
    data, removed, info = clone_correct(data, method, clone_list,
                                        clone_threshold, auto_clone)
    n_units = len(data['names'])
    unit = 'genets' if removed else 'individuals'
    print('{0} samples -> {1} {2} ({3} removed); {4}'.format(
        n_total, n_units, unit, len(removed), info))
    if n_units < 3:
        sys.exit('Error: need at least 3 {0} for cluster exploration '
                 '(have {1}).'.format(unit, n_units))

    print('\n###2 - Distance matrix and UPGMA linkage')
    sim, dist, n_no_overlap = compute_distance_matrix(data, method)
    if n_no_overlap:
        sys.stderr.write('Warning: {0} pair(s) share no genotyped sites; '
                         'set to maximum distance.\n'.format(n_no_overlap))
    linkage_matrix = build_linkage(dist)
    print('{0} {1}, {2} loci, {3} pair(s) with no shared sites'.format(
        n_units, unit, data['n_loci'], n_no_overlap))

    print('\n###3 - K evaluation (cut-off & separation gap per K)')
    records, k_stop, best_k = evaluate_k(linkage_matrix, n_units, max_k,
                                         min_cluster_size, min_gap)
    selected_k = best_k
    if force_k:
        selected_k = int(force_k)
        if selected_k not in [r['K'] for r in records]:
            sys.stderr.write('Warning: --k {0} not available; using best K {1}.'
                             '\n'.format(selected_k, best_k))
            selected_k = best_k
    selected = select_record(records, selected_k)
    display_max = max(k_stop, selected_k)
    print_k_table(records, k_stop, best_k, selected_k)
    has_structure = any(r['rel_gap'] >= min_gap and r['min_size'] >=
                        min_cluster_size for r in records)
    if has_structure:
        print('Best-supported K = {0}; selected K = {1} (well-supported up to '
              'K = {2})'.format(best_k, selected_k, k_stop))
    else:
        print('No well-supported split (no K reached the min-gap/min-cluster '
              'criteria); showing K = {0} as the weakest default.'.format(
                  selected_k))

    print('\n###4 - Cluster membership at selected K = {0}'.format(selected_k))
    labels = selected['labels']
    clusters = sorted(set(int(x) for x in labels))
    for cluster in clusters:
        members = sorted(nm for i, nm in enumerate(data['names'])
                         if int(labels[i]) == cluster)
        print('C{0} ({1}): {2}'.format(cluster, len(members),
                                       ', '.join(members)))
    if pop_filename:
        indivs_pops = vcf_clone_detect.get_pop_assignments_from_popfile(
            pop_filename)
        print('\nCluster x population cross-tabulation:')
        for cluster in clusters:
            pops = {}
            for i, nm in enumerate(data['names']):
                if int(labels[i]) == cluster and nm in indivs_pops:
                    pops[indivs_pops[nm]] = pops.get(indivs_pops[nm], 0) + 1
            summary = ', '.join('{0}={1}'.format(p, c)
                                for p, c in sorted(pops.items()))
            print('C{0}: {1}'.format(cluster, summary or 'NA'))

    print('\n###5 - Differentiation summary (selected K = {0})'.format(
        selected_k))
    stats = cluster_allele_stats(data['dosage'], labels)
    priv = private_alleles(stats)
    fixed_diff, _ = pairwise_matrices(stats)
    fst = pairwise_fst(stats)
    het_obs, poly = cluster_diversity(data['state_code'],
                                      data['state_is_hom'], labels, stats)
    sizes = cluster_size_map(labels, clusters)
    print('Per-cluster summary:')
    print('{0:>6}  {1:>5}  {2:>8}  {3:>8}  {4:>6}  {5:>7}'.format(
        'clus', 'n', 'private', 'fixed', 'Ho', '%poly'))
    for ci, c in enumerate(clusters):
        print('{0:>6}  {1:>5}  {2:>8}  {3:>8}  {4:>6.3f}  {5:>7.2f}'.format(
            'C{0}'.format(c), sizes[c], priv[ci][0], priv[ci][1],
            het_obs[ci], poly[ci]))
    print('\nPairwise Fst (Hudson):')
    for a, ca in enumerate(clusters):
        cells = []
        for b in range(len(clusters)):
            if a == b:
                cells.append('{0:>7}'.format('-'))
            else:
                cells.append('{0:>7}'.format('{0:.3f}'.format(fst[a, b])
                             if fst[a, b] == fst[a, b] else 'NA'))
        print('C{0:>4} '.format(ca) + ' '.join(cells))
    print('Fixed differences:')
    for a, ca in enumerate(clusters):
        cells = ['{0:>7}'.format('-' if a == b else int(fixed_diff[a, b]))
                 for b in range(len(clusters))]
        print('C{0:>4} '.format(ca) + ' '.join(cells))

    print('\n###6 - Output files')
    csv_filename, pdf_filename = derive_outputs(output_filename, vcf_filename,
                                                pdf_output)
    write_assignment_csv(data['names'],
                         [r for r in records if r['K'] <= display_max],
                         csv_filename)
    if make_pdf:
        geno = (data['state_code'] >= 0).sum(axis=0)
        perc_genotyped = {nm: 100.0 * geno[i] / data['n_loci']
                          for i, nm in enumerate(data['names'])}
        print('PDF report: {0}'.format(pdf_filename))
        try:
            write_pdf_report(sim, dist, linkage_matrix, data['names'], records,
                             display_max, selected, perc_genotyped, data,
                             method, tree_mode, ordination, pdf_filename)
        except Exception as error:
            sys.stderr.write('Warning: PDF report failed ({0})\n'.format(error))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-p', '--pop', dest='pop_filename', metavar='pop_file',
                        help='optional text file (tsv/csv) with individuals and '
                             'populations (used only for a cluster x population '
                             'cross-tabulation, not for clustering)')
    parser.add_argument('-o', '--output', dest='output_filename',
                        metavar='cluster_file',
                        help='output file (csv) for per-sample cluster '
                             'assignments across K (default: derived from vcf)')
    parser.add_argument('-m', '--method', dest='method', default=DEFAULT_METHOD,
                        choices=vcf_clone_detect.METHODS,
                        help='similarity measure (default: ibs; dosage is the '
                             'documented alternative)')
    parser.add_argument('--auto-clone', dest='auto_clone', action='store_true',
                        help='clone-correct using an auto-inferred clone '
                             'threshold (default: clones are NOT removed)')
    parser.add_argument('--clone-list', dest='clone_list', metavar='file',
                        default=None,
                        help='external list of samples to drop for clone-'
                             'correction (e.g. vcf_clone_detect "to remove" '
                             'output)')
    parser.add_argument('--clone-threshold', dest='clone_threshold',
                        metavar='pct', default=None,
                        help='clone-correct using this manual similarity %% '
                             'threshold above which individuals are clones')
    parser.add_argument('--max-k', dest='max_k', type=int, default=DEF_MAX_K,
                        metavar='K', help='maximum K to evaluate (default: '
                        '{0})'.format(DEF_MAX_K))
    parser.add_argument('--min-cluster-size', dest='min_cluster_size', type=int,
                        default=DEF_MIN_CLUSTER_SIZE, metavar='N',
                        help='stop splitting once a cluster would fall below '
                             'this size (default: {0})'.format(
                                 DEF_MIN_CLUSTER_SIZE))
    parser.add_argument('--min-gap', dest='min_gap', type=float,
                        default=DEF_MIN_GAP, metavar='frac',
                        help='minimum relative separation gap for a split to '
                             'count as clean (default: {0})'.format(DEF_MIN_GAP))
    parser.add_argument('--k', dest='force_k', type=int, default=None,
                        metavar='K', help='force which K drives the '
                        'differentiation panels (default: best-supported K)')
    parser.add_argument('--tree', dest='tree_mode', default='upgma',
                        choices=('upgma', 'nj'),
                        help='tree to draw (default: upgma; nj is display only, '
                             'clusters still come from upgma)')
    parser.add_argument('--ordination', dest='ordination', default='pca',
                        choices=('pca', 'pcoa'),
                        help='ordination for the scatter panel (default: pca)')
    parser.add_argument('--pdf-output', dest='pdf_output', default=None,
                        metavar='pdf_file', help='filename for the PDF report')
    parser.add_argument('--no-pdf', dest='no_pdf', action='store_true',
                        help='do not generate the PDF report (text only)')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.output_filename,
         args.method, args.max_k, args.min_cluster_size, args.min_gap,
         args.force_k, args.tree_mode, args.ordination, args.clone_list,
         args.clone_threshold, args.auto_clone, make_pdf=not args.no_pdf,
         pdf_output=args.pdf_output)
