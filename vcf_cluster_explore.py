#!/usr/bin/env python
"""
Explores the deeper genetic clusters (lineages) in a `.vcf` dataset in a way
that is NOT biased by the presence of clones. Where `vcf_clone_detect.py` asks
"which samples are near-identical (clones)", this script asks "what are the
genetic groups, how many of them are there, and how differentiated are they".

Clones are non-independent samples that distort per-group allele frequencies and
therefore bias essentially every population-genetic measure (Fst, private and
fixed-private alleles, heterozygosity, ordination, and even the clustering
itself). By default the script does NOT remove clones - it runs on all
samples as given (use this when the input is already clone-corrected, or to
inspect the raw structure). Optional CLONE-CORRECTION reduces each clonal genet
to a single representative ramet (the one with the least missing data) and runs
the entire analysis - similarity, tree, K-evaluation and all differentiation
statistics - on that clone-corrected (genet) set. Enable it with `--auto-clone`
(clonal genets detected internally by reusing `vcf_clone_detect.py` at an auto-
inferred threshold), `--clone-threshold PCT` (manual threshold), or
`--clone-list FILE` (external list of samples to drop, e.g. the "samples to
remove" output of `vcf_clone_detect.py`).

The script (1) computes pairwise genetic similarities (`--method`, default
`ibs`; `dosage` is the documented alternative), (2) builds a UPGMA / average-
linkage tree, (3) evaluates K = 2 .. `--max-k`, reporting for each K the
genetic-similarity cut-off that splits the tree into K groups, the merge-height
separation gap, and the silhouette width, and picks the best K by silhouette
(robust to between-cluster GD overlap, unlike the raw gap), (4) assigns every
sample to a cluster at each K, and (5) produces a two-page PDF. PAGE 1 is a
full-page-width tree with readable tip labels, the per-K cluster-assignment
columns, a % genotyped bar and - when a popfile (or `--pops-from-sample-id`) is
given - one dot-strip panel per metadata track (population from column 2 and,
if present, a lineage/species/region from column 3; or the 2nd/3rd/4th `_`-
delimited fields of the sample name) so each sample's category memberships line
up with its tip; the selected-K column is boxed, and tracks can be renamed with
`--fields`. PAGE 2 is a fixed A4 page (for supplementary docs), a row-based
grid: (row 1) a fanned circular tree with cluster-coloured branches/tips, then
ordination axes 1-vs-2 and 2-vs-3 with hulls; (row 2) the metric-vs-K support
curve and per-sample silhouette; (row 3, with tracks) one stacked bar per field
showing how the clusters distribute across that field's categories; (row 4)
white-yellow-red differentiation heatmaps - shared loci (or shared genotyped
SNPs without a `.loci` file), pairwise private alleles excluding singletons, and
fixed differences (all counted over sites with >= 2 samples per cluster).

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
__copyright__ = 'Copyright (C) 2026 Pim Bongaerts'
__license__ = 'GPL'


DEFAULT_METHOD = 'ibs'
DEF_MAX_K = 10
DEF_MIN_CLUSTER_SIZE = 2
MAX_CATEGORIES = 20            # annotation tracks with more values are dropped
NAME_TRIM = 32                 # tree tip labels are trimmed to this many chars
TRACK_COLORS = ['0.0', '0.45', '0.7']   # one grey per field (black -> light)
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
        match_mat, both_mat, geno_count = vcf_clone_detect.compute_pair_matrices(
            data, method)
        comparisons = vcf_clone_detect.comparisons_for_indices(
            names, match_mat, both_mat, geno_count, {}, range(len(names)))
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
def silhouette_samples(dist, labels):
    """ Per-sample silhouette width from a precomputed dissimilarity matrix.
    s(i) compares how much better sample i fits its own cluster than the nearest
    other cluster; robust to between-cluster distance overlap (unlike the raw
    merge-height gap). Returns an array (one value per sample, range ~[-1, 1]). """
    clusters = np.unique(labels)
    n_clusters = len(clusters)
    n = len(labels)
    if n_clusters < 2:
        return np.zeros(n)
    mean_to = np.zeros((n, n_clusters))      # mean dist from sample to cluster c
    counts = np.zeros(n_clusters)
    own = np.zeros(n, dtype=int)
    for ci, cluster in enumerate(clusters):
        member = labels == cluster
        counts[ci] = member.sum()
        mean_to[:, ci] = dist[:, member].mean(axis=1)
        own[member] = ci
    own_count = counts[own]
    own_mean = mean_to[np.arange(n), own]
    # a(i): mean within-cluster distance, excluding self (self-distance is 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        a = np.where(own_count > 1, own_mean * own_count / (own_count - 1), 0.0)
    # b(i): smallest mean distance to any other cluster
    other = mean_to.copy()
    other[np.arange(n), own] = np.inf
    b = other.min(axis=1)
    denom = np.maximum(a, b)
    return np.where(denom > 0, (b - a) / denom, 0.0)


def mean_silhouette(dist, labels):
    """ Mean silhouette width (higher = better-separated clusters). """
    return float(np.mean(silhouette_samples(dist, labels)))


def calinski_harabasz(dist, labels):
    """ Calinski-Harabasz index from a precomputed distance matrix (via the
    Huygens identity on squared distances). A second, variance-based K-selection
    index that corroborates the silhouette; higher = better-separated. """
    clusters = np.unique(labels)
    k = len(clusters)
    n = len(labels)
    if k < 2 or n <= k:
        return 0.0
    sq = dist ** 2
    total_ss = sq.sum() / (2.0 * n)
    within_ss = 0.0
    for cluster in clusters:
        member = labels == cluster
        nc = int(member.sum())
        if nc > 0:
            within_ss += sq[np.ix_(member, member)].sum() / (2.0 * nc)
    between_ss = total_ss - within_ss
    if within_ss <= 0:
        return 0.0
    return float((between_ss / (k - 1)) / (within_ss / (n - k)))


def evaluate_k(linkage_matrix, dist, n, max_k, min_cluster_size):
    """ Evaluate every K from 2 to max_k: the similarity cut-off that yields K
    clusters, the merge-height separation gap, and the silhouette width. Returns
    (records, best_k). The whole range is scanned (no early break); best_k is the
    K with the highest silhouette among splits whose clusters all meet
    `min_cluster_size` (so that adding singletons at high K does not win).

    Indexing: `Z[:,2]` are the n-1 merge heights ascending. K clusters occur for
    a cut height in (Z[n-1-K, 2], Z[n-K, 2]); the gap supporting K is the width
    of that interval. """
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
                        'silhouette': mean_silhouette(dist, labels),
                        'ch': calinski_harabasz(dist, labels),
                        'min_size': int(sizes.min()),
                        'n_clusters': int((sizes > 0).sum())})

    valid = [r for r in records if r['min_size'] >= min_cluster_size]
    pool = valid if valid else records
    best_k = max(pool, key=lambda r: r['silhouette'])['K']
    return records, best_k


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


def pairwise_private_no_singletons(stats):
    """ Pairwise private-allele counts EXCLUDING singletons: priv[a, b] = alleles
    present in cluster a on >= 2 allele copies (so not a singleton) yet entirely
    absent from cluster b (which has data at the site). KxK, NaN-free ints. """
    n_clusters = stats['n_clusters']
    altcnt, ncall = stats['altcnt'], stats['ncall']
    alt_present, ref_present = stats['alt_present'], stats['ref_present']
    has_data = stats['has_data']
    refcnt = 2.0 * ncall - altcnt
    alt_ok = altcnt >= 2                      # >= 2 copies -> not a singleton
    ref_ok = refcnt >= 2
    priv = np.zeros((n_clusters, n_clusters), dtype=int)
    for a in range(n_clusters):
        for b in range(n_clusters):
            if a == b:
                continue
            priv[a, b] = int(
                (alt_ok[a] & ~alt_present[b] & has_data[b]).sum() +
                (ref_ok[a] & ~ref_present[b] & has_data[b]).sum())
    return priv


def pairwise_fixed_diff_min2(stats):
    """ Pairwise fixed-difference counts requiring >= 2 genotyped samples in
    BOTH clusters at a site (so a difference is never called off a single
    sample). Symmetric KxK int matrix. """
    n_clusters = stats['n_clusters']
    freq, ncall = stats['freq'], stats['ncall']
    fixed = np.zeros((n_clusters, n_clusters), dtype=int)
    for a in range(n_clusters):
        for b in range(a + 1, n_clusters):
            valid = ((ncall[a] >= 2) & (ncall[b] >= 2)
                     & np.isfinite(freq[a]) & np.isfinite(freq[b]))
            fdiff = int((valid & (np.abs(freq[a] - freq[b]) == 1.0)).sum())
            fixed[a, b] = fixed[b, a] = fdiff
    return fixed


def pairwise_shared_snps(stats):
    """ Number of SNP sites genotyped in >= 2 samples of both clusters, for every
    cluster pair (the SNP-based fallback for shared loci when no `.loci` file is
    given; the >= 2 rule matches the other differentiation panels). Symmetric
    KxK int matrix (NaN-free). """
    n_clusters = stats['n_clusters']
    ncall = stats['ncall']
    shared = np.zeros((n_clusters, n_clusters), dtype=int)
    for a in range(n_clusters):
        for b in range(n_clusters):
            if a != b:
                shared[a, b] = int(((ncall[a] >= 2) & (ncall[b] >= 2)).sum())
    return shared


def shared_loci_matrix(presence, labels, names):
    """ Number of loci recovered in >= 2 samples of both clusters, for every
    cluster pair (from an ipyrad `.loci` presence list; the >= 2 rule matches
    the other differentiation panels). Symmetric KxK int. """
    name_cluster = {nm: int(labels[i]) for i, nm in enumerate(names)}
    clusters = sorted(set(name_cluster.values()))
    cidx = {c: i for i, c in enumerate(clusters)}
    k = len(clusters)
    counts = np.zeros((k, len(presence)), dtype=int)     # samples/cluster/locus
    for li, locus in enumerate(presence):
        for nm in locus:
            counts[cidx[name_cluster[nm]], li] += 1
    present2 = counts >= 2
    shared = np.zeros((k, k), dtype=int)
    for a in range(k):
        for b in range(k):
            if a != b:
                shared[a, b] = int((present2[a] & present2[b]).sum())
    return shared


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


def pairwise_dxy(stats):
    """ Absolute divergence dxy between every cluster pair: mean over sites (with
    data in both) of the expected proportion of pairwise allelic differences
    between a sequence drawn from each cluster, `pa*(1-pb) + pb*(1-pa)`. Unlike
    Fst (relative), dxy is not deflated by low within-cluster diversity. Returns
    KxK (NaN diag). """
    n_clusters = stats['n_clusters']
    freq, has_data = stats['freq'], stats['has_data']
    dxy = np.full((n_clusters, n_clusters), np.nan)
    for a in range(n_clusters):
        for b in range(a + 1, n_clusters):
            pa, pb = freq[a], freq[b]
            use = has_data[a] & has_data[b] & np.isfinite(pa) & np.isfinite(pb)
            if use.any():
                value = float((pa[use] * (1 - pb[use])
                               + pb[use] * (1 - pa[use])).mean())
            else:
                value = np.nan
            dxy[a, b] = dxy[b, a] = value
    return dxy


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


def parse_loci_presence(filename, names):
    """ Parse an ipyrad `.loci` file into per-locus sample-presence sets.

    Format: each locus is a block of `samplename<whitespace>sequence` lines,
    terminated by a `//...|locus_id|` separator line. A sample line means that
    sample recovered the locus. Returns a list (one set of sample names per
    locus, restricted to `names`) and the count of `.loci` samples that matched
    the VCF. """
    name_set = set(names)
    presence = []
    current = set()
    matched = set()
    with open(filename, 'r') as handle:
        for line in handle:
            if line.startswith('//'):
                if current:
                    presence.append(current)
                current = set()
                continue
            sample = line.split(None, 1)[0] if line.strip() else ''
            if sample in name_set:
                current.add(sample)
                matched.add(sample)
    if current:
        presence.append(current)
    return presence, len(matched)


def locus_sharing_stats(presence, labels, names):
    """ From per-locus sample presence + cluster labels, compute per-cluster
    private loci (recovered in >=1 sample of only that cluster) and a pairwise
    shared-loci Jaccard matrix. Returns (clusters, private_counts, jaccard). """
    name_cluster = {nm: int(labels[i]) for i, nm in enumerate(names)}
    clusters = sorted(set(name_cluster.values()))
    cidx = {c: i for i, c in enumerate(clusters)}
    k = len(clusters)
    present_any = np.zeros((k, len(presence)), dtype=bool)   # cluster x locus
    for li, locus in enumerate(presence):
        for nm in locus:
            present_any[cidx[name_cluster[nm]], li] = True
    n_present = present_any.sum(axis=0)
    private = np.zeros(k, dtype=int)
    for ci in range(k):
        private[ci] = int((present_any[ci] & (n_present == 1)).sum())
    jaccard = np.full((k, k), np.nan)
    for a in range(k):
        for b in range(k):
            union = (present_any[a] | present_any[b]).sum()
            inter = (present_any[a] & present_any[b]).sum()
            jaccard[a, b] = inter / union if union > 0 else np.nan
    return clusters, private, jaccard


def ordinate(dosage, dist, kind, n_comp=3):
    """ Ordination of samples on the first `n_comp` axes. Returns
    (coords[n, n_comp], pct_var[n_comp]). """
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
        coords = u_mat[:, :n_comp] * sing[:n_comp]
        total = (sing ** 2).sum()
        pct = (100.0 * sing[:n_comp] ** 2 / total if total > 0
               else np.zeros(n_comp))
        return coords, pct
    # PCoA on the distance matrix (coherent with the tree, but D is non-Euclidean)
    n = dist.shape[0]
    centering = np.eye(n) - 1.0 / n
    gram = -0.5 * centering @ (dist ** 2) @ centering
    eigval, eigvec = np.linalg.eigh(gram)
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    coords = eigvec[:, :n_comp] * np.sqrt(np.clip(eigval[:n_comp], 0, None))
    positive = np.clip(eigval, 0, None).sum()
    pct = (100.0 * np.clip(eigval[:n_comp], 0, None) / positive
           if positive > 0 else np.zeros(n_comp))
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
    """ A list of n categorical colours. """
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


def assign_hierarchical_colors(records):
    """ Descent-consistent cluster colours across K. Because UPGMA cuts are
    strictly nested, a cluster either persists or splits as K grows; we let a
    colour follow a lineage so a new colour only appears at a genuine split (and
    is never reused for an unrelated cluster). Returns {K: {cluster_id: rgba}}.

    Per K (ascending), each cluster is a frozenset of member indices. Processing
    new clusters largest-first: a cluster identical to a previous one keeps its
    colour; otherwise its parent is the previous cluster that is its superset and
    the first (largest) child to claim that parent's colour inherits it, while
    later children draw the next unused palette colour. """
    max_clusters = max(r['n_clusters'] for r in records)
    palette = cluster_palette(max_clusters)
    ordered = sorted(records, key=lambda r: r['K'])
    colors_by_k = {}
    prev = {}                       # frozenset(members) -> rgba (previous K)
    next_color = 0
    for record in ordered:
        labels = record['labels']
        members = {}
        for i, lab in enumerate(labels):
            members.setdefault(int(lab), set()).add(i)
        parts = {cid: frozenset(s) for cid, s in members.items()}
        current = {}                # frozenset -> rgba (this K)
        claimed = set()             # previous framesets whose colour was taken
        cid_color = {}
        # Largest clusters first so the dominant child inherits a parent colour
        for cid in sorted(parts, key=lambda c: (-len(parts[c]), min(parts[c]))):
            part = parts[cid]
            if part in prev:                         # cluster persisted
                colour = prev[part]
            else:
                parent = next((p for p in prev if part <= p
                               and p not in claimed), None)
                if parent is not None:               # first child inherits
                    colour = prev[parent]
                    claimed.add(parent)
                else:                                # new lineage -> new colour
                    colour = palette[next_color % len(palette)]
                    next_color += 1
            current[part] = colour
            cid_color[cid] = colour
        colors_by_k[record['K']] = cid_color
        prev = current
    return colors_by_k


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


def draw_circular_tree(ax, tree, sample_cluster, cluster_color, title,
                       span_deg=300.0):
    """ Draw a Bio.Phylo tree as a FANNED radial dendrogram: radius = branch-
    length distance from the root (centre), tips spread over a `span_deg` fan
    (an open wedge, not a full circle, so the deep splits open up), with an extra
    angular gap inserted between clusters so genetically differentiated groups
    separate visually. Branches that sit entirely within one cluster are drawn in
    that cluster's colour (grey where a branch still joins several clusters); tip
    dots are coloured by cluster to match the PCA panels. """
    tips = tree.get_terminals()
    n = len(tips)
    if n < 2:
        ax.axis('off')
        return
    # slot positions along the fan, with a gap between adjacent clusters
    gap = max(1.0, n / 40.0)
    slots, pos, prev = {}, 0.0, None
    for tip in tips:
        cid = sample_cluster.get(tip.name)
        if prev is not None and cid != prev:
            pos += gap
        slots[tip] = pos
        pos += 1.0
        prev = cid
    total = max(pos - 1.0, 1.0)
    span = np.radians(span_deg)
    a0 = np.pi / 2.0 + span / 2.0             # fan centred on top, opens downward
    tip_angle = {tip: a0 - span * (slots[tip] / total) for tip in tips}
    depths = tree.depths()                     # branch-length distance from root
    # post-order (reversed level order): mean tip angle + cluster set per clade
    order = list(tree.find_clades(order='level'))
    ang_sum, ang_cnt, clset = {}, {}, {}
    for clade in reversed(order):
        if clade.is_terminal():
            ang_sum[clade], ang_cnt[clade] = tip_angle[clade], 1
            clset[clade] = {sample_cluster.get(clade.name)}
        else:
            ang_sum[clade] = sum(ang_sum[c] for c in clade.clades)
            ang_cnt[clade] = sum(ang_cnt[c] for c in clade.clades)
            clset[clade] = set().union(*(clset[c] for c in clade.clades))
    angle = {clade: ang_sum[clade] / ang_cnt[clade] for clade in ang_sum}

    def edge_color(clade):
        cids = clset[clade]
        if len(cids) == 1:
            return cluster_color.get(next(iter(cids)), '0.35')
        return '0.4'

    for clade in order:
        if clade.is_terminal():
            continue
        r = depths[clade]
        child_angles = [angle[c] for c in clade.clades]
        thetas = np.linspace(min(child_angles), max(child_angles), 40)
        ax.plot(r * np.cos(thetas), r * np.sin(thetas),
                color=edge_color(clade), lw=0.6, zorder=1)
        for child, a in zip(clade.clades, child_angles):
            rc = depths[child]
            ax.plot([r * np.cos(a), rc * np.cos(a)],
                    [r * np.sin(a), rc * np.sin(a)], color=edge_color(child),
                    lw=0.6, zorder=2)
    xs, ys, cols = [], [], []
    for tip in tips:
        a, r = tip_angle[tip], depths[tip]
        xs.append(r * np.cos(a))
        ys.append(r * np.sin(a))
        cols.append(cluster_color.get(sample_cluster.get(tip.name), '0.5'))
    ax.scatter(xs, ys, s=14, color=cols, edgecolor='none', zorder=3)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=8)


def make_ax_in(fig, fig_w, fig_h):
    """ Return an `ax_in(x, ytop, w, h)` that places an axes by inches-from-the-
    top-left on `fig` (whose size is fig_w x fig_h inches). """
    def ax_in(x_in, ytop_in, w_in, h_in):
        return fig.add_axes([x_in / fig_w,
                             (fig_h - ytop_in - h_in) / fig_h,
                             w_in / fig_w, h_in / fig_h])
    return ax_in


def despine(ax):
    """ Hide the top and right spines of an axes. """
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)


def draw_heatmap(fig, ax, mat, cmap, fmt, title, n_clusters, cluster_labels):
    """ Draw a KxK cluster-pair heatmap (NaN diagonal) with cell annotations.
    Annotation text colour follows each cell's background luminance, so it stays
    legible on any colormap (including the light end of white->yellow->red). """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    disp = np.array(mat, dtype=float)
    np.fill_diagonal(disp, np.nan)
    image = ax.imshow(disp, cmap=cmap, aspect='auto')
    ax.set_xticks(range(n_clusters))
    ax.set_yticks(range(n_clusters))
    ax.set_xticklabels(cluster_labels, fontsize=6, rotation=90)
    ax.set_yticklabels(cluster_labels, fontsize=6)
    finite = disp[np.isfinite(disp)]
    lo = float(finite.min()) if finite.size else 0.0
    hi = float(finite.max()) if finite.size else 1.0
    cmap_obj = plt.get_cmap(cmap)
    norm = Normalize(vmin=lo, vmax=hi)
    for a in range(n_clusters):
        for b in range(n_clusters):
            if a != b and np.isfinite(disp[a, b]):
                r, g, bl, _ = cmap_obj(norm(disp[a, b]))
                lum = 0.299 * r + 0.587 * g + 0.114 * bl
                ax.text(b, a, fmt.format(mat[a][b]), ha='center', va='center',
                        fontsize=5, color='black' if lum > 0.5 else 'white')
    ax.set_title(title, fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)


# --------------------------------------------------------------------------- #
#  Sample metadata annotation tracks (population / lineage, drawn by the tree)
# --------------------------------------------------------------------------- #
def _make_track(name, assign, names):
    """ Package a sample->value mapping into an annotation track (with the
    ordered set of distinct values), restricted to `names`. """
    values = sorted(set(assign[nm] for nm in names if nm in assign))
    return {'name': name, 'assign': assign, 'values': values}


def build_category_tracks(names, pop_filename, pops_from_sample_id,
                          field_names=None):
    """ Build the metadata annotation tracks drawn beside the tree, each a
    sample->category mapping visualised as one dot per sample on an x-axis of
    that category's distinct values.

    With `--pops-from-sample-id` three tracks come from the 2nd, 3rd and 4th
    `_`-delimited fields of each sample name (the 1st field is skipped). From a
    popfile: a `population` track (column 2) and, when a 3rd column is present, a
    `lineage / region` track (column 3). `field_names` (e.g. from `--fields
    location,depth`) overrides the default track titles positionally (name 1 =
    field 2 / population, name 2 = field 3 / lineage, ...). Returns a (possibly
    empty) list of tracks; tracks with no usable values, or with more than
    MAX_CATEGORIES distinct values (e.g. a per-sample id field), are dropped
    (with a note). """
    field_names = field_names or []

    def title_for(pos, default):
        return field_names[pos] if pos < len(field_names) else default

    tracks = []
    if pops_from_sample_id:
        for pos, (default, idx) in enumerate((('id field 2', 1),
                                              ('id field 3', 2),
                                              ('id field 4', 3))):
            assign = {}
            for nm in names:
                parts = nm.split('_')
                if len(parts) > idx:
                    assign[nm] = parts[idx]
            if assign:
                tracks.append(_make_track(title_for(pos, default), assign,
                                          names))
    elif pop_filename:
        indivs_pops, indivs_lineages = \
            vcf_clone_detect.get_assignments_from_popfile(pop_filename)
        if any(nm in indivs_pops for nm in names):
            tracks.append(_make_track(title_for(0, 'population'), indivs_pops,
                                      names))
        if any(nm in indivs_lineages for nm in names):
            tracks.append(_make_track(title_for(1, 'lineage / region'),
                                      indivs_lineages, names))
    kept = []
    for track in tracks:
        n_values = len(track['values'])
        if n_values == 0:
            continue
        if n_values > MAX_CATEGORIES:
            sys.stderr.write('Note: annotation track `{0}` has {1} distinct '
                             'values (> MAX_CATEGORIES={2}); not drawn.\n'.format(
                                 track['name'], n_values, MAX_CATEGORIES))
            continue
        kept.append(track)
    return kept


def draw_category_track(ax, track, tip_y, ylim, dot_size, color):
    """ Draw one annotation track aligned to the tree tips: for every sample a
    dot at the x-position of its category value (categories along the x-axis,
    samples along the shared tree y-axis). All dots of a track share one colour
    (the value is read off the x-axis, not the colour). """
    values = track['values']
    vidx = {v: i for i, v in enumerate(values)}
    for i in range(len(values)):
        ax.axvline(i, color='0.92', lw=0.5, zorder=0)
    xs, ys = [], []
    for nm, y in tip_y.items():
        value = track['assign'].get(nm)
        if value is None:
            continue
        xs.append(vidx[value])
        ys.append(y)
    ax.scatter(xs, ys, s=dot_size, color=color, edgecolor='none', zorder=3)
    ax.set_xlim(-0.5, len(values) - 0.5)
    ax.set_ylim(ylim)
    ax.set_yticks([])
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(values, fontsize=6, rotation=90)
    ax.tick_params(length=0)
    ax.set_title(track['name'], fontsize=7)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)


# --------------------------------------------------------------------------- #
#  PDF report (page 1: full-width tree + annotations; page 2: analysis panels)
# --------------------------------------------------------------------------- #
def write_tree_page(pdf, dist, linkage_matrix, names, display, selected,
                    perc_genotyped, method, tree_mode, colors_by_k, sel_color,
                    category_tracks):
    """ PAGE 1 - a full-page-width tree with readable tip labels, followed (left
    to right) by the per-K cluster-assignment columns, a %genotyped bar and one
    dot-strip panel per metadata track (population / lineage or the sample-id
    fields). Returns the tip y-positions dict, or None if tips can't be found. """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    import Bio.Phylo

    name_idx = {nm: i for i, nm in enumerate(names)}
    n_tips = len(names)
    n_cols = len(display)

    left_margin, right_margin = 0.6, 0.6
    top_margin, bottom_margin = 0.6, 0.6
    tree_w = 7.5
    kcol_w = max(0.9, 0.22 * n_cols)
    bar_w = 1.0
    track_gap = 0.15
    track_ws = [max(0.7, 0.16 * len(t['values'])) for t in category_tracks]
    per_tip = 0.16
    dot_size = 16

    tree_h = max(3.0, n_tips * per_tip)
    fig_w = (left_margin + tree_w + 0.1 + kcol_w + 0.15 + bar_w
             + sum(w + track_gap for w in track_ws) + right_margin)
    fig_h = top_margin + tree_h + bottom_margin
    font_size = max(5.0, min(9.0, per_tip * 72 * 0.7))

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax_in = make_ax_in(fig, fig_w, fig_h)

    # --- tree ----------------------------------------------------------------
    ax_tree = ax_in(left_margin, top_margin, tree_w, tree_h)
    if tree_mode == 'nj':
        tree = build_nj_tree(dist, names)
    else:
        tree = build_upgma_tree(linkage_matrix, names)
    sys.stderr.write('Drawing {0} tree ({1} tips)...\n'.format(tree_mode, n_tips))
    Bio.Phylo.draw(tree, axes=ax_tree, do_show=False,
                   label_func=lambda c: '  ' + c.name if c.name else '')
    ax_tree.set_ylabel('')
    ax_tree.set_xlabel('Genetic distance (1 - {0} similarity)'.format(method),
                       fontsize=8)
    ax_tree.tick_params(labelsize=7)

    tip_y = {}
    max_disp = 1
    for text in ax_tree.texts:
        label = text.get_text().strip()
        if label in name_idx:
            tip_y[label] = text.get_position()[1]
            text.set_fontsize(font_size)
            text.set_color(sel_color[label])
            # trim long names so the labels stay inside the tree box
            shown = label if len(label) <= NAME_TRIM else label[:NAME_TRIM] + '…'
            text.set_text('  ' + shown)
            max_disp = max(max_disp, len(shown))
    if not tip_y:
        sys.stderr.write('Warning: could not recover tip positions.\n')
        plt.close(fig)
        return None
    ys = sorted(tip_y.values())
    row_pitch = (ys[1] - ys[0]) if len(ys) > 1 else 1.0
    # Reserve horizontal room for the (trimmed) tip labels so they don't spill
    # out of the tree box into the columns: extend xlim by the label width.
    xmin, xmax = ax_tree.get_xlim()
    data_per_inch = (xmax - xmin) / tree_w if tree_w > 0 else 1.0
    label_inch = (max_disp + 2) * font_size * 0.6 / 72.0
    ax_tree.set_xlim(xmin, xmax + label_inch * data_per_inch)
    ylim = ax_tree.get_ylim()
    if tree_mode == 'nj':
        title = '{0} tree (display); clusters = UPGMA'.format(tree_mode.upper())
    else:
        title = 'UPGMA tree (clusters = tree cuts)'
    ax_tree.set_title(title, fontsize=9)

    # --- per-K cluster-assignment columns (descent-consistent colours) -------
    ax_k = ax_in(left_margin + tree_w + 0.1, top_margin, kcol_w, tree_h)
    rects, colors = [], []
    for col, record in enumerate(display):
        cid_color = colors_by_k[record['K']]
        for nm, y in tip_y.items():
            cid = int(record['labels'][name_idx[nm]])
            rects.append(Rectangle((col, y - row_pitch * 0.45), 1,
                                   row_pitch * 0.9))
            colors.append(cid_color[cid])
    ax_k.add_collection(PatchCollection(rects, facecolors=colors,
                                        edgecolors='none'))
    ax_k.set_xlim(0, max(n_cols, 1))
    ax_k.set_ylim(ylim)
    ax_k.set_yticks([])
    ax_k.set_xticks([c + 0.5 for c in range(n_cols)])
    ax_k.set_xticklabels(['K={0}'.format(r['K']) for r in display],
                         fontsize=6, rotation=90)
    ax_k.tick_params(length=0)
    ax_k.set_title('cluster @ K', fontsize=7)
    # Highlight the selected-K column with a thick black border + bold label
    sel_col = next((i for i, r in enumerate(display)
                    if r['K'] == selected['K']), None)
    if sel_col is not None:
        y0, y1 = min(ylim), max(ylim)
        ax_k.add_patch(Rectangle((sel_col, y0), 1, y1 - y0, fill=False,
                                 edgecolor='black', lw=2.5, zorder=5))
        ax_k.get_xticklabels()[sel_col].set_fontweight('bold')

    # --- % genotyped bar -----------------------------------------------------
    ax_bar = ax_in(left_margin + tree_w + 0.1 + kcol_w + 0.15, top_margin,
                   bar_w, tree_h)
    bar_names = [nm for nm in names if nm in tip_y]
    bar_y = [tip_y[nm] for nm in bar_names]
    bar_v = [perc_genotyped[nm] for nm in bar_names]
    ax_bar.barh(bar_y, bar_v, height=row_pitch * 0.8,
                color=[sel_color[nm] for nm in bar_names], edgecolor='none')
    ax_bar.set_ylim(ylim)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(max(0, min(bar_v) - 5), 100)
    ax_bar.set_title('% genotyped', fontsize=7)
    ax_bar.tick_params(labelsize=6)
    despine(ax_bar)

    # --- metadata annotation tracks (population / lineage / sample-id) -------
    # one grey per field (cycled) so a track reads as a single category axis
    x_track = left_margin + tree_w + 0.1 + kcol_w + 0.15 + bar_w
    for i, (track, tw) in enumerate(zip(category_tracks, track_ws)):
        x_track += track_gap
        ax_t = ax_in(x_track, top_margin, tw, tree_h)
        color = TRACK_COLORS[i % len(TRACK_COLORS)]
        draw_category_track(ax_t, track, tip_y, ylim, dot_size, color)
        x_track += tw

    pdf.savefig(fig)
    plt.close(fig)
    return tip_y


def write_analysis_page(pdf, dist, linkage_matrix, names, display, selected,
                        data, method, ordination, tree_mode, colors_by_k,
                        category_tracks, loci_presence):
    """ PAGE 2 (fixed A4) - analysis panels in a row-based grid:
      ROW 1: fanned circular tree | ordination axes 1v2 | ordination axes 2v3
             (all coloured by cluster, so differentiated groups pop out).
      ROW 2: metric-vs-K support curve | per-sample silhouette.
      ROW 3 (only with annotation tracks): one stacked bar per field showing how
             the clusters distribute across that field's categories (one column
             per field).
      ROW 4 (differentiation): shared loci (or shared genotyped SNPs when no
             `.loci`) | pairwise private alleles excluding singletons | fixed
             differences (>= 2 samples per cluster). """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from scipy.spatial import ConvexHull

    n_tips = len(names)
    name_idx = {nm: i for i, nm in enumerate(names)}
    labels = selected['labels']
    stats = cluster_allele_stats(data['dosage'], labels)
    clusters = stats['clusters']
    n_clusters = stats['n_clusters']
    sel_cid_color = colors_by_k[selected['K']]
    sample_cluster = {nm: int(labels[name_idx[nm]]) for nm in names}
    cluster_labels = ['C{0}'.format(c) for c in clusters]

    ks = [r['K'] for r in display]
    sils = [r['silhouette'] for r in display]
    chs = [r['ch'] for r in display]
    sample_sil = silhouette_samples(dist, labels)
    coords, pct = ordinate(data['dosage'], dist, ordination, n_comp=3)
    if tree_mode == 'nj':
        tree = build_nj_tree(dist, names)
    else:
        tree = build_upgma_tree(linkage_matrix, names)

    def make_ordination(ix, iy):    # scatter of axis ix vs iy + convex hulls
        def draw(ax):
            for c in clusters:
                pts = coords[labels == c]
                col = sel_cid_color[c]
                ax.scatter(pts[:, ix], pts[:, iy], s=14, color=col,
                           edgecolor='none', label='C{0}'.format(c))
                if pts.shape[0] >= 3:
                    try:
                        hull = ConvexHull(pts[:, [ix, iy]])
                        ax.add_patch(Polygon(pts[:, [ix, iy]][hull.vertices],
                                             closed=True, facecolor=col,
                                             alpha=0.15, edgecolor=col, lw=0.8))
                    except Exception:
                        pass
            ax.set_xlabel('{0}{1} ({2:.1f}%)'.format(ordination.upper(), ix + 1,
                                                     pct[ix]), fontsize=7)
            ax.set_ylabel('{0}{1} ({2:.1f}%)'.format(ordination.upper(), iy + 1,
                                                     pct[iy]), fontsize=7)
            ax.set_title('Ordination ({0}) {1} vs {2}'.format(
                ordination.upper(), ix + 1, iy + 1), fontsize=8)
            ax.legend(fontsize=5, ncol=2, framealpha=0.9)
            ax.tick_params(labelsize=6)
            despine(ax)
        return draw

    def draw_circular(ax):
        tlabel = 'NJ' if tree_mode == 'nj' else 'UPGMA'
        draw_circular_tree(ax, tree, sample_cluster, sel_cid_color,
                           'Circular {0} tree (fan)'.format(tlabel))

    def draw_a1(ax):        # metric-vs-K curve (silhouette + Calinski-Harabasz)
        ax.plot(ks, sils, '-o', color='#1f6f6f', markersize=4,
                label='silhouette')
        peak = max(display, key=lambda r: r['silhouette'])['K']
        ax.plot([peak], [dict(zip(ks, sils))[peak]], '*', color='#1f6f6f',
                markersize=12, zorder=5)
        ax.axvline(selected['K'], color='red', ls='--', lw=1,
                   label='selected K={0}'.format(selected['K']))
        ax.set_xlabel('K (number of clusters)', fontsize=7)
        ax.set_ylabel('silhouette width', fontsize=7, color='#1f6f6f')
        ax.set_xticks(ks)
        ax.tick_params(labelsize=6)
        axb = ax.twinx()
        axb.plot(ks, chs, '-s', color='0.6', markersize=3,
                 label='Calinski-Harabasz')
        axb.set_ylabel('Calinski-Harabasz', fontsize=7, color='0.5')
        axb.tick_params(labelsize=6)
        ax.set_title('Support for K (peak = best silhouette)', fontsize=8)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = axb.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=5, loc='best', framealpha=0.9)

    def draw_a2(ax):        # per-sample silhouette at the selected K
        cluster_gap = max(1, n_tips // 100)
        ypos = 0
        yticks, yticklabels = [], []
        for c in clusters:
            idx = np.where(labels == c)[0]
            vals = np.sort(sample_sil[idx])[::-1]
            yr = np.arange(ypos, ypos + len(vals))
            ax.barh(yr, vals, height=1.0, color=sel_cid_color[c],
                    edgecolor='none')
            yticks.append(ypos + len(vals) / 2.0)
            yticklabels.append('C{0}'.format(c))
            ypos += len(vals) + cluster_gap
        mean_s = float(np.mean(sample_sil))
        ax.axvline(mean_s, color='red', ls='--', lw=1,
                   label='mean {0:.3f}'.format(mean_s))
        ax.axvline(0, color='0.5', lw=0.6)
        ax.set_ylim(-1, ypos)
        ax.invert_yaxis()
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=6)
        ax.set_xlabel('silhouette width', fontsize=7)
        ax.set_title('Per-sample silhouette (K={0})'.format(selected['K']),
                     fontsize=8)
        ax.legend(fontsize=5, loc='lower right', framealpha=0.9)
        ax.tick_params(labelsize=6)
        despine(ax)

    def make_field_dist(track):     # stacked cluster counts across a field
        values = track['values']
        vidx = {v: i for i, v in enumerate(values)}
        cidx = {c: i for i, c in enumerate(clusters)}
        counts = np.zeros((n_clusters, len(values)))
        for i, nm in enumerate(names):
            value = track['assign'].get(nm)
            if value is not None:
                counts[cidx[int(labels[i])], vidx[value]] += 1

        def draw(ax):
            x = np.arange(len(values))
            bottom = np.zeros(len(values))
            for ci, c in enumerate(clusters):
                ax.bar(x, counts[ci], bottom=bottom, color=sel_cid_color[c],
                       edgecolor='white', linewidth=0.2, label='C{0}'.format(c))
                bottom += counts[ci]
            ax.set_xticks(x)
            ax.set_xticklabels(values, fontsize=6, rotation=90)
            ax.set_ylabel('# samples', fontsize=7)
            ax.set_title('Clusters by {0}'.format(track['name']), fontsize=8)
            ax.tick_params(labelsize=6)
            ax.legend(fontsize=4, ncol=2, framealpha=0.9)
            despine(ax)
        return draw

    # --- differentiation matrices (row 4) ------------------------------------
    if loci_presence is not None:
        shared = shared_loci_matrix(loci_presence, labels, names)
        shared_title = 'Shared loci'
    else:
        shared = pairwise_shared_snps(stats)
        shared_title = 'Shared genotyped SNPs'
    priv_pair = pairwise_private_no_singletons(stats)
    fixed_diff = pairwise_fixed_diff_min2(stats)

    def draw_shared(ax):
        draw_heatmap(fig, ax, shared, 'YlOrRd', '{0:d}', shared_title,
                     n_clusters, cluster_labels)

    def draw_priv_pair(ax):
        draw_heatmap(fig, ax, priv_pair, 'YlOrRd', '{0:d}',
                     'Private alleles (pairwise, no singletons)',
                     n_clusters, cluster_labels)

    def draw_fixed(ax):
        draw_heatmap(fig, ax, fixed_diff, 'YlOrRd', '{0:d}',
                     'Fixed differences (>=2 samples/cluster)',
                     n_clusters, cluster_labels)

    # --- assemble rows (each row is a list of panel callables) ---------------
    rows = [[draw_circular, make_ordination(0, 1), make_ordination(1, 2)],
            [draw_a1, draw_a2]]
    if category_tracks:
        rows.append([make_field_dist(t) for t in category_tracks])
    rows.append([draw_shared, draw_priv_pair, draw_fixed])

    # Fixed A4 portrait canvas so page 2 drops straight into supplementary docs
    fig_w, fig_h = 8.27, 11.69
    margin, col_gap, row_gap = 0.5, 0.6, 0.7
    n_rows = len(rows)
    content_w = fig_w - 2 * margin
    content_h = fig_h - 2 * margin
    row_h = (content_h - (n_rows - 1) * row_gap) / n_rows
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax_in = make_ax_in(fig, fig_w, fig_h)
    for ri, row in enumerate(rows):
        ncols = len(row)
        cw = (content_w - (ncols - 1) * col_gap) / ncols
        ytop = margin + ri * (row_h + row_gap)
        for ci, draw in enumerate(row):
            draw(ax_in(margin + ci * (cw + col_gap), ytop, cw, row_h))
    pdf.savefig(fig)
    plt.close(fig)


def write_pdf_report(sim, dist, linkage_matrix, names, records, col_max,
                     selected, perc_genotyped, data, method, tree_mode,
                     ordination, pdf_filename, category_tracks=None,
                     loci_presence=None):
    """ Two-page PDF. PAGE 1: a full-page-width tree with readable tip labels,
    the per-K cluster-assignment columns (selected K boxed), a %genotyped bar and
    one dot-strip panel per metadata track (population / lineage from a popfile,
    or the sample-id fields). PAGE 2 (row-based grid): ordination + circular tree;
    metric-vs-K curve + per-sample silhouette; per-field cluster-distribution
    stacked bars; and differentiation heatmaps (shared loci or genotyped SNPs,
    pairwise private alleles excl. singletons, fixed differences). Cluster colours
    are descent-consistent across K (assign_hierarchical_colors), shared by both
    pages. """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    category_tracks = category_tracks or []
    name_idx = {nm: i for i, nm in enumerate(names)}
    display = [r for r in records if r['K'] <= col_max]

    # Descent-consistent colours, shared by K-columns, tips, bar and all panels
    colors_by_k = assign_hierarchical_colors(display)
    labels = selected['labels']
    sel_cid_color = colors_by_k[selected['K']]
    sel_color = {nm: sel_cid_color[int(labels[name_idx[nm]])] for nm in names}

    with PdfPages(pdf_filename) as pdf:
        tip_y = write_tree_page(pdf, dist, linkage_matrix, names, display,
                                selected, perc_genotyped, method, tree_mode,
                                colors_by_k, sel_color, category_tracks)
        if tip_y is None:
            plt.close('all')
            return
        write_analysis_page(pdf, dist, linkage_matrix, names, display, selected,
                            data, method, ordination, tree_mode, colors_by_k,
                            category_tracks, loci_presence)
    sys.stderr.write('PDF report written to {0}\n'.format(pdf_filename))


# --------------------------------------------------------------------------- #
#  Text-output helpers
# --------------------------------------------------------------------------- #
def print_k_table(records, best_k, selected_k):
    """ ###3 - K evaluation table. """
    print('{0:>3}  {1:>10}  {2:>8}  {3:>8}  {4:>10}  {5:>9}  {6:>10}  {7}'.format(
        'K', 'cutoff(%)', 'abs_gap', 'rel_gap', 'silhouette', 'min_clus',
        'n_clusters', 'flag'))
    for record in records:
        flags = []
        if record['K'] == best_k:
            flags.append('best')
        if record['K'] == selected_k:
            flags.append('selected')
        print('{0:>3}  {1:>10.2f}  {2:>8.4f}  {3:>8.3f}  {4:>10.3f}  {5:>9}  '
              '{6:>10}  {7}'.format(record['K'], record['cutoff'],
                                    record['abs_gap'], record['rel_gap'],
                                    record['silhouette'], record['min_size'],
                                    record['n_clusters'], ', '.join(flags)))


def cluster_size_map(labels, clusters):
    return {c: int((labels == c).sum()) for c in clusters}


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main(vcf_filename, pop_filename, output_filename, method, max_k,
         min_cluster_size, force_k, tree_mode, ordination,
         clone_list, clone_threshold, auto_clone, loci_filename,
         pops_from_sample_id, field_names, make_pdf, pdf_output):

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
    unit = 'genets' if removed else 'samples'
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

    print('\n###3 - K evaluation (per-K cut-off, gap and silhouette)')
    records, best_k = evaluate_k(linkage_matrix, dist, n_units, max_k,
                                 min_cluster_size)
    selected_k = best_k
    if force_k:
        selected_k = int(force_k)
        if selected_k not in [r['K'] for r in records]:
            sys.stderr.write('Warning: --k {0} not available; using best K {1}.'
                             '\n'.format(selected_k, best_k))
            selected_k = best_k
    selected = select_record(records, selected_k)
    display_max = max(r['K'] for r in records)
    print_k_table(records, best_k, selected_k)
    best_sil = select_record(records, best_k)['silhouette']
    print('Best K = {0} (silhouette {1:.3f}); selected K = {2}. Silhouette is '
          'the K-selection metric (higher = better-separated clusters); with '
          'overlapping lineages it stays modest. Use --k to pick another K.'
          .format(best_k, best_sil, selected_k))

    print('\n###4 - Cluster membership at selected K = {0}'.format(selected_k))
    labels = selected['labels']
    clusters = sorted(set(int(x) for x in labels))
    for cluster in clusters:
        members = sorted(nm for i, nm in enumerate(data['names'])
                         if int(labels[i]) == cluster)
        print('C{0} ({1}): {2}'.format(cluster, len(members),
                                       ', '.join(members)))
    if pop_filename:
        indivs_pops, _ = vcf_clone_detect.get_assignments_from_popfile(
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
    dxy = pairwise_dxy(stats)
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
    print('Pairwise dxy (absolute divergence):')
    for a, ca in enumerate(clusters):
        cells = []
        for b in range(len(clusters)):
            if a == b:
                cells.append('{0:>7}'.format('-'))
            else:
                cells.append('{0:>7}'.format('{0:.4f}'.format(dxy[a, b])
                             if dxy[a, b] == dxy[a, b] else 'NA'))
        print('C{0:>4} '.format(ca) + ' '.join(cells))
    print('Fixed differences:')
    for a, ca in enumerate(clusters):
        cells = ['{0:>7}'.format('-' if a == b else int(fixed_diff[a, b]))
                 for b in range(len(clusters))]
        print('C{0:>4} '.format(ca) + ' '.join(cells))

    # Optional ipyrad .loci shared/unique loci
    loci_presence = None
    if loci_filename:
        loci_presence, n_matched = parse_loci_presence(loci_filename,
                                                       data['names'])
        if n_matched < 3:
            sys.stderr.write('Warning: only {0} `.loci` sample(s) match the VCF;'
                             ' skipping locus analysis.\n'.format(n_matched))
            loci_presence = None
        else:
            lcl, lpriv, _ = locus_sharing_stats(loci_presence, labels,
                                                data['names'])
            print('Shared/unique loci ({0} loci, {1} samples matched):'.format(
                len(loci_presence), n_matched))
            for ci, c in enumerate(lcl):
                print('  C{0}: {1} private loci'.format(c, int(lpriv[ci])))

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
        category_tracks = build_category_tracks(data['names'], pop_filename,
                                                pops_from_sample_id, field_names)
        if category_tracks:
            print('Annotation tracks: {0}'.format(
                ', '.join('{0} ({1})'.format(t['name'], len(t['values']))
                          for t in category_tracks)))
        print('PDF report: {0}'.format(pdf_filename))
        try:
            write_pdf_report(sim, dist, linkage_matrix, data['names'], records,
                             display_max, selected, perc_genotyped, data,
                             method, tree_mode, ordination, pdf_filename,
                             category_tracks=category_tracks,
                             loci_presence=loci_presence)
        except Exception as error:
            sys.stderr.write('Warning: PDF report failed ({0})\n'.format(error))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-p', '--pop', dest='pop_filename', metavar='pop_file',
                        help='optional text file (tsv/csv): col 1 = sample, col 2 '
                             '= population, optional col 3 = lineage/species/'
                             'region. Drawn as dot-strip annotation tracks beside '
                             'the tree and cross-tabulated against clusters (not '
                             'used for clustering itself)')
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
                             'threshold above which samples are clones')
    parser.add_argument('--max-k', dest='max_k', type=int, default=DEF_MAX_K,
                        metavar='K', help='maximum K to evaluate (default: '
                        '{0})'.format(DEF_MAX_K))
    parser.add_argument('--min-cluster-size', dest='min_cluster_size', type=int,
                        default=DEF_MIN_CLUSTER_SIZE, metavar='N',
                        help='clusters smaller than this exclude a K from being '
                             'chosen as best (default: {0})'.format(
                                 DEF_MIN_CLUSTER_SIZE))
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
    parser.add_argument('--loci', dest='loci_filename', metavar='loci_file',
                        default=None,
                        help='optional ipyrad `.loci` file; adds a shared/unique '
                             'loci panel (private loci per cluster + pairwise '
                             'Jaccard of recovered loci)')
    parser.add_argument('--pops-from-sample-id', dest='pops_from_sample_id',
                        action='store_true',
                        help='derive up to three annotation tracks (drawn beside '
                             'the tree) from each sample name, splitting on "_" '
                             'and using the 2nd, 3rd and 4th fields (the 1st is '
                             'skipped); used instead of a popfile. A field with '
                             'more than {0} distinct values (e.g. a per-sample '
                             'id) is dropped'.format(MAX_CATEGORIES))
    parser.add_argument('--fields', dest='field_names', metavar='names',
                        default=None,
                        help='comma-separated names for the annotation tracks '
                             'in order (e.g. "location,depth"); overrides the '
                             'default track titles')
    parser.add_argument('--pdf-output', dest='pdf_output', default=None,
                        metavar='pdf_file', help='filename for the PDF report')
    parser.add_argument('--no-pdf', dest='no_pdf', action='store_true',
                        help='do not generate the PDF report (text only)')
    args = parser.parse_args()
    field_names = ([s.strip() for s in args.field_names.split(',')]
                   if args.field_names else None)
    main(args.vcf_filename, args.pop_filename, args.output_filename,
         args.method, args.max_k, args.min_cluster_size,
         args.force_k, args.tree_mode, args.ordination, args.clone_list,
         args.clone_threshold, args.auto_clone, args.loci_filename,
         args.pops_from_sample_id, field_names, make_pdf=not args.no_pdf,
         pdf_output=args.pdf_output)
