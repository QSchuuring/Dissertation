#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from polyleven import levenshtein
from sklearn.cluster import DBSCAN
from utils import *  # expects: load_fasta, encode_seqs, write_fasta, consensusVoting, lprint

# Optional imports are deferred until needed to avoid requiring plotting libs by default
plt = None
sns = None

def plot_landscape(x, y, df, query_, plot_type, args):
    global plt, sns
    if plt is None or sns is None:
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401

    plt.figure(figsize=(5, 5))
    tmp = df.loc[df.dbscan_label == -1]
    plt.scatter(tmp[x], tmp[y], color='lightgray', marker='x', label='unclustered')

    tmp = df.loc[df.dbscan_label > 9]
    plt.scatter(tmp[x], tmp[y], color='black', label='other clusters')

    tmp = df[(df.dbscan_label >= 0) & (df.dbscan_label <= 9)]
    sns.scatterplot(x=x, y=y, hue='dbscan_label', data=tmp, palette='tab10', linewidth=0)

    plt.scatter(query_[x], query_[y], color='red', marker='*', s=150, label='Ref Seq')
    plt.legend(bbox_to_anchor=(1, 1), frameon=False)

    plt.xlabel(x)
    plt.ylabel(y)
    plt.tight_layout()
    os.makedirs(args.o, exist_ok=True)
    plt.savefig(os.path.join(args.o, f"{args.keyword}_{plot_type}.pdf"), bbox_inches='tight')


def parse_uniform_sizes(expr: str):
    """
    Parse a sizes expression like:
      "10:30:2,70,100:150:10"
    -> [10,12,14,16,18,20,22,24,26,28,30,70,100,110,120,130,140,150]
    Accepts plain integers or ranges start:end[:step] (inclusive end).
    """
    sizes = []
    if not expr:
        return sizes
    for token in expr.split(','):
        token = token.strip()
        if not token:
            continue
        if ':' in token:
            parts = token.split(':')
            if len(parts) not in (2, 3):
                raise ValueError(f"Bad range token: {token}")
            start = int(parts[0])
            end = int(parts[1])
            step = int(parts[2]) if len(parts) == 3 else 1
            if step == 0:
                raise ValueError("Step cannot be 0")
            # inclusive end
            sizes.extend(list(range(start, end + (1 if step > 0 else -1), step)))
        else:
            sizes.append(int(token))
    # de-dup and sort
    return sorted(set(sizes))


if __name__ == '__main__':

    p = argparse.ArgumentParser(description="""
Cluster sequences in an MSA using DBSCAN and write .a3m files for each cluster.
Assumes first sequence in the input FASTA/A3M is the query sequence.

H Wayment-Steele, 2022 (lightly extended for uniform subsampling ranges)
    """.strip())

    p.add_argument("keyword", help="Keyword used to name outputs/logs.")
    p.add_argument("-i", required=True, help="Input FASTA/A3M alignment (query must be first).")
    p.add_argument("-o", required=True, help="Output directory.")
    p.add_argument("--n_controls", default=10, type=int,
                   help="Number of uniform subsample replicates to generate per size (default 10).")
    p.add_argument("--verbose", action='store_true', help="Print cluster info to stdout.")

    # DBSCAN / scan options
    p.add_argument('--scan', action='store_true',
                   help='Scan eps on 1/4 of sequences (shuffled) to pick value.')
    p.add_argument('--eps_val', type=float,
                   help="Fixed eps value for DBSCAN (overrides --scan).")
    p.add_argument('--resample', action='store_true',
                   help='If set, resample the original MSA with replacement before clustering.')
    p.add_argument("--gap_cutoff", type=float, default=0.25,
                   help='Remove sequences with > this fraction gaps (default 0.25).')
    p.add_argument('--min_eps', type=float, default=3.0,
                   help='Min eps to scan (default 3.0).')
    p.add_argument('--max_eps', type=float, default=20.0,
                   help='Max eps to scan (default 20.0).')
    p.add_argument('--eps_step', type=float, default=0.5,
                   help='Step for eps scan (default 0.5).')
    p.add_argument('--min_samples', type=int, default=3,
                   help='min_samples for DBSCAN (default 3; do not go lower).')

    # Dimensionality reduction plots
    p.add_argument('--run_PCA', action='store_true',
                   help='Run PCA on one-hot sequences and save plot + metadata.')
    p.add_argument('--run_TSNE', action='store_true',
                   help='Run t-SNE on one-hot sequences and save plot + metadata.')

    # Uniform subsampling controls
    p.add_argument('--uniform_sizes', type=str, default="10:30:2,70,100:150:10",
                   help=('Comma-separated sizes / ranges for uniform subsamples. '
                         'Format: "start:end[:step]" or single ints. '
                         'Default: "10:30:2,70,100:150:10" (U10–U30 step 2, U70, U100–U150 step 10).'))
    p.add_argument('--no_uniform', action='store_true',
                   help='Disable writing any uniform subsamples (only clustered outputs).')

    args = p.parse_args()

    # Lazy import for DR if requested
    if args.run_PCA:
        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401
    if args.run_TSNE:
        from sklearn.manifold import TSNE
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401

    os.makedirs(args.o, exist_ok=True)
    log_f = open(f"{args.keyword}.log", 'w')

    # Load and clean MSA
    IDs, seqs = load_fasta(args.i)
    # keep only aligned uppercase and gaps
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]

    df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs})

    # First row is the query
    query_ = df.iloc[:1].copy()
    df = df.iloc[1:].copy()

    if args.resample:
        df = df.sample(frac=1.0, replace=True, random_state=None).reset_index(drop=True)

    L = len(df.sequence.iloc[0])
    N = len(df)

    # Filter by gap fraction
    df['frac_gaps'] = [s.count('-') / L for s in df['sequence']]
    former_len = len(df)
    df = df.loc[df.frac_gaps < args.gap_cutoff].copy()
    new_len = len(df)

    lprint(args.keyword, log_f)
    lprint(f"{former_len - new_len} seqs removed for > {int(args.gap_cutoff * 100)}% gaps, {new_len} remaining", log_f)

    # One-hot encode sequences (for DBSCAN and optional DR plots)
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)

    # ---- eps selection ----
    if args.eps_val is None and args.scan:
        lprint('eps\tn_clusters\tn_not_clustered', log_f)
        eps_test_vals = np.arange(args.min_eps, args.max_eps + args.eps_step, args.eps_step)
        n_clusters_obs = []

        for eps in eps_test_vals:
            # 1/4 sample for quick scan
            sample = df.sample(frac=0.25, random_state=None)
            sample_ohe = encode_seqs(sample.sequence.tolist(), max_len=L)
            clustering = DBSCAN(eps=float(eps), min_samples=args.min_samples).fit(sample_ohe)
            labels = clustering.labels_
            n_clust = len(set(labels)) - (1 if -1 in labels else 0)
            n_not = int(np.sum(labels == -1))
            lprint(f'{eps:.2f}\t{n_clust}\t{n_not}', log_f)
            n_clusters_obs.append(n_clust)
            if eps > 10 and n_clust == 1:
                break

        eps_to_select = float(eps_test_vals[int(np.argmax(n_clusters_obs))])
    else:
        eps_to_select = float(args.eps_val) if args.eps_val is not None else 7.0  # reasonable default if no scan

    # ---- DBSCAN on full set ----
    clustering = DBSCAN(eps=eps_to_select, min_samples=args.min_samples).fit(ohe_seqs)
    df['dbscan_label'] = clustering.labels_

    lprint(f'Selected eps={eps_to_select:.2f}', log_f)
    lprint(f'{len(df)} total seqs', log_f)

    clusters = sorted([int(x) for x in df.dbscan_label.unique() if x >= 0])
    unclustered = int((df.dbscan_label == -1).sum())
    lprint(f'{len(clusters)} clusters, {unclustered} of {len(df)} not clustered ({unclustered/len(df):.2f})', log_f)

    # Distances to query
    avg_id_uncl = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L
                           for x in df.loc[df.dbscan_label == -1, 'sequence'].tolist()]) if unclustered else float('nan')
    lprint(f'avg identity to query of unclustered: {avg_id_uncl:.2f}', log_f)

    avg_id_cl = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L
                         for x in df.loc[df.dbscan_label != -1, 'sequence'].tolist()]) if (len(df) - unclustered) else float('nan')
    lprint(f'avg identity to query of clustered: {avg_id_cl:.2f}', log_f)

    # ---- Write clustered MSAs and gather metadata ----
    cluster_metadata = []
    for clust in clusters:
        tmp = df.loc[df.dbscan_label == clust].copy()

        cs = consensusVoting(tmp.sequence.tolist())
        avg_dist_to_cs = np.mean([1 - levenshtein(x, cs) / L for x in tmp.sequence.tolist()])
        avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in tmp.sequence.tolist()])

        if args.verbose:
            print(f'Cluster {clust} consensus seq, {len(tmp)} seqs:')
            print(cs)
            print('#########################################')
            for _, row in tmp.iterrows():
                print(row['SequenceName'], row['sequence'])
            print('#########################################')

        out_df = pd.concat([query_, tmp], axis=0)
        out_path = os.path.join(args.o, f"{args.keyword}_{clust:03d}.a3m")
        write_fasta(out_df.SequenceName.tolist(), out_df.sequence.tolist(), outfile=out_path)

        cluster_metadata.append({
            'cluster_ind': clust,
            'consensusSeq': cs,
            'avg_lev_dist': f'{avg_dist_to_cs:.3f}',
            'avg_dist_to_query': f'{avg_dist_to_query:.3f}',
            'size': len(out_df)
        })

    # ---- Uniform subsampling exactly as in the dissertation ----
    if not args.no_uniform:
        target_sizes = parse_uniform_sizes(args.uniform_sizes)
        for size in target_sizes:
            if len(df) < size:
                lprint(f"skipping U{size}: requested size exceeds available sequences ({len(df)} < {size})", log_f)
                continue

            lprint(f"writing {args.n_controls} size-{size} uniformly sampled clusters", log_f)
            for i in range(args.n_controls):
                tmp = df.sample(n=size, replace=args.resample, random_state=None)
                out_df = pd.concat([query_, tmp], axis=0)
                outfile = os.path.join(args.o, f"{args.keyword}_U{size}-{i:03d}.a3m")
                write_fasta(out_df.SequenceName.tolist(), out_df.sequence.tolist(), outfile=outfile)

    # ---- Optional DR plots ----
    if args.run_PCA:
        lprint('Running PCA ...', log_f)
        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401

        ohe_vecs = encode_seqs(df.sequence.tolist(), max_len=L)
        mdl = PCA()
        embedding = mdl.fit_transform(ohe_vecs)
        query_embedding = mdl.transform(encode_seqs(query_.sequence.tolist(), max_len=L))

        df['PC 1'] = embedding[:, 0]
        df['PC 2'] = embedding[:, 1]
        query_['PC 1'] = query_embedding[:, 0]
        query_['PC 2'] = query_embedding[:, 1]

        plot_landscape('PC 1', 'PC 2', df, query_, 'PCA', args)
        lprint('Saved PCA plot to ' + os.path.join(args.o, f"{args.keyword}_PCA.pdf"), log_f)

    if args.run_TSNE:
        lprint('Running TSNE ...', log_f)
        from sklearn.manifold import TSNE
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401

        # t-SNE has no .transform, so fit on all + query together
        ohe_vecs = encode_seqs(df.sequence.tolist() + [query_['sequence'].iloc[0]], max_len=L)
        mdl = TSNE()
        embedding = mdl.fit_transform(ohe_vecs)

        df['TSNE 1'] = embedding[:-1, 0]
        df['TSNE 2'] = embedding[:-1, 1]
        query_['TSNE 1'] = embedding[-1:, 0]
        query_['TSNE 2'] = embedding[-1:, 1]

        plot_landscape('TSNE 1', 'TSNE 2', df, query_, 'TSNE', args)
        lprint('Saved TSNE plot to ' + os.path.join(args.o, f"{args.keyword}_TSNE.pdf"), log_f)

    # ---- Write metadata ----
    clust_outfile = os.path.join(args.o, f"{args.keyword}_clustering_assignments.tsv")
    df.to_csv(clust_outfile, index=False, sep='\t')
    lprint(f'wrote clustering data to {clust_outfile}', log_f)

    metad_outfile = os.path.join(args.o, f"{args.keyword}_cluster_metadata.tsv")
    metad_df = pd.DataFrame.from_records(cluster_metadata)
    metad_df.to_csv(metad_outfile, index=False, sep='\t')
    lprint(f'wrote cluster metadata to {metad_outfile}', log_f)

    lprint(f"Saved this output to {args.keyword}.log", log_f)
    log_f.close()
