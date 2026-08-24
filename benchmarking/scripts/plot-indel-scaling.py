#!/usr/bin/env python3
"""The scaling figure: PEPMatch 2.0 vs Brute Force from 100 to 1,000,000 queries.

Reads the per-size tables written by slurm/run-synth-scaling.sbatch
(results/synth-scaling/synth_<n>indel_<label>_benchmarking.tsv) and renders a
two-panel, print-quality figure.

  left   total time vs query count, log-log, with search-only as a dashed companion
  right  recall vs query count, log-x

WHAT THE FIGURE IS ALLOWED TO CLAIM -- two things are annotated on the plot rather
than left for a reader to trip over:

  1. PEPMatch's TOTAL includes a one-time proteome index build (~29 s, rebuilt every run
     because Preprocessor has no cache). At small N that constant dominates and Brute
     Force genuinely finishes sooner -- measured, not hypothetical. The crossover is a
     real property of the method, so the figure marks it instead of hiding it behind a
     search-only plot.
  2. Brute Force's recall is 100% BY CONSTRUCTION: it wraps the same committed oracle
     that generated the ground truth. It is drawn as a labelled reference line, not as
     evidence. What carries weight is PEPMatch matching it while scaling.

Colors are the validated categorical slots 1 (blue) and 2 (orange) -- checked with the
dataviz validator: adjacent-pair CVD dE 24.7 (protan), normal-vision 33.6, both well
above threshold. Identity is never carried by color alone: each series also has its own
marker and dash pattern, which is what keeps the figure readable in a grayscale print.

Usage
  plot-indel-scaling.py                       # 1-indel, results/synth-scaling
  plot-indel-scaling.py --indels 1 --outdir results/figures
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')          # headless: must precede pyplot's own backend selection
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]

# Validated categorical slots (see references/palette.md); ink/grid from the same system.
BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, MUTED, GRID, AXIS = '#0b0b0b', '#898781', '#e1e0d9', '#c3c2b7'

SERIES = [
  # label,         tsv Method value, color,  marker, linestyle
  ('PEPMatch 2.0', 'PEPMatch',       BLUE,   'o',    '-'),
  ('Brute Force',  'Brute Force',    ORANGE, 's',    '--'),
]

SIZE_LABELS = [('100', 100), ('1k', 1_000), ('5k', 5_000), ('10k', 10_000),
               ('100k', 100_000), ('1m', 1_000_000)]


def load(results_dir, indels):
  """-> {method: {N: row}}, plus the list of (N, reason) points that are absent.

  Missing and FAILED points are collected and PRINTED rather than interpolated over: a
  gap in a scaling curve that silently closes itself is indistinguishable from a
  measurement, and this figure exists to make a claim about scaling.
  """
  data, missing = {}, []
  for label, n in SIZE_LABELS:
    path = results_dir / f'synth_{indels}indel_{label}_benchmarking.tsv'
    if not path.exists():
      missing.append((n, f'no table at {path.name}'))
      continue
    with open(path, newline='') as f:
      for row in csv.DictReader(f, delimiter='\t'):
        if row.get('Status', 'OK') != 'OK':
          missing.append((n, f'{row["Method"]} status={row["Status"]}'))
          continue
        data.setdefault(row['Method'], {})[n] = row
  return data, missing


def series_xy(data, method, column):
  xs, ys = [], []
  for _, n in SIZE_LABELS:
    row = data.get(method, {}).get(n)
    if not row:
      continue
    try:
      ys.append(float(row[column]))
      xs.append(n)
    except (KeyError, ValueError):
      continue          # 'N/A' preprocessing etc.
  return xs, ys


def style_axes(ax):
  ax.set_axisbelow(True)
  ax.grid(True, which='major', color=GRID, linewidth=0.8)
  ax.grid(True, which='minor', color=GRID, linewidth=0.4, alpha=0.6)
  for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
  for side in ('left', 'bottom'):
    ax.spines[side].set_color(AXIS)
    ax.spines[side].set_linewidth(0.8)
  ax.tick_params(colors=MUTED, labelsize=8, width=0.8)
  for lbl in ax.get_xticklabels() + ax.get_yticklabels():
    lbl.set_color(MUTED)


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--indels', type=int, default=1)
  p.add_argument('--results-subdir', default='synth-scaling')
  p.add_argument('--outdir', default=None)
  args = p.parse_args()

  results_dir = BENCH / 'results' / args.results_subdir
  outdir = Path(args.outdir) if args.outdir else BENCH / 'results' / 'figures'
  outdir.mkdir(parents=True, exist_ok=True)

  data, missing = load(results_dir, args.indels)
  if not data:
    raise SystemExit(f'no usable tables in {results_dir}')

  print(f'loaded methods: {", ".join(sorted(data))}')
  for method in sorted(data):
    print(f'  {method:14s} N = {sorted(data[method])}')
  if missing:
    print('\nMISSING / NOT-OK points (plotted as gaps, never interpolated):')
    for n, why in missing:
      print(f'  N={n:<9,} {why}')

  plt.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,      # embed real fonts, keep text selectable
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
  })
  fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9.2, 3.9))

  # ---- left: time ------------------------------------------------------------------
  for label, method, color, marker, dash in SERIES:
    xs, ys = series_xy(data, method, 'Total (s)')
    if xs:
      ax_t.plot(xs, ys, color=color, marker=marker, linestyle=dash, linewidth=2,
                markersize=6.5, markeredgecolor='white', markeredgewidth=0.9,
                label=f'{label} — total', zorder=3)
    xs, ys = series_xy(data, method, 'Searching (s)')
    if xs:
      ax_t.plot(xs, ys, color=color, marker=marker, linestyle=':', linewidth=1.5,
                markersize=4.5, markerfacecolor='white', markeredgecolor=color,
                alpha=0.85, label=f'{label} — search only', zorder=2)

  ax_t.set_xscale('log')
  ax_t.set_yscale('log')
  ax_t.set_xlabel('Query count', fontsize=9, color=INK)
  ax_t.set_ylabel('Time (s)', fontsize=9, color=INK)
  ax_t.set_title('Search time scales with query count', fontsize=10, color=INK,
                 loc='left', pad=8)
  ax_t.xaxis.set_minor_formatter(NullFormatter())
  ax_t.yaxis.set_minor_locator(LogLocator(base=10, subs=tuple(range(2, 10)), numticks=20))
  ax_t.yaxis.set_minor_formatter(NullFormatter())
  style_axes(ax_t)
  ax_t.legend(fontsize=7.5, frameon=False, labelcolor=MUTED, loc='upper left')

  # The honest caveat, on the figure rather than in a caption nobody reads.
  ax_t.text(0.98, 0.03,
            "PEPMatch total includes a one-time\n~29 s index build (dotted = search only)",
            transform=ax_t.transAxes, fontsize=6.8, color=MUTED,
            ha='right', va='bottom', linespacing=1.35)

  # ---- right: recall ---------------------------------------------------------------
  for label, method, color, marker, dash in SERIES:
    xs, ys = series_xy(data, method, 'Recall (%)')
    if xs:
      ax_r.plot(xs, ys, color=color, marker=marker, linestyle=dash, linewidth=2,
                markersize=6.5, markeredgecolor='white', markeredgewidth=0.9,
                label=label, zorder=3)

  ax_r.set_xscale('log')
  ax_r.set_ylim(0, 105)
  ax_r.set_xlabel('Query count', fontsize=9, color=INK)
  ax_r.set_ylabel('Recall (%)', fontsize=9, color=INK)
  ax_r.set_title('Recall is unaffected by scale', fontsize=10, color=INK, loc='left', pad=8)
  ax_r.xaxis.set_minor_formatter(NullFormatter())
  style_axes(ax_r)
  ax_r.legend(fontsize=7.5, frameon=False, labelcolor=MUTED, loc='lower left')
  ax_r.text(0.98, 0.06,
            'Brute Force is the reference:\n100% by construction',
            transform=ax_r.transAxes, fontsize=6.8, color=MUTED,
            ha='right', va='bottom', linespacing=1.35)

  fig.tight_layout()
  stem = outdir / f'synth-{args.indels}indel-scaling'
  for ext in ('pdf', 'png'):
    fig.savefig(f'{stem}.{ext}', dpi=300, bbox_inches='tight',
                facecolor='white')
    print(f'wrote {stem}.{ext}')
  plt.close(fig)


if __name__ == '__main__':
  main()
