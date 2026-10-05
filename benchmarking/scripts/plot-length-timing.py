#!/usr/bin/env python3
"""Figure 5: PEPMatch 2.0 search time per query, as a function of query length.

Reads the tables written by time-per-length.py and draws two panels, one per indel count
(1 left, 2 right), on a SHARED log y-axis so the two indel counts compare directly.

  PEPMatch 2.0, per-length k   solid line, filled circles   (each length at its own k)
  PEPMatch 2.0, fixed k        dashed line, open triangles  (k=4 / k=3, Figure 4's k)
  Brute Force (naive)          dotted line, squares         (64-thread equivalent)

The y value is the ENGINE marginal cost, (T_10k - T_1k) / 9,000 queries, in ms. The
harness timer (Matcher.match, Figure 4's definition) is not plotted: it adds
single-threaded Python per hit and per call, which is not the algorithm. The seed length
k is printed under each per-length point, and alternate k steps are lightly shaded, so a
reader can see the curve step where k changes.

Colors are Figure 4's: blue = PEPMatch 2.0, orange = Brute Force. Identity is never
color alone (marker + dash differ for every series), so a grayscale print still reads.

  plot-length-timing.py --results results/per-length --outdir results/figures
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]

BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, MUTED, GRID, AXIS = '#0b0b0b', '#898781', '#e1e0d9', '#c3c2b7'
BAND = '#f4f3ee'


def read_tsv(path):
  with open(path) as f:
    return list(csv.DictReader(f, delimiter='\t'))


def style_axes(ax):
  ax.set_axisbelow(True)
  ax.grid(True, which='major', axis='y', color=GRID, linewidth=0.8)
  ax.grid(True, which='minor', axis='y', color=GRID, linewidth=0.4, alpha=0.6)
  for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
  for side in ('left', 'bottom'):
    ax.spines[side].set_color(AXIS)
    ax.spines[side].set_linewidth(0.8)
  ax.tick_params(colors=MUTED, labelsize=8, width=0.8)


def main():
  p = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
  p.add_argument('--results', default=str(BENCH / 'results' / 'per-length'))
  p.add_argument('--outdir', default=str(BENCH / 'results' / 'figures'))
  args = p.parse_args()
  results, outdir = Path(args.results), Path(args.outdir)
  outdir.mkdir(parents=True, exist_ok=True)

  summary = read_tsv(results / 'per-length-summary.tsv')
  bf_path = results / 'per-length-bf.tsv'
  bf = read_tsv(bf_path) if bf_path.exists() else []

  bad = [r for r in summary if r['recall_min'] and float(r['recall_min']) != 100.0]
  if bad:
    raise SystemExit(f'Refusing to plot: {len(bad)} sets below 100% recall, '
                     f'e.g. d={bad[0]["indels"]} L={bad[0]["length"]}')

  plt.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
  })
  indels = sorted({int(r['indels']) for r in summary})
  fig, axes = plt.subplots(1, len(indels), figsize=(4.6 * len(indels), 3.9),
                           sharey=True, squeeze=False)

  lows = []
  for ax, d in zip(axes[0], indels):
    def series(arm):
      rows = sorted((r for r in summary if int(r['indels']) == d and r['arm'] == arm
                     and r['engine_marginal_ms_per_query']),
                    key=lambda r: int(r['length']))
      return ([int(r['length']) for r in rows],
              [float(r['engine_marginal_ms_per_query']) for r in rows],
              [int(r['k']) for r in rows])

    xs, ys, ks = series('per-length')

    # Shade alternate k steps so the step structure is visible without reading labels.
    if xs:
      start, shade = xs[0], False
      for i in range(1, len(xs) + 1):
        if i == len(xs) or ks[i] != ks[i - 1]:
          if shade:
            ax.axvspan(start - 0.5, xs[i - 1] + 0.5, color=BAND, zorder=0, lw=0)
          shade = not shade
          if i < len(xs):
            start = xs[i]

    fx, fy, fk = series('fixed-k')
    if fx:
      fixed = max(set(fk), key=fk.count)
      ax.plot(fx, fy, color=BLUE, marker='^', markerfacecolor='white', linestyle='--',
              linewidth=1.6, markersize=6, label=f'PEPMatch 2.0, fixed k = {fixed}',
              zorder=3)
    if xs:
      ax.plot(xs, ys, color=BLUE, marker='o', linestyle='-', linewidth=2, markersize=6.5,
              markeredgecolor='white', markeredgewidth=0.9,
              label='PEPMatch 2.0, k per length', zorder=4)
      for x, y, k in zip(xs, ys, ks):
        # Below the point: the fixed-k point at the same length is never faster (a
        # shorter seed only adds candidates), so above is where it would collide.
        ax.annotate(f'{k}', (x, y), textcoords='offset points', xytext=(0, -12),
                    ha='center', fontsize=6.5, color=MUTED)

    brows = sorted((r for r in bf if int(r['indels']) == d), key=lambda r: int(r['length']))
    if brows:
      ax.plot([int(r['length']) for r in brows],
              [float(r['ms_per_query_at_threads']) for r in brows],
              color=ORANGE, marker='s', linestyle=':', linewidth=2, markersize=6,
              markeredgecolor='white', markeredgewidth=0.9,
              label='Brute Force (naive)', zorder=4)

    ax.set_yscale('log')
    lows += ys + fy + [float(r['ms_per_query_at_threads']) for r in brows]
    ax.set_xlim(7.4, 25.6)
    ax.set_xticks(range(8, 26, 2))
    ax.set_xlabel('Query length (residues)', fontsize=9, color=INK)
    ax.set_title(f'{d} indel{"s" if d > 1 else ""}', fontsize=10, color=INK, loc='left',
                 pad=8)
    style_axes(ax)
    ax.legend(fontsize=7.2, frameon=False, labelcolor=MUTED, loc='upper right')

  if lows:
    # Shared y-axis: set ONCE from the lowest point in either panel, or the last panel's
    # limit would clip the other panel's k labels.
    axes[0][0].set_ylim(bottom=min(lows) / 3)
  axes[0][0].set_ylabel('Search time per query (ms)', fontsize=9, color=INK)
  info = dict(line.rstrip('\n').split('\t', 1) for line in open(results / 'run-info.txt'))
  fig.text(0.995, 0.005, 'Numbers under points: seed length k. Engine-only marginal time, '
           f'{info["threads"]} threads.', fontsize=6.5, color=MUTED, ha='right', va='bottom')
  fig.tight_layout(rect=(0, 0.03, 1, 1))
  stem = outdir / 'per-length-timing'
  for ext in ('pdf', 'png'):
    fig.savefig(f'{stem}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    print(f'wrote {stem}.{ext}')
  plt.close(fig)


if __name__ == '__main__':
  main()
