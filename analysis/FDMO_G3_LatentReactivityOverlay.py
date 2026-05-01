#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Aidan Cosgrove <aidancos@umich.edu>"

"""
FDMO_G3_LatentReactivityOverlay.py

FDMO Gen3 2D latent space with reactivity overlay for all three generations.

  - Training sequences    : very faint gray background
  - Library seqs inactive : pastel pink
  - Library seqs active   : pink -> purple -> dark blue gradient by reactivity

Uses posterior mean (mu) directly from the G3 VAE training pkl.
"""

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

################################################
# Paths
################################################

LATENT_PKL = os.path.expanduser(
    '~/KJ_G3Paper/Models/FDMO_G3_latent2_2026-04-10_19-49-31_6594459/'
    'FDMO_G3_NoCV_vae_d2_layer4_w0.0005_b256_l0.0005_beta1_300epoch_seed19_latent.pkl'
)
REACTIVITY_CSV = os.path.expanduser(
    '~/KJ_G3Paper/KJ_G3PAPER/master_reactivity_clean.csv'
)
PLATE_XLSX = os.path.expanduser(
    '~/KJ_G3Paper/KJ_G3PAPER/MasterReactivity.xlsx'
)
OUT_DIR = os.path.expanduser('~/KJ_G3Paper/figures/latent_space_overlays')
os.makedirs(OUT_DIR, exist_ok=True)

################################################
# Color scheme
################################################

INACTIVE_COLOR = '#FFD6DE'
INACTIVE_EDGE  = '#F0A0B0'
TRAIN_COLOR    = '#E8E8E8'
CMAP = mcolors.LinearSegmentedColormap.from_list(
    'reactivity',
    ['#FFB3C6', '#C77DFF', '#3A0CA3', '#115472'], N=256
)

ROW_LABELS = list('ABCDEFGH')

################################################
# Load latent space
################################################

print("Loading latent space...")
with open(LATENT_PKL, 'rb') as f:
    latent = pickle.load(f)

all_keys  = latent['key'] + latent['t_key']
all_mu    = np.vstack([latent['mu'], latent['t_mu']])
key_to_mu = {k: all_mu[i] for i, k in enumerate(all_keys)}
latent_set = set(all_keys)

# Separate training background from library sequences
anc_keys   = set(k for k in all_keys if 'Anc' in str(k) or 'AncFDMO' in str(k))
train_keys = [k for k in all_keys if k not in anc_keys]
train_mu   = np.array([key_to_mu[k] for k in train_keys])

print(f"  Total sequences  : {len(all_keys):,}")
print(f"  Training (bg)    : {len(train_keys):,}")
print(f"  Ancestor/library : {len(anc_keys):,}")

################################################
# Load plate layouts from MasterReactivity.xlsx
################################################

print("\nLoading plate layouts...")
xl = pd.read_excel(PLATE_XLSX, sheet_name=None)
print(f"  Sheets: {list(xl.keys())}")

def parse_layout(sheet_name):
    """Extract {(row_label, col_num): enzyme_name} from a plate layout sheet."""
    if sheet_name not in xl:
        print(f"  WARNING: sheet '{sheet_name}' not found")
        return {}
    df      = xl[sheet_name]
    mapping = {}
    # find header row with column numbers
    header_row = None
    for i, row in df.iterrows():
        nums = [v for v in row if isinstance(v, (int, float)) and 1 <= v <= 12]
        if len(nums) >= 10:
            header_row = row
            break

    for i, row in df.iterrows():
        row_label = row.iloc[0]
        if row_label not in ROW_LABELS:
            continue
        for j in range(1, len(row)):
            val = row.iloc[j]
            if pd.notna(val) and isinstance(val, str) and len(val.strip()) > 1:
                name = val.strip()
                if name.lower() in ('blk', 'blank', 'empty', ''):
                    continue
                col_num = j  # fallback to column index
                if header_row is not None:
                    try:
                        col_num = int(float(header_row.iloc[j]))
                    except (ValueError, TypeError):
                        col_num = j
                mapping[(row_label, col_num)] = name
    return mapping

# Try common sheet name variants
def find_sheet(candidates):
    for c in candidates:
        if c in xl:
            return c
    # case-insensitive fallback
    for c in candidates:
        for s in xl:
            if c.lower() == s.lower():
                return s
    return None

sheet_g1 = find_sheet(['Plate_layout_g1', 'Gen 1', 'gen1', 'G1'])
sheet_g2 = find_sheet(['Plate_layout_g2', 'Gen 2', 'gen2', 'G2'])
sheet_g3 = find_sheet(['Plate_layout_g3', 'Gen 3', 'gen3', 'G3',
                        'Gen 3 pt1', 'gen3p1'])

print(f"  Gen1 sheet: {sheet_g1}")
print(f"  Gen2 sheet: {sheet_g2}")
print(f"  Gen3 sheet: {sheet_g3}")

layout_g1 = parse_layout(sheet_g1) if sheet_g1 else {}
layout_g2 = parse_layout(sheet_g2) if sheet_g2 else {}
layout_g3 = parse_layout(sheet_g3) if sheet_g3 else {}

print(f"  Gen1 wells mapped: {len(layout_g1)}")
print(f"  Gen2 wells mapped: {len(layout_g2)}")
print(f"  Gen3 wells mapped: {len(layout_g3)}")

################################################
# Build enzyme -> latent key mapping
################################################

def find_latent_key(enzyme_name):
    """Try several ID variants to find a match in the latent space."""
    candidates = [
        enzyme_name,
        enzyme_name.replace(' ', ''),
        enzyme_name.replace(' ', '_'),
        enzyme_name + '_v2',
        enzyme_name + '_v1',
        enzyme_name.replace(' ', '') + '_v2',
        enzyme_name.replace(' ', '') + '_v1',
        # capitalise first letter of suffix (e.g. _hypo. -> _Hypo.)
        enzyme_name[0].upper() + enzyme_name[1:],
    ]
    for c in candidates:
        if c in latent_set:
            return c
    # case-insensitive fallback
    lower_map = {k.lower(): k for k in latent_set}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None

enzyme_to_latent = {}
for enzyme in set(layout_g1.values()):
    key = find_latent_key(enzyme)
    if key:
        enzyme_to_latent[('gen1', enzyme)] = key

for enzyme in set(layout_g2.values()):
    key = find_latent_key(enzyme)
    if key:
        enzyme_to_latent[('gen2', enzyme)] = key

for enzyme in set(layout_g3.values()):
    key = find_latent_key(enzyme)
    if key:
        enzyme_to_latent[('gen3', enzyme)] = key

all_lib_keys = set(enzyme_to_latent.values())

print(f"\n  Gen1 in latent space: {sum(1 for (g,_) in enzyme_to_latent if g=='gen1')}")
print(f"  Gen2 in latent space: {sum(1 for (g,_) in enzyme_to_latent if g=='gen2')}")
print(f"  Gen3 in latent space: {sum(1 for (g,_) in enzyme_to_latent if g=='gen3')}")
print(f"  Total unique         : {len(all_lib_keys)}")

################################################
# Load reactivity
################################################

print("\nLoading reactivity data...")
react_df = pd.read_csv(REACTIVITY_CSV)
print(f"  Substrates: {sorted(react_df['substrate'].unique())}")

def get_reactivity(substrate_name):
    """
    Returns {latent_key: max_reactivity} for all library sequences
    for a given substrate. Matches directly on enzyme name column.
    Zero for inactive/missing.
    """
    sub    = react_df[react_df['substrate'] == substrate_name]
    result = {k: 0.0 for k in all_lib_keys}

    for (gen, enzyme), lkey in enzyme_to_latent.items():
        match = sub[
            (sub['generation'] == gen) &
            (sub['enzyme'] == enzyme)
        ]
        if len(match) > 0:
            result[lkey] = max(result[lkey],
                               float(match['reactivity'].max()))
    return result


################################################
# Substrates to plot — use exact names from master_reactivity.csv
################################################

SUBSTRATES = [
    ('2-Ph-indole',              '2-Ph-Indole'),
    ('3-Indole-acetamide',       '3-Indole-Acetamide'),  # was 3-Me-Indole in old CSV
    ('Tryptophol',               'Tryptophol'),
    ('Gramine',                  'Gramine'),
    ('Melatonin',                'Melatonin'),
    ('Methyl-quinoline',         'Methyl-Quinoline'),
    ('nPr-benzaldehyde',         'nPr-Benzaldehyde'),
    ('p-OH-benzaldehyde',        'p-OH-Benzaldehyde'),
    ('Phenox',                   'Phenox'),
    ('Aniline',                  'Aniline'),
    ('TropB-native',             'TropB Native'),
    ('Tetralone',                'Tetralone'),
]

################################################
# Plot function
################################################

def plot_panel(ax, anc_react, name, vmax, show_cbar=False, cbar_ax=None):
    # Training background
    ax.scatter(train_mu[:, 0], train_mu[:, 1],
               c=TRAIN_COLOR, s=1.5, alpha=0.08,
               linewidths=0, rasterized=True, zorder=1)

    inactive = [k for k in all_lib_keys if anc_react.get(k, 0) <= 0]
    active   = [(k, anc_react[k]) for k in all_lib_keys if anc_react.get(k, 0) > 0]

    if inactive:
        anc_inactive   = [k for k in inactive if k in key_to_mu and 'Anc' in str(k)]
        extant_inactive = [k for k in inactive if k in key_to_mu and 'Anc' not in str(k)]
        if anc_inactive:
            coords = np.array([key_to_mu[k] for k in anc_inactive])
            ax.scatter(coords[:, 0], coords[:, 1],
                       c=INACTIVE_COLOR, s=20, alpha=0.85, marker='o',
                       edgecolors=INACTIVE_EDGE, linewidths=0.3, zorder=2)
        if extant_inactive:
            coords = np.array([key_to_mu[k] for k in extant_inactive])
            ax.scatter(coords[:, 0], coords[:, 1],
                       c=INACTIVE_COLOR, s=25, alpha=0.85, marker='^',
                       edgecolors=INACTIVE_EDGE, linewidths=0.3, zorder=2)

    sc = None
    if active:
        anc_active    = [(k, v) for k, v in active if k in key_to_mu and 'Anc' in str(k)]
        extant_active = [(k, v) for k, v in active if k in key_to_mu and 'Anc' not in str(k)]
        sc = None
        if anc_active:
            keys, scores = zip(*anc_active)
            coords = np.array([key_to_mu[k] for k in keys])
            sc = ax.scatter(coords[:, 0], coords[:, 1],
                            c=np.array(scores), cmap=CMAP, vmin=0, vmax=vmax,
                            s=55, alpha=0.95, marker='o',
                            edgecolors='white', linewidths=0.4, zorder=3)
        if extant_active:
            keys, scores = zip(*extant_active)
            coords = np.array([key_to_mu[k] for k in keys])
            sc_ext = ax.scatter(coords[:, 0], coords[:, 1],
                                c=np.array(scores), cmap=CMAP, vmin=0, vmax=vmax,
                                s=65, alpha=0.95, marker='^',
                                edgecolors='white', linewidths=0.4, zorder=3)
            if sc is None:
                sc = sc_ext
        if show_cbar and cbar_ax is not None and sc is not None:
            cb = plt.colorbar(sc, cax=cbar_ax)
            cb.set_label('Reactivity (prod/IS)', fontsize=8)
            cb.ax.tick_params(labelsize=7)

    n_active = len(active)
    ax.set_title(name, fontsize=8.5, fontweight='bold', pad=14)
    ax.text(0.5, 1.02,
            f'({n_active} active / {len(all_lib_keys)} total)',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=7, fontstyle='italic', color='#888888')
    ax.set_xlabel('$Z_1$', fontsize=7, color='#888888')
    ax.set_ylabel('$Z_2$', fontsize=7, color='#888888')
    ax.tick_params(labelsize=6, colors='#888888')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#bbbbbb')
    ax.spines['bottom'].set_color('#bbbbbb')
    return sc

################################################
# Build figures
################################################

print("\nComputing reactivity per substrate...")
all_reacts = {raw: get_reactivity(raw) for raw, _ in SUBSTRATES}
all_scores = [v for r in all_reacts.values() for v in r.values() if v > 0]
global_vmax = np.percentile(all_scores, 95) if all_scores else 1.0
print(f"  Global vmax (95th pct): {global_vmax:.4f}")
for raw, label in SUBSTRATES:
    n_active = sum(1 for v in all_reacts[raw].values() if v > 0)
    print(f"  {label:<25} : {n_active} active")

# ── Grid figure ───────────────────────────────────────────────────────
print("\nBuilding grid figure...")
ncols, nrows = 4, 3
fig = plt.figure(figsize=(ncols * 4.5, nrows * 4.2 + 1.0))
gs  = GridSpec(nrows, ncols + 1,
               width_ratios=[1] * ncols + [0.04],
               hspace=0.65, wspace=0.35,
               top=0.93, bottom=0.08)
cbar_ax = fig.add_subplot(gs[:, -1])

last_sc = None
for idx, (raw, name) in enumerate(SUBSTRATES):
    ax = fig.add_subplot(gs[idx // ncols, idx % ncols])
    sc = plot_panel(ax, all_reacts[raw], name, global_vmax,
                    show_cbar=(idx == len(SUBSTRATES) - 1),
                    cbar_ax=cbar_ax)
    if sc is not None:
        last_sc = sc

if last_sc is not None:
    cb = plt.colorbar(last_sc, cax=cbar_ax)
    cb.set_label('Reactivity (prod/IS)', fontsize=8)
    cb.ax.tick_params(labelsize=7)

legend_handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=TRAIN_COLOR,
           markersize=7, label='Training sequences', alpha=0.6,
           markeredgecolor='#cccccc', markeredgewidth=0.5),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=INACTIVE_COLOR,
           markeredgecolor=INACTIVE_EDGE, markersize=9,
           label='Inactive ancestor'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor=INACTIVE_COLOR,
           markeredgecolor=INACTIVE_EDGE, markersize=9,
           label='Inactive extant'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#115472',
           markeredgecolor='white', markersize=9,
           label='Active ancestor'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='#115472',
           markeredgecolor='white', markersize=9,
           label='Active extant'),
]
fig.legend(handles=legend_handles,
           loc='lower center', ncol=5, fontsize=9,
           frameon=True, framealpha=0.9,
           bbox_to_anchor=(0.47, -0.01),
           borderpad=0.8, handletextpad=0.5)

fig.suptitle('FDMO Library — Latent Space Reactivity (Gen1 + Gen2 + Gen3)',
             fontsize=14, fontweight='bold', y=1.02)

out_grid = os.path.join(OUT_DIR, 'FDMO_G3_latent_reactivity_grid.png')
plt.savefig(out_grid, dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print(f"  Saved: {out_grid}")

# ── Individual panels ─────────────────────────────────────────────────
print("Building individual panels...")
for raw, name in SUBSTRATES:
    fig2, ax2 = plt.subplots(figsize=(5.5, 4.8))
    cbar_ax2  = fig2.add_axes([0.92, 0.15, 0.025, 0.65])
    sc2 = plot_panel(ax2, all_reacts[raw], name, global_vmax,
                     show_cbar=True, cbar_ax=cbar_ax2)
    if sc2 is None:
        # no active seqs — still draw colorbar scale
        import matplotlib.cm as cm
        sm = plt.cm.ScalarMappable(cmap=CMAP,
                                   norm=mcolors.Normalize(0, global_vmax))
        sm.set_array([])
        cb2 = plt.colorbar(sm, cax=cbar_ax2)
        cb2.set_label('Reactivity (prod/IS)', fontsize=8)
        cb2.ax.tick_params(labelsize=7)
    plt.subplots_adjust(right=0.88, top=0.88)
    safe = name.replace(' ','_').replace('/','_').replace('+','and')
    out_ind = os.path.join(OUT_DIR, f'FDMO_G3_latent_{safe}.png')
    plt.savefig(out_ind, dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: FDMO_G3_latent_{safe}.png")

print(f"\nAll figures saved to: {OUT_DIR}")
