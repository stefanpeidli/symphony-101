import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns
import yaml
import PyComplexHeatmap as pch

# make figures paper quality
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['font.family'] = 'Arial'

import warnings
warnings.filterwarnings("ignore")

with open('../configuration/config.yaml', "r") as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
cell_line_dict=config['cell_lines']

cell_lines = ['H1299', 'H1975', 'A549']  # exclude SYM4 and SYM5 in determining thresholds
all_guide_rankings = pd.concat({cell_line: pd.read_csv(f'../results/mageck/summaries/{cell_line}/gene_SYM101_counts_mageck.sgrna_summary.txt', sep='\t', index_col=0) for cell_line in cell_line_dict.values()}, axis=1)
# aggregate
all_guide_rankings[('summary', 'max_score')] = all_guide_rankings.loc[:, (slice(None), 'score')].max(axis=1)
all_guide_rankings[('summary', 'mean_score')] = all_guide_rankings.loc[:, (slice(None), 'score')].mean(axis=1)
all_guide_rankings[('summary', 'min_score')] = all_guide_rankings.loc[:, (slice(None), 'score')].min(axis=1)
all_guide_rankings[('summary', 'max_LFC')] = all_guide_rankings.loc[:, (slice(None), 'LFC')].max(axis=1)
all_guide_rankings[('summary', 'mean_LFC')] = all_guide_rankings.loc[:, (slice(None), 'LFC')].mean(axis=1)
all_guide_rankings[('summary', 'min_LFC')] = all_guide_rankings.loc[:, (slice(None), 'LFC')].min(axis=1)
all_guide_rankings[('summary', 'absmax_LFC')] = np.abs(all_guide_rankings.loc[:, (slice(None), 'LFC')]).max(axis=1)
all_guide_rankings[('summary', 'is_control')] = all_guide_rankings.index.str.contains('Control').astype(str)

data = all_guide_rankings.loc[:, 'summary']
data_control = data[(data.is_control=='True')]
y_key = 'max_score'
p = 50
with sns.axes_style('whitegrid'):
    fig, ax = pl.subplots(figsize=(6,4))
sns.scatterplot(data=data, x='absmax_LFC', y=y_key, edgecolor=None, label='Targeting guide')
sns.scatterplot(data=data_control, x='absmax_LFC', y=y_key, color='red', label='NT guide')
tresh = np.percentile(data_control[y_key], p)
pl.axhline(tresh, color='r', linestyle='--', label=f'{p}th percentile of NT guides')
percent_removed = ((data[y_key]<tresh) & (data.is_control=='False')).mean()
pl.text(0.5, 0.25, f'{percent_removed:.1%} guides removed', transform=pl.gca().transAxes)
pl.title('MAGeCK sgRNA test results\nMax across lines, All Time Points')
pl.ylabel('MAGeCK score (Max across cell lines)')
pl.xlabel('Max Absolute LFC across cell lines')
pl.legend()
pl.savefig(snakemake.output['filtering_qc_plot'], dpi=300, bbox_inches='tight')
pl.close()

# annotate
all_guide_rankings[('summary', 'is_working')] = all_guide_rankings.loc[:, 'summary'].max_score>tresh
# How many guides are marked per target?
data = all_guide_rankings.loc[:, 'summary']
# add infos
df = pd.read_csv(snakemake.input['counts'], index_col=0)
for col in ['Gene', 'Transcript_ID', 'transcript_biotype', 'Region']:
    data[col] = df.loc[data.index, col]
# counting
sdf = data.groupby(['Transcript_ID', 'transcript_biotype'], observed=False).sum().is_working.reset_index()
sdf = sdf[(sdf.transcript_biotype.isin(['protein_coding', 'lncRNA', 'NT_Control']))]
sdf.rename(columns={'is_working':'good_guides_count'}, inplace=True)
pl.figure(figsize=(8,3))
sns.barplot(data=sdf.value_counts(['transcript_biotype', 'good_guides_count']).reset_index(),
            x='good_guides_count', y=0, hue='transcript_biotype')
pl.xlabel('Number of good / non-filtered guides')
pl.ylabel('Number of transcripts')
pl.axvline(2.5, color='k')
pl.savefig(snakemake.output['filtering_hist_plot'], dpi=300, bbox_inches='tight')
pl.close()

# filtering would remove too many guides for these target, so we enforce keeping at least 3 guides
N = 3
key = 'keep_flag'
too_few = sdf.loc[sdf.good_guides_count<N, 'Transcript_ID']
print(len(too_few))
data = all_guide_rankings.loc[:, 'summary']
# add infos
df = pd.read_csv(snakemake.input['counts'], index_col=0)
for col in ['Gene', 'Transcript_ID', 'transcript_biotype', 'Region']:
    data[col] = df.loc[data.index, col]
keep = data[data.Transcript_ID.isin(too_few)].groupby('Transcript_ID').max_score.nlargest(N).index.get_level_values(1)
all_guide_rankings.loc[:, ('summary', key)] = all_guide_rankings.loc[:, ('summary', 'is_working')].copy()
all_guide_rankings.loc[df.loc[all_guide_rankings.index].Transcript_ID.isin(too_few), ('summary', key)] = False
all_guide_rankings.loc[keep, ('summary', key)] = True


# visualize in heatmap
df = pd.read_csv(snakemake.input['counts'], index_col=0)
sample_names = np.sort(df.columns[12:])
pools = ['Sym-pool-1-3', 'Sym-pool-1-1', 'Sym-pool-1-2']
pools = np.sort(pools)
sample_names_wo_pools = [x for x in sample_names if not any([p in x for p in pools])]
controls = [x for x in df.index if 'control' in x.lower()]
# First I normalize for seq-depth by bringing all NT guides to the same count
target = df.loc[controls, sample_names_wo_pools].mean().mean()  # mean of NTs in control
ist = df.loc[controls, sample_names].mean()  # mean of NTs in all
norm_factor = target / ist  # make all NT mean counts equal to those in control
ndf = df.loc[:, sample_names] * norm_factor  # normalize
# Foldchange to Pool after seq-depth normalization based on NT depth
fc = (ndf.loc[:, sample_names].T / ndf.loc[:, pools].mean(axis=1)).T  # normalize to control (pools)
l2fc = np.log2(fc)
l2fc[l2fc==-np.infty] = None
df['mean_l2fc'] = l2fc.mean(axis=1)
biotype = df.sort_values(['transcript_biotype', 'mean_l2fc']).transcript_biotype
biotype = biotype[~biotype.isin(['protein_coding_CDS_not_defined', 'retained_intron'])]
is_working = all_guide_rankings.loc[:, 'summary'].loc[data.index].sort_values(by=[key, 'max_score'])[key].astype(str).map({'True': 'Good guides', 'False': 'Bad guides'})
cell_line_dict_ = {f'0{x}': y for x, y in cell_line_dict.items()}
cell_line = pd.Series({x: x.split('-')[1] for x in sample_names}).replace(cell_line_dict_)
replicate = pd.Series({x: x.split('-')[2] for x in sample_names})
time_point = pd.Series({x: x.split('-')[3].split('_')[0] for x in sample_names})
col_order = pd.concat([cell_line, replicate, time_point], axis=1).sort_values([0, 2, 1]).index  # sort by time before replicate
# select and reorder
data = l2fc.loc[:, col_order]
top_ha = pch.HeatmapAnnotation(
    cell_line=pch.anno_simple(cell_line.loc[col_order], add_text=True, height=7, cmap='Pastel1'),
    time_point=pch.anno_simple(time_point.loc[col_order], add_text=True, height=7, cmap='Set1'),
    replicate=pch.anno_simple(replicate.loc[col_order], add_text=False),
    verbose=0, axis=1, plot_legend=False
)
left_ha = pch.HeatmapAnnotation(
    biotype=pch.anno_simple(biotype, legend=True, height=7, add_text=False, cmap='Pastel2'),
    Filtered=pch.anno_simple(is_working, legend=True, height=7, add_text=True, cmap='Set1'),
    verbose=0, axis=0, plot_legend=True
)
pl.figure(figsize=(15, 10))
cm = pch.ClusterMapPlotter(
    data=data, 
    left_annotation=left_ha, top_annotation=top_ha,
    row_cluster=False, col_cluster=False,
    row_split=is_working,
    col_split=cell_line.loc[col_order],
    col_split_gap=2, row_split_gap=1,
    label='beta', row_dendrogram=False, col_dendrogram=False,
    legend_gap=7, rasterized=True, cmap='rocket', vmin=-3, vmax=0,
    legend_kws=dict(label='L2FC'),
    xlabel="Samples", ylabel='Guides', legend_hpad=0,
    verbose=0
)
pl.suptitle('L2FC Matrix for SYMPHONY pooled screen 1', y=0.92)
for plot in snakemake.output['filtered_heatmap']:
    pl.savefig(plot, dpi=300, bbox_inches='tight')
pl.close()

# Save filtered count object
df = pd.read_csv(snakemake.input['counts'], index_col=0)
good_guides = all_guide_rankings.index[all_guide_rankings.loc[:, 'summary'].keep_flag]
tdf = df.loc[good_guides]
tdf.to_csv(snakemake.output['filtered_counts'])
all_guide_rankings.to_csv(snakemake.output['filtered_tab'], sep='\t')