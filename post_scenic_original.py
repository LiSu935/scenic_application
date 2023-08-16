# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scanpy.plotting._tools.scatterplots import plot_scatter
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='pyscenic pipeline after pyscenic step4')
parser.add_argument('--prefix', type=str, default='SULI', help='prefix: prefix for the output files')
parser.add_argument('--path_ad_file_from_seurat_integrated', type=str, default="/storage/htc/joshilab/Su_Li/Spencerlab/data/c1_integrated.h5ad", help='SCT mtx from seurat containing all features')
parser.add_argument('--wdir', type=str, default="/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/", help='working dir of the scenic, major dir.')

args = parser.parse_args()
print(args)

# here to set a prefix for ease
prefix = args.prefix

# set a working directory

wdir = args.wdir
wdir = wdir+"results/"+prefix+"/"
os.chdir( wdir )

if not os.path.exists('top5_regulator_module_regulon_list'):
   os.makedirs('top5_regulator_module_regulon_list')


# =====No need if already QC by Seurat================================================= #
# ===================================================================================== #
# path to unfiltered loom file (this will be created in the optional steps below)
#f_loom_path_unfilt = "pbmc10k_unfiltered.loom" # test dataset, n=500 cells
# ===================================================================================== #

# path to loop file that has been QC and clustering by Seurat:
#f_anndata_path_input = "/storage/htc/joshilab/Su_Li/StowersHSC/scenic_application/seurat_scenicInput/hsc_slim.h5ad"

# path to unfiltered loom file (this will be created in the optional steps below)
#f_loom_path_unfilt = prefix+"_unfiltered.loom"

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "hsc_integrated.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
# replaced with h5ad file converted from Seurat:
f_anndata_path = args.path_ad_file_from_seurat_integrated

# path to pyscenic output
f_pyscenic_output = prefix+"_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = prefix+'_scenic_integrated-output.loom'

reg_path = prefix+"reg.csv"
adj_path = prefix+"adj.csv"



sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=150)

# scenic output
lf = lp.connect( f_final_loom, mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).head()

print(lf.ra.keys())
print(lf.ca.keys())

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
    [
        
        pd.DataFrame( lf.ca.ClusterID, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.Louvain_clusters_Scanpy, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
    ],
    axis=1
)
cellAnnot.columns = [

 'ClusterID',
 'Louvain_clusters_Scanpy',
 'Percent_mito',
 'nGene',
 'nUMI']

cellAnnot

# capture embeddings:
dr = [
    pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
]
dr_names = [
    meta['embeddings'][0]['name'].replace(" ","_")
]

# add other embeddings
drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )

for i in range( len(drx.columns) ):
    dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
    dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )

# rename columns:
for i,x in enumerate( dr ):
    x.columns = ['X','Y']

lf.close()

## Display a motifs table with motif logos

# View the motifs table along with motif logos

# helper functions (not yet integrated into pySCENIC):

from pyscenic.utils import load_motifs
import operator as op
from IPython.display import HTML, display

BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"

def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', 200)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)

df_motifs = load_motifs(reg_path)

selected_motifs = ['ATF4','TCF3','EBF1']
df_motifs_sel = df_motifs.iloc[ [ True if x in selected_motifs else False for x in df_motifs.index.get_level_values('TF') ] ,:]

#display_logos(df_motifs.head())
display_logos( df_motifs_sel.sort_values([('Enrichment','NES')], ascending=False).head(9))


## Regulon specificity scores (RSS) across predicted cell types

from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

### Calculate RSS

rss_louvain = regulon_specificity_scores( auc_mtx, cellAnnot['Louvain_clusters_Scanpy'] )
rss_louvain

cats = sorted( list(set(cellAnnot['Louvain_clusters_Scanpy'])), key=int )

fig = plt.figure(figsize=(15, 12))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_louvain.T[c]
    ax = fig.add_subplot(4,5,num)
    plot_rss(rss_louvain, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig(prefix+"_cluster-RSS-top5.tiff", dpi=150, bbox_inches = "tight")

# Select the top 5 regulons from each cell type
topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_louvain.T[c].sort_values(ascending=False)[:5].index)
    )
topreg = list(set(topreg))

topreg_o = []
for i,c in enumerate(cats):
    topreg_o.extend(
        list(rss_louvain.T[c].sort_values(ascending=False)[:5].index)
    )

print("This is the original topreg_o list:")
print(topreg_o) 

# Generate a Z-score for each regulon to enable comparison between regulons
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
#auc_mtx_Z.sort_index(inplace=True)

# Generate a heatmap
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in cellAnnot['Louvain_clusters_Scanpy'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats, size=1.0)

sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')    
g.ax_heatmap.set_xlabel('')  
plt.savefig(prefix+"_cluster-heatmap-top5.tiff", dpi=600, bbox_inches = "tight")


### Further exploration of modules directly from the network inference output

adjacencies = pd.read_csv(adj_path, index_col=False, sep=',')

from pyscenic.utils import modules_from_adjacencies
modules = list(modules_from_adjacencies(adjacencies, exprMat))



def write_tf_module_regulon_list(tf):
    tf_mods = [ x for x in modules if x.transcription_factor==tf ]
    for i,mod in enumerate( tf_mods ):
        with open( "./top5_regulator_module_regulon_list/"+tf+'_module_'+str(i)+'.txt', 'w') as f:
            for item in mod.genes:
                f.write("%s\n" % item)
            
    with open( "./top5_regulator_module_regulon_list/"+tf+'_regulon.txt', 'w') as f:
        for item in regulons[tf+'_(+)']:
            f.write("%s\n" % item)

            
            
for tf in topreg:
    tf = tf[:-4]
    write_tf_module_regulon_list(tf)
    
            
