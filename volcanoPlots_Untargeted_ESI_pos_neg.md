#Volcano plots generation for Figures 1, S4, S5, S6#



import pandas as pd
import seaborn as sns 
from matplotlib import pyplot as plt,cm,lines
from matplotlib import colormaps
list(colormaps)
import numpy as np
CB_color_cycle = {
    'blue':    '#377eb8',
    'orange':  '#ff7f00',
    'green':   '#4daf4a',
    'pink':    '#f781bf',
    'brown':   '#a65628',
    'purple':  '#984ea3',
    'gray':    '#999999',
    'red':     '#e41a1c',
    'yellow':  '#dede00',
} 

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
from matplotlib import cm
import pandas as pd

# Example dataframe (replace this with your actual data loading step)
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")
# Define key column names
padj = 'Adj. P-value: (Bb_pos_12022024_Newanalysis) / (CM_pos_12022024_Newanalysis)'
l2fc = 'Log2 Fold Change: (Bb_pos_12022024_Newanalysis) / (CM_pos_12022024_Newanalysis)'

# Rename columns for clarity
df = df.rename(columns={l2fc: 'l2fc'})

# Compute log-adjusted p-values and significance
df['logpadj'] = -np.log10(df[padj])
df['significance'] = (df[padj] < 0.05) & (abs(df.l2fc) >= 2)

# Filter for valid names
name_idx = df['Name'].dropna().index
new_df = df.loc[name_idx]
new_df['color'] = new_df.significance.apply(lambda x: 'DarkGray' if not x else 'gainsboro')

# Drop unnecessary columns
columnstodrop = [c for c in new_df.columns if ':' in c] + ['Tags', 'Checked']
sig_df = new_df[(new_df[padj] < 0.05) & (new_df['l2fc'] >= 2)].drop(columnstodrop, axis=1)
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0])

# Assign colors to lipid classes (customized color palette)
class_colors = {
    'TG': '#FF00FF',      # dark blue
    'PG': '#4B0082',      # indigo
    'PE': '#4B0082',      # violet
    'PC': '#4B0082',      # purple
    'MGDG': '#FF00FF',    # magenta
    'MePC': '#4B0082',    # pink
    'LPC': '#4B0082',     # salmon
    'LBPA': '#4B0082',    # orange     
    'DG': '#FF00FF',      # gold
    'AcHexCmE': "#2287AE", # yellow
    'ZyE': "#2287AE",      # forest green
    'WE': "#2287AE",       # dark turquoise
    'SPH': '#8B0000',     # aquamarine
    'SM': '#8B0000',      # light sea green
    'MG': '#FF00FF',      # magenta
    'Hex1Cer': '#8B0000', # green-yellow
    'ChE': "#2287AE",     # firebrick red
    'Car': "#2287AE",     # crimson
    'CarE': "#2287AE",    # dark red
    'BisMeLPA': "#4B0082", # slate gray
    'AcHexChE': "#2287AE", # silver
    'Co': '#2287AE',       # orange red
    'AcCa': "#2287AE",     # slate blue
    'PI': '#4B0082',       # lime green
    'Cer': '#8B0000',      # saddle brown
    'SiE': "#2287AE"       # deep pink
}
sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 5))

sns.scatterplot(data=new_df, x='l2fc', y='logpadj', ax=ax, color='gainsboro')
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.8)

# Add threshold lines
ax.hlines(-np.log10(0.05), -10, 10, colors='black', linestyles='dashed', linewidth=0.8)
ax.set_xlim(-10, 10)
ax.set_xticks([-10, -5, 0, 5, 10])
#ax.vlines([-2, 2], 0, 10, colors='black', linestyles='dashed', linewidth=0.8)

# Labels and title
#ax.set_xlabel('$log_2$(fold change)')#
#ax.set_ylabel('$-log_{10}$(adjusted p-value)')#
#ax.set_title(r'$\mathbf{\mathit{Borrelia\ burgdorferi}}$ vs. complete media (positive mode)', fontsize=10)

# Custom legend
#handles = [lines.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=lipid) for lipid, color in class_colors.items()]
#ax.legend(handles=handles, title='Lipid Class', frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')
ax.legend().remove()
plt.tight_layout()
plt.savefig('03052026_BbvsCM_Pos.svg', dpi=600)
plt.show()






'''''''''''''''''''''''''''''''''''''
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
from matplotlib import cm
import pandas as pd

# Example dataframe (replace this with your actual data loading step)

# Define key column names
padj_col = 'Adj. P-value: (Bb_pos_12022024_Newanalysis) / (SM_pos_12022024_Newanalyis)'
l2fc_col = 'Log2 Fold Change: (Bb_pos_12022024_Newanalysis) / (SM_pos_12022024_Newanalyis)'

# Rename columns for clarity
df = df.rename(columns={l2fc: 'l2fc'})

# Compute log-adjusted p-values and significance
df['logpadj'] = -np.log10(df[padj])
df['significance'] = (df[padj] < 0.05) & (abs(df.l2fc) >= 2)

# Filter for valid names
name_idx = df['Name'].dropna().index
new_df = df.loc[name_idx]
new_df['color'] = new_df.significance.apply(lambda x: 'DarkGray' if not x else 'gainsboro')

# Drop unnecessary columns
columnstodrop = [c for c in new_df.columns if ':' in c] + ['Tags', 'Checked']
sig_df = new_df[(new_df[padj] < 0.05) & (new_df['l2fc'] >= 2)].drop(columnstodrop, axis=1)
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0])

# Assign colors to lipid classes (customized color palette)
class_colors = {
    'TG': '#FF00FF',      # dark blue
    'PG': '#4B0082',      # indigo
    'PE': '#4B0082',      # violet
    'PC': '#4B0082',      # purple
    'MGDG': '#FF00FF',    # magenta
    'MePC': '#4B0082',    # pink
    'LPC': '#4B0082',     # salmon
    'LBPA': '#4B0082',    # orange
    'DG': '#FF00FF',      # gold
    'AcHexCmE':"#2287AE", # yellow
}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 5))

sns.scatterplot(data=new_df, x='l2fc', y='logpadj', ax=ax, color='gainsboro')
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.8)

# Add threshold lines
ax.hlines(-np.log10(0.05), -10, 10, colors='black', linestyles='dashed', linewidth=0.8)
ax.set_xlim(-10, 10)
ax.set_xticks([-10, -5, 0, 5, 10])
#ax.vlines([-2, 2], 0, 10, colors='black', linestyles='dashed', linewidth=0.8)

# Labels and title
#ax.set_xlabel('$log_2$(fold change)')#
#ax.set_ylabel('$-log_{10}$(adjusted p-value)')#
#ax.set_title(r'$\mathbf{\mathit{Borrelia\ burgdorferi}}$ vs. complete media (positive mode)', fontsize=10)

# Custom legend
#handles = [lines.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=lipid) for lipid, color in class_colors.items()]
#ax.legend(handles=handles, title='Lipid Class', frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')
ax.legend().remove()
plt.tight_layout()
plt.savefig('03052026_BbvsSM_Pos.svg', dpi=600)
plt.show()







''''''''''''''''''''''''''''''''''''''''''''''''''''
#Complete media vs borrelia (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (Bb_neg_newanalysis_12022024) / (CM_neg_newanalysis_12022024)'
l2fc_col = 'Log2 Fold Change: (Bb_neg_newanalysis_12022024) / (CM_neg_newanalysis_12022024)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]
#Assign colours
# === Assign colors to lipid classes ===
lipid_colors = {
    'TG': '#FF00FF',       # magenta
    'PG': '#4B0082',       # indigo
    'PE': '#4B0082',       # violet
    'PC': '#4B0082',       # purple
    'MGDG': '#FF00FF',     # magenta
    'MePC': '#4B0082',     # pink
    'LPC': '#4B0082',      # salmon
    'LBPA': '#4B0082',     # orange
    'DG': '#FF00FF',       # gold
    'AcHexCmE': '#2287AE', # yellow
    'ZyE': '#2287AE',      # forest green
    'WE': '#2287AE',       # dark turquoise
    'SPH': '#8B0000',      # aquamarine
    'SM': '#8B0000',       # light sea green
    'MG': '#FF00FF',       # magenta
    'Hex1Cer': '#8B0000',  # green-yellow
    'ChE': '#2287AE',      # firebrick red
    'Car': '#2287AE',      # crimson
    'CarE': '#2287AE',     # dark red
    'BisMeLPA': '#4B0082', # slate gray
    'AcHexChE': '#2287AE', # silver
    'Co': '#2287AE',       # orange red
    'AcCa': '#2287AE',     # slate blue
    'PI': '#4B0082',       # lime green
    'Cer': '#8B0000',      # saddle brown
    'SiE': '#2287AE',      # deep pink
    'SQDG': '#FF00FF',     # magenta
    'PS': '#4B0082',
    'FA': "#820016",
    'LdMePE': "#4B0082"       # dark turquoise
}

# Map colors to lipid classes in the significant DataFrame
sig_df['color'] = sig_df['lipid_class'].map(lipid_colors)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 5))

# Scatter plots for significant and non-significant data
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='gainsboro')
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', 
                palette=lipid_colors, s=60, edgecolor='white', linewidth=0.8)

# Add threshold lines for significance
ax.hlines(-np.log10(0.05), -10, 10, colors='black', linestyles='dashed', linewidth=0.8)
ax.set_xlim(-10, 10)
ax.set_xticks([-10, -5, 0, 5, 10])

# Labels and title
#ax.set_xlabel('$log_2$(fold change)')
#ax.set_ylabel('$-log_{10}$(adjusted p-value)')
#ax.set_title(r'$\mathbf{\mathit{Borrelia\ burgdorferi}}$ vs. complete media (positive mode)', fontsize=10)#

# Custom legend
#handles = [lines.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                         #markersize=8, label=lipid) for lipid, color in lipid_colors.items()]
#ax.legend(handles=handles, title='Lipid Class', frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')

ax.legend().remove()  # Not needed as we don't add a legend
plt.tight_layout()
plt.tight_layout()
plt.savefig('BbvsCM_NEGATIVE_filtered_0305.svg',dpi=600)






'''''''''''''''''''''''''''''''''''''''''''''''''''''
#spent media vs borrelia (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (Bb_neg_newanalysis_12022024) / (SM_neg_newanalysis_12022024)'
l2fc_col = 'Log2 Fold Change: (Bb_neg_newanalysis_12022024) / (SM_neg_newanalysis_12022024)'
# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

#Assign colours
# === Assign colors to lipid classes ===
lipid_colors = {
    'TG': '#FF00FF',       # magenta
    'PG': '#4B0082',       # indigo
    'PE': '#4B0082',       # violet
    'PC': '#4B0082',       # purple
    'MGDG': '#FF00FF',     # magenta
    'MePC': '#4B0082',     # pink
    'LPC': '#4B0082',      # salmon
    'LBPA': '#4B0082',     # orange
    'DG': '#FF00FF',       # gold
    'AcHexCmE': '#2287AE', # yellow
    'ZyE': '#2287AE',      # forest green
    'WE': '#2287AE',       # dark turquoise
    'SPH': '#8B0000',      # aquamarine
    'SM': '#8B0000',       # light sea green
    'MG': '#FF00FF',       # magenta
    'Hex1Cer': '#8B0000',  # green-yellow
    'ChE': '#2287AE',      # firebrick red
    'Car': '#2287AE',      # crimson
    'CarE': '#2287AE',     # dark red
    'BisMeLPA': '#4B0082', # slate gray
    'AcHexChE': '#2287AE', # silver
    'Co': '#2287AE',       # orange red
    'AcCa': '#2287AE',     # slate blue
    'PI': '#4B0082',       # lime green
    'Cer': '#8B0000',      # saddle brown
    'SiE': '#2287AE',      # deep pink
    'SQDG': '#FF00FF',     # magenta
    'PS': '#4B0082',
    'FA': "#820016",
    'LdMePE': "#4B0082"       # dark turquoise
}

# Map colors to lipid classes in the significant DataFrame
sig_df['color'] = sig_df['lipid_class'].map(lipid_colors)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 5))

# Scatter plots for significant and non-significant data
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='gainsboro')
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', 
                palette=lipid_colors, s=60, edgecolor='white', linewidth=0.8)

# Add threshold lines for significance
ax.hlines(-np.log10(0.05), -10, 10, colors='black', linestyles='dashed', linewidth=0.8)
ax.set_xlim(-10, 10)
ax.set_xticks([-10, -5, 0, 5, 10])

# Labels and title
#ax.set_xlabel('$log_2$(fold change)')
#ax.set_ylabel('$-log_{10}$(adjusted p-value)')
#ax.set_title(r'$\mathbf{\mathit{Borrelia\ burgdorferi}}$ vs. complete media (positive mode)', fontsize=10)#

# Custom legend
#handles = [lines.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                         #markersize=8, label=lipid) for lipid, color in lipid_colors.items()]
#ax.legend(handles=handles, title='Lipid Class', frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')

ax.legend().remove()  # Not needed as we don't add a legend
plt.tight_layout()
plt.tight_layout()
plt.savefig('BbvsSPM_NEGATIVE_filtered_0305.svg',dpi=600)



'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# Load data
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")

# Define key columns
padj_col = 'Adj. P-value: (Rabbitserum_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)'
l2fc_col = 'Log2 Fold Change: (Rabbitserum_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)'

# Preprocess data
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (df['l2fc'] >= 2)

# Filter valid entries and add lipid classes
valid_df = df[df['Name'].notna()].copy()
valid_df['lipid_class'] = valid_df['Name'].str.split('(').str[0]

# Create significance-filtered dataframe
sig_df = valid_df[(valid_df[padj_col] < 0.05) & 
                 (valid_df['l2fc'] >= 2)].copy()

# Define color palette
class_colors = {
    'TG': '#00008B', 'PG': '#4B0082', 'PE': '#8A2BE2', 'PC': '#800080',
    'MGDG': '#FF00FF', 'MePC': '#FF69B4', 'LPC': '#FA8072', 'LBPA': '#FFA500',
    'DG': '#FFD700', 'AcHexCmE': '#FFFF00', 'ZyE': '#228B22', 'WE': '#00CED1',
    'SPH': '#7FFFD4', 'SM': '#20B2AA', 'MG': '#008000', 'Hex1Cer': '#ADFF2F',
    'ChE': '#B22222', 'Car': '#DC143C', 'CarE': '#8B0000', 'BisMeLPA': '#708090',
    'AcHexChE': '#C0C0C0', 'Co': '#FF4500', 'AcCa': '#6A5ACD', 'PI': '#32CD32',
    'Cer': '#8B4513', 'SiE': '#FF1493', 'SQDG': "#138B87", 'PS': '#00FF7F',
    'FA': '#FF6347', 'LdMePE': "#667F93", 'ST': '#FFB6C1', 'LPE': '#FF1493',
    'PA': '#ADFF2F', 'PEt': "#864B36", 'OAHFA': "#634282", 'dMePE': '#4682B4',
    'CmE': '#FF8C00'
}

# Create dynamic color mapping based on present lipids
present_classes = sorted(sig_df['lipid_class'].unique())
class_colors_filtered = {k: v for k, v in class_colors.items() 
                        if k in present_classes}

# Create plot
fig, ax = plt.subplots(figsize=(10, 10))

# Background points
sns.scatterplot(
    data=valid_df,
    x='l2fc',
    y='logpadj',
    color='gainsboro',
    ax=ax,
    label='Non-significant'
)

# Significant points
sns.scatterplot(
    data=sig_df,
    x='l2fc',
    y='logpadj',
    hue='lipid_class',
    palette=class_colors_filtered,
    s=60,
    edgecolor='white',
    linewidth=0.4,
    ax=ax
)

# Threshold lines
ax.hlines(-np.log10(0.05), -10, 10, 
         colors='black', linestyles='dashed', linewidth=0.8)
ax.vlines([2], 0, 10,  # Only positive L2FC threshold
         colors='black', linestyles='dashed', linewidth=0.8)

# Custom legend
legend_elements = [
    lines.Line2D([0], [0],
                 marker='o',
                 color='w',
                 markerfacecolor=color,
                 markersize=10,
                 label=lipid)
    for lipid, color in class_colors_filtered.items()
]

# Add legend with criteria in title
ax.legend(
    handles=legend_elements,
    title='Lipid Classes\n(padj ≤ 0.05, L2FC ≥ 2)',
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    handletextpad=0.2,
    fontsize=9
)

# Axis labels and title
ax.set_xlabel('$log_2$(fold change)', fontsize=12)
ax.set_ylabel('$-log_{10}$(adjusted p-value)', fontsize=12)
ax.set_title('Rabbit Serum vs. Water Control (Positive Mode)', fontsize=14)

plt.tight_layout()
plt.savefig(
    'rabbit_serum_vs_water_control_positivemode.svg',
    bbox_inches='tight',
    dpi=600
)
plt.show()





'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#complete media vs water control positive mode

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# Load data
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")

# Define key columns
padj_col = 'Adj. P-value: (CM_pos_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)'
l2fc_col = 'Log2 Fold Change: (CM_pos_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)'

# Preprocess data
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (df['l2fc'] >= 2)

# Filter valid entries and add lipid classes
valid_df = df[df['Name'].notna()].copy()
valid_df['lipid_class'] = valid_df['Name'].str.split('(').str[0]

# Create significance-filtered dataframe
sig_df = valid_df[(valid_df[padj_col] < 0.05) & 
                 (valid_df['l2fc'] >= 2)].copy()

# Define color palette (expanded)
class_colors = {
    'TG': '#00008B', 'PG': '#4B0082', 'PE': '#8A2BE2', 'PC': '#800080',
    'MGDG': '#FF00FF', 'MePC': '#FF69B4', 'LPC': '#FA8072', 'LBPA': '#FFA500',
    'DG': '#FFD700', 'AcHexCmE': '#FFFF00', 'ZyE': '#228B22', 'WE': '#00CED1',
    'SPH': '#7FFFD4', 'SM': '#20B2AA', 'MG': '#008000', 'Hex1Cer': '#ADFF2F',
    'ChE': '#B22222', 'Car': '#DC143C', 'CarE': '#8B0000', 'BisMeLPA': '#708090',
    'AcHexChE': '#C0C0C0', 'Co': '#FF4500', 'AcCa': '#6A5ACD', 'PI': '#32CD32',
    'Cer': '#8B4513', 'SiE': '#FF1493', 'SQDG': "#138B87", 'PS': '#00FF7F',
    'FA': '#FF6347', 'LdMePE': '#4682B4', 'ST': '#FFB6C1', 'LPE': '#FF1493',
    'PA': '#ADFF2F', 'PEt': '#FF4500', 'OAHFA': '#8A2BE2', 'dMePE': '#4682B4',
    'CmE': '#FF8C00'
}

# Create dynamic color mapping based on present lipids
present_classes = sorted(sig_df['lipid_class'].unique())
class_colors_filtered = {k: v for k, v in class_colors.items() 
                        if k in present_classes}

# Create plot
fig, ax = plt.subplots(figsize=(10, 10))

# Background points
sns.scatterplot(
    data=valid_df,
    x='l2fc',
    y='logpadj',
    color='gainsboro',
    ax=ax,
    label='Non-significant'
)

# Significant points
sns.scatterplot(
    data=sig_df,
    x='l2fc',
    y='logpadj',
    hue='lipid_class',
    palette=class_colors_filtered,
    s=60,
    edgecolor='white',
    linewidth=0.4,
    ax=ax
)

# Threshold lines
ax.hlines(-np.log10(0.05), -10, 10, 
         colors='black', linestyles='dashed', linewidth=0.8)
ax.vlines([-2, 2], 0, 10, 
         colors='black', linestyles='dashed', linewidth=0.8)

# Create custom legend
legend_elements = [
    lines.Line2D([0], [0],
                 marker='o',
                 color='w',
                 markerfacecolor=color,
                 markersize=10,
                 label=lipid)
    for lipid, color in class_colors_filtered.items()
]

# Add legend with sorting
ax.legend(
    handles=legend_elements,
    title='Lipid Classes\n(padj ≤ 0.05, L2FC ≥ 2)',
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    handletextpad=0.2,
    fontsize=9
)

# Axis labels
ax.set_xlabel('Log2 Fold Change', fontsize=12)
ax.set_ylabel('-log10(Adj. P-value)', fontsize=12)
ax.set_xlim(-2, max(valid_df['l2fc'].max() + 1, 3))
ax.set_ylim(0, max(valid_df['logpadj'].max() + 1, 5))

plt.tight_layout()
plt.savefig(
    'complete_media_vs_water_control_positivemode_filtered_12222025.svg',
    bbox_inches='tight',
    dpi=600
)
plt.show()





''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#spent media vs water control positive mode
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# =========================
# Load data
# =========================
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")

# =========================
# Clean & detect columns
# =========================
df.columns = df.columns.str.strip()

padj_col = [c for c in df.columns if 'Adj. P-value' in c and 'SM_pos' in c][0]
l2fc_col = [c for c in df.columns if 'Log2 Fold Change' in c and 'SM_pos' in c][0]

df = df.rename(columns={
    padj_col: 'padj',
    l2fc_col: 'l2fc'
})

# =========================
# Preprocess
# =========================
df = df[df['Name'].notna()].copy()
df['logpadj'] = -np.log10(df['padj'])
df['lipid_class'] = df['Name'].str.split('(').str[0]

# =========================
# Enriched lipids only
# =========================
sig_df = df[
    (df['padj'] <= 0.05) &
    (df['l2fc'] >= 2)
].copy()

# =========================
# Color palette
# Define color palette (expanded)
class_colors = {
    'TG': '#00008B', 'PG': '#4B0082', 'PE': '#8A2BE2', 'PC': '#800080',
    'MGDG': '#FF00FF', 'MePC': '#FF69B4', 'LPC': '#FA8072', 'LBPA': '#FFA500',
    'DG': '#FFD700', 'AcHexCmE': '#FFFF00', 'ZyE': '#228B22', 'WE': '#00CED1',
    'SPH': '#7FFFD4', 'SM': '#20B2AA', 'MG': '#008000', 'Hex1Cer': '#ADFF2F',
    'ChE': '#B22222', 'Car': '#DC143C', 'CarE': '#8B0000', 'BisMeLPA': '#708090',
    'AcHexChE': '#C0C0C0', 'Co': '#FF4500', 'AcCa': '#6A5ACD', 'PI': '#32CD32',
    'Cer': '#8B4513', 'SiE': '#FF1493', 'SQDG': "#138B87", 'PS': '#00FF7F',
    'FA': '#FF6347', 'LdMePE': '#4682B4', 'ST': '#FFB6C1', 'LPE': '#FF1493',
    'PA': '#ADFF2F', 'PEt': '#FF4500', 'OAHFA': '#8A2BE2', 'dMePE': '#4682B4',
    'CmE': '#FF8C00'
}
# =========================

enriched_classes = sorted(sig_df['lipid_class'].unique())
class_colors_filtered = {k: v for k, v in class_colors.items()
                         if k in enriched_classes}

# =========================
# Plot
# =========================
fig, ax = plt.subplots(figsize=(10, 10))

sns.scatterplot(
    data=df,
    x='l2fc',
    y='logpadj',
    color='gainsboro',
    s=25,
    ax=ax
)

sns.scatterplot(
    data=sig_df,
    x='l2fc',
    y='logpadj',
    hue='lipid_class',
    palette=class_colors_filtered,
    s=70,
    edgecolor='white',
    linewidth=0.5,
    ax=ax
)

ax.axhline(-np.log10(0.05), linestyle='--', color='black', linewidth=0.8)
ax.axvline(2, linestyle='--', color='black', linewidth=0.8)

legend_elements = [
    lines.Line2D(
        [0], [0],
        marker='o',
        color='w',
        markerfacecolor=class_colors_filtered[lipid],
        markersize=9,
        label=lipid
    )
    for lipid in enriched_classes
]

ax.legend(
    handles=legend_elements,
    title='Enriched Lipid Classes\n(padj ≤ 0.05, L2FC ≥ 2)',
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    fontsize=9
)

ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-log10 (Adj. P-value)')

plt.tight_layout()
plt.savefig(
    'SM_vs_Water_positive_mode_volcano_enriched_only.svg',
    dpi=600,
    bbox_inches='tight'
)
plt.show()




''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Complete media vs water control (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (CM_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'
l2fc_col = 'Log2 Fold Change: (CM_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

#Assign colours
class_colors = {
'TG': '#00008B',       # dark blue
'PG': '#4B0082',       # indigo
'PE': '#8A2BE2',       # violet
'PC': '#800080',       # purple
'MGDG': '#FF00FF',     # magenta
'MePC': '#FF69B4',     # pink
'LPC': '#FA8072',      # salmon
'LBPA': '#FFA500',     # orange
'DG': '#FFD700',       # gold
'AcHexCmE': '#FFFF00', # yellow
'ZyE': '#228B22',      # forest green
'WE': '#00CED1',       # dark turquoise
'SPH': '#7FFFD4',      # aquamarine
'SM': '#20B2AA',       # light sea green
'MG': '#008000',       # green
'Hex1Cer': '#ADFF2F',  # green-yellow
'ChE': '#B22222',      # firebrick red
'Car': '#DC143C',      # crimson
'CarE': '#8B0000',     # dark red
'BisMeLPA': '#708090', # slate gray
'AcHexChE': '#C0C0C0', # silver
'Co': '#FF4500',       # orange red (added for missing key)
'AcCa': '#6A5ACD',     # slate blue (added for missing key)
'PI': '#32CD32',       # lime green (added for missing key)
'Cer': '#8B4513',      # saddle brown (added for missing key)
'SiE': '#FF1493',      # deep pink (added for missing key)
'SQDG': "#138B87",     # saddle brown (added for SQDG)
'PS': '#00FF7F',       # spring green (added for PS)
'FA': '#FF6347',       # tomato (added for FA)
'LdMePE': '#4682B4',   # steel blue (added for LdMePE)
'ST': '#FFB6C1',       # light pink (added for ST)
'LPE': '#FFD700',      # gold (added for LPE)
'PA': '#ADFF2F',       # green-yellow (added for PA)
'PEt': '#FF4500',      # orange red (added for PEt)
'OAHFA': "#695B76",    # violet (added for OAHFA)
'dMePE': '#4682B4',    # steel blue (added for dMePE)
'CmE': '#FF8C00' 
}
# Filter class_colors to include only enriched lipid classes
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(10,10))

# Background: all data in gray
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='lightgray', s=20)

# Highlight: significant lipid classes (positive log2 fold change and p-value < 0.05)
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.4)
# === Custom legend ===
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]
ax.legend(handles=handles, title='Lipid Class', frameon=False,
          bbox_to_anchor=(1.02, 1), loc='upper left')
# === Add significance threshold lines ===
ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(2, color='black', linestyle='dashed', linewidth=0.8)
plt.tight_layout()
plt.savefig('complete_media_versus_watercontrol(negative mode)12232025.svg',
            bbox_inches='tight', dpi=600)



''''''''''''''''''''''            '''''''''''''''''''
#rabbit serum vs water control (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (RS_NEG_NEWANALYSIS_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'
l2fc_col = 'Log2 Fold Change: (RS_NEG_NEWANALYSIS_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

#Assign colours
class_colors = {
'TG': '#00008B',       # dark blue
'PG': '#4B0082',       # indigo
'PE': '#8A2BE2',       # violet
'PC': '#800080',       # purple
'MGDG': '#FF00FF',     # magenta
'MePC': '#FF69B4',     # pink
'LPC': '#FA8072',      # salmon
'LBPA': '#FFA500',     # orange
'DG': '#FFD700',       # gold
'AcHexCmE': '#FFFF00', # yellow
'ZyE': '#228B22',      # forest green
'WE': '#00CED1',       # dark turquoise
'SPH': '#7FFFD4',      # aquamarine
'SM': '#20B2AA',       # light sea green
'MG': '#008000',       # green
'Hex1Cer': '#ADFF2F',  # green-yellow
'ChE': '#B22222',      # firebrick red
'Car': '#DC143C',      # crimson
'CarE': '#8B0000',     # dark red
'BisMeLPA': '#708090', # slate gray
'AcHexChE': '#C0C0C0', # silver
'Co': '#FF4500',       # orange red (added for missing key)
'AcCa': '#6A5ACD',     # slate blue (added for missing key)
'PI': '#32CD32',       # lime green (added for missing key)
'Cer': '#8B4513',      # saddle brown (added for missing key)
'SiE': '#FF1493',      # deep pink (added for missing key)
'SQDG': "#138B87",     # saddle brown (added for SQDG)
'PS': '#00FF7F',       # spring green (added for PS)
'FA': '#FF6347',       # tomato (added for FA)
'LdMePE': '#4682B4',   # steel blue (added for LdMePE)
'ST': '#FFB6C1',       # light pink (added for ST)
'LPE': '#FFD700',      # gold (added for LPE)
'PA': '#ADFF2F',       # green-yellow (added for PA)
'PEt': '#FF4500',      # orange red (added for PEt)
'OAHFA': "#695B76",    # violet (added for OAHFA)
'dMePE': '#4682B4',    # steel blue (added for dMePE)
'CmE': '#FF8C00' 
}
# Filter class_colors to include only enriched lipid classes
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(10, 10))

# Background: all data in gray
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='lightgray', s=20)

# Highlight: significant lipid classes (positive log2 fold change and p-value < 0.05)
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.4)

# === Labels and title ===
ax.set_xlabel('$log_2$(fold change)', fontsize=14)
ax.set_ylabel('$-log_{10}$(adjusted p-value)', fontsize=14)
ax.set_title('rabbit serum vs. water control (negative mode)', fontsize=10)

# === Add significance threshold lines ===
ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(2, color='black', linestyle='dashed', linewidth=0.8)

# === Custom legend ===
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]
ax.legend(handles=handles, title='Lipid Class', frameon=False,
          bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig('rabbitserum_versus_watercontrol(negative mode)12232025.svg',
            bbox_inches='tight', dpi=600)







            

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#spent media vs water control (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (SM_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'
l2fc_col = 'Log2 Fold Change: (SM_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

#Assign colours
class_colors = {
'TG': '#00008B',       # dark blue
'PG': '#4B0082',       # indigo
'PE': '#8A2BE2',       # violet
'PC': '#800080',       # purple
'MGDG': '#FF00FF',     # magenta
'MePC': '#FF69B4',     # pink
'LPC': '#FA8072',      # salmon
'LBPA': '#FFA500',     # orange
'DG': '#FFD700',       # gold
'AcHexCmE': '#FFFF00', # yellow
'ZyE': '#228B22',      # forest green
'WE': '#00CED1',       # dark turquoise
'SPH': '#7FFFD4',      # aquamarine
'SM': '#20B2AA',       # light sea green
'MG': '#008000',       # green
'Hex1Cer': '#ADFF2F',  # green-yellow
'ChE': '#B22222',      # firebrick red
'Car': '#DC143C',      # crimson
'CarE': '#8B0000',     # dark red
'BisMeLPA': '#708090', # slate gray
'AcHexChE': '#C0C0C0', # silver
'Co': '#FF4500',       # orange red (added for missing key)
'AcCa': '#6A5ACD',     # slate blue (added for missing key)
'PI': '#32CD32',       # lime green (added for missing key)
'Cer': '#8B4513',      # saddle brown (added for missing key)
'SiE': '#FF1493',      # deep pink (added for missing key)
'SQDG': "#138B87",     # saddle brown (added for SQDG)
'PS': '#00FF7F',       # spring green (added for PS)
'FA': '#FF6347',       # tomato (added for FA)
'LdMePE': '#4682B4',   # steel blue (added for LdMePE)
'ST': '#FFB6C1',       # light pink (added for ST)
'LPE': '#FFD700',      # gold (added for LPE)
'PA': '#ADFF2F',       # green-yellow (added for PA)
'PEt': '#FF4500',      # orange red (added for PEt)
'OAHFA': "#695B76",    # violet (added for OAHFA)
'dMePE': '#4682B4',    # steel blue (added for dMePE)
'CmE': '#FF8C00' 
}
# Filter class_colors to include only enriched lipid classes
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(10,10))

# Background: all data in gray
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='lightgray', s=20)

# Highlight: significant lipid classes (positive log2 fold change and p-value < 0.05)
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.4)


# === Custom legend ===
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]
ax.legend(handles=handles, title='Lipid Class', frameon=False,
          bbox_to_anchor=(1.02, 1), loc='upper left')

# === Add significance threshold lines ===
ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(2, color='black', linestyle='dashed', linewidth=0.8)
plt.tight_layout()
plt.savefig('SPENT_media_versus_watercontrol(negative mode)12232025.svg',
            bbox_inches='tight', dpi=600)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



#complete media vs rabbit serum
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (CM_pos_12022024_Newanalysis) / (Rabbitserum_12022024_Newanalysis)'
l2fc_col = 'Log2 Fold Change: (CM_pos_12022024_Newanalysis) / (Rabbitserum_12022024_Newanalysis)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

# === Assign colors ===
class_colors = {
    'TG': '#00008B',       # dark blue
    'PG': '#4B0082',       # indigo
    'PE': '#8A2BE2',       # violet
    'PC': '#800080',       # purple
    'MGDG': '#FF00FF',     # magenta
    'MePC': '#FF69B4',     # pink
    'LPC': '#FA8072',      # salmon
    'LBPA': '#FFA500',     # orange
    'DG': '#FFD700',       # gold
    'AcHexCmE': '#ADD8E6', # yellow
    'ZyE': '#228B22',      # forest green
    'WE': '#00CED1',       # dark turquoise
    'SPH': '#7FFFD4',      # aquamarine
    'SM': '#20B2AA',       # light sea green
    'MG': '#008000',       # green
    'Hex1Cer': '#ADFF2F',  # green-yellow
    'ChE': '#B22222',      # firebrick red
    'Car': '#DC143C',      # crimson
    'CarE': '#8B0000',     # dark red
    'BisMeLPA': '#708090', # slate gray
    'AcHexChE': '#C0C0C0', # silver
    'Co': '#FF4500',       # orange red (added for missing key)
    'AcCa': '#6A5ACD',     # slate blue (added for missing key)
    'PI': '#32CD32',       # lime green (added for missing key)
    'Cer': '#8B4513',      # saddle brown (added for missing key)
    'SiE': '#FF1493'       # deep pink (added for missing key)
}

# Filter class_colors to include only enriched lipid classes
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(6, 5))

# Background: all data in gray
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='lightgray', s=20)

# Highlight: significant lipid classes (positive log2 fold change and p-value < 0.05)
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.4)

# === Add significance threshold lines ===
ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(2, color='black', linestyle='dashed', linewidth=0.8)

# === Labels and title ===
ax.set_xlabel('$log_2$(fold change)', fontsize=10)
ax.set_ylabel('$-log_{10}$(adjusted p-value)', fontsize=10)
ax.set_title ('complete media vs rabbit serum (positive mode))', fontsize=10)

# === Custom legend ===
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]
ax.legend(handles=handles, title='Lipid Class', frameon=False,
          bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig('cm vs rs _positive_mode_colored_volcano_filtered_12232025.svg',
            bbox_inches='tight', dpi=600)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# complete media vs rabbit serum + excel file 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# =========================
# Load data
# =========================
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")
df.columns = df.columns.str.strip()

# =========================
# Define columns
# =========================
padj_col = 'Adj. P-value: (CM_pos_12022024_Newanalysis) / (Rabbitserum_12022024_Newanalysis)'
l2fc_col = 'Log2 Fold Change: (CM_pos_12022024_Newanalysis) / (Rabbitserum_12022024_Newanalysis)'

# =========================
# Preprocess
# =========================
df = df[df['Name'].notna()].copy()
df['l2fc'] = df[l2fc_col]
df['logpadj'] = -np.log10(df[padj_col])

df['lipid_class'] = df['Name'].str.split('(').str[0].str.strip()

# =========================
# CM-enriched lipids (RIGHT SIDE)
# =========================
cm_enriched = df[
    (df[padj_col] <= 0.05) &
    (df['l2fc'] >= 2)
].copy()

# =========================
# EXPORT TO EXCEL  ✅
# =========================
output_cols = [
    'Name',
    'lipid_class',
    'l2fc',
    padj_col
]

cm_enriched[output_cols].sort_values(
    by='l2fc', ascending=False
).to_excel(
    'CM_enriched_lipids_vs_RabbitSerum_posmode.xlsx',
    index=False
)

print(f"CM-enriched lipids exported: {cm_enriched.shape[0]}")

# =========================
# Color palette
# =========================
class_colors = {
    'TG': '#00008B',
    'PG': '#4B0082',
    'PE': '#8A2BE2',
    'PC': '#800080',
    'MGDG': '#FF00FF',
    'MePC': '#FF69B4',
    'LPC': '#FA8072',
    'LBPA': '#FFA500',
    'DG': '#FFD700',
    'AcHexCmE': '#ADD8E6',
    'ZyE': '#228B22',
    'WE': '#00CED1',
    'SPH': '#7FFFD4',
    'SM': '#20B2AA',
    'MG': '#008000',
    'Hex1Cer': '#ADFF2F',
    'ChE': '#B22222',
    'Car': '#DC143C',
    'CarE': '#8B0000',
    'BisMeLPA': '#708090',
    'AcHexChE': '#C0C0C0',
    'Co': '#FF4500',
    'AcCa': '#6A5ACD',
    'PI': '#32CD32',
    'Cer': '#8B4513',
    'SiE': '#FF1493'
}

# keep only enriched classes
enriched_classes = cm_enriched['lipid_class'].unique()
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

# =========================
# Volcano plot
# =========================
fig, ax = plt.subplots(figsize=(6, 5))

# background
sns.scatterplot(
    data=df, x='l2fc', y='logpadj',
    color='lightgray', s=20, ax=ax
)

# CM-enriched points
sns.scatterplot(
    data=cm_enriched,
    x='l2fc', y='logpadj',
    hue='lipid_class',
    palette=class_colors,
    s=60, edgecolor='white', linewidth=0.4,
    ax=ax
)

# thresholds
ax.axhline(-np.log10(0.05), linestyle='dashed', color='black', linewidth=0.8)
ax.axvline(2, linestyle='dashed', color='black', linewidth=0.8)

ax.set_xlabel('$log_2$(fold change)', fontsize=10)
ax.set_ylabel('$-log_{10}$(adjusted p-value)', fontsize=10)
ax.set_title('Complete media vs rabbit serum (positive mode)', fontsize=10)

# legend
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]

ax.legend(
    handles=handles,
    title='Lipid Class',
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc='upper left'
)

plt.tight_layout()
plt.savefig(
    'CM_vs_RS_positive_mode_volcano_CM_enriched.svg',
    dpi=600, bbox_inches='tight'
)
plt.show()
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''





#complte media vs rabbit serum (negative mode)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.lines as lines
import pandas as pd

# === Load data ===
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")

# === Define correct columns ===
padj_col = 'Adj. P-value: (CM_neg_newanalysis_12022024) / (RS_NEG_NEWANALYSIS_12022024)'
l2fc_col = 'Log2 Fold Change: (CM_neg_newanalysis_12022024) / (RS_NEG_NEWANALYSIS_12022024)'

# === Preprocess data ===
df = df.rename(columns={l2fc_col: 'l2fc'})
df['logpadj'] = -np.log10(df[padj_col])
df['significance'] = (df[padj_col] < 0.05) & (abs(df['l2fc']) >= 2.3)
df = df[df['Name'].notna()]

# === Identify significant subset and lipid classes ===
sig_df = df[df['significance']].copy()
sig_df['lipid_class'] = sig_df['Name'].apply(lambda x: x.split('(')[0].strip())

# === Filter for positive log2 fold change and p-value < 0.05 ===
sig_df = sig_df[(sig_df['l2fc'] > 0) & (sig_df[padj_col] < 0.05)]

# === Exclude lipid classes not enriched on the right side ===
# Get lipid classes with at least one significant point in the positive log2 fold change region
enriched_classes = sig_df['lipid_class'].unique()

# Filter sig_df to include only enriched lipid classes
sig_df = sig_df[sig_df['lipid_class'].isin(enriched_classes)]

#Assign colours
class_colors = {
'TG': '#00008B',       # dark blue
'PG': '#4B0082',       # indigo
'PE': '#8A2BE2',       # violet
'PC': '#800080',       # purple
'MGDG': '#FF00FF',     # magenta
'MePC': '#FF69B4',     # pink
'LPC': '#FA8072',      # salmon
'LBPA': '#FFA500',     # orange
'DG': '#FFD700',       # gold
'AcHexCmE': '#FFFF00', # yellow
'ZyE': '#228B22',      # forest green
'WE': '#00CED1',       # dark turquoise
'SPH': '#7FFFD4',      # aquamarine
'SM': '#20B2AA',       # light sea green
'MG': '#008000',       # green
'Hex1Cer': '#ADFF2F',  # green-yellow
'ChE': '#B22222',      # firebrick red
'Car': '#DC143C',      # crimson
'CarE': '#8B0000',     # dark red
'BisMeLPA': '#708090', # slate gray
'AcHexChE': '#C0C0C0', # silver
'Co': '#FF4500',       # orange red (added for missing key)
'AcCa': '#6A5ACD',     # slate blue (added for missing key)
'PI': '#32CD32',       # lime green (added for missing key)
'Cer': '#8B4513',      # saddle brown (added for missing key)
'SiE': '#FF1493',      # deep pink (added for missing key)
'SQDG': "#138B87",     # saddle brown (added for SQDG)
'PS': '#00FF7F',       # spring green (added for PS)
'FA': '#FF6347',       # tomato (added for FA)
'LdMePE': '#4682B4',   # steel blue (added for LdMePE)
'ST': '#FFB6C1',       # light pink (added for ST)
'LPE': '#FFD700',      # gold (added for LPE)
'PA': '#ADFF2F',       # green-yellow (added for PA)
'PEt': '#FF4500',      # orange red (added for PEt)
'OAHFA': '#8A2BE2',    # violet (added for OAHFA)
'dMePE': '#4682B4',    # steel blue (added for dMePE)
'CmE': '#FF8C00' 
}
# Filter class_colors to include only enriched lipid classes
class_colors = {k: v for k, v in class_colors.items() if k in enriched_classes}

sig_df['color'] = sig_df['lipid_class'].map(class_colors)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(6,5))

# Background: all data in gray
sns.scatterplot(data=df, x='l2fc', y='logpadj', ax=ax, color='lightgray', s=20)

# Highlight: significant lipid classes (positive log2 fold change and p-value < 0.05)
sns.scatterplot(data=sig_df, x='l2fc', y='logpadj', ax=ax, hue='lipid_class', palette=class_colors, s=60, edgecolor='white', linewidth=0.4)

# === Add significance threshold lines ===
ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(2, color='black', linestyle='dashed', linewidth=0.8)

# === Labels and title ===
ax.set_xlabel('$log_2$(fold change)', fontsize=14)
ax.set_ylabel('$-log_{10}$(adjusted p-value)', fontsize=14)
ax.set_title('complete media vs. rabbit serum (negative mode)', fontsize=10)

# === Custom legend ===
handles = [
    lines.Line2D([0], [0], marker='o', color='w',
                 markerfacecolor=color, markersize=8, label=lipid)
    for lipid, color in class_colors.items()
]
ax.legend(handles=handles, title='Lipid Class', frameon=False,
          bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig('complete media vs. rabbit serum (negative mode)12232025.svg',
            bbox_inches='tight', dpi=600)

