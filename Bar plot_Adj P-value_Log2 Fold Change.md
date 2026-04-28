Bar plot (Figs 1A and B: Adj. P-value/Log2 Fold Change)

import pandas as pd
import seaborn as sns 
from matplotlib import pyplot as plt 
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

#borrelia vs Internal control (negative mode)#
df = pd.read_excel("Dremoved_onesheet_Compounds_negmode_12022024.xlsx")
padj = 'Adj. P-value: (Bb_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)'
l2fc = "Log2 Fold Change: (Bb_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)"
df = df.rename(columns={'Log2 Fold Change: (Bb_neg_newanalysis_12022024) / (WATER_NEG_NEWANALYSIS_12022024)': 'l2fc'})
df[padj]
df.columns
Name = "Name"
df["logpadj"]= -np.log10(df[padj])
df["significance"]=(df[padj]<0.05)&(abs(df['l2fc'])>=2)
fig,_ax=plt.subplots(figsize=[5,4])
name_idx=df[Name].dropna().index
new_dfbbNEG=df.loc[name_idx]
temp_df = new_dfbbNEG.sort_values(by='logpadj',ascending=False).head(7)

columnstodrop=[column for column in new_dfbbNEG.columns if ":" in column] + ['Tags','Checked']

sig_df=new_dfbbNEG[(new_dfbbNEG[padj]<0.05)&(new_dfbbNEG['l2fc']>=2)].drop_duplicates(subset="Name").drop(columnstodrop, axis=1)
sig_df.Name.apply(lambda x:x.split("(")[0])
sig_df["lipid_class"]=sig_df.Name.apply(lambda x:x.split("(")[0])

names_to_drop = sig_df.Name.apply(lambda x: False if 'D' in x.split('(')[1] else True)
sig_df = sig_df[names_to_drop]
display(sig_df[['lipid_class','l2fc','logpadj']].groupby('lipid_class').mean().T)
to_plot = sig_df[['lipid_class','l2fc','logpadj']].groupby("lipid_class").mean().sort_values(by="logpadj",ascending=False).reset_index()
to_plot['logpadj']= to_plot['logpadj'].round(2)
_palette=sns.color_palette("flare",as_cmap=True,n_colors=len(to_plot.lipid_class.unique()))
fig,_ax = plt.subplots(figsize=[5,3],layout='constrained')
sns.barplot(data=to_plot,x='lipid_class',y='l2fc',hue='logpadj',palette=_palette,ax=_ax,errorbar='se')
_ax.tick_params(axis='both',labelsize=9)
_ax.tick_params(axis='x',rotation=90)
_ax.legend(frameon=False, fontsize=9, bbox_to_anchor=(1.0,0.9), title= "-log$_{10}$(adjpval)", title_fontsize=9)
plt.savefig('bb vs Internal control (neg mode_barplot_l2fc_pvalue_11082025).svg',dpi=600)






''''''''''''''''''''''''''''
#borrelia vs Internal control (positive mode)#
df = pd.read_excel("Compounds_posmode_12022024_D5TG_REM.xlsx")
padj = 'Adj. P-value: (Bb_pos_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)'
l2fc = "Log2 Fold Change: (Bb_pos_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)"
df = df.rename(columns={'Log2 Fold Change: (Bb_pos_12022024_Newanalysis) / (Water_pos_12022024_Newanalysis)': 'l2fc'})
df[padj]
df.columns
Name = "Name"
df["logpadj"]= -np.log10(df[padj])
df["significance"]=(df[padj]<0.05)&(abs(df['l2fc'])>=2)
fig,_ax=plt.subplots(figsize=[5,4])
name_idx=df[Name].dropna().index
new_dfbbNEG=df.loc[name_idx]
temp_df = new_dfbbNEG.sort_values(by='logpadj',ascending=False).head(7)

columnstodrop=[column for column in new_dfbbNEG.columns if ":" in column] + ['Tags','Checked']

sig_df=new_dfbbNEG[(new_dfbbNEG[padj]<0.05)&(new_dfbbNEG['l2fc']>=2)].drop_duplicates(subset="Name").drop(columnstodrop, axis=1)
sig_df.Name.apply(lambda x:x.split("(")[0])
sig_df["lipid_class"]=sig_df.Name.apply(lambda x:x.split("(")[0])

names_to_drop = sig_df.Name.apply(lambda x: False if 'D' in x.split('(')[1] else True)
sig_df = sig_df[names_to_drop]
display(sig_df[['lipid_class','l2fc','logpadj']].groupby('lipid_class').mean().T)
to_plot = sig_df[['lipid_class','l2fc','logpadj']].groupby("lipid_class").mean().sort_values(by="logpadj",ascending=False).reset_index()
to_plot['logpadj']= to_plot['logpadj'].round(2)
_palette=sns.color_palette("flare",as_cmap=True,n_colors=len(to_plot.lipid_class.unique()))
fig,_ax = plt.subplots(figsize=[5,3],layout='constrained')
sns.barplot(data=to_plot,x='lipid_class',y='l2fc',hue='logpadj',palette=_palette,ax=_ax,errorbar='se')
_ax.tick_params(axis='both',labelsize=9)
_ax.tick_params(axis='x',rotation=90)
_ax.legend(frameon=False, fontsize=9, bbox_to_anchor=(1.0,0.9), title= "-log$_{10}$(adjpval)", title_fontsize=9)
plt.savefig('bb vs Internal control (pos mode_barplot_l2fc_pvalue_11082025).svg',dpi=600)

