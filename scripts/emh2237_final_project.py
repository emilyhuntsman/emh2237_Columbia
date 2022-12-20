# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.cluster import KMeans
#%% import data (S1) from Tian paper
tian_human = pd.read_excel('tian_sd01.xlsx',sheet_name='D. Human ECM signatures',skiprows=1)

# calculating means of each category from MS data
tian_human.loc[:,'Normal'] = tian_human[['N1', 'N2']].mean(axis=1)
tian_human.loc[:,'PanIN'] = tian_human[['P1', 'P2']].mean(axis=1)
tian_human.loc[:,'Pancreatitis'] = tian_human[['PT1', 'PT2']].mean(axis=1)
tian_human.loc[:,'PDAC'] = tian_human[[ 'WD1', 'WD2', 'WD3', 'MD1', 'MD2','MD3', 'MD4', 'PD1', 'PD2', 'PD3', 'PD4']].mean(axis=1)

# calculating means of each abundance category
tian_human.loc[:,'Normal Abundance'] = tian_human.loc[:,'N1.1']
tian_human.loc[:,'PanIN Abundance'] = tian_human[['P1.1', 'P2.1']].mean(axis=1)
tian_human.loc[:,'Pancreatitis Abundance'] = tian_human.loc[:,'PT1.1']
tian_human.loc[:,'PDAC Abundance'] = tian_human[[ 'WD1.1', 'MD1.1', 'MD2.1', 'PD1.1', 'PD2.1']].mean(axis=1)

# attempting to calculate logFC from raw abundance, did not match the log2(x/normal) conditions from paper
tian_human.loc[:,'Normal vs PanIN log'] = np.log2(tian_human.loc[:,'PanIN Abundance'] / tian_human.loc[:,'Normal Abundance'])
tian_human.loc[:,'Normal vs PDAC log'] = np.log2(tian_human.loc[:,'PDAC Abundance'] / tian_human.loc[:,'Normal Abundance'])
tian_human.loc[:,'Normal vs Pancreatitis log'] = np.log2(tian_human.loc[:,'Pancreatitis Abundance'] / tian_human.loc[:,'Normal Abundance'])

# calculate fold change, already in logspace so subtract
tian_human.loc[:,'Normal vs PanIN'] = tian_human.loc[:,'PanIN'] - tian_human.loc[:,'Normal']
tian_human.loc[:,'Normal vs PDAC'] = tian_human.loc[:,'PDAC']- tian_human.loc[:,'Normal']
tian_human.loc[:,'Normal vs Pancreatitis'] = tian_human.loc[:,'Pancreatitis'] - tian_human.loc[:,'Normal']

# %% PanIN vs Normal
for index,row in tian_human.iterrows():
    if tian_human.loc[index,'Overrepresented in PanIN (vs. normal corrected P<0.1, fold>1.5) Related to Figure 2A, D'] == 'yes':
        if  tian_human.loc[index,'Overrepresented in PDAC (vs. normal corrected P<0.1, fold>1.5) Related to Figure 2B, D'] == 'yes':
            val = 'Overrepresented in PanIN (and PDAC)'
        else:
            val = 'Overrepresented in PanIN only'
    else:
        val = 'Other'
    tian_human.loc[index,'key'] = val 

plt.figure(figsize=(4.25, 5))
color_dict = {'Other':'lightgrey','Overrepresented in PanIN (and PDAC)':'crimson','Overrepresented in PanIN only':'gold'}
groups = tian_human.groupby('key')
for name, group in groups:
    plt.plot(group['Normal vs PanIN'], group['Corrected t test P values (normal vs. PanIN) yellow: P<0.1'], marker='o', linestyle='', markersize=3, label=name, color=color_dict[name])

plt.xlabel('Normal vs PanIN', fontsize = 12)
plt.ylabel('Adjusted P-value', fontsize = 12)
plt.legend(fontsize = 8, loc= 'upper left', frameon=False)
plt.yscale(value='log')
ax = plt.gca()
ax.set_xlim([-3, 3])
ax.set_ylim([.01, 1])
ax.invert_yaxis()
plt.hlines(y=[.1],xmin=[-3],xmax=[3],linestyles='dotted')

# %% PDAC vs Normal
for index,row in tian_human.iterrows():
    if (row['Corrected t test P values (normal vs. PDAC) yellow: P<0.1']<.1):
        if tian_human.loc[index,'Overrepresented in PDAC (vs. normal corrected P<0.1, fold>1.5) Related to Figure 2B, D'] == 'yes':
            if  tian_human.loc[index,'Overrepresented in PanIN (vs. normal corrected P<0.1, fold>1.5) Related to Figure 2A, D'] == 'yes':
                val = 'Overrepresented in PDAC (and PanIN)'
            else:
                val = 'Overrepresented in PDAC only'
        else:
            val = 'Other'
    else:
        val = 'Other'
    tian_human.loc[index,'key'] = val 

plt.figure(figsize=(4.25, 5))
color_dict = {'Other':'lightgrey','Overrepresented in PDAC (and PanIN)':'crimson','Overrepresented in PDAC only':'gold'}
groups = tian_human.groupby('key')
for name, group in groups:
    plt.plot(group['Normal vs PDAC'], group['Corrected t test P values (normal vs. PDAC) yellow: P<0.1'], marker='o', linestyle='', markersize=3, label=name, color=color_dict[name])

plt.xlabel('Normal vs PDAC', fontsize = 12)
plt.ylabel('Adjusted P-value', fontsize = 12)
plt.legend(fontsize = 8, loc= 'upper left', frameon=False)
plt.yscale(value='log')
ax = plt.gca()
ax.set_xlim([-3.5, 3.5])
ax.set_ylim([.00001, 1])
ax.invert_yaxis()
plt.hlines(y=[.1],xmin=[-3.5],xmax=[3.5],linestyles='dotted')

# %% Pancreatitis vs Normal
plt.figure(figsize=(4.25, 5))
for index,row in tian_human.iterrows():
    if (row['Corrected t test p values (normal vs. pancreatitis) yellow: P<0.1']<.1):
        if tian_human.loc[index,'Overrepresented in pancreatits (vs. normal corrected P<0.1, fold>1.5) Related to Figure 2C, D'] == 'yes':
            val = 'Overrepresented in Pancreatitis'
        else:
            val = 'Other'
    else:
        val = 'Other'
    tian_human.loc[index,'key'] = val 

groups = tian_human.groupby('key')
color_dict = {'Other':'lightgrey','Overrepresented in Pancreatitis':'crimson'}
for name, group in groups:
    plt.plot(group['Normal vs Pancreatitis'], group['Corrected t test p values (normal vs. pancreatitis) yellow: P<0.1'], marker='o', linestyle='', markersize=3, label=name, color=color_dict[name])

plt.xlabel('Normal vs Pancreatitis', fontsize = 12)
plt.ylabel('Adjusted P-value', fontsize = 12)
plt.legend(fontsize = 8, loc= 'upper left', frameon=False)
ax = plt.gca()
ax.set_xlim([-3.5, 3.5])
ax.set_ylim([.01, 1])
ax.invert_yaxis()
plt.hlines(y=[.1],xmin=[-3.5],xmax=[3.5],linestyles='dotted')
plt.yscale(value='log')

# %% PDAC Normal Abundance
plt.figure() # smaller figure for higher abundances
tian_pdac_abun = tian_human[(tian_human['PDAC Abundance']>= 10000000000) & (tian_human['Normal Abundance']>= 10000000000)]

groups = tian_pdac_abun.groupby('Matrisome Category')
for name, group in groups:
    plt.plot(group['Normal Abundance'], group['PDAC Abundance'], marker='o', linestyle='', markersize=3, label=name)

for index,row in tian_pdac_abun.iterrows():
    plt.gca().text(row['Normal Abundance'], row['PDAC Abundance'], str(row['Entrez Gene Symbol']),fontsize='5')

plt.gca().plot([0,1],[0,1], transform=plt.gca().transAxes)
plt.yscale(value='log')
plt.xscale(value='log')
plt.xlabel('Protein abundance in normal pancreas', fontsize = 12)
plt.ylabel('Protein abundance in PDAC', fontsize = 12)
plt.title('Human')
plt.legend(fontsize='10')

# more populated figure for lower abundances
plt.figure()
tian_pdac_abun = tian_human[(tian_human['PDAC Abundance'] < 10000000000) & (tian_human['Normal Abundance'] < 10000000000)]

groups = tian_pdac_abun.groupby('Matrisome Category')
for name, group in groups:
    plt.plot(group['Normal Abundance'], group['PDAC Abundance'], marker='o', linestyle='', markersize=2, label=name)

for index,row in tian_pdac_abun.iterrows():
    if row['Entrez Gene Symbol'] in ['S100A6','S100A10','S100A4','S100A8','S100A9','S100A11']: # label select proteins from paper figure
        plt.gca().text(row['Normal Abundance'], row['PDAC Abundance'], str(row['Entrez Gene Symbol']),fontsize='5')

plt.gca().plot([0,1],[0,1], transform=plt.gca().transAxes)
plt.yscale(value='log')
plt.xscale(value='log')
plt.xlabel('Protein abundance in normal pancreas', fontsize = 12)
plt.ylabel('Protein abundance in PDAC', fontsize = 12)
plt.title('Human')
plt.legend(fontsize='10')
# %% 3D clustering attempt
tian_human_nonan = tian_human.dropna(subset=['Normal vs PanIN', 'Normal vs PDAC','Normal vs Pancreatitis'])
kmeans = KMeans(n_clusters=4)
y = kmeans.fit_predict(tian_human_nonan[['Normal vs PanIN', 'Normal vs PDAC','Normal vs Pancreatitis']])
tian_human_nonan['Cluster'] = y
print(tian_human_nonan)
# cluster on fold changes
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
groups = tian_human_nonan.groupby('Cluster')
for name, group in groups:
    ax.scatter(group['Normal vs PanIN'], group['Normal vs PDAC'],group['Normal vs Pancreatitis'], marker='o', label=name)
ax.set_xlabel('Normal vs PanIN')
ax.set_ylabel('Normal vs PDAC')
ax.set_zlabel('Normal vs Pancreatitis')
plt.show()
# matrisome category
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
groups = tian_human_nonan.groupby('Matrisome Category')
for name, group in groups:
    ax.scatter(group['Normal vs PanIN'], group['Normal vs PDAC'],group['Normal vs Pancreatitis'], marker='o', label=name)
ax.set_xlabel('Normal vs PanIN')
ax.set_ylabel('Normal vs PDAC')
ax.set_zlabel('Normal vs Pancreatitis')
# %% 2D clustering PanIN and PDAC
tian_human_nonan = tian_human.dropna(subset=['Normal vs PanIN', 'Normal vs PDAC'])
kmeans = KMeans(n_clusters=7)
y = kmeans.fit_predict(tian_human_nonan[['Normal vs PanIN', 'Normal vs PDAC']])
tian_human_nonan['Cluster'] = y
print(tian_human_nonan)
# clusters
fig = plt.figure()
ax = fig.add_subplot()
groups = tian_human_nonan.groupby('Cluster')
for name, group in groups:
    ax.scatter(group['Normal vs PanIN'], group['Normal vs PDAC'], marker='o', label=name)
ax.set_xlabel('Normal vs PanIN')
ax.set_ylabel('Normal vs PDAC')
plt.show()
# matrisome category
fig = plt.figure()
ax = fig.add_subplot()
groups = tian_human_nonan.groupby('Matrisome Category')
for name, group in groups:
    ax.scatter(group['Normal vs PanIN'], group['Normal vs PDAC'], marker='o', label=name)
ax.set_xlabel('Normal vs PanIN')
ax.set_ylabel('Normal vs PDAC')