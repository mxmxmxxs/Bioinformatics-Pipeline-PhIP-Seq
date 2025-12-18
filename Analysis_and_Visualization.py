# %%
#Import txt files and transpose them

import os
import pandas as pd #needed for reading csv
 
#Get all txt files in subfolder

subfolder = '/Users/mima/Desktop/Masterthesis/000_SpecialCourse_KU/Code/data/' #specify path to folder with txt files

txt_files = [f for f in os.listdir(subfolder) if f.endswith('.txt')] #list of all txt files in folder

#Iterate and transpose each file

all_dfs = []

for file in txt_files: #loop through all txt files in the folder 

    filepath = os.path.join(subfolder, file)

    df = pd.read_csv(filepath, sep=";", header=None, usecols=[2]) 

    df_transposed = df.T

    df_transposed.index = [file] * len(df_transposed)

    all_dfs.append(df_transposed) #store each transposed dataframe in a list
 
#Concatenate all dataframes

result = pd.concat(all_dfs, axis=0)

print(result)

result.to_csv("transposed_output.csv", sep=";")  # Saves to combined CSV file in working directory
 

# %%
#not needed but checks if csv was saved correctly

transposed_df=pd.read_csv("transposed_output.csv", sep=";")
transposed_df


# %%
#!pip3 install openpyxl #uncomment if openpyxl is not installed
import openpyxl #needed for reading excel files

#Import and read the VirScan_Info.xlsx file
virscan_info = pd.read_excel("/Users/mima/Desktop/Masterthesis/000_SpecialCourse_KU/Code/data/VirScan_Info.xlsx")
print(virscan_info)

# %% [markdown]
# 

# %%
#Sort by Species column alphabetically which is needed to group by species later
virscan_info_sorted = virscan_info.sort_values(by='Species', ascending=True).reset_index(drop=True)

#Save sorted dataframe as a new csv file
virscan_info_sorted.to_csv("VirScan_Info_sorted_by_species.csv", index=False)

print("Sorted VirScan_Info saved as 'VirScan_Info_sorted_by_species.csv'")
print(virscan_info_sorted)

# %%
#Number and list of unique proteins in order to see how many different proteins are in the VirScan_Info file
unique_proteins_list = virscan_info['Protein names'].unique().tolist()
print(f"Number of unique proteins: {len(unique_proteins_list)}")
print(f"Unique protein names:\n{unique_proteins_list}")

# %%
#Number and list of unique viruses
unique_viruses_list = virscan_info_sorted['Species'].unique().tolist()
print(f"Number of unique viruses: {len(unique_viruses_list)}")
print(f"Unique virus names:\n{unique_viruses_list}")

# %%
#Count rows for each unique protein name (inorder of appearance) along with Species, Entry, Gene names which is needed to sum up proteins later
protein_counts = virscan_info.groupby(['Species','Entry','Gene names','Protein names'], sort=False).size().reset_index(name='Count')
print(protein_counts)

protein_counts.to_csv("protein_counts.csv", sep=";")  # Save to CSV filex

# %%
#Count rows for each unique virus name using the sorted dataframe
virus_counts = virscan_info_sorted.groupby(['Species'], sort=False).size().reset_index(name='Count')
print(virus_counts)

virus_counts.to_csv("virus_counts.csv", sep=";")  # Save to CSV file

# %%
#Counts based on protein name and counted rows, saves results as new csv file

transposed_df = pd.read_csv("transposed_output.csv", sep=";", index_col=0)
protein_counts = pd.read_csv("protein_counts.csv", sep=";")

#Clean up index names (remove "_idxstats.txt")
transposed_df.index = transposed_df.index.str.replace("_idxstats.txt", "")

#Replace index with numbers (for each of the 16 samples)
transposed_df.index = [f"Sample_{i}" for i in range(1, len(transposed_df) + 1)]

#Get the count for each unique protein combination
protein_groups = protein_counts.groupby(['Protein names', 'Species', 'Gene names']).size()

#Create a new dataframe to store summed up results
result_df_protein = pd.DataFrame(index=transposed_df.index)

#Sum up columns for each unique protein
col_index = 0
for idx, (protein_info, count) in enumerate(protein_counts.iterrows()):
    protein_name = protein_counts.loc[idx, 'Protein names']
    species = protein_counts.loc[idx, 'Species']
    gene_name = protein_counts.loc[idx, 'Gene names']

    #Sum up the next 'Count' columns
    cols_to_sum = transposed_df.iloc[:, col_index:col_index + protein_counts.loc[idx, 'Count']]
    result_df_protein[f"{protein_name}_{species}_{gene_name}"] = cols_to_sum.sum(axis=1)

    col_index += protein_counts.loc[idx, 'Count']

#Save results
result_df_protein.to_csv("protein_summed_output.csv", sep=";")
print(result_df_protein)

result_df_protein

# %%
#Counts based on virus name and counted rows, saves results as new csv file

transposed_df = pd.read_csv("transposed_output.csv", sep=";", index_col=0)
virus_counts = pd.read_csv("virus_counts.csv", sep=";")

#Clean up index names (remove "_idxstats.txt")
transposed_df.index = transposed_df.index.str.replace("_idxstats.txt", "")

#Replace index with numbers (for each of the 16 samples)
transposed_df.index = [f"Sample_{i}" for i in range(1, len(transposed_df) + 1)]

#Get the count for each unique virus combination
virus_groups = virus_counts.groupby(['Species']).size()

#Create a new dataframe to store summed up results
result_df_virus = pd.DataFrame(index=transposed_df.index)

#Sum up columns for each unique virus
col_index = 0
for idx, (virus_info, count) in enumerate(virus_counts.iterrows()):
    species = virus_counts.loc[idx, 'Species']

    #Sum up the next 'Count' columns
    cols_to_sum = transposed_df.iloc[:, col_index:col_index + virus_counts.loc[idx, 'Count']]
    result_df_virus[f"{species}"] = cols_to_sum.sum(axis=1)

    col_index += virus_counts.loc[idx, 'Count']

#Save results
result_df_virus.to_csv("virus_summed_output.csv", sep=";")
print(result_df_virus)

result_df_virus

# %%
#Create heatmap for top 20 proteins and saves it as png file

#!pip3 install matplotlib seaborn #remove hashtag if not installed
import matplotlib.pyplot as plt #needed for plotting
import seaborn as sns #needed for heatmap
from matplotlib.colors import LinearSegmentedColormap #needed for custom colormap

#Read protein results
result_df_protein = pd.read_csv("protein_summed_output.csv", sep=";", index_col=0)

#Transpose the dataframe to switch axes
result_df_transposed_protein = result_df_protein.T
result_df_transposed_protein.sort_values(by=result_df_transposed_protein.columns.tolist(), axis=0, ascending=False, inplace=True) #sort rows by total counts
results_top20_protein = result_df_transposed_protein.head(20)

#Create custom colormap with white for zero values
colors = ['white', '#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026']
n_bins = 100
cmap = LinearSegmentedColormap.from_list('YlOrRd_white', colors, N=n_bins)

#Shorten protein names on y-axis (maximum 50 characters) because they were wayyyy too long
results_top20_protein.index = [name[:50] + '...' if len(name) > 50 else name for name in results_top20_protein.index]

#Create heatmap
plt.figure(figsize=(14, 8))
sns.heatmap(results_top20_protein, cmap=cmap, cbar_kws={'label': 'Count'}, vmin=0)
plt.title('Protein Expression Heatmap - Top 20', fontsize=14, fontweight='bold')
plt.xlabel('Samples', fontsize=12, fontweight='bold')
plt.ylabel('Proteins', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=9, rotation=0)
plt.tight_layout()
plt.savefig('protein_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

# %%
#Create heatmap for top 20 viruses and saves it as png file

#Read virus results
result_df_virus = pd.read_csv("virus_summed_output.csv", sep=";", index_col=0)

#Transpose the dataframe to switch axes
result_df_transposed_virus = result_df_virus.T
result_df_transposed_virus.sort_values(by=result_df_transposed_virus.columns.tolist(), axis=0, ascending=False, inplace=True) #sort rows by total counts
results_virus_top20 = result_df_transposed_virus.head(20)

#Create custom colormap with white for zero values
colors = ['white', '#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026']
n_bins = 100
cmap = LinearSegmentedColormap.from_list('YlOrRd_white', colors, N=n_bins)

#Create heatmap
plt.figure(figsize=(14, 8))
sns.heatmap(results_virus_top20, cmap=cmap, cbar_kws={'label': 'Count'}, vmin=0)
plt.title('Virus Expression Heatmap - Top 20', fontsize=14, fontweight='bold')
plt.xlabel('Samples', fontsize=12, fontweight='bold')
plt.ylabel('Viruses', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=9, rotation=0)
plt.tight_layout()
plt.savefig('virus_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()


