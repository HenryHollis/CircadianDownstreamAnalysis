#!/usr/bin/env python
# coding: utf-8



import json
import requests
import sys
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap

# using argparse to handle command line arguments, demands outputFilename positional argument
parser = argparse.ArgumentParser(description="Python Interface to EnrichR",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-l', "--lib", nargs='+', help="EnrichR library (https://maayanlab.cloud/Enrichr/#libraries)", default = ["KEGG_2021_Human", \
"GO_Biological_Process_2023", "MSigDB_Hallmark_2020", "WikiPathway_2021_Human", "TRRUST_Transcription_Factors_2019", "GO_Cellular_Component_2023",\
"TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "TargetScan_microRNA_2017", "miRTarBase_2017"])
# possible other libs: "GO_Biological_Process_2023", "MSigDB_Hallmark_2020", "WikiPathway_2021_Human"
parser.add_argument('-g',"--genes_list_file", help="Path to gene list")
parser.add_argument('-b',"--background_list_file", help="Path to Background")
args = parser.parse_args()
config = vars(args)
print("Arguments: ")
print(config)

# list of enrichR libraries, by default, Kegg, GOBP, Hallmarks and wikipathways
libs = args.lib

# #Open file selector
# def open_file_dialog():
#     Tk().withdraw()  # Hide the root window
#     filename = askopenfilename()  # Open the file select dialog
#     return filename

# Get genes from file selected by user
# print("Please select file containing list of genes (Symbols)")
gene_list_path =  args.genes_list_file   #open_file_dialog()
# print("Path to genes list:", gene_list_path)
parent_folder_path = os.path.dirname(os.path.dirname(gene_list_path))
file_name = os.path.splitext(os.path.basename(gene_list_path))[0]
file_name = str.replace(file_name, "_", "")

print()

#Get Background from file selected by user
# print("Please select file containing list of BACKGROUND genes (Symbols)")
background_path =  args.background_list_file #open_file_dialog()
# print("Path to Background list:", background_path)
background_name = os.path.splitext(os.path.basename(background_path))[0]
background_name = str.replace(background_name, "_", "")

def parse_csv_file(file_path):
    try:
        dataframe = pd.read_csv(file_path)  # Read the CSV file into a DataFrame
        return dataframe
    except FileNotFoundError:
        print("File not found.")
        return None
    except pd.errors.EmptyDataError:
        print("Empty file.")
        return None
    except pd.errors.ParserError:
        print("Error parsing the file.")
        return None

#call function to parse csv file containing gene list
gene_list_df = parse_csv_file(gene_list_path)
# if gene_list_df is not None:
#     print("Background DataFrame (nrow: {}):".format(gene_list_df.shape[0]))
#     print(gene_list_df.head())
    
#call function to parse csv file containing background list

background_df = parse_csv_file(background_path)
# if background_df is not None:
#     print("Background DataFrame (nrow: {}):".format(background_df.shape[0]))
#     print(background_df.head(3))



# cast these into lists instead of pandas objects
genes =list(gene_list_df['Gene_Symbols'])
background = list(background_df['Gene_Symbols'])

####### get user and background response info ######
base_url = "https://maayanlab.cloud/speedrichr"
description = "sample gene set with background"
res = requests.post(
    base_url+'/api/addList',
    files=dict(
      list=(None, '\n'.join(genes)),
      description=(None, description),
    )
  )
if res.ok:
	userlist_response = res.json()
else:
    print("ERROR with requests.post")

res = requests.post(
        base_url+'/api/addbackground',
        data=dict(background='\n'.join(background)),
  )

if res.ok:
        background_response = res.json()
################################################
def plot_results(out, filename):
    fig, ax = plt.subplots(figsize = (0.1, 2.75))
    cmap = mpl.cm.bwr_r
    norm = mpl.colors.Normalize(vmin = 0, vmax = .2)

    mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)

    cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
    cbl = cbl.ax.invert_xaxis()

    plt.figure(figsize = (8,6))


    ax = sns.barplot(data = out.head(20), x = 'gene count', y = 'Term name', palette = mapper.to_rgba(out.head(20)["Adjusted p-value"]))

    ax.set_yticklabels([textwrap.fill(e, 100) for e in out.head(20)['Term name']])
    ax.figure.colorbar(mapper)
    plt.title("EnrichR Results: {}".format(filename))
    plt.subplots_adjust(left=0.4, right=0.99, top=0.9, bottom=0.1)

    plt.savefig((filename + "_barplot.png") )


#Run through libraries and download enrichment results to outfiles
dir_path = os.path.join(parent_folder_path, "EnrichR_results", ("EnrichR_" + file_name + "_BCKGRND-" + background_name ))
if not os.path.exists(dir_path):
    # If it doesn't exist, create the directory
    os.makedirs(dir_path)
    print(f'Directory "{dir_path}" created.')
    os.chdir(dir_path)
    out_data = []

    try:
        for lib in libs:
            res = requests.post(
                    base_url+'/api/backgroundenrich',
                    data=dict(
                    userListId=userlist_response["userListId"],
                    backgroundid=background_response["backgroundid"],
                    backgroundType=lib
                    )
                )
            if res.ok:
                results = res.json()
            else:
                print("Error with {}".format(lib))

            out = pd.DataFrame(results[lib], columns = ["Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"])
            out['gene count'] = out['Overlapping genes'].apply(lambda x: len(x))
            output_path = (file_name + "_LIBRY-"+ str.replace(lib, "_", "") + "_BCKGRND-" + background_name + ".csv" )
            plot_results(out,  (file_name + "_LIBRY-"+ str.replace(lib, "_", "")  ) )
            out_data.append(out)
            out.to_csv(output_path, index = False)
            print("Saved {} results to file".format(lib))
    except:
        pass
    cmap = mpl.cm.bwr_r
    norm = mpl.colors.Normalize(vmin = 0, vmax = .2)

    mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)

    # Create a new figure and axes for the combined plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()

    if len(out_data) > 4:
        out_data = out_data[0:4]

    for i, out in enumerate(out_data):
        ax = sns.barplot(data = out.head(20), x = 'gene count', y = 'Term name', palette = mapper.to_rgba(out.head(20)["Adjusted p-value"]), ax=axes[i])
        ax.set_yticklabels([textwrap.shorten(e, 35, placeholder="...") for e in out.head(20)['Term name']])
        axes[i].set_title(libs[i])
        
    axes[0].figure.colorbar(mapper)
    plt.tight_layout(h_pad= 1, w_pad = .3)
    plt.subplots_adjust(left=.2, right=0.98, top=0.9, bottom=0.1)
    plt.savefig((file_name + "_barplot_AllLibs.png") )
else:
    print(f'Directory "{dir_path}" already exists. Skipping.')
