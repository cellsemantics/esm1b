#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:49:20 2023

@author: abhishekh
"""
import pandas as pd
import re
import numpy as np

def return_groupwise_amino_acid(filepath: str, sheet_name: str):
    
    
    """
    The function takes a file path and sheet name to extract score data, and outputs a DataFrame where 
    rows represent mutant codons, columns indicate mutation positions, and cell values indicate the 
    mean score for each mutation at the respective position and amino acid, after removing the stop codon.
    
    Args:
        filepath (str): Path to the metasheet.
        sheet_name (str): The name of the sheet in the metasheet where the score of interest is located.  
    
    Returns:
        pandas.DataFrame: A DataFrame where each row index indicates the mutant amino acid, each column 
        index indicates the position of mutation, and the value in each cell indicates the mean score 
        value for the mutation at a position with the mutant amino acid present at the DataFrame row index. 
    """
    
    df = pd.read_excel(filepath, sheet_name = sheet_name) # Read the score sheet in dataframe
    df = df[df.loc[:, "aa"]!="Stop"] # exclude the stop codon 
    original_order = df['aa'].unique()
    df['aa'] = pd.Categorical(df['aa'], categories=original_order, ordered=True)
    grouped_amino_acid= df.groupby('aa').mean().reset_index() 
    grouped_amino_acid.index = grouped_amino_acid["aa"]
    grouped_amino_acid.drop(["aa"], inplace=True, axis=1)
    grouped_amino_acid.columns = range(len(grouped_amino_acid.columns))
    return grouped_amino_acid

def pre_process_score(data:pd.DataFrame()):
    
    """
    Processes fitness scores data.
    
    Parameters:
        data (pandas DataFrame): Contains fitness scores with columns 'Amino acid', 'Codon', and 'aa'.
    
    Returns:
        df_tmp (pandas DataFrame): Contains all point mutations.
        complete_amino_acid (str): Represents the reference amino acid sequence after cleaning.
    """
    
    data = data.copy()
    data = data.iloc[:, 2:]  # drop column 'Amino acid', 'Codon'
    #### remove any stop codon from data
    data = data[data.loc[:, "aa"]!="Stop"]
    #### group the fitness w.r.t aa (amino acid) by mean value
    grouped_amino_acid = data.groupby('aa').mean().reset_index()
    #### get the list of amino acid in the ref sequence
    ref_seq_amino_acid = grouped_amino_acid.columns.to_list()
    ref_seq_amino_acid.remove("aa")


    lst = list()
    df_tmp = pd.DataFrame(columns=["mutant"])

    def remove_numbers_special_chars(s):
        return re.sub('[^A-Za-z]', '', s)  # Keep only alphabetic characters


    # Remove numbers and special characters from each string
    cleaned_list = [remove_numbers_special_chars(s) for s in ref_seq_amino_acid]
    
    complete_amino_acid = ''.join(cleaned_list) # this is the actual reference amino acid

    mutated_codon_list = list(np.unique(grouped_amino_acid["aa"])) # list of mutated codon


    for i in range(len(cleaned_list)):
        for j in range(len(mutated_codon_list)):

            lst.append(cleaned_list[i]+ str(i) + mutated_codon_list[j])
            
    df_tmp["mutant"] = lst
#     df_tmp["score"] = np.array(grouped_amino_acid.iloc[:, 1:]).T.ravel()
    return df_tmp, complete_amino_acid # #df_tmp : mutant dataframe and complete_amino_acid: actual referce amino acid sequence

def return_ztransformed_dataframe(data:pd.DataFrame()):

    """
    Z-Transforms a pandas dataframe for enhanced comparability and comprehension.
    
    Parameters:
        data (pandas dataframe): The dataframe to be transformed.

    Returns:
        z_transform (pandas dataframe): The transformed dataframe after applying the Z-transformation.
    """

    data = data.copy()
    z_transform = (data - np.nanmean(np.array(data)))/np.nanstd(np.array(data))
    return z_transform

def custom_figure_axis(ax, fontsize=10, show_ticks = True):

    """
    Customize the appearance of matplotlib axis for a figure.

    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis to be customized.
        fontsize (int, optional): Font size for axis labels and ticks. Default is 10.
        show_ticks (bool, optional): Whether to display ticks and labels. Default is True.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The customized axis.
    """

    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.2)
    ax.spines['left'].set_linewidth(0.2)
    ax.tick_params(axis='x', labelsize=fontsize, rotation=90)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.tick_params(axis='both', which='both', width=0.5)
    ax.xaxis.label.set_fontsize(fontsize)
    ax.yaxis.label.set_fontsize(fontsize)
    
    if show_ticks==False:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
    return ax

def man_whiteney(group1, group2):   

    """
    Perform a one-sided Mann-Whitney U test and return the p-value.

    This function compares two independent groups and tests if 'group1' tends to have larger values than 'group2'.

    Parameters:
        group1 (array-like): The data of the first group.
        group2 (array-like): The data of the second group.

    Returns:
        float: The p-value from the one-sided Mann-Whitney U test.

    Notes:
        The null hypothesis is that the distribution of 'group1' is not greater than 'group2'.
    """


    from scipy.stats import mannwhitneyu
    statistic, p_value = mannwhitneyu(group1, group2, alternative='greater')
    return p_value



def box_plot_esm_vs_all_three_score(data: pd.DataFrame(), ax):

    """
    Generate 3 custom box plots based on median cutoffs for fitness, functional score, and transformation score.

    Parameters:
        data (pd.DataFrame): Input dataframe containing various scores.
        ax (matplotlib.axes._subplots.AxesSubplot): The axis for plotting.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The customized axis after plotting.
    """

    import seaborn as sns
    import scipy.stats as stats

    data = data.copy()

    mean_funtional_cutoff = np.round(data.loc[:, "mean_funtional"].median(), 3)
    weighted_mean_fitness_cutoff =    np.round(data.loc[:, "weighted_mean_fitness"].median(), 3)
    transformation_score_cutoff =        np.round(data.loc[:, "transformation score"].median(), 3)

    lst_mean_funtional = list()
    lst_weighted_mean_fitness = list()
    lst_transformation_score = list()

    lst_weighted_mean_fitness.append(data[data.loc[:, "weighted_mean_fitness"]<=weighted_mean_fitness_cutoff]["esm"])
    lst_weighted_mean_fitness.append(data[data.loc[:, "weighted_mean_fitness"]>weighted_mean_fitness_cutoff]["esm"])

    lst_mean_funtional.append(data[data.loc[:, "mean_funtional"]<=mean_funtional_cutoff]["esm"])
    lst_mean_funtional.append(data[data.loc[:, "mean_funtional"]>mean_funtional_cutoff]["esm"])

    lst_transformation_score.append(data[data.loc[:,  'transformation score']<=transformation_score_cutoff]["esm"])
    lst_transformation_score.append(data[data.loc[:,  'transformation score']>transformation_score_cutoff]["esm"])

    sns.boxplot(lst_weighted_mean_fitness, ax=ax[0], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    sns.boxplot(lst_mean_funtional, ax=ax[1], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    sns.boxplot(lst_transformation_score, ax=ax[2], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    ax[0].set_ylabel("esm score", fontsize=5)
    ax[0].set_xlabel("weighted mean fitness")
    ax[1].set_xlabel("mean functional score")
    ax[2].set_xlabel("Transformation score")

    ax[0].set_xticklabels(["<=" + str(weighted_mean_fitness_cutoff), ">" + str(weighted_mean_fitness_cutoff)], fontsize=5)
    ax[1].set_xticklabels(["<=" + str(mean_funtional_cutoff), ">" + str(mean_funtional_cutoff)], fontsize=5)
    ax[2].set_xticklabels(["<=" + str(transformation_score_cutoff), ">" + str(transformation_score_cutoff)], fontsize=5)

    offset = 0.01

    formatted_p = "{:.2e}".format(man_whiteney(lst_weighted_mean_fitness[1], lst_weighted_mean_fitness[0]))
    ax[0].text((max(ax[0].get_xlim()) - offset), (max(ax[0].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')
    # print(formatted_p)
    formatted_p = "{:.2e}".format(man_whiteney(lst_mean_funtional[1], lst_mean_funtional[0]))
    ax[1].text((max(ax[1].get_xlim()) - offset), (max(ax[1].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')
    # print(formatted_p)

    formatted_p = "{:.2e}".format(man_whiteney(lst_transformation_score[1], lst_transformation_score[0]))
    ax[2].text((max(ax[2].get_xlim()) - offset), (max(ax[2].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')

    for i in range(3):
        ax[i] =  custom_figure_axis(ax[i], fontsize=5, show_ticks = True)

    return ax


# """Plot three subplots for the column name vs generation number for mutator, non mutator amd all"""
def return_mutator_non_mutator_column_name_wise_graph_together(data, mutator_list, non_mutator_list, column_name, fontsize=10):
    
    """
    Plot three subplots for the specified column against generation number for mutator, non-mutator, and the entire population.

    Parameters:
        data (pd.DataFrame): The input DataFrame containing the data.
        mutator_list (list): List of labels for mutator data points.
        non_mutator_list (list): List of labels for non-mutator data points.
        column_name (str): The column to plot against generation number.
        fontsize (int, optional): Font size for plot labels. Defaults to 10.

    Returns:
        ax (numpy.ndarray): Array of subplot axes.

    Notes:
        - Requires 'label' and 'generation_number' columns in the input DataFrame.

    """
    
    
    
    
    
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    import seaborn as sns
    data = data.copy()
    df_non_mutator_full_mutation_data = data[data['label'].isin(non_mutator_list)]
    df_mutator_full_mutation_data = data[data['label'].isin(mutator_list)]
    
    fig, ax = plt.subplots(3, 1 , dpi = 600, figsize = (2.1, 2.1), sharex=True);
    sns.lineplot(data = df_mutator_full_mutation_data, x = "generation_number", y=column_name, color="red", lw = 0.25, estimator ="median", ax = ax[0]);
    sns.lineplot(data = df_non_mutator_full_mutation_data, x = "generation_number", y=column_name, color="green", lw = 0.25, estimator ="median", ax = ax[1]);
    sns.lineplot(data = data, x = "generation_number", y=column_name, color="blue", lw = 0.25, estimator ="median", ax = ax[2]);



    ax[0].set_title(column_name + " vs generation for mutator", fontsize=fontsize);
    ax[1].set_title(column_name + " vs generation for non mutator", fontsize=fontsize);
    ax[2].set_title(column_name + " vs generation for all population", fontsize=fontsize);
    
    # Apply the custom tick formatter
    formatter = FuncFormatter(format_ticks)
    ax[2].xaxis.set_major_formatter(formatter)
    for i in range(3):
        ax[i]=custom_figure_axis(ax[i], fontsize=fontsize, show_ticks = True);
        ax[i].set_ylabel(column_name)  ;
        ax[i].legend(ncol=2, fontsize=3, frameon=False)
    plt.subplots_adjust(hspace=0.2)
    plt.subplots_adjust(wspace=0.2)  # You can adjust the value as needed

    plt.tight_layout();
    
    return ax

"""Convert xtick label to the 1K scale"""

def format_ticks(x:int, pos:int):

    """
    Format tick labels for a plot axis.

    Parameters:
        x (int): The tick value.
        pos (int): The position of the tick.

    Returns:
        str: The formatted tick label.

    Behavior:
        If x is greater than or equal to 1000, the function returns the value in 'k' format.
        Otherwise, the function returns the original value.

    """
    if x >= 1000:
        return f'{x//1000}k'
    else:
        return x


def return_combined_fitness_esm_data(fitness_dataframe, esm_dataframe):

    """
    Combines fitness data and ESM score data by calculating the median fitness and
    median ESM score for each generation. It then merges the dataframes based on the generation number.

    Parameters:
        fitness_dataframe (pd.DataFrame): DataFrame containing fitness data with columns 'Generation' and 'Fitness'.
        esm_dataframe (pd.DataFrame): DataFrame containing ESM score data with columns 'generation_number' and 'esm_score'.

    Returns:
        pd.DataFrame: Combined DataFrame with columns 'generation_number', 'median fitness', and 'median esm_score'.
    """

    fitness_median_generation_wise = fitness_dataframe.groupby(["Generation"])["Fitness"].agg(['median', 'std'])
    fitness_median_generation_wise = fitness_median_generation_wise.reset_index()
    fitness_median_generation_wise.columns = ["generation_number", "median fitness", "std fitness"]

    esm_median_generation_wise = esm_dataframe.groupby(["generation_number"])["esm_score"].agg(['median', 'std'])
    esm_median_generation_wise = esm_median_generation_wise.reset_index()
    esm_median_generation_wise.columns = ["generation_number", "median esm_score", "std esm_scor"]

    combined_fitness_esm = pd.merge(fitness_median_generation_wise, esm_median_generation_wise)

    return combined_fitness_esm


def return_generation_grouped_dataframe(data, gen_cut_off, cut_off_string1, cut_off_string2):

    """
    Group and aggregate data by generation and mutator status, and create a new column based on generation cut-offs.

    Parameters:
        data (pd.DataFrame): The input DataFrame containing the data.
        gen_cut_off (int): The generation number used as a cut-off.
        cut_off_string1 (str): Label for data points below or equal to the generation cut-off.
        cut_off_string2 (str): Label for data points above the generation cut-off.

    Returns:
        pd.DataFrame: A DataFrame with grouped and aggregated data, including the new grouping column.

    Notes:
        - The input DataFrame 'data' should have columns 'generation_number' and 'mutator'.
    """
    
    data = data.copy()
    
    esm_score_median = data.groupby(["generation_number", "mutator"])["esm_score"].median()
    esm_score_median_reset_index = esm_score_median.reset_index()
    esm_score_median_reset_index['group_gen'] = esm_score_median_reset_index['generation_number'].apply(lambda x: cut_off_string1 if x <= gen_cut_off else cut_off_string2)
    
    
    return esm_score_median_reset_index
    