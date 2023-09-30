 # Importing Libraries
import streamlit as st
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from IPython.display import set_matplotlib_formats
from decimal import Decimal
import re
import os
from lib.display import *
from lib.tools import *

def extract_compounds(uploaded_file, selected_sheet):
    """
    Extract compounds from the specified Excel file and sheet.
    
    Parameters:
    - uploaded_file: The uploaded Excel file.
    - selected_sheet: The sheet name to extract data from.
    
    Returns:
    - A flattened numpy array of the compounds.
    """
    # Read specific range for compounds from the selected sheet
    compounds_D = pd.read_excel(uploaded_file, usecols="D", skiprows=5, nrows=16, header=None, sheet_name=selected_sheet).fillna("NONE")
    compounds_O = pd.read_excel(uploaded_file, usecols="O", skiprows=5, nrows=16, header=None, sheet_name=selected_sheet).fillna("NONE")

    # Combine the dataframes
    combined_compounds = pd.concat([compounds_D, compounds_O], axis=1)
    
    # Convert the combined dataframe to a numpy array and flatten
    compounds_np_array = combined_compounds.to_numpy()
    flattened_array = np.concatenate([compounds_np_array[:, 0], compounds_np_array[:, 1]])
    
    return combined_compounds, flattened_array


def extract_head_data(uploaded_file, selected_sheet):
    """
    Extract assay and title texts from the specified Excel file and sheet.
    
    Parameters:
    - uploaded_file: The uploaded Excel file.
    - selected_sheet: The sheet name to extract data from.
    
    Returns:
    - Tuple of (assay_text, title_text).
    """
    assay = pd.read_excel(uploaded_file, usecols="A", skiprows=0, nrows=1, header=None, sheet_name=selected_sheet).fillna("NONE")
    title = pd.read_excel(uploaded_file, usecols="A", skiprows=1, nrows=1, header=None, sheet_name=selected_sheet).fillna("NONE")
    
    assay_text = assay.iloc[0, 0]
    title_text = title.iloc[0, 0]
    
    return assay_text, title_text


def extract_concentrations(uploaded_file, selected_sheet):
    """
    Extract data from columns C through M and N through X, for rows 25 to 40, from the specified Excel file and sheet.
    
    Parameters:
    - uploaded_file: The uploaded Excel file.
    - selected_sheet: The sheet name to extract data from.
    
    Returns:
    - A numpy array containing the extracted data.
    """
    # Define the columns to extract
    columns_CM = [chr(i) for i in range(ord('C'), ord('M')+1)]
    columns_NX = [chr(i) for i in range(ord('N'), ord('X')+1)]

    # Extract data for CM columns
    data_list_CM = []
    for col in columns_CM:
        try:
            data = pd.read_excel(uploaded_file, usecols=col, skiprows=24, nrows=16, header=None, sheet_name=selected_sheet)
            data_list_CM.append(data.squeeze())  # Use squeeze to convert single columns to Series
        except ValueError:
            continue
    combined_data_CM = pd.concat(data_list_CM, axis=1)
    
    # Extract data for NX columns
    data_list_NX = []
    for col in columns_NX:
        try:
            data = pd.read_excel(uploaded_file, usecols=col, skiprows=24, nrows=16, header=None, sheet_name=selected_sheet)
            data_list_NX.append(data.squeeze())  # Use squeeze to convert single columns to Series
        except ValueError:
            continue
    combined_data_NX = pd.concat(data_list_NX, axis=1)

    # Combine the data
    all_combined = np.c_[combined_data_CM.T.to_numpy(), combined_data_NX.T.to_numpy()]
    
    return all_combined


def extract_experiment(uploaded_file, selected_sheet, start_row):
    """
    Extract data from columns C through M and N through X, for rows 25 to 40, from the specified Excel file and sheet.
    
    Parameters:
    - uploaded_file: The uploaded Excel file.
    - selected_sheet: The sheet name to extract data from.
    
    Returns:
    - A numpy array containing the extracted data.
    """
    # Define the columns to extract
    columns_CM = [chr(i) for i in range(ord('C'), ord('M')+1)]
    columns_NX = [chr(i) for i in range(ord('N'), ord('X')+1)]

    # Extract data for CM columns
    data_list_CM = []
    for col in columns_CM:
        try:
            data = pd.read_excel(uploaded_file, usecols=col, skiprows=start_row, nrows=16, header=None, sheet_name=selected_sheet)
            data_list_CM.append(data.squeeze())  # Use squeeze to convert single columns to Series
        except ValueError:
            continue
    combined_data_CM = pd.concat(data_list_CM, axis=1)
    
    # Extract data for NX columns
    data_list_NX = []
    for col in columns_NX:
        try:
            data = pd.read_excel(uploaded_file, usecols=col, skiprows=start_row, nrows=16, header=None, sheet_name=selected_sheet)
            data_list_NX.append(data.squeeze())  # Use squeeze to convert single columns to Series
        except ValueError:
            continue
    combined_data_NX = pd.concat(data_list_NX, axis=1)

    # Combine the data
    all_combined = np.c_[combined_data_CM.T.to_numpy(), combined_data_NX.T.to_numpy()]
    
    return all_combined


# Get y-label:
def extract_ylabel(uploaded_file, selected_sheet):
    # Read the Excel file into a DataFrame
    df = pd.read_excel(uploaded_file, sheet_name=selected_sheet,  header=None)
    
    # Extract the value from cell C42
    y_label = df.at[41, 2] # because pandas uses 0-based indexing
    
    return y_label
