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
from lib.extract import *


def highlight_values(series):
    """Function to highlight specific values."""
    colors_list = []
    indices = []
    color_mapping = {
        'none': 'color: red',
        'mock': 'color: lightblue',
        'mek': 'color: orange'
    }
    
    for idx, val in enumerate(series):
        val_lower = val.lower()
        if val_lower in color_mapping:
            colors_list.append(color_mapping[val_lower])
            indices.append(idx)
        else:
            colors_list.append('color: black')
            
    return colors_list, indices


def reshape_dataframe(df):
    """Reshape a long dataframe into two columns for better display."""
    half_length = len(df) // 2
    first_half = df.iloc[:half_length].reset_index(drop=True)
    second_half = df.iloc[half_length:].reset_index(drop=True)
    reshaped_data = pd.concat([first_half, second_half], axis=1)
    reshaped_data.columns = ["Compounds 1-16", "Compounds 17-32"]
    return reshaped_data



def style_table(styler, css_color, index_column):
    """
    Apply a specific color to a specific column in a DataFrame.

    Parameters:
    - styler (pd.io.formats.style.Styler): The DataFrame Styler object to update.
    - css_color (str): The CSS color property like 'color: red'.
    - index_column (int): The index of the column to apply the color to (starting from 0).

    Returns:
    - updated Styler object
    """
    # Extract the color value from the CSS property
    color = css_color.split(':')[-1].strip()

    col_name = styler.data.columns[index_column]
    return styler.applymap(lambda x: f'background-color: {color}', subset=pd.IndexSlice[:, col_name])



# Define the function to be fitted
def variable_slope_log_inhibitor_response(x, a, b, c, d):
    y = a + ((b - a) / (1 + np.power(10, (c - x) * d)))
    return y


# Get initial IC50 for the fit
def compute_x_at_ymid(x, y):
    # Compute Ymid
    y_mid = np.mean([np.max(y), np.min(y)])
    
    # Interpolate the data
    f = interp1d(y, x)
    
    # Find the X value at Ymid
    x_mid = f(y_mid)
    
    return x_mid


def filter_compounds(arr):
    arr = np.array(arr,  dtype=str)
    arr_lower = np.char.lower(arr)
    
    # Find indices of 'MEK' and 'MOCK'
    mek_indices  = np.where(arr_lower == 'mek')[0]
    mock_indices = np.where(arr_lower == 'mock')[0]
    none_indices = np.where(arr_lower == 'none')[0]


    
    # Create mask for elements to keep
    mask = ~np.isin(arr_lower, ['mek', 'mock', 'none'])
    
    # Apply mask to get filtered array
    filtered_arr = arr[mask]


    return mek_indices, mock_indices, none_indices, filtered_arr

