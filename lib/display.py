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
from lib.extract import *
from lib.tools import *

def display_compounds(df):
    """
    Reshape, style, and display the compounds dataframe in Streamlit.
    
    Parameters:
    - df: The combined compounds dataframe.
    
    Returns:
    - Tuple of (colored_indices, combined_colors).
    """
    # Melt the dataframe
    melted_df = df.melt(value_name="Compounds").drop(columns="variable")

    # Reset index to start from 1
    melted_df.index = melted_df.index + 1

    # Reshape the dataframe for display
    reshaped_df = reshape_dataframe(melted_df)

    # Apply styling and collect indices
    colors_col1, indices_col1 = highlight_values(reshaped_df["Compounds 1-16"])
    colors_col2, indices_col2 = highlight_values(reshaped_df["Compounds 17-32"])
    
    colored_indices = indices_col1 + [i + 16 for i in indices_col2]
    combined_colors = colors_col1 + colors_col2
    print(combined_colors)

    styled_df = reshaped_df.style.apply(lambda x: colors_col1, axis=0, subset=["Compounds 1-16"]) \
                                .apply(lambda x: colors_col2, axis=0, subset=["Compounds 17-32"])
    
    # Display the styled dataframe in Streamlit
    st.markdown("## Extracted Compounds")
    st.markdown(styled_df.to_html(escape=False), unsafe_allow_html=True)

    # To check the indices
    return colored_indices, combined_colors



def display_head_data(assay_text, title_text):
    """
    Display the assay and title texts in Streamlit with specified formatting.
    
    Parameters:
    - assay_text: The extracted assay text.
    - title_text: The extracted title text.
    
    Returns:
    - None (directly displays the texts in Streamlit).
    """
    st.markdown(f"# **{assay_text}**")
    st.markdown(f'<p style="font-size:24px">{title_text}</p>', unsafe_allow_html=True)



def display_concentrations(data, colored_indices, combined_colors):
    """
    Display the extracted data in a table format in Streamlit with colored columns.
    
    Parameters:
    - data: The numpy array containing the extracted data.
    - colored_indices: List of column indices to be colored.
    - combined_colors: List of colors corresponding to each column.
    
    Returns:
    - None (directly displays the data in Streamlit).
    """
    df = pd.DataFrame(data)  # Convert the numpy array to a dataframe
    styler = df.style

    for idx, color in zip(colored_indices, combined_colors):
        styler = style_table(styler, combined_colors[idx], idx)
    
    st.markdown("## Concentration uM")
    st.markdown(styler.to_html(escape=False), unsafe_allow_html=True)




def display_experiment(df, number):
    """
    Display the extracted data in a table format in Streamlit.
    
    Parameters:
    - df: The dataframe containing the extracted data.
    
    Returns:
    - None (directly displays the data in Streamlit).
    """
    st.markdown("## Experiment "+str(number))
    st.table(df)

