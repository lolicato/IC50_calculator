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
import subprocess
import re
import os
from lib.display import *
from lib.extract import *
from lib.tools import *


# Set the page configuration
st.set_page_config(
    page_title="Compounds Extraction",
    layout="centered",
    initial_sidebar_state="collapsed",
    page_icon="🧪"
)

plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.figsize'] = 12, 8
plt.rcParams['axes.labelsize'] = 50
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['font.size'] = 30
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['xtick.major.size'] = 15 # Set the length of x-axis ticks to 5
plt.rcParams['ytick.major.size'] = 15 # Set the length of y-axis ticks to 5
plt.rcParams["xtick.major.width"] = 2 # Set the length of x-axis ticks to 5
plt.rcParams["ytick.major.width"] = 2 # Set the length of y-axis ticks to 5


def main():

    st.title("IC50 Calculator")

    uploaded_file = st.file_uploader("Choose an Excel file", type="xlsx")

    if uploaded_file:
        # Check the available sheets in the Excel file
        xls = pd.ExcelFile(uploaded_file)
        sheet_names = xls.sheet_names

        # If multiple sheets, prompt user to select a sheet
        if len(sheet_names) > 1:
            selected_sheet = st.selectbox("Choose a sheet to process", sheet_names)
        else:
            selected_sheet = sheet_names[0]

        # Let user select the number of experiments
        n_experiments = st.selectbox("Select the number of experiments:", [1, 2, 3, 4, 5, 6])
        
        # Update start_row based on n_experiments
        start_row_dict = {
            1: [44],
            2: [44, 62],
            3: [44, 62, 80],
            4: [44, 62, 80, 98],
            5: [44, 62, 80, 98, 116],
            6: [44, 62, 80, 98, 116, 134],

        }
        start_row = start_row_dict[n_experiments]

        # Display the "Confirm Selection" button
        confirm_button = st.button("Confirm Selection")

        # Only execute the following if the confirm button is pressed
        if confirm_button:



    # EXTRACTS ###########
            # Extract the assay and title texts using the extract_text_data function
            assay_text, title_text = extract_head_data(uploaded_file, selected_sheet)

            # Extract compounds using the function
            compounds = extract_compounds(uploaded_file, selected_sheet)
        
            # Extract concentrations using the function
            concentrations = extract_concentrations(uploaded_file, selected_sheet)

            # Extract Y label using the function
            y_label = extract_ylabel(uploaded_file, selected_sheet)

            # Extract experimental data using the function
            experiments_data = []
            for idx in range(n_experiments):
                exp_data = extract_experiment(uploaded_file, selected_sheet, start_row= start_row[idx])
                experiments_data.append(exp_data)





    # DISPLAY ###########

            # Display the extracted texts using the display_text_data function
            display_head_data(assay_text, title_text)

            # Use the function to display the compounds in Streamlit
            colored_indices, combined_colors = display_compounds(compounds[0])


            # Display the extracted data using the display function
            display_concentrations(concentrations, colored_indices, combined_colors)


            for idx, exp_data in enumerate(experiments_data, 1):
                display_experiment(exp_data, number=idx)


    # Plotting #########

            control_idx, mocks_idx, none_idx, filtered_compounds = filter_compounds(compounds[1])
            compounds = np.array(compounds[1])

            x_axis                 = np.array(concentrations)
            experiments_data       = np.array(experiments_data)

            print(x_axis.shape)
            print(experiments_data.shape)
            print(mocks_idx)
            
            #print(experiments_data.shape)
            #print(experiments_data[0,:,0])

            CNAME_LIST   = []
            IC50_LIST    = []
            SLOPE_LIST   = []
            PDFNAME_LIST = []



            for i in range(x_axis.shape[1]):

                if i in mocks_idx or i in none_idx:
                    continue

                else:
                    xdata = np.log10(x_axis[:,i]*1e-6)
                    ydata = np.mean(experiments_data[:, :, i], axis=0)
                    yerr  = np.std(experiments_data[:, :, i], axis=0)

                    # Use curve_fit to fit the function to the data
                    p0         = (min(ydata),max(ydata), compute_x_at_ymid(xdata, ydata), -1)  # initial guess for parameters a, b, c, and d
                    popt, pcov = curve_fit(variable_slope_log_inhibitor_response, xdata, ydata, p0, method="lm", maxfev=100000)
                    IC50       = 10 ** popt[2]
                    
                    # Appending to list
                    IC50_LIST.append(str(np.round(IC50,8)))
                    SLOPE_LIST.append(str(np.round(popt[3],8)))
                    CNAME_LIST.append(compounds[idx])
                    
                    # Print the fitted parameters
                    print('Bottom  =', popt[0])
                    print('Top     =', popt[1])
                    print('LogIC50 =', popt[2])
                    print('IC50    =', 10 ** popt[2])
                    print('Slope   =', popt[3])


                    # Generate the fitted curve using the fitted parameters
                    xfit = np.linspace(min(xdata), max(xdata), num=1000)
                    yfit = variable_slope_log_inhibitor_response(xfit, *popt)



                    # Plot the data
                    fig, ax = plt.subplots()
                    plt.plot(xdata, ydata, 'ro', markersize=15, color="red", label="Response")
                    plt.errorbar(xdata, ydata, yerr=yerr, linestyle='None', capsize=6, capthick=2, elinewidth=1, color="red")

                    # Plot the fitted curve
                    plt.plot(xfit, yfit, 'b-', color="red", lw=3)

                    # Plot Mocks
                    #If there are more than 2 mocks only the first and the last one will be plotted

                    # Mock A
                    m1      = np.mean(experiments_data[:, :, mocks_idx[0]], axis=0)
                    m1_err  = np.std (experiments_data[:, :, mocks_idx[0]], axis=0)
                    
                    
                    # Mock P right
                    m3      = np.mean(experiments_data[:, :, mocks_idx[-1]], axis=0)
                    m3_err  = np.std (experiments_data[:, :, mocks_idx[-1]], axis=0)

                    plt.plot(xdata, m1, 'ro', markersize=15, color="blue", label="mock A", zorder=0, alpha=0.5)
                    plt.errorbar(xdata, m1, yerr=m1_err, linestyle='None', capsize=6, capthick=2, elinewidth=1, color="blue", zorder=0, alpha=0.5)
                    
                    #plt.plot(xdata, m2, 'ro', markersize=15, color="brown", label="mock P", zorder=0, alpha=0.5)
                    #plt.errorbar(xdata, m2, yerr=m2_err, linestyle='None', capsize=6, capthick=2, elinewidth=1, color="brown", zorder=0, alpha=0.5)
                    
                    plt.plot(xdata, m3, 'ro', markersize=15, color="black", label="mock P right", zorder=0, alpha=0.5)
                    plt.errorbar(xdata, m3, yerr=m3_err, linestyle='None', capsize=6, capthick=2, elinewidth=1, color="black", zorder=0, alpha=0.5)

                    plt.title(compounds[i], pad=20, size=25, weight = 'bold')

                    # Add a legend and axis labels
                    plt.xlabel('Compound [$log_{10}(M)$]',weight = 'bold', size=30)
                    plt.ylabel(y_label,weight = 'bold', size=30)
                    plt.ylim(0,200)
                    plt.xlim(-9,-3)
                    plt.xticks(weight = 'bold', size=30)
                    plt.yticks(weight = 'bold', size=25)

                    # Text
                    plt.text(-6, 160, "IC50="+str( '%.2E' % Decimal(IC50*1e6))+" $\mu  M$", weight = 'bold')

                    plt.legend(frameon=False, ncol=3, loc="upper left")

                    #plt.savefig("PDF/"+compounds[i]+".pdf", format="pdf", dpi=300, bbox_inches='tight')
                    
                    PDFNAME_LIST.append(compounds[i]+".pdf")



                    # Display the plot in Streamlit
                    st.pyplot(fig)
                    plt.close()  # Clear the current figure
                    st.write("")  # Add an empty line as padding
                    st.write("---")  # Draw a line for better separation (optional)
                    st.write("")  # Add another empty line as padding





if __name__ == "__main__":
    main()
