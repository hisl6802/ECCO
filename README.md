# Clustering toolbox
## Design
The Clustering Toolbox available here was designed as a user interface that allows for clustering analysis of untargeted metabolomics data. This user interface contains 23 functionalities spanning from hierarchical clustering to clustering comparision to KEGG pathways analysis. Most functionalities within this UI contain user options allow for customization of analyses. Furthermore, we provide a start-up window that allows the user to define the number of threads (e.g., processes they would like to use during analysis allowing for rapid clustering analysis. 

## Ensemble Clustering combined with Clustering Optimization (ECCO)
Ensemble clustering combined with Clustering Optimization or ECCO is the main functionality within the clustering toolbox. This form of ensemble clustering is discussed in our pre-print [pre-print](https://www.biorxiv.org/content/10.1101/2022.11.03.515009v1.abstract) and can be adapted for your preferences. 


-Please see the Example Files directory for examples of input files for all functionalities expect the Selected Clusters Figure functionality

### Make sure to have python installed on the command line and that pip is installed for ease of implementation.
- If you do not have python installed, install at [python.org](https://www.python.org/)


# Set-Up

## 1) Download the source code. 
`git clone https://github.com/hisl6802/ECCO`
or by downloading a zip file from the Code <> button above.

## 2) Set up a virtual environment (example of how to establish with conda below):
- It is not required that you run this in a virtual environment but highly recommended
- `conda create ECCO_env`
- `conda activate ECCO_env`

## 3) Install the packages necessary to run the UI
`pip install <package>`

Make sure to download the following packages before getting started with this UI.
- [numpy](https://numpy.org/)  
- [scipy](https://www.scipy.org/)
- [pandas](https://pandas.pydata.org/)
- [tkinter](https://docs.python.org/3/library/tkinter.html)
- [sklearn](https://scikit-learn.org/stable/index.html)
- [matplotlib](https://matplotlib.org/3.2.1/index.html)
- [fpdf](https://pyfpdf.readthedocs.io/en/latest/#:~:text=%20FPDF%20for%20Python%20%201%20Main%20features.,priority%20technical%20support%2C%20you%20can%20contact...%20More%20)
- [xlrd](https://pypi.org/project/xlrd/)
- [seaborn](https://seaborn.pydata.org/index.html)
- [mplcursors](https://pypi.org/project/mplcursors/)
- openpyxl (needed for reading in excel sheets)

## 4) After downloading the requisite packages needed. Open a shell in IDLE, open the file listed below and run the following program:

`ECCO_UI.py`

## 5) Or this program can be run from the command line by running the following command from the directory containing ECCO_UI.py

`/path/to/directory> python3 ECCO_UI.py`

## 6) Troubleshooting
I recommend that you run in the IDLE or on the command line as any missing modules will be listed in the intialization error. 


## 7) Help/Issue Reporting
Please reach out to me at brady.hislop@student.montana.edu with any issues, or troubles getting things set up. 
 
