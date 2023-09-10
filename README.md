# Overlap calculator
This script can calculate the overlap between two pi systems as described in my PhD thesis. There are multiple ways to use this app.

1. A live version should exist at [https://overlap.streamlit.app/]

2. You can use the GUI version offline by following these instructions
 - Installation/start:
   1.  Unzip the file in a folder (Example A:\unzip_folder)
   2.  The easiest way to install is using miniconda or Anaconda.
      *   Install the python requirements with: "conda env create -f  A:\unzip_folder\environment.yml".
      *   Activate the environment with conda activate overlap (This needs to be done every time you start a console)
   3.  Alternatively, you can install the packages manually.
       Requirements are: numpy, pandas, shapely, matplotlib and streamlit
   4.  Start the program with: "streamlit run A:\unzip_folder\overlap_streamlit.py". It might ask you for a mail
       address, which is skippable. Afterwards a browser tab should open automatically if it does not, the console has
       the correct web address.
3. Starting the command line program with "python overlap.py"

## References
If you have used these scripts in your research cite the following two references
 - P. N. Ruth, Method Development for Benchmarking Key Interactions in Quantum Crystallography, p. 39, http://dx.doi.org/10.53846/goediss-9798

 - To be added when on github

