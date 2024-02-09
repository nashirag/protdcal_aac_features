# protdcal_aac_features
Quick little script to generate Protdcal and amino acid composition (AAC) features for an input file of peptides. 

Friendly to those unfamiliar with python + command line!

## Requirements
### Python 3+
First ensure you have python installed, version 3 or later. 
To check your python version in command line:
```
python --version
```
If you don't have python, [install it here](https://www.python.org/downloads/)!
### pip
You'll need pip to install the required python packages. Check if you have pip installed first in command line:
```
pip --version
```
To install pip in command line:
```
python -m ensurepip --upgrade
```
### Pandas, NumPy, and ProPythia
Now that you have pip and python installed, let's install the libraries we need. Run the following in command line:
```
pip install pandas
pip install numpy
pip install propythia
```
Once these are good to go, you can carry on with running the application!


## Basic Usage
Download the **generate_features folder**. 

Once in the generate_features folder within command line, use the following to run:
```
python generate_features.py -p "<column>" -i <input file> -o <output file>
```
  **column** - the name of the column in your file that contains the peptide sequences (be sure to put it in quotes!)

  **input file** - the location + name of the file containing your peptide sequences

  **output file** - where you want to save the peptide sequences to (include a file name and file extension if desired, i.e. features.xlsx, features.csv, or features.tsv)
