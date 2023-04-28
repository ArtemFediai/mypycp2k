# CP2K Functional Names

This folder contains a collection of functional names organized by CP2K version, category, and functional type. The organization of functional names has been pre-generated using a Python script `get_page.py` from this folder to ensure the data corresponds to the specific cp2k version. The included files and folders are already generated to avoid security risks associated with accessing third-party websites on the fly and to account for potential internet unavailability.

## Organization

The functional names are organized into the following categories:

1. `c`: Correlation functionals
2. `x`: Exchange functionals
3. `xc`: Exchange-correlation functionals

For each CP2K version, a folder is created, e.g., `2022.1`, containing an `all.txt` file with all functional names for that version. The `c.txt`, `x.txt`, and `xc.txt` files are saved directly in the version folder, containing the functional names for correlation, exchange, and exchange-correlation functionals, respectively.

Furthermore, functional names are saved in more specific subcategory text files based on the type of functionals:

1. `lda`: Local Density Approximation (LDA) functionals
2. `gga`: Generalized Gradient Approximation (GGA) functionals
3. `mgga`: Meta-Generalized Gradient Approximation (MGGA) functionals
4. `hyb`: Hybrid functionals

These subcategory files are stored within the `c`, `x`, and `xc` subfolders for each CP2K version.

## Usage

The `mypycp2k` package is shipped with these pre-generated files and folders. CP2K Engine or users can easily access the functional names according to their category (correlation, exchange, or exchange-correlation), functional type (LDA, GGA, MGGA, or Hybrid), and CP2K version without the need for internet access or worrying about security risks associated with fetching data from third-party websites.

Feel free to browse the folders and files to find the functional names you are interested in for a specific CP2K version, category, and functional type.


## For developers.
If you are from year 2024+ and mypycp2k was not updated for some reason, add the current cp2k version and its identifier in the versions.csv file.
