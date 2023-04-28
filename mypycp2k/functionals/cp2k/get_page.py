"""
This script reads the input file versions.csv, which contains CP2K version names and their corresponding version identifiers. For each version, the script fetches the HTML content of the XC_FUNCTIONAL webpage using the specified URL format. It then extracts the functional names from the webpage content and saves them into text files, organized by category and subcategory.

The script creates a folder for each version, e.g., 2022.1, and saves an all.txt file containing all functional names for that version. The functional names are organized into three main categories:

    c: Correlation functionals
    x: Exchange functionals
    xc: Exchange-correlation functionals

It creates subfolders c, x, and xc for each version, and saves the corresponding functional names into c.txt, x.txt, and xc.txt files directly in the version folder.

Additionally, the script saves functional names in more specific subcategory text files based on the type of functionals:

    lda: Local Density Approximation (LDA) functionals
    gga: Generalized Gradient Approximation (GGA) functionals
    mgga: Meta-Generalized Gradient Approximation (MGGA) functionals
    hyb: Hybrid functionals

These subcategory files are stored within the c, x, and xc subfolders.

This organization enables users to easily find and access functional names according to their category (correlation, exchange, or exchange-correlation), functional type (LDA, GGA, MGGA, or Hybrid), and CP2K version.
"""


import requests
import re
import os
import csv

def save_names_to_files(names, folder):
    os.makedirs(folder, exist_ok=True)

    with open(os.path.join(folder, 'all.txt'), 'w') as f:
        for name in names:
            f.write(name + '\n')

    categories = {
        'xc': '_XC_',
        'c': '_C_',
        'x': '_X_',
    }

    category_names = {}
    for cat, keyword in categories.items():
        os.makedirs(os.path.join(folder, cat), exist_ok=True)
        with open(os.path.join(folder, f'{cat}.txt'), 'w') as f:
            category_names[cat] = [name for name in names if keyword in name]
            for name in category_names[cat]:
                f.write(name + '\n')

    subcategories = {
        'lda': 'LDA_',
        'gga': 'GGA_',
        'mgga': 'MGGA_',
        'hyb': 'HYB_',
    }

    for cat in categories.keys():
        for subcat, prefix in subcategories.items():
            with open(os.path.join(folder, cat, f'{subcat}.txt'), 'w') as f:
                for name in category_names[cat]:
                    if name.startswith(prefix):
                        f.write(name + '\n')


def get_names_from_url(url):
    response = requests.get(url)
    html_content = response.text
    names = re.findall(r'<a href="XC_FUNCTIONAL/[^"]+">(.*?)<\/a>', html_content)
    return names


with open('versions.csv', newline='') as csvfile:
    versions_reader = csv.DictReader(csvfile)
    for row in versions_reader:
        version = row['version']
        name = row['name']
        url = f"https://manual.cp2k.org/{name}/CP2K_INPUT/ATOM/METHOD/XC/XC_FUNCTIONAL.html"
        names = get_names_from_url(url)
        save_names_to_files(names, version)

