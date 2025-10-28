#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import yaml
import sys

# Load data files
set1184 = pd.read_csv('set1184_downstairs.dat')

with open('isotope_data.yaml', 'r') as f:
    content = f.read()

# Split and parse each isotope section
isotopes = [yaml.safe_load('Isotope:' + section) 
            for section in content.split('Isotope:')[1:]]

for isotope in isotopes:
    half_life = isotope['Half-life']
    if isinstance(half_life, str):
        # Remove 'years' or 'year' and any extra whitespace
        half_life = half_life.replace('years', '').replace('year', '').strip()
        isotope['Half-life'] = float(half_life)

# Create a dictionary with isotope names as keys
isotope_dict = {isotope['Isotope']: isotope for isotope in isotopes}

# Easy access examples:
co60_energies = isotope_dict['Co-60']['Energy (keV)']
ba133_halflife = isotope_dict['Ba-133']['Half-life']

# Get all half-lives
all_half_lives = [data['Half-life'] for data in isotope_dict.values()]

# Get all energies (nested list)
all_energies = [data['Energy (keV)'] for data in isotope_dict.values()]

# Create half-life dictionary
keys = ['60-Co', '133-Ba', '137-Cs', '241-Am']
half_life_dict = dict(zip(keys, all_half_lives))

# Create nuclide activity dictionary from dataframe
nuclide_activity = {}
for index, row in set1184.iterrows():
    nuclide_activity[row.iloc[0]] = row.iloc[3]

# Strip whitespace from keys
nuclide_activity = {key.strip(): value for key, value in nuclide_activity.items()}

# Define half-life calculation function
def calculate_half_life(nuclide, elap_time=45.88):
    curie_amt = nuclide_activity[nuclide] * (1/2)**(elap_time/half_life_dict[nuclide])
    photon_amt = curie_amt * 37000
    print(f'{curie_amt:.2f} is the current activity in uCi')
    print(f'{photon_amt:.2f} is the photons per second given the activity')
    print(f'Half life of {nuclide} is {half_life_dict[nuclide]} years')


def main():
    if len(sys.argv) > 1:
        nuclide = sys.argv[1]
        calculate_half_life(nuclide)
    else:
        print("Usage: python your_script.py <nuclide>")
        print("Example: python your_script.py 60-Co")


if __name__ == "__main__":
    main()