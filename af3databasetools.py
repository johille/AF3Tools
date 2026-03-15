#!/usr/bin/env python3
"""
AF3DataBaseTools - A tool to make a folder structure usable for the Foldseek Custom Database Command
"""

import argparse
import os

def extract_cif(input_folder):
    cif_paths = []
    for folder in os.listdir(input_folder):
        for file in os.listdir(input_folder + "/" + folder):
            if ".cif" in file:
                cif_paths.append(input_folder + "/" + folder + "/" + file)
    return cif_paths

def create_db_folder(cif_paths, output_folder):
    os.mkdir(output_folder)
    for file in cif_paths:
        cmd = f'cp "{file}" "{output_folder}"'
        os.system(cmd)
    pass

def main():
    argparser = argparse.ArgumentParser(description="AF3JsonTools - make a folder structure usable for the Foldseek Custom Database Command")
    argparser.add_argument("input_folder", help="Path to the Alphafold Outpuit")
    argparser.add_argument("output_folder", help="Path to the created db_folder")
    args = argparser.parse_args()

    cif_paths = extract_cif(args.input_folder)
    create_db_folder(cif_paths, args.output_folder)

if __name__ == "__main__":
    main()
