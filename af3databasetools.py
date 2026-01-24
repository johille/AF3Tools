#!/usr/bin/env python3
"""
AF3DataBaseTools - A tool to make a folder structure usable for the Foldseek Custom Database Command
"""

import argparse

def extract_cif():
    pass

def create_db_folder():
    pass

def main():
    "Main function to facilitate command line arguments."
    argparser = argparse.ArgumentParser(description="AF3JsonTools - make a folder structure usable for the Foldseek Custom Database Command")
    argparser.add_argument("input_folder", help="Path to the Alphafold Outpuit")
    argparser.add_argument("output_folder", help="Path to the created db_folder")
    args = argparser.parse_args()

if __name__ == "__main__":
    main()
