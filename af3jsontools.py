#!/usr/bin/env python3
"""
AF3JsonTools - A tool to prepare and sanitize JSON files for AF3 processing.
"""

import json
import argparse

def trim_cutoff(limit_aa_length, aa_sequence):
    "omit sequences in the provided json file for a certain cutoff"
    pass

def main():
    "Main function to facilitate command line arguments."
    argparser = argparse.ArgumentParser(description="AF3JsonTools - Prepare and sanitize JSON files for AF3 processing.")
    argparser.add_argument("input_file", help="Path to the input JSON file.")
    argparser.add_argument("output_file", help="Path to the output sanitized JSON file.")
    args = argparser.parse_args()
    

if __name__ == "__main__":
    main()
