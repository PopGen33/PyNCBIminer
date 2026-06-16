# This script is an attempt at automatting the process of retrieving new sequences 
# from NCBI and adding them to a phylogeny
# NOTE: Does not work with modern pandas; needs pandas=2 (or earlier)

import os
from datetime import date
from pathlib import Path
from configobj import ConfigObj
import argparse
from Bio import SeqIO

## Imports from pyNCBIminer; Did this by copying the latest version to a new directory
# It seems that the things I'm importing can find the things *they* import easier this way
# TODO: figure out a way to submit this as a pull request to pyNCBIminer

from iterated_blast import iterated_blast_main
from my_entrez import format_entrez_query, entrez_count, entrez_summary

## pyNCBIminer imports end here


if __name__ == "__main__":
    print("Retrieving new sequences from NCBI...")
    
    # Parse arguments (get config file path)
    parser = argparse.ArgumentParser(description="Automates adding new sequences to alignments and supermatrix.")
    parser.add_argument("--config", 
                        help="Path to the config file",
                        required=True)
    args = parser.parse_args()

    # Read the config file
    config = ConfigObj(args.config)

    # get today's date
    today = date.today().strftime('%Y-%m-%d')

    # get organisms list from config
    organisms = config['DEFAULT']['organisms'].strip().split(", ") # Orginal uses splitlines(), but comma delimited in fine in the config file

    # Get working direcotry from config and create it if it doesn't exist
    working_directory = Path(config['DEFAULT']['working_directory'])
    # Near verbatim from pyNCBIminer_00_main.py to create working directory paths
    if not os.path.exists(Path(working_directory)):
        os.makedirs(Path(working_directory))
    if not os.path.exists(Path(working_directory) / Path("parameters")):
        os.makedirs(Path(working_directory) / Path("parameters"))
    if not os.path.exists(Path(working_directory) / Path("parameters") / Path("ref_seq")):
        os.makedirs(Path(working_directory) / Path("parameters") / Path("ref_seq"))
    if not os.path.exists(Path(working_directory) / Path("parameters") / Path("ref_msa")):
        os.makedirs(Path(working_directory) / Path("parameters") / Path("ref_msa"))
    if not os.path.exists(Path(working_directory) / Path("tmp_files")):
        os.makedirs(Path(working_directory) / Path("tmp_files"))
    if not os.path.exists(Path(working_directory) / Path("results")):
        os.makedirs(Path(working_directory) / Path("results"))
    # end verbatim from pyNCBIminer_00_main.py
    
    # Iterate over loci in config
    for locus in config['loci']:
        print(f"Processing locus: {locus}")

        # Retrieve BLAST parameters for this locus
        params = config['loci'][locus]['blast']

        # Copy initial quaries to their expected location in the working directory
        # so pyNCBIminer can find them
        initial_queries = list(SeqIO.parse(params['initial_queries_path'], "fasta"))
        SeqIO.write(initial_queries, Path(working_directory) / Path("parameters") / Path("initial_queries.fasta"), "fasta")

        # Get blast_round; near verbatim from pyNCBIminer_00_main.py
        queries_file_list = os.listdir(Path(working_directory) / Path("parameters") / Path("ref_seq"))
        if len(queries_file_list) > 0:
            round_list = []
            for queries_file in queries_file_list:
                round_list.append(int(Path(queries_file).stem.split("_")[-1]))
            round_list.sort()
            blast_round = round_list[-1]
        else:
            blast_round = 1
        # end verbatim from pyNCBIminer_00_main.py

        # Call the iterated blast function from pyNCBIminer
        iterated_blast_main(
            wd = working_directory,
            organisms = organisms,
            count = entrez_count(
                config['DEFAULT']['entrez_email'], 
                organisms,
                params['entrez_qualifier'], 
                config['DEFAULT']['last_check'], 
                today
                ),
            expect = float(params['expect_value']),
            gapcosts = '2 1',
            #gapcosts = params['gap_costs'],
            word_size = int(params['word_size']),
            nucl_reward = int(params['nucl_reward']),
            nucl_penalty = int(params['nucl_penalty']),
            max_length = int(params['max_length']),
            key_annotations = params['key_annotations'],
            exclude_sources = params['exclude_sources'],
            ref_number = 5, # Hard-coded in original pyNCBIminer code; could be added to config?
            date_from = config['DEFAULT']['last_check'],
            date_to = today,
            entrez_email = config['DEFAULT']['entrez_email'],
        )
    