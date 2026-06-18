# This script is an attempt at automatting the process of retrieving new sequences 
# from NCBI and adding them to a phylogeny
# NOTE: Does not work with modern pandas; needs pandas=2 (or earlier)

import os
from datetime import date
from pathlib import Path
import shutil
from configobj import ConfigObj
import argparse
from Bio import SeqIO

## Imports from pyNCBIminer; Did this by copying the latest version to a new directory
# It seems that the things I'm importing can find the things *they* import easier this way
# TODO: figure out a way to submit this as a pull request to pyNCBIminer

from iterated_blast import iterated_blast_main
from my_entrez import entrez_count
from miner_filter import Miner_filter
from my_filter import rename_results, combine_keep_records, put_filtered_seq_together

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
    # TODO: Write the configspec file so that the config file can be validated

    # get today's date
    today = date.today().strftime('%Y-%m-%d')

    # get organisms list from config
    organisms = config['DEFAULT']['organisms'].strip().split(", ") # Orginal uses splitlines(), but comma delimited in fine in the config file

    # Get working direcotry from config and create it if it doesn't exist
    working_directory = Path(config['DEFAULT']['working_directory'])

    # Make working directory
    os.makedirs(Path(working_directory), exist_ok=True)
    
    # Sequence Retrieval Module
    # Iterate over loci in config
    for locus in config['loci']:
        print(f"Processing locus: {locus}")

        # Retrieve BLAST parameters for this locus
        params = config['loci'][locus]['blast']

        # Build working directories in PyNCBIminer format
        locus_working_directory = working_directory / Path(locus)
        os.makedirs(locus_working_directory / Path("parameters"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("parameters") / Path("ref_seq"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("parameters") / Path("ref_msa"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("tmp_files"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("results"), exist_ok=True)

        # Copy initial quaries to their expected location in the working directory
        # so pyNCBIminer can find them
        initial_queries = list(SeqIO.parse(params['initial_queries_path'], "fasta"))
        SeqIO.write(initial_queries, locus_working_directory / Path("parameters") / Path("initial_queries.fasta"), "fasta")

        # Get blast_round; near verbatim from pyNCBIminer_00_main.py
        queries_file_list = os.listdir(locus_working_directory / Path("parameters") / Path("ref_seq"))
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
            wd = locus_working_directory,
            organisms = organisms,
            count = entrez_count(
                config['DEFAULT']['entrez_email'], 
                organisms,
                params['entrez_qualifier'], 
                config['DEFAULT']['last_check'], 
                today
                ),
            expect = float(params['expect_value']),
            gapcosts = params['gap_costs'],
            word_size = int(params['word_size']),
            nucl_reward = int(params['nucl_reward']),
            nucl_penalty = int(params['nucl_penalty']),
            max_length = int(params['max_length']),
            key_annotations = params['key_annotations'].split("|"), # my config file uses pipe to delimit the key annotations
            exclude_sources = params['exclude_sources'].split("|"), # my config file uses pipe to delimit the exclude sources
            ref_number = 5, # Hard-coded in original pyNCBIminer code; could be added to config?
            date_from = config['DEFAULT']['last_check'],
            date_to = today,
            entrez_email = config['DEFAULT']['entrez_email'],
        )

    # Superatrix Construction Module
    # Sequence filtering; we'll byupass call_miner_filter and just make the Miner object ourselves
    for locus in config['loci']:
        print(f"Filtering sequences for locus: {locus}")
        # Get filtering parameters for this locus; if not present, get default parameters
        try:
            params = config['loci'][locus]['filtering']
        except KeyError:
            params = config['DEFAULT']['filtering']
        # Locus working directory
        locus_working_directory = working_directory / Path(locus)
        # if not_controlled directory exists, move to not_controlled_old (and delete old not_controlled_old if it exists)
        not_controlled_dir = locus_working_directory / Path("results") / Path("not_controlled")
        if (not_controlled_dir.parent.resolve() / Path("not_controlled_old")).exists():
                shutil.rmtree(not_controlled_dir.parent.resolve() / Path("not_controlled_old"))
        if not_controlled_dir.exists():
            not_controlled_dir.rename(not_controlled_dir.parent.resolve() / Path("not_controlled_old"))

        # Create Miner_filter object
        my_miner_filter = Miner_filter(locus_working_directory, locus_working_directory)

        # Extended segments refinement
        if params['extended_segments_refinement'].lower() == 'true':
            print("Extended segments refinement for locus: %s" % locus)
            my_miner_filter.control_extension(length_ratio=0.6, 
                                              max_subset_size=200, 
                                              gappyness_threshold=float(params['extension_gappyness_threshold'])
            )
        
        # Species-level sequence selection
        if params['species-level_sequences_selection'].lower() == 'true':
            print("Reducing dataset to a single sequence per species for locus: %s" % locus)
            # Dunno what this does, but it's in the original pyNCBIminer code
            rename_results(locus_working_directory)
            # Name correction below uses tNRS which is a plant database from what I've read; may be able to use a more general alternative like the GBIF taxonomy API if this proves necessary
            my_miner_filter.reduce_dataset(
                name_correction = params['name_correction'].lower() == 'true',  # for name correction using trns (rtrns)
                subsp = True, 
                var = True, 
                f = True,
                sp = True, 
                cf = True, 
                aff = True, 
                x = True,
                consensus_value = True, # I think this is read from the "abnormal index" boolean in the GUI, but it's *really* hard to trace; mystery parameter, honestly
                length_threshold = int(params['length_threshold']),
                ignore_gap = True # Read from nowhere and always true?
            )
    # below is bit lazy, but should work so long as nothing else edits the working diretory
    wd_list = [working_directory / Path(locus) for locus in config['loci']] 
    combine_keep_records(wd_list)
    put_filtered_seq_together(wd_list)    

    # Finally, update the last_check date in the config file to today's date
    print("last_check in config file updated to today's date: %s" % today)
    config['DEFAULT']['last_check'] = today
    config.write()

    