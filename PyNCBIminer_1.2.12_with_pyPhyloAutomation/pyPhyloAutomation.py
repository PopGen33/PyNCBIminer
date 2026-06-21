# This script is an attempt at automatting the process of retrieving new sequences 
# from NCBI and adding them to a phylogeny
# NOTE: Does not work with modern pandas; needs pandas=2 (or earlier)

import filecmp
import os
from datetime import date
from pathlib import Path
import shutil
from configobj import ConfigObj
import argparse
from Bio import SeqIO
import subprocess
import filecmp

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
    #today = '2015-01-01' # TESTING ~*~*~*~*~*~*~*~
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
        blast_params = config['loci'][locus]['blast']

        # Build working directories in PyNCBIminer format
        locus_working_directory = working_directory / Path(locus)
        os.makedirs(locus_working_directory / Path("parameters"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("parameters") / Path("ref_seq"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("parameters") / Path("ref_msa"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("tmp_files"), exist_ok=True)
        os.makedirs(locus_working_directory / Path("results"), exist_ok=True)

        # Copy initial quaries to their expected location in the working directory
        # so pyNCBIminer can find them
        initial_queries = list(SeqIO.parse(blast_params['initial_queries_path'], "fasta"))
        SeqIO.write(initial_queries, locus_working_directory / Path("parameters") / Path("initial_queries.fasta"), "fasta")

        # Modified from pyNCBIminer_00_main.py
        # The original code was to finish runs that hadn't completed
        # I'm abusing it a bit to run *new* BLAST rounds after iterative BLAST has already stopped
        # by passing a high blast_round number than was run. This does mean that BLAST runs will accumulate over 
        # time. May have to add some cleanup for that
        queries_file_list = os.listdir(locus_working_directory / Path("parameters") / Path("ref_seq"))
        if len(queries_file_list) > 0:
            round_list = []
            for queries_file in queries_file_list:
                round_list.append(int(Path(queries_file).stem.split("_")[-1]))
            blast_round = max(round_list) + 1 # add one more round than is present in the queries file list to force a new round of BLAST
        else:
            blast_round = 1
        # end modified from pyNCBIminer_00_main.py

        # Call the iterated blast function from pyNCBIminer
        iterated_blast_main(
            wd = locus_working_directory,
            organisms = organisms,
            count = entrez_count(
                config['DEFAULT']['entrez_email'], 
                organisms,
                blast_params['entrez_qualifier'], 
                config['DEFAULT']['last_check'], 
                today
                ),
            expect = float(blast_params['expect_value']),
            gapcosts = blast_params['gap_costs'],
            word_size = int(blast_params['word_size']),
            nucl_reward = int(blast_params['nucl_reward']),
            nucl_penalty = int(blast_params['nucl_penalty']),
            max_length = int(blast_params['max_length']),
            key_annotations = blast_params['key_annotations'].split("|"), # my config file uses pipe to delimit the key annotations
            exclude_sources = blast_params['exclude_sources'].split("|"), # my config file uses pipe to delimit the exclude sources
            ref_number = 5, # Hard-coded in original pyNCBIminer code; could be added to config?
            date_from = config['DEFAULT']['last_check'],
            date_to = today,
            entrez_email = config['DEFAULT']['entrez_email'],
            blast_round = blast_round
        )

    # Superatrix Construction Module
    # Sequence filtering; we'll bypass call_miner_filter and just make the Miner object ourselves
    for locus in config['loci']:
        print(f"Filtering sequences for locus: {locus}")
        # Get filtering parameters for this locus; if not present, get default parameters
        try:
            filter_params = config['loci'][locus]['filtering']
        except KeyError:
            filter_params = config['DEFAULT']['filtering']
        # Locus working directory
        locus_working_directory = working_directory / Path(locus)
        # if not_controlled directory exists, move to not_controlled_old (and delete old not_controlled_old if it exists)
        not_controlled_dir = locus_working_directory / Path("results") / Path("not_controlled")
        if (not_controlled_dir.parent.resolve() / Path("not_controlled_old")).exists():
                shutil.rmtree(not_controlled_dir.parent.resolve() / Path("not_controlled_old"))
        if not_controlled_dir.exists():
            not_controlled_dir.rename(not_controlled_dir.parent.resolve() / Path("not_controlled_old"))
        
        # Stop filtering if there are no new sequences
        # TODO: Find a better way to do this; this doesn't entirely work...
        # Check for the 'history_backup' directory; if it exists, get the most recent directory and blast_results_checked_seq_info_modified.txt
        # Compare that to the copy in the locus directory; if they are the same, skip the filtering step (blast results haven't changed)
        # if os.path.exists(locus_working_directory / Path("results") / Path("history_backup")) and len(os.listdir(locus_working_directory / Path("results") / Path("history_backup"))) > 0:
        #     most_recent_backup = max((locus_working_directory / Path("results") / Path("history_backup")).iterdir(), key=os.path.getmtime)
        #     last_blast_results_checked_seq_info = most_recent_backup / Path("blast_results_checked_seq_info_modified.txt")
        #     current_blast_results_checked_seq_info = locus_working_directory / Path("results") / Path("blast_results_checked_seq_info_modified.txt")
        #     if filecmp.cmp(last_blast_results_checked_seq_info, current_blast_results_checked_seq_info, shallow=False):
        #         print(f"No changes in blast results for locus {locus} since last filtering; skipping filtering step")
        #         continue

        # Create Miner_filter object
        my_miner_filter = Miner_filter(locus_working_directory, locus_working_directory)

        # Extended segments refinement
        if filter_params['extended_segments_refinement'].lower() == 'true':
            print("Extended segments refinement for locus: %s" % locus)
            my_miner_filter.control_extension(length_ratio=0.6, 
                                              max_subset_size=200, 
                                              gappyness_threshold=float(filter_params['extension_gappyness_threshold'])
            )
        
        # Species-level sequence selection
        if filter_params['species-level_sequences_selection'].lower() == 'true':
            print("Reducing dataset to a single sequence per species for locus: %s" % locus)
            # Dunno what this does, but it's in the original pyNCBIminer code
            rename_results(locus_working_directory)
            # Name correction below uses tNRS which is a plant database from what I've read; may be able to use a more general alternative like the GBIF taxonomy API if this proves necessary
            my_miner_filter.reduce_dataset(
                name_correction = filter_params['name_correction'].lower() == 'true',  # for name correction using trns (rtrns)
                subsp = True, 
                var = True, 
                f = True,
                sp = True, 
                cf = True, 
                aff = True, 
                x = True,
                consensus_value = True, # I think this is read from the "abnormal index" boolean in the GUI, but it's *really* hard to trace; mystery parameter, honestly
                length_threshold = int(filter_params['length_threshold']),
                ignore_gap = True # Read from nowhere and always true?
            )
    # below is bit lazy, but should work so long as nothing else edits the working diretory
    wd_list = [working_directory / Path(locus) for locus in config['loci']] 
    combine_keep_records(wd_list)
    put_filtered_seq_together(wd_list)

    # Alignment with MAFFT
    # Results from previous steps are in working_directory/filtered_seqs
    os.makedirs(working_directory / Path("filtered_seqs_aligned"), exist_ok=True)
    for locus in config['loci']:
        print(f"Aligning sequences for locus: {locus}")
        input_fasta = working_directory / Path("filtered_seqs") / Path(f"{locus}.fasta")
        output_fasta = working_directory / Path("filtered_seqs_aligned") / Path(f"{locus}.fasta")

        # Get filtering parameters for this locus; if not present, get default parameters
        try:
            mafft_params = config['loci'][locus]['mafft']
        except KeyError:
            mafft_params = config['DEFAULT']['mafft']

        # Determine algorithm; put it inside a list so it unpacks nicely even if it's a single argument
        algo = mafft_params['mafft_algorithm'].lower()
        if algo == 'auto':
            algo = ['--auto']
        elif algo == 'linsi':
            algo = ['--localpair']
        elif algo == 'ginsi':
            algo = ['--globalpair']
        elif algo == 'einsi':
            algo = ['--ep', '0', '--genafpair']
        elif algo == 'qinsi':
            # changes the base call to mafft-qinsi
            pass
        elif algo == 'xinsi':
            # changes the base call to mafft-xinsi
            pass
        else:
            raise ValueError(f"Invalid MAFFT algorithm specified for locus {locus}: {mafft_params['mafft_algorithm']}. Options are: auto, linsi, ginsi, einsi")
        
        # Build mafft command
        if algo in ['qinsi', 'xinsi']:
            # qinsi and xinsi are their own executable, so have to change the base call
            # These are for RNA structure-aware alignment
            # *** This requires compiling MAFFT and its extensions from source and having them available in PATH
            mafft_command = [
                f"mafft-{algo}",
                "--maxiterate", mafft_params['max_iterations'],
                "--thread", mafft_params['thread']
                ]
        else:
            mafft_command = [
                "mafft",
                *algo,
                "--maxiterate", mafft_params['max_iterations'],
                "--thread", mafft_params['thread']
            ]
        # add --dash if specified in config; this is a MAFFT option that includes structural information for protein alignments; might be useful
        # Only usable if alignments are amino acid sequences, so removed for now.
        # if mafft_params['dash'].lower() == 'true':
        #     mafft_command.append("--dash")
        # mafft_command.append(str(input_fasta))

        # Run mafft
        with open(output_fasta, "w") as output_handle:
            subprocess.run(mafft_command, stdout=output_handle)


    # Finally, update the last_check date in the config file to today's date
    print("last_check in config file updated to today's date: %s" % today)
    config['DEFAULT']['last_check'] = today
    config.write()
    