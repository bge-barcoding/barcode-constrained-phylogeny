import logging
import tempfile
import argparse
import sqlite3
import os

from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
from subprocess import run

logging.basicConfig()
logger = logging.getLogger('msa_hmm')


def correct_revcom(hmmfile, seqfile):
    """
    Read a fasta file per record, align record using hmmer with and without reverse complementing the sequence,
    check which one to keep(highest "*" count returned from with align_score()),
    and write all sequences together to a FASTA file.
    :param hmmfile: HMM model file
    :param seqfile: name of FASTA file input
    :return:
    """
    logger.info(f'Aligning sequences in FASTA file {seqfile}')

    # Read each sequence record from fasta file
    sequences = []
    log_reverse = 0
    for record in SeqIO.parse(seqfile, "fasta"):
        # Run align_score with original sequence
        count_1, alignment_original = align_score(record, hmmfile)
        logger.debug(f'Forward alignment score {count_1} for {record.id}')

        # Reverse complement the sequence and align
        record.seq = record.seq.reverse_complement()
        count_2, alignment_reverse = align_score(record, hmmfile)
        logger.debug(f'Reverse complemented alignment score {count_2} for {record.id}')

        # Check if the first count is higher
        if count_1 >= count_2:

            # Reverse the reverse complement
            record.seq = record.seq.reverse_complement()
            sequences.append(record)
        else:

            # Add alignment to list
            sequences.append(record)

            # Log counter
            log_reverse += 1

    logger.info(f'Corrected {log_reverse} reverse complemented sequences out of {len(sequences)}')
    return sequences


def align_write(sequences, outfile, conn):
    """
    Aligns the reverse complement-corrected sequences. Remaps barcode_id to process_id using a database
    lookup by way of the provided conn. Writes output to outfile.
    :param sequences:
    :param outfile:
    :param conn:
    :return:
    """

    # Save the sequence to the temporary file
    with open(f'{outfile}.tmp1', mode='w+') as temp_fasta:
        for seq in sequences:
            temp_fasta.write(f'>{seq.id}\n')
            temp_fasta.write(f'{seq.seq}\n')

    # Align with hmmalign, capture and parse output as phylip
    run(['hmmalign', '--trim', '-o', f'{outfile}.tmp2', '--outformat', 'phylip', hmmfile, f'{outfile}.tmp1'])
    aligned = read_alignment(f'{outfile}.tmp2', 'phylip')
    os.remove(f'{outfile}.tmp1')
    os.remove(f'{outfile}.tmp2')

    # Map barcode_id to process_id and write to outfile
    with open(f'{outfile}', mode='w') as output:
        for seq in aligned:
            bid = seq.id.split('|')[0]

            process_id = conn.execute(f'SELECT processid FROM barcode WHERE barcode_id={bid}').fetchone()
            opentol_id = conn.execute(f'SELECT taxon.opentol_id FROM taxon JOIN barcode ON taxon.taxon_id = barcode.taxon_id WHERE barcode.barcode_id={bid}').fetchone()
            output.write(f'>{process_id[0]}|{opentol_id[0]}\n')
            output.write(f'{seq.seq}\n')

    logger.info(f'Wrote all aligned sequences to {outfile}')


def align_score(record, hmmfile):
    """
    Uses a Hidden Markov Model to align a sequence using hmmalign. This gives an stockholm file as output.
    Returns the number of "*" in the stockholm file. Reads the output file and saves its contents in the
    variable 'alignment'.
    :param record: a sequence from the fasta input
    :param hmmfile: HMM model file
    :return: int number of "*" in a stockholm file, content of the hmmalignment output file
    """
    # Open two temporary files to store a fasta sequence and a stockholm sequence
    with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta, tempfile.NamedTemporaryFile(mode='w+') as temp_stockholm:

        # Save the sequence to the temporary file
        SeqIO.write(record, temp_fasta.name, 'fasta')

        # Run hmm align, read the aligned sequence
        run(['hmmalign', '--trim', '-o', temp_stockholm.name, hmmfile, temp_fasta.name])
        alignment = read_alignment(temp_stockholm.name, "stockholm")
        seq = alignment[0]

    # Return probability colum of the stockholm file
    quality_string = alignment.column_annotations['posterior_probability']

    count = 0

    # Count . and * characters
    dot_count = quality_string.count('.')
    star_count = quality_string.count('*')

    # Give value 0 to . and value 10 to *
    count += dot_count * 0
    count += star_count * 10

    # Add all numbers
    digit_sum = sum(int(char) for char in quality_string if char.isdigit())

    # Add numbers to the values calculated from . and *
    count += digit_sum

    # Calculate average count to account for gaps (0's)
    average_count = count / len(quality_string)
    return average_count, seq


if __name__ == '__main__':
    hmmfile = snakemake.params.hmm  # noqa: F821
    database = snakemake.input[0] # noqa: F821
    in_file = snakemake.input[1]  # noqa: F821
    out_file = snakemake.output[0]  # noqa: F821
    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {database}")
    connection = sqlite3.connect(database)

    # Check the HMM file
    logger.info(f"Going to use Hidden Markov Model from {hmmfile}")

    # Check the input file
    logger.info(f"Going to align sequences from input file {in_file}")

    # Announce the output file

    logger.info(f"Will write alignment to output file {out_file}")

    # Do the reverse complement correction
    seqs = correct_revcom(hmmfile, in_file)

    # Write alignment
    align_write(seqs, out_file, connection)
