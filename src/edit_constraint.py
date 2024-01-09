from Bio import SeqIO


# def replace_newick(alignment_input, newick_input, newick_output):
#     """
#     :param dict: To retrieve the values from the dict (consisting of ott and barcode id)
#     :return:
#     """
#     headers_dict = {}
#     for record in SeqIO.parse(alignment_input, "fasta"):
#         opentol_id = record.description.split('|')[0]
#         headers_dict.setdefault(opentol_id, []).append(record.description.split('|')[1])
#     with open(newick_input, "r") as input:
#         with open(newick_output, "w+") as output:
#             ott = input.readline()
#             print(ott)
#             ott = ott.replace(":0", "")
#             ott = ott.replace(";", "")
#             for key in headers_dict.keys():
#                 print(key)
#                 if len(headers_dict[key]) > 1:
#                     ott = ott.replace(key, "(" + str(headers_dict[key]) + ")")
#                 else:
#                     ott = ott.replace(key, str(headers_dict[key]))
#                 ott = ott.replace("[", "")
#                 ott = ott.replace("]", "")
#             output.write(ott)
def replace_newick(alignment_input, newick_input, newick_output):
    """
    Replace Open Tree of Life (OTOL) IDs in a Newick tree with groups of sequences from the same species.

    :param alignment_input: Path to the input alignment file (fasta format).
    :param newick_input: Path to the input Newick tree file.
    :param newick_output: Path to the output Newick tree file with replaced IDs.
    """
    headers_dict = {}

    # Create a dictionary with OTOL IDs as keys and lists of sequence IDs as values
    for record in SeqIO.parse(alignment_input, "fasta"):
        opentol_id = record.description.split('|')[1]
        headers_dict.setdefault(opentol_id, []).append(record.description)

    # Read the input Newick tree and replace OTOL IDs with groups of sequences
    with open(newick_input, "r") as input_file:
        newick_tree = input_file.readline().strip()
        for key in headers_dict.keys():
            if len(headers_dict[key]) > 1:
                # Group sequences within the same species
                newick_tree = newick_tree.replace(key, "(" + ",".join(headers_dict[key]) + ")")
            else:
                # Replace with a single sequence ID
                newick_tree = newick_tree.replace(key, str(headers_dict[key][0]))

    # Write the modified Newick tree to the output file
    with open(newick_output, "w+") as output_file:
        output_file.write(newick_tree.rstrip(";").replace(":0", "").replace(";", ""))


if __name__ == '__main__':
    alignment_input = snakemake.input[0] # noqa: F821
    newick_input = snakemake.input[1] # noqa: F821
    newick_output = snakemake.output[0] # noqa: F821

    replace_newick(alignment_input, newick_input, newick_output)
