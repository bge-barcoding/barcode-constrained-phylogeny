## Scripts

### [bold_data_dump.py](bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.

### [map_opentol.py](map_opentol.py)
- Uses the [Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 

### [family_fasta.py](family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{family}.fasta'

### [create_MSA.py](create_MSA.py)
- Uses masce to create a Multiple Sequence Alignment as output files:'{family_NT}.fasta' and '{family_AA}.fasta'.

### [replace_alignment_ids.py](replace_alignment_ids.py)
- Replaces the FASTA headers in alignments from '>{barcode_id}' to '>{opentol_id}\_{barcode_id}' 

### [edit_constraint.py](edit_constraint.py)
- Edits the constraint trees to accomodate for multiple barcodes from one species and prepares the newick files for raxml.
