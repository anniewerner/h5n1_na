"""
This script will read in the NA alignment file, translate the sequence to amino acids,
find the beginning of the 370-, 400-, and 430-loops, and indicate the sequences at the
three loops that make up the secondary sialic acid binding site (sbs).
This produces 3 json files, which yields the actual sequence at each of the three loops
in the sbs. These can  be used in augur export so that the amino acid sequences for
each of the three loops show up in auspice as a color by.
"""

"""This script was heavily influenced by the work of the Moncla Laboratory's script
for annotating HA furin cleavage sites and their sequences for HPAI A/H5N1 HA
sequences. We acknowledge the advancements made by the Moncla group and thank
them for sharing their tools such that they are adaptable to other functionalities,
such as the script below."""

import Bio
from Bio import SeqIO
import json

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--alignment', type=str, help='alignment file output by rule augur align')
parser.add_argument('--loop370_sequence', type=str, help='name of output json file that annotates the 370 loop sequence for each tip')
parser.add_argument('--loop400_sequence', type=str, help='name of output json file that annotates the 400 loop sequence for each tip')
parser.add_argument('--loop430_sequence', type=str, help='name of output json file that annotates the 430 loop sequence for each tip')

args = parser.parse_args()
alignment = args.alignment
loop370_sequence_json = args.loop370_sequence
loop400_sequence_json = args.loop400_sequence
loop430_sequence_json = args.loop430_sequence

def translate_nucleotide_to_aa(input_nucleotide_str):
    # convert to sequence object 
    sequence = Bio.Seq.Seq(input_nucleotide_str)
    # translate
    aa_sequence = str(sequence.translate())
    return(aa_sequence)


def return_loop370_start_position(aa_sequence):
    
    loop370_begin = "WIGR"
    
    # if .find does not match part of the string, it will return a value of -1 
    start_pos_loop370_aa = aa_sequence.find(loop370_begin)+4
    start_pos_loop370_nt = start_pos_loop370_aa*3
    return(start_pos_loop370_nt)


def output_loop370_site_aa_sequence(loop370_nt_start, nt_sequence):

    """In the following function, we start from the 4 aa downstream
    of the conserved 'WIGR' motif (positions 358-361), or 12 nts
    downstream. We will then collect 24 nts to cover the 8-aa range
    that makes up the 370 loop. Since avian H5N1 NA sequences often
    have an additional 5 aa preceding the starting methionine,
    this method helps correct for sequencing errors that alter
    sequence numbering."""

    # find the 24 nucleotides downstream of the 370-loop start
    loop370_site_nts = []
    current_site = loop370_nt_start
    # starting at the current base, move downstream until we accumulate 24 nts that aren't Ns
    while len(loop370_site_nts) < 24:
        nt_site = nt_sequence[current_site+1]
        if nt_site != "N":
            loop370_site_nts.append(nt_site)
        current_site += 1
            
    # make into list, translate it, and join into a string
    loop_370_nts = "".join(loop370_site_nts[::+1])
    loop370_sequence = Bio.Seq.Seq(loop_370_nts)
    loop370_site = str(loop370_sequence.translate())
    
    return(loop370_site)

def return_loop400_start_position(aa_sequence):

    loop400_begin = "GYSG"
    
    # if .find does not match part of the string, it will return a value of -1
    start_pos_loop400_aa = aa_sequence.find(loop400_begin)-2
    start_pos_loop400_nt = start_pos_loop400_aa*3
    return(start_pos_loop400_nt)
    
def output_loop400_site_aa_sequence(loop400_nt_start, nt_sequence):

    """In the following function, we start from the 2 aa upstream
    of the conserved 'GYSG' motif (positions 401-404), or 6 nts
    upstream. We will then collect 21 nts to cover the 7-aa range that
    makes up the 400-loop. This is needed to correct for the inclusion
    of 5 additional amino acids in a large portion of H5N1 NA sequences
    prior to the starting methionine."""
    
    # find the 21 nucleotides downstream of the 400-loop start
    loop400_site_nts = []
    current_site = loop400_nt_start
    # starting at the current base, move downstream until we accumulate 21 nts that aren't Ns
    while len(loop400_site_nts) < 21:
        nt_site = nt_sequence[current_site+1]
        if nt_site != "N":
            loop400_site_nts.append(nt_site)
        current_site +=1
    
    # make into list, translate it, and join into a string
    loop_400_nts = "".join(loop400_site_nts[::+1])
    loop400_sequence = Bio.Seq.Seq(loop_400_nts)
    loop400_site = str(loop400_sequence.translate())
    
    return(loop400_site)

def return_loop430_start_position(aa_sequence):

    # We will use the "GYSG" motif used for finding the 400-loop to orient to position 430.
    loop430_begin = "GYSG"

    # if .find for the 400-loop returns a negative value, it will not be able to find the 430-loop either
    # we will move 25 aa upstream of the "GYSG" motif so we begin at position 429
    start_pos_loop430_aa = aa_sequence.find(loop430_begin)+25
    start_pos_loop430_nt = start_pos_loop430_aa*3
    return(start_pos_loop430_nt)

def output_loop430_site_aa_sequence(loop430_nt_start, nt_sequence):

    """In the following function, we start 25 aa downstream of the
    conserved 'GYSG' motif (positions 401-404), or 75 nts downstream.
    We will then collect 15 nts to cover the 5-aa range that makes
    up the 430-loop. This is needed to correct for the inclusion
    of 5 additional amino acids in a large portion of H5N1 NA sequences
    prior to the starting methionine."""
    
    # find the 15 nucleotides downstream of the 430-loop start
    loop430_site_nts = []
    current_site = loop430_nt_start
    # starting at the current base, move downstream until we accumulate 15 nts that aren't Ns
    while len(loop430_site_nts) < 15:
        nt_site = nt_sequence[current_site+1]
        if nt_site != "N":
            loop430_site_nts.append(nt_site)
        current_site +=1
        
    # make into list, translate it, and join into a string
    loop_430_nts = "".join(loop430_site_nts[::+1])
    loop430_sequence = Bio.Seq.Seq(loop_430_nts)
    loop430_site = str(loop430_sequence.translate())
    
    return(loop430_site)
    

def output_sbs_site_jsons(alignment, output_json1, output_json2, output_json3):
    
    output_dict_370seq = {"nodes":{}}
    output_dict_400seq = {"nodes":{}}
    output_dict_430seq = {"nodes":{}}
    
    with open(output_json1, "w") as outfile: 
        outfile.write("")
    with open(output_json2, "w") as outfile:
        outfile.write("")
    with open(output_json3, "w") as outfile:
        outfile.write("")

    for seq in SeqIO.parse(alignment, "fasta"):

        strain_name = seq.description

        # convert gaps to Ns to avoid translation errors
        nt_sequence = str(seq.seq).upper().replace("-","N")

        # translate to amino acids and find the start of 370-, 400-, and 430-loops in aa and nt coordinates
        aa_sequence = translate_nucleotide_to_aa(nt_sequence)

        start_pos_loop370_nt = return_loop370_start_position(aa_sequence)
        start_pos_loop400_nt = return_loop400_start_position(aa_sequence)
        start_pos_loop430_nt = return_loop430_start_position(aa_sequence)
        
        # if .find returns a null, it outputs a value of -1. This means that there is not adequate sequence data at the 370-loop, the start_pos will be a negative number. If it is positive, we can annotate. We will apply the same ideology to the 400- and 430-loops below.
        if start_pos_loop370_nt > 4:

            # collect nts and translate the 370-loop sequence
            loop370_site = output_loop370_site_aa_sequence(start_pos_loop370_nt, nt_sequence)
        
        else:
            loop370_site = "NNNNNNNN"
        
        output_dict_370seq["nodes"][strain_name] = {"loop370_site_sequence":loop370_site.replace("X","-")}
        
        # now for the 400-loop...
        if start_pos_loop400_nt > -3:
        
            # collect nts and translate the 400-loop sequence
            loop400_site = output_loop400_site_aa_sequence(start_pos_loop400_nt, nt_sequence)
        
        else:
            loop400_site = "NNNNNNN"
            
        output_dict_400seq["nodes"][strain_name] = {"loop400_site_sequence":loop400_site.replace("X","-")}
        
        # and finally for the 430-loop...
        if start_pos_loop430_nt > 0:
        
            # collect nts and translate the 430-loop sequence
            loop430_site = output_loop430_site_aa_sequence(start_pos_loop430_nt, nt_sequence)
        
        else:
            loop430_site = "NNNNN"
        
        output_dict_430seq["nodes"][strain_name] = {"loop430_site_sequence":loop430_site.replace("X","-")}
         
        
    f = open(output_json1, "w")
    json.dump(output_dict_370seq, f)
    f.close()
    
    f = open(output_json2, "w")
    json.dump(output_dict_400seq, f)
    f.close()
    
    f = open(output_json3, "w")
    json.dump(output_dict_430seq, f)
    f.close()

output_sbs_site_jsons(alignment, loop370_sequence_json, loop400_sequence_json, loop430_sequence_json)
