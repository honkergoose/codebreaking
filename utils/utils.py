# A module of utilities for solving ROSALIND problems and learning a bit of biology
from typing import List, Iterator, NamedTuple

class GeneticString(NamedTuple):
    name: str
    description: str
    data: str

def overlap_graph_connection(gs1: GeneticString, gs2: GeneticString, k: int) -> set:
    "Given two genetic strings return entries in an adjacency list format"
    adj_set = set()
    # import pdb; pdb.set_trace()
    if gs1.data == gs2.data:
        # String are the same. Do not connect.
        return adj_set
    
    if gs1.data[len(gs1.data)-k:len(gs1.data)] == gs2.data[:k]:
        # gs1 points to gs2
        adj_set.add(f"{gs1.name} {gs2.name}")

    if gs2.data[len(gs2.data)-k:len(gs2.data)] == gs1.data[:k]:
        # gs2 points to gs1
        adj_set.add(f"{gs2.name} {gs1.name}")
    return adj_set
    

PROFILE_MATRIX_MAP = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}
PROFILE_MATRIX_MAP_REV = {v: k for k, v in PROFILE_MATRIX_MAP.items()}
def build_consensus_string(pm: List[List[int]]) -> str:
    # Build the consensus string from a profile matrix
    # Profile matrix is in the format of A, C, G, T
    consensus = ""
    for col in range(len(pm[0])):
        maxidx = 0
        for row in range(len(pm)):
            if pm[row][col] >= pm[maxidx][col]:
                maxidx = row
        consensus += PROFILE_MATRIX_MAP_REV[maxidx]
    return consensus


def add_to_profile_matrix(dna: str, pm: List[List[int]]) -> List[List[int]]:
    # Given a DNA string and a profile matrix, return the profile matrix with the added string
    # Profile matrix is in the format of A, C, G, T
    for i in range(len(dna)):
        pm[PROFILE_MATRIX_MAP[dna[i]]][i] += 1
    return pm

def shared_dna_motif(dna1: str, dna2: str) -> List[str]:
    # Tired; unfinished
    # The looping structure will be a bit compilicated but I see it in my head
    # this will solve: https://rosalind.info/problems/lcsm/
    # as the follow up is a function that just goes through each subsequent
    # da string and just checks each motif from this func with dna_motif
    # Returns the list of shared motifs between the two dna strings
    # you then pick the largest and return one of them
    motifs = []
    k = i = 0
    while i < len(dna1):
        while k < len(dna2):
            if dna1[i] == dna2[k]:
                # A matching character was found. A motif
    return [1]


# Looks for a motif in a DNA strand
# A motif is a common shared interval of DNA
# For example the Alu repeat is a repeat in the human genome
# Knowing these allows you to correlate between organisms genes or sets of genes and correlate their purpose
def dna_motif(dna, motif) -> List[int]:
    # Returns positions of all found instances of the substring or motif
    # returns an empty list if no instance is found
    positions = []
    k = i = 0
    while i < len(dna):
        if dna[i] == motif[k]:
            # A matching character is found
            # If we're at the end of the motif add its position
            if k == len(motif) - 1:
                # go back to its start and record that position
                positions.append(i - k)
                # reset motif's index and take i back
                i = i - k + 1
                k = 0
                continue
            # If we're not at the end of the motif increment the index
            k += 1
            i += 1
        else:
            # No match is found.
            # reset the index of the motif
            k = 0
            i += 1
    return positions
    

# I'm tired and should re-read this one
# https://rosalind.info/problems/hamm/
def hamming_distance(dna1, dna2):
    distance = 0
    for i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            distance += 1
    return distance

# You can determined what species something is by the GC content in its genome
# Not sure why but for example eukaryotic genomes have 50%
# prokaryotes have greater than 50%.. So right there you can determine
# if something is eukaryotic or prokaryotes just by genome GC content
def gc_count(dna: str) -> float:
    count = 0
    for n in dna:
        if n in ["G", "C"]:
            count += 1
    return count

def dna_to_rna(dna: str) -> str:
    # Transcribes dna to rna
    return "".join([n if n != "T" else "U" for n in dna])

# mRNA made from DNA or messanger RNA
# it serves as a template for transation into protien
# Why uracil instead of thymine?
# they're similar with what looks like a <s>single additional carbon</s> methyl group
# for thymine?
# internet says "Uracil is energetically less expensive for the production of thymine"
# nice blog -- https://earthlingnature.wordpress.com/2012/09/29/why-thymine-instead-of-uracil/
def count_nucleotides(dna: str, nucleotide: str) -> int:
    count = 0
    for n in dna:
        if n == nucleotide:
            count += 1
    return count



def parse_fasta(source: Iterator[str]) -> Iterator[GeneticString]:
    """
    Given an iterator of lines of text in FASTA format (like a .fas or .fasta file object)
    Parse the lines into GeneticStrings assuming FASTA format
    """
    name = description = data = ""
    for line in source:
        if line[0] == ">":
            if data:
                # yield genetic string we just finished parsing before parsing a new one
                yield GeneticString(name, description, data)
                # reset data
                data = ""
            # Name is first work after >
            name = line.split(" ")[0][1:]
            # Description is everything after the first word
            description = " ".join(line.split(" ")[1:]).strip("\n")
        else:
            # Append this line to the genetic data
            data += line.strip("\n")
                
    # yield last genetic string
    yield GeneticString(name, description, data)
  
def find_reverse_palindromes(dna: str, min_p: int, max_p: int) -> List[List[int]]:
    """
    min_p: minimum sequence length to be considered a reverse palindrome
    max_p: maximum sequence length to be considered a reverse palindrome
    
    returns list[list[int, int]] where the tuple of ints
    is the start and end of the slice that can be 
    inspected by taking the original dna string and slicing it dna[x:y] where x and y are a tuple
    """
    rp_slices = []
    low = 0
    high = min_p
    # sliding window checking if each segment of length between min_p and max_p is a reverse palindrome
    while low + min_p <= len(dna):
        # while our smallest window is within the bounds of the DNA strand
        while high - low <= max_p and high <= len(dna):
            # while our window is still expanding AND is within the bounds of the DNA strand
            if (reverse_palindrome(dna[low:high])):
                rp_slices.append([low, high])
            # increase our window size
            high += 1
        # move our window start point forward in the dna strand
        low += 1
        # reset high to reset our window to minimum palindrome length
        high = low + min_p
    return rp_slices

def reverse_palindrome(dna) -> bool:
    return dna == reverse_compliment(dna)

# Reverse complements are interesting.
# Its hard to wrap your brain around the shape
# and without a grasp of oCHEM its hard to get the WHY
# Need to re-read: https://rosalind.info/problems/revc/
def reverse_compliment(dna):
    return "".join(list(reversed(dna_compliment(dna))))

# Maps DNA nucleotides to their "complement".
# Complementarity describes a relationship between two structures following the lock-and-key principle. (Fisher's theory)
# The lock-and-key is a model for enzyme-substrate interaction. The two structures have complementary geometric shapes.
# Enzymes are highly specific they must bind to a specific substrate before they can catalyze a chemical reaction.
# An Enzyme is a protien that acts as a biological cayalyst. Catalyst is a substance that increases the rate of chemical reaction.
# Protiens are...? BioJunk in a 3D structure.
# Nucleotides are organic molecules consisting of nucleoside and a phosphate.
COMPLIMENT_MAP = {
    "A" : "T",
    "T" : "A",
    "C" : "G",
    "G" : "C"
}
def dna_compliment(dna: str) -> str:
    return "".join([COMPLIMENT_MAP[s] for s in dna])