# Unit tests for utils.py

from utils import utils

def test_dna_compliment():
    dna = "ATCG"
    compliment = "TAGC"
    assert utils.dna_compliment(dna) == compliment
    
def test_reverse_compliment():
    dna = "ATCG"
    r_compliment = "CGAT"
    assert utils.reverse_compliment(dna) == r_compliment
    
def test_reverse_palindrome():
    rp_dna = "GCATGC"
    not_rp_dna = "ATGCGGGTCTATAT"
    assert utils.reverse_palindrome(rp_dna)
    assert not utils.reverse_palindrome(not_rp_dna)
    
def test_find_reverse_palindromes():
    dna = "TCAATGCATGCGGGTCTATATGCAT"
    min_p = 4
    max_p = 12
    answer = [
        [3, 9], [4, 8], [5, 11], [6, 10], [16, 20], [17, 21], [19, 25], [20, 24]
    ]
    assert utils.find_reverse_palindromes(dna, min_p, max_p) == answer
    
def test_parse_fasta():
    text = """>Taxon1 The first one
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Taxon2 The second one
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Taxon3 The third one
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
    lines = text.split("\n")
    results = [genetic_string for genetic_string in utils.parse_fasta(lines)]
    assert results[0].name == "Taxon1"
    assert results[0].description == "The first one"
    assert results[0].data == "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"
    assert results[1].name == "Taxon2"
    assert results[1].description == "The second one"
    assert results[1].data == "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"
    assert results[2].name == "Taxon3"
    assert results[2].description == "The third one"
    assert results[2].data == "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
        
def test_dna_to_rna():
    dna = "GATGGAACTTGACTACGTAAATT"
    rna = "GAUGGAACUUGACUACGUAAAUU"
    assert utils.dna_to_rna(dna) == rna

def test_dna_motif():
    dna = "GATATATGCATATACTT"
    motif = "ATAT"
    output = [1, 3, 9]
    assert utils.dna_motif(dna, motif) == output

def test_build_consensus_string():
    consensus_string = "ATGCAACT"
    profile_matrix = [
        [5, 1, 0, 0, 5, 5, 0, 0],
        [0, 0, 1, 4, 2, 0, 6, 1],
        [1, 1, 6, 3, 0, 1, 0, 0],
        [1, 5, 0, 0, 0, 1, 1, 6]
    ]
    assert utils.build_consensus_string(profile_matrix) == consensus_string

def test_add_adjacency_list():
    # This test is pretty lame but it does what I want
    text = """>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT"""
    iter_input = text.split("\n")

    expected_set = set(["Rosalind_0498 Rosalind_2391"])
    actual_set = set()
    for gs1 in utils.parse_fasta(iter_input):
        for gs2 in utils.parse_fasta(iter_input):
            # Sets seemed more efficent but I don't think this is. I am just re-creating hashsets again
            # and again. There's probably a better way to handle duplicates in an adjacency list
            # like maybe having a dictionary but I don't really feel like doing that right now
            actual_set = actual_set.union(utils.overlap_graph_connection(gs1, gs2, 3))
    assert expected_set == actual_set