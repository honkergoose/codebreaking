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
        