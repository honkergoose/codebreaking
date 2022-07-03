# A module of utilities for solving ROSALIND problems


def find_rp(dna, mi, mx):
  if len(dna) < mi:
    raise "hell"

  res = []

  lw = 0
  hi = mi
  
  while lw + mi < len(dna): #  | 0 4 | 0 5 | 0 6 | ... | 0 12 | 1 5 | ...
    while hi - lw <= mx and hi <= len(dna):
      print(dna[lw:hi])
      print(rc(dna[lw:hi]) + " -- ")
      if rp(dna[lw:hi]):
        res.append([lw+1, hi-lw])
      hi += 1
    lw += 1
    hi = lw + mi
  return res
  

def rp(dna):
  return dna == rc(dna)

def rc(dna):
  return "".join(list(reversed(compliment(dna))))

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

if __name__ == "__main__":
  dna = "TCAATGCATGCGGGTCTATATGCAT"
  print(find_rp(dna, 4, 12))
  print("---")
  print(dna)
