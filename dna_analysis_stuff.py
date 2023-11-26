#DNA STUFF
# code starts here

class DNA:
  aminoAcidDict = {'Phe': ['UUU', 'UUC'],
                   'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                   'Ile': ['AUU', 'AUC', 'AUA'],
                   'Met': ['AUG'],
                   'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
                   'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                   'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
                   'Thr': ['ACU', 'ACC', 'ACA', 'ACG'],
                   'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
                   'Tyr': ['UAU', 'UAC'],
                   'Stop': ['UAA', 'UAG', 'UGA'],
                   'His': ['CAU', 'CAC'],
                   'Gln': ['CAA', 'CAG'],
                   'Asn': ['AAU', 'AAC'],
                   'Lys': ['AAA', 'AAG'],
                   'Asp': ['GAU', 'GAC'],
                   'Glu': ['GAA', 'GAG'],
                   'Cys': ['UGU', 'UGC'],
                   'Trp': ['UGG'],
                   'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                   'Gly': ['GGU', 'GGC', 'GGA', 'GGG']
                   }

  def __init__(self, dnaStrand):
    self.dnaStrand = dnaStrand.upper()
    self.dnaPair = self.createDnaPair(dnaStrand)
    self.mrnaStrand = self.transcribeDnaToRna(dnaStrand)
    self.aminoAcids = self.translateAminoAcids(self.mrnaStrand)
    self.gcContentValue = self.gcContent(dnaStrand)
    
  @staticmethod
  def createDnaPair(dnaStrand):
    dnaStrand = dnaStrand.lower()
    dnaStrand = dnaStrand.replace('a', 'T')
    dnaStrand = dnaStrand.replace('t', 'A')
    dnaStrand = dnaStrand.replace('c', 'G')
    dnaStrand = dnaStrand.replace('g', 'C')
    return dnaStrand

  @staticmethod
  def transcribeDnaToRna(dnaStrand):
    dnaStrand = dnaStrand.lower()
    dnaStrand = dnaStrand.replace('a', 'U')
    dnaStrand = dnaStrand.replace('t', 'A')
    dnaStrand = dnaStrand.replace('c', 'G')
    dnaStrand = dnaStrand.replace('g', 'C')
    return dnaStrand

  @staticmethod
  def translateAminoAcids(mrnaStrand):
    aminoAcids = ''
    for i in range(0, len(mrnaStrand) - 2, 3):
      kmer = mrnaStrand[i:i+3]
      for key, value in DNA.aminoAcidDict.items():
        if kmer in value:
          aminoAcids += key + ' '
    return aminoAcids.strip()

  @staticmethod
  def gcContent(dnaStrand):
    return round((dnaStrand.count('G') + dnaStrand.count('C')) / len(dnaStrand), 1) * 100
  

  


# tests
# d = DNA('TACGCATTAATT')
# print(d.dnaStrand)
# print(d.dnaPair)
# print(d.mrnaStrand)
# print(d.aminoAcids)
# print(d.gcContentValue)
