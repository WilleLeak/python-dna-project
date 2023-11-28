#DNA STUFF
# code starts here
from collections import Counter

class DNA:
    # dict for the amino acid creation
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

    def __init__(self, dnaStrand, kmerLen = 3):
        self.dnaStrand = dnaStrand.upper()
        self.dnaPair = self.createDnaPair(dnaStrand)
        self.mrnaStrand = self.transcribeDnaToRna(dnaStrand)
        self.aminoAcids = self.translateAminoAcids(self.mrnaStrand)
        self.gcContentValue = self.gcContent(dnaStrand)
        self.kmerInfo = self.createKmerInfo(dnaStrand, kmerLen, False)

    # returns the complementary sequence  
    @staticmethod
    def createDnaPair(dnaStrand):
        dnaStrand = dnaStrand.lower()
        dnaStrand = dnaStrand.replace('a', 'T')
        dnaStrand = dnaStrand.replace('t', 'A')
        dnaStrand = dnaStrand.replace('c', 'G')
        dnaStrand = dnaStrand.replace('g', 'C')
        return dnaStrand

    # returns the mrna sequence transcribed from the dna sequence
    @staticmethod
    def transcribeDnaToRna(dnaStrand):
        dnaStrand = dnaStrand.lower()
        dnaStrand = dnaStrand.replace('a', 'U')
        dnaStrand = dnaStrand.replace('t', 'A')
        dnaStrand = dnaStrand.replace('c', 'G')
        dnaStrand = dnaStrand.replace('g', 'C')
        return dnaStrand

    # returns amino acids for the mrna sequence
    @staticmethod
    def translateAminoAcids(mrnaStrand):
        aminoAcids = ''
        for i in range(0, len(mrnaStrand) - 2, 3):
            kmer = mrnaStrand[i:i+3]
            for key, value in DNA.aminoAcidDict.items():
                if kmer in value:
                    aminoAcids += key + ' '
        return aminoAcids.strip()

    # returns gc content as a percent
    @staticmethod
    def gcContent(dnaStrand):
        dnaStrand = dnaStrand.upper()
        return round(((dnaStrand.count('G') + dnaStrand.count('C')) / len(dnaStrand)) * 100, 1)

    # kmers of n length, reverse complement kmer, canonical kmer => returns array with those in order
    @staticmethod
    def createKmerInfo(dnaStrand, kmerLen, forGUI = False):
        kmers = []
        for i in range(len(dnaStrand) - kmerLen + 1):
            kmers.append(dnaStrand[i:i + kmerLen])
        
        reverseComplementKmers = []
        for kmer in kmers:
            reverseComplementKmers.append(DNA.createDnaPair(kmer))
        
        canonicalKmers = []
        for i in range(len(kmers)):
            kmer = kmers[i]
            reverseKmer = reverseComplementKmers[i]
            canonicalKmer = min(kmer, reverseKmer, key = lambda x: x.lower())
            canonicalKmers.append(canonicalKmer)

        values = [kmers, reverseComplementKmers, canonicalKmers]
        guiFriendlyValues = list(map(list, zip(*values)))
        if forGUI: return guiFriendlyValues 
        else: return values

    @staticmethod
    def generateExtraKmerInfo(kmerInfo, forGUI = False):
      kmers = kmerInfo[0]
      kmerTotals = list(dict(Counter(kmers)).values())
      distinctKmers = [1] * len(kmers)
      uniqueKmers = [1 if num == 1 else 0 for num in kmerTotals]
      extraKmerInfo = [kmers, kmerTotals, distinctKmers, uniqueKmers]

      canonKmers = kmerInfo[2]
      kmerTotals = list(dict(Counter(canonKmers)).values())
      distinctKmers = [1] * len(canonKmers)
      uniqueKmers = [1 if num == 1 else 0 for num in kmerTotals]
      extraCanonicalKmerInfo = [kmers, kmerTotals, distinctKmers, uniqueKmers]
      
      if forGUI:
        guiKmerInfo = [list(x) for x in zip(extraKmerInfo[0], extraKmerInfo[1], extraKmerInfo[2], extraKmerInfo[3])]
        guiCanonicalKmerInfo = [list(x) for x in zip(extraCanonicalKmerInfo[0], extraCanonicalKmerInfo[1], extraCanonicalKmerInfo[2], extraCanonicalKmerInfo[3])]
        return [guiKmerInfo, guiCanonicalKmerInfo]
      
      return [extraKmerInfo, extraCanonicalKmerInfo]
    
    # @staticmethod
    # def translateCodon(self, codon):
    #   for aminoAcid, codons in self.aminoAcidDict.items():
    #       print(codon)
    #       if codon in codons:
    #         return aminoAcid          
            
      
#   def translateReadingFrames(self):
#     mrnaStrand = self.mrnaStrand
#     frames = []
#     for frame in range(3):
#       proteinSequence = ''
#       for i in range(frame, len(mrnaStrand), 3):
#         codon = mrnaStrand[i:i+3]
#         proteinSequence += self.translateCodon(codon)
#       frames.append(proteinSequence)
      
#     reverseComplement = self.translateDnaToRna(self.dnaPair).reverse()
#     print(reverseComplement)
      
#     for frame in range(3):
#       proteinSequence = ''
#       for i in range(frame, len(mrnaStrand), 3):
#         codon = ''
  

  


# tests

# d = DNA('cgcaaacg', 3)
# # print(d.dnaStrand)
# # print(d.dnaPair) 
# # print(d.mrnaStrand)
# # print(d.aminoAcids)
# print(d.gcContentValue)
# # print(DNA.createKmerInfo(d.dnaStrand, 3))
# #print(DNA.generateExtraKmerInfo(d.kmerInfo, True))
