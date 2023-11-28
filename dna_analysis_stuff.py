#DNA STUFF
# code starts here

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

    def __init__(self, dnaStrand, kmerLen):
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
        return round((dnaStrand.count('G') + dnaStrand.count('C')) / len(dnaStrand), 1) * 100

    @staticmethod
    def translateCodon(self, codon):
        for aminoAcid, codons in self.aminoAcidDict.items():
            print(codon)
            if codon in codons:
                return aminoAcid

    # kmers of n length, reverse complement kmer, canonical kmer => returns array with those in order
    @staticmethod
    def createKmerInfo(dnaStrand, kmerLen, forGUI):
        kmers = []
        for i in range(len(dnaStrand) - kmerLen):
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

# d = DNA('TACGCATTAATT')
# print(d.dnaStrand)
# print(d.dnaPair) 
# print(d.mrnaStrand)
# print(d.aminoAcids)
# print(d.gcContentValue)
# print(DNA.createKmerInfo(d.dnaStrand, 3))
