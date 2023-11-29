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
        self.extraKmerInfo = self.generateExtraKmerInfo(self.kmerInfo, False)
        self.readingFrames = self.findReadingFrames(self.mrnaStrand)
        self.openReadingFrames = self.getOpenReadingFrames(self.readingFrames)

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
        dnaStrand = dnaStrand.replace('T', 'U')
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
            reverseComplementKmers.append(DNA.createDnaPair(kmer)[::-1])
        
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
    
    # loops through the codons and amino acid dict to create a list of amino acids
    @staticmethod
    def translateCodons(codons):
      proteinSequence = []
      for codon in codons:
        for aminoAcid, codonList in DNA.aminoAcidDict.items():
            if codon in codonList:
              proteinSequence.append(aminoAcid)
              break # save time by not running through whole dict
      return proteinSequence
    
    # returns all six reading frames, 3 that are forwards and 3 that are backwards        
    @staticmethod  
    def findReadingFrames(mrnaStrand):
      proteinSequences = []
      
      # forward frames
      for readingFrame in range(3):
        codons = [mrnaStrand[i:i + 3] for i in range(readingFrame, len(mrnaStrand), 3)] # skips 0, 1, 2 chars for each of the three forward frames
        proteinSequence = DNA.translateCodons(codons)
        proteinSequences.append(proteinSequence)
      
      # backwards frames  
      mrnaStrand = mrnaStrand.lower()
      mrnaStrand = mrnaStrand.replace('a', 'U')
      mrnaStrand = mrnaStrand.replace('u', 'A')
      mrnaStrand = mrnaStrand.replace('c', 'G')
      mrnaStrand = mrnaStrand.replace('g', 'C')
      # reverse the string bc last three reading frames backwards for some reason
      mrnaStrand = mrnaStrand[::-1] 
      
      for readingFrame in range(3):
        codons = [mrnaStrand[i:i + 3] for i in range(readingFrame, len(mrnaStrand), 3)] # skips 0, 1, 2 chars for each of the three forward frames
        proteinSequence = DNA.translateCodons(codons)
        proteinSequences.append(proteinSequence)
        
      return proteinSequences
        
    # this method is a mouthful => finds reading frames that have Met and Stop in sequential order 
    # and returns a list of strings that are amino acids from Met to Stop (Met and Stop cant touch or so I think)    
    @staticmethod  
    def getOpenReadingFrames(readingFrames):
      listOfOpenReadingFrames = []
      for proteinSequence in readingFrames:
        # met and stop cant touch but stop has to be after met
        if 'Met' in proteinSequence and 'Stop' in proteinSequence[proteinSequence.index('Met') + 2:]:
          openReadingFrame = []
          # met is the start index for adding to the list
          startIndex = proteinSequence.index('Met') 
          while startIndex < len(proteinSequence) and proteinSequence[startIndex] != 'Stop':
            openReadingFrame.append(proteinSequence[startIndex])
            startIndex += 1
          # add stop bc the loop stops at stop
          openReadingFrame.append('Stop')
          # append to grand list
          listOfOpenReadingFrames.append(openReadingFrame)
      
      # this is all for gui => it makes it a single list of strings so its easier for me to display     
      openReadingFrames = []
      for openReadingFrame in listOfOpenReadingFrames:
        openReadingFrames.append(' '.join(openReadingFrame).strip())
          
      return openReadingFrames

  

  


# tests

# d = DNA('CGCTACGTCTTACGCTGGAGCTCTCATGGATCGGTTCGGTAGGGCTCGATCACATCGCTAGCCAT', 3)
# # print(d.dnaStrand)
# # print(d.dnaPair) 
# # print(d.mrnaStrand)
# # print(d.aminoAcids)
# # print(d.gcContentValue)
# # print(DNA.createKmerInfo(d.dnaStrand, 3))
# # print(DNA.generateExtraKmerInfo(d.kmerInfo, True))
# #print(DNA.translateReadingFrames(d.mrnaStrand))
# print(d.openReadingFrames)
