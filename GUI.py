import PySimpleGUI as sg
import dna_analysis_stuff as dna
import random

# ==================HELPER METHODS==================================

# returns color of the base => i got these colors by googling and seeing the first colors that came up
def baseToColor(base):
  colorMap = {'A': 'green', 'T': 'red', 'G': 'yellow', 'C': 'lightblue'}
  return colorMap.get(base, 'white')

def drawDNA(graph, dnaObj):
    graph.erase()
    
    # set x, y values for left side dna pairs
    x = 25
    y = 30

    baseSpace = 5  # arbitrary number just to space out pairs so they aren't hugging
    minWidth = 10
    availableSpace = graphSizeY - (len(dnaObj.dnaStrand) * baseSpace) - 100
    baseHeight = max(minWidth, availableSpace / (len(dnaObj.dnaStrand)))

    # draws the left side dna base pairs
    for base in dnaObj.dnaStrand:
      # get color for each base
      color = baseToColor(base)
      graph.draw_rectangle((x, y), (x + 30, y + baseHeight), line_color = 'black', fill_color = color)
      # only draw letters if base height supports it
      if baseHeight > 15:
        centerX = x + 15
        centerY = y + (baseHeight / 2)
        graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = ('Arial Bold', 12)) # labels center of rectangle with base pair
      y += baseSpace + baseHeight
    
    # set x, y values for right side blocks
    x = 55
    y = 30
    
    # draws the right side dna base pairs
    for base in dnaObj.dnaPair:
      # get color for each base
      color = baseToColor(base)
      graph.draw_rectangle((x, y), (x + 30, y + baseHeight), line_color = 'black', fill_color = color)
      # only draw letters if base height supports it
      if baseHeight > 15:
        centerX = x + 15
        centerY = y + (baseHeight / 2)
        graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = ('Arial Bold', 12)) # labels center of rectangle with base pair
      y += baseSpace + baseHeight
    
    # draws the two sides to the ladder  
    graph.draw_rectangle((20, 25), (25, (baseHeight + baseSpace) * len(dnaObj.dnaStrand) + 30), fill_color = 'black')
    graph.draw_rectangle((85, 25), (90, (baseHeight + baseSpace) * len(dnaObj.dnaStrand) + 30), fill_color = 'black')
      
  
def analyzeDNA(dnaInput, kmerLen):
  # creates dna object and then pulls data
  dnaObj = dna.DNA(dnaInput, kmerLen)
  window['dna_desc'].update(visible = True)
  window['analysis_output'].update(f'The GC content for the DNA sequence is {dnaObj.gcContentValue}%\n\nThe transcribed mRNA is\n{dnaObj.mrnaStrand}\n\nThe amino acids for the DNA sequence are:\n{dnaObj.aminoAcids}')
  
  # dna visualization happens here
  drawDNA(window['vis'], dnaObj)
  # fills out table with info
  fillTable(dnaObj, kmerLen)
  # fills out reading frame info
  fillReadingFrames(dnaObj)
    
    
def fillTable(dnaObj, kmerLen): 
  # calls kmerInfo method to get gui friendly kmers and then updates data and visibility of table
  values = dna.DNA.createKmerInfo(dnaObj.dnaStrand, kmerLen, True)
  window['kmer_table'].update(values = values, visible = True)
  window['kmer_desc'].update(visible = True)
  
  kmerInfo = dna.DNA.generateExtraKmerInfo(dnaObj.kmerInfo, True)

  window['kmer_counts_table'].update(values = kmerInfo[0], visible = True)
  window['kmer_count_desc'].update(visible = True)
  
  window['canon_kmer_counts_table'].update(values = kmerInfo[1], visible = True)
  window['canon_kmer_count_desc'].update(visible = True)
    
def fillReadingFrames(dnaObj):
  readingFrames = dnaObj.readingFrames
  openReadingFrames = dnaObj.openReadingFrames
  
  text = ''  
  text += f'FRAME +1: {" ".join(readingFrames[0])}\n'
  text += f'FRAME +2: {" ".join(readingFrames[1])}\n'
  text += f'FRAME +3: {" ".join(readingFrames[2])}\n'
  text += f'FRAME -1: {" ".join(readingFrames[3])}\n'
  text += f'FRAME -2: {" ".join(readingFrames[4])}\n'
  text += f'FRAME -3: {" ".join(readingFrames[5])}'
    
  window['reading_frames_disp'].update(text, visible = True)
  window['reading_frames_desc'].update(visible = True)
  
  text = ''
  for openReadingFrame in openReadingFrames:
    text += openReadingFrame + '\n'
    
  window['open_reading_frames_disp'].update(text, visible = True)
  window['open_reading_frames_desc'].update(visible = True)  

# =================END OF HELPER METHODS============================

# ======================GLOBAL VARIABLES============================

# window size
windowSizeX = 800
windowSizeY = 1000

# graph size
# graphSizeX = windowSizeX - 100
# graphSizeY = 90

graphSizeX = 100
graphSizeY = windowSizeY - 10

# fonts
headerFont = ('Arial Bold', 15)
textFont = ('Arial Bold', 12)
errorFont = ('Arial Bold', 20)


# =================END OF GLOBAL VARIABLES==========================


# ================GUI LAYOUT========================================
sg.theme('DarkGrey')

# layer with all of the input fields
inputLayout = [
  [sg.Text('Enter DNA sequence and kmer length', font = headerFont)],
  [sg.InputText(key = 'sequence_input', default_text = 'Enter DNA sequence or randomize', size = (55, 10), enable_events = True, font = textFont), sg.InputText(key = 'kmer_len_input', default_text = '3', font = textFont, size = (3, 1))]
]

# layer with all buttons 
builderLayout = [
  [sg.Text('Or build your own DNA sequence!', font = headerFont)],
  [sg.Button('A', key = 'a_builder'), sg.Button('T', key = 't_builder'), sg.Button('G', key = 'g_builder'), sg.Button('C', key = 'c_builder'), sg.Button('Delete', key = 'builder_delete'), sg.Button('Clear', key = 'builder_clear')],
  [sg.Button('Analyze'), sg.Button('Randomize', key = 'builder_random')]
]

# layer with output box
outputLayout = [
  [sg.Multiline('', size = (54, 8), key = 'analysis_output', font = headerFont, no_scrollbar = False, disabled = True)]
]

# layer with image of dna strand (ngl my favorite part)
visualizationLayer = [
  [sg.Text('DNA Strand', key = 'dna_desc', justification = 'center', font = headerFont ,visible = False)],
  [sg.Graph(canvas_size = (graphSizeX, graphSizeY), key = 'vis', graph_bottom_left = (0, graphSizeY), graph_top_right = (graphSizeX, 0))]
]

# tab for list of kmers
kmerTableListTab = [
  [sg.Text('kmer information', key = 'kmer_desc', font = headerFont, visible = False)],
  [sg.Table(values = [], key = 'kmer_table', headings = ['  kmers  ', 'reverse complement kmer', 'canonical kmer'], font = textFont, display_row_numbers = True, alternating_row_color = 'darkgreen', visible = False)]
]

# tab for kmer information
kmerTableCountsTab = [
  [sg.Text('kmer count information', key = 'kmer_count_desc', font = headerFont, visible = False)],
  [sg.Table(values = [], key = 'kmer_counts_table', headings = ['  kmer  ', 'total count', 'distinct count', 'unique count'], font = textFont, display_row_numbers = True, alternating_row_color = 'darkgreen', visible = False)]
]

# tab for canonical kmer information
canonicalKmerTableCountsTab = [
  [sg.Text('canonical kmer count information', key = 'canon_kmer_count_desc', font = headerFont, visible = False)],
  [sg.Table(values = [], key = 'canon_kmer_counts_table', headings = ['  kmer  ', 'total count', 'distinct count', 'unique count'], font = textFont, display_row_numbers = True, alternating_row_color = 'darkgreen', visible = False)]
]

# tab for paragraph of kmer information
kmerInformationTab = [
  [sg.Text('Understanding kmers', font = headerFont)],
  [sg.Multiline('A kmer is a sequence of nucleotides. In order to get a sequence of kmers, you get the first k chars and then slide a \'window\' over by one char, thus creating a overlap with only one new char. You repeat this until the end of the sequence.\n\nSince DNA is double stranded, it is usually beneficial to take a reverse complement kmer which is the same as the complementary sequence but reversed. A canonical kmer is the lexicographically smaller kmer of the standard kmer and reverse complement kmer.\n\nDistinct kmers are only counted once, even if they appear multiple times. Unique kmers are all kmers that appear only once.', font = textFont, size = (65, 11), no_scrollbar = False, disabled = True)]
]

# tabgroup made of tabs
kmerInfoLayer = [
  [sg.TabGroup([
    [sg.Tab('kmers', kmerTableListTab)],
    [sg.Tab('kmer info', kmerTableCountsTab)],
    [sg.Tab('canonical kmer info', canonicalKmerTableCountsTab)],
    [sg.Tab('how to read kmers', kmerInformationTab)]]
               , key = 'kmer_tabs')
  ]
]

# tab for all six reading frames
readingFramesTab = [
  [sg.Text('Six Reading frames for this DNA sequence', font = headerFont, key = 'reading_frames_desc', visible = False)],
  [sg.Multiline('', font = textFont, visible = False, key = 'reading_frames_disp', size = (65, 11), disabled = True)]
]

# tab for only open reading frames
openReadingFramesTab = [
    [sg.Text('Open reading frames for this DNA sequence', font = headerFont, key = 'open_reading_frames_desc', visible = False)],
    [sg.Multiline('', font = textFont, visible = False, key = 'open_reading_frames_disp', size = (65, 11), disabled = True)]
]

readingFramesDescTab = [
    [sg.Text('What are reading frames?', font = headerFont)],
    [sg.Multiline('DNA is made of four base pairs, (A)denine, (T)hymine, (C)ytosine, (G)uanine. These base pairs are then transcribed to mRNA, where the base pairs bond with a sugar and phosphate and T becomes (U)racil. This new combination is known as a nucleotide, and a group of three nucleotides is called a codon which codes for an amino acid.\n\nA reading frame is generated from mRNA, they are created by taking an mRNA sequence and splitting it up into chunks of three and then figuring out the amino acids for each codon. This is done times, each time sliding the start of the first codon over once to the right. In order to get the last three reading frames, the reverse complement sequence is found and then the same process is done.\n\nFor example:\nAUCGCA has frames:\nAUC-GCA\nA-UCG-CA\nAU-CGC-A\nUGC-GAU\nU-GCG-AU\nUG-CGA-U    Note: anything less than length three is not counted as a codon\n\nIf a reading frame has a start codon (Met) and a stop codon (UAA, UAG, UGA) in sequential order, the amino acid sequence between the start and stop codon is an open reading frame. For simplicity\'s sake, the example provided was short and did not show a start and stop codon. Proteins are formed from open reading frames.', font = textFont, size = (65, 11), no_scrollbar = False, disabled = True)]
]

readingFramesLayer = [
  [sg.TabGroup([
    [sg.Tab('reading frames', readingFramesTab)],
    [sg.Tab('open reading frames', openReadingFramesTab)],
    [sg.Tab('reading frames information', readingFramesDescTab)]]
               , key = 'reading_frames_tabs')
  ]
]

leftColumnLayout = [
  [sg.Column(layout = inputLayout)],
  [sg.Column(layout = builderLayout)],
  [sg.Column(layout = outputLayout)],
  [sg.Column(layout = kmerInfoLayer)],
  [sg.Column(layout = readingFramesLayer)]
]


# main layout => lmao its short bc its all stuffed into the left column
mainLayout = [
  [sg.Column(leftColumnLayout, vertical_alignment = 'top'), sg.Column(visualizationLayer, vertical_alignment = 'top')]
]

 # ================END OF GUI LAYOUT================================


# =================EVENT LOOP=======================================

# generates window
window = sg.Window('DNA ANALYSIS', mainLayout, size = (windowSizeX, windowSizeY))

# event loops that continually runs until exit  
while True:
  event, values = window.read()

  # closes the tab 
  if event == sg.WIN_CLOSED or event == 'Exit':
    break
  
  # randomize a sequence of dna
  if event == 'builder_random':
    length = random.randint(20, 34)
    length = (length // 3) * 3
    randomSequence = 'ATG' # need a start codon for each sequence
    randomSequence += ''.join(random.choice('ATGC') for i in range(length))
    randomSequence += random.choice(['TAA', 'TAG', 'TGA']) #  need a stop codon for each sequence
    window['sequence_input'].update(randomSequence)
    kmerLen = values['kmer_len_input']
    
    # if kmer not number create error
    if not kmerLen.isnumeric():
        sg.popup_error('kmer length must be an integer! Defaulting to 3. Please try again.', title = 'ERROR!', font = errorFont, keep_on_top = True)
        window['kmer_len_input'].update('3') # default to 3
        continue
    
    # if kmer longer than dna sequence thats a problem => error
    if int(kmerLen) > len(randomSequence):
        sg.popup_error('kmer length cannot be longer than DNA sequence! Defaulting to 3. Please try again.', title = 'ERROR!', font = errorFont, keep_on_top = True)
        window['kmer_len_input'].update('3') # default to 3
        continue
    
    analyzeDNA(randomSequence, int(kmerLen))
    
  # build your own sequence => i love how well this works
  if event.startswith(('a_builder', 't_builder', 'g_builder', 'c_builder')):
    base = event.split('_')[0]
    updatedText = values['sequence_input'] + base.upper()
    window['sequence_input'].update(updatedText)
    
  # deletes the last input  
  if event == 'builder_delete':
    updatedText = values['sequence_input'][0:len(values['sequence_input']) - 1]
    window['sequence_input'].update(updatedText)
    
  # wipes everything off the screen
  if event == 'builder_clear':
      window['sequence_input'].update('')
      window['analysis_output'].update('')
      window['vis'].erase()
      window['dna_desc'].update(visible = False)
      window['kmer_table'].update(visible = False)
      window['kmer_desc'].update(visible = False)
      window['kmer_counts_table'].update(visible = False)
      window['kmer_count_desc'].update(visible = False)
      window['canon_kmer_counts_table'].update(visible = False)
      window['canon_kmer_count_desc'].update(visible = False)
      window['reading_frames_desc'].update(visible = False)
      window['reading_frames_disp'].update(visible = False)
      window['open_reading_frames_desc'].update(visible = False)
      window['open_reading_frames_disp'].update(visible = False)

  # does dna analysis
  if event == 'Analyze':
    dnaInput = values['sequence_input']
    kmerLen = values['kmer_len_input']
    
    # this is needed so gui does not crash => not sure why it crashes when text field empty and you click analyze
    if dnaInput == '':
      continue
  
    # if kmer not number create error
    if not kmerLen.isnumeric():
        sg.popup_error('kmer length must be an integer! Defaulting to 3. Please try again.', title = 'ERROR!', font = errorFont, keep_on_top = True)
        window['kmer_len_input'].update('3')
        continue
      
    # if kmer longer than dna sequence thats a problem => error
    if int(kmerLen) > len(dnaInput):
        sg.popup_error('kmer length cannot be longer than DNA sequence! Defaulting to 3. Please try again.', title = 'ERROR!', font = errorFont, keep_on_top = True)
        window['kmer_len_input'].update('3')
        continue
    
    # causes error popup and wipes input field if extra letters are in the field
    if not all(base in 'atcgATGC' for base in dnaInput):
      sg.popup_error('Can only use characters \'A\' \'T\' \'G\' \'C\' in DNA sequence. Please try again.', title = 'ERROR!', font = errorFont, keep_on_top = True)
      window['sequence_input'].update('')
      continue
    
    analyzeDNA(dnaInput, int(kmerLen))

window.close()
# bye bye window
# ==============END OF EVENT LOOP===================================