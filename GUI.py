import PySimpleGUI as sg
import dna_analysis_stuff as dna
import random

# ==================HELPER METHODS==================================

# returns color of the base => i got these colors by googling and seeing the first colors that came up
def baseToColor(base):
  colorMap = {'A': 'green', 'T': 'red', 'G': 'yellow', 'C': 'lightblue'}
  return colorMap.get(base, 'white')

# first erases graph and then draws rectangles as base pairs
def drawDNA(graph, dnaObj):
  graph.erase()
  x = 15
  y = 50
  
  baseSpace = 5 # arbitrary number just to space out pairs so they arent hugging
  minWidth = 10
  availibleSpace = graphSizeX - (len(dnaObj.dnaStrand) * baseSpace) - 15
  baseWidth = max(minWidth, availibleSpace / (len(dnaObj.dnaStrand)))
  
  for base in dnaObj.dnaStrand:
    color = baseToColor(base)
    graph.draw_rectangle((x, y), (x + baseWidth, y - 30), line_color = 'black', fill_color = color) # draws colored rectangle
    
    # dont draw letters if bases too small
    if baseWidth > 15: 
      centerX = x + (baseWidth / 2)
      centerY = y - 15
      graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = 'Any 15') # labels center of rectangle with base pair
    x += baseWidth + baseSpace
    
  x = 15
  y = 80
    
  for base in dnaObj.dnaPair:
    color = baseToColor(base)
    graph.draw_rectangle((x, y), (x + baseWidth, y - 30), line_color = 'black', fill_color = color) # draws colored rectangle
    
    # dont draw letters if bases too small
    if baseWidth > 15:
      centerX = x + (baseWidth / 2)
      centerY = y - 15
      graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = 'Any 12') # labels center of rectangle with base pair
    x += baseWidth + baseSpace
  
  # the dna sticks idk what theyre called
  graph.draw_rectangle((10, 20), (15 + (baseWidth + baseSpace) * len(dnaObj.dnaStrand), 15), line_color = 'black', fill_color = 'black')
  graph.draw_rectangle((10, 80), (15 + (baseWidth + baseSpace) * len(dnaObj.dnaStrand), 85), line_color = 'black', fill_color = 'black')
  
def analyzeDNA(dnaInput):
    # creates dna object and then pulls data
    dnaObj = dna.DNA(dnaInput)
    window['dna_desc'].update(visible = True)
    window['analysis_output'].update(f'The GC content for the DNA sequence is {dnaObj.gcContentValue}%\n\nThe transcribed mRNA is\n{dnaObj.mrnaStrand}\n\nThe amino acids for the DNA sequence are:\n{dnaObj.aminoAcids}')
    
    # dna visualization happens here
    
    drawDNA(window['vis'], dnaObj)
    fillTable(dnaObj, int(values['kmer_len_input']))
    
def fillTable(dnaObj, kmerLen):
    values = dna.DNA.createKmerInfo(dnaObj.dnaStrand, kmerLen)
    window['kmer_table'].update(values = values, visible = True)
    window['kmer_desc'].update(visible = True)

# =================END OF HELPER METHODS============================

# ======================GLOBAL VARIABLES============================

# window size
windowSizeX = 800
windowSizeY = 700

# graph size
graphSizeX = windowSizeX - 100
graphSizeY = 90

# =================END OF GLOBAL VARIABLES==========================


# ================GUI LAYOUT========================================
sg.theme('DarkTeal4')

inputLayout = [
  [sg.Text('Enter DNA sequence and kmer length', font = ('Arial Bold', 15))],
  [sg.InputText(key = 'sequence_input', default_text = 'Enter DNA sequence or randomize', size = (50, 10), enable_events = True, font = ('Arial Bold', 12)), sg.InputText(key = 'kmer_len_input', default_text = '3', font = ('Arial Bold', 12), size = (3, 1))]
]

builderLayout = [
  [sg.Text('Or build your own DNA sequence!', font = ('Arial Bold', 15))],
  [sg.Button('A', key = 'a_builder'), sg.Button('T', key = 't_builder'), sg.Button('G', key = 'g_builder'), sg.Button('C', key = 'c_builder'), sg.Button('Delete', key = 'builder_delete'), sg.Button('Clear', key = 'builder_clear')],
  [sg.Button('Randomize', key = 'builder_random')],
  [sg.Button('Analyze')]
]

outputLayout = [
  [sg.Multiline('', size = (50, 10), key = 'analysis_output', font = ('Arial Bold', 15), no_scrollbar = True, disabled = True)]
]

visualizationLayer = [
  [sg.Text('Visualization of DNA Strand', key = 'dna_desc', justification = 'center', font = ('Arial Bold', 15) ,visible = False)],
  [sg.Graph(canvas_size = (graphSizeX, graphSizeY), key = 'vis', graph_bottom_left = (0,0), graph_top_right = (graphSizeX, graphSizeY))]
]

kmerInfoLayer = [
    [sg.Text('Kmer Information', key = 'kmer_desc', font = ('Arial Bold', 20), visible = False)],
    [sg.Table(values = [], key = 'kmer_table', headings = ['  kmers  ', 'reverse complement kmer', ' canonical kmer '],font = ('Arial Bold', 12), visible = False)]
]

# layout that smushes it all together
mainLayout = [
  [sg.Column(layout = inputLayout)],
  [sg.Column(layout = builderLayout)],
  [sg.Column(layout = outputLayout)],
  [sg.Column(layout = visualizationLayer)],
  [sg.Column(layout = kmerInfoLayer)]
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
    randomSequence = 'TAC' # need a start codon for each sequence
    randomSequence += ''.join(random.choice('ATGC') for i in range(length))
    randomSequence += random.choice(['ATT', 'ATC', 'ACT']) #  need a stop codon for each sequence
    window['sequence_input'].update(randomSequence)
    analyzeDNA(randomSequence)
    
  # build your own sequence => i love how well this works
  if event.startswith(('a_builder', 't_builder', 'g_builder', 'c_builder')):
    base = event.split('_')[0]
    updatedText = values['sequence_input'] + base.upper()
    window['sequence_input'].update(updatedText)
    
  # deletes the last input  
  if event == 'builder_delete':
    updatedText = values['sequence_input'][0:len(values['sequence_input']) - 1]
    window['sequence_input'].update(updatedText)
    
  # clears the text box and graph
  if event == 'builder_clear':
      window['sequence_input'].update('')
      window['analysis_output'].update('')
      window['vis'].erase()
      window['dna_desc'].update(visible = False)

  # does dna analysis
  if event == 'Analyze':
    dnaInput = values['sequence_input']
    kmerLen = values['kmer_len_input']
    
    # this is needed so gui does not crash => not sure why it crashes when text field empty and you click analyze
    if dnaInput == '':
      window['vis'].erase() # clear the graph just because
      continue
  
    if not kmerLen.isnumeric():
        sg.popup_error('kmer length must be an integer! Defaulting to 3. Please try again.', title = 'ERROR!', font = ('Arial Bold', 20), keep_on_top = True)
        window['kmer_len_input'].update('3')
        continue
    
    if int(kmerLen) > len(dnaInput):
        sg.popup_error('kmer length cannot be longer than DNA sequence! Please try again.', title = 'ERROR!', font = ('Arial Bold', 20), keep_on_top = True)
        continue
    
    # causes error popup and wipes input field if extra letters are in the field
    if not all(base in 'atcgATGC' for base in dnaInput):
      sg.popup_error('Can only use characters \'A\' \'T\' \'G\' \'C\' in DNA sequence. Please try again.', title = 'ERROR!', font = ('Arial Bold', 20), keep_on_top = True)
      window['sequence_input'].update('')
      continue
    
    analyzeDNA(dnaInput)


window.close()

# ==============END OF EVENT LOOP===================================
