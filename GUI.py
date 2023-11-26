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
      graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = 'Any 12') # labels center of rectangle with base pair
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
    window['analysis_output'].update(f'The GC content for the DNA sequence is {dnaObj.gcContentValue}%\n\nThe transcribed mRNA is\n{dnaObj.mrnaStrand}\n\nThe amino acids for the DNA sequence are:\n{dnaObj.aminoAcids}')
    
    # dna visualization happens here
    
    drawDNA(window['vis'], dnaObj)

# =================END OF HELPER METHODS============================

# ======================GLOBAL VARIABLES============================

# window size
windowSizeX = 800
windowSizeY = 500

# graph size
graphSizeX = windowSizeX - 100
graphSizeY = 90

# =================END OF GLOBAL VARIABLES==========================


# ================GUI LAYOUT========================================
sg.theme('DarkTeal4')

inputLayout = [
  [sg.Text('Enter DNA sequence')],
  [sg.InputText(key = 'sequence_input', size = (50, 10))]
]

builderLayout = [
  [sg.Text('Or build your own DNA sequence!')],
  [sg.Button('A', key = 'a_builder'), sg.Button('T', key = 't_builder'), sg.Button('G', key = 'g_builder'), sg.Button('C', key = 'c_builder'), sg.Button('Delete', key = 'builder_delete')],
  [sg.Button('Randomize', key = 'builder_random')],
  [sg.Button('Analyze')]
]

outputLayout = [
  [sg.Multiline('', size = (50, 10), key = 'analysis_output',expand_x = True, no_scrollbar = True)]
]

visualizationLayer = [
  [sg.Text('Visualization of DNA Strand', justification = 'center')],
  [sg.Graph(canvas_size = (graphSizeX, graphSizeY), key = 'vis', graph_bottom_left = (0,0), graph_top_right = (graphSizeX, graphSizeY))]
]

# layout that smushes it all together
mainLayout = [
  [sg.Column(layout = inputLayout)],
  [sg.Column(layout = builderLayout)],
  [sg.Column(layout = outputLayout)],
  [sg.Column(layout = visualizationLayer)]
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
    randomSequence = 'TAC'
    randomSequence += ''.join(random.choice('ATGC') for i in range(length))
    randomSequence += random.choice(['ATT', 'ATC', 'ACT'])
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

  # does dna analysis
  if event == 'Analyze':
    dnaInput = values['sequence_input']
    
    # this is needed so gui does not crash => not sure why it crashes when text field empty and you click analyze
    if dnaInput == '':
      window['vis'].erase() # clear the graph just because
      continue
    
    # causes error popup and wipes input field if extra letters are in the field
    if not all(base in 'atcgATGC' for base in dnaInput):
      sg.popup_error('Can only use characters \'A T G C\' in DNA sequence. Please try again.', title = 'ERROR!')
      window['sequence_input'].update('')
      continue
    
    analyzeDNA(dnaInput)


window.close()

# ==============END OF EVENT LOOP===================================
