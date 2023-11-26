#GUI

import PySimpleGUI as sg
import dna_analysis_stuff as dna
import random
import traceback

# ==================HELPER METHODS==================================

# returns color of the base => i got these colors by googling and seeing the first colors that came up
def baseToColor(base):
  colorMap = {'A': 'green', 'T': 'red', 'G': 'Yellow', 'C': 'blue'}
  return colorMap.get(base, 'white')

# first erases graph and then draws rectangles as base pairs
def drawDNA(graph, dnaObj):
  graph.erase()
  x = 15
  y = 50
  
  minWidth = 10
  baseWidth = max(minWidth, 380 / (len(dnaObj.dnaStrand)))
  baseSpace = 5 # arbitrary number just to space out pairs so they arent hugging
  
  for base in dnaObj.dnaStrand:
    color = baseToColor(base)
    graph.draw_rectangle((x, y), (x + baseWidth, y - 30), line_color = 'black', fill_color = color) # draws colored rectangle
    centerX = x + (baseWidth / 2)
    centerY = y - 15
    graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = 'Any 12') # labels center of rectangle with base pair
    x += baseWidth + baseSpace
    
  x = 15
  y = 80
    
  for base in dnaObj.dnaPair:
    color = baseToColor(base)
    graph.draw_rectangle((x, y), (x + baseWidth, y - 30), line_color = 'black', fill_color = color) # draws colored rectangle
    centerX = x + (baseWidth / 2)
    centerY = y - 15
    graph.draw_text(text = base, location = (centerX, centerY), color = 'black', font = 'Any 12') # labels center of rectangle with base pair
    x += baseWidth + baseSpace
  
  # the dna sticks idk what theyre called
  graph.draw_rectangle((10, 20), (15 + (baseWidth + baseSpace) * len(dnaObj.dnaStrand), 15), line_color = 'black', fill_color = 'black')
  graph.draw_rectangle((10, 80), (15 + (baseWidth + baseSpace) * len(dnaObj.dnaStrand), 85), line_color = 'black', fill_color = 'black')

# =================END OF HELPER METHODS============================

sg.theme('DarkTeal9')

inputLayout = [
  [sg.Text('Enter DNA sequence')],
  [sg.InputText(key = 'sequence_input')]
]

builderLayout = [
  [sg.Text('Or build your own DNA sequence!')],
  [sg.Button('A', key = 'a_builder'), sg.Button('T', key = 't_builder'), sg.Button('G', key = 'g_builder'), sg.Button('C', key = 'c_builder'), sg.Button('Delete', key = 'builder_delete')],
  [sg.Button('Randomize', key = 'builder_random')],
  [sg.Button('Analyze')]
]

outputLayout = [
  [sg.Output(size = (50, 10), key = 'analysis_output')]
]

visualizationLayer = [
  [sg.Text('Visualization of DNA sequence', justification = 'center')],
  [sg.Graph(canvas_size = (400, 100), key = 'vis', graph_bottom_left = (0,0), graph_top_right = (500, 90))]
]

# layout that smushes it all together
mainLayout = [
  [sg.Column(layout = inputLayout)],
  [sg.Column(layout = builderLayout)],
  [sg.Column(layout = outputLayout)],
  [sg.Column(layout = visualizationLayer)]
]

# generates window
window = sg.Window('DNA ANALYSIS', mainLayout, size = (500, 500))

# event loops that continually runs until exit  
while True:
  event, values = window.read()

  # closes the tab 
  if event == sg.WIN_CLOSED or event == 'Exit':
    break
  
  # randomize a sequence of dna
  if event == 'builder_random':
    length = random.randint(10, 20)
    randomSequence = ''.join(random.choice('ATGC') for i in range(length))
    window['sequence_input'].update(randomSequence)
    
  # build your own sequence => i love how well this works
  if event.startswith(('a_builder', 't_builder', 'g_builder', 'c_builder')):
    letter = event.split('_')[0]
    updatedText = values['sequence_input'] + letter.upper()
    window['sequence_input'].update(updatedText)
    
  # deletes the last input  
  if event == 'builder_delete':
    updatedText = values['sequence_input'][0:len(values['sequence_input']) - 1]
    window['sequence_input'].update(updatedText)

  # does dna analysis
  if event == 'Analyze':
    dna_input = values['sequence_input']
    
    # this is needed so gui does not crash => not sure why it crashes when text field empty and you click analyze
    if dna_input == '':
      window['vis'].erase() # clear the graph just because
      continue
    
    # causes error popup and wipes input field if extra letters are in the field
    if not all(base in 'atcgATGC' for base in dna_input):
      sg.popup_error('Can only use characters \'A T G C\' in DNA sequence. Please try again.', title = 'ERROR!')
      window['sequence_input'].update('')
      continue
    
    # creates dna object and then pulls data
    dnaObj = dna.DNA(dna_input)
    window['analysis_output'].update(f'The GC content for the DNA sequence is {dnaObj.gcContentValue}%\n\nThe transcribed mRNA is\n{dnaObj.mrnaStrand}\n\nThe amino acids for the DNA sequence are:\n{dnaObj.aminoAcids}')
    
    # dna visualization happens here
    try:
      drawDNA(window['vis'], dnaObj)
    except Exception as e:
      sg.popup_error(f'An error occurred: {str(e)}', title='ERROR!')
      traceback.print_exc()

window.close()

