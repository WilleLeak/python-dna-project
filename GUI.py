#GUI

import PySimpleGUI as sg
import dna_analysis_stuff as dna

layout = [
    [sg.Text('Enter DNA sequence')],
    [sg.InputText(key = 'sequence_input')],
    [sg.Button('Analyze')],
    [sg.Output(size = (50, 10), key = 'analysis_output')]
]

window = sg.Window('DNA ANALYSIS GUI', layout, size = (500, 500))

while True:
  event, values = window.read()

  if event == sg.WIN_CLOSED or event == 'Exit':
    break

  if event == 'Analyze':
    dna_input = values['sequence_input']
    dna_object_instance = dna.DNA(dna_input)
    
    print(f'the mrna strand is {dna_object_instance.mrnaStand}')

window.close()

