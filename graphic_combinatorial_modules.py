#!/usr/bin/env python
"""
Created on Mon Feb 13 16:07:29 2017

@author: sarahguiziou

To generate the graphical biological design

"""
import dnaplotlib as dpl
import recombinase_graph as rec
import matplotlib.pyplot as plt
from matplotlib import gridspec

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'


# Color maps 
col_map = {}
col_map['red']     = (0.95, 0.30, 0.25)
col_map['green']   = (0.30, 0.75, 0.30)
col_map['blue']    = (0.38, 0.65, 0.87)
col_map['orange']  = (1.00, 0.75, 0.17)
col_map['purple']  = (0.55, 0.35, 0.64)

#col_map_matrix=[(0,0,1), (1,0,0), (0,1,0), (1, 0.75,0), (0.75,0,1), (0,0.75,0.75), (1,1,0), (0.9,0,0.8)]
col_map_matrix=[col_map['blue'], col_map['red'], col_map['green'], col_map['purple'], col_map['orange'], (0,0.75,0.75), (1,1,0), (0.9,0,0.8),(0.5, 0.5, 0) , (0, 0.5, 0.5) , (0.5, 0, 0.5)]

# Function to calculate darker color
def dark (col, fac=1.4):
	return (col[0]/fac, col[1]/fac, col[2]/fac)

"""
BOOLEAN COMPUTATION DEVICE

Generate an image corresponding to the design of the computation device.

"""

def computation_device(inp, strain, int_list, title, name_directory, directory_number):

    # Create the DNAplotlib renderer
    dr = dpl.DNARenderer()

    # Use default renderers and append our custom ones for recombinases
    reg_renderers = dr.std_reg_renderers()
    reg_renderers['Connection'] = rec.flip_arrow
    part_renderers = dr.SBOL_part_renderers()
    part_renderers['RecombinaseSite'] = rec.sbol_recombinase1
    part_renderers['RecombinaseSite2'] = rec.sbol_recombinase2

    # Create the construct programmably to plot
    PF = {'type':'Promoter', 'name':'prom', 'fwd':True}
    out_off = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic', 'label':'GOI'}}
    term = {'type':'Terminator', 'name':'term', 'fwd':True}
    sp={'type':'EmptySpace', 'name':'s', 'fwd':True, 'opts':{'x_extent':1}}

    ##initalization of list
    # list of attB sites
    rec_matrix_attB=[]
    # list of attP sites
    rec_matrix_attP=[]
    # list with the full design
    design=[]
    
    #number of input in this module
    nb_input=int(inp[0])
    
    # for all nb_input add the corresponding attB and attP site in fwd direction to the corresponding list
    for X in range(0,nb_input):
        
        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attB.append(siteB)
        siteP=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attP.append(siteP)

    # number of ZERO correspond to number of total variables (0) - number of ONE (1)
    ZERO=int(inp[0])-int(inp[1])
    
    # variables that correspond to ZERO, in reverse order
    for Z in reversed(range(0, ZERO)):
        # append to the design the attB site
        design.append(rec_matrix_attB[Z])
    # append to the design the promoter in fwd direction        
    design.append(PF)
    # variable which correspond to ZERO        
    for Z in range(0,ZERO):
        # append attP site
        design.append(rec_matrix_attP[Z])
    # variable which correspond to ONE
    for O in range(ZERO, int(inp[0])):
        # add attB site+term+attP site
        design.append(rec_matrix_attB[O])
        design.append(term)
        design.append(rec_matrix_attP[O])
    # append output gene
    design.append(out_off)
    
    ##DESIGN OF INTEGRASE CASSETTE
    
    design_int=[]
    i=0
    name_input='abcdefghijklmnopqrstuvwxyz'
    #name_input='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for integrase in int_list:
        i+=1
        if integrase!='0':
            design_int.append({'type':'Promoter', 'name':'prom', 'fwd':True, 'opts':{'label':'P'+str(name_input[i-1]), 'label_x_offset':-2, 'label_y_offset':-6}})
            int_gene = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map_matrix[int(integrase)-1], 'label_x_offset':-2, 'label_y_offset':-0.5, 'label':'int'+str(integrase),'label_style':'italic'}}
            design_int.append(int_gene)
            design_int.append({'type':'Terminator', 'name':'term', 'fwd':True, 'opts':{'color':(0.5,0.5,0.5)}})
            design_int.append(sp)
            design_int.append(sp)
            design_int.append(sp)
            design_int.append(sp)

    # Create the figure
    fig = plt.figure(figsize=(3,1))
    gs = gridspec.GridSpec(2, 1)
    ax_dna1 = plt.subplot(gs[1])
    ax_dna2 = plt.subplot(gs[0])

    # Redender the DNA to axis
    start, end = dr.renderDNA(ax_dna1, design, part_renderers)
    ax_dna1.set_xlim([start, end])
    ax_dna1.set_ylim([-18,18])
    ax_dna1.set_aspect('equal')
    ax_dna1.set_xticks([])
    ax_dna1.set_yticks([])
    ax_dna1.axis('off')
    
    start, end = dr.renderDNA(ax_dna2, design_int, part_renderers)
    ax_dna2.set_xlim([start, end])
    ax_dna2.set_ylim([-18,18])
    ax_dna2.set_aspect('equal')
    ax_dna2.set_xticks([])
    ax_dna2.set_yticks([])
    ax_dna2.axis('off')
    

    # Update subplot spacing
    plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)
    
    #hfont={'fontname':'Helevetica Neue'}
    
    #plt.title('For strain '+strain+', the computation device is:', fontsize=4)
    
    #plt.suptitle(subtitle, x = 0.5, y= 0.3, fontsize=5, **hfont)
    
    # Save the figure 
    fig.savefig(name_directory+'/'+directory_number+'_'+title+'.png', dpi=300)

    # Clear the plotting cache
    plt.close('all')
    
"""TEST FUNCTION"""

computation_device('10', 'hello', ['0', '1'], 'hello', '../results', 'hello')

""""
LEGEND

Generate the legend of the images
Which site corresponds to which integrase

"""

def legend(nb_input, name_directory, directory_number):


    # Create the DNAplotlib renderer
    dr = dpl.DNARenderer()

    # Use default renderers and append our custom ones for recombinases
    reg_renderers = dr.std_reg_renderers()
    reg_renderers['Connection'] = rec.flip_arrow
    part_renderers = dr.SBOL_part_renderers()
    part_renderers['RecombinaseSite'] = rec.sbol_recombinase1
    part_renderers['RecombinaseSite2'] = rec.sbol_recombinase2

    # list with the full design
    design=[]
    subtitle=''

    # for all nb_input add the corresponding attB and attP site in fwd direction to the corresponding list
    for X in range(0,nb_input):
        
        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        design.append(siteB)
        subtitle+='int'+str(X+1)+'    '
        

    # Create the figure
    fig = plt.figure(figsize=(3,1))
    gs = gridspec.GridSpec(2, 1)
    ax_dna1 = plt.subplot(gs[0])

    # Redender the DNA to axis
    start, end = dr.renderDNA(ax_dna1, design, part_renderers)
    #start, end = dr.renderDNA(design)
    ax_dna1.set_xlim([start, end])
    ax_dna1.set_ylim([-18,18])
    ax_dna1.set_aspect('equal')
    ax_dna1.set_xticks([])
    ax_dna1.set_yticks([])
    ax_dna1.axis('off')

    # Update subplot spacing
    plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)
    
    #plt.title('legend', fontsize=7)
    
    plt.suptitle(subtitle, x = 0.5, y= 0.9, fontsize=6)
   
    # Save the figure 
    fig.savefig(name_directory+'/'+directory_number+'_legend.png', dpi=300)

    # Clear the plotting cache
    plt.close('all')

"""TEST LEGEND"""
    
#legend(2,'results/','007')
