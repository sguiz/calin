#!/usr/bin/env python
"""
	Recombinase NOT-gate
"""

import math
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches

__author__  = 'Bryan Der <bder@mit.edu>, Voigt Lab, MIT\n\
               Thomas Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'MIT'
__version__ = '1.0'

def sbol_recombinase1 (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	""" SBOL recombinase site renderer - forward direction
	"""
	# Default parameters
	color = (0,0,0)
	color2 = (0,0,0)
	start_pad = 0.0
	end_pad = 0.0
	x_extent = 6.0
	y_extent = 6.0
	linestyle = '-'
	# Update default parameters if provided
	if opts != None:
		if 'start_pad' in list(opts.keys()):
			start_pad = opts['start_pad']
		if 'end_pad' in list(opts.keys()):
			end_pad = opts['end_pad']
		if 'x_extent' in list(opts.keys()):
			x_extent = opts['x_extent']
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'linestyle' in list(opts.keys()):
			linestyle = opts['linestyle']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'scale' in list(opts.keys()):
			scale = opts['scale']
		if 'color' in list(opts.keys()):
			color = opts['color']
		if 'color2' in list(opts.keys()):
			color2 = opts['color2']
	# Check direction add start padding
	final_end = end
	final_start = prev_end
	y_lower = -1 * y_extent/2
	y_upper = y_extent/2
	if start > end:
		start = prev_end+end_pad+x_extent+linewidth
		end = prev_end+end_pad
		final_end = start+start_pad
		color = color2
	else:
		start = prev_end+start_pad+linewidth
		end = start+x_extent
		final_end = end+end_pad
	# Draw the site
	p1 = Polygon([(start, y_lower), 
		          (start, y_upper),
		          (end,0)],
		          edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
		          path_effects=[Stroke(joinstyle="miter")])		
	ax.add_patch(p1)
	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_recombinase2 (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	""" SBOL recombinase site renderer - reverse direction
	"""
	# Default parameters
	color = (0,0,0)
	color2 = (0,0,0)
	start_pad = 0.0
	end_pad = 0.0
	x_extent = 6.0
	y_extent = 6.0
	linestyle = '-'
	# Update default parameters if provided
	if opts != None:
		if 'start_pad' in list(opts.keys()):
			start_pad = opts['start_pad']
		if 'end_pad' in list(opts.keys()):
			end_pad = opts['end_pad']
		if 'x_extent' in list(opts.keys()):
			x_extent = opts['x_extent']
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'linestyle' in list(opts.keys()):
			linestyle = opts['linestyle']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'scale' in list(opts.keys()):
			scale = opts['scale']
		if 'color' in list(opts.keys()):
			color = opts['color']
		if 'color2' in list(opts.keys()):
			color2 = opts['color2']
		else:
			if 'color' in list(opts.keys()):
				r2 = float(color[0]) / 2
				g2 = float(color[1]) / 2
				b2 = float(color[2]) / 2
				color2 = (r2,g2,b2)
	# Check direction add start padding
	final_end = end
	final_start = prev_end
	y_lower = -1 * y_extent/2
	y_upper = y_extent/2
	if start > end:
		start = prev_end+end_pad+x_extent+linewidth
		end = prev_end+end_pad
		final_end = start+start_pad
		temp = color
		color = color2
		color2 = temp
	else:
		start = prev_end+start_pad+linewidth
		end = start+x_extent
		final_end = end+end_pad
	# Draw the site
	p1 = Polygon([(start, y_lower), 
		         (start, y_upper),
		          (end,0)],
		          edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
		          path_effects=[Stroke(joinstyle="miter")]) 
	midpoint = (end + start) / 2
	hypotenuse = math.sqrt( (y_extent/2)**2 + (x_extent)**2 )
	hypotenuse2 = hypotenuse / 2
	cosineA = (y_extent/2) / hypotenuse
	f = hypotenuse2 * cosineA
	p2 = Polygon([(midpoint, -1*f), 
		          (midpoint, f),
		          (end,0)],
		          edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=12, 
		          path_effects=[Stroke(joinstyle="miter")]) 	
	ax.add_patch(p1)
	ax.add_patch(p2)	
	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def flip_arrow (ax, type, num, from_part, to_part, scale, linewidth, arc_height_index, opts):
	""" Regulation arcs for recombinase sites
	"""
	# Default parameters
	color = (0.0,0.0,0.0)
	arcHeightStart = 10
	arcHeightEnd = 10
	# Update default parameters if provided
	if opts != None:
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'color' in list(opts.keys()):
			color = opts['color']
		if 'arc_height_start' in list(opts.keys()):
			arcHeightStart = opts['arc_height_start']
		if 'arc_height_end' in list(opts.keys()):
			arcHeightEnd = opts['arc_height_end']
	start = (from_part['start'] + from_part['end']) / 2
	end   = (to_part['start']   + to_part['end']) / 2
	# Check direction and draw arc
	if start > end:
		arcHeightStart = -arcHeightStart
		arcHeightEnd = -arcHeightEnd
	ax.annotate('', (end, arcHeightEnd), (start, arcHeightStart), ha="right", va="center", size=8, arrowprops=dict(arrowstyle='->',connectionstyle="arc3,rad=-.4",lw=linewidth, color=color))

# Color maps 
col_map = {}
col_map['red']     = (0.95, 0.30, 0.25)
col_map['green']   = (0.38, 0.82, 0.32)
col_map['blue']    = (0.38, 0.65, 0.87)
col_map['orange']  = (1.00, 0.75, 0.17)
col_map['purple']  = (0.55, 0.35, 0.64)

#col_map_matrix=[(0,0,1), (1,0,0), (0,1,0), (1, 0.75,0), (0.75,0,1), (0,0.75,0.75), (1,1,0), (0.9,0,0.8)]
col_map_matrix=[col_map['blue'], col_map['red'], col_map['green'], col_map['purple'], col_map['orange']]

# Function to calculate darker colour
def dark (col, fac=1.4):
	return (col[0]/fac, col[1]/fac, col[2]/fac)



def main(nb_input, list_integrase, list_gene_state, title, name_directory, dirc_number):


    # Global line width
    lw = 1.0

    # Create the DNAplotlib renderer
    dr = dpl.DNARenderer()

    # Use default renderers and append our custom ones for recombinases
    reg_renderers = dr.std_reg_renderers()
    reg_renderers['Connection'] = flip_arrow
    part_renderers = dr.SBOL_part_renderers()
    part_renderers['RecombinaseSite'] = sbol_recombinase1
    part_renderers['RecombinaseSite2'] = sbol_recombinase2

    # Create the construct programmably to plot
    PF = {'type':'Promoter', 'name':'prom', 'fwd':True}
    PR={'type':'Promoter', 'name':'prom', 'fwd':False}
    out_fwd = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.80,0.8,0.8), 'label':' ', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
    out_off = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.4,0.4,0.4), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
    out_rev = {'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':' ', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
    term = {'type':'Terminator', 'name':'term', 'fwd':True}

    ##initalization of list
    # list of attB fwd sites
    rec_matrix_attB_fwd=[]
    #list of attP fwd sites
    rec_matrix_attP_fwd=[]
    # list of attB rev sites
    rec_matrix_attB_rev=[]
    # list of attP fwd sites
    rec_matrix_attP_rev=[]
    # list with the full design
    design=[]
    
    # for all nb_input add the corresponding attB and attP site in fwd direction to the corresponding list
    for X in range(0,nb_input):
        
        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[list_integrase[X]-1], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attB_fwd.append(siteB)
        
        siteP=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[list_integrase[X]-1], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attP_fwd.append(siteP)

        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':False,   'opts':{'color2':col_map_matrix[list_integrase[X]-1], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attB_rev.append(siteB)
        
        siteP=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':False,   'opts':{'color2':col_map_matrix[list_integrase[X]-1], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attP_rev.append(siteP)

    # construction for 2 inputs
    if nb_input==2:
        
        if list_gene_state[2]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[2]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        
        design.append(rec_matrix_attB_fwd[1])
        design.append(PF)
        design.append(rec_matrix_attB_fwd[0])
        
        if list_gene_state[0]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[0]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))

        if list_gene_state[1]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[1]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
            
        design.append(rec_matrix_attP_fwd[1])
        design.append(rec_matrix_attP_rev[0])

    # construction for 3 inputs
    elif nb_input==3:
        
        if list_gene_state[3]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[3]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        
        design.append(rec_matrix_attB_fwd[2])
        
        if list_gene_state[0]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[0]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
            
        design.append(rec_matrix_attB_fwd[0])
        design.append(rec_matrix_attB_rev[1])
        design.append(PR)
        design.append(rec_matrix_attP_rev[0])
        
        if list_gene_state[1]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[1]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        if list_gene_state[2]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[2]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
            
        design.append(rec_matrix_attP_fwd[2])
        design.append(rec_matrix_attP_rev[1])

    # construction for 4 inputs
    elif nb_input==4:

        if list_gene_state[4]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[4]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        
        design.append(rec_matrix_attB_fwd[3])
        design.append(rec_matrix_attB_fwd[0])
        design.append(rec_matrix_attB_rev[1])

        if list_gene_state[1]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[1]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        if list_gene_state[0]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[0]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
            
        design.append(rec_matrix_attP_rev[0])
        design.append(rec_matrix_attB_rev[2])
        design.append(PR)
        design.append(rec_matrix_attP_rev[1])
        
        if list_gene_state[2]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[2]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
        if list_gene_state[3]!=0:
            design.append(({'type':'CDS', 'name':'cds', 'fwd':False, 'opts':{'color':(0.8,0.8,0.8), 'label':str(list_gene_state[3]), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}))
            
        design.append(rec_matrix_attP_fwd[3])
        design.append(rec_matrix_attP_rev[2])


    # Create the figure
    fig = plt.figure(figsize=(3,1))
    gs = gridspec.GridSpec(3, 1)
    ax_dna1 = plt.subplot(gs[1])

    # Redender the DNA to axis
    start, end = dr.renderDNA(ax_dna1, design, part_renderers)
    ax_dna1.set_xlim([start, end])
    ax_dna1.set_ylim([-18,18])
    ax_dna1.set_aspect('equal')
    ax_dna1.set_xticks([])
    ax_dna1.set_yticks([])
    ax_dna1.axis('off')

    # Update subplot spacing
    plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)
    
    plt.title(title, fontsize=7)
    
#    plt.suptitle(subtitle, x = 0.5, y= 0.3, fontsize=5)
    
#        
#    rect1 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[0])
#    rect2 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[1])
#    rect3 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[2])
#    rect4 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[3])
#    rect5 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[4])
#    rect6 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[5])
#    rect7 = Line2D([], [], marker="s", markersize=2.5, linewidth=0, color=col_map_matrix[6])
#
#    plt.legend((rect1, rect2, rect3, rect4, rect5, rect6, rect7), ('Int1', 'Int2', 'Int3', 'Int4', 'Int5', 'Int6', 'Int7'), fontsize=3, loc='lower right')
#    
    # Save the figure 
    #fig.savefig(name_directory+'/'+name_directory+'_'+title+'.pdf', transparent=True)
    fig.savefig(name_directory+'/'+dirc_number+'_'+title+'.png', dpi=300)

    # Clear the plotting cache
    plt.close('all')
    
########TEST#########
#main(2, (2,1), [0,2,0], 'hello','../results/000344','000344')
