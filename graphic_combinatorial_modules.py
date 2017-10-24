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
col_map_matrix=[col_map['blue'], col_map['red'], col_map['green'], col_map['purple'], col_map['orange'], (0,0.75,0.75), (1,1,0), (0.9,0,0.8)]

# Function to calculate darker colour
def dark (col, fac=1.4):
	return (col[0]/fac, col[1]/fac, col[2]/fac)


def computation_device(inp, strain, subtitle, title, name_directory, directory_number):


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
    sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':1}}
    PF = {'type':'Promoter', 'name':'prom', 'fwd':True}
    PR={'type':'Promoter', 'name':'prom', 'fwd':False}
    out_on = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['green'], 'label':'OUTPUT', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
    out_off = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':(0.8,0.8,0.8), 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
    term = {'type':'Terminator', 'name':'term', 'fwd':True}

    ##initalization of list
    # list of attB sites
    rec_matrix_attB=[]
    # list of attP sites
    rec_matrix_attP=[]
    # list with the full design
    design=[]
    
    ## if we have a design with a xor or nxor
    if ('R' in inp) or ('N' in inp):
        #correspond to the number of input in this module, 1-nb of input in xor/nxor
        # 2- number in the other module
        nb_input=int(inp[1])+int(inp[2])
    
    else:
        #number of input in this module
        nb_input=int(inp[0])
    
    # for all nb_input add the corresponding attB and attP site in fwd direction to the corresponding list
    for X in range(0,nb_input):
        
        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attB.append(siteB)
        siteP=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        rec_matrix_attP.append(siteP)

    # if xor or nxor in the construct    
    if ('R' in inp) or ('N' in inp):
    
        # number of variable in the xor or nxor
        nb_xor=int(inp[1])
        # number of variable otherwise
        nb_other=int(inp[2])
        
        # for variable in the xor or nxor
        for S in range(0, nb_xor):
            # add to the design the attB site
            design.append(rec_matrix_attB[S])
        # if xor, add a reverse promoter
        if 'R' in inp:
            design.append(PR)
        # if nxor add a fwd promoter
        else:
            design.append(PF)
        #for variable in the xor or nxor add the attP site in reverse orientation in the reverse order of variables
        for S in reversed(range(0, nb_xor)):
            design.append({'type':'RecombinaseSite',  'name':'Site1',  'fwd':False,   'opts':{'color2':col_map_matrix[S], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        
        # number of ZERO in the second part of the construct
        # inp[3] correspond to the number of ONE in the construct
        ZERO=nb_other-int(inp[3])
        # variable which correspond to ZERO   
        for Z in range(nb_xor, ZERO+nb_xor):
            #add at the beginning of the design the attB site
            design.insert(0, rec_matrix_attB[Z])
        # variable which correspond to ZERO        
        for Z in range(nb_xor, ZERO+nb_xor):
            # add at the end of the design the attP site
            design.append(rec_matrix_attP[Z])
        # variable which correspond to ONE
        for O in range(ZERO+nb_xor, nb_xor+nb_other):
            # add attB site+term+attP site
            design.append(rec_matrix_attB[O])
            design.append(term)
            design.append(rec_matrix_attP[O])
        
    ## if no xor neither nxor   
    else:
        # number of ZERO correspond to number of total variables (0) - number of ONE (1)
        ZERO=int(inp[0])-int(inp[1])
        
        # variable which correspond to ZERO, in reverse order
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
    
    plt.title('For strain '+strain+', the computation device is:', fontsize=7)
    
    plt.suptitle(subtitle, x = 0.5, y= 0.3, fontsize=5)
    
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
    #fig.savefig(name_directory+'/'+title+'.pdf', transparent=True)
    fig.savefig(name_directory+'/'+directory_number+'_'+title+'.png', dpi=300)

    # Clear the plotting cache
    plt.close('all')
    

#computation_device('70', 'hello', 'Var2 - int1, Var3 - int2, Var5 - int3','hello')


def legend(nb_input, name_directory, directory_number):


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

    # list with the full design
    design=[]
    subtitle=''

    # for all nb_input add the corresponding attB and attP site in fwd direction to the corresponding list
    for X in range(0,nb_input):
        
        siteB=({'type':'RecombinaseSite',  'name':'Site1',  'fwd':True,   'opts':{'color':col_map_matrix[X], 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}})
        design.append(siteB)
        subtitle+='Integrase'+str(X+1)+'   '
        

    # Create the figure
    fig = plt.figure(figsize=(3,1))
    gs = gridspec.GridSpec(3, 1)
    ax_dna1 = plt.subplot(gs[1])

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
    
    plt.title('legend', fontsize=7)
    
    plt.suptitle(subtitle, x = 0.5, y= 0.3, fontsize=5)
    
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
    #fig.savefig(name_directory+'/'+title+'.pdf', transparent=True)
    fig.savefig(name_directory+'/'+directory_number+'_legend.png', dpi=300)

    # Clear the plotting cache
    plt.close('all')
    
#legend(2,'results/','007')
