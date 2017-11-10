# -*- coding: utf-8 -*-

#!/usr/bin/env python
"""
	Recombinase NOT-gate - PART OF THE CODE
"""

import math
from matplotlib.patheffects import Stroke
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch

__author__  = 'Bryan Der <bder@mit.edu>, Voigt Lab, MIT\n\
               Thomas Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT\n\
               Sarah Guiziou <guiziou@cbs.cnrs.fr, Bonnet lab, MIT'
__license__ = 'MIT'
__version__ = '2.0'


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
