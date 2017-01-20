#!/usr/bin/env python
"""
Generate a heatmap from a pairwise comparison matrix

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

import sys
import numpy as np
import toyplot
import toyplot.browser
#import toyplot.png

HEADER_CHAR = '#'

def main(matrix_filename, width, height):
    # Get matrix labels and number of cols in file
    matrix_file = open(matrix_filename, 'r')
    llabels = []
    nrows = -1
    for line in matrix_file:
        nrows += 1
        # Header
        if nrows == 0:
            tlabels = line.split()
            ncols = len(tlabels)
        # Data
        else:
            llabels.append(line.split()[0])
    matrix_file.close()
    
    # Read file into numpy matrix
    lxp = np.loadtxt(matrix_filename, skiprows = 1, 
                     usecols=range(1, ncols + 1))

    # Plot heatmap
    print(lxp.max())
    canvas = toyplot.Canvas(width = width, height = height)
    colormap = toyplot.color.LinearMap(toyplot.color.brewer("Reds"), 
                                       domain_min = 0.1,
                                       domain_max = 0.2)
    tlocator = toyplot.locator.Explicit(range(0, ncols), tlabels)    
    llocator = toyplot.locator.Explicit(range(0, nrows), llabels)

    canvas.matrix((lxp, colormap), 
                  rshow=True, 
                  bshow=True,
                  tlocator = tlocator,
                  llocator = llocator)

    # Render heatmap
    plot_filename = '{0}.png'.format(matrix_filename.split('.')[0])
    #toyplot.png.render(canvas, plot_filename)
    toyplot.browser.show(canvas)
    
if __name__ == "__main__":
   main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))