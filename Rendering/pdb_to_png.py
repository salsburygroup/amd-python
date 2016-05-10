#!/usr/bin/env python
#Modified from https://www.biostars.org/p/42474/
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pymol

pymol.finish_launching()
pymol.cmd.set('ray_trace_frames', 1)

##
# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]

# Load Structures

pymol.cmd.load(spath, sname)
pymol.cmd.disable("all")
pymol.cmd.enable(sname)
pymol.cmd.show_as('cartoon')
pymol.cmd.bg_color('white')
pymol.cmd.util.chainbow("all")

#Prepare name for image file
outname = os.path.dirname(spath) + '/' + os.path.splitext(sname)[0] + '.png'
pymol.cmd.png(outname)

# Get out!
pymol.cmd.quit()
