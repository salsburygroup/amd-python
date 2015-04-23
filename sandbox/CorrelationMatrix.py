#! /usr/bin/env python

import MDAnalysis, MDAnalysis.core.Timeseries

u = MDAnalysis.Universe('/Users/melvrl13/Documents/RyanM/F10/f10.psf', '/Users/melvrl13/Desktop/foldingPaperWorking/weighted3200.dcd')

ag = u.selectAtoms('all')

