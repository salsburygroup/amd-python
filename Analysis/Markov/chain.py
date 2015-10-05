#! /usr/bin/env python

def from_vmd(cluster_file):
    """
    ============
    Split frames
    ============

    Usage: %from_vmd(cluster.dat)
    Creates a Markov Chain from vmd clustering output

    inputs:
        -cluster: path to raw VMD output.
        -trjectory: path to trajectory as string
        -output_prefix: prefix for output frames as string. If a prefix is not specified, 'Frame' will be used.

    :Author: Ryan Melvin
    :Date: Fri Oct  2 17:39:22 EDT 2015

    Sources: 

    Dependencies:
        - re
        - pandas

    """
    import re
    import pandas as pd

    with open (cluster_file) as file:
        all_clusters = [result for result in re.findall('\{(.*?)\}',file.read(),re.S)] #re.S makes dot anything, even \n
    
    for line in all_clusters:
        cluster = map(int, all_clusters[line].split(' '))
        dftemp = pd.dataframe(columns = [frame, cluster_number])

    

