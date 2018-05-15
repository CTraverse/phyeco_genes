# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:48:36 2017

@author: chuck
"""

import pandas as pd

# This script sorts a master file used in an existing data pipeline, the output here i used in phyeco_gene_grabber.py

phyeco = '/home/chuck/Documents/midas/phyeco.map'

markers = pd.read_csv(phyeco, header=0, sep='\t')

i = 0

markers['marker_value'] = markers['marker_id'].str[1:7]

markers.sorted = markers.sort(['species_id', 'marker_value'], ascending=[1,1])

del markers.sorted['marker_value']
print(markers.sorted)

markers.sorted.to_csv('/home/chuck/Documents/midas/phyeco.map.sorted', sep='\t', index=False)


