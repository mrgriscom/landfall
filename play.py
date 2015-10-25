import shapely
import shapely.wkt
import shapely.ops
import shapely.validation
from shapely.geometry import Polygon, MultiPolygon
import sys
import csv
import itertools
import collections
import tempfile
import json

COAST_DATA = '/home/drew/tmp/landbig/land_polygons.csv'

def load(path):
    csv.field_size_limit(sys.maxsize)
    with open(path) as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            poly = shapely.wkt.loads(row['WKT'])
            yield poly

for i, k in enumerate(load(COAST_DATA)):
    import pdb;pdb.set_trace()
    
    

def quadrant_overlap(poly, quadrant):

    for each ring in poly:
        for each segment (p0, p1) in ring:
            if crosses quadrant edge:
                partial overlap
                determine the crossing indexes for quicker search at next level

    no overlap: either fully inside or outside
    test centerpoint of quadrant is interior or exterior



test_centerpoint:
  determine all crossings of coordindate along a certain axis
  if number of crossings with other coord greater than
    odd: inside
    even: outside
