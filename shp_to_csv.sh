#!/bin/bash

SHP=$1    # path to .shp file
OUTDIR=$2 # directory to place output (will be created if not exist)

ogr2ogr -f csv -lco GEOMETRY=AS_WKT -progress $OUTDIR $SHP
