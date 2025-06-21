#!/bin/bash

previous_run=/home/aasim.jan/Projects/GW200105_like_dark_sirens/runs/eccentricity_0.1/incl_pi6/rundir
latest_grid=`ls ${previous_run}/overlap-grid-* | sort | tail -n 1`
echo "Latest grid is ${latest_grid}, using that to generate the new grid"

echo "Puffing once"
util_ParameterPuffball.py --parameter mc --parameter eta --fmin 19.0 --fref 19.0 --parameter eccentricity --parameter meanPerAno --parameter chieff_aligned --downselect-parameter eta --downselect-parameter-range  '[0.07,0.21]' --puff-factor 2.5 --force-away 0.05 --downselect-parameter eccentricity --downselect-parameter-range  '[0.0,0.28]'  --downselect-parameter chi1 --downselect-parameter-range  '[0,0.99]' --downselect-parameter chi2 --downselect-parameter-range  '[0,0.05]'  --downselect-parameter mc --downselect-parameter-range  '[3.6,3.64]' --inj-file ${latest_grid} --inj-file-out puffball-1

echo "Puffing twice"
util_ParameterPuffball.py --parameter mc --parameter eta --fmin 19.0 --fref 19.0 --parameter eccentricity --parameter meanPerAno --parameter chieff_aligned --downselect-parameter eta --downselect-parameter-range  '[0.07,0.21]' --puff-factor 3 --force-away 0.05 --downselect-parameter eccentricity --downselect-parameter-range  '[0.0,0.28]'  --downselect-parameter chi1 --downselect-parameter-range  '[0,0.99]' --downselect-parameter chi2 --downselect-parameter-range  '[0,0.05]'  --downselect-parameter mc --downselect-parameter-range  '[3.6,3.64]' --inj-file ${latest_grid} --inj-file-out puffball-2

echo "Combining grids"
ligolw_add ${latest_grid} puffball-1.xml.gz puffball-2.xml.gz -o overlap-grid.xml.gz

echo "Removing intermediate puffball grids"
rm puffball-1.xml.gz puffball-2.xml.gz 
