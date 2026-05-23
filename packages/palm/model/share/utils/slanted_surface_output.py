#!/usr/bin/env python

import vtk
import numpy as np
from netCDF4 import Dataset
import json
import os
import sys
from vtk.util import numpy_support

""" 
DESCRIPTION
python3 ./slanted_surface_output.py SURFACE_OUTPUT_file

The script will process all available variables in SURFACE_OUTPUT_file and all avalable timesteps. 
Working both for AV and inst files.
The outputs are .vtu file located in paraview_visualization folder, that in created in SURFACE_OUTPUT_file directory.

Paraview:
Start Paraview
Open .vtu file (or whole time series), open as "XLM Polydata Reader"
Hit Apply
Visualization should appear

NOTE: Ignore all warning in python and in Paraview. If not working let me {buresm@cs.cas.cz} know.
"""

ncfile = sys.argv[1]
output_dir = os.path.join(os.path.dirname(os.path.abspath(ncfile)), 'paraview_visualization')
if not os.path.isdir(output_dir):
    print('Creating output folder {}'.format(output_dir))
    os.makedirs(output_dir)

nc = Dataset(ncfile, 'r')
print('Loading data')
vertices = nc.variables['cct_vertices'][:]
faces = nc.variables['cct_faces'][:]
numPoints = faces.shape[0]

cells = vtk.vtkCellArray()
qp = vtk.vtkPoints()
data_mask = np.ones(numPoints).astype(np.bool)

print('Creating vertices')
for idx, vec in enumerate(vertices):
    qp.InsertNextPoint(vec[::-1])

print('Start preparing polygons')
for idx in range(numPoints):
    poly = faces[idx, :]

    if np.sum(poly != -1) < 3:
        print(np.sum(poly != -1), poly)
        data_mask[idx] = False
        continue

    pol = vtk.vtkPolygon()
    pol.GetPointIds().SetNumberOfIds(np.sum(poly != -1) + 1)
    pol.GetPointIds().SetId(0, poly[0]-1)
    pol.GetPointIds().SetId(1, poly[1]-1)
    pol.GetPointIds().SetId(2, poly[2]-1)
    if poly[3] == -1:
        pol.GetPointIds().SetId(3, poly[0]-1)
        cells.InsertNextCell(pol)
        continue
    else:
        pol.GetPointIds().SetId(3, poly[3] - 1)

    if poly[4] == -1:
        pol.GetPointIds().SetId(4, poly[0]-1)
        cells.InsertNextCell(pol)
        continue
    else:
        pol.GetPointIds().SetId(4, poly[4] - 1)

    if poly[5] == -1:
        pol.GetPointIds().SetId(5, poly[0]-1)
        cells.InsertNextCell(pol)
        continue
    else:
        pol.GetPointIds().SetId(5, poly[5] - 1)

    if poly[6] == -1:
        pol.GetPointIds().SetId(6, poly[0]-1)
        cells.InsertNextCell(pol)
        continue
    else:
        pol.GetPointIds().SetId(6, poly[6] - 1)

print('Done with Polygon preparation')
variables = nc.VAR_LIST.split(';')[1:-1]

p_times = nc.variables['time'][:]
if p_times.mask.size == 1:
    p_times_mask = np.ones(p_times.size, dtype=bool)
else:
    p_times_mask = ~p_times.mask
times = p_times[p_times_mask].data.squeeze()
if times.size == 1:
    times = [times, ]

# TODO: for fast testing
# times = times[0:2]
# variables = variables[0:2]

for var in variables:
    print('Variable: {}'.format(var))
    for it, time in enumerate(times):
        print('Time: {}, {}/{}'.format(time, it, len(times)))
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(qp)
        polydata.SetPolys(cells)
        data = nc.variables[var][it, :]
        data = data[data_mask]

        VTK_data = numpy_support.numpy_to_vtk(num_array=data.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        polydata.GetCellData().SetScalars(VTK_data)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(os.path.join(output_dir, 'slanted_visualization_' + var + '_{:03d}'.format(it) + '.vtu'))
        writer.SetInputData(polydata)
        writer.Write()

print('DONE, all visualization are located here: {}'.format(output_dir))
nc.close()