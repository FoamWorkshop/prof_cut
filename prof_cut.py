#!/usr/bin/python
from __future__ import division

import sys
sys.path.append('/home/adam/Documents/00.projects/02.python/cncfc')
import cncfclib

import os
import argparse
import sys
import dxfgrabber
import numpy as np
from numpy import NaN, pi
from cncfclib import *
import collections
from scipy import interpolate


pt=collections.namedtuple('pt',['x', 'y', 'R', 'th'])
def angl(x):
    out = x

    if x<0:
        out = 2*pi + x
    return out
vangl=np.vectorize(angl)

def radius_segment(P1,P2,P3):
    P4 = np.abs( (P2[:,1]-P1[:,1]) * P3[:,0] - (P2[:,0]-P1[:,0]) * P3[:,1] + P2[:,0]*P1[:,1] - P2[:,1]*P1[:,0]) / np.sqrt( (P2[:,1] - P1[:,1])**2 + (P2[:,0]-P1[:,0])**2 )
    return np.vstack(P4)

def angle_segment(P1,P2,P3):
    k = ((P2[:,1]-P1[:,1]) * (P3[:,0]-P1[:,0]) - (P2[:,0]-P1[:,0]) * (P3[:,1]-P1[:,1])) / ((P2[:,1]-P1[:,1])**2 + (P2[:,0]-P1[:,0])**2)
    P4=np.vstack([P3[:,0] - k * (P2[:,1]-P1[:,1]), P3[:,1] + k * (P2[:,0]-P1[:,0])]).T
    angl=np.arctan2(P4[:,1] - P3[:,1], P4[:,0] - P3[:,0])
    return np.vstack(angl)

def coords2file(name,coords_XU, coords_YV):
    pref_1='xyuv_'
    pref_2='r_'

    f = open(name, 'w')
    for XU, YV in zip(coords_XU,coords_YV):
        f.write('{0:.3f} {1:.3f}\n'.format(XU, YV))

    f.close()

def Rcoords2file(name,coords_R):
    pref_1='xyuv_'
    pref_2='r_'

    f = open(name, 'w')
    for R in coords_R:
        f.write('{0:.3f}\n'.format(R))

    f.close()

def list_entities(dxf):
    dxf_summary = [shape.dxftype for shape in dxf.entities]
    print('{0:<10}: {1}'.format('LINES', dxf_summary.count('LINE')))
    print('{0:<10}: {1}'.format('ARCS', dxf_summary.count('ARC')))

def knots_dict(knots_list):
    return [[i, var] for i, var in enumerate(list(set(knots_list)))]

def dxf_read(files, layer_name, tol):

    line_count = 0
    circle_count = 0
    p_set= []
    c_set= []

    for shape in dxf.entities:
        if shape.layer == layer_name:
            if shape.dxftype == 'LINE':
                line_count += 1
                p_st = shape.start
                p_en = shape.end
                p_set.append( pt(round(p_st[0],tol), round(p_st[1],tol),NaN, NaN))
                p_set.append( pt(round(p_en[0],tol), round(p_en[1],tol),NaN, NaN))

            if shape.dxftype == 'CIRCLE':
                circle_count += 1
                p_cnt = shape.center
                c_set.append(pt(round(p_cnt[0],tol), round(p_cnt[1],tol),NaN, NaN))

    # hrd_element_list.append(elements_list[0])

    return (p_set, c_set, [line_count, circle_count])


#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'
dflt_layer_list= 'all'

dflt_dec_acc = 4  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 0.1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
prof_segment_list = []
prof_R_list = []
prof_th_list = []

#*********************************************************************PROGRAM
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
--------------------------
lofted shape rules for dxf
--------------------------

the lofted shape is defined by:

1. profile layers:
[prefix_XX_YYYY] where: XX   - part numeber; YYYY - profile number.
Each profile is build with n LINES and a CIRCLE. The other shapes are neglected. he profile is builf with n LINES. The lines have one end coincident (O_POINT) and the other is coincident with the cut contour of the profile. There is one distiguinshed line (0_LINE) which end is marked with a CIRCLE center. This line is treated as the first and the other lines are ordered counterclockwise. All profiles must be build with THE SAME NUMBER of LINES.

2. spin layer:
[prefix_XX_spin] where: XX   - part numeber
The spin is built with n CIRCLES which define JOINTS. The other shapes are neglected. Y position of the CIRCLE center points define coresponding positions of the profile O_POINTS. O_POINTs are assigned to JOINTS as the ordered profile names (LAYER NAMES) along the Y axis. Number of joints must correspond to the number of the profiles.

    ''',
                                usage='%(prog)s [-i dxf_filename] [-l lofted_shape_prefix]')



parser.add_argument('-i', '--input', type=str, required=True, help='input filename')
parser.add_argument('-l', '--profile_prefix', type=str, required=True, help='profile data prefix. prefix_xxxx - profile; prefix_spin - loft spin data')
parser.add_argument('-a', '--accuracy', type=int, default=dflt_dec_acc, help='decimal accuracy, default: 3')
parser.add_argument('-n_sect', '--nb_of_sections', type=int, default=dflt_n_arc, help='number of interpolation sections along the spin, default: number of profiles')
parser.add_argument('-interp', '--interp_meth', type=str, default=dflt_l_arc, help='interpolation method between profiles:\n, linear -default\n ')
# parser.add_argument('-cw', '--collection_dir', type=int,default=dflt_path_dir, help='closed path collection dir')

args = parser.parse_args()

dxf_list = args.input
layer_list = args.profile_prefix
dec_acc = args.accuracy
aaaa=args.interp_meth
ffff=args.nb_of_sections

# n_arc = args.arc_seg_num
# l_arc = args.arc_seg_len
# path_dir = args.collection_dir

dir_path = os.getcwd()
dxf_files = [i for i in os.listdir(dir_path) if i.endswith('.dxf')]

if not dxf_list:
    print('the input file is not specified, use: -i "file name"')
elif not (dxf_list in dxf_files):
    print('the current dir does not include requested dxf files')
elif not layer_list:
    print('the profile prefix is not specified, use: -l "profile prefix"')

else:
    prof_pref= layer_list
    prof_spin= layer_list+'_spin'

    files_dxf_member = dxf_list

    print("\n{0:24} {1:8}".format('extracting from: ',files_dxf_member))
    print('{0}'.format('-' * 13*5))

    case_name = os.path.splitext(files_dxf_member)
    dxf = dxfgrabber.readfile(files_dxf_member, {"assure_3d_coords": True})
    dxf_layers = dxf.layers

    prof_layer_name_list = [var.name for var in dxf_layers if prof_pref in var.name]
    prof_spin_layer_name_list = [var.name for var in dxf_layers if prof_spin in var.name]

    for layer_name in sorted(prof_layer_name_list):
        p_set, c_set, shape_count = dxf_read(dxf, layer_name, dec_acc)

        p_cnt = [x for x, y in collections.Counter(p_set).items() if y >  1]
        p_sec = [x for x, y in collections.Counter(p_set).items() if y == 1]

        if len(c_set) == 0:
            print('no circle defining the start point')
        elif len(c_set) >1:
            print('more than one circle defining the start point:\n {}'.format(c_set))
        elif not(c_set[0] in p_sec):
            print('circle center position does not match any section node')
        elif len(p_cnt) != 1:
                print('error, not single center point: {}'.format(p_cnt))
        elif len(p_set)/2 != len(p_sec):
            print('error some extra sections are present')
        else:

            print('{}'.format(layer_name))

            O = p_cnt[0]

            p_cnt_arr=np.array(p_cnt[0])
            p_sec_arr=np.array(p_sec)
            p_0_arr=np.ones((np.size(p_sec_arr[:,0]) , 2)) * np.array([O.x, O.y])

            p_0 = collections.namedtuple('p_0',['x', 'y', 'R', 'th'])
            p_0.x = c_set[0].x
            p_0.y = c_set[0].y

            p_0.R = np.sqrt((p_0.x - O.x)**2 + (p_0.y - O.y)**2)
            p_0.th = np.arctan2(p_0.y - O.y, p_0.x - O.x)
            p_sec_arr[:,2]=np.sqrt((p_sec_arr[:,0] - O.x)**2 + (p_sec_arr[:,1] - O.y)**2)
            p_sec_arr[:,3]=vangl(np.arctan2(p_sec_arr[:,1] - O.y, p_sec_arr[:,0] - O.x) - p_0.th)
            p_sec_arr=np.vstack(sorted(p_sec_arr, key=lambda a_entry: a_entry[3]))
            segment = np.hstack([p_sec_arr[:,:2],np.roll(p_sec_arr[:,:2],-1, axis=0)])
            # prof_R_list.append( np.hstack([radius_segment(segment[:,:2], segment[:,2:4], p_0_arr[:,:2])]) )
            prof_R_list.append(  radius_segment(segment[:,:2], segment[:,2:4], p_0_arr[:,:2]) )
            prof_th_list.append( angle_segment( segment[:,:2], segment[:,2:4], p_0_arr[:,:2]) *180/pi )

    for layer_name in sorted(prof_spin_layer_name_list):
        dummy_1, c_set, dummy_2 = dxf_read(dxf, layer_name, dec_acc)
        z_set= np.hstack([var.y for var in c_set])
        z_set= np.sort(z_set)
    print('R list \n {}'.format(np.hstack(prof_R_list).T))
    print('th list \n {}'.format(np.hstack(prof_th_list).T))
    print('spin sections \n {}'.format(z_set))

    n_sect = 10
    yv = np.hstack((np.linspace(z_set[0],z_set[-1],n_sect),z_set))
    yv = np.unique(yv)
    yv = np.sort(yv)

    R_val =np.hstack(prof_R_list).T
    th_val=np.hstack(prof_th_list).T

    print('R val')
    for i in range(len(R_val[0,:])):
        path_R = interpolate.interp1d(z_set, R_val[:,i])
        print(path_R(yv))
        print(yv)
        name='{0}_xyuv_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
        coords2file(name, path_R(yv), yv)

    for i in range(len(th_val[0,:])):
        path_th = interpolate.interp1d(z_set, th_val[:,i])
        name='{0}_r_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
        print(path_th(yv))
        Rcoords2file(name, path_th(yv))

print "\n end of program. thank you!"
