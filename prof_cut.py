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


pt=collections.namedtuple('pt',['x', 'y', 'R', 'th'])
def angl(x):
    out = x

    if x<0:
        out = 2*pi + x
    return out
vangl=np.vectorize(angl)

def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]

def radius_segment_intersection():
    # // first convert line to normalized unit vector
    # double dx = x2 - x1;
    # double dy = y2 - y1;
    # double mag = sqrt(dx*dx + dy*dy);
    # dx /= mag;
    # dy /= mag;
    #
    # // translate the point and get the dot product
    # double lambda = (dx * (x3 - x1)) + (dy * (y3 - y1));
    # x4 = (dx * lambda) + x1;2
    # y4 = (dy * lambda) + y1;
    return 0

def radius_segment(P1,P2,P3):
    return np.abs(-np.linalg.det([P2-P1, P3])-np.linalg.det([P2,P1]))/np.sqrt(np.sum(np.square(P2-P1)))



def sub_points(p1, p2):
    vect = []
    p1 = [x for x in p1[0]]
    p2 = [x for x in p2[0]]

    if len(p1) == len(p2):
        for i, n in enumerate(p2):
            vect.append(n - p1[i])
        return vect
    return len(p1) * [None]


def knots_rank_find(knots_rank, rank):
    knots = [x[0] for x in knots_rank if x[1] == rank]
    if len(knots) > 0:
        return knots
    else:
        return [None]


def knots_rank_list(el_kt_list, sorted_knots, skip_knot):
    knots_rank = []
    for var in sorted_knots:
        if var[0] == skip_knot:
            knots_rank.append([var[0], None])
        else:
            knots_rank.append([var[0], [x for a in el_kt_list for x in a].count(var[0])])
    return knots_rank


def knot2coord(sorted_knots, knot):
    for var in sorted_knots:
        if var[0] == knot:
            return var[1]
    return None


def knots2file(name, io_path, sorted_knots):
    f = open(name, 'w')

    for var in io_path:
        coord = knot2coord(sorted_knots, var[0])
        f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    f.close()


def list_entities(dxf):
    dxf_summary = [shape.dxftype for shape in dxf.entities]
    print('{0:<10}: {1}'.format('LINES', dxf_summary.count('LINE')))
    print('{0:<10}: {1}'.format('ARCS', dxf_summary.count('ARC')))


def knots_dict(knots_list):
    return [[i, var] for i, var in enumerate(list(set(knots_list)))]


def elements_coords2knots(el_list, kt_list):
    el_kt_list = []
    for el in el_list:
        for kt in kt_list:
            if kt[1] == el[0]:
                p1 = kt[0]
            if kt[1] == el[1]:
                p2 = kt[0]
        el_kt_list.append([p1, p2])
    return el_kt_list


def elements_knots2coords(el_list, kt_list):
    el_coord_list = []
    for el in el_list:
        for kt in kt_list:
            if kt[0] == el[0]:
                p1 = kt[1]
            if kt[0] == el[1]:
                p2 = kt[1]
        el_coord_list.append([p1, p2])
    return el_coord_list


def knots_rank_list_summary(knots_rank):
    print('{0:<16}: {1}'.format('IO knots', [x[1] for x in knots_rank].count(1),
                                [x[0] for x in knots_rank if x[1] == 1]))
    print('{0:<16}: {1}'.format('master knots', [x[1] for x in knots_rank].count(3),
                                [x[0] for x in knots_rank if x[1] == 3]))
    print('{0:<16}: {1}'.format('chain knots', [x[1] for x in knots_rank].count(2),
                                [x[0] for x in knots_rank if x[1] == 2]))


def paths_summary(io_path, ct_path):
    print("-----IN PATHS-----")
    for var in io_path:
        print("{0} {1}".format(var[0], var[1]))

    print("-----CT PATHS-----")
    for var in ct_path:
        print("{0} {1}".format(var[0], var[1]))

    print("-----OUT PATHS-----")
    for var in reversed(io_path):
        print("{0} {1}".format(var[1], var[0]))


def find_path(crit, el_kt_list, sorted_knots, excl_knot):
    path = []

    knots_rank = knots_rank_list(el_kt_list, sorted_knots, excl_knot)

    curr_knot = knots_rank_find(knots_rank, 1)
    last_knot = knots_rank_find(knots_rank, 3)

    curr_element=[]
    print("{0:8}").format('OK'),
    while not ((curr_element is None) or curr_knot[0]==last_knot[0]):
        # print '\rpool size: {0}'.format(len(path)),

        curr_element=next((element for element in el_kt_list if curr_knot[0] in element), None)
        if not (curr_element is None):

            if curr_element[0] == curr_knot[0]:
                curr_knot=[curr_element[1]]
                path.append(curr_element)
            else:
                curr_knot=[curr_element[0]]
                path.append(curr_element[::-1])

            el_kt_list.remove(curr_element)

    if crit == 1:
        path.append([path[-1][1], path[0][0]])
    # print '\n'
    return path


def find_l_el(read_dir, el_kt_list, sorted_knots, master_knot):
    # find all elements including master_knot and put into el_index list
    el_index = [i for i, element in enumerate(
        el_kt_list) if master_knot in element]/home/user/python-libs

    seg1, seg2 = elements_knots2coords(
        [el_kt_list[i] for i in el_index], sorted_knots)

    if cw_order(seg1, seg2) == read_dir:
        cur_ind = el_index[1]
    else:
        cur_ind = el_index[0]

    last_el = el_kt_list.pop(cur_ind)
    excl_knot = [x for x in last_el if x != master_knot]  # take the other knot

    return (last_el, excl_knot)


def cw_order(seg1, seg2):
    common_el = [x for x in list(set(seg1) & set(seg2))]
    u = sub_points(common_el, list(set(seg1) - set(common_el)))
    v = sub_points(common_el, list(set(seg2) - set(common_el)))
    if cross_prod(u, v) > 0:
        return False
    else:
        return True


def print_list(data_list, common_text):
    for var in data_list:
        print('{0} {1}'.format(common_text, var))


def dxf_read(files, layer_name, tol, n_arc, l_arc):

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
prof_pref= 'prof'
dflt_dec_acc = 4  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 0.1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction

#*********************************************************************PROGRAM
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i', '--input', nargs='+', help='input filenames')
parser.add_argument('-l', '--layer', nargs='+', help='input layers')
parser.add_argument('-a', '--accuracy', type=int, default=dflt_dec_acc, help='decimal accuracy, default: 3')
parser.add_argument('-narc', '--arc_seg_num', type=int, default=dflt_n_arc, help='arc segments number, default: 10')
parser.add_argument('-larc', '--arc_seg_len', type=int, default=dflt_l_arc, help='minimal arc segment length, default: 1')
parser.add_argument('-cw', '--collection_dir', type=int,default=dflt_path_dir, help='closed path collection dir')

args = parser.parse_args()

dxf_list = args.input
layer_list = args.layer
dec_acc = args.accuracy
n_arc = args.arc_seg_num
l_arc = args.arc_seg_len
path_dir = args.collection_dir

dir_path = os.getcwd()
dxf_files = [i for i in os.listdir(dir_path) if i.endswith('.dxf')]

if dxf_list:
    print dxf_files
    files_dxf = [i for i in dxf_files if i in dxf_list]
else:
    files_dxf = dxf_files

if not files_dxf:
    print 'dir does not include requested dxf files'

# else, execute the program
else:

    print('SETTINGS:')
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'decimal accuracy', dec_acc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'arc segments count', n_arc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'minimal arc segment length', l_arc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'closed path collection dir', path_dir))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'files', files_dxf))


    for i, files_dxf_member in enumerate(files_dxf):
        print("\n{0:24} {1:8}".format('extracting from: ',files_dxf_member))
        print('\n{:12}|{:12}|{:12}|{:12}|{:12}|'.format('layer', 'max rank', 'rank 1','start mark.','offs. ang.'))
        print('{0}'.format('-' * 13*5))

        case_name = os.path.splitext(files_dxf_member)
        dxf = dxfgrabber.readfile(files_dxf_member, {"assure_3d_coords": True})
        dxf_layers = dxf.layers

        if layer_list:
            layer_name_list= [var.name for var in dxf_layers if var in layer_list]
        else:
            layer_name_list = [var.name for var in dxf_layers if prof_pref in var.name]

        for layer_name in sorted(layer_name_list):
            p_set, c_set, shape_count = dxf_read(dxf, layer_name, dec_acc, n_arc, l_arc)
            print('{}'.format('point data'))
            print(np.array(p_set))
            print('{}'.format('start data'))
            print(np.array(c_set))
            print(np.array(list(collections.OrderedDict.fromkeys(p_set))))
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
                print('start point: \n{}'.format(c_set[0]))
                print('center point: \n{}'.format(p_cnt[0]))
                print('section points: \n{}'.format(p_sec))
                O = p_cnt[0]

                p_cnt_arr=np.array(p_cnt[0])
                p_sec_arr=np.array(p_sec)

                # p_ang0_arr = np.array(c_set[0])
                # p_ang0_arr[:,2]=np.sqrt((p_ang0_arr[:,0] - O.x)**2 + (p_ang0_arr[:,1] - O.y)**2)
                # p_ang0_arr[:,3]=np.arccos((p_ang0_arr[:,0] - O.x) / (p_ang0_arr[:,2]))
                p_0 = collections.namedtuple('p_0',['x', 'y', 'R', 'th'])
                p_0.x = c_set[0].x
                p_0.y = c_set[0].y
                p_0.R = np.sqrt((p_0.x - O.x)**2 + (p_0.y - O.y)**2)
                p_0.th = np.arctan2(p_0.y - O.y, p_0.x - O.x)
                print('pppppp')
                print(p_0.th)
                p_sec_arr[:,2]=np.sqrt((p_sec_arr[:,0] - O.x)**2 + (p_sec_arr[:,1] - O.y)**2)
                p_sec_arr[:,3]=vangl(np.arctan2(p_sec_arr[:,1] - O.y, p_sec_arr[:,0] - O.x) - p_0.th)
                p_sec_arr=np.vstack(sorted(p_sec_arr, key=lambda a_entry: a_entry[3]))
                print(p_sec_arr)
            # print('{:12}|{:12}|{:12}|{:12}'.format(layer_name, max([x[1] for x in knots_rank]), [x[1] for x in knots_rank].count(1), len(start_point_list))),


print "\n end of program. thank you!"
print(radius_segment(np.array([1,1]),np.array([2,2]),np.array([0,1])))
