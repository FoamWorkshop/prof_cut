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
from cncfclib import *


def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]


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


def dxf_read(files, layer_name, dec_acc, n_arc, l_arc):
    tol = dec_acc
    knots_list = []
    elements_list = []
    hrd_knots_list=[]
    hrd_element_list=[]
    start_point_list = []
    line_count = 0
    arc_count = 0
    circle_count = 0
#    list_entities(dxf)
# ,
#                                                                       [x[1] for x in knots_rank].count(3), ' '
    for shape in dxf.entities:
        if shape.layer == layer_name:
            if shape.dxftype == 'LINE':
                line_count += 1
                p1 = tuple(round(x, tol) for x in shape.start)
                p2 = tuple(round(x, tol) for x in shape.end)
                if p1!=p2:
                    knots_list.append(p1)
                    knots_list.append(p2)
                    elements_list.append([p1, p2])

            if shape.dxftype == 'CIRCLE':
                circle_count += 1
                O = tuple(round(x, tol) for x in shape.center)
                # print('CIRCLE {0}: center: {1}'.format(circle_count, O))
                start_point_list.append(0)

            if shape.dxftype == 'ARC':
                arc_count += 1
                ARC_knots_list = []
                n = n_arc  # number of segments
                min_len = l_arc
                O = shape.center
                R = shape.radius
                angl_1 = shape.startangle * np.pi / 180
                angl_2 = shape.endangle * np.pi / 180

                if angl_2 >= angl_1:
                    angl_list = np.linspace(angl_1, angl_2, n)
                else:
                    angl_list = np.linspace(angl_1, angl_2 + 2 * np.pi, n)

                arc_len = R * np.absolute(angl_2 - angl_1)

                if arc_len / n < min_len:
                    n = max(int(arc_len / min_len), 3)

                for angl in angl_list:
                    ARC_knots_list.append(
                        (round(O[0] + R * np.cos(angl), tol), round(O[1] + R * np.sin(angl), tol), O[2]))

                for i in range(n - 1):
                    elements_list.append(ARC_knots_list[i:i + 2])

                knots_list.extend(ARC_knots_list)

    #hrd_knots_list.append(knots_list[0])
    hrd_element_list.append(elements_list[0])

    for var in elements_list:
        tmp=[var for hrd_var in hrd_element_list if (var[0] in hrd_var) and (var[1] in hrd_var)]
        if  not len(tmp):
            hrd_element_list.append(var)

    # print 'removed elemts: ', len(elements_list)-len(hrd_element_list)
    return (knots_list, hrd_element_list, start_point_list, [line_count, arc_count])


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
parser.add_argument('-a', '--accuracy', type=int,
                    default=dflt_dec_acc, help='decimal accuracy, default: 3')
parser.add_argument('-narc', '--arc_seg_num', type=int,
                    default=dflt_n_arc, help='arc segments number, default: 10')
parser.add_argument('-larc', '--arc_seg_len', type=int,
                    default=dflt_l_arc, help='minimal arc segment length, default: 1')
parser.add_argument('-cw', '--collection_dir', type=int,
                    default=dflt_path_dir, help='closed path collection dir')

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
            knots_list, elements_list, start_point_list, shape_count = dxf_read(
                dxf, layer_name, dec_acc, n_arc, l_arc)

            sorted_knots = knots_dict(knots_list)
            el_kt_list = elements_coords2knots(elements_list, sorted_knots)

            knots_rank = knots_rank_list(el_kt_list, sorted_knots, None)
            master_knot = knots_rank_find(knots_rank, 3)
            IO_knot = knots_rank_find(knots_rank, 1)

            print('{:12}|{:12}|{:12}|{:12}'.format(layer_name, max([x[1] for x in knots_rank]), [x[1] for x in knots_rank].count(1), len(start_point_list))),

            if len(IO_knot) != 1 or len(master_knot) != 1 or IO_knot[0] == None or master_knot[0] == None:
                 print('{0:^20}|'.format('SKIPPED'))
            #
            #     for var in IO_knot:
            #         print("IO knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))
            #
            #     for var in master_knot:
            #         print("master knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))
            #
            #     var =10
            #     print("master knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))

            else:


                io_path = find_path(2, el_kt_list, sorted_knots, None)  # IO path

                # print('{0:3}: {1:4d}|'.format('i/o', len(io_path))),

                last_el, excl_knot = find_l_el(
                    path_dir, el_kt_list, sorted_knots, master_knot[0])

                ct_path = find_path(
                    1, el_kt_list, sorted_knots, excl_knot[0])  # loop path
                # print('{0:3}: {1:4d}|'.format('ct', len(ct_path)))

#{{{
                i_file_name = '{1}{2}.{3}'.format(
                    case_name[0], layer_name, '1', 'knt')
                knots2file(i_file_name, io_path, sorted_knots)
#                print('i_path saved to: {0}'.format(i_file_name))

                o_file_name = '{1}{2}.{3}'.format(
                    case_name[0], layer_name, '3', 'knt')
                knots2file(o_file_name, [var[::-1]
                                         for var in io_path[::-1]], sorted_knots)
#                print('o_path saved to: {0}'.format(o_file_name))

                ct_file_name = '{1}{2}.{3}'.format(
                    case_name[0], layer_name, '2', 'knt')

                knots2file(ct_file_name, ct_path, sorted_knots)
                print("\n")
print "\n end of program. thank you!"
