#!/Users/yamamoriyuu/anaconda3/bin/python

import argparse as arg
import numpy as np
import mdtraj as md

parser = arg.ArgumentParser(description='calc. 3J const for trajectory.')
parser.add_argument('xtc', metavar='xtc', type=str, help='name of xtc')
parser.add_argument('pdb', metavar='pdb', type=str, help='name of pdb')
parser.add_argument('basename', metavar='basename', type=str, help='basename of output')
parser.add_argument('pairs', metavar='pairs', type=str, help='filename of pairs to calc 3J const')
parser.add_argument('weight', metavar='weight', type=str, help='filename of weights')

args     = parser.parse_args()
xtc      = args.xtc
pdb      = args.pdb
basename = args.basename

t = md.load(xtc, top=pdb)

list, J3_HN_HA = md.compute_J3_HN_HA(t)

pairs = []
f = open(args.pairs, 'r')
l = f.readlines()
f.close()
for item in l:
    pairs.append(item)

w = []
f = open(args.weight, 'r')
l = f.readlines()
f.close()
for item in l:
    i=item.split()
    w.append(float(i[1]))

J3_HN_HA_ave=np.average(J3_HN_HA, axis=0, weights=w)

file_ave = "%s_wave.txt" % basename
f = open(file_ave,'w')
for i, v in enumerate(J3_HN_HA_ave):
    f.write("%5.1f %s_%s_%s_%s\n" % (v,t.topology.atom(list[i][0]), t.topology.atom(list[i][1]), t.topology.atom(list[i][2]), t.topology.atom(list[i][3])))
f.close()


