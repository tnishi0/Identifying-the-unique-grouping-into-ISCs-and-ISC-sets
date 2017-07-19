import networkIRR_subgroups_independent
import random
import itertools
import sys

from time import time
reload(networkIRR_subgroups_independent)
from networkIRR_subgroups_independent import *

if len(sys.argv) != 2:
	print("./sage all_CS_patterns.sage name_of_adjacency_matrix_file")
	sys.exit()

#filename = raw_input("Enter name of adjacency matrix file: ")
filename = str(sys.argv[1])
filename2 = filename.split()
filename = filename2[0]

ins = open(filename, "r")

_Aij32 = []
for line in ins:
	number_strings = line.split()
	numbers = [int(n) for n in number_strings]
	if(len(numbers)>0):
		_Aij32.append(numbers)

Aij32 = np.array(_Aij32)

data32=NetworkIRR(Aij32)
nnodes = data32._nnodes

data32.get_orbits()
AutG = data32.get_automorphism_group()

gen_orbits = [] 
for i in range(len(AutG.gens())):
	gen_orbits.append(PermutationGroup([AutG.gens()[i]]).orbits())

supp = []
for i in range(len(gen_orbits)):
	_supp = []
	for j in range(len(gen_orbits[i])):
		if len(gen_orbits[i][j]) > 1:
			for k in range(len(gen_orbits[i][j])):
				_supp.append(gen_orbits[i][j][k])
		
	if len(_supp) > 1:
		supp.append(_supp)

Hcheck = [0]*len(supp)

H = []
for i in range(len(supp)):
	if Hcheck[i]==0:
		Hcheck[i] = 1
		_H = [i]
		minnewelement = 0
		maxnewelement = 1
		while True:
			for j in range(minnewelement, maxnewelement): # for added generators
				for k in range(len(supp)): # searching common elements again
					if Hcheck[k] == 0:
						support_intersect_check = 0
						for l in range(len(supp[_H[j]])):
							for m in range(len(supp[k])):
								if supp[_H[j]][l] == supp[k][m]:
									support_intersect_check = 1
									break

							if support_intersect_check == 1:
								break
				
						if support_intersect_check == 1:
							Hcheck[k] = 1
							_H.append(k)

			minnewelement = maxnewelement
			maxnewelement = len(_H)
				
			if maxnewelement == minnewelement:
				break

		H.append(_H)

Horbits = []
for i in range(len(H)):
	gen_subH = []
	for j in range(len(H[i])):
		gen_subH.append(AutG.gens()[H[i][j]])

	_Horbits = []
	for j in range(len(PermutationGroup(gen_subH).orbits())):
		if len(PermutationGroup(gen_subH).orbits()[j])>1:
			_Horbits.append(PermutationGroup(gen_subH).orbits()[j])
	Horbits.append(_Horbits)

for i in range(len(Horbits)):
	for j in range(len(Horbits[i])):
		Horbits[i][j] = sorted(Horbits[i][j], key=int)

Hsupports = []
for i in range(len(Horbits)):
	_supp = []
	for j in range(len(Horbits[i])):
		for k in range(len(Horbits[i][j])):
			_supp.append(Horbits[i][j][k])
	
	Hsupports.append(_supp)

for i in range(len(Hsupports)):
	Hsupports[i] = sorted(Hsupports[i], key=int)

Hsupports_first_nodes = []

for i in range(len(Hsupports)):
	Hsupports_first_nodes.append(Hsupports[i][0]) 

reorder_H = []
for i in range(len(Hsupports)):
	reorder_H.append(i)

while True:
	num_swap = 0
	for i in range(len(Hsupports_first_nodes)-1):
		if Hsupports_first_nodes[i] > Hsupports_first_nodes[i+1]:
			swap = Hsupports_first_nodes[i]
			Hsupports_first_nodes[i] = Hsupports_first_nodes[i+1]
			Hsupports_first_nodes[i+1] = swap
			Hswap = reorder_H[i]
			reorder_H[i] = reorder_H[i+1]
			reorder_H[i+1] = Hswap
			num_swap += 1

	if num_swap == 0:
		break

print "\n\n****** Geometric decomposition ******\n"

if len(Hsupports)>0:
	print "G = H1",
	for i in range(1, len(Hsupports)):
		print("x H%d"%(i+1)),

print "\n"
print "support of subgroup components:\n"

for i in range(len(Hsupports)):
	print ("H%d:"%(i+1)),
	for j in range(len(Hsupports[reorder_H[i]])):
		print Hsupports[reorder_H[i]][j]+1, 
	print ""

Hgenerators = []

for i in range(len(H)):
	_Hgenerators = []
	for j in range(len(H[i])):
		_Hgenerators.append(AutG.gens()[H[i][j]])
	
	Hgenerators.append(_Hgenerators)

H_belong_node = [-1]*nnodes

for j in range(len(H)):
	for k in range(len(Hsupports[j])):
		H_belong_node[Hsupports[j][k]] = j
	
print "\n\n"
subgroup_NTorbitsets = []
for i in range(len(H)):
	H_nu_subgroup_NTorbitsets = []
	print ""
	print("Possible CS patterns for H%d"%(i+1))
	H_nu_subgroups = PermutationGroup(Hgenerators[reorder_H[i]]).subgroups()
	for j in range(len(H_nu_subgroups)):
		subgroup_orbits = H_nu_subgroups[j].orbits()
		subgroup_NTorbits = []
		for k in range(len(subgroup_orbits)):
			if len(subgroup_orbits[k]) > 1:
				subgroup_NTorbits.append(sorted(subgroup_orbits[k], key=int))
 
		same_NTorbitset_exist_check = 0
		for k in range(len(H_nu_subgroup_NTorbitsets)):
			if len(H_nu_subgroup_NTorbitsets[k]) != len(subgroup_NTorbits):
				continue                                                       
			same_NTorbitset_check = 0       
			for l in range(len(subgroup_NTorbits)):            
				same_NTorbit_exist_check = 0
				for m in range(len(H_nu_subgroup_NTorbitsets[k])):
					if subgroup_NTorbits[l] == H_nu_subgroup_NTorbitsets[k][m]:
						same_NTorbit_exist_check = 1
						break

				if same_NTorbit_exist_check == 0: 
					same_NTorbitset_check = 1
					break

				if same_NTorbit_exist_check == 1:
					continue

			if same_NTorbitset_check == 0:
				same_NTorbitset_exist_check = 1
				break

		if same_NTorbitset_exist_check == 0:
			if len(subgroup_NTorbits) == 0:
				continue
			H_nu_subgroup_NTorbitsets.append(subgroup_NTorbits)

	subgroup_NTorbitsets.append(H_nu_subgroup_NTorbitsets)

	num_option = 1
	for j in range(len(H_nu_subgroup_NTorbitsets)):
		print "(",num_option,") ", 
		for k in range(len(H_nu_subgroup_NTorbitsets[j])):
			for l in range(len(H_nu_subgroup_NTorbitsets[j][k])):
				if l == 0:
					print ("[%d"%(H_nu_subgroup_NTorbitsets[j][k][l]+1)),
				elif l == len(H_nu_subgroup_NTorbitsets[j][k])-1:
					print ("%d]"%(H_nu_subgroup_NTorbitsets[j][k][l]+1)),
				else:
					print H_nu_subgroup_NTorbitsets[j][k][l]+1,
		print ""
		num_option += 1
