#import networkIRR_subgroups_independent
#import random
#import itertools
#from time import time
#reload(networkIRR_subgroups_independent)
#from networkIRR_subgroups_independent import *
import numpy as np
import sys

def gram_schmidt(original_basis):
	original_basis = np.array(original_basis)
	new_basis = []
	for i in range(len(original_basis)):
		_new_basis = []
		for j in range(len(original_basis[i])):
			_new_basis.append(original_basis[i][j])
		for j in range(len(new_basis)):
			_new_basis -= np.dot(new_basis[j], original_basis[i].transpose())*new_basis[j]
		_new_basis /= np.linalg.norm(_new_basis)
		new_basis.append(_new_basis)

	return np.array(new_basis)

if len(sys.argv) is not 2:
	print("python grouping_clusters.py name_of_adjacency_matrix_file")
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

nnodes = len(Aij32)

###### sgorbits include isolated clusters ######

num_NTorbits = input("Enter the number of nontrivial clusters: ")
NTsgorbits = []
for i in range(num_NTorbits):
	print ("cluster #%d (example: 3 5 6): "%(i+1)),
	_NTsgorbits = map(int, raw_input().split())
	for j in range(len(_NTsgorbits)):
		_NTsgorbits[j] = _NTsgorbits[j]-1
	NTsgorbits.append(_NTsgorbits)

for i in range(len(NTsgorbits)):
	NTsgorbits[i] = sorted(NTsgorbits[i], key=int)

orbit_check = [0]*nnodes
for i in range(len(NTsgorbits)):
	for j in range(len(NTsgorbits[i])):
		orbit_check[NTsgorbits[i][j]] = 1

sgorbits = []
for i in range(len(NTsgorbits)):
	sgorbits.append(NTsgorbits[i])

for i in range(nnodes):
	if orbit_check[i] == 0:
		sgorbits.append([i])

print ("\n\nAll clusters(%d): "%(len(sgorbits))),

for i in range(len(sgorbits)):
	if len(sgorbits[i])>1:
		for j in range(len(sgorbits[i])):
			if j == 0:
				print ("[%d"%(sgorbits[i][j]+1)),
			elif j == len(sgorbits[i])-1:
				print ("%d]"%(sgorbits[i][j]+1)),
			else:
				print sgorbits[i][j]+1,
	
	else:
		print ("[%d]"%(sgorbits[i][0]+1)),
	
T = []
for i in range(nnodes):
	_T = [0]*nnodes
	T.append(_T)

parallel_coord = 0 
transverse_coord = len(sgorbits)
for i in range(len(sgorbits)):
	VR = []
	for j in range(len(sgorbits[i])):
		_VR = []
		for k in range(len(sgorbits[i])):
			if k == 0:
				_VR.append(1)
			elif k == j+1:
				_VR.append(1)
			else:
				_VR.append(0)

		VR.append(_VR)

#	Tblock = matrix(RDF, VR)
	Q, R = np.linalg.qr(VR)
	Tblock = Q.transpose()

	for j in range(len(sgorbits[i])):
		T[parallel_coord][sgorbits[i][j]] = Tblock[0][j]
	parallel_coord += 1

	for j in range(1, len(sgorbits[i])):
		for k in range(len(sgorbits[i])):
			T[transverse_coord][sgorbits[i][k]] = Tblock[j][k]	
		transverse_coord += 1

InvT = np.transpose(T)
B = np.dot(T, np.dot(Aij32, InvT))

transverse_coord = []
n_transverse_coord = len(sgorbits)
for i in range(len(sgorbits)):
	_transverse_coord = []
	for j in range(len(sgorbits[i])-1):
		_transverse_coord.append(n_transverse_coord)
		n_transverse_coord += 1
	transverse_coord.append(_transverse_coord)
	
ISC_ISCset_check = [0]*len(transverse_coord)
ISC_ISCset = []

for i in range(len(transverse_coord)):
	if ISC_ISCset_check[i] == 0:
		ISC_ISCset_check[i] = 1
		_ISC_ISCset = [i]
		minnewelement = 0
		maxnewelement = 1
		while True:
			for j in range(minnewelement, maxnewelement):
				for k in range(len(transverse_coord)):
					if ISC_ISCset_check[k] == 0:
						Intertwined_check = 0
						for l in range(len(transverse_coord[_ISC_ISCset[j]])):
							for m in range(len(transverse_coord[k])):
								if np.abs(B[transverse_coord[_ISC_ISCset[j]][l]][transverse_coord[k][m]]) > 1e-12:
									Intertwined_check = 1
									break

							if Intertwined_check == 1:
								break

						if Intertwined_check == 1:
							ISC_ISCset_check[k] = 1
							_ISC_ISCset.append(k)

			minnewelement = maxnewelement
			maxnewelement = len(_ISC_ISCset)

			if maxnewelement == minnewelement:
				break

		ISC_ISCset.append(_ISC_ISCset)	

ISC_ISCset_cluster = []

for i in range(len(ISC_ISCset)):
	_ISC_ISCset = []
	for j in range(len(ISC_ISCset[i])):
		_ISC_ISCset.append(sgorbits[ISC_ISCset[i][j]])

	ISC_ISCset_cluster.append(_ISC_ISCset)

print "\n\nISCs:"

num_ISCs = 0
for i in range(len (ISC_ISCset_cluster)):
	if len(ISC_ISCset_cluster[i]) == 1 and len(ISC_ISCset_cluster[i][0]) > 1:
		num_ISCs += 1
		for l in range(len(ISC_ISCset_cluster[i][0])):
			if l == 0:
				print ("[%d"%(ISC_ISCset_cluster[i][0][l]+1)),
			elif l == len(ISC_ISCset_cluster[i][0])-1:
				print ("%d]"%(ISC_ISCset_cluster[i][0][l]+1)),
			else:
				print ISC_ISCset_cluster[i][0][l]+1,
		print ""
		
if num_ISCs == 0:
	print "none"

print "\n\nISC sets:"

num_ISCset = 0
for i in range(len(ISC_ISCset_cluster)):
	if len(ISC_ISCset_cluster[i]) > 1:
		num_ISCset += 1
		for j in range(len(ISC_ISCset_cluster[i])):
			if j == 0:
				for k in range(len(ISC_ISCset_cluster[i][j])):
					if k == 0:
						print ("{[%d"%(ISC_ISCset_cluster[i][j][k]+1)),
					elif k == len(ISC_ISCset_cluster[i][j])-1:
						 print ("%d]"%(ISC_ISCset_cluster[i][j][k]+1)),
					else:
						print ISC_ISCset_cluster[i][j][k]+1,

			elif j == len(ISC_ISCset_cluster[i])-1:
				for k in range(len(ISC_ISCset_cluster[i][j])):
					if k == 0:
						print ("[%d"%(ISC_ISCset_cluster[i][j][k]+1)),
					elif k == len(ISC_ISCset_cluster[i][j])-1:
						print ("%d]}"%(ISC_ISCset_cluster[i][j][k]+1)),
					else:
						print ISC_ISCset_cluster[i][j][k]+1,

			else:
				for k in range(len(ISC_ISCset_cluster[i][j])):
					if k == 0:
						print ("[%d"%(ISC_ISCset_cluster[i][j][k]+1)),
					elif k == len(ISC_ISCset_cluster[i][j])-1:
						print ("%d]"%(ISC_ISCset_cluster[i][j][k]+1)),
					else:
						print ISC_ISCset_cluster[i][j][k]+1,
		print ""

if num_ISCset == 0:
	print "none"

U = []
for i in range(nnodes):
	_U = []
	for j in range(nnodes):
		if j == i:
			_U.append(1)
		if j is not i:
			_U.append(0)
	U.append(_U)

for i in range(len(ISC_ISCset)):
	if len(ISC_ISCset[i]) == 1 and len(sgorbits[ISC_ISCset[i][0]]) > 1:
		Bblock = []	
		Cm = ISC_ISCset[i][0]
		for j in range(len(transverse_coord[Cm])):
			_Bblock = []
			for k in range(len(transverse_coord[Cm])):
				_Bblock.append(B[transverse_coord[Cm][j]][transverse_coord[Cm][k]])
			Bblock.append(_Bblock)

		VR = np.linalg.eig(Bblock)
		eigenvectors = VR[1]
		eigenvectors = eigenvectors.transpose()

		for j in range(len(eigenvectors)):
			eigenvectors[j] = eigenvectors[j]/np.linalg.norm(eigenvectors[j])

		check = [0]*len(eigenvectors)
		degen_eigenvalues = []
		degen_eigenvectors = []
		n_not_check = len(eigenvectors)
		while n_not_check > 0:
			n_not_check = 0
			for j in range(len(eigenvectors)):
				if check[j] == 0:
					_degen_eigenvectors = []
					_degen_eigenvalues = []
					degen_check_eigenvector = j
					check[degen_check_eigenvector] = 1
					__eigenvector = []
					for k in range(len(eigenvectors[degen_check_eigenvector])):
						__eigenvector.append(eigenvectors[degen_check_eigenvector][k])
					_degen_eigenvectors.append(__eigenvector)
					_degen_eigenvalues.append(VR[0][degen_check_eigenvector])
					n_not_check += 1
					for k in range(degen_check_eigenvector+1, len(eigenvectors)):
						if np.abs(VR[0][degen_check_eigenvector] - VR[0][k]) < 1e-12:
							__eigenvector = []
							for l in range(len(eigenvectors[k])):
								__eigenvector.append(eigenvectors[k][l])
							_degen_eigenvectors.append(__eigenvector)
							_degen_eigenvalues.append(VR[0][k])
							check[k] = 1
							n_not_check += 1

					degen_eigenvectors.append(_degen_eigenvectors)
					degen_eigenvalues.append(_degen_eigenvalues)

		for j in range(len(degen_eigenvectors)):
			degen_eigenvectors[j] = gram_schmidt(degen_eigenvectors[j])

		eigenvectors = []
		for j in range(len(degen_eigenvectors)):
			for k in range(len(degen_eigenvectors[j])):
				eigenvectors.append(degen_eigenvectors[j][k])

		Ublock = []
		for j in range(len(transverse_coord[Cm])):
			_Ublock = []
			for k in range(len(transverse_coord[Cm])):
				_Ublock.append(eigenvectors[j][k])
			Ublock.append(_Ublock)

		for j in range(len(transverse_coord[Cm])):
			for k in range(len(transverse_coord[Cm])):
				U[transverse_coord[Cm][j]][transverse_coord[Cm][k]] = Ublock[j][k]

T = np.dot(U, T)
InvT = np.transpose(T)
B = np.dot(T, np.dot(Aij32, InvT))

for i in range(len(sgorbits)):
	for j in range(nnodes):
		if np.abs(T[i][j]) > 1e-12:
			if T[i][j] < 0:
				T[i][j] = -T[i][j]


for i in range(len(T)):
	for j in range(len(T)):
		if i is not j:
			norm = 0
			for k in range(len(T)):
				norm += T[i][k]*T[j][k]

			if np.abs(norm) > 1e-12:
				print "U is not orthogonalized"
				raise Exception

print "\n\nTranspose of U:"

print "       ",
for i in range(nnodes):
	wch = str(i+1)
	printstring = ""
	printstring = printstring + wch + " "*(6-len(wch))
	print printstring,
print ""

printstring = "    "
printstring = printstring + "-------"*nnodes
print printstring

for i in range(nnodes):
	wch = str(i+1)
	print "   |"
	printstring = ""
	printstring = printstring + wch + " "*(2-len(wch))
	print printstring, "|",

	for j in range(nnodes):
		if np.abs(T[j][i]) < 1e-12:
				print("   -- "),
		else:
			if T[j][i] < 0:
				print(" %.2f"%(T[j][i])),
			else:
				print("  %.2f"%(T[j][i])),
	print ""



B_out_block = []
for i in range(nnodes):
	_B_out_block = []
	for j in range(nnodes):
		_B_out_block.append(0)
	B_out_block.append(_B_out_block)

for i in range(len(ISC_ISCset)):
	for j in range(len(ISC_ISCset)):
		if i is not j:
			for k in range(len(ISC_ISCset[i])):
				for l in range(len(ISC_ISCset[j])):
					Cm = ISC_ISCset[i][k]
					Cn = ISC_ISCset[j][l]
					for m in range(len(transverse_coord[Cm])):
						for n in range(len(transverse_coord[Cn])):
							B_out_block[transverse_coord[Cm][m]][transverse_coord[Cn][n]] = 1

for i in range(len(sgorbits)):
	for j in range(len(sgorbits), nnodes):
		B_out_block[i][j] = 1
		B_out_block[j][i] = 1

filename1 = "U.txt"
fp = file(filename1, "w")
for i in range(nnodes):
	for j in range(nnodes):
		if np.abs(T[i][j]) < 1e-12:
			fp.write("  0.00000"), 
		else:
			if T[i][j] < 0:
				fp.write(" %.5f"%(T[i][j])),
			else:
				fp.write("  %.5f"%(T[i][j])),
	fp.write("\n")
fp.close()

filename2 = "B.txt"
fp2 = file(filename2, "w")
for i in range(nnodes):
	for j in range(nnodes):
		if np.abs(B[i][j]) < 1e-12:
			fp2.write("  0.00000"),
		else:
			if B[i][j] < 0:
				fp2.write(" %.5f"%(B[i][j])),
			else:
				fp2.write("  %.5f"%(B[i][j])),
	fp2.write("\n")
fp2.close()

print "\n\n\n\nB:"
print "       ", 
for i in range(nnodes):
	wch = str(i+1)
	printstring = ""
	printstring = printstring + wch + " "*(6-len(wch))
	print printstring,
print ""

printstring = "    "
printstring = printstring + "-------"*nnodes
print printstring

for i in range(nnodes):
	wch = str(i+1)
	print "   |"
	printstring = ""
	printstring = printstring + wch + " "*(2-len(wch))
	print printstring, "|",

	for j in range(nnodes):
		if np.abs(B[i][j]) < 1e-12:
			if B_out_block[i][j] == 0:
				print("  ===="),
			elif B_out_block[i][j] == 1:
				print("   -- "),
		else:
			if B[i][j] < 0:
				print(" %.2f"%(B[i][j])),
			else:
				print("  %.2f"%(B[i][j])),
	print ""

#check the block form"

for i in range(nnodes):
	for j in range(nnodes):
		if B_out_block[i][j] == 1:
			if np.abs(B[i][j]) > 1e-12:
				print "*** Outside of block is nonzero! ***"
				print "*** Please check you input correct set of clusters and network ***"
				raise Exception

#check the diagonalization of each block for ISC"

for i in range(len(ISC_ISCset)):
	if len(ISC_ISCset[i]) == 1:
		Cm = ISC_ISCset[i][0]
		for j in range(len(transverse_coord[Cm])):
			for k in range(nnodes):
				if transverse_coord[Cm][j] != k:
					if np.abs(B[transverse_coord[Cm][j]][k]) > 1e-12:
						print "*** off diagonal element is nonzero ***"
						print "*** Please check you input correct set of clusters and network ***"
						raise Exception
					
					if np.abs(B[k][transverse_coord[Cm][j]]) > 1e-12:
						print "*** off diagonal element is nonzero ***"
						print "*** Please check you input correct set of clusters and network ***"
						raise Exception
