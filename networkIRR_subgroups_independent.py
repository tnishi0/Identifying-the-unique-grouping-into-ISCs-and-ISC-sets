from sage.matrix.constructor import *
#from sage.matrix.matrix import is_Matrix
import scipy as sp
import numpy.random as rn
import numpy.linalg as la
import numpy as np
import sage.all as sg

class GroupIRR:
	def __init__(self, Pgroup, adjmat):
		self._reset_data(Pgroup, adjmat)
		
	
	def _reset_data(self, Pgroup, adjmat): 
		self._adjmat=adjmat
		self._group=Pgroup
		self._nnodes=len(adjmat)
		self._orbits=None
		self._permutation_matrices=None
		self._conjugacy_classes=None
		self._conjugacy_classes_matrices=None
		self._character_table=None
		self._IRR_degeneracies=None
		self._projection_operators=None
		self._transformation_operator=None
		self._condition_check=None
		self._condition_check1=None
		self._condition_check2=None
	
	
	def get_permutation_matrices(self):
		if self._permutation_matrices==None:
			self._permutation_matrices=[]

		for element in self._group.list():
			self._permutation_matrices.append(np.array(element.matrix()))

		return list(self._permutation_matrices)

	
	def get_orbits(self):
		if self._orbits==None:
			self._orbits = self._group.orbits()

		return self._orbits

	
	def get_character_table(self):
		if self._character_table==None:
			self._character_table=self._group.character_table()

		return self._character_table


	def get_conjugacy_classes(self):
		if self._conjugacy_classes==None:
			self._conjugacy_classes=self._group.conjugacy_classes()

		return self._conjugacy_classes

	
	def get_conjugacy_classes_matrices(self):
		if self._conjugacy_classes==None:
			self.get_conjugacy_classes()

		if self._conjugacy_classes_matrices==None:
			self._conjugacy_classes_matrices=[]

			for conjclass in self._conjugacy_classes:
				sublist=[]
				clist = sg.ConjugacyClass(self._group, conjclass.representative()).list()

				for element in clist:
					sublist.append(np.array(element.matrix()))

				self._conjugacy_classes_matrices.append(sublist)

		return list(self._conjugacy_classes_matrices)
		

	def get_numIRRs(self):
		characters=self.get_character_table()
		numIRRs=len(characters[0])
		return numIRRs

					
	def get_IRR_degeneracies(self):
		if self._IRR_degeneracies==None:
			characters = self.get_character_table()
			numIRRs = len(characters[0])	
			group_order = self._group.order()
			
			self._IRR_degeneracies=[]
			matrices = self.get_conjugacy_classes_matrices()			

			for i in range(numIRRs):
				total=0
				for j in range(numIRRs):
					total = total + (float(len(matrices[j]))/group_order)*np.conj(np.complex(characters[i][j]))*np.trace(matrices[j][0])

				self._IRR_degeneracies.append(round(np.real(total)))

		return list(self._IRR_degeneracies)

		
	def get_projection_operator(self, j):
		degen=self.get_IRR_degeneracies()
		if self._projection_operators==None:
			self._projection_operators=[None]*len(degen)

		if self._projection_operators[j] == None:
			IRR_dimension=self.get_character_table()[j,0]
			group_order=self._group.order()
			characters=self.get_character_table()[j]
			matrices=self.get_conjugacy_classes_matrices()
			result=np.zeros((self._nnodes, self._nnodes))

			for i in range(len(characters)):
				for mat in matrices[i]:
					result=result+mat*np.conj(np.complex(characters[i]))

			self._projection_operators[j]=result*np.float(IRR_dimension)/np.float(group_order)

		return self._projection_operators[j]

	
	def get_transformation_operator(self):
		epsilon=1e-12
		if self._transformation_operator==None:
			result=[]
			degens=self.get_IRR_degeneracies()
			total=0
			for j in range(len(degens)):
				IRR_dimension = self.get_character_table()[j,0]
				
				if degens[j] > 1e-3:
					P = self.get_projection_operator(j)
					VR = la.eig(P)				
					
					eigenvectors = VR[1]
					eigenvectors = eigenvectors.transpose()

					R1 = int(IRR_dimension)*degens[j]
					R2 = 0
					for i, w in enumerate(VR[0]):
						if np.abs(w-1.0)<epsilon:
							R2 = R2+1
							result.append(eigenvectors[i])

#					U, W, VH = la.svd(P)
#					print "U", U, "W", W, "VH", VH
#					R1 = int(IRR_dimension)*degens[j]
#					R2 = 0
#					total=total+R1
#					print j,"-th projector eigenvectors"
#					for i, w in enumerate(W):
#						if np.abs(w-1.0)<epsilon:
#							R2 = R2+1
#							result.append(VH[i])
#							print VH[i],

					if R1!=R2:
						print "Warning!"
						print "Found, ",R2," singular vectors"
						print "there should be", R1
						raise Exception

			self._transformation_operator=np.array(result)
		return self._transformation_operator.copy()
	
	
	def condition_check(self):
		dimension_check = 0
		for i in range(len(self.get_IRR_degeneracies())):
			dimension_check += int(self.get_character_table()[i,0])*self.get_IRR_degeneracies()[i]

		if dimension_check != self._nnodes:
			print "dimension_error!\n"
			raise Exception

		if self._condition_check==None:
			characters=self.get_character_table()
			num_trivial_IRR = 0
			for i in range(len(characters[0])):				
				_characters = []
				for j in range(len(characters[0])):
					_characters.append(complex(characters[i,j]))

				_characters = np.array(_characters)
				if(max(np.abs(_characters-1.0))<1e-10):
					trivial_IRR = i
					num_trivial_IRR += 1
			
			if num_trivial_IRR > 1 or num_trivial_IRR < 1:
				print "# trivial_IRR is not 1"
				raise Exception

#			print "trivial_IRR", trivial_IRR

			degens=self.get_IRR_degeneracies()
			self._condition_check1 = True
			for i in range(len(degens)): # for all IRRs
				if i is not trivial_IRR and degens[i] > 1.5:
					self._condition_check1 = False
					break

			if self._condition_check1 == True:
				self._condition_check = True
				self._condition_check2 = False
				return self._condition_check, self._condition_check1, self._condition_check2

			# If there exists at least one IRR whose multiplicity is larger than one
			tmatrix = [] 
			dim_sub = [] # dimension of subspace for IRR
			for i in range(len(degens)):
				P = self.get_projection_operator(i)
				epsilon = 1e-12
				IRR_dimension = self.get_character_table()[i,0]
				R1 = int(IRR_dimension)*degens[i]
				R2 = 0
				
				if(degens[i] > 0.1):
					VR = la.eig(P)
					eigenvectors = VR[1]
					eigenvectors = eigenvectors.transpose()
					for j, w in enumerate(VR[0]):
						if np.abs(w-1.0)<epsilon:
							R2 += 1
							tmatrix.append(eigenvectors[j])
						
				if R1 != R2:
					print "number of eigenvectors error in",i,"-th IRR"
					raise Exception
				
				dim_sub.append(R2)
			
			Inv_tmatrix = la.inv(tmatrix)
			bmatrix = np.dot(tmatrix, np.dot(self._adjmat, Inv_tmatrix))
			check_row = 0
			for i in range(len(degens)):
				if i != trivial_IRR:
					eigenvalues = [] 
					neigenvalues = 0
					for j in range(check_row, check_row + dim_sub[i]):
						eigenvalues.append(bmatrix[j][j])
						neigenvalues += 1

					maxdeveigenvalue = 0
					for j in range(neigenvalues):
						for k in range(neigenvalues):
							if j != k:
								if np.abs(eigenvalues[j]-eigenvalues[k]) > maxdeveigenvalue:
									maxdeveigenvalue = np.abs(eigenvalues[j] - eigenvalues[k])

					if maxdeveigenvalue > 1e-8:
						self._condition_check = False
						self._condition_check2 = False	
						return self._condition_check, self._condition_check1, self._condition_check2
					
					for j in range(check_row, check_row + dim_sub[i]):
						for k in range(check_row, check_row + dim_sub[i]):			
							if j != k and np.abs(bmatrix[j][k]) > 1e-10:
								self._condition_check = False
								self._condition_check2 = False
								return self._condition_check, self._condition_check1, self._condition_check2
					
				check_row += dim_sub[i]	
			
			if check_row != self._nnodes:
				print "condition2 is not checked exactly\n"
				raise Exception
								 		
			self._condition_check = True
			self._condition_check2 = True
			return self._condition_check, self._condition_check1, self._condition_check2		
												

class NetworkIRR:
	def __init__(self, adjmat=None):
		self._reset_data()
		self._set_adjacency_matrix(adjmat)


	def _reset_data(self):
		self._adjmat =None
		self._nnodes=0
		self._group=None
		self._subgroups=None
		self._orbits=None
		self._permutation_matrices=None
		self._conjugacy_classes=None
		self._conjugacy_classes_matrices=None
		self._character_table=None
		self._IRR_degeneracies=None
		self._projection_operators=None
		self._transformation_operator=None
	
	
	def _set_adjacency_matrix(self,adjmat):
		self._reset_data();
#		self._adjmat=np.array(adjmat.copy())
		self._adjmat=matrix(adjmat.copy())
		self._nnodes=self._adjmat.nrows()
	    
	
	
	def get_adjacency_matrix(self):
		#return self._adjmat.copy() 
		return self._adjmat

	
	
	def get_automorphism_group(self):
		if self._group==None:
			self._group,self._orbits = (sg.Graph(sg.Matrix(self._adjmat))).automorphism_group(orbits=True)
 
 		return self._group       
#		return sg.copy(self._group)
        
	
	
	def get_automorphism_group_matrices(self):
		if self._group==None:
			self.get_automorphism_group()
        
		if self._permutation_matrices==None:
			self._permutation_matrices=[]
            
		for element in self._group.list():
			self._permutation_matrices.append(np.array(element.matrix()))
                
		return list(self._permutation_matrices)
        
	

	def get_automorphism_group_subgroups(self):
		if self._group==None:
			self.get_automorphism_group()
			
		if self._subgroups==None:
			self._subgroups = self._group.subgroups()	

		return self._subgroups
	
	def get_orbits(self):
		if self._orbits==None:
			self._group, self._orbits = sg.Graph(self.get_adjacency_matrix()).automorphism_group(orbits=True)

		return sg.copy(self._orbits)

	
	
	def get_character_table(self):
		if self._character_table==None:
			self._character_table=self.get_automorphism_group().character_table()

		return sg.copy(self._character_table)
    
	
	
	def get_conjugacy_classes(self):
		if self._conjugacy_classes==None:
			if self._group==None:
				self.get_automorphism_group()
                
			self._conjugacy_classes=self._group.conjugacy_classes()
            
		return sg.copy(self._conjugacy_classes)
        
        
        
	def get_conjugacy_classes_matrices(self):
		if self._conjugacy_classes==None:
			self.get_conjugacy_classes()
        
		if self._conjugacy_classes_matrices==None:
			self._conjugacy_classes_matrices=[]
            
			for conjclass in self._conjugacy_classes:
				sublist=[]
                #This line makes no sense, but conjclass.list()
                #returns an error
				clist=sg.ConjugacyClass(self._group,conjclass.representative()).list()
                
				for element in clist:
					sublist.append(np.array(element.matrix()))
                    
				self._conjugacy_classes_matrices.append(sublist)
            
		return list(self._conjugacy_classes_matrices)

	
	
	def get_numIRRs(self):
		characters=self.get_character_table()
		numIRRs=len(characters[0])
		return numIRRs
		
	
	
	def get_IRR_degeneracies(self):
		if self._IRR_degeneracies==None:
			characters=self.get_character_table()
			numIRRs=len(characters[0])
			group_order=self.get_automorphism_group().order()
            
			self._IRR_degeneracies=[]
			matricies=self.get_conjugacy_classes_matrices()
			
			for i in range(numIRRs): #Loop over IRRs
				total=0
				for j in range(numIRRs): #loop over classes
					total=total +(float(len(matricies[j]))/group_order)* np.conj(np.complex(characters[i][j]))* np.trace(matricies[j][0])
#					print float(len(matricies[j])), group_order, np.conj(np.complex(characters[i][j])), np.trace(matricies[j][0])
					
				self._IRR_degeneracies.append(round(np.real(total)))
                
		return list(self._IRR_degeneracies)
        
	
	
	def get_projection_operator(self, j):
		degen=self.get_IRR_degeneracies()
		if self._projection_operators==None:
			self._projection_operators=[None]*len(degen)
        
		if self._projection_operators[j] == None:
            #the dimension of the IRR is the character of the identity.
			IRR_dimension=self.get_character_table()[j,0]
			group_order=self.get_automorphism_group().order()
			characters=self.get_character_table()[j]
			matricies=self.get_conjugacy_classes_matrices()
			result=np.zeros((self._nnodes,self._nnodes))
            
			for i in range(len(characters)):
				for mat in matricies[i]:
					result=result+mat*complex(characters[i])
                    
			self._projection_operators[j]=result*np.float(IRR_dimension)/np.float(group_order)
        
		return self._projection_operators[j]

	
	
	def get_transformation_operator(self):
		epsilon=1e-12
		if self._transformation_operator==None:
			result=[]
			degens=self.get_IRR_degeneracies()
			total=0
			for j in range(len(degens)): #for all kinds of irreducible representations
				IRR_dimension=self.get_character_table()[j,0]
				
				if degens[j] > 1e-3:
					P=self.get_projection_operator(j)
					U,W,VH = la.svd(P)
                    
					#print 'U', U, 'W', W, 'VH', VH
                    
                    #Uncomment this output if you want to see it.
                    #print "Representation ", j, " with degeneracy ", degens[j]
                    #print "W=",W
                    #print "dimension=",IRR_dimension
                    
					R1=int(IRR_dimension)*degens[j] #expected number of IRR coordinates of $j$-th irreducible representation
					R2=0
					total=total+R1
					for i,w in enumerate(W):
						if np.abs(w-1.0)<epsilon:
							R2=R2+1
							result.append(VH[i]) #list result forms a transformation matrix T
                        
					if R1!=R2:
						print "Warning!"
						print "Found, ",R2," singular vectors"
						print "there should be", R1
						raise Exception
                        
                    #print "Rank=",R1,R2
                    #print "P=",np.real(P)
                    
            #print "sum of dimension*degeneracy", total
			self._transformation_operator=np.array(result)
		return self._transformation_operator.copy()
