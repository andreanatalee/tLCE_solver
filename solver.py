#!/usr/bin/sage -python
# -*- coding: utf8 -*-

CYELLOW = '\33[33m'
CRED = '\033[91m'
CEND = '\033[0m'

import sys
import argparse
from itertools import combinations
import multiprocessing as mp


# SageMath imports
from sage.all import (
	randint,
    matrix,
    identity_matrix,
    random_matrix,
    zero_matrix,
    FiniteField,
    ZZ,
)

# -------------------------------
def arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-n", "--code_length", type=int, help="code length", required=True)
    parser.add_argument("-k", "--code_dimension", type=int, help="code dimension", required=False)
    parser.add_argument("-q", "--prime", type=int, help="Field characteristic", required=True)
    parser.add_argument("-t", "--samples_number", type=int, help="number of samples", required=True)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    options = parser.parse_args(args)
    return options

# -------------------------------
def is_monomial(input_matrix, n):
	for i in range(0, n, 1):
		if len(input_matrix.nonzero_positions_in_row(i)) != 1:
			return False
	return True

# ---------------------------------
def sample_invertible_matrix(k, q):
	while True:
		# Generate a random n x n matrix with entries in Zq
		A = random_matrix(FiniteField(q), k, k)
		# Check if the matrix is invertible
		if A.is_invertible():
			return A

# ----------------------------------
# sample permutation with fisher-yates
def sample_permutation_matrix(n, q):
	P = zero_matrix(FiniteField(q), n, n)
	a = [i for i in range(0, n)]
	for i in range(n-1, 0, -1):
		j = randint(0, i)
		tmp = a[i]
		a[i] = a[j]
		a[j] = tmp

	for i in range(0, n):
		P[i,a[i]] = 1

	return P


# -------------------------------
def sample_monomial_matrix(n, q):
	M = zero_matrix(FiniteField(q), n, n)
	for i in range(0, n):
		M[i,i] = randint(1, q-1)

	P = sample_permutation_matrix(n, q)
	return M * P


# ------------------------------
def oracle_call_lce(M, n, k, q):
	F = FiniteField(q)
	A_rand = random_matrix(F, k, n - k)
	A_code = identity_matrix(F, k).augment(A_rand)
	B_code = A_code * M
	W = B_code[0:k,0:k]
	while W.rank() < k :
		A_rand = random_matrix(F, k, n - k)
		A_code = identity_matrix(F, k).augment(A_rand)
		B_code = A_code * M
		W = B_code[0:k,0:k]
	FSB_code = (W.inverse())*B_code
	B = FSB_code[0:k,k:n]
	return [A_rand,B]

# --------------------------------------------------------------
def get_linear_system(prefix, g, g_, suffix, embed=lambda x: x):
	return embed(prefix.augment(g, subdivide=False)).tensor_product(embed((-g_.transpose()).augment(suffix, subdivide=False)), subdivide=False)


# ------------------------------------------------------------------------
def get_linear_system_transpose(prefix, g, g_, suffix, embed=lambda x: x):
	return embed((-g.transpose()).augment(suffix, subdivide=False)).tensor_product(embed(prefix.augment(g_, subdivide=False)), subdivide=False)


# -----------------------------------------------------
def get_reduced_system(n, k, q, system_, L):                # L is a list of guesses
	b = zero_matrix(FiniteField(q), system_.nrows(), 1)
	to_delete = []
	for l in range(len(L)):
		row = L[l][0]
		column = L[l][1]
		to_delete = to_delete + [(i * n + j) for i in range(0, n) for j in range(0, n) if ((i == row) or (j == column))] # xor

		if l == 0 :
			b = b + system_[:,(row * n + column)]
		if l > 0 :
			to_delete.remove(row * n + column)

	reduced_system = system_.delete_columns(to_delete)
	return reduced_system, b

#-------------------------------------------------

def task_2(n,k,q,t,system,L, survivals):
	R,b = get_reduced_system(n,k,q,system,L)
	l = len(L)
	#d = (t*k - n + l)*(t*(n-k)-(t-1)*(n-l)) - (l - 1)
	#if (R.augment(b)).rank() <= (R.nrows() - d):			# enhanced Rouche-Capelli (needs testing)
	#	return L
	if R.rank() == (R.augment(b)).rank():					# classic Rouche-Capelli
		return L
	return 0


def nested_loops(depth, current_depth=0,current_combination=[]):
	if current_depth == depth:
		L = [ [I[r],current_combination[r]] for r in range(len(W))]
		temp = task_2(n,k,q,t,system,L,survivals)
		if temp:
			survivals_temp += L
	for i in range(len(W[current_depth])):
		print(current_depth)
		nested_loops(depth, current_depth + 1, current_combination + [W[current_depth][i]])


def algorithm_2(n,k,q,t,system, I, survivals):			# Algorithm 2
	W = [ survivals[i] for i in I ]
	survivals_temp = []

	def nested_loops(depth, current_depth=0,current_combination=[]):		# subroutine for exploring all the compbinations
		if current_depth == depth:
			L = [ [I[r],current_combination[r]] for r in range(len(W))]
			# print(L)
			temp = task_2(n,k,q,t,system,L,survivals)
			if temp:
					for pair in L:
						survivals_temp.append(pair)
				#print(L)
			return
		for i in range(len(W[current_depth])):
			nested_loops(depth, current_depth + 1, current_combination + [W[current_depth][i]])

	nested_loops(len(I))
	# print(survivals_temp)

	# if len(I) == 1:
	for i in I:								# we leave only the variables that passed test 2
		to_remove = []
		for j in survivals[i]:
			if [i,j] in survivals_temp:
				h = 0
			else:
				to_remove += [j]
				# print([i,j])
		for h in to_remove:
			survivals[i].remove(h)
	return survivals

#	if len(I) > 1:
#		for i in I:								# we leave only the variables that passed test 2
#			to_remove = []
#			for j in survivals[i]:
#				if [i,j] in survivals_temp[0]:
#					h = 0
#				else:
#					to_remove += [j]
#					# print([i,j])
#			for h in to_remove:
#				survivals[i].remove(h)
#		return survivals


def count_vars(survivals):
	counter = 0
	for y in survivals:
		counter += len(y)
	return counter
		

def algorithm_3(n,k,q,t,system):
	survivals = [ [i for i in range(n)] for j in range(n) ]
	l = n - k*t + 1
	for i in range(n - l + 1):
		I = [i + h for h in range(l)]
		algorithm_2(n,k,q,t,system,I,survivals)
	return survivals


def survivals_to_entries(survivals):
	entries= []
	for i in range(n):
		for j in survivals[i]:
			entries += [ [i,j] ]
	return entries

#------------------------------------------------------------------------------------

def main(n,k,q,t):
	l = n - k*t + 1
	d = (t*k - n + l)*(t*(n-k)-(t-1)*(n-l)) - (l - 1)

	if l > n - k :
		print(f'Not enough samples. The algorithm needs more samples to recover the matrix.\n')
		return
	
	if n^l > q^d :
		print(f'The prime q is too small for the algorithm to be efficient.\n')
		return
	
	prefix = identity_matrix(FiniteField(q), k)
	suffix = identity_matrix(FiniteField(q), n - k)
    
	def test_algorithm(n,k,q):
		Q = sample_monomial_matrix(n,q)
		print(f'\nSecret monomial matrix, Q:\n{Q}\n')

		G_i = []				# Lists of codes
		G_i_prime = []
		#H_i_prime = []
		system = matrix(FiniteField(q), 0, n**2)
		print(f'The samples are:\n')
		for i in range(t):
			M, M_ = oracle_call_lce(Q,n,k,q)
			G = identity_matrix(FiniteField(q), k).augment(M, subdivide=False)
			G_ = identity_matrix(FiniteField(q), k).augment(M_, subdivide=False)
			print(f'\nG{i}:\n{G}\n')
			print(f'\nG{i}^:\n{G_}\n')
			#H_ = (-M_.transpose()).augment(identity_matrix(GF(q), n-k))
			G_i += [G]
			G_i_prime += [G_]
			#H_i_prime += [H_]
			system = system.stack(get_linear_system(prefix, M, M_, suffix))
		
		survivals = algorithm_3(n,k,q,t,system)
		if count_vars(survivals) <= system.rank():								# if the number of variables is less than the number of equation
			columns = [i*n + j for i in range(n) for j in survivals[i]]		    # we can recover the secret
			# print(survivals)
			# print(columns)
			vec = system[:, columns].right_kernel_matrix()[0]
			# print(system[:, columns].rank())
			# print(vec)
			Monomial = zero_matrix(FiniteField(q), n,n)
			entries = survivals_to_entries(survivals)
			# print(entries)
			for h in range(len(entries)):
				entry = entries[h]
				Monomial[entry[0],entry[1]] = vec[h]
			return Monomial, G_i, G_i_prime

	M, G_i, G_i_prime = test_algorithm(n,k,q)

	for i in range(t):
		if G_i_prime[i] != (G_i[i]*M).rref() :
			print(G_i_prime[i], (G_i_prime[i]*M).rref())
			print(f'\n Algorithm failed. The matrix found does not solve the instance\n')
			print(M)
			return
	
	print(f'The retrieved matrix is \n{M}')
	

#---------------------------------------------------

if __name__ == '__main__':

	n = arguments(sys.argv[1:]).code_length
	k = arguments(sys.argv[1:]).code_dimension
	q = arguments(sys.argv[1:]).prime
	t = arguments(sys.argv[1:]).samples_number
	main(n, k, q, t)	
