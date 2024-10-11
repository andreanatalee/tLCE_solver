#!/usr/bin/sage -python
# -*- coding: utf8 -*-

CYELLOW = '\33[33m'
CRED = '\033[91m'
CEND = '\033[0m'

import sys
import argparse
from itertools import combinations
import multiprocessing as mp
import time
import random
from tqdm import tqdm


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
	parser.add_argument("-k", "--code_dimension", type=int, help="code dimension", required=True)
	parser.add_argument("-q", "--prime", type=int, help="Field characteristic", required=True)
	parser.add_argument("-t", "--samples_number", type=int, help="number of samples", required=True)
	parser.add_argument("-m", "--iterations", type=int, help="number of instances tested", required=True)

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

#---------------------------------

def sample_different_entries(n,l):
	list =[]
	a = [i for i in range(0, n)]
	for i in range(n-1, 0, -1):
		j = randint(0, i)
		tmp = a[i]
		a[i] = a[j]
		a[j] = tmp

	for i in range(0, l):
		list += [a[i]]

	return list

def sample_L_random(n,l):
	L = []
	a = sample_different_entries(n,l)
	b = sample_different_entries(n,l)
	for i in range(l):
		L += [[a[i],b[i]]]
  
	return L
        

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

def task_2(n,k,q,t,system,L):
	R,b = get_reduced_system(n,k,q,system,L)
	l = len(L)
	#d = (t*k - n + l)*(t*(n-k)-(t-1)*(n-l)) - (l - 1)
	#if (R.augment(b)).rank() <= (R.nrows() - d):			# enhanced Rouche-Capelli (needs testing)
	#	return L
	if R.rank() == (R.augment(b)).rank():					# classic Rouche-Capelli
		return 1
	return 0

#-------------------------------------------------

def test(n,k,q,t,m):
	success_count = 0
	l = n - k*t + 1
	d = (t*k - n + l)*(t*(n-k)-(t-1)*(n-l)) - (l - 1)

	prefix = identity_matrix(FiniteField(q), k)
	suffix = identity_matrix(FiniteField(q), n - k)

	Q = sample_monomial_matrix(n,q)
	# print(f'\nSecret monomial matrix, Q:\n{Q}\n')

	while (True):
		G_i = []				# Lists of codes
		G_i_prime = []
		#H_i_prime = []
		system = matrix(FiniteField(q), 0, n**2)
		# print(f'The samples are:\n')
		for i in range(t):
			M, M_ = oracle_call_lce(Q,n,k,q)
			G = identity_matrix(FiniteField(q), k).augment(M, subdivide=False)
			G_ = identity_matrix(FiniteField(q), k).augment(M_, subdivide=False)
			# print(f'\nG{i}:\n{G}\n')
			# print(f'\nG{i}^:\n{G_}\n')
			#H_ = (-M_.transpose()).augment(identity_matrix(GF(q), n-k))
			G_i += [G]
			G_i_prime += [G_]
			#H_i_prime += [H_]
			system = system.stack(get_linear_system(prefix, M, M_, suffix))
		if (system.rank() == system.nrows()):
			break

	for i in range(m):
		flag = 1
		while flag:
#			L = sample_L_random(n,l)
			list = sample_different_entries(n,l)
			L = [ [j,list[j]] for j in range(l) ]
			for pair in L:
				if not Q[pair[0],pair[1]]:
					flag = 0
		if task_2(n,k,q,t,system,L):
			success_count += 1
			#print(L)

	return success_count

def main(n,k,q,t,m):
	false_positives = []
	l = n - k*t + 1
	d = (t*k - n + l)*(t*(n-k)-(t-1)*(n-l)) - (l - 1)
 
	# Loop to run the algorithm 10 times
	for i in tqdm(range(10)):
		false_positives += [test(n,k,q,t,m)]

	# Calculate average number of false positives
	average = sum(false_positives) / 10

	# Print the results
	print(f"Average False Positives Found: {average}")
	print(f"Expected: {m/(q**d) - 0.}")

#---------------------------------------------------

if __name__ == '__main__':

	n = arguments(sys.argv[1:]).code_length
	k = arguments(sys.argv[1:]).code_dimension
	q = arguments(sys.argv[1:]).prime
	t = arguments(sys.argv[1:]).samples_number
	m = arguments(sys.argv[1:]).iterations
	main(n, k, q, t, m)
