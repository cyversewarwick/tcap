import numpy as np
import scipy as sp
import scipy.stats as sps
import pandas as pd
import argparse
import multiprocessing as mp
import ctypes
import sys
import time

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='CSV expression file. Gene names in first column, time points in first row. Single profile (average your replicates).')
	parser.add_argument('--Pool', dest='pool', type=int, default=1, help='Number of processes to use when precomputing the Qian matrix and later on in AP clustering. Default: 1 (no parallelisation)')
	parser.add_argument('--Iterations', dest='iters', type=int, default=1000, help='Number of Affinity Propagation steps to take in the clustering. Default: 1000')
	parser.add_argument('--Converge', dest='conv', type=int, default=1000, help='If the cluster affinity of all genes remains unchanged for this many steps, the clustering is deemed completed ahead of time. Default: 100')
	parser.add_argument('--Lambda', dest='damp', type=float, default=0.9, help='Dampening of the A and R matrix values to prevent numerical oscillations. Default: 0.9')
	parser.add_argument('--SelfSim', dest='s', type=float, default=0, help='Self similarity score override in the Qian similarity matrix to prevent degenerate cases. If 0, then set to median of matrix. Default: 0')
	args = parser.parse_args()
	return args

def qian(X,Y):
	'''A function that computes the Qian similarity measure between two expression profiles.'''
	holderplus = np.zeros((len(X)+1,len(X)+1))
	holderminus = np.zeros((len(X)+1,len(X)+1))
	for i in np.arange(1,len(X)+1):
		for j in np.arange(1,len(X)+1):
			holderplus[i,j] = np.max([0, holderplus[i-1,j-1] + X[i-1]*Y[j-1]])
			holderminus[i,j] = np.max([0, holderminus[i-1,j-1] - X[i-1]*Y[j-1]])
	return np.max([np.max(holderplus),np.max(holderminus)])

def ap_qian(X,Y):
	'''A function that computes the AP-corrected Qian similarity measure between two expression profiles.
	Don't forget to set ddof=1 in zscoring the data or this will be useless.'''
	val = qian(X,Y)
	return val-len(X)+1

def single_psistar(i, psistar, input):
	'''Computes a single row of the psistar thing.'''
	for j in np.arange(i,input.shape[0]):
		if i==j:
			#self-similarity. will be dealt with later
			continue
		psistar[i,j] = ap_qian(input[i,:], input[j,:])
		#the matrix is symmetrical
		psistar[j,i] = psistar[i,j]
	return psistar

def pool_single_psistar(i):
	'''The same as above, wrapped in a parallel-happy manner'''
	#no need for a return as this is on a shared memory array
	single_psistar(i, psistar, input)

def _psistar_pool_init(_psistar, _input):
	global psistar, input
	psistar = _psistar
	input = _input

def get_psistar(input, args):
	'''A function that computes the psistar similarity matrix using the Qian measure.'''
	if args.pool > 1:
		#shared memory life
		psistar_base = mp.Array(ctypes.c_double, input.shape[0]*input.shape[0])
		psistar = np.ctypeslib.as_array(psistar_base.get_obj())
		psistar = psistar.reshape(input.shape[0], input.shape[0])
		p = mp.Pool(args.pool, _psistar_pool_init, (psistar, input))
		p.map(pool_single_psistar, np.arange(input.shape[0]))
	else:
		psistar = np.zeros((input.shape[0], input.shape[0]))
		for i in np.arange(input.shape[0]):
			psistar = single_psistar(i, psistar, input)
	#now, time for self-similarity
	if args.s == 0:
		#slot in the median psistar score
		mask = np.ones(psistar.shape, dtype=bool)
		np.fill_diagonal(mask, 0)
		s1 = np.median(psistar[mask])
		np.fill_diagonal(psistar, s1)
	else:
		#just put in the provided value
		np.fill_diagonal(psistar, args.s)
	return psistar

def single_r(i, a, r, psistar, args):
	'''Updates a single row of the R matrix'''
	for j in np.arange(r.shape[0]):
		old = r[i,j]
		holder = np.delete(a[i,:]+psistar[i,:],j)
		r[i,j] = psistar[i,j] - np.max(holder)
		#damping. screw first step skip, Matlab ain't having any of it
		r[i,j] = args.damp*old + (1-args.damp)*r[i,j]
	return r

def pool_single_r(i):
	'''The above, for parallel purposes'''
	#no need for a return as this is on a shared memory array
	single_r(i, a, r, psistar, args)

def single_a(i, a, r, psistar, args):
	'''Updates a single row of the A matrix'''	
	for j in np.arange(r.shape[0]):
		old = a[i,j]
		#for i==j this will point at the same index and that's okay per formulas
		holder = np.delete(r[:,j],(i,j))
		holder[holder<0] = 0
		if i==j:
			a[i,j] = np.sum(holder)
		else:
			a[i,j] = np.min((0, r[j,j]+np.sum(holder)))
		#damping
		a[i,j] = args.damp*old + (1-args.damp)*a[i,j]
	return a

def pool_single_a(i):
	'''The above, for parallel purposes'''
	#no need for a return as this is on a shared memory array
	single_a(i, a, r, psistar, args)

def _ar_pool_init(_a, _r, _psistar, _args):
	global a, r, psistar, args
	a = _a
	r = _r
	psistar = _psistar
	args = _args

def ap_step(psistar, a, r, args, p=None):
	'''A function that computes the updated a, r, e for a single AP step.
	
	Input:
		* psistar - Qian similarity matrix
		* a, r - as in the TCAP paper, from the previous iteration, for damping
		* args - command line arguments
		* p - pool, if parallel. None'd by default'''
	if args.pool > 1:
		#step one - r
		p.map(pool_single_r, np.arange(r.shape[0]))
		#step two - a
		p.map(pool_single_a, np.arange(r.shape[0]))
	else:
		#step one - r
		for i in np.arange(r.shape[0]):
			r = single_r(i, a, r, psistar, args)
		#step two - a
		for i in np.arange(r.shape[0]):
			a = single_a(i, a, r, psistar, args)
	#step three - e
	e = np.zeros(r.shape[0])
	for i in np.arange(r.shape[0]):
		e[i] = np.argmax(a[i,:]+r[i,:])
	#output
	return (a, r, e)

def tcap_cluster(input, args):
	'''A function that takes standardised data on input, along with command line arguments,
	and returns the e(i) vector that contains cluster centers upon either hitting convergence
	or maxing out the iteration count.
	
	Input:
		* input - the .values of the imported CSV file
		* args - the command line arguments'''	
	#first, we need the similarity matrix
	t0 = time.time()
	psistar = get_psistar(input, args)
	t1 = time.time()
	print(t1-t0)
	#allocate variables
	t0 = time.time()
	if args.pool > 1:
		#shared memory life
		a_base = mp.Array(ctypes.c_double, psistar.shape[0]*psistar.shape[0])
		r_base = mp.Array(ctypes.c_double, psistar.shape[0]*psistar.shape[0])
		a = np.ctypeslib.as_array(a_base.get_obj())
		a = a.reshape(psistar.shape[0], psistar.shape[0])
		r = np.ctypeslib.as_array(r_base.get_obj())
		r = r.reshape(psistar.shape[0], psistar.shape[0])
		#initialise one pool and use it always and forever
		p = mp.Pool(args.pool, _ar_pool_init, (a, r, psistar, args))
	else:
		p = None
		a = np.zeros(psistar.shape)
		r = np.zeros(psistar.shape)
	e = np.zeros(input.shape[0])
	#counter thing to assess premature convergence
	counter = 1
	#do steps
	for i in np.arange(args.iters):
		e_old = e
		(a, r, e) = ap_step(psistar, a, r, args, p)
		#assess convergence
		if list(e_old) == list(e):
			counter += 1
			if counter == args.conv:
				#we have hit premature convergence
				break
		else:
			#reset
			counter = 1
	#do as it says on the tin - return e
	t1 = time.time()
	print(t1-t0)
	return e

def main():
	#grab command line arguments
	args = parse_args()
	#grab the data and zscore it. this is not a may. this is needed for the Qian
	data = pd.read_csv(args.input, header=None, sep='\t')
	#data = pd.read_csv(args.input, header=0, index_col=0)
	data[:][:] = sps.mstats.zscore(data,axis=1,ddof=1)
	#run the clustering proper
	#t0 = time.clock()
	e = tcap_cluster(data.values, args)
	#t1 = time.clock()
	#print(t1-t0)
	#print(e)

if __name__ == "__main__":
	main()