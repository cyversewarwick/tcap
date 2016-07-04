#iPlant compatibility
import matplotlib
matplotlib.use('Agg')

import numpy as np
import scipy.stats as spstats
import pandas as pd
import argparse
import multiprocessing as mp
import ctypes
import sys
import os
import time
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='CSV expression file. Gene names in first column, time points in first row. Single profile (average your replicates).')
	parser.add_argument('--Pool', dest='pool', type=int, default=0, help='Number of processes to use when precomputing the Qian matrix. Default: 0 (automated parallelisation)')
	parser.add_argument('--Iterations', dest='iters', type=int, default=1000, help='Number of Affinity Propagation steps to take in the clustering. Default: 1000')
	parser.add_argument('--Converge', dest='conv', type=int, default=100, help='If the cluster affinity of all genes remains unchanged for this many steps, the clustering is deemed completed ahead of time. Default: 100')
	parser.add_argument('--Lambda', dest='damp', type=float, default=0.9, help='Dampening of the A and R matrix values to prevent numerical oscillations. Default: 0.9')
	parser.add_argument('--SelfSim', dest='s', type=float, default=0, help='Self similarity score override in the Qian similarity matrix to prevent degenerate cases. If 0, then set to median of matrix. Default: 0')
	parser.add_argument('--NoHL', dest='hl', action='store_false', help='Flag. If provided, the centroid expression profile will not be highlighted in red on the plots.')
	args = parser.parse_args()
	return args

def print_time(t1, t0, jobname):
	'''Helper function to print the time taken by procedures.'''
	m, s = divmod(round(t1-t0), 60)
	h, m = divmod(m, 60)
	holderstring = 'Took %d:%02d:%02d. '+jobname+' complete.\n'
	sys.stdout.write(holderstring % (h, m, s))

def qian(X,Y):
	'''A function that computes the AP-Qian similarity measure between two expression profiles.'''
	holderplus = np.zeros((len(X)+1,len(X)+1))
	holderminus = np.zeros((len(X)+1,len(X)+1))
	out = 0
	#range() is 2x faster than np.arange() apparently
	for i in range(len(X)):
		for j in range(len(X)):
			#don't ask. this is faster. np.max() is glacial
			h1 = holderplus[i,j] + X[i]*Y[j]
			h2 = holderminus[i,j] - X[i]*Y[j]
			if h1>0:
				holderplus[i+1,j+1] = h1
				if h1>out:
					out=h1
			if h2>0:
				holderminus[i+1,j+1] = h2
				if h2>out:
					out=h2
	return out-len(X)+1

def single_psistar(i, psistar, data):
	'''Computes a single gene's worth of Qian similarity matrix'''
	for j in np.arange(i+1, data.shape[0]):
		val = qian(data[i,:], data[j,:])
		psistar[i,j] = val
		psistar[j,i] = val
	#no need for a return as this is done on a shared memory array

def pool_single_psistar(i):
	'''The same as above, wrapped in a parallel-happy manner'''
	#no need for a return as this is on a shared memory array
	single_psistar(i, psistar, data)

def _psistar_pool_init(_psistar, _data):
	global psistar, data
	psistar = _psistar
	data = _data

def get_psistar(data, args):
	'''Computes the AP-Qian similarity matrix'''
	#set up shared memory regardless of parallelisation or not
	psistar_base = mp.Array(ctypes.c_double, data.shape[0]*data.shape[0])
	psistar = np.ctypeslib.as_array(psistar_base.get_obj())
	psistar = psistar.reshape(data.shape[0], data.shape[0])
	if args.pool > 1:
		p = mp.Pool(args.pool, _psistar_pool_init, (psistar, data))
		p.map(pool_single_psistar, np.arange(data.shape[0]))
	else:
		for i in np.arange(data.shape[0]-1):
			single_psistar(i, psistar, data)
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

def rstep(a,r,psistar,args):
	'''A full vectorised single R step'''
	old = np.copy(r)
	apsi = a+psistar
	ind1 = np.arange(apsi.shape[0])
	ind2 = np.argmax(apsi,axis=1)
	val1 = apsi[ind1, ind2]
	apsi[ind1, ind2] = float('-inf')
	val2 = np.max(apsi, axis=1)
	r = psistar - val1[:,None]
	r[ind1, ind2] = psistar[ind1, ind2] - val2
	r = args.damp*old + (1-args.damp)*r
	return r

def astep(a,r,psistar,args):
	'''A full vectorised single A step'''
	old = np.copy(a)
	rp = np.copy(r)
	rp[rp<0] = 0
	np.fill_diagonal(rp, r.diagonal())
	a = np.sum(rp,axis=0)[None,:] - rp
	ad = np.copy(a.diagonal())
	a[a>0] = 0
	np.fill_diagonal(a, ad)
	a = args.damp*old + (1-args.damp)*a
	return a

def tcap_cluster(data, args):
	'''A function that takes standardised data on input, along with command line arguments,
	and returns the e(i) vector that contains cluster centers upon either hitting convergence
	or maxing out the iteration count.
	
	Input:
		* input - the .values of the imported CSV file
		* args - the command line arguments'''	
	#first, we need the similarity matrix
	sys.stdout.write('Commencing Qian similarity matrix computation...\n')
	t0 = time.time()
	psistar = get_psistar(data, args)
	t1 = time.time()
	print_time(t1, t0, 'Similarity matrix computation')
	#allocate variables
	sys.stdout.write('Commencing Affinity Propagation clustering...\n')
	t0 = time.time()
	r = np.zeros(psistar.shape)
	a = np.zeros(psistar.shape)
	#well, this one's not shared
	e = np.zeros(data.shape[0])
	#counter thing to assess premature convergence
	counter = 1
	for i in np.arange(args.iters):
		e_old = e
		#step one - r
		r=rstep(a,r,psistar,args)
		#step two - a
		a=astep(a,r,psistar,args)
		#assess convergence
		e = np.argmax(a+r,axis=1)
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
	print_time(t1, t0, 'Clustering')
	return e

def make_output(data, e, args):
	'''Create all the relevant plots/exports.
	
	Input:
		* data - the parsed CSV file
		* e - the resulting clustering vector from tcap_cluster()
		* args - the command line arguments'''
	#get gene list and time points
	genes = np.array(data.index)
	tps = np.array([float(x) for x in data.columns])
	#get centroids
	centroids = np.unique(e)
	#set up plot output folder
	if not os.path.exists('plots'):
		os.makedirs('plots')
	#cluster output file
	writer = open('clusters.txt','w')
	for i in range(len(centroids)):
		clusgenes = genes[e==centroids[i]]
		writer.write('---Cluster '+str(i+1)+' ---\n'+'\n'.join(clusgenes)+'\n')
		plotname = 'plots/Cluster'+str(i+1).zfill(int(np.floor(np.log10(len(centroids)))+1))+'.eps'
		plt.figure(figsize=(10,6))
		with sns.axes_style("ticks"):
			for gene in clusgenes:
				plt.plot(tps,data.loc[gene,:],color='#03012d',linewidth=1)
			if args.hl:
				centroid = genes[centroids[i]]
				plt.suptitle('Cluster '+str(i+1)+ ' (Centroid '+centroid+')')
				plt.plot(tps,data.loc[centroid,:],color='#d9544d',linewidth=2)
			else:
				plt.suptitle('Cluster '+str(i+1))
			#this is large enough to support viewing all the time points most of the time
			if len(tps)<30:
				plt.xticks(tps)
			plt.xlim([min(tps),max(tps)])
			plt.xlabel('Time')
			plt.ylabel('Standardised Expression')
			plt.tight_layout()
			sns.despine()
			plt.savefig(os.path.normcase(plotname), format='eps', dpi=1000, bbox_inches='tight')
	writer.close()

def main():
	#grab command line arguments
	args = parse_args()
	#assess parallelisation
	if args.pool==0:
		args.pool=mp.cpu_count()
	#grab the data and zscore it. this is not a may. this is needed for the Qian
	#data = pd.read_csv(args.input, header=None, sep='\t')
	sys.stdout.write('Parsing expression file...\n')
	data = pd.read_csv(args.input, header=0, index_col=0)
	data[:][:] = spstats.mstats.zscore(data,axis=1,ddof=1)
	#run the clustering proper
	e = tcap_cluster(data.values, args)
	#create output
	make_output(data, e, args)

if __name__ == "__main__":
	main()