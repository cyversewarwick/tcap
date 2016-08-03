import argparse
import numpy as np
import pandas as pd
import sys
import math

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='CSV expression file. Gene names in first column, time points (\'Time\') in first or second row, condition (\'Treatment\') in first row if present.')
	parser.add_argument('--GP2S', dest='gp2s', type=argparse.FileType('r'), help='scores.txt as produced by GP2S')
	parser.add_argument('--GeneList', dest='list', type=argparse.FileType('r'), help='An optional gene list (such as TFs) to filter the expression data. One line per gene identifier.')
	parser.add_argument('--ScoreThreshold', dest='thresh', type=float, default=5, help='GP2S score significance threshold to deem a gene differentially expressed. Default: 5')
	parser.add_argument('--RemoveCondition', dest='cond', type=str, default=None, help='If provided, remove the data from a given condition (such as control in preparation for treated expression clustering). Default: None (no condition removal)')
	parser.add_argument('--AverageReps', dest='reps', action='store_true', help='Flag. If provided, the individual replicates will be averaged to return one value per time point')
	args = parser.parse_args()
	return args

def check_header(header):
	#for some reason, spacebars get imported as floaty nans because reasons
	if isinstance(header, float):
		if math.isnan(header):
			return True
	else:
		#if it's not a floaty nan, proceed as previously
		if header.upper() in ['TREATMENT','CONDITION','TIME','']:
			return True
	#is not header anymore
	return False

def load_data(args):
	data = pd.read_csv(args.input,index_col=None,header=None).values
	headind = []
	i=0
	while check_header(data[i,0]):
		headind.append(i)
		i+=1
	header = data[headind,:]
	data = np.delete(data,headind,axis=0)
	for i in range(data.shape[0]):
		for j in range(data.shape[1]):
			if j==0:
				data[i,j] = data[i,j].upper()
			else:
				data[i,j] = float(data[i,j])
	return (data, header)

def filter_gp2s(data, args):
	gp2s = pd.read_csv(args.gp2s,index_col=None,header=None,sep='\t').values
	for i in range(gp2s.shape[0]):
		gp2s[i,0] = gp2s[i,0].upper()
		gp2s[i,1] = float(gp2s[i,1])
	degs = gp2s[gp2s[:,1]>=args.thresh,0]
	del_inds = []
	for i in range(data.shape[0]):
		if data[i,0] not in degs:
			del_inds.append(i)
	if del_inds:
		data = np.delete(data,del_inds,axis=0)
	else:
		sys.stdout.write('No genes deemed not differentially expressed by GP2S.\n')
	return data

def filter_list(data,args):
	genes = args.list.readlines()
	for i in range(len(genes)):
		genes[i] = genes[i].strip().upper()
	del_inds = []
	for i in range(data.shape[0]):
		if data[i,0] not in genes:
			del_inds.append(i)
	if del_inds:
		data = np.delete(data,del_inds,axis=0)
	else:
		sys.stdout.write('No genes deleted when compared against the provided gene list.\n')
	return data

def filter_condition(data, header, args):
	header2 = list(header[0,:])
	args.cond = args.cond.upper()
	cond_del_ind = []
	for i in range(1,len(header2)):
		header2[i] = header2[i].upper()
		if header2[i] == args.cond:
			cond_del_ind.append(i)
	if cond_del_ind:
		header = np.delete(header,cond_del_ind,axis=1)
		data = np.delete(data,cond_del_ind,axis=1)
	else:
		sys.stdout.write('No instances of the provided condition name were found in the first line of the header.\n')
	#automatically kick the treatment line if applicable, i.e. if only one treatment left
	if len(np.unique(header[0,1:]))==1:
		sys.stdout.write('Dataset now features one condition only. Removing treatment line.\n')
		header = np.delete(header,0,axis=0)
	return (data, header)

def average_reps(data, header):
	data2 = np.delete(data,np.arange(1,data.shape[1]),axis=1)
	header2 = np.delete(header,np.arange(1,header.shape[1]),axis=1)
	free = [True]*(data.shape[1]-1)
	for i in range(1,data.shape[1]):
		if free[i-1]:
			#found a new thing that wasn't part of a rep group before
			minimask = [i]
			free[i-1] = False
			for j in range(i+1,data.shape[1]):
				if free[j-1] and all(header[:,i]==header[:,j]):
					#found something that's part of the same rep group
					minimask.append(j)
					free[j-1] = False
			header2 = np.hstack((header2, header[:,[i]]))
			newcol = np.mean(data[:,minimask],axis=1)
			newcol.shape = (newcol.shape[0],1)
			data2 = np.hstack((data2, newcol))
	return (data2, header2)

def write_data(data, header):
	#turn header into strings again (in case it's relevant)
	for i in range(header.shape[0]):
		for j in range(header.shape[1]):
			header[i,j] = str(header[i,j])
			#continue empty line thing
			if header[i,j] == 'nan':
				header[i,j] = ''
	#turn expression into strings again
	for i in range(data.shape[0]):
		for j in range(1,data.shape[1]):
			data[i,j] = str(data[i,j])
	output = np.vstack((header, data))
	toggle = ''
	with open('FilteredOutput.csv','w') as fid:
		for i in range(output.shape[0]):
			fid.write(toggle+','.join(output[i,:]))
			toggle='\n'

def main():
	#grab the arguments
	args = parse_args()
	#quick sanity check - do we actually have anything to do?
	if (not args.gp2s) and (not args.list) and (not args.cond) and (not args.reps):
		sys.stderr.write('No task given. No task performed. Aborting.\n')
		sys.exit(1)
	#if we're here, then there's at least something to do.
	#load data first
	(data, header) = load_data(args)
	#filter GP2S scores
	if args.gp2s:
		data = filter_gp2s(data, args)
	#filter gene list
	if args.list:
		data = filter_list(data, args)
	#filter out condition
	if args.cond:
		(data, header) = filter_condition(data, header, args)
	#average out replicates
	if args.reps:
		(data, header) = average_reps(data, header)
	#export the final thing
	write_data(data, header)
	
if __name__ == "__main__":
	main()