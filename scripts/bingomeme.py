import numpy as np

clusters = {}
keys = []

with open('clusters.txt','r') as fid:
	for line in fid:
		if line[0:3] == '---':
			#we have a new key
			holder = line.split()
			keys.append(holder[1])
			clusters[holder[1]] = []
		else:
			if not keys:
				#weird x at the start
				continue
			clusters[keys[-1]].append(line.strip())

#buffer with zeros to avoid dumb
bingokeys = []
digits = int(np.floor(np.log10(len(keys)))+1)
for key in keys:
	bingokeys.append(key.zfill(digits))

with open('functional_analysis_inputs/bingo.txt','w') as fid1:
	with open('functional_analysis_inputs/meme.txt','w') as fid2:
		for i in np.arange(len(keys)):
			fid1.write('>Cluster'+bingokeys[i]+'\n'+'\n'.join(clusters[keys[i]])+'\nbatch\n')
			holder = [keys[i]+'\t'+gene+'\n' for gene in clusters[keys[i]]]
			fid2.write(''.join(holder))