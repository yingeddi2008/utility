import sys
import re
import itertools
import os
import glob
from itertools import chain
import argparse


def rline(lnin):
	ln=lnin.rstrip().strip(";")
	line=re.split("; |;",ln)
	return line

def bacteria_megan(meganout,pattern,level='g'):
	na=meganout.replace(pattern,"")
	unkn = []
	bacdict = {}
	with open(meganout) as m:
		for ln1, ln2 in itertools.izip_longest(*[m]*2):
			ln1=rline(ln1)
			ln2=rline(ln2)
			if len(ln1) < 3 and len(ln2) < 3:
				unkn.append(ln1[0])
				unkn.append(ln2[0])	
			if ln1[2::2] == ln2[2::2] and 'd__Bacteria' in ln2:
				tax = ln1[2::2]
				tmp = dict(x.split("__") for x in tax)
				if level in tmp and tmp[level] not in bacdict:
					bacdict[tmp[level]] = [ln1[0],ln2[0]]
				if level in tmp and tmp[level] in bacdict:
					bacdict[tmp[level]].append(ln1[0])
					bacdict[tmp[level]].append(ln2[0]) 
	###### write out file for unknowns #######
	try:
		os.remove(na+".unknown")
	except OSError:
		pass
	unknout = open(na+".unknown","w")
	unknout.write("\n".join(unkn)+"\n")
	unknout.close()
	###### write out bacteria file #######
	try:
            	os.remove(na+".bac")
        except OSError:
                pass
	bacout = open(na+".bac","w")
	for k,v in bacdict.items():
		for r in v:
			bacout.write("\t".join([na,k,r])+"\n") 
	bacout.close()
	return bacdict

def mergeDict(dict1, dict2):
	dict3 = defaultdict(list)
	for k, v in chain(dict1.items(), dict2.items()):
		dict3[k].append(v)
	return dict3

def writeBac(bacteriaName,meganMapList,minR=100):
	temp=[d[bacteriaName] for d in meganMapList if bacteriaName in d]
	read = [y for x in temp for y in x]
	if len(read) >= minR:
		try:
			os.remove(bacteriaName+".list")
		except OSError:
			pass
		###### write bacteria list #####
		baclist = open(bacteriaName+".list","w")
		baclist.write("\n".join(read)+"\n")
		baclist.close()
	
def main():
	parser = argparse.ArgumentParser(description=
	'''	A simple script to parse through several megan blast2lca outputs and extract bacteria reads names for downstream assembly programs. Outputs are:
		1. bacteria read list for each sample in three columns in the order of sampleName, BacteriaGenusName, and ReadName;
		2. unknown read list for each sample;
		3. abundant bacteria read list from all samples.''',
		epilog="No warranty comes with this script. Author: ying.eddi2008@gmail.com. \nAny suggestions or bugs report are welcomed.",add_help=False,formatter_class=argparse.RawTextHelpFormatter)
	##### Required arguments #####
	required = parser.add_argument_group('required arguments')
	required.add_argument("-p","--pattern", help="input file pattern to be parsed",type=str,required=True)
	##### Optional arguments #####
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument("-m","--minR",default=1000,help="minimum number of reads to write out in a bacteria taxonomy level [default: 1000]",type=int)
	optional.add_argument("-l","--level",default='g',help="at what taxonomy level to cluster, choices are g stands for genus, f stands for family, o stands for order, c stands for class [default: g]",choices=['g','f','o','c'],type=str)
	optional.add_argument("-h","--help",help="show this help message and exit",action="help")
	args = parser.parse_args()
	
	##### start ######
	samlist = glob.glob(args.pattern)
	if len(samlist) < 1:
		sys.exit("No samples to scan through, there might be a problem with your input pattern, please correct!")
	tor = map(bacteria_megan,samlist,itertools.repeat(args.pattern,len(samlist)),itertools.repeat(args.level,len(samlist)))
	bacnames = [d.keys() for d in tor ]
	fbac = list(set([y for x in bacnames for y in x]))
	map(writeBac,fbac,itertools.repeat(tor,len(fbac)),itertools.repeat(args.minR,len(fbac)))

if __name__ == "__main__":
	main()
