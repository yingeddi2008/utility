import sys
import os
import re
import argparse
from itertools import count,islice
from scipy.stats import poisson

def readVCF(vcffile):
	vlist = []
	with open(vcffile) as v:
		for line in v:
			if not line.startswith("#"):
				pos = line.rstrip().split("\t")[1]
				vlist.append(int(pos))
	return vlist

def generateCSV(vlist,genomelength,stepsize,windowsize,outcsv):
	try:
		os.remove(outcsv)
	except OSError:
		pass
	out=open(outcsv,"aw")
	ave=float(len(vlist))/genomelength*windowsize
	for i in islice(count(0),1,genomelength,stepsize):
		window = range(i,i+windowsize)
		cnt=sum([1 for x in vlist if x in window])
		lowerp=poisson.cdf(cnt,ave)
		upperp=1-lowerp
		out.write(','.join([str(i),str(i+windowsize-1),str(cnt),str(ave),str(lowerp),str(upperp)])+"\n")
	out.close()
	
def main():
	parser = argparse.ArgumentParser(description=
	''' Take an input of VCF file generated from GATK SNP analysis pipeline against a bacteria genome, and output the frequency of SNP occurrence and Poisson probability.''',
		epilog="No warranty comes with this script. Author: ying.eddi2008@gmail.com. \nAny suggestions or bugs report are welcomed.",add_help=False,formatter_class=argparse.RawTextHelpFormatter)
	##### Required arguments #####
	required = parser.add_argument_group('required arguments')
	required.add_argument("-i","--infile", help="input vcffile to be parsed",type=str,required=True)
	required.add_argument("-g","--glen", help="bacteria genome length",type=int,required=True)
	##### Optional arguments #####
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument("-o","--outfile",help="output file name [default: <infile>.csv]",type=str)
	optional.add_argument("-s","--ssize",default=500,help="step size [default: 500]",type=int)
	optional.add_argument("-w","--wsize",default=1000,help="windown size [default: 1000]",type=int)
	optional.add_argument("-h","--help",help="show this help message and exit",action="help")
	args = parser.parse_args()
	##### Output file name ####
	if not args.outfile:
		args.outfile=args.infile.replace(".vcf",".csv")
	##### Procedure #####
	SNPlist = readVCF(args.infile)
	generateCSV(SNPlist,args.glen,args.ssize,args.wsize,args.outfile)
	
		
if __name__ == "__main__":
	main()
	