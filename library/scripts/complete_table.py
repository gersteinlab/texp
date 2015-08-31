#!/usr/bin/python
import sys, getopt
from collections import defaultdict

def main(argv):
	
	file1 = ""
	file2 = ""
	infinite_defaultdict = lambda: defaultdict(infinite_defaultdict)

	try:
		opts, args = getopt.getopt(argv,"h1:2:",["list1=","list2="])	
	except getopt.GetoptError:
		print 'test.py -1 <file1> -2 <file2>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -1 <file1> -2 <file2>'
			sys.exit()
		elif opt in ("-1", "--list1"):
			file1 = arg
			name1 = file1.rsplit(".",1)[0]
		elif opt in ("-2", "--list2"):
			file2 = arg
			name2 = file2.rsplit(".",1)[0]
			
#LTR_Subfamily LTR_Ref_bases
#LTR7A 0.0000506129243799997
#LTR5 0.0000812283159704486
#LTR109A2 0.000122376759182875
#LTR22C 0.000127363421303424
	element_dict = []
	f = open(file1,"r")
	for line in f:
		tokens = line.split()
		element_dict.append(tokens[0])
	element_dict.pop(0)

#LTR06*LTR06 29.4
#LTR06*THE1D-int 0.01
#LTR1*LTR1 104.96
	means = defaultdict(dict)

	f2 = open(file2,"r")
	for line in f2:
		tokens = line.split()
		elements = tokens[0].split("*")
		means[elements[0]][elements[1]] = tokens[1]

	str_dict = defaultdict()
	str_title = "LTR_Subfamily"
	for key,value in means.iteritems():
		str_title += " "+key+"_Transcript"
		for element in element_dict:
			if ( element not in str_dict ):
				str_dict[element] = element

			if ( element not in  means[key].keys() ):
				str_dict[element]+=" 0"
			else:
				str_dict[element]+=" "+means[key][element]

	print str_title
	for element in element_dict:
		print str_dict[element]

def dummy():
	print dummy

if __name__ == "__main__":
	main(sys.argv[1:])

