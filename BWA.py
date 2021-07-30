#!/usr/bin/python3

#taking inputs as arguments:

import sys
import os
from subprocess import call
import subprocess
import argparse
import matplotlib.pyplot as plt


inputFile = sys.argv[1] 	 # 1st Argument in the terminal is used as our input file
outFile = sys.argv[2]		 # 2nd Argument in the terminal is used as our output file
query = sys.argv[3]		 # 3rd Argument in the terminal is used as our Query Seq that we find 						in our total seq
kmer_length = sys.argv[4]	 # 4th Argument in the terminal is used as our Kmer Length

#checking fasta extension:

if not inputFile.endswith('.fasta'):			#Try to rename the extension to .fasta
    print("Extension of the file is not fasta")
    sys.exit("Run the program again with correct file")

#checking fasta format:

f=open(inputFile, "r")
g=open(outFile, "a")
line=f.readlines()
x=line[0]
if not x[:1]=='>':		#Checks weather the file starts with > or not
    sys.exit("Format of file is not fasta")

#Summary of the genome:
print("")
print("File is checked and the Extension is in fasta")
print("The format is also in fasta")
seq = ' '
seqs = []
print("")
print("Calculating the total length of the genome ")
print("")
new_file = inputFile
with open(new_file, 'r') as f:
    for line in f:
        if line.startswith('>'):	#line was previously defined that read lines in the inputfile
            seqs.append(seq)
            seq = ' '
        else:
            seq += line.strip()
            dna_length = len(seq)
print("Length is: " + str(dna_length))	# Give us the total length of the DNA calculated from above loop
print("")
print("calculating no of As Ts Gs Cs..")
print("")
A = seq.count("A")	#it counts occurance of A in our variable seq and saves it in another variable A
print("Number of As=:" + str(A))
T = seq.count("T")	#it counts occurance of T in our variable seq and saves it in another variable T
print("Number of Ts=:" + str(T))
C = seq.count("C")	#it counts occurance of C in our variable seq and saves it in another variable C
print("Number of Cs=:" + str(C))
G = seq.count("G")	#it counts occurance of G in our variable seq and saves it in another variable G
print("Number of Gs=:" + str(G))
print("")
print("calculating % of As Ts Gs Cs..")
length = len(seq)	# This was our Dna seq Length
print("")
print("%age of A: " + str((A * 100 / length))) #Calculates percentage by (obtained/total)*100
print("%age of T: " + str((T * 100 / length))) #          //
print("%age of C: " + str((C * 100 / length))) #          //
print("%age of G: " + str((G * 100 / length))) #          //
print("")



print("Calculating Nucleotides Content") # we can also use the above variables that already stores the data of occurance 						of A, T, G, C to calculate the AT & GC content  
countAT=0
countGC=0
rem=0
with open(inputFile, 'r+') as infile:
    next(infile)	#Skips the first line that is the header line of our seq
    for x in infile.read():
        if x=="A" or x=="T":
            countAT+=1
        if x=="G" or x=="C":
            countGC+=1

g.write("AT count: ")	#Writes this data in the file
g.write(str(countAT))	#          //
g.write("\n")		#          //
g.write("GC count: ")	#          //
g.write(str(countGC))	#          //
g.write("\n")		#          //

def atgc_perc(): #This function gives us the AT-GC % comparison by showing us a pie chart graph
	activities = ['AT', 'GC'] 
  
	slices = [countAT , countGC]  
	colors = ['r', 'y'] 
  
	plt.pie(slices, labels = activities, colors=colors,  
        startangle=90, shadow = True, explode = (0, 0), 
        radius = 1.2, autopct = '%1.1f%%')  
	plt.legend()  
	plt.show() 


atgc_perc() #Calling the above function

#converting whole sequence into a string:

with open(inputFile, 'r') as myfile:
    next(myfile) # skips the header line 
    seq=myfile.read().replace('\n','')
    seq=seq.translate({ord(i): None for i in 'qwertyuiopasdfghjklzxcvbnm:;|\SDFHJKLZXVBNMQWERYUIOP< >,.?/{}1234567890@#$%^&*()_-+=!`~[]'}) #does not read the above letters to be converted into the string

#burrows wheeler transform
sequence=seq
seq=seq+"#" 
a=len(seq)
List1=[]

for x in range(a):
    seq=seq+seq[0]
    seq=seq[:0]+seq[1:]
    List1.append(seq)

List1.sort()
orig_str=""
jumb_str=""

for x in List1:
    orig_str=orig_str+x[0]
    jumb_str=jumb_str+x[-1]

#creating list with relative indexes

def process(jumb_str):
    word_dict={}
    jumb_list=[0]*len(jumb_str)

    for i,c in enumerate(jumb_str):
        if c not in word_dict:
            word_dict[c]=0
        else:
            word_dict[c]+=1
        jumb_list[i]=word_dict[c]
    return jumb_list

#burrows wheeler aligner

def reverse_search(orig_str,jumb_str,match_str,jumb_list):
    os_m_range=None

    for i in range(len(match_str)-1,-1,-1):
        match_value=match_str[i]

        if not os_m_range:
            start=orig_str.find(match_value)
            end=start+1
            while end<len(orig_str) and orig_str[end]==match_value:
                end=end+1
            os_m_range=(start,end)

        else:
            start=os_m_range[0]
            finish=os_m_range[1]
            while start<len(jumb_str) and jumb_str[start]!=match_value:
                start=start+1
            assert start<len(jumb_str), "Exception: Query doesn't match anywhere in the sequence"
            count=1
            end=start+1
            while end<finish:
                if jumb_str[end]==match_value:
                    count+=1
                end+=1
            js_m_range=(start,count)

            os_start=orig_str.find(match_value)+jumb_list[start]
            os_m_range=(os_start, os_start+count)
  
    return os_m_range

def character_position(str1,index):
    List=[]
    for x in range(index[0],index[1]):
        a=str1[x]
        count=0
        for i,c in enumerate(str1):
            if i<x:
                if c==a:
                    count+=1
        List.append(count)
    return List

def original_index(main_string,value1,a):
    List=[]
    for x in value1:
        count=0
        for i,c in enumerate(main_string):
            if c==a:
                if count==x:
                    List.append(i)
                count+=1
    return List

def kmer(Tseq, k):
    kfreq={}
    for i in range(0, len(Tseq) -k +1):
        kmer=Tseq[i:i+k]
        if kmer in kfreq:
            kfreq[kmer]+=1
        else:
            kfreq[kmer]=1
    return kfreq

k=int(kmer_length) #stores the input we gave as argument and saves it as a int in variable K
rf=kmer(sequence,k) #Calls the kmer function giving it 2 arguments as well 1= seq to be used 2= length of kmers
ListRF=[[rf[s], s] for s in rf] #saves the return of kmer seq into a list
ListRF.sort(reverse=True)     # reverse sorts the list
g.write("Kmer: ")	#saves the data in the file
g.write(str(ListRF))	#           //
g.write("\n")		#           //

Dikhao = len(ListRF)  #takes the list length and stores in variable to create a graph for kmers
numBins=[ i for i in range(Dikhao)]	#loop to set range acc to length of our list
bins=[name[0] for name in ListRF]
occurances = [name[1] for name in ListRF]
plt.bar(numBins, occurances, align='center')
plt.xticks(numBins, bins)
plt.xticks(rotation = 90)
plt.show()			#displays the graph on the screen

jumb_list = process(jumb_str)
match_string=query
z=match_string[0]
value=reverse_search(orig_str,jumb_str,match_string,jumb_list)
value1=character_position(orig_str,value)
value2=original_index(sequence,value1,z)
g.write("Query is found at Indexes:")
g.write(str(value2))

def prodigal():
	file = inputFile
	def Average(length):
		return sum(length)/ len(length)
	print("------------------------STEP 8------------------------")
	print('working on prodigal')
	seq=''
	seq_array=[]
	
	prodigal_file = file
	with open(prodigal_file,'r') as f: 
		for line in f:
			if line.startswith('>'):
				seq_array.append(seq)
				seq = '' 
			else:
				seq =seq+line.strip()
				dna_length = len(seq)
		print("length of genome is:" + str(dna_length))
	print("prediction of genes is in process of given fasta file")		
	pro ='prodigal -i ' + file + ' -a protein -c -d nucleotide -f gff -s data > output'
	prodigal = [pro]
	for pro in prodigal:
		subprocess.call(pro, shell = True)
	f=open('output','r')
	start=[]
	end=[]
	length=[]
	for line in f:
    		if line.startswith("#"):
        		continue
    		else:
        		new_line=line.split("\t")
        		new_line[3]=int(new_line[3])
        		start.append(new_line[3])
        		new_line[4]=int(new_line[4])
        		end.append(new_line[4])

#start end array for beginning and end of gene can be print by following code
#	print('start array:\n')
#	print(start)
#	print('end array:\n')
#	print(end)

	for i in range(0,len(start)):
   		 gene=(end[i]-start[i])
   		 length.append(gene)
	print('array with length of genes:\n')
	print(length)
	print('max length of genes:\n')
	print(max(length))
	print('min length of genes:\n')
	print(min(length))
	average=Average(length)
	print('average of genes:')
	print(round(average,2))
	print(' base pairs forming coding regions ')
	print(sum(length))
	seqq=''
	seqss=[]
	new_file = inputFile
	with open(new_file,'r') as f: 
		for line in f:
			if line.startswith('>'):
				seqss.append(seqq)
				seqq = '' 
			else:
				seqq +=line.strip()
				dna_length = len(seqq)
		#print("length is " + str(dna_length))

	print('base pairs forming non-coding regions')
	print(int(dna_length)-sum(length))


	print("% of coding region")
	cp=(float(sum(length)*100/int(dna_length)))
	print(cp)
	nc=100-cp
	print("% of non-coding region")
	print(nc)



	plt.title('GENE LENGTH BAR CHART')
	plt.ylabel("Number of Genes")
	plt.xlabel("Gene Length")
	plt.hist(length,bins=[100,200,300,400,500,600,700,800,900,1000,1100,1200,1300],rwidth=0.95,color='b')
	fname=input(str('Enter name for the graph\n'))

	plt.savefig(fname, transparent=False, bb0x_inches='tight')

	plt.show()
prodigal()
