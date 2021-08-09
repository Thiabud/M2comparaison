# -*- coding: utf-8 -*-
#!/usr/bin/python3


import argparse
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of genomic fasta
# - output filename
# - output sequence type dna/rna/protein

#GESTION DES ARGUMENTS 
parser = argparse.ArgumentParser(description='')
parser.add_argument("-r", "--ref", dest="ref" , type=str, help="sequence fasta d'exon") 
parser.add_argument("-n", "--new", dest="new" , type=str, help="new alignement") 
#parser.add_argument("-o", "--out", dest="out_repo" , type=str, help="fichier de sorti") 
args = parser.parse_args()
new = args.new
ref = args.ref
finalresult= []
present=[]
totalref=0
totalnew=0
count=0
tmp=0
exon=[]
lineoftmp=""
with open(ref) as reffile:
    for refline in reffile:
        if refline[0]=='>' :
            longeur=int(len(refline.split('_'))-1)
            numberofcds=int(refline.split('_')[longeur][0])
            #print(numberofcds)
            #print(tmp)
            if (int(tmp) >= numberofcds):
                tmp = refline.split('_')[longeur][0]
                lineoftmp=refline
            else:
                exon.append(lineoftmp[1:-1])
                tmp = int(refline.split('_')[longeur][0])
                lineoftmp=refline
exon = exon[1:]
print("exon",exon)
number =0
gff=[]
gene=""
with open(new) as newfile:
    for newline in newfile :
        if (newline.split('\t')[2] == "gene"):
            #print(newline.split('\t')[2])
            gff.append(gene+":cds_"+str(number))
            number=0
            gene=(newline.split('=')[2]).split(' ')[0]
        else :
            number+=1
print(gff)
plusdansgff = 0
plusdansgfflist = {}
plusdansgfftotal=0

plusdansexon = 0
plusdansexonlist = {}
plusdansexontotal=0

pareille = 0
pareillelist = {}
pareillenumber=0
for lastexon in exon:
    posexon=lastexon.split('_')[1]
    for lastgff in gff :
        posgff=lastgff.split('_')[1]
        if posexon == posgff:
            longeur=int(len(lastexon.split('_')))-1
            exonnumber = int(lastexon.split('_')[longeur])
            gffnumber = int(lastgff.split('_')[longeur])
            if gffnumber < exonnumber :
                plusdansexon +=1
                plusdansexonlist[lastexon] = exonnumber-gffnumber
                plusdansexontotal += exonnumber-gffnumber
            if gffnumber > exonnumber :
                plusdansgff += 1
                plusdansgfflist[lastexon] = gffnumber- exonnumber
                plusdansgfftotal += gffnumber - exonnumber
            if gffnumber == exonnumber :
                pareille +=1
                pareillelist[lastexon] = exonnumber-gffnumber
                pareillenumber += exonnumber-gffnumber
print("pareille",pareille)
print("pareillelist",pareillelist)
print("pareillenumber",pareillenumber)

print("plusdansgff",plusdansgff)
print("plusdansgfflist",plusdansgfflist)
print("plusdansgffnumber",plusdansgfftotal)

print("plusdansexon",plusdansexon)
print("plusdansexonlist",plusdansexonlist)
print("plusdansexontotal",plusdansexontotal)