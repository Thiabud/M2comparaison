# -*- coding: utf-8 -*-
#!/usr/bin/python3
import os
import sys
import subprocess
import shutil
import argparse
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import *
from collections import Counter
#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of genomic fasta
# - output filename
# - output sequence type dna/rna/protein

#GESTION DES ARGUMENTS 
parser = argparse.ArgumentParser(description='')
parser.add_argument("-r", "--ref", dest="ref" , type=str, help="alignement of reference") 
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

#Ouverture du GFF de reference et extraction de toute les ligne correspondant a des genes
fileref=open("ref.txt","w")
with open(ref) as reffile:
    for refline in reffile:
        if (refline.split('\t')[2])=="gene":
            fileref.write(refline)
fileref.close()
#Ouverture du new GFF et extraction de toute les ligne correspondant a des genes
filenew=open("new.txt","w")
with open(new) as newfile:
    for newline in newfile:
        if (newline.split('\t')[2]=="gene"):
            filenew.write(newline)
filenew.close()
allref=[]
allnew=[]

#recuperation de tout les identifiant pour 
with open("ref.txt") as reffile:
    for refline in reffile:
        id1=refline.split('\t')[8]
        allref.append(id1)

with open("new.txt") as newfile:
    for newline in newfile:
        #print(newline)
        #print(1)
        id2=newline.split('\t')[8]
        allnew.append(id2)
count2=0

#ANALYSE comparant les genes du gff de ref et du nouveau gff sortant un tableau  comprenannt :
#l'identifiant pour blast et mmseq, si le nouveau modele est compris, comprant ou cheavauche la ref, la taille du modele pour la ref et new, et le pourcentage de couverage entre les deux 
with open("ref.txt") as reffile:
    totalnew=0
    for refline in reffile:
        totalnew=0
        count+=1
        count2=0
        with open("new.txt") as newfile:
            deb1=int(refline.split('\t')[3])
            fin1=int(refline.split('\t')[4])
            id1=refline.split('\t')[8]
            taille1ref =fin1-deb1
            totalref+=taille1ref
            for newline in newfile:
                count2+=1
                if ((refline.split('\t')[0])==(newline.split('\t')[0])):
                    res=[]
                    id2=newline.split('\t')[8]
                    deb2=int(newline.split('\t')[3])
                    fin2=int(newline.split('\t')[4])
                    taille2mm = (fin2-deb2)
                    totalnew+=taille2mm
                    if (deb2<=deb1) :
                        res.append(0)
                    else :
                        res.append(1)
                    
                    if (fin2<=fin1) :
                        res.append(0)
                    else:
                        res.append(1)
                    if((res==[1,0])):
                        deb=deb2-deb1
                        fin=fin2-fin1
                        finalresult.append(["idblast",id1,"idmmseqs",id2,"compris","taille2mm",fin2-deb2,"taille1ref",fin1-deb1, "prct",(taille2mm/taille1ref)*100])
                    if (res==[0,1]):
                        finalresult.append(["idblast",id1,"idmmseqs",id2,"comprant","taille2mm",fin2-deb2,"taille1ref",fin1-deb1, "prct",(taille1ref/taille2mm)*100])
                    if (res==[0,0]):
                        if fin2>=deb1 :
                            deb=deb2-deb1
                            fin=fin2-fin1
                            finalresult.append(["idblast",id1,"idmmseqs",id2,"chevauchementgauche","taille2mm",fin2-deb2,"taille1ref",fin1-deb1,"couverture",((fin2-deb1)/taille1ref)*100])
                    if (res==[1,1]):
                        if deb2<=fin1 :
                            deb=deb2-deb1
                            fin=fin2-fin1
                            finalresult.append(["idblast",id1,"idmmseqs",id2,"chevauchementdroit","taille2mm",fin2-deb2,"taille1ref",fin1-deb1,"couverture",((fin1-deb2)/taille1ref)*100])
#Ecriture du tableau 
file = open('resultanalyse.csv', 'w+', newline ='\n') 
with file:     
    write = csv.writer(file) 
    write.writerows(finalresult) 
#Ecrire d'un tableau comprenant uniquement les modèles de gene different     
different=[]    
for element in finalresult:
    #print(element[10])
    if (element[10] != 100):
        different.append(element)
file = open('different.csv', 'w+', newline ='\n') 
with file:    
    write = csv.writer(file) 
    write.writerows(different) 
#Ecriture d'une liste permettant de faire la distribution du nombre de gene par pourcentage dde similitude par pas de 5 
#print(len(finalresult))
#print("totalref",totalref)
#print("totalnew",totalnew)
list999=[]
totalnumber=0
allcoverage=[]
distribution={}
for i in range(0,105,5):
    distribution[i]=0
for element in finalresult :
    allcoverage.append(element[10])
    totalnumber+=1
    for key in distribution.keys():
        if((element[10]>=key) and (element[10]<key+5)):
            distribution[key]+=1
    if (element[10]>=99.9):
        list999.append(element)
#print(distribution)
with open('mycsvfile.csv', 'w') as f: 
    w = csv.DictWriter(f, distribution.keys())
    w.writeheader()
    w.writerow(distribution)

#diagramm de venn des region 99,9% identique
identifier=[]
for element in finalresult :
    identifier.append(element[1])
uniqueidentifier= set(identifier)
numberofuniqueidentifier= len(uniqueidentifier)
print("______________________",numberofuniqueidentifier)
numberof999=len(list999)
print(numberof999)


venn2(subsets = (count-numberof999, count2-numberof999, numberof999), set_labels = ('tblastn', 'mmseqs2'), set_colors=('red', 'skyblue'), alpha = 0.7)
#venn2_circles(subsets = (count-numberof999, count2-numberof999, numberof999),linestyle='dashed', linewidth=2, color='k')
plt.show()


print("allref", len(allref))
print("allnew", len(allnew))
#creation d'une liste de tout les identifiant présent dans les résulats
inres=[]
for element in finalresult:
    inres.append(element[1])
#création d'une liste comprenant toute les regions d'interet du gff de ref non préssent dans le resultats final
refnotinres=[]    
for element in allref:
    if element not in inres:
       # print(element)
        refnotinres.append(element)
#création d'une liste comprenant toute les regions d'interet du gff new non préssent dans le resultats final
newnotinres=[]
print(len(inres))
print(len(allnew))
for element in allnew:
    if element not in inres:
        #print(element)
        newnotinres.append(element)
print("newnotinres",len(newnotinres))
print("refnotinres",len(refnotinres))

a_set = set(allref)
number_of_unique_values = len(a_set)
print("----------------------------",number_of_unique_values)



#### Extract all problematic model in gff format for each tool ###

#Extract from the ref gff the problematic annotation in a new gff
passed=False
bool2=False
with open("refproblem.gff", 'w+') as file :  
    with open(ref) as reffile:
        for line in reffile:
            for element in different:
                if element[1]==line.split('\t')[8]:
                    file.write(line)
                    passed=True  
                    bool2=False 
                if bool2==True and line.split('\t')[2] == "gene":
                    passed=False 
                    bool2==False        
                if line.split('\t')[2] != "gene" and passed==True:
                    file.write(line)
                    bool2=True
                    break

#Extract from the new gff the problematic annotation in a new gff
passed=False
bool2=False
with open("newproblem.gff", 'w+') as file :  
    with open(new) as newfile:
        for line in newfile:
            for element in different:
                if element[3]==line.split('\t')[8]:
                    file.write(line)
                    passed=True  
                    bool2=False 
                if bool2==True and line.split('\t')[2] == "gene":
                    passed=False 
                    bool2==False   
                if line.split('\t')[2] != "gene" and passed==True:
                    file.write(line)
                    bool2=True
                    break

###Extract one file for each chromosome##
shutil.rmtree('./chr/', ignore_errors=True)
#os.rmdir('./chr/')
os.mkdir('./chr/') 
with open("refproblem.gff") as file:
    for line in file:
        chr=line.split('\t')[0]
        with open('./chr/'+"ref_"+chr+"_problem.gff", "a") as f:
            #print(line)
            f.write(line)
with open("newproblem.gff") as file2:
    for line in file2:
        chr=line.split('\t')[0]
        with open('./chr/'+"new"+chr+"_problem.gff", "a") as f:
            #print(line)
            f.write(line)

print(totalnew)
print(totalref)
#print(distribution)

#count the exon in gff
numberofcdsinref="cut -f3 "+ ref +" | grep -c ""CDS"""
numberofcdsinref = subprocess.check_output(numberofcdsinref, shell=True)
numberofcdsinref=int(numberofcdsinref)
print("numberofcdsinref",numberofcdsinref)
numberofcdsinnew="cut -f3 "+ new +" | grep -c ""CDS"""
numberofcdsinnew = subprocess.check_output(numberofcdsinnew, shell=True)
numberofcdsinnew=int(numberofcdsinnew)
print("numberofcdsinnew",numberofcdsinnew)
nuclforef=0
with open(ref) as f:
    for line in f :
        if line.split('\t')[2]=="CDS":
            nuclforef+= ((int(line.split('\t')[4])-int(line.split('\t')[3])))
print("nuclforef",nuclforef)
nuclfornew=0
with open(new) as f:
    for line in f :
        if line.split('\t')[2]=="CDS":
            nuclfornew+= ((int(line.split('\t')[4])-int(line.split('\t')[3])))
print("nuclfornew",nuclfornew)