# -*- coding: utf-8 -*-

import argparse
import csv
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
with open(ref) as reffile:
    for refline in reffile:
        with open(new) as newfile:
            deb1=int(refline.split('\t')[3])
            fin1=int(refline.split('\t')[4])
            id1=refline.split('\t')[8]
            taille1ref =fin1-deb1
            totalref+=taille1ref
            totalnew=0
            for newline in newfile:
                if ((refline.split('\t')[0])==(newline.split('\t')[0])):
                    res=[]
                    deb2=int(newline.split('\t')[3])
                    print("deb",deb2)
                    fin2=int(newline.split('\t')[4])
                    print("fin",fin2)
                    taille2mm = (fin2-deb2)
                    print("taille2mm",taille2mm)
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
                        finalresult.append([id1,"compris","taille2mm",fin2-deb2,"taille1ref",fin1-deb1, "prct",(taille2mm/taille1ref)*100])
                    if (res==[0,1]):
                        finalresult.append([id1,"comprant","taille2mm",fin2-deb2,"taille1ref",fin1-deb1, "prct",(taille2mm/taille1ref)*100])
                    if (res==[0,0]):
                        if fin2>=deb1 :
                            deb=deb2-deb1
                            fin=fin2-fin1
                            finalresult.append([id1,"chevauchementgauche","taille2mm",fin2-deb2,"taille1ref",fin1-deb1,"couverture",((fin2-deb1)/taille1ref)*100])
                    if (res==[1,1]):
                        if deb2<=fin1 :
                            deb=deb2-deb1
                            fin=fin2-fin1
                            finalresult.append([id1,"chevauchementdroit","taille2mm",fin2-deb2,"taille1ref",fin1-deb1,"couverture",((fin1-deb2)/taille1ref)*100])


    print(len(finalresult))
file = open('resultanalyse.csv', 'w+', newline ='\n') 
with file:     
    write = csv.writer(file) 
    write.writerows(finalresult) 
print(len(finalresult))
print("totalref",totalref)
print("totalnew",totalnew)