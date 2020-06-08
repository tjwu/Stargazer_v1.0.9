#/usr/bin/env python

# Importing libraries:
import pandas as pd
import numpy as np
import os
import re
import argparse

pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Parsing arguments:
parser = argparse.ArgumentParser()
parser.add_argument('--dirname','-dir', action='store', dest='dirname', required=True, help='Sample directory name')
parser.add_argument('--ploidyfile','-p', action='store', dest='ploidy', required=True, help='Sample ploidy fle')
parser.add_argument('--outfile','-o', action='store', dest='outfile', default="out", help='Sample report output name')
parser.add_argument('--aouvar', action='store', dest='aouvar', default="~/ReformatSGResults/aou_variants.txt", help='Path to aou_variants.txt')
parser.add_argument('--starshared', action='store', dest='starshared', default="~/ReformatSGResults/starall_shared.txt", help='Path to starall_shared.txt')
parser.add_argument('--starlook', action='store', dest='starlook', default="~/ReformatSGResults/star_lookup.txt", help='Path to star_lookup.txt')
parser.add_argument('--starfunc', action='store', dest='starfunc', default="~/ReformatSGResults/star_function.txt", help='Path to star_function.txt')

args = parser.parse_args()
sampdir=args.dirname
ploidyfile=args.ploidy
outfile=args.outfile
aouvar=args.aouvar
starshared=args.starshared
starlook=args.starlook
starfunc=args.starfunc

# Open file for writing print output
outF=open("ReformatSGResults.out","w")

# Get filenames from directory:
filenames=[]
for file in os.listdir(sampdir):
    if file.endswith(".txt"):
        filenames.append(os.path.join(sampdir, file))

filenames=sorted(set(filenames))
print("These are the files that will be read:", filenames, file = outF)

# Read ploidy        
plfile=open(ploidyfile,"r")
ploidyline=plfile.readline()
print(ploidyline,file=outF)
sex=ploidyline.split("\t")
sex=sex.pop().strip()
print(sex,file=outF)
noX=sex.count("X")
print(noX,file=outF)
plfile.close()

# Read AoU variants of interest file:
aou_variants=pd.read_csv(aouvar,sep="\t")

# Read shared star alleles table:
starshared=pd.read_csv(starshared,sep='\t')

# Read star lookup table (for dpyd & g6pd):
starlookup=pd.read_csv(starlook,sep="\t")

# Read star alleles function table:
starfunc=pd.read_csv(starfunc,sep='\t')

# Function for renaming star alleles for dpyd & g6pd:
def renamestar(star,gene):
    return(starlookup.loc[(starlookup['gene']==gene) & (starlookup['name']==star),"other_names"].to_string(index=False).strip(" "))

# Rename star alleles in the data :
def redodata(gdata,gene):
    for i,row in gdata.iterrows():
        gdata.at[i,'hap1_main']=renamestar(row['hap1_main'],gene)
        gdata.at[i,'hap2_main']=renamestar(row['hap2_main'],gene)
        tmp=gdata.at[i,'dip_cand'].split(",")
        tmp1=[]
        for star in tmp:
            tmp1.append(renamestar(star,gene))
        gdata.at[i,'dip_cand']=";".join(tmp1)
    return(gdata)

# Function that splits the diplotype candidates by comma:
def diplocand(dipcand,gene):
    if gene=="g6pd":
         return(sorted(set(dipcand.split(';'))))
    else:
         return(sorted(set(dipcand.split(","))))

# Function that returns a diplotype with two haplotypes:
def return_diplotype(hap1,hap2):
    return("/".join([hap1,hap2]))

# Function that removes the star alleles that AoU not of our interest:
def onlyaou(cands,gene):
    #print(cands,gene)
    #print(aou_variants.loc[aou_variants['gene']==gene.upper(),"star allele"].squeeze())
    aou=aou_variants.loc[aou_variants['gene']==gene.upper(),"star allele"].squeeze().split(";")
    print(aou,file=outF)
    return(list(sorted(set(aou).intersection(cands))))

# Function that converts the star alleles that are not of interest to *1 | G6PD B:
def convertaou(cand,gene):
    print(cand,file=outF)
    aou=aou_variants.loc[aou_variants['gene']==gene.upper(),"star allele"].squeeze().split(";")
    x=cand in aou
    if not x and gene!="g6pd":
        return("*1")
    elif not x and gene=="g6pd":
        return("B")
    else:
        return(cand)

# Function to check and remove star alleles with shared variants:
flatten = lambda l: [item for sublist in l for item in sublist]

def check_share(cands,gene):
    ret=[]
    n=0
    #print(gene,len(cands),cands)
    while n<len(cands):
        cand=cands[n]
        if gene!="g6pd":
            shares=starshared.loc[(starshared["gene"]==gene) & (starshared["name"]==cand),"shared star alleles"]
        else:
            shares=starshared.loc[(starshared["gene"]==gene) & (starshared["other_names"]==cand),"shared star alleles"]
        if len(shares)>0:
             shares=shares.to_string(index=False).split(',')
        #print(shares)
        if len(shares)==0:
            ret.append(cand)
            cands.remove(cand)
            continue
        elif len(shares)==1 and shares[0].strip() not in cands:
            ret.append(cand)
            cands.remove(cand)
            continue
        else:
            print(cand,file=outF)
            #shares=shares.to_string(index=False).split(',')
            print("Shared:",shares,file=outF)
            checkbetwn=flatten([[cand],shares])
            print("Checkbetwn:",checkbetwn,file=outF)
            numvar=[]
            for c in checkbetwn:
                #print(c)
                c=c.strip(" ")
                if gene!="g6pd":
                     numvar.append(int(starshared.loc[(starshared["gene"]==gene) & (starshared["name"]==c),"no of variant"].to_string(index=False)))
                else:
                     numvar.append(int(starshared.loc[(starshared["gene"]==gene) & (starshared["other_names"]==c),"no of variant"].to_string(index=False)))
#
            print("No. of var:",numvar,file=outF)
            res=checkbetwn[numvar.index(max(numvar))]
            #print(res)
            if res in checkbetwn:
                checkbetwn.remove(res)
            elif res.strip(" ") in checkbetwn:
                res=res.strip(" ")
                checkbetwn.remove(res)
            for c in checkbetwn:
               #print(gene,cands,c)
               c=c.strip(" ")
               if c in cands:
                  cands.remove(c)
           # print(gene,cands)
            if len(cands)==1:
               #print(len(cands),res)
               ret.append(res)
               n=n+1
            elif len(cands)==2:
             #  print(len(cands),res)
               cand=cands[n+1]
               shares=starshared.loc[(starshared["gene"]==gene) & (starshared["name"]==cand),"shared star alleles"]
               if len(shares)==0:
              #      print(len(cands),res)
                    ret.append(res)
                    cands.remove(res)
               else:
                    n=n+1
            elif len(cands)==0:
               ret.append(res.strip(" "))
            else:
               n=n+1
    return(list(set(ret)))

# Confirm diplotypes are same or return the accurate diplotype:
def check_diplotypes(dip1,dip2):
    if dip1==dip2:
        return dip1
    else:
        return dip2

# Reorder scores:
def reorder_scores(scores):
#    print(scores)
    scoreorder=[]
    orderscore={'0':0,
                '0.1':1,
                '0.5':2,
                '1.5':3,
                '1':4,
                "unknown":5}
    for score in scores:
        scoreorder.append(orderscore[score])
    return(np.argsort(scoreorder).tolist())

# Check function of star-allele and return the most likely
def check_function(x,gene):
    score=[]
    for cand in x:
        if gene!="g6pd":
            score.append(starfunc.loc[(starfunc["gene"]==gene) & (starfunc["name"]==cand),"score"].to_string(index=False).strip())
        else:
            score.append(starfunc.loc[(starfunc["gene"]==gene) & (starfunc["other_names"]==cand),"score"])
    ord=reorder_scores(score)
    hapsret=[x[ord[0]],x[ord[1]]]
    otherret=[]
    for i in range(2,len(ord)):
         otherret.append(x[ord[i]])
    return(hapsret,otherret)

# Return the correct result 
def get_correct_result(row,gene):
    hap1=convertaou(row['hap1_main'],gene)
    hap2=convertaou(row['hap2_main'],gene)
    print("----",gene,hap1,hap2,file=outF)
    a=return_diplotype(hap1=hap1, hap2=hap2)
    print("dip1:",a,file=outF)
    dipcand=onlyaou(diplocand(row['dip_cand'],gene=gene),gene=gene)
    print(dipcand,file=outF)
    x=check_share(dipcand,gene=gene)
    print("x:",x,file=outF)
    print("dip2 candidates:",x,file=outF)
    if len(x)>2:
        x1=check_function(x,gene)
        x=x1[0]
        b=return_diplotype(hap1=x[0],hap2=x[1])
        result=check_diplotypes(a,b)
        result=result+" (More star-alleles:"+str(x1[1])+")"
    elif len(x)>1:
        b=return_diplotype(hap1=x[0],hap2=x[1])
        result=check_diplotypes(a,b)
    elif len(x)==1 and hap1=="*1" or hap1=="B" or hap2=="*1" or hap2=="B" and hap2!=x[0]:
        #print(x,hap1,hap2)
        if hap1=="*1" or hap1=="B":
            b=return_diplotype(hap1=hap1,hap2=x[0])
        else:
            b=return_diplotype(hap1=hap2,hap2=x[0])
        result=check_diplotypes(a,b)
    elif len(x)==1 and hap1=="*1" or hap1=="B" or hap2=="*1" or hap2=="B":
        if hap1=="*1" or hap1=="B":
            b=return_diplotype(hap1=hap1,hap2=x[0])
        else:
            b=return_diplotype(hap1=hap2,hap2=x[0])
        result=check_diplotypes(a,b)
    else:
        result=a
    print("Result:",result,file=outF)
    return(result)

def g6pddiphap(diplotype,noX,sex):
    #print(noX,sex)
    haps=diplotype.split("/")
    if noX==2 and sex=="XX":
    #   print(diplotype)
       return(diplotype)
    elif noX==1 and sex=="XY":
       if haps[0]==haps[1]:
          return(haps[0])
       elif haps[0]=="B":
          return(haps[1])
       elif haps[1]=="B":
          return(haps[0])
       else:
          print("check G6PD diplotype!")
          return("N/A")
    else:
       print("check ploidy!")
       print("Ploidy file contents:\n\t",ploidyline)
       return(diplotype+" (Sex:"+sex+")")

# Compile result report:
result={}
for file in filenames:
    data=pd.read_csv(file,sep='\t')
    print(data.head(),file=outF)
    gene=file.replace(sampdir,"")
    gene=re.sub("\/?","",gene)
    gene=re.sub("(_.*)?\.stargazer-genotype\.txt","",gene).lower()
    print(gene,file=outF)
    if gene=="g6pd" or gene=="dpyd":
        redodata(data,gene)
    for i,row in data.iterrows():
        if gene!="g6pd":
            dipcand=onlyaou(row['dip_cand'].split(','),gene)
            print(dipcand,file=outF)
        else:
            dipcand=onlyaou(row['dip_cand'].split(';'),gene)
            print(dipcand,file=outF)
        print("len dipcand, gene:", len(dipcand),gene,file=outF)
        if len(dipcand)==1 and ((row['hap1_main'] == dipcand[0]) or (row['hap2_main']==dipcand[0])):
            hap1=convertaou(row['hap1_main'],gene)
            hap2=convertaou(row['hap2_main'],gene)
            result[gene]=return_diplotype(hap1,hap2)
        elif len(dipcand)==1 and ((row['hap1_main'] != dipcand[0]) and (row['hap2_main']!=dipcand[0])):
            result[gene]=get_correct_result(row,gene)
        elif len(dipcand)==0 and gene!="g6pd":
            result[gene]="*1/*1"
        elif len(dipcand)==0 and gene=="g6pd":
            result[gene]="B/B"
        else:
            result[gene]=get_correct_result(row,gene)
        if gene=="g6pd":
            result[gene]=g6pddiphap(result[gene],noX,sex)
            
        
result_df=pd.DataFrame.from_dict(result,orient='index',columns=["diplotype"])
print(result_df,file=outF)


# Write result report:
result_df.to_csv(outfile, sep='\t')

outF.close()
os.remove("ReformatSGResults.out")
