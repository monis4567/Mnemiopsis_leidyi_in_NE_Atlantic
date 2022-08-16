#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################################################################################################
# python code that uses biopython to import accession numbers and species names and museum 
# voucher specimen numbers and
# family names and NCBI-taxon-Id number by searching for a higher taxon level and a selection of 
# genetic markers
# This code is specifically developed with the aim of mining the NCBI GenBank for the gene fragments: 
# mtDNA-Co1, nDNA-zic1
# , nDNA-bmp4, nDNA-H3 from the order of bony fishes: 'Myctophiformes'
# The second aim of this code is to make it possible to see which publications and authors 
# originally where associated
# with the NCBI accession number record deposited for each gene marker obtained. 
# The idea behind this is that it will allow
# for inferring whether the species identity assigned to the fish can be considered more or less 
# likely to be a trusted
# species ID. 
# If the species ID associated with the accession number stems from a study performed by researchers with 
# very little taxonomical expertise, this code should allow for identifying such untrusted 
# studies and thereby make it
# possible to avoid accession numbers where species ID potentially can have been incorrectly assigned

# A third aim with this code is to make it possible to build a local custom BLAST database that 
# only includes selected
# genes (mainly mitochondrial) and only from short list of species. The idea behind this setup is 
# that it will make it possible to
# prepare a local and limited BLAST database for the species of marine fishes in a specific geographical region 
# - e.g. the NE Atlantic Ocean.
# Restricting such a local BLAST database to only include already known occurrences, will reduce 
# the risk of getting BLAST matches
# with irrelevant species - e.g. if you have sequences that stems from samples collected in the 
# NE Atlantic, there is no need
# to match it up against an Indo-Pacific endemic or some other irrelevant tropical species.
# This will off course make it impossible to identify seqeunce reads that originates from a species 
# not included in the list
# but it should make BLAST searches faster, and return more reliable hits. 
# An additional BLAST search can then be made subsequently 
# on a global BLAST database for the more odd sequence reads, by choosing the remote blast option.

# The code can be run directly in a bash terminal by calling the code with python3 and specifying at least 
# 2 input variables.
# The first input variable must be taxon -  e.g. 'Myctophiformes' the second variable can be an mtDNA gene
# like 12s or co1. The gene fragments are separated by a comma.
# The line below is an example of how this python3 code can be run:
# python3 py01_retrieve_acc_from_taxanomy_and_gene01.py Tripterygiidae co1
# The line below is a second example of how this python3 code can be run:
# python3 py01_retrieve_acc_from_taxanomy_and_gene01.py Tripterygiidae co1,12s

# Note that this code has not been tested exhaustively. The code works for fetching COX1, 
# and 12s and 16s and the control region of the mitochondrial genome for various families of bony fish
# but the code is not yet tested on other taxon groups or other mtDNA gene fragments.
#
########################################################################################################################
# Before you try out this python code it might be worth installing ipython3, and TKinter and biopython for python3
# Use the commands below in a unix terminal and see if this can make ipython3 work

########################################################################################################################

# $ sudo apt install ipython3
# $ sudo apt-get install python-support
# $ sudo apt-get install python3-tk
# $ pip3 install biopython

# you will also need to have pandas and numpy available for python3
# install pandas for python3
# https://stackoverflow.com/questions/38768996/how-to-install-pandas-for-python-3
# $ sudo apt-get install python3-pip
# $ sudo -H pip3 install pandas

########################################################################################################################
# if in the remote directory, then start out by loading the required module
#module load Python/3.7.3-gimkl-2018b


#_______________________________________________________________________________________________________________________
# This python code is setup to be run from a bash terminal.
# you can start this code by running the code with two input variables.
# The first variable is the taxonomic group for which you are about to obtain sequences
# The second varaible is the genetic region you want to obtain gene fragments for. 
# Note the alternative spellings of the gene fragment is separated by commas

# Here is one example of how this code can be executed:
# python3 retrieve_acc_from_taxanomy_and_gene10.py Tripterygiidae 16s
# python3 retrieve_acc_from_taxanomy_and_gene09.py Enneapterygius 12s

# Here is a second example of how this code can be executed:
# python3 retrieve_acc_from_taxanomy_and_gene09.py Tripterygiidae co1
# python3 retrieve_acc_from_taxanomy_and_gene09.py Tripterygiidae cr
#_______________________________________________________________________________________________________________________

#______________________________________________________________________________________________________________________________
#  section 01 start : import modules
#______________________________________________________________________________________________________________________________


## import the biopython module
from Bio import Entrez
#also import the part that can write out a fasta files
from Bio import SeqIO
#import a module that can set a pause on the queries send to NCBI
import time
#import a module that will enable you to write out a csv file
import csv
#from urllib.request
from urllib.error import HTTPError
##report a valid email address
# a valid email address is required by biopython
Entrez.email = "swknudsen@snm.ku.dk"
#Entrez.email = "sknu003@aucklanduni.ac.nz"
##import the module that will allow you to replace regex
import re
#import a module that will enable you to delete files
import os
#https://stackoverflow.com/questions/14155669/call-python-script-from-bash-with-argument
import sys
# install pandas for python3
# https://stackoverflow.com/questions/38768996/how-to-install-pandas-for-python-3
#sudo apt-get install python3-pip
#sudo -H pip3 install pandas
#import the modules
import numpy as np
import pandas as pd

#______________________________________________________________________________________________________________________________
#  section 01 end : import modules
#______________________________________________________________________________________________________________________________
#______________________________________________________________________________________________________________________________
#  section 02 start : define input variables to use when running code from bash
#______________________________________________________________________________________________________________________________

# #______________________________________________________________________________________________________________________________
#print sys.argv[0] # prints python_script.py

var01 = sys.argv[1] # prints var1
var02 = sys.argv[2] # prints var2
#var03 = sys.argv[3] # prints var2
#var04 = sys.argv[4] # prints var2
#var05 = sys.argv[5] # prints var2
#var06 = sys.argv[6] # prints var2
#var07 = sys.argv[7] # prints var2

#var01 = 'Mnemiopsis'
#var01 = 'Sparidae'
#var02 = 'its1'

# #https://www.geeksforgeeks.org/python-program-convert-string-list/
# Define function to convert string to list 
# Notice the inout needs to be separated by comma
def convert_to_list(string): 
    li = list(string.split(",")) 
    return li 
  
# # use function to convert string to list on variable
ls_organisms = convert_to_list(var01)
ls_prod02 = convert_to_list(var02)
#ls_prod02 = convert_to_list(var03)
#ls_prod03 = convert_to_list(var04)
#ls_prod04 = convert_to_list(var05)
#ls_prod05 = convert_to_list(var06)
#ls_prod06 = convert_to_list(var07)

# the scetion above enables you to run this code like this:
# python3 py01_retrieve_acc_from_taxanomy_and_gene01.py Tripterygiidae 12s,12S,s-rRNA 

# or like this
# python3 py01_retrieve_acc_from_taxanomy_and_gene01.py Tripterygiidae co1,coI,COI,CO1,'cytochrome c oxidase subunit I'
# python3 py01_retrieve_acc_from_taxanomy_and_gene01.py Tripterygiidae co1,coI,COI,CO1,'cytochrome c oxidase subunit I', 'cytochrome c oxidase subunit 1',COX1,COXI,'cytochrome oxidase subunit I', 'cytochrome oxidase subunit 1'
#______________________________________________________________________________________________________________________________
#  section 02 end : define input variables to use when running code from bash
#______________________________________________________________________________________________________________________________


for genefragm in ls_prod02:
    ls_prod01 = convert_to_list(genefragm)

    #______________________________________________________________________________________________________________________________
    #  section 04 start : define a variable you will use to search for species - note this is a taxon level higher than species
    #______________________________________________________________________________________________________________________________
    
    # define a variable you will use to search for species - note this is a taxon level higher than species
    # but this can be changed to generic level or even 'Class'. But be careful not using 'Vertebrata' as NCBI GenBank has 
    # an enormuous number of sequences from vertebrates
    #ls_organisms=["Tetraodontiformes"]
    #ls_organisms=["Myctophiformes"]
    #ls_organisms=["Syngnathus"]
    #ls_organisms=["Salmonidae"]
    #ls_organisms=["Labrus bergylta", "Limanda limanda"]
    #ls_organisms=["Labrus bergylta", "Limanda limanda", "Syngnathus", "Epinephelus awoara"]
    #ls_organisms=["Lampadena urophaos"]
    #ls_organisms=["Clariidae"]
    
    accno="AP006025" #Enneapterygius complete mt genome, except for dloop
    #accno="LC340220" #is only 12s
    
    accno="KR381771"
    accno="MT584161"
    accno="MT584159"
    accno="HM147271"


    #define ls_genes that the code is supposed to retreive
    #ls_genes = ["bmp4", "zic1", "H3"]
    #ls_genes = ["H3"]
    #ls_genes = ["12s", "12S", "s-rRNA"]
    #it appears that 'product' might be a better match with a query
    #ls_products = ["12S"]
    
    # some 12s regions are stored in full mitochondria, and are deposited as '12s' or '12S' or
    # as 's-rRNA'. The check if the product is in the GBQaulifier record distinguishes between
    # '12s' and '12S'
    # make lists of possible product names
    # lists with mitochondrial gene fragments
    ls_products_16srRNAmtDNA = ["16s", "16S", "l-rRNA","16S ribosomal RNA"]
    ls_products_12srRNAmtDNA = ["12s", "12S", "s-rRNA","12S ribosomal RNA","srna","srRNA"]
    ls_products_COX1mtDNA =["co1","coI","Cox1","cox1","COI","CO1","cytochrome c oxidase subunit I","cytochrome oxidase I","cytochrome c oxidase subunit 1","COX1","COXI","cytochrome oxidase subunit I", "cytochrome oxidase subunit 1"]
    ls_products_dloopmtDNA=["cr","CR","control region","d loop", "D-loop", "d-loop"]
    # lists with nuclear gene fragments
    ls_products_ets2nDNA=["ets2","ETS2","ETS proto-oncogene 2", "Ets2"]
    ls_products_bmp4nDNA=["bone morphogenetic protein 4","bmp4"]
    ls_products_zic1nDNA=["zinc finger 1", "zic1"]
    ls_products_h3nDNA=["histone H3", "h3","H3"]
    ls_products_its1_2nDNA=["internal transcribed spacer 1", "internal transcribed spacer 2","Internal transcribed spacer 1", "Internal transcribed spacer 2", "ITS1","ITS2","its1","its2"]
    
    #ls_prod01=["h3"]
    # check what fragment was used as input variable, and get the appropriate list of possible fragment names to search for
    for inputprod in ls_prod01:
        if inputprod in ls_products_12srRNAmtDNA:
            ls_products = ls_products_12srRNAmtDNA
            prod02 = "12srRNAmtDNA"
        elif inputprod in ls_products_16srRNAmtDNA:
            ls_products = ls_products_16srRNAmtDNA
            prod02 = "16srRNAmtDNA"
        elif inputprod in ls_products_COX1mtDNA:
            ls_products = ls_products_COX1mtDNA
            prod02 = "COX1mtDNA"
        elif inputprod in ls_products_dloopmtDNA:
            ls_products = ls_products_dloopmtDNA
            prod02 = "dloopmtDNA"
        elif inputprod in ls_products_ets2nDNA:
            ls_products = ls_products_ets2nDNA
            prod02 = "ets2nDNA"
        elif inputprod in ls_products_bmp4nDNA:
            ls_products = ls_products_bmp4nDNA
            prod02 = "bmp4nDNA"
        elif inputprod in ls_products_zic1nDNA:
            ls_products = ls_products_zic1nDNA
            prod02 = "zic1nDNA"
        elif inputprod in ls_products_h3nDNA:
            ls_products = ls_products_h3nDNA
            prod02 = "h3nDNA"
        elif inputprod in ls_products_its1_2nDNA:
            ls_products = ls_products_its1_2nDNA
            prod02 = "its1its2"
        else:
            print("gene fragment name provided as input variable is not recognized by code")
    
    #make an empty list to fill with organism names
    ls_organisms02=[]
    #______________________________________________________________________________________________________________________________
    #  section 04 end : define a variable you will use to search for species - note this is a taxon level higher than species
    #______________________________________________________________________________________________________________________________
    
    #______________________________________________________________________________________________________________________________
    #  section 03 start : define output files and delete previous version
    #______________________________________________________________________________________________________________________________
    #put filenames from previous outputs in to a variable
    fileP_05_01_accn_spnm = 'output01_accn_spnm_'+var01+'_'+prod02+'.txt'
    fileP_05_02_accn_ncbitaxhierarchy = 'output01_accn_ncbitaxhierarchy_'+var01+'_'+prod02+'.txt'
    fileP_05_03_accn_fam = 'output01_accn_fam_'+var01+'_'+prod02+'.txt'
    fileP_05_04_accn_taxno = 'output01_accn_taxno_'+var01+'_'+prod02+'.txt'
    fileP_05_05_accn_vouch = 'output01_accn_vouch_'+var01+'_'+prod02+'.txt'
    fileP_05_06_accn_ntseq = 'output01_accn_ntseq_'+var01+'_'+prod02+'.txt'
    fileP_05_07_accn_publication = 'output01_accn_publication_'+var01+'_'+prod02+'.csv'
    
    #check in an 'if-else' loop if the previous files should be deleted or not
    # this is to make it possible to avoid deleting files, or start deleting files
    # with just commenting out a single line instead of commenting out the entire block 
    # below if the iteration is started again with the intention
    # of appending to files that already have acceession numbers and species and
    # information
    
    #comment out this line , or  comment out the next line
    #crit_for_del="append_to_previous_list"
    crit_for_del="do_not_append_to_previous_list"
    # the iteration that deletes the files
    if crit_for_del == "do_not_append_to_previous_list":
        ## check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_01_accn_spnm):
            os.remove(fileP_05_01_accn_spnm)
        else:
            print("Cannot delete",fileP_05_01_accn_spnm)
        ## check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_02_accn_ncbitaxhierarchy):
            os.remove(fileP_05_02_accn_ncbitaxhierarchy)
        else:
            print("Cannot delete",fileP_05_02_accn_ncbitaxhierarchy)
        ## check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_03_accn_fam):
            os.remove(fileP_05_03_accn_fam)
        else:
            print("Cannot delete",fileP_05_03_accn_fam)
        ## check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_04_accn_taxno):
            os.remove(fileP_05_04_accn_taxno)
        else:
            print("Cannot delete",fileP_05_04_accn_taxno)
        ## check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_05_accn_vouch):
            os.remove(fileP_05_05_accn_vouch)
        else:
            print("Cannot delete",fileP_05_05_accn_vouch)
        #    # check if file exists and delete it in case it exists
        if os.path.exists(fileP_05_06_accn_ntseq):
            os.remove(fileP_05_06_accn_ntseq)
        else:
            print("Cannot delete",fileP_05_06_accn_ntseq)
    else:
        pass
    #______________________________________________________________________________________________________________________________
    #  section 03 end : define output files and delete previous version
    #______________________________________________________________________________________________________________________________
    
    
    #______________________________________________________________________________________________________________________________
    #  section 05 start : Use section to iterate over a list of species names
    #______________________________________________________________________________________________________________________________
    #flnm_lsspcs="lst_spc_marine_fish_Kenya_long_lst03.txt"
    # flnm_lsspcs = ls_organisms
    # #https://runestone.academy/runestone/books/published/thinkcspy/Files/Iteratingoverlinesinafile.html
    # #iterate over lines in the file and place the elements in the empty list
    # #open read the file
    # qbfile = open(flnm_lsspcs, "r")
    # for aline in qbfile:
    #     # Replace the end of line in each entry , see this website: https://stackoverflow.com/questions/3094659/editing-elements-in-a-list-in-python
    #     aline2=aline.replace('\n', '')
    #     #check if the accession number already is in the list, in case it is, then report that it is already in the list 
    #     if aline2 in ls_organisms02:
    #         print(aline2,"is already in list of species and is not appended to list 'ls_organisms02'")
    #     else:
    #         #if the accession number is not in the list then print the number and append it to the list
    #         print(aline2)
    #     # append to the empty list while iterating : https://stackoverflow.com/questions/41452819/list-append-in-for-loop
    #         ls_organisms02.append(aline2)
    #         #print(aline2)
    # #close read the file
    # qbfile.close()
    # #count the number of elements in the list
    # vr_elem_in_ls_organisms02=len(ls_organisms02)
    # print("the list 'ls_organisms02' has",vr_elem_in_ls_organisms02," species")
    
    #______________________________________________________________________________________________________________________________
    #  section 05 end : Use section to iterate over a list of species names
    #______________________________________________________________________________________________________________________________
    
    #______________________________________________________________________________________________________________________________
    #  section 06 start : Use section to continue iterate over a list of species names that you have started fetching but not finished
    #______________________________________________________________________________________________________________________________
    
    # #make an empty list to fill with organism names
    # ls_organisms03=[]
    # #import the regex module to be able to check if the line starts with
    # import re
    # #f re.match(r'^Thu', somestring):
    #     # do stuff
    
    # # If a previous search on GenBank stopped at some letter
    # # then use the if part below to modify the 'ls_organisms02' list
    # # By setting the capital letter to search for to something else than 'A-Z'
    # # e.g. setting it to 'T-Z' will result in a new list with only strings 
    # # starting with T-Z.
    # # With this the list of species names can be adjusted to exclude the species
    # # for which a previous tryout with this code already have obtained accession
    # # information
    
    # # Iterate each element in list 
    # # and add them in variale total 
    # for ele in ls_organisms02: 
    #     if re.match(r'^[T-Z]+', ele):
    #         # append to the empty list while iterating : https://stackoverflow.com/questions/41452819/list-append-in-for-loop
    #         ls_organisms03.append(ele)
    #         #print(ele) 
    #         #list1.remove(ele) 
    #     else:
    #         pass
    
    # #replace the previous list with the newly prepared list
    # #ls_organisms=ls_organisms02
    # #ls_organisms=ls_organisms03
    #______________________________________________________________________________________________________________________________
    #  section 06 end : Use section to continue iterate over a list of species names that you have started fetching but not finished
    #______________________________________________________________________________________________________________________________
    
    #______________________________________________________________________________________________________________________________
    #  section 07 start : Iterate over elements in list of organisms
    #______________________________________________________________________________________________________________________________
    
    #make empty dictionaries
    dct_specs = {}
    dct_species = {}
    dct_spls_genes = {}
    dct_ls_genesp = {}
    dct_pubdoi ={}
    dct_pubauthor ={}
    dct_pubjournal={}
    dct_pubtitle={}
    dct_tax ={}
    dct_fam={}
    dct_product={}
    dct_speciesgb = {}
    dct_speciesgb02 = {}
    dct_vouchergb={}
    dct_taxongb={}
    dct_taxno={}
    dct_gene = {}
    dct_location = {}
    dct_latlonpos = {}
    
    # 
    #___________________
    #for every element in the list with ls_organisms do ... -  i.e. iterate over all elements in organsims
    for org in ls_organisms:
        #for gene in ls_genes:
        for product in ls_products: #-  i.e. iterate over all elements in ls_products
            #query= org+"[organism] AND "+gene+"[gene]"
            #query= org+"[organism] AND "+product+"[product]"
            #query= org+" AND "+product+" AND (\"1\"[SLEN] : \"20000\"[SLEN])"
            query= org+" AND "+product+" AND 1:20000[Sequence Length]"
            #query= org+"[organism] AND "+gene+"[gene] AND "+accno+"[accn]"
            #query= product+"[product] AND "+accno+"[accn]"
            #query= product+" AND "+accno
            res = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
            rec = Entrez.read(res)
            res = Entrez.efetch(db="nucleotide", id=rec["IdList"],  retmode = "xml")
            #setting rettype to "gb" will return the full genbank record, and this includes all extra information, including details on voucher specimens
            # see this question about setting try / except for avoiding HTTPError : https://www.biostars.org/p/126143/
            time.sleep(8) # to make sure not many requests go per second to ncbi           
            #try in try test 01
            try:
                resgb = Entrez.efetch(db="nucleotide", id=rec["IdList"],  rettype="gb", retmode="xml")
            #except in try test 01
            except HTTPError:
                time.sleep(80)
                resgb = Entrez.efetch(db="nucleotide", id=rec["IdList"],  rettype="gb", retmode="xml")
            
    #http://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    
            # test if the record is in the 'Entrez.read(res)', this will be empty
            # if the search did not return any hits
            #try in try test 02
            try:
                #iterate over all the records, provided the 'try' works
                # for test 01 -start
                for record in Entrez.read(res):
                    #print(len(record))
                    print("reading record for ",accn)
                    # for the record in each search try to get the organism name and accession number
                    speciesName = record["GBSeq_organism"]
                    accn = record["GBSeq_accession-version"]
                    #split the accession number by the point
                    accn = accn.split('.')[0]
                    #dct_gene[accn] = gene
                    dct_product[accn] = prod02
                    tax = record["GBSeq_taxonomy"]
                    family=tax.split(";")[-2].strip()
                    #assert accn not in dct_species
                    if accn in dct_species:
                        pass
                        #print(accn,"is already in dct_species")
                    else:
                        dct_species[accn] = speciesName
                        #prepare the line to be written to the file as a string
                        lninfl=f"{accn}\t{speciesName}\n"
                        #Open a file and use 'a' to append to a file
                        f = open(fileP_05_01_accn_spnm, 'a')
                        #write the line in the string
                        f.write(lninfl) 
                        #close the file again
                        f.close()
    
                        #dct_doi[accn] = doi
                        #dct_author[accn] = author
                        dct_tax[accn] = tax
                        # make the line a string before it is written
                        lninfl=f"{accn}\t{tax}\n"
                        #use 'a' to append to a file
                        f = open(fileP_05_02_accn_ncbitaxhierarchy, 'a')
                        f.write(lninfl) 
                        f.close()
    
                        dct_fam[accn] = family
                        #use 'a' to append to a file
                        lninfl=f"{accn}\t{family}\n"
                        f = open(fileP_05_03_accn_fam, 'a')
                        f.write(lninfl) 
                        f.close()
            # as some searches will return nothing, the 'try' will not work
            # return 'empty record' for these searches
            #except in try test 02
            except:
                pass
                #print("empty search record")
                # for test 01 -end
                
            try:
                # for test 02 -start
                for recordgb in Entrez.read(resgb):
                    accngb = recordgb["GBSeq_accession-version"]
                    accn = accngb.split('.')[0]
                    # check 'accn' is not already in the other dictionary 'dct_speciesgb'
                    #assert accn not in dct_speciesgb
                    if accn in dct_speciesgb:
                        pass
                        #print(accn,"is already in 'dct_speciesgb'")
                    else:
                        #append the organism name to the dictionary of species names
                        dct_speciesgb[accn] = recordgb["GBSeq_organism"]
                        #append detailes about the publication to dictionaries with accession number
                        # as a key
                        try:
                            dct_pubdoi[accn] = recordgb['GBSeq_references'][0]['GBReference_xref'][0]['GBXref_id']
                        except:
                            dct_pubdoi[accn] = "NA"
                        try:
                            dct_pubauthor[accn] = recordgb['GBSeq_references'][0]['GBReference_authors']
                        except:
                            dct_pubauthor[accn] = "NA"
                        try:
                            dct_pubjournal[accn] = recordgb['GBSeq_references'][0]['GBReference_journal']
                        except:
                            dct_pubjournal[accn] = "NA"
                        try:
                            dct_pubtitle[accn] = recordgb['GBSeq_references'][0]['GBReference_title']
                        except:
                            dct_pubtitle[accn] = "NA"
                        # dct_pubtitle[accn] = recordgb['GBSeq_references'][0]['GBReference_title']
                        try:
                            dct_location[accn] = recordgb['GBSeq_feature-table'][0]['GBFeature_quals'][4]['GBQualifier_value']
                        except:
                            dct_location[accn] = "NA"

                        try:
                            dct_latlonpos[accn] = recordgb['GBSeq_feature-table'][0]['GBFeature_quals'][5]['GBQualifier_value']
                        except:
                            dct_latlonpos[accn] = "NA"
                        dct_speciesgb02[accn] = recordgb
                        # loop over 'qual' in all the 'recordgb'-elements
                        for qual in recordgb['GBSeq_feature-table'][0]['GBFeature_quals']:
                            #get the qual value and put it in the variable 'name'
                            name = qual['GBQualifier_name']
                            #get the qual value and put it in the variable 'value'
                            value = qual['GBQualifier_value']                        
                            #check if the 'name' matches 'specimen_voucher'
                            if name == 'specimen_voucher':
                                #and if it matches put it in the dictionary 'dct_vouchergb'
                                dct_vouchergb[accn] = value
                                vouch01 = str(value)
                                #use 'a' to append to a file
                                lninfl=f"{accn}\t{vouch01}\n"
                                f = open(fileP_05_05_accn_vouch, 'a')
                                f.write(lninfl) 
                                f.close()
                            #also check if any of the 'name' matches 'db_xref'
                            elif name == 'db_xref':
                                #and if it matches and starts with 'taxon'
                                if value.startswith('taxon:'):
                                #then add this value to the dictionary 'dct_taxno'
                                    dct_taxno[accn] = value
                                    val01 = str(value)
                                    val02 = re.sub(r" ","_",val01)
                                    #use 'a' to append to a file
                                    lninfl=f"{accn}\t{val02}\n"
                                    f = open(fileP_05_04_accn_taxno, 'a')
                                    f.write(lninfl) 
                                    f.close()
                                else:
                                    pass
                                    #print(value,"does not appear to be a taxon")
                            else:
                                
                                pass
                                #print(name,"for name is not voucher or db_xref")
            except:
                pass
                #print("not in search 02")
    #______________________________________________________________________________________________________________________________
    #  section 07 end : Iterate over elements in list of organisms
    #______________________________________________________________________________________________________________________________
    
    #______________________________________________________________________________________________________________________________
    #  section 08 start : Iterate over elements in dictionary with fetched gebnbank records
    #______________________________________________________________________________________________________________________________
    #define functions to use in the elif-ladder below:
    #___________________
    # make an inner function to test in an elid ladder
    # see this question : https://stackoverflow.com/questions/2069662/how-to-exit-an-if-clause
    # define a function to check if product is in recordgb_path, if it is then append the sequence to the dictionary
    # from the correct seqpath01
    # if the product is not matched then return to the outer elif-ladder
    def ch_recbg_seqpath01(prod,recordgb_path,dcts,accnnumb,seqpath01):
                #if recordgb_path.startswith(prod):
                if prod in recordgb_path:
                    dcts[accnnumb] = seqpath01
                return
                #else:
                #    print("prod did not match")
                #    return
    
    #_____________________________________________________________________________________________________________________________
    dct_seqstrandgb={}
    #make empty dictionaries for recordgb_paths
    dct_GBft1_GBfeq0_GBQv = {}
    dct_GBft1_GBfeq2_GBQv = {}
    dct_GBft2_GBfeq0_GBQv = {}
    dct_GBft3_GBfeq0_GBQv = {}
    dct_GBft3_GBfeq1_GBQv = {}
    dct_GBft4_GBfeq0_GBQv = {}
    dct_GBft5_GBfeq0_GBQv = {}
    dct_GBft6_GBfeq0_GBQv = {}
    dct_GBft6_GBfeq1_GBQv = {}
    dct_GBft8_GBfeq1_GBQv = {}
    #dct_GBft9_GBfeq1_GBQv = {}
    dct_GBft9_GBfeq0_GBQv = {}
    dct_GBft12_GBfeq0_GBQv = {} #for CO1
    dct_GBft13_GBfeq0_GBQv = {} #for CO1
    dct_GBft15_GBfeq0_GBQv = {} #for CO1
    dct_GBft16_GBfeq0_GBQv = {} #for CO1
    dct_GBft18_GBfeq0_GBQv = {}
    dct_GBft20_GBfeq0_GBQv = {}
    dct_GBft28_GBfeq0_GBQv = {}
    dct_GBft30_GBfeq0_GBQv = {}
    dct_GBft33_GBfeq0_GBQv = {}
    dct_GBft34_GBfeq0_GBQv = {}
    dct_GBft38_GBfeq0_GBQv = {} #for CO1 
    dct_GBft40_GBfeq0_GBQv = {} #for CO1
    dct_GBft42_GBfeq0_GBQv = {} #for CO1
    dct_GBft44_GBfeq0_GBQv = {} #for CO1
    dct_GBft46_GBfeq0_GBQv = {} #for CO1
    dct_GBft1_GBfek = {}
    dct_GBft2_GBfek = {}
    dct_GBft40_GBfek = {} #for CO1
    dct_GBft42_GBfek = {} #for CO1
    dct_GBft44_GBfek = {} #for CO1
    dct_GBft50_GBfek = {} # for dloop
    dct_GBft51_GBfek = {} # for dloop
    dct_GBft52_GBfek = {} # for dloop
    dct_GBft75_GBfek = {} # for dloop
    dct_GBseq1ft1qv={}
    dct_GBseq2ft2qv={}
    dct_GBseq2ft1qv={}
    dct_GBseq3ft1qv={}
    dct_GBseq4ft1qv={}
    dct_GBseq5ft1qv={}
    dct_GBseq4ft3qv={}
    dct_GBseq6ft2qv={}
    dct_GBseq6ft3qv={}
    dct_GBseq8ft3qv={}
    dct_GBseq={}
    #_____________________________________________________________________________________________________________________________-
    #_____________________________________________________________________________________________________________________________-
    #iterate over both keys and values in dictionary
    #https://stackoverflow.com/questions/3294889/iterating-over-dictionaries-using-for-loops
    # example
    #for k,v in d.items():
    #    print(k, 'corresponds to', v)
    for accn02,recordgb02 in dct_speciesgb02.items():
        accn=accn02
        #print to screen how far the code has come
        print("fetching sequence for ",accn)
        try:
            dct_GBft1_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft1_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft1_GBfeq2_GBQv[accn] = recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][2]['GBQualifier_value']
        except:
            dct_GBft1_GBfeq2_GBQv[accn] = "NA"
        try:
            dct_GBft2_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft2_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft3_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][3]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft3_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft3_GBfeq1_GBQv[accn] = recordgb02['GBSeq_feature-table'][3]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBft3_GBfeq1_GBQv[accn] = "NA"
        try:
            dct_GBft4_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft4_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft5_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][5]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft5_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft6_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft6_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft6_GBfeq1_GBQv[accn] = recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBft6_GBfeq1_GBQv[accn] = "NA"
        try:
            dct_GBft8_GBfeq1_GBQv[accn] = recordgb02['GBSeq_feature-table'][8]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBft8_GBfeq1_GBQv[accn] = "NA"
        
        try:
            dct_GBft9_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][9]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft9_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft12_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][12]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft12_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft13_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][13]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft13_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft15_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][15]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft15_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft16_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][16]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft16_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft18_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][18]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft18_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft20_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][20]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft20_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft28_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][28]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft28_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft30_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][30]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft30_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft34_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][34]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft34_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft33_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][33]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft33_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft38_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][38]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft38_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft40_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][40]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft40_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft42_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][42]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft42_GBfeq0_GBQv[accn] = "NA"
        try:
            dct_GBft44_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][44]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft44_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft46_GBfeq0_GBQv[accn] = recordgb02['GBSeq_feature-table'][46]['GBFeature_quals'][0]['GBQualifier_value']
        except:
            dct_GBft46_GBfeq0_GBQv[accn] = "NA"

        try:
            dct_GBft1_GBfek[accn] = recordgb02['GBSeq_feature-table'][1]['GBFeature_key']
        except:
            dct_GBft1_GBfek[accn] = "NA"
        try:
            dct_GBft2_GBfek[accn] = recordgb02['GBSeq_feature-table'][2]['GBFeature_key']
        except:
            dct_GBft2_GBfek[accn] = "NA"
        
        try:
            dct_GBft40_GBfek[accn] = recordgb02['GBSeq_feature-table'][40]['GBFeature_key']
        except:
            dct_GBft40_GBfek[accn] = "NA"

        try:
            dct_GBft44_GBfek[accn] = recordgb02['GBSeq_feature-table'][44]['GBFeature_key']
        except:
            dct_GBft44_GBfek[accn] = "NA"

        try:
            dct_GBft50_GBfek[accn] = recordgb02['GBSeq_feature-table'][50]['GBFeature_key']
        except:
            dct_GBft50_GBfek[accn] = "NA"
        try:
            dct_GBft51_GBfek[accn] = recordgb02['GBSeq_feature-table'][51]['GBFeature_key']
        except:
            dct_GBft51_GBfek[accn] = "NA"
        try:
            dct_GBft52_GBfek[accn] = recordgb02['GBSeq_feature-table'][52]['GBFeature_key']
        except:
            dct_GBft52_GBfek[accn] = "NA"
        try:
            dct_GBft75_GBfek[accn] = recordgb02['GBSeq_feature-table'][75]['GBFeature_key']
        except:
            dct_GBft75_GBfek[accn] = "NA"
    #_____________________________________________________________________________________________________________________________-
    #_____________________________________________________________________________________________________________________________-
        try:
            dct_GBseq1ft1qv[accn] = recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBseq1ft1qv[accn] = "NA"            
        try:
            dct_GBseq2ft2qv[accn] = recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][2]['GBQualifier_value']
        except:
            dct_GBseq2ft2qv[accn] = "NA"
        try:
            dct_GBseq2ft1qv[accn] = recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBseq2ft1qv[accn] = "NA"
        try:
            dct_GBseq3ft1qv[accn] = recordgb02['GBSeq_feature-table'][3]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBseq3ft1qv[accn] = "NA"
        try:
            dct_GBseq4ft1qv[accn] = recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBseq4ft1qv[accn] = "NA"
        try:
            dct_GBseq5ft1qv[accn] = recordgb02['GBSeq_feature-table'][5]['GBFeature_quals'][1]['GBQualifier_value']
        except:
            dct_GBseq5ft1qv[accn] = "NA"
        try:
            dct_GBseq4ft3qv[accn] = recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][3]['GBQualifier_value']
        except:
            dct_GBseq4ft3qv[accn] = "NA"
        try:
            dct_GBseq6ft2qv[accn] = recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][2]['GBQualifier_value']
        except:
            dct_GBseq6ft2qv[accn] = "NA"
        try:
            dct_GBseq6ft3qv[accn] = recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][3]['GBQualifier_value']
        except:
            dct_GBseq6ft3qv[accn] = "NA"
        try:
            dct_GBseq8ft3qv[accn] = recordgb02['GBSeq_feature-table'][8]['GBFeature_quals'][3]['GBQualifier_value']
        except:
            dct_GBseq8ft3qv[accn] = "NA"
        try:
            dct_GBseq[accn] = recordgb02['GBSeq_sequence']
        except:
            dct_GBseq[accn] = "NA"
    #_____________________________________________________________________________________________________________________________-
    #_____________________________________________________________________________________________________________________________-
    #_____________________________________________________________________________________________________________________________-
            
        for product in ls_products:
            #print(product)
    
            #start - check if the recordgb exists, and if it does use the functions defined above - notice the functions returns to the outer loop if false
            if (product in dct_GBft1_GBfeq0_GBQv[accn]) and (dct_GBseq1ft1qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][1]['GBQualifier_value'])
                #print("in 1")
            elif (product in dct_GBft1_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"): 
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_sequence'])
                #print("in 2")
            elif (product in dct_GBft1_GBfeq0_GBQv[accn]) and (dct_GBseq2ft2qv[accn] != "NA"): 
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][2]['GBQualifier_value'])
                #print("in 2")
            elif (product in dct_GBft1_GBfeq2_GBQv[accn]) and (dct_GBseq[accn] != "NA"): 
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][1]['GBFeature_quals'][2]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_sequence'])
                #print("in 3")
            elif (product in dct_GBft2_GBfeq0_GBQv[accn]) and (dct_GBseq2ft2qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][2]['GBQualifier_value'])
                #print("in 4")
            elif (product in dct_GBft2_GBfeq0_GBQv[accn]) and (dct_GBseq2ft1qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][1]['GBQualifier_value'])
                #print("in 4")


            elif (product in dct_GBft2_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][2]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_sequence'])
                #print("in 4")


            elif (product in dct_GBft3_GBfeq0_GBQv[accn]) and (dct_GBseq3ft1qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][3]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][3]['GBFeature_quals'][1]['GBQualifier_value'])
                #print("in 5")
            elif (product in dct_GBft4_GBfeq0_GBQv[accn]) and (dct_GBseq4ft1qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][1]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft5_GBfeq0_GBQv[accn]) and (dct_GBseq5ft1qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][5]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][5]['GBFeature_quals'][1]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft4_GBfeq0_GBQv[accn]) and (dct_GBseq4ft3qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][4]['GBFeature_quals'][3]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft6_GBfeq0_GBQv[accn]) and (dct_GBseq6ft2qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][2]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft6_GBfeq1_GBQv[accn]) and (dct_GBseq6ft3qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][1]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][3]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft8_GBfeq1_GBQv[accn]) and (dct_GBseq8ft3qv[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][8]['GBFeature_quals'][1]['GBQualifier_value'],dct_seqstrandgb,accn,recordgb02['GBSeq_feature-table'][8]['GBFeature_quals'][3]['GBQualifier_value'])
                #print("in 6")
            elif (product in dct_GBft1_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][1]['GBFeature_key'],dct_seqstrandgb,accn,recordgb02['GBSeq_sequence'])
                #print("in 7")
            elif (product in dct_GBft2_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][2]['GBFeature_key'],dct_seqstrandgb,accn,recordgb02['GBSeq_sequence'])
                #print("in 8")

            elif (product in dct_GBft6_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                # For 'Clariidae' and 'cox1' some of the intervals appear to be stored in the 7th 'GBSeq_feature-table' instead of the 9th 'GBSeq_feature-table'
                # the try except part first looks for the 9th part, and if this outside the range of the list index, 
                #then instead looks for the 7th part, if this also fails then look in th 6th part
                try:
                    int_GBftt6_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][9]['GBFeature_intervals'][0]['GBInterval_from'])
                except:
                    try:
                        int_GBftt6_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][7]['GBFeature_intervals'][0]['GBInterval_from'])
                    except:
                        int_GBftt6_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][6]['GBFeature_intervals'][0]['GBInterval_from'])
                try:
                    int_GBftt6_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][9]['GBFeature_intervals'][0]['GBInterval_to'])
                except:
                    try:
                        int_GBftt6_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][7]['GBFeature_intervals'][0]['GBInterval_to'])
                    except:
                        int_GBftt6_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][6]['GBFeature_intervals'][0]['GBInterval_to'])
                try:
                    #subtract one because Python starts with zero
                    int_recGBf = (int_GBftt6_GBfti0_GBif-1)
                    int_recGBt =(int_GBftt6_GBfti0_GBit-1)
                    # get the sequence
                    str_GBseq01 = recordgb02['GBSeq_sequence']
                    #slice the string
                    str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                    ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][6]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
                #For 'Clariidae' and 'cox1' the 'try' above will most likely not work for an accession number where the 'GBSeq_feature-table' is outside the list
                # index range. In this case no sequence will be assigned, but instead NA will be added for the sequence
                # instead I have assigned 'NA' to 'str_GBseq02'
                except:
                    int_recGBf = "0"
                    int_recGBt = "0"
                    str_GBseq02 ="NA"


            elif (product in dct_GBft9_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt9_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][9]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt9_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][9]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt9_GBfti0_GBif-1)
                int_recGBt =(int_GBftt9_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][9]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)


            
            elif (product in dct_GBft12_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt12_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][12]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt12_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][12]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt12_GBfti0_GBif-1)
                int_recGBt =(int_GBftt12_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][12]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft13_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt13_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][13]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt13_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][13]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt13_GBfti0_GBif-1)
                int_recGBt =(int_GBftt13_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][13]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft15_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt15_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][15]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt15_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][15]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt15_GBfti0_GBif-1)
                int_recGBt =(int_GBftt15_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][15]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft16_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt16_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][16]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt16_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][16]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt16_GBfti0_GBif-1)
                int_recGBt =(int_GBftt16_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][16]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)


            elif (product in dct_GBft18_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt18_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][18]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt18_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][18]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt18_GBfti0_GBif-1)
                int_recGBt =(int_GBftt18_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][18]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
            

            elif (product in dct_GBft20_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt20_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][20]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt20_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][20]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt20_GBfti0_GBif-1)
                int_recGBt =(int_GBftt20_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][20]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
            
            elif (product in dct_GBft28_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt28_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][28]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt28_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][28]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt28_GBfti0_GBif-1)
                int_recGBt =(int_GBftt28_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][28]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
            
            elif (product in dct_GBft30_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt30_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][30]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt30_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][30]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt30_GBfti0_GBif-1)
                int_recGBt =(int_GBftt30_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][30]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft33_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt33_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][33]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt33_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][33]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt33_GBfti0_GBif-1)
                int_recGBt =(int_GBftt33_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][33]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
            
            elif (product in dct_GBft34_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt34_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][34]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt34_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][34]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt34_GBfti0_GBif-1)
                int_recGBt =(int_GBftt34_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][34]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)


            elif (product in dct_GBft38_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt38_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][38]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt38_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][38]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt38_GBfti0_GBif-1)
                int_recGBt =(int_GBftt38_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][38]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)


            elif (product in dct_GBft40_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt40_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][40]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt40_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][40]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt40_GBfti0_GBif-1)
                int_recGBt =(int_GBftt40_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][40]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft42_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt42_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][42]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt42_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][42]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt42_GBfti0_GBif-1)
                int_recGBt =(int_GBftt42_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][42]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)

            elif (product in dct_GBft44_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt44_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][44]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt44_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][44]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt44_GBfti0_GBif-1)
                int_recGBt =(int_GBftt44_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][44]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
                        
            elif (product in dct_GBft46_GBfeq0_GBQv[accn]) and (dct_GBseq[accn] != "NA"):
                #to get CO1 from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt46_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][46]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt46_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][46]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt46_GBfti0_GBif-1)
                int_recGBt =(int_GBftt46_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][46]['GBFeature_quals'][0]['GBQualifier_value'],dct_seqstrandgb,accn,str_GBseq02)
             
            elif (product in dct_GBft50_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                #to get dloop from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt50_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][50]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt50_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][50]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt50_GBfti0_GBif-1)
                int_recGBt =(int_GBftt50_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][50]['GBFeature_key'],dct_seqstrandgb,accn,str_GBseq02)
                #print("in 9")        
            elif (product in dct_GBft51_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                #to get dloop from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt51_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][51]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt51_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][51]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt51_GBfti0_GBif-1)
                int_recGBt =(int_GBftt51_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][51]['GBFeature_key'],dct_seqstrandgb,accn,str_GBseq02)
                #print("in 9")
            elif (product in dct_GBft52_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                #to get dloop from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt52_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][52]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt52_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][52]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt52_GBfti0_GBif-1)
                int_recGBt =(int_GBftt52_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][52]['GBFeature_key'],dct_seqstrandgb,accn,str_GBseq02)
                #print("in 9")
            elif (product in dct_GBft75_GBfek[accn]) and (dct_GBseq[accn] != "NA"):
                #to get dloop from mt genome - also get the 'from' and 'to' to get the position in the mt-genome, and turn the numbers in to integers with 'int'
                int_GBftt75_GBfti0_GBif = int(recordgb02['GBSeq_feature-table'][75]['GBFeature_intervals'][0]['GBInterval_from'])
                int_GBftt75_GBfti0_GBit = int(recordgb02['GBSeq_feature-table'][75]['GBFeature_intervals'][0]['GBInterval_to'])
                #subtract one because Python starts with zero
                int_recGBf = (int_GBftt75_GBfti0_GBif-1)
                int_recGBt =(int_GBftt75_GBfti0_GBit-1)
                # get the sequence
                str_GBseq01 = recordgb02['GBSeq_sequence']
                #slice the string
                str_GBseq02 = str_GBseq01[int_recGBf:int_recGBt]
                ch_recbg_seqpath01(product,recordgb02['GBSeq_feature-table'][75]['GBFeature_key'],dct_seqstrandgb,accn,str_GBseq02)
                #print("in 9")
            else:
                pass
                #print("no match")
    
    #________
    #define a function to get keys in a dictionary: https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
    def getListkeys(dict): 
        return dict.keys() 
    #________
    #use the function on your dictionaries
    dctk_accno_speciesgb02=getListkeys(dct_speciesgb02)
    dctk_accno_seqstrandgb=getListkeys(dct_seqstrandgb)
    #check what kind of object you have with 'type'
    type(dctk_accno_seqstrandgb)
    # it is a 'dict_keys' object
    # turn it in to a list
    ls_spg02=list(dctk_accno_speciesgb02)
    ls_seqst=list(dctk_accno_seqstrandgb)
    #find the difference between the two lists
    #https://stackoverflow.com/questions/3462143/get-difference-between-two-lists
    ls_dfspgb_seqst=list(set(ls_spg02) - set(ls_seqst))
    #print a comment to screen
    print("ensure the accession numbers below are not needed for the search made")
    print("if only '[]' appears then all records read have also been fetched")
    print(ls_dfspgb_seqst)
    
    
    #______________________________________________________________________________________________________________________________
    #  section 08 end : Iterate over elements in dictionary with fetched gebnbank records
    #______________________________________________________________________________________________________________________________
    
    #______________________________________________________________________________________________________________________________
    #  section 09 start : Build dictionaries from the files and dictionaries to prepare a csv table with all information and sequences
    #______________________________________________________________________________________________________________________________
    
    #define a function to get keys in a dictionary: https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
    def getListkeys(dict): 
        return dict.keys() 
    #use the function on your dictionaries
    ls_accno_seqstrandgb=getListkeys(dct_seqstrandgb)
    ls_accno_product=getListkeys(dct_product)
    #Get the difference between the two lists:
    #this will help identify which accession numbers the code failed to get sequences to match to
    #https://stackoverflow.com/questions/3462143/get-difference-between-two-lists
    list(set(ls_accno_product) - set(ls_accno_seqstrandgb))
    
    # get value for specific key
    #str_seq_AF354976=dct_seqstrandgb['AF354976']
    # get the lenght of this text string
    #len(str_seq_AF354976)
    
    #define a function to get values from a dictionary
    def getListvalues(dict): 
        return dict.values()
    #get the values from a dictionary and place in a list
    ls_val_dct_seqstrandgb=getListvalues(dct_seqstrandgb)
    
    #make an empty dictionary to fill with sequence lengths
    dct_lenseq01={}
    #iterate over both keys and values in dictionary
    #https://stackoverflow.com/questions/3294889/iterating-over-dictionaries-using-for-loops
    #for k,v in d.items():
    #    print(k, 'corresponds to', v)
    #using the example on the 'dct_seqstrandgb' dictionary prepared earlier
    # iterating over both keys and values in the items in the dictionary,
    # where items cover both keys and values
    for accno,ntseq in list(dct_seqstrandgb.items()):
        #get the lenght of the value in the dictionary -  i.e. the length of the nt seq
        ls_lenseq01=len(ntseq)
        #place this length in the empty dictionary prepared
        dct_lenseq01[accno]=ls_lenseq01
        #print the length obtained for each value
        #print(ls_lenseq01)
    #note that some sequences are very short, but this appears to be because they are partial
    #dct_seqstrandgb['HM773107']
    # looking up 'HM773107' on NCBI Nucleotide database showed that this sequence is mainly mtDNA
    # control region, and that only 8 nts are from the '12s' region. So the loop above appears to work.
    
    
    # install pandas for python3
    # https://stackoverflow.com/questions/38768996/how-to-install-pandas-for-python-3
    #sudo apt-get install python3-pip
    #sudo -H pip3 install pandas
    
    #import the modules
    import numpy as np
    import pandas as pd
    #following the examples on this website
    #https://pandas-docs.github.io/pandas-docs-travis/getting_started/dsintro.html
    # make data frame from a list of dicts
    ls_dcts_data02=[dct_species,dct_fam,dct_seqstrandgb,dct_lenseq01,dct_taxno,dct_vouchergb,dct_product, dct_pubjournal, dct_pubauthor, dct_pubtitle,dct_pubdoi, dct_location, dct_latlonpos]
    #ls_dcts_data02=[dct_species,dct_fam,dct_seqstrandgb,dct_lenseq01,dct_taxno,dct_vouchergb,dct_product]
    #ls_dcts_data02=[dct_species,dct_fam,dct_seqstrandgb,dct_lenseq01,dct_taxno]
    # make it a data frame with pandas
    df_data02=pd.DataFrame(ls_dcts_data02)
    # transpose this dataframe
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.transpose.html
    df_data03 = df_data02.T
    #change column names on panda data frame. See this website: https://stackoverflow.com/questions/11346283/renaming-columns-in-pandas
    df_data03.columns = ['speciesgenus', 'family', 'ntseq','length_of_ntseq','taxno','voucher','geneproduct', 'pubjournal', 'pubauthor', 'pubtitle','pubdoi','location','latlonpos']
    #df_data03.columns = ['speciesgenus', 'family', 'ntseq','length_of_ntseq','taxno','voucher','geneproduct']
    # Get rows from pandas with max values, 
    # I used the second answer , which sorts by the the value of 'length_of_ntseq', and then 
    # drop any duplicates in 'speciesgenus'
    #https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-value-in-groups-using-groupby
    df_data04=df_data03.sort_values('length_of_ntseq', ascending=False).drop_duplicates(['speciesgenus'])
    
    #Remove specific columns of the pandas dataframe:
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html
    df_data05=df_data04.drop(columns=[ 'family', 'ntseq','length_of_ntseq','taxno'])
    # turn the data frame with only two columns in to a dictionary
    # see this website: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    df_data06=df_data05.T
    #convert the dataframe to a dictionary - with parameters 'index' as here: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    df_data07=df_data06.to_dict('index')
    #as the above creates a nested dictionary, access the nested dictionary entitled 'speciesgenus'
    dct_accno_specs_for_custom_BLASTdb=df_data07['speciesgenus']
    
    #Remove specific columns of the pandas dataframe:
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html
    df_data08=df_data04.drop(columns=['speciesgenus','family', 'length_of_ntseq','taxno'], labels=None)
    #transpose data frame to be able to make it a dictionary
    df_data09=df_data08.T
    #convert the dataframe to a dictionary - with parameters 'index' as here: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    df_data10=df_data09.to_dict('index')
    #as the above creates a nested dictionary, access the nested dictionary entitled 'speciesgenus'
    dct_accno_ntseq_for_custom_BLASTdb=df_data10['ntseq']
    # turn the data frame with only two columns in to a dictionary
    #Remove specific columns of the pandas dataframe:
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html
    df_data11=df_data04.drop(columns=['speciesgenus','family', 'ntseq','length_of_ntseq'])
    # turn the data frame with only two columns in to a dictionary
    # see this website: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    # see this website: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    df_data12=df_data11.T
    #convert the dataframe to a dictionary - with parameters 'index' as here: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
    df_data14=df_data12.to_dict('index')
    #as the above creates a nested dictionary, access the nested dictionary entitled 'speciesgenus'
    dct_accno_taxno_for_custom_BLASTdb=df_data14['taxno']
    
    #The data frame prepared can be written as a csv file
    #get module to write csv files
    #https://stackoverflow.com/questions/16923281/writing-a-pandas-dataframe-to-csv-file
    import csv
    #define a file name that the resulting data frame is to be written out to
    file_name="output01_df_w_ntseq_for_customBLASTdb_"+var01+"_"+prod02+".csv"
    #write data frame to csv
    df_data04.to_csv(file_name, sep=';', encoding='utf-8')
    
    #write the nucleotide sequence dictionary to a file
    # see this example: https://stackoverflow.com/questions/36965507/writing-a-dictionary-to-a-text-file
    #Open a file and use 'a' to append to a file, use 'w' to write
    with open(fileP_05_06_accn_ntseq, 'w') as f:
        for key, val in dct_seqstrandgb.items():
            f.write(str(key)+'\t'+str(val)+'\n')
    #close the file again
    f.close()

    #write a file with accession numbers matching publications
    #make data frame from a list of dicts
    ls_dcts_data15=[dct_pubjournal, dct_pubauthor, dct_pubtitle,dct_pubdoi, dct_location,dct_latlonpos]
    df_data16=pd.DataFrame(ls_dcts_data15)
    # transpose this dataframe
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.transpose.html
    df_data17 = df_data16.T
    #change column names on panda data frame. See this website: https://stackoverflow.com/questions/11346283/renaming-columns-in-pandas
    df_data17.columns = ['pubjournal', 'pubauthor', 'pubtitle','pubdoi','location','latlonpos']
    #The data frame prepared can be written as a csv file
    #get module to write csv files
    #https://stackoverflow.com/questions/16923281/writing-a-pandas-dataframe-to-csv-file
    import csv
    #define a file name that the resulting data frame is to be written out to
    file_name=fileP_05_07_accn_publication
    #write data frame to csv
    df_data17.to_csv(file_name, sep=';', encoding='utf-8')


    #______________________________________________________________________________________________________________________________
    #  section 09 end : Build dictionaries from the files and dictionaries to prepare a csv table with all information and sequences
    #______________________________________________________________________________________________________________________________



#______________________________________________________________________________________________________________________________
#  section 10 start : Use the csv file and dictionaries to build a local BLAST database
#______________________________________________________________________________________________________________________________

# # The code below is set up with the goal of making input files to match the format described here:
# #https://www.ncbi.nlm.nih.gov/books/NBK279688/

# #like this:
# #$ cat test.fsa 
# #>seq1
# #MSFSTKPLDMATWPDFAALVERHNGVWGGCWCMAFHAKGSGAVGNREAKEARVREGSTHAALVFDGSACVGWCQFGPTGE
# #LPRIKHLRAYEDGQAVLPDWRITCFFSDKAFRGKGVAAAALAGALAEIGRLGGGTVESYPEDAQGRTVAGAFLHNGTLAM
# #>seq2
# #MKAIDLKAEEKKRLIEGIQDFFYEERNEEIGIIAAEKALDFFLSGVGKLIYNKALDESKIWFSRRLEDISLDYELLYK
# #>seq3 
# #MTLAAAAQSATWTFIDGDWYEGNVAILGPRSHAMWLGTSVFDGARWFEGVAPDLELHAARVNASAIALGLAPNMTPEQIV
# #GLTWDGLKKFDGKTAVYIRPMYWAEHGGYMGVPADPASTRFCLCLYESPMISPTGFSVTVSPFRRPTIETMPTNAKAGCL
# #YPNNGRAILEAKARGFDNALVLDMLGNVAETGSSNIFLVKDGHVLTPAPNGTFLSGITRSRTMTLLGDYGFRTTEKTLSV
# #RDFLEADEIFSTGNHSKVVPITRIEGRDLQPGPVAKKARELYWDWAHSASVG
# #>seq4
# #MRSFFHHVAAADPASFGVAQRVLTIPIKRAHIEVTHHLTKAEVDALIAAPNPRTSRGRRDRTFLLFLARTGARVSEATGV
# #NANDLQLERSHPQVLLRGKGRRDRVIPIPQDLARALTALLAEHGIANHEPRPIFIGARQERLTRFGATHIVRRAAAQAVT
# #IKPALAHKPISPHIFRHSLAMKLLQSGVDLLTIQAWLGHAQVATTHRYAAADVEMMRKGLEKAGVSGDLGLRFRPNDAVL
# #QLLTSI
# #>seq5
# #MTISRVCGSRTEAMLTNGQEIAMTSILKSTGAVALLLLYTLTANATSLMISPSSIERVAPDRAAVFHLRNQMDRPISIKV
# #RVFRWSQKGGVEKLEPTGDVVASPISAQLSPNGNRAVRVVRVSKEPLRSEEGYRVVIDEADPTRNTPEAESLSARHVLPV
# #LFRPPDVLGPEIELSLTRSDGWLMLVVENKGASRLRRSDVTLAQGSAGIARREGFVGYVLPGLTRHWRVGREDSYSGGIV
# #TVSANSSGGAIGEQLVVSGR
# #>seq6
# #TTLLLQVPIGWGVLHQGGALVVLGFAIAHWRGFVGTYTRDTAIEMRD
# #
# #An additional (optional) file mapping the identifiers to taxids (a number identifying a taxonomic node) may be used to associate each sequence with a taxonomic node.
# #
# #$ cat test_map.txt
# #seq1 68287
# #seq2 2382161
# #seq3 68287
# #seq4 382
# #seq5 382
# #seq6 382


# #https://www.biostars.org/p/314630/
# #Get the module that will allow you to write out a fasta file
# from Bio import SeqIO

# #following the question here
# #https://www.biostars.org/p/314630/
# #convert the dictionary in to a file that the 'SeqIO.convert' part can read

# #first write the dictionary to a tab delimited file
# # get the csv module - this will be useful for writing out the tab delimited file that has tax ids matching accession numbers
# import csv
# #write the new dictionary to a csv file
# #https://www.quora.com/How-do-I-write-a-dictionary-to-a-file-in-Python
# with open('dct_accno_ntseq_for_custom_BLASTdb01.csv', 'w') as f:
#     for key, value in dct_accno_ntseq_for_custom_BLASTdb.items():
#         f.write('%s\t%s\n' % (key, value))

# #use the tab delimited file prepared above as input, and make the output format 'fasta' format
# #also note the delimeter 'tab' is called 'tab' . Not the regex '\t'.
# SeqIO.convert('dct_accno_ntseq_for_custom_BLASTdb01.csv', 'tab', 'dct_accno_ntseqs_for_custom_BLASTdb02.fas', 'fasta')
# #the resulting 'dct_seqstrandgb.output4.fas' file is now a fasta file with accession numbers as headers and the matching sequences

# #import the module that will allow you to replace regex
# import re

# #make an empty dictionary
# dct_taxno02 = {}
# #iterate over both keys and values for the items( items covers both keys and values in a dictionary) 
# #and do different things fot the keys and the values
# for key, value in dct_accno_taxno_for_custom_BLASTdb.items():
#     #make the value a string to make it possible to replace regex on the value - 
#     #here I match 'taxon:' and replace only 'taxon:'
#     val01 = str(value)
#     val02 = re.sub(r"taxon:","",val01)
#     # I used print() during the check of this loop to see if I got the part replaced that I aimed to replace
#     #print(val02)
#     # for the keys assign them to the dictionary, and make them match the modified value
#     dct_taxno02[key] = val02
# # write out the new dictionary as a tab delimited file - this should match the file format described here : https://www.ncbi.nlm.nih.gov/books/NBK279688/
# # - this will be a tab separated file with acc numbers as keys and NCBItaxID numbers names as values
# with open('dct_accno_taxno_for_custom_BLASTdb01.txt', 'w') as f:
#     for key, value in dct_taxno02.items():
#         f.write('%s\t%s\n' % (key, value))

# # make another dictionary to hold acc numbers that match species names
# dct_specsdict02 = {}
# #iterate over both keys and values for the items( items covers both keys and values in a dictionary) 
# # and do different things on the values
# for key, value in dct_accno_specs_for_custom_BLASTdb.items():
#     #make the value a string to make it possible to replace regex on the value - here I match ' ' and replace with '_'
#     val01 = str(value)
#     val02 = re.sub(r" ","_",val01)
#     # I used print() during the check of this loop to see if I got the part replaced that I aimed to replace
#     #print(val02)
#     # for the keys assign them to the dictionary, and make them match the modified value
#     dct_specsdict02[key] = val02

# # write out the a second new dictionary - this will be a tab separated file with acc numbers as 
# #keys and species names as values
# with open('dct_accno_specs_for_custom_BLASTdb.txt', 'w') as f:
#     for key, value in dct_specsdict02.items():
#         f.write('%s\t%s\n' % (key, value))

# # write out a list to a txt file
# # this was to be able to try out the code above on a list of species
# # Use the list 'ls_organisms'
# #https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python
# with open('ls_organisms01.txt', 'w') as f:
#     for item in ls_organisms:
#         f.write("%s\n" % item)
# ######################################################################################################################################################
# ######################################################################################################################################################
# #  #In a unix terminal you should now be able to first load the BLAST module
# #module load BLAST

# #  # To check what programs are available when you run the command 'module load BLAST'  . You can check if you can find the programs :  'makeblastdb' and 'blastdbcp'. These can be found by checking with:
# #ls $(dirname $(which blastn))
# #  # ie: list the files which are in the same directory as blastn.

# #  #Also compare with the example provided here https://www.ncbi.nlm.nih.gov/books/NBK279688/
# #  #that makes use of this line :
# #makeblastdb -in test.fsa -parse_seqids -blastdb_version 5 -taxid_map test_map.txt -title "Cookbook demo" -dbtype prot
# #  #to build a local database
# #  
# #  #instead try this line -  notice that the 'dct_seqstrandgb.output4.fas' replaces the 'test.fsa' and that 'accno_to_taxid.txt' replaces 'test_map.txt'
# #  #also notice that '-dbtype nucl' is nucl not 'prot' because I am building a nucleotide mathcing database, not a protein matching database
# # 
# #  makeblastdb -in /home/sknu003/uoa00029_runs/local_BLAST_DB/dct_seqstrandgb.output4.fas -parse_seqids -blastdb_version 5 -taxid_map /home/sknu003/uoa00029_runs/local_BLAST_DB/accno_to_taxidno.txt -title "Syngn_Limand_DB_v01" -dbtype nucl

# # Here is second example with some other file names for input, but the idea is the same as above.

# #  makeblastdb -in /home/sknu003/uoa00029_runs/local_BLAST_DB/dct_accno_ntseqs_for_custom_BLASTdb02.fas -parse_seqids -blastdb_version 5 -taxid_map /home/sknu003/uoa00029_runs/local_BLAST_DB/dct_accno_taxno_for_custom_BLASTdb01.txt -title "custDB_Syngnathid_Limanda_v02" -dbtype nucl# 

# #  #also see this website: https://www.biostars.org/p/76551/
# #  #IMPORTANT! Note that including the '-taxid_map' argument to a file that is able to match the taxID number with an accession number will be required 
# #  #for the later blast search query against your local database. This argument and file is required to get the species names with sscinames (or ssciname)  argument in the blast search
# #  #But before you can start blasting against your own local database you will need a copy of this file : ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# #  #download it with the wget command like this:
# #  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# #  #Then uncompress the file like this:
# #  tar -zxvf taxdb.tar.gz
# #  #afterwards I tried making this small fasta file:
# #  #>partofAF414201 CCCGGGTACTACGAGCACCAGCTTAAAACCCAAAGGACTTGGCGGTGCTTTAGATCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCCCGTTAAACCTCACCTTTTCTTGTTTTCCCCGCCTATATACCGCCGTCGTCAGCTTACCCTGTG#AAGGACCAATAGTAAGCAAAATTGGCACAGCCCAAAACGCCAGGTCGAGGTGTAGCGTATGGAAAGGGAAGAAATGGGCTACATTCAATATAACAATGCATACGGATGGTGGCCTGAAATAACCACCTGAAGGAGGATTTAGCAGTAAGCAGGA#AATAGAGAGTC
# #  #and named it : 'partofAF414201.fas'
# #  #And since I already have loaded the BLAST module
# #  module load BLAST
# #  #I tried this:
# # blastn -db dct_seqstrandgb.output4.fas -query partofAF414201.fas -out results_partofAF414201.fas.out
# # 
# # 
# #  #following some of the instructions here:
# #  #https://www.biostars.org/p/76551/
# #  #I discovered that after I have downloaded the ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz file with wget, and after I had uncompressed the file with 
# #  tar -zxvf taxdb.tar.gz
# #  #I was able to include some of the settings in the '-outfmt' as in the example on this webpage: https://www.biostars.org/p/76551/
# # 
# #  blastn -db dct_seqstrandgb.output4.fas -query partfrom_AY141392.fas -out results_partofAF414201.fas.out -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle'
# # 
# #  #In this query the input file 'partfrom_AY141392.fas' only hold one single nt-sequence already obtained from NCBI -  i.e. from 'Labrus bergylta', with the '-outfmt' argument I was able to get the file : 'results_partofAF414201.fas.out' which then look like this :
# # 
# #  #partfrom_AY141392       gb|AF414201|    0.0     652     0       AF414201        56723   Labrus bergylta ballan wrasse   N/A
# #  #partfrom_AY141392       gb|AY141392|    0.0     652     0       AY141392        56723   Labrus bergylta ballan wrasse   N/A
# #
# #  #which is what I expected.
# # 
# #  #I also tried following the instructions here: https://www.ncbi.nlm.nih.gov/books/NBK279680/
# #  #and here:  https://www.biostars.org/p/304824/
# #  #This of course threw an unforeseen error. Looking up this error lead me here: http://www.metagenomics.wiki/tools/blast/error-too-many-positional-arguments
# #  #this line should work:
# #  blastn -db nt -query partofAF414201.fas -out results.out -remote
# # 
# #  #but this one does not work :
# #  blastn –db nt –query partofAF414201.fas –out results.out -remote
# #  #The lines ('–' versus '-') are different -  how terribly annoying!
# #  #Another weblog : https://danielbruzzese.wordpress.com/2018/03/26/installing-and-querying-a-local-ncbi-nucleotide-database-nt/


#______________________________________________________________________________________________________________________________
#  section 10 end : Use the csv file and dictionaries to build a local BLAST database
#______________________________________________________________________________________________________________________________











####