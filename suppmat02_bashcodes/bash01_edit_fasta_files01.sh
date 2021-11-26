#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)
# get base directory
WD00="/home/hal9000/MS_Mnemiopsis/"
#WD00="/Users/steenknudsen/Documents/Documents/MS_DNA_og_Liv_citizen_science/"
#get dir with the present bash code inside
WD02="suppmat02_bashcodes"
#define output directory
WD04="suppmat04_edited_fasta_fl"
# #define input dirs
WD01="suppmat01_inp_files"
# WD05="supma05_qpcr_csv_data_from_students_tests"
# in Finder on Macintosh connect to
# https://webfile.science.ku.dk/webdav

#make strings with new dirs
WD00_WD02=$(echo "$WD00""$WD02")
WD00_WD04=$(echo "$WD00""$WD04")
WD00_WD01=$(echo "$WD00""$WD01")
# WD00_WD05=$(echo "$WD00""$WD05")
# WD00_WD04=$(echo "$WD00""$WD04")
echo "$WD00_WD02"
echo "$WD00_WD04"

# OUTDIR1 is the dir where the final merged csv file is to be stored
OUTDIR1="${WD00_WD04}"
INPDIR1="${WD00_WD01}"
#define outfile
OUTFILE1="algn_Mnelei_18s_10.fas"
#define inputfile
INPFILE1="algn_Mnelei_18s_09.fas"
#remove the old versions of the output directory
rm -rf "${OUTDIR1}"
#make new versions of the in- and output directory
mkdir "${OUTDIR1}"

cd "$WD00_WD01"
cat "$INPFILE1" | \
	# use sed to replace
	sed "s/_consensus_sequence//g" | \
	# use sed to replace underscores at end of line
	sed "s/_$//g" | \
	sed "s/>Mnemiopsis_leidyi/>Mnelei/g" | \
	head -22 | cut -c1-80 

