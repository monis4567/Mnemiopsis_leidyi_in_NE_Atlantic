#!/bin/bash
# -*- coding: utf-8 -*-
echo $BASH_VERSION

##get the present directory
WD=$(pwd)

INF01="NCBI_deposti_seq_Mnemiopsis_leydii_18s.fasta"
INF01="algn_Mnelei_18s_18.fas.aligned.fasta"
INF01="Bankit_NCBI_submiss_file_Mnelei_18s.fasta"
#INF01="Mnemiopsis_18s_2024apr15_v01.fas"
#remove previous versions of feature table
rm -rf *Mne*feature_table*txt
# loop over all fas files in the current directory
for file in  "$INF01"
do
    nf=$(echo $file | sed 's/\.fasta/_feature_table_03.txt/g')
    echo $nf

    # read individual fasta file line by line
    while read line
    do
        if [[ $line == ">"* ]] # if line starts with '>'
        then
            # Use awk to retain first element in line using '_' as field seperator 
            # and write to the new file name
            #echo $line | awk 'BEGIN { FS = "_"; OFS="_" } ; {print $1,$2,$3,$4,$5}' |\
            echo $line | awk 'BEGIN { FS = "_" } ; {print $1}' |\
            sed -E 's/>/>Feature /g' >> $nf
        else # if line doesn't start with >, then simply just echo the line
            SEQL=$(echo $line | wc -c)
            #echo "$SEQL" 
            echo "" >> $nf
            printf "<1    >"$SEQL"    gene
                                    gene          18Sr
            <1    >"$SEQL"    mRNA
                                    product       18Sr ribosomal RNA

                    <1    >"$SEQL"    gene
                                    gene          its1
            <1    "$SEQL"    CDS
                                    product       internal transcribed spacer 1
                                    product       internal transcribed spacer 1
                                    codon_start   2
            <1    >"$SEQL"    mRNA
                                    product       internal transcribed spacer 1
                   	<1    >"$SEQL"    gene
                                    gene          5.8Sr
            <1    >"$SEQL"    mRNA
                                    product       5.8S ribosomal RNA

                    <1    >"$SEQL"    gene
                                    gene          its2
            <1    "$SEQL"    CDS
                                    product       internal transcribed spacer 2
                                    product       internal transcribed spacer 2
                                    codon_start   1
            <1    >"$SEQL"    mRNA
                                    product       internal transcribed spacer 1
                    
                    <1    >"$SEQL"    gene
                                    gene          28Sr ribosomal RNA
            
            <1    >"$SEQL"    mRNA
                                    product       28Sr
                                    " >> $nf
            echo "" >> $nf
        fi
    #end the iteration over the file, and define the input file to do the
    # iteration over
 	done < ${file}
done


# cat ${nf} | sed 's/\t/    /g' > tmp.txt
# mv tmp.txt ${nf}

head -30 $nf

# Example
# https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
# >Feature Seq1
# <1    >1050    gene
#                         gene          ATH1
# <1    1009    CDS
#                         product       acid trehalase
#                         product       Athlp
#                         codon_start   2
# <1    >1050    mRNA
#                         product       acid trehalase

# >Feature Seq2
# 2626  2590    tRNA
# 2570  2535
#                         product       tRNA-Phe

# >Feature Seq3
# 1080  1210  CDS
# 1275  1315
#                         product       actin
#                         note          alternatively spliced
# 1055  1210  mRNA
# 1275  1340
#                         product       actin
# 1055  1340  gene
#                         gene          ACT
# 1055  1079  5'UTR
# 1316  1340  3'UTR
