from Bio import SeqIO, Entrez
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

#take in NCBI accesion number for a nucleotide record and download the genbank file
def download_genbank(accession):
    Entrez.email = "bhua@gritstone.com" #  enter your email here

    handle=Entrez.efetch(db='nucleotide',id=accession,rettype='gbwithparts',retmode='text')
    f = open(str(accession)+'.gb', 'w') # save to disk
    f.write(handle.read())
    handle.close()
    f.close()

#take parse genbank file for     
def extract_RNA(accession):
    #read in genbank file
    record = SeqIO.parse(open(accession+'.gb'), "genbank").__next__()
    
    #take in SeqRecord object and extract gene name, protein id, protein product, -30:+30 sequence, +1:+60 sequence, and +31:+90 sequence and write out to csv
    h = open(accession+'_RNA_info.csv', "w")
    h.write('feature_index,gene,protein_id,product,-30:+30,+1:+60,+31:+90\n')
    
    #initialize index
    count=0
    
    for feature in record.features:
        if feature.type == "CDS":
            h.write("%s," % str(count))
            h.write("%s," % feature.qualifiers["gene"][0])
            try:
                h.write("%s," % feature.qualifiers["protein_id"][0])
            except:
                h.write(",")
            
            if 'product' in feature.qualifiers.keys():
                h.write("%s," % feature.qualifiers["product"][0].replace(",",""))
            else:
                h.write(",")
            
            #extract sequences around start codon using coordinates dependent of stand on which gene lies
            if feature.strand == 1:
                h.write("%s," % record.seq[feature.location.nofuzzy_start-30:feature.location.nofuzzy_start+30])
                h.write("%s," % record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_start+60])
                h.write("%s\n" % record.seq[feature.location.nofuzzy_start+31:feature.location.nofuzzy_start+91])
            else:
                h.write("%s," % record.seq[feature.location.nofuzzy_end-30:feature.location.nofuzzy_end+30].reverse_complement())
                h.write("%s," % record.seq[feature.location.nofuzzy_end-60:feature.location.nofuzzy_end].reverse_complement())
                h.write("%s\n" % record.seq[feature.location.nofuzzy_end-90:feature.location.nofuzzy_end-30].reverse_complement())
            count+=1
        else:
            count+=1
    h.close()


#clean RNA and iBAQ files and merge
def merge_iBAQ(accession):
    #clean up RNA_info file
    RNA_df=pd.read_csv(accession+'_RNA_info.csv')
    RNA_df.drop('protein_id',axis=1,inplace=True)
    RNA_df.drop('feature_index',axis=1,inplace=True)
    RNA_df.drop('product',axis=1,inplace=True)
    RNA_df.dropna(inplace=True)
    RNA_df.drop_duplicates(inplace=True)

    RNA_df.to_csv(accession+'_RNA_info_cleaned.csv',index=False)

    #read in and merge HEK-293 iBAQ data
    HEKiBAQ_df=pd.read_csv('./iBAQ_input/HEK_iBAQ_avg_cleaned.csv')
    joined_genes=HEKiBAQ_df.merge(RNA_df,how='left',left_on='gene_names', right_on='gene')
    joined_genes.dropna(inplace=True)
    joined_genes.to_csv('HEK_iBAQ_with_RNA_info.csv',index=False)
    
    return joined_genes
    
#generate fasta of all RNA snippets
def write_fasta(joined_genes):
    with open(f'-30:+30.fasta','w') as file:
        count=1
        for entry in joined_genes['-30:+30'].tolist():
            file.write(f'>seq{count}\n{entry}\n')
            count+=1

    with open(f'+1:+60.fasta','w') as file:
        count=1
        for entry in joined_genes['+1:+60'].tolist():
            file.write(f'>seq{count}\n{entry}\n')
            count+=1

    with open(f'+31:+90.fasta','w') as file:
        count=1
        for entry in joined_genes['+31:+90'].tolist():
            file.write(f'>seq{count}\n{entry}\n')
            count+=1
            
#read in fasta and calculate MFE with RNAfold
def calculate_MFE(joined_genes,fasta,col_name):
    path=os.getcwd()
    #peform RNAfold on fasta file
    output=os.popen(f'cd && Library/ViennaRNA2/bin/RNAfold<{path}/{fasta}').read() #change this to the path where the RNAfold executable is installed
    
    #reformat output from stdout and convert to list
    reformat=output.split('\n')
    MFE_text_list=[]
    for i in range(2,len(reformat),3):
        MFE_text_list.append(reformat[i])
    MFE_floats=[]
    for i in MFE_text_list:
        MFE=''
        for char in i[::-1]:
            if char==')':
                pass
            elif char!='(':
                MFE+=char
            elif char=='(':
                break
        MFE_floats.append(float(MFE[::-1]))
    joined_genes[col_name]=MFE_floats
    joined_genes.to_csv('HEK_iBAQ_with_RNA_info_and_MFE.csv',index=False)
    return joined_genes

#generate plots
def generate_plots(joined_genes):
    #make new output folder for figures if it doesn't already exist
    if os.path.exists('./plots')==False:
        os.mkdir('plots')
    
    regions=['MFE_-30:+30','MFE_+1:+60','MFE_+31:+90']
    for region in regions:
        #generate scatterplot
        plt.figure()
        fig=sns.scatterplot(x=region,y='iBAQ_avg',data=joined_genes,alpha=0.1)
        plt.savefig('./plots/'+region+'_scatter.jpg')
    
        #generate kdeplot
        plt.figure()
        fig2=sns.kdeplot(x=region,data=joined_genes,shade=True)
        plt.savefig('./plots/'+region+'_kdeplot.jpg')
