import click
from RNAiBAQ_functions import *

@click.command
@click.argument('accession')

def main(accession):
    
    """
    This script takes an NCBI accession ID then downloads the entry as a genbank and saves it to local storage, then parses the genbank to extract the sequences of the -30:+30, +1:+60, and +31:+90 regions of the start codon. It then joins the iBAQ mass spec average values based on gene name. Finally, the lowest minimum free energy is calculated for each RNA snippet using RNAfold and the data is outputed into a final csv file.
    """    
    
    download_genbank(accession)
    extract_RNA(accession)
    joined_genes = merge_iBAQ(accession)
    write_fasta(joined_genes)
    joined_genes=calculate_MFE(joined_genes,'-30:+30.fasta','MFE_-30:+30')
    joined_genes=calculate_MFE(joined_genes,'+1:+60.fasta','MFE_+1:+60')
    joined_genes=calculate_MFE(joined_genes,'+31:+90.fasta','MFE_+31:+90')
    generate_plots(joined_genes)
    
if __name__ == '__main__':
    main()