# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 16:34:13 2016

@author: sarahguiziou
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'

"""
DESIGN DNA SEQUENCE

Generate a gb file with the DNA sequence of a computational device.
The computational device implemented in a specific strain will permit the partial implementation
of a specific Boolean function.

Input:
1: inp, define the type of computational device to generate. It is composed of 2 numbers,
with the first number that corresponds to the number of integrase needed for this device,
and the second number corresponds to the number of IMPLY function.

2: nbr_strain; integer that corresponds to the ID number of this strain.
3: output; binary number that corresponds to the input Boolean function.
4: mypath: the path for the result directory
5: directory_name: the name of the result directory

"""

def design_DNAsequence(inp, nbr_strain, output, path, name_directory):
  
    ###############################################################################
    # SEQUENCE DEFINITION
    
    name_int=['Bxb1','Tp901','Int5','Int7','Int4','Int3']
    
    # define attB sites
    BF_site=[]
    # Bxb1=1
    attB_Bxb1='TCGGCCGGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCATCCGGGC'
    BF_site.append(attB_Bxb1)
    # Tp901=2
    attB_Tp901='ATGCCAACACAATTAACATCTCAATCAAGGTAAATGCTTTTTGCTTTTTTTGC'
    BF_site.append(attB_Tp901)
    #Int5=3
    attB_Int5='GAGCGCCGGATCAGGGAGTGGACGGCCTGGGAGCGCTACACGCTGTGGCTGCGGTCGGTGC'
    BF_site.append(attB_Int5)
    #int7=4
    attB_Int7='AGACGAGAAACGTTCCGTCCGTCTGGGTCAGTTGGGCAAAGTTGATGACCGGGTCGTCCGTT'
    BF_site.append(attB_Int7)
    #int4=5
    attB_int4='ttccaaagagcgcccaacgcgacctgaaatttgaataagactgctgcttgtgtaaaggcgatgatt'
    BF_site.append(attB_int4)
    #int3=6
    attB_int3='gtttgtaaaggagactgataatggcatgtacaactatactcgtcggtaaaaaggcatcttat'
    BF_site.append(attB_int3)
      
    # define attP sites  
    PF_site=[]
    # Bxb1=1
    attP_Bxb1='TCGTGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCC'
    PF_site.append(attP_Bxb1)
    # Tp901=2
    attP_Tp901='GCGAGTTTTTATTTCGTTTATTTCAATTAAGGTAACTAAAAAACTCCTTT'
    PF_site.append(attP_Tp901)
    #Int5=3
    attP_Int5='CCCTAATACGCAAGTCGATAACTCTCCTGGGAGCGTTGACAACTTGCGCACCCTGATCTG'
    PF_site.append(attP_Int5)
    #Int7=4
    attP_Int7='GTGTTATAAACCTGTGTGAGAGTTAAGTTTACATGCCTAACCTTAACTTTTACGCAGGTTCAGCTT'
    PF_site.append(attP_Int7)
    #int4=5
    attP_Int4='caaaaattacaaagttttcaacccttgatttgaattagcggtcaaataatttgtaattcgttt'
    PF_site.append(attP_Int4)
    #int3=6
    attP_Int3='atggataaaaaaatacagcgtttttcatgtacaactatactagttgtagtgcctaaataatgctt'
    PF_site.append(attP_Int3)    
               
    # define terminators
    
    name_ter=['ECK120033737', 'B0015', 'L3S2P21', 'L3S3P21', 'L3S2P11', 'L3S3P22']
    ter=[]
    
    ECK120033737='GGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGG'
    T1=ECK120033737
    ter.append(T1)
    
    B0015='ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata'
    T2=B0015
    ter.append(T2)
    
    L3S2P21='CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC'
    T3=L3S2P21
    ter.append(T3)
    
    L3S3P21='CCAATTATTGAAGGCCTCCCTAACGGGGGGCCTTTTTTTGTTTCTGGTCTCCC'
    T4=L3S3P21
    ter.append(T4)
    
    L3S2P11='CTCGGTACCAAATTCCAGAAAAGAGACGCTTTCGAGCGTCTTTTTTCGTTTTGGTCC'
    T5=L3S2P11
    ter.append(T5)

    L3S3P22='CCAATTATTGAAGGCCGCTAACGCGGCCTTTTTTTGTTTCTGGTCTCCC'
    T6=L3S3P22
    ter.append(T6)
    
    # define spacers, promoters, BCD, RiboJ...
    # P2
    P7f='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGA'
    RiboJ='AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA'
    BCD2='GGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCTAAGGAGGTTTTCTA'
    sfGFP='atgtcaaaaggagaagaactttttacaggtgtagtacctatcttggttgaattggatggtgatgttaacggtcacaaattttctgtacgtggtgaaggtgaaggtgatgcaactaacggtaaattgacacttaaattcatttgtacaactggaaaacttcctgttccttggcctactcttgttacaacattgacatatggagtacaatgtttttcacgttatcctgatcatatgaaacgtcacgatttttttaaatctgctatgccagaaggttatgtacaagaacgtacaatttcatttaaagatgacggaacatataaaacacgtgctgaagtaaaattcgaaggtgacactcttgttaatcgtatcgaattgaaaggaatcgatttcaaagaagatggtaacattttgggacacaaacttgaatacaacttcaactctcataatgtttatatcacagctgacaaacaaaaaaacggtattaaagctaattttaaaattcgtcacaatgttgaagatggatctgttcaattggctgatcattatcaacaaaatacaccaatcggagacggaccagtattgcttccagataaccactacctttctactcaatcagttctttcaaaagatcctaacgaaaaacgtgaccatatggtacttcttgaatttgttacagcagcaggtatcactcacggtatggacgaactttataaataa'
    sp0='CTCGGATACCCTTACTCTGTTGAAAACGAATAGATAGGTT'
    spN='ATTATTGACCACTTCCGAGTAGAATCGTGCTTCAGTAAGA'
    L3S3P00='CCAATTATTGAAGGGGAGCGGGAAACCGCTCCCCTTTTTTTGTTTCTGGTCTCCC'
    #sp20_1='TAGTTGCGTCTCAGGGACCC'
    #sp20_2='TAAGTGGCAATCCCGCCTGA'
    #sp20_3='AAACCCGTCGCAGTATCCCT'
    #sp20_4='ACTCAGGTCTGCCGTAAGGG'
        
    ##############################################################################
    # SEQUENCES OF THE DESIGN
    
    # initialization of list sequence
    # it will be incremented with the DNA sequence and name of the feature that composed the full final DNA sequence.
    sequence=[]
    
    # append spacer sp0
    sequence.append([sp0,'spacer sp0'])
    
    # number of ZERO correspond to number of NOT elements: total variables (0) - number of ONE (1)
    ZERO=int(inp[0])-int(inp[1])
    
    # loop to generate all attB site for the NOT elements
    # in reverse order between the number of NOT elements to 0
    # Z corresponds to the indice of integrase to implement the NOT element
    for Z in reversed(range(0, ZERO)):
        # append to the design the attB site of NOT elements
        sequence.append([BF_site[Z],'attB '+name_int[Z]])
        
    # append to the design the promoter P7 in fwd direction        
    sequence.append([P7f,'promoter P7'])
    
    # loop to generate all attP site for the NOT elements
    # between 0 to the number of NOT elements
    # Z corresponds to the indice of integrase to implement the NOT element      
    for Z in range(0,ZERO):
        # append attP site
        sequence.append([PF_site[Z],'attP '+name_int[Z]])

    # loop to generate all IMPLY elements
    # O corresponds to the indice of integrase to implement the IMPLY element
    for O in range(ZERO, int(inp[0])):
        # add attB site+term+attP site
        sequence.append([BF_site[O],'attB '+name_int[O]])
        sequence.append([ter[O],name_ter[O]])
        sequence.append([BF_site[O],'attP '+name_int[O]])
        
    # append RiboJ, BCD2, sfGFP, spacer N at the end of the sequence
    sequence.append([RiboJ, 'RiboJ'])
    sequence.append([BCD2, 'BCD2'])
    sequence.append([sfGFP, 'sfGFP'])
    sequence.append([spN, 'spN'])
    sequence.append([L3S3P00, 'Term L3S3P00'])
    
    # initialization of the DNA sequence
    DNA_seq=''
    
    # loop in sequence list to generate the string corresponding to the DNA sequence of the device.
    for seq in sequence:
        DNA_seq+=seq[0]

    # creation of the formated DNA sequence for the genbank file
    seq_final=Seq(DNA_seq, IUPAC.unambiguous_dna)
    record = SeqRecord(seq_final,
                   id='NA', # random accession number
                   name='Strain'+nbr_strain+'_'+inp,
                   description='Sequence of the computational device in the strain'+nbr_strain+' to implement the Boolean function required')
    
    # initialization of the variable len_seq 
    len_seq=0
    
    # loop in sequence list to generate the genbank feature of the sequence
    for feat in sequence:
        
        feature = SeqFeature(FeatureLocation(start=len_seq, end=len_seq+len(feat[0])), id=feat[1], type=feat[1])
        record.features.append(feature)
        len_seq+=len(feat[0]) 

    # Save as GenBank file
    output_file = open(path+'/'+name_directory+'_Strain'+nbr_strain+'.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')
    output_file.close()

   
"""TEST of the design DNA function"""
    
#design_DNAsequence('63','3','0000101010','../results/0/','0')    
