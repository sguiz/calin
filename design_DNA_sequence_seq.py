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
import string

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

def design_DNAsequence(nb_input, list_integrase, list_gene_state, nbr_strain, output, path, name_directory):
  
    ###############################################################################
    # SEQUENCE DEFINITION
  
    # To reverse complement sequences
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)  
    
    
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
    attB_int4='TTCCAAAGAGCGCCCAACGCGACCTGAAATTTGAATAAGACTGCTGCTTGTGTAAAGGCGATGATT'
    BF_site.append(attB_int4)
    #int3=6
    attB_int3='GTTTGTAAAGGAGACTGATAATGGCATGTACAACTATACTCGTCGGTAAAAAGGCATCTTAT'
    BF_site.append(attB_int3)
    
    # attB reverse site
    BR_site=[]
    
    for a in range(len(BF_site)):
        BR_site.append(BF_site[a].translate(tab)[::-1])
        
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
    attP_Int4='CAAAAATTACAAAGTTTTCAACCCTTGATTTGAATTAGCGGTCAAATAATTTGTAATTCGTTT'
    PF_site.append(attP_Int4)
    #int3=6
    attP_Int3='ATGGATAAAAAAATACAGCGTTTTTCATGTACAACTATACTAGTTGTAGTGCCTAAATAATGCTT'
    PF_site.append(attP_Int3) 
    
    # attP reverse site
    PR_site=[]
    
    for a in range(len(PF_site)):
        PR_site.append(PF_site[a].translate(tab)[::-1])
    
    # define spacers, promoters, BCD, RiboJ...
    # Promoter
    P6f='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGA'
    P6r=P6f.translate(tab)[::-1]
    
    # Ribozyme
    RiboJf='AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA'
    RiboJr=RiboJf.translate(tab)[::-1]
    BydvJf='AGGGTGTCTCAAGGTGCGTACCTTGACTGATGAGTCCGAAAGGACGAAACACCCCTCTACAAATAATTTTGTTTAA'
    BydvJr=BydvJf.translate(tab)[::-1]
    ElvJf='AGCCCCATAGGGTGGTGTGTACCACCCCTGATGAGTCCAAAAGGACGAAATGGGGCCTCTACAAATAATTTTGTTTAA'
    ElvJr=ElvJf.translate(tab)[::-1]
    AraJf='AGTGGTCGTGATCTGAAACTCGATCACCTGATGAGCTCAAGGCAGAGCGAAACCACCTCTACAAATAATTTTGTTTAA'
    AraJr=AraJf.translate(tab)[::-1]
    
    # RBS
    B0034f='AAAGAGGAGAAA'
    B0034r=B0034f.translate(tab)[::-1]
    spacer_RBSf='TACTAG'
    spacer_RBSr=spacer_RBSf.translate(tab)[::-1]
    
    # gene sequence
    sfGFP='ATGCGTAAAGGCGAAGAGCTGTTCACTGGTGTCGTCCCTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCATATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATGATAA'
    BFP='ATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGACGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAA'
    mKate2='ATGTCAGAATTAATTAAAGAAAATATGCACATGAAATTATATATGGAAGGTACTGTCAACAATCATCATTTCAAATGCACATCCGAAGGTGAAGGTAAACCATATGAAGGCACACAAACAATGCGCATCAAAGCAGTTGAAGGTGGACCCCTGCCCTTTGCGTTTGACATTCTCGCAACGAGCTTTATGTACGGGTCTAAAACTTTTATCAATCACACCCAAGGCATTCCTGACTTTTTTAAACAGTCCTTTCCTGAAGGCTTTACCTGGGAACGTGTAACAACTTATGAAGATGGCGGTGTACTTACAGCAACTCAAGATACGAGTTTACAAGATGGCTGTCTGATTTACAATGTTAAAATCCGTGGCGTAAATTTCCCGAGTAACGGACCCGTAATGCAAAAAAAAACTCTTGGTTGGGAAGCATCAACAGAAACCTTATATCCTGCGGACGGTGGCTTAGAAGGACGCGCAGACATGGCACTGAAATTAGTTGGAGGCGGTCATTTAATCTGCAACCTGAAAACAACCTATCGTTCCAAAAAACCCGCTAAAAACCTTAAAATGCCTGGAGTATACTATGTTGATCGTCGCTTAGAACGTATTAAAGAAGCTGATAAAGAAACCTACGTTGAACAACATGAAGTAGCCGTAGCCCGTTATTGTGACCTTCCGTCGAAATTAGGACATCGTTGATAA'
    LacZ='ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGA'
    
    out_f=[sfGFP, BFP, mKate2, LacZ]
    out_name=['sfGFP', 'BFP', 'mKate2', 'LacZ']
    # attP reverse site
    out_r=[]
    
    for a in range(len(out_f)):
        out_r.append(out_f[a].translate(tab)[::-1])
        
    # spacers
    sp0f='CTCGGATACCCTTACTCTGTTGAAAACGAATAGATAGGTT'
    sp4f='AGGCAACTGAAACGATTCGGATCCTGTATTACTATTCTTA'
    sp5f='ACTTTATCTGAGAATAGTCAATCTTCGGAAATCCCAGGTG'
    sp6f='CCGTCTCAGAATCGGCCGTGAACAATAAAATAGTTTCGGT'
    sp6r=sp6f.translate(tab)[::-1]
    sp7f='TAATAAAAGGTCCCGTCTGAACTTACTGTGAATTCGACTA'
    sp7r=sp7f.translate(tab)[::-1]
    spNf='ATTATTGACCACTTCCGAGTAGAATCGTGCTTCAGTAAGA'

    sp20_1f='TAGTTGCGTCTCAGGGACCC'
    sp20_2f='TAAGTGGCAATCCCGCCTGA'
    sp20_3f='AAACCCGTCGCAGTATCCCT'
    sp20_4f='ACTCAGGTCTGCCGTAAGGG'
    
    # terminator
    L3S3P00f='CCAATTATTGAAGGGGAGCGGGAAACCGCTCCCCTTTTTTTGTTTCTGGTCTCCC'
    L3S3P00r=L3S3P00f.translate(tab)[::-1]
    L3S2P21f='CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC'
    L3S2P21r=L3S2P21f.translate(tab)[::-1]
    B0014f='TCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATATACTAGAGAGAGAATATAAAAAGCCAGATTATTAATCCGGCTTTTTTATTATTT'
    B0014r=B0014f.translate(tab)[::-1]
    J61048f='ccggcttatcggtcagtttcacctgatttacgtaaaaacccgcttcggcgggtttttgcttttggaggggcagaaagatgaatgactgtccacgacgctatacccaaaagaaa'
    J61048r=J61048f.translate(tab)[::-1]
        
    ##############################################################################
    # SEQUENCES OF THE DESIGN
    
    # initialization of list sequence
    # it will be incremented with the DNA sequence and name of the feature that composed the full final DNA sequence.
    sequence=[]
    
    # append spacer sp0
    sequence.append([sp0f,'spacer sp0',1])

    # construction for 2 inputs
    if nb_input==2:
        
        if list_gene_state[2]!=0:
            sequence.append([out_r[list_gene_state[2]-1], out_name[list_gene_state[2]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])
        
        sequence.append([BF_site[1], 'attB '+name_int[1],1])
        sequence.append([P6f, 'P6', 1])
        sequence.append([BF_site[0], 'attB '+name_int[0],1])
        
        if list_gene_state[0]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[0]-1], out_name[list_gene_state[0]-1], 1])
            
        sequence.append([sp4f, 'spacer 4', 1])
        sequence.append([B0014f, 'B0014',1])
        sequence.append([sp5f, 'spacer 5', 1])

        if list_gene_state[1]!=0:
            sequence.append([out_f[list_gene_state[1]-1], out_name[list_gene_state[1]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([PF_site[1], 'attP '+name_int[1],1])
        sequence.append([PR_site[0], 'attP '+name_int[0],-1])

    # construction for 3 inputs
    elif nb_input==3:
        
        if list_gene_state[3]!=0:
            sequence.append([out_r[list_gene_state[3]-1], out_name[list_gene_state[3]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([BF_site[2], 'attB '+name_int[2],1])
        sequence.append([sp7r, 'spacer 7', -1])
        sequence.append([L3S2P21r, 'L3S2P21',-1])
        sequence.append([sp6r, 'spacer 6', -1])
        
        if list_gene_state[0]!=0:
            sequence.append([out_r[list_gene_state[0]-1], out_name[list_gene_state[0]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([BF_site[0], 'attB '+name_int[0],1])
        sequence.append([BR_site[1], 'attB '+name_int[1],-1])
        sequence.append([P6r, 'P6',-1])
        sequence.append([PR_site[0], 'attP '+name_int[0],-1])

        if list_gene_state[1]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[1]-1], out_name[list_gene_state[1]-1], 1])

        sequence.append([sp4f, 'spacer 4', 1])
        sequence.append([B0014f, 'B0014',1])
        sequence.append([sp5f, 'spacer 5', 1])

        
        if list_gene_state[2]!=0:
            sequence.append([out_f[list_gene_state[2]-1], out_name[list_gene_state[2]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([PF_site[2], 'attP '+name_int[2],1])
        sequence.append([PR_site[1], 'attP '+name_int[1],-1])

    # construction for 4 inputs
    elif nb_input==4:

        if list_gene_state[4]!=0:
            sequence.append([out_r[list_gene_state[4]-1], out_name[list_gene_state[4]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([BF_site[3], 'attB '+name_int[3],1])
        sequence.append([BF_site[0], 'attB '+name_int[0],1])
        sequence.append([BR_site[1], 'attB '+name_int[1],-1])

        if list_gene_state[1]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[1]-1], out_name[list_gene_state[1]-1], 1])
        
        sequence.append([sp7r, 'spacer 7', -1])
        sequence.append([L3S2P21r, 'L3S2P21',-1])
        sequence.append([sp6r, 'spacer 6', -1])
        
        if list_gene_state[0]!=0:
            sequence.append([out_r[list_gene_state[0]-1], out_name[list_gene_state[0]-1], -1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([PR_site[0], 'attP '+name_int[0],-1])
        sequence.append([BR_site[2], 'attB '+name_int[2],-1])
        sequence.append([P6r, 'P6',-1])
        sequence.append([PR_site[1], 'attP '+name_int[1],-1])
        
        if list_gene_state[2]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[2]-1], out_name[list_gene_state[2]-1], 1])
        
        sequence.append([sp4f, 'spacer 4', 1])
        sequence.append([B0014f, 'B0014',1])
        sequence.append([sp5f, 'spacer 5', 1])
        
        
        if list_gene_state[3]!=0:
            sequence.append([out_r[list_gene_state[3]-1], out_name[list_gene_state[3]-1],-1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([PF_site[3], 'attP '+name_int[3],1])
        sequence.append([PR_site[2], 'attP '+name_int[2],-1])

    # construction for 5 inputs
    elif nb_input==5:

        if list_gene_state[5]!=0:
            sequence.append([out_r[list_gene_state[5]-1], out_name[list_gene_state[5]-1],-1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([BF_site[4], 'attB '+name_int[4],1])
        sequence.append([BF_site[1], 'attB '+name_int[1],1])
        sequence.append([BR_site[2], name_int[2],-1])

        if list_gene_state[2]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[2]-1], out_name[list_gene_state[2]-1],1])
        
        sequence.append([sp7r, 'spacer 7', -1])
        sequence.append([L3S2P21r, 'L3S2P21',-1])
        sequence.append([sp6r, 'spacer 6', -1])
        
        if list_gene_state[1]!=0:
            sequence.append([out_r[list_gene_state[1]-1], out_name[list_gene_state[1]-1],-1])
            sequence.append([spacer_RBSr, 'spacer RBS', -1])
            sequence.append([B0034r, 'B0034', -1])

        sequence.append([BF_site[0], 'attB '+name_int[0],1])
        sequence.append([PF_site[1], 'attP '+name_int[1],1])
        sequence.append([J61048r, 'J61048', -1]) 
        sequence.append([PR_site[0], 'attP '+name_int[0],-1])

        sequence.append([BR_site[3], 'attB '+name_int[3],-1])
        sequence.append([P6r, 'P6',-1])
        sequence.append([PR_site[2], 'attP '+name_int[2],-1])
        
        if list_gene_state[3]!=0:
            sequence.append([B0034f, 'B0034', 1])
            sequence.append([spacer_RBSf, 'spacer RBS', 1])
            sequence.append([out_f[list_gene_state[3]-1], out_name[list_gene_state[3]-1],1])
        
        sequence.append([sp4f, 'spacer 4', 1])
        sequence.append([B0014f, 'B0014', 1])
        sequence.append([sp5f, 'spacer 5', 1])
        
        if list_gene_state[4]!=0:
            sequence.append([out_r[list_gene_state[4]-1], out_name[list_gene_state[4]-1],-1])

        sequence.append([PF_site[4], 'attP '+name_int[4],1])
        sequence.append([PR_site[3], 'attP '+name_int[3],-1])
    
    sequence.append([spNf, 'spacer N', 1])
        
##################################################### GENBANK FORMAT

    # initialization of the DNA sequence
    DNA_seq=''
    
    # loop in sequence list to generate the string corresponding to the DNA sequence of the device.
    for seq in sequence:
        DNA_seq+=seq[0]

    # creation of the formated DNA sequence for the genbank file
    seq_final=Seq(DNA_seq, IUPAC.unambiguous_dna)
    record = SeqRecord(seq_final,
                   id='NA', # random accession number
                   name='Strain'+nbr_strain,
                   description='Sequence of the computational device in the strain'+nbr_strain+' to implement the Boolean function required')
    
    # initialization of the variable len_seq 
    len_seq=0
    
    # loop in sequence list to generate the genbank feature of the sequence
    for feat in sequence:
        
        feature = SeqFeature(FeatureLocation(start=len_seq, end=len_seq+len(feat[0])), id=feat[1], type=feat[1], strand=feat[2])
        record.features.append(feature)
        len_seq+=len(feat[0]) 

    # Save as GenBank file
    output_file = open(path+'/'+name_directory+'_Strain'+nbr_strain+'.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')
    output_file.close()

   
"""TEST of the design DNA function"""    
#design_DNAsequence(3,(0,1,2),[4,3,2,1],'2','0000101010','../results/0','0')    
