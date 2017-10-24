# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 16:34:13 2016

@author: sarahguiziou
"""

import string

def design_DNAsequence(inp):
    
    ###############################################################################
    # SEQUENCE DEFINITION
    
    # To reverse complement sequences
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    
    
    # define attB sites
    BF_site=[]
    # Bxb1=1
    BF1='TCGGCCGGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCATCCGGGC'
    BF_site.append(BF1)
    # Tp901=2
    BF2='ATGCCAACACAATTAACATCTCAATCAAGGTAAATGCTTTTTGCTTTTTTTGC'
    BF_site.append(BF2)
    #Int7=3
    BF3='AGACGAGAAACGTTCCGTCCGTCTGGGTCAGTTGGGCAAAGTTGATGACCGGGTCGTCCGTT'
    BF_site.append(BF3)
    #int5=4
    BF4='GAGCGCCGGATCAGGGAGTGGACGGCCTGGGAGCGCTACACGCTGTGGCTGCGGTCGGTGC'
    BF_site.append(BF4)
    
    # attB reverse site
    BR_site=[]
    
    for a in range(len(BF_site)):
        BR_site.append(BF_site[a].translate(tab)[::-1])
    
    
    # define attP sites
    
    PF_site=[]
    # Bxb1=1
    PF1='TCGTGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCC'
    PF_site.append(PF1)
    # Tp901=2
    PF2='GCGAGTTTTTATTTCGTTTATTTCAATTAAGGTAACTAAAAAACTCCTTT'
    PF_site.append(PF2)
    #Int7=3
    PF3='GTGTTATAAACCTGTGTGAGAGTTAAGTTTACATGCCTAACCTTAACTTTTACGCAGGTTCAGCTT'
    PF_site.append(PF3)
    #Int5=4
    PF4='CCCTAATACGCAAGTCGATAACTCTCCTGGGAGCGTTGACAACTTGCGCACCCTGATCTG'
    PF_site.append(PF4)
    
    # attB reverse site
    PR_site=[]
    
    for a in range(len(PF_site)):
        PR_site.append(PF_site[a].translate(tab)[::-1])
        
    # define terminators
    
    ter=[]
    
    ECK120033737='GGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGG'
    T1=ECK120033737
    ter.append(T1)
    
    ECK120029600='TTGAGAAGAGAAAAGAAAACCGCCGATCCTGTCCACCGCATTACTGCAAGGTAGTGGACAAGACCGGCGGTCTTAAGTTTTTTGGCTGAA'
    T2=ECK120029600
    ter.append(T2)
    
    L3S2P21='CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC'
    T3=L3S2P21
    ter.append(T3)
    
    L3S3P21='CCAATTATTGAAGGCCTCCCTAACGGGGGGCCTTTTTTTGTTTCTGGTCTCCC'
    T4=L3S3P21
    ter.append(T4)
    
    
    # define spacers, promoters, BCD, RiboJ...
    # P2
    PF='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGA'
    PR=PF.translate(tab)[::-1]
    RiboJ='AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA'
    BCD2='GGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCTAAGGAGGTTTTCTA'
    sfGFP='atgtcaaaaggagaagaactttttacaggtgtagtacctatcttggttgaattggatggtgatgttaacggtcacaaattttctgtacgtggtgaaggtgaaggtgatgcaactaacggtaaattgacacttaaattcatttgtacaactggaaaacttcctgttccttggcctactcttgttacaacattgacatatggagtacaatgtttttcacgttatcctgatcatatgaaacgtcacgatttttttaaatctgctatgccagaaggttatgtacaagaacgtacaatttcatttaaagatgacggaacatataaaacacgtgctgaagtaaaattcgaaggtgacactcttgttaatcgtatcgaattgaaaggaatcgatttcaaagaagatggtaacattttgggacacaaacttgaatacaacttcaactctcataatgtttatatcacagctgacaaacaaaaaaacggtattaaagctaattttaaaattcgtcacaatgttgaagatggatctgttcaattggctgatcattatcaacaaaatacaccaatcggagacggaccagtattgcttccagataaccactacctttctactcaatcagttctttcaaaagatcctaacgaaaaacgtgaccatatggtacttcttgaatttgttacagcagcaggtatcactcacggtatggacgaactttataaataa'
    sp0='CTCGGATACCCTTACTCTGTTGAAAACGAATAGATAGGTT'
    sp20_1='TAGTTGCGTCTCAGGGACCC'
    sp20_2='TAAGTGGCAATCCCGCCTGA'
    sp20_3='AAACCCGTCGCAGTATCCCT'
    sp20_4='ACTCAGGTCTGCCGTAAGGG'
    
    # asymetric terminators
    ter=[]
    
    B0015='ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata'
    ter.append(B0015)
    
    J61048='ccggcttatcggtcagtttcacctgatttacgtaaaaacccgcttcggcgggtttttgcttttggaggggcagaaagatgaatgactgtccacgacgctatacccaaaagaaa'
    ter.append(J61048)
    
    ECK120015170='ACAATTTTCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGCTAAAA'
    ter.append(ECK120015170)
    
    ECK120010855='GTAACAACGGAAACCGGCCATTGCGCCGGTTTTTTTTGGCCT'
    ter.append(ECK120010855)
    
    ##############################################################################
    # SEQUENCES OF THE DESIGN
    
    ## output file, definition of directory and open file
    #output_file='/Users/sarahguiziou/Desktop/seq.txt'
    #f=open(output_file, 'w')
    
    design=''
    # if xor or nxor in the construct    
    if ('R' in inp) or ('N' in inp):
    
        # number of variable in the xor or nxor
        nb_xor=int(inp[1])
        # number of variable otherwise
        nb_other=int(inp[2])
        
        # for variable in the xor or nxor
        for S in range(0, nb_xor):
            # add to the design the attB site
            design+=BF_site[S]
        # if xor, add a reverse promoter
        if 'R' in inp:
            design+=PR
        # if nxor add a fwd promoter
        else:
            design+=PF
        #for variable in the xor or nxor add the attP site in reverse orientation in the reverse order of variables
        for S in reversed(range(0, nb_xor)):
            design+=PR_site[S]
        
        # number of ZERO in the second part of the construct
        # inp[3] correspond to the number of ONE in the construct
        ZERO=nb_other-int(inp[3])
        # variable which correspond to ZERO   
        for Z in range(nb_xor, ZERO+nb_xor):
            #add at the beginning of the design the attB site
            design=BF_site[Z]+design
        # variable which correspond to ZERO        
        for Z in range(nb_xor, ZERO+nb_xor):
            # add at the end of the design the attP site
            design+=PF_site[Z]
        # variable which correspond to ONE
        for O in range(ZERO+nb_xor, nb_xor+nb_other):
            # add attB site+term+attP site
            design+=BF_site[O]
            design+=ter[O]
            design+=PF_site[O]
        
    ## if no xor neither nxor   
    else:
        # number of ZERO correspond to number of total variables (0) - number of ONE (1)
        ZERO=int(inp[0])-int(inp[1])
        
        # variable which correspond to ZERO, in reverse order
        for Z in reversed(range(0, ZERO)):
            # append to the design the attB site
            design+=BF_site[Z]
        # append to the design the promoter in fwd direction        
        design+=PF
        # variable which correspond to ZERO        
        for Z in range(0,ZERO):
            # append attP site
            design+=PF_site[Z]
        # variable which correspond to ONE
        for O in range(ZERO, int(inp[0])):
            # add attB site+term+attP site
            design+=BF_site[O]
            design+=ter[O]
            design+=PF_site[O]
    # append output gene
    design+=RiboJ+BCD2+sfGFP
    #f.close()
    
    return design
    

"""TEST"""
    
#print design_DNAsequence('32')    
