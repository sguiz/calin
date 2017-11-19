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
import pandas as pd

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'

    
"""
CREATE DICO SEQ

Extracte from a csv file, a dictionnary of DNA sequences.

"""    
    
def create_dico_seq(name_file):

    list_seq=pd.read_csv(name_file, sep=";",usecols=[0,1,2], names=['id','name','seq'])
    dico_seq={}
    
    for X in range(0,len(list_seq)):
        dico_seq[list_seq['id'][X]]=[list_seq['name'][X], list_seq['seq'][X]]
    
    return dico_seq

#print create_dico_seq('test.csv')
    
"""
GB CREATE

From the sequence list (composed of the DNA seq, name and orientation)
Generation of a gb file

"""

def gb_create(sequence, nb_strain, inp, path, name_directory):

    # initialization of the DNA sequence
    DNA_seq=''
    
    # loop in sequence list to generate the string corresponding to the DNA sequence of the device.
    for seq in sequence:
        DNA_seq+=seq[1]

    # creation of the formated DNA sequence for the genbank file
    seq_final=Seq(DNA_seq, IUPAC.unambiguous_dna)
    
    record = SeqRecord(seq_final,
                   id='NA', # random accession number
                   name='Strain'+nb_strain+'_'+inp,
                   description='Sequence of the computational device in the strain'+nb_strain+' to implement the Boolean function required')
       
    # initialization of the variable len_seq 
    len_seq=0
    
    # loop in sequence list to generate the genbank feature of the sequence
    for feat in sequence:
        
        feature = SeqFeature(FeatureLocation(start=len_seq, end=len_seq+len(feat[1])), id=feat[0], type=feat[0], strand=feat[2])
        record.features.append(feature)
        len_seq+=len(feat[1]) 

    # Save as GenBank file
    output_file = open(path+'/'+name_directory+'_Strain'+nb_strain+'.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')
    output_file.close()
    
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

def design_DNAsequence(inp, nbr_strain, output, path, name_directory, seq_file):
  
    ##############################################################################
    # SEQUENCES OF THE DESIGN
    
    # initialization of list sequence
    dico=create_dico_seq(seq_file)
    
    # it will be incremented with the DNA sequence and name of the feature that composed the full final DNA sequence.
    sequence=[]
    
    # append spacer sp0
    sequence.append([dico['sp0'][0],dico['sp0'][1], 1])
    
    # number of ZERO correspond to number of NOT elements: total variables (0) - number of ONE (1)
    ZERO=int(inp[0])-int(inp[1])
    
    # loop to generate all attB site for the NOT elements
    # in reverse order between the number of NOT elements to 0
    # Z corresponds to the indice of integrase to implement the NOT element
    for Z in reversed(range(0, ZERO)):
        # append to the design the attB site of NOT elements
        sequence.append([dico['BF'+str(Z+1)][0],dico['BF'+str(Z+1)][1], 1])
        
    # append to the design the promoter P7 in fwd direction        
    sequence.append([dico['Prom'][0],dico['Prom'][1],1])
    
    # loop to generate all attP site for the NOT elements
    # between 0 to the number of NOT elements
    # Z corresponds to the indice of integrase to implement the NOT element      
    for Z in range(0,ZERO):
        # append attP site
        sequence.append([dico['PF'+str(Z+1)][0],dico['PF'+str(Z+1)][1], 1])

    # loop to generate all IMPLY elements
    # O corresponds to the indice of integrase to implement the IMPLY element
    for O in range(ZERO, int(inp[0])):
        # add attB site+term+attP site
        sequence.append([dico['BF'+str(O+1)][0],dico['BF'+str(O+1)][1], 1])
        sequence.append([dico['T'+str(O+1)][0], dico['T'+str(O+1)][1], 1])
        sequence.append([dico['PF'+str(O+1)][0], dico['PF'+str(O+1)][1], 1])
        
    # append RiboJ, BCD2, sfGFP, spacer N at the end of the sequence
    sequence.append([dico['RiboJ'][0], dico['RiboJ'][1], 1])
    sequence.append([dico['RBS'][0], dico['RBS'][1], 1])
    sequence.append([dico['GOI1'][0], dico['GOI1'][1], 1])
    sequence.append([dico['Tend'][0], dico['Tend'][1], 1])
    sequence.append([dico['spN'][0], dico['spN'][1], 1])

    gb_create(sequence, nbr_strain, inp, path, name_directory)
   
"""TEST of the design DNA function"""
    
#design_DNAsequence('62','1','0000101010','../results/0','0', 'test.csv')    
