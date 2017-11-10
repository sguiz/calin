# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 12:13:18 2017

@author: sarahguiziou
"""
import math
import itertools
import graphic_sequential_modules as gsm
import graphic_combinatorial_modules as gcm
import design_DNA_sequence_seq as DNAseq
import design_DNA_sequence_comb as DNAcomb
import os
import copy
import random

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'

"""
CONSTRUCTION STATE MATRIX

Input: number of inputs

Output:
Generate the sequential state matrix for the defined number of inputs
The size of the matrix is equal to (factorial(N), N+1) and composed of
the state index. Each line corresponds to a different lineage, a specific
order of occurence of inputs, and each row corresponds to a number of
input present (from 0 input present to N input presents)

"""

def construction_state_matrix(nb_input):
    
    # Number of lineage equal to factorial of the number of inputs
    nb_lineage=math.factorial(nb_input)
    
    # generate a matrix of the size of the state matrix fill with zero.
    matrix=[[0 for X in range(0, nb_input+1)] for Y in range(0, nb_lineage)]
    
    # intialisation of the state indice
    indice=1
    
    # from 1 input used to N inputs, which correspond to the different column
    for nb_input_used in range(1, nb_input+1):
        
        # number of time that one number should be used, which depend
        # to the number of inputs used and the number of inputs.
        frequence=math.factorial(nb_input-nb_input_used)
        
        # the lineage indice, from 0 to nb_lineage
        lineage=0
        
        # while the lineage indice is different to the number of lineage
        while lineage!=nb_lineage:
            # for a number of time correspond to the frequence, we add 
            #the indice to the matrice at the specific place
            for X in range(0, frequence):
                
                matrix[lineage][nb_input_used]=indice
                
                #increment the lineage indice
                lineage+=1
                
            # increment the state indice
            indice+=1
    
    # return the matrice of interest       
    return matrix


"""
TEST Construction state matrix
"""
#print construction_state_matrix(3)

"""
LINEAGE TO ORDER OF INPUTS

Inputs:
1: number of inputs
2: lineage index

This function  permits to convert the number of the lineage
to a list composed of the inputs in the order of occurence.

Output:
Order of occurence of inputs: list of number in the corresponding order,
from the first input (left) to the last one (right).
"""

def lineage_to_order_inputs(nb_input, lineage):
    
    #list of inputs
    liste=[X for X in range(1, nb_input+1)]
    
    # permutation of the list of inputs to obtain all orders of inputs for all lineages.
    all_lineage_order=list(itertools.permutations(liste))
    
    # order of inputs obtains from the list with all permutations. 
    order_inputs=all_lineage_order[lineage]

    # return the order of inputs
    return order_inputs

"""
TEST Lineage to order inputs
ex
>>>print lineage_to_order_inputs(4, 5)
(1,4,3,2)

"""   
#print lineage_to_order_inputs(4, 5)    

"""
EXTRACTION OF CONSTRUCTIONS

This function permits from a number of inputs and an input sequential function
to obtain for each devices needed, the corresponding lineage and the GOI position
with a gene.
Input function is of the same size than the state matrix and is composed of the output state,
either 0 for no output, or any integer for any type of inputs. 

Input example: 
f=[[1,1,1,1], [1,1,0,0],[1,1,1,0], [1, 1, 0, 1], [1,0,0,0], [1,0,0,0]]

Output example: 
[[0, [1, 1, 1, 1]], [3, [0, 1, 0, 1]], [2, [0, 0, 1, 0]]]
The first number corresponds to the lineage number and the second one to 
which gene is present from the state 0 to N state...

"""

def extraction_of_constructions(nb_input, function):
    
    # initialization of imp, that will correspond to the different devices to implement
    # list with for each devices, the indice of the lineage and a list with the presence or absence of each genes
    imp=[]    
    
    # state matrix from function: construction_state_matrix
    matrix=construction_state_matrix(nb_input)
    
    # to determine the max id of output gene
    max_output=0
    
    # LOOP TO SIMPLIFY THE INPUT FUNCTION, to obtain the reduced set of lineage to implement
    # Scan the input function from the right to the left.
    # If state is not equal to Zero, we store the lineage for implementation
    # and simplify the input function to not implement twice the same states.
    
    # priority: the column in the right
    for column in reversed(range(0, nb_input+1)):
        
        # index of the lineage
        id_lineage=0
        
        # for the lineage in the function (which correspond to the line)
        for lineage in function:
            
            # if the value of the specific lineage and specific column.
            # is equal to 1.     
            if lineage[column]!=0:
                
                # construction with the specific lineage.
                imp.append([id_lineage, list(lineage)])
                
                # indice of the different value in the lineage.
                id_value=column+1
                
                # loop into the different column of the selected lineage
                while id_value!=0:
                    id_value+=-1
                    
                    # value corresponds to the value of the state of 
                    #the corresponded lineage
                    value=lineage[id_value]
                    
                    if value>max_output:
                        max_output=copy.deepcopy(value)
                    
                    # if the value is equal to 1
                    #if value==1:
                    if value!=0:
                        
                        # take the id corresponding to this state in the matrix
                        id_correspondance=matrix[id_lineage][id_value]
                        
                        # loop into all lineages for this column
                        # SIMPLIFICATION OF THE MATRIX
                        for X in range(0, len(matrix)):

                            # if the id of this cell corresponds to the previous one
                            # we set the value to 0
                            if matrix[X][id_value]==id_correspondance:
                                function[X][id_value]=0
                
                # set the lineage to 0
                function[id_lineage]=[0 for X in range(0,nb_input+1)]
            
            # increment the lineage
            id_lineage+=1
    
    return imp, max_output
    
""" 
TEST extraction of constructions
ex:
>>> print extraction_of_constructions(3,[[1,1,1,1], [1,1,0,0],[1,0,0,0], [1, 0, 2, 0], [1,0,0,0], [1,0,0,0]])
[[0, [1,1,1,1]], [3,[0,0,2,0]]]
""" 

#f=[[1,1,1,1], [1,1,0,0],[1,0,0,0], [1,0,2,0], [1,0,4,0], [1,0,8,0]]
#print extraction_of_constructions(3, f)

"""
STT FROM OUTPUT STRING
Generation of the sequential truth table (input function) from a string of outputs
for the different sequential states. 

Inputs:
1- string of outputs: from the last state to the first one
ex: '12034' for 2 inputs, output 1 for state 5, output 4 for state
2- number of inputs

Output:
Sequential truth table matrix
"""

def STT_from_output_string(bin_output, nb_input):

    matrix_eq=construction_state_matrix(nb_input)

    nb_lineage=math.factorial(nb_input)
    
    # generate a matrix of the size of the equivalent matrix fill with zero.
    matrix=[[0 for X in range(0, nb_input+1)] for Y in range(0, nb_lineage)]

    n=0
    
    for state in reversed(bin_output):
        st=int(state)
        
        if st!=0:
            coor_lineage=0
            coor_node=0
            
            for lineage in matrix_eq:
                coor_node=0
                for node in lineage:
                    if node==n:
                        matrix[coor_lineage][coor_node]=st
                    coor_node+=1
                coor_lineage+=1
        n+=1
    

    return matrix

""" 
TEST STT from output string

>>>STT_from_output_string('1231141110000110', 3)
[[0, 1, 0, 4], [0, 1, 0, 1], [0, 1, 0, 1], [0, 1, 1, 3], [0, 0, 1, 2], [0, 0, 1, 1]]
"""
        
#print STT_from_output_string('1231141110000110', 3)

"""
MAIN function

Input: 
1: number of inputs
2: input string corrsponding to the sequential program to implement
3: directory

Generation of the output files in the directory; a text file, png files with 
the biological construtions and gb file with the DNA sequences

Use the functions:
STT_from_output_string
extraction_of_constructions
graphic_sequential_modules

"""
 
def main(nb_input, f, directory_name):

    # number of input in integer
    nb_input=int(nb_input)

    # generate the STT using the following function
    function=STT_from_output_string(f, nb_input)
    
    # copy of the STT
    STT=copy.deepcopy(function)
    
    # check the len of the directory name and complete with zero if lower than 6
    while len(directory_name)<6:
        directory_name='0'+directory_name
    
    # directory path to store the output files    
    mypath='../results/'+directory_name
    os.makedirs(mypath)
    
    name_output_file=mypath+'/'+directory_name+'_output.txt'
    of=open(name_output_file, 'w')

    # From the STT, extraction of the caracteristic for the required devices
    imp, max_output=extraction_of_constructions(nb_input, function)
    
    strain=0
    
    for construction in imp:
        strain+=1
        
        list_integrase=lineage_to_order_inputs(nb_input, construction[0])

        title='Strain'+str(strain)
        gsm.computational_device(nb_input, list_integrase, construction[1], title, mypath, directory_name)
        gsm.legend(nb_input, mypath, directory_name)
        
        if max_output<5:
            DNAseq.design_DNAsequence(nb_input, list_integrase, construction[1], str(strain), '', mypath, directory_name)
   
    if nb_input==5:
        if f[325]=='1':
            strain+=1
            gcm.computation_device('50', str(strain), ['1','2','3','4','5'], 'Strain'+str(strain), mypath, directory_name)
            
            if max_output<5:
                DNAcomb.design_DNAsequence('50', str(strain), '', mypath, directory_name)
            
    #of.write('For the following sequential program:\n'+str(STT)+'\n\n')
    #of.write(str(strain)+' strains are needed.\n\n')
    #of.close()
    
    of.write("This "+str(nb_input)+'-input history-dependent function, \n'+str(STT)+', \nis implemented with '+str(strain)+' strain(s).\n\n')
    
    if max_output<5:
        of.write('The DNA sequence is generated using <i>E. coli</i> promoters, terminators and GOI sequences.\n')
        name_int=['Bxb1', 'Tp901', 'Int5', 'Int7', 'Int3', 'Int4']
        name_output=['GFP', 'RFP', 'BFP', 'LacZ']
        
        
        texte_int=''
        texte_int1=''
        texte_int2=''
        texte_out1=''
        texte_out2=''
    
        for name in range(0, nb_input):
            if name==nb_input-1:
	        texte_int1+=(str(name+1))
	        texte_int2+=(name_int[name]+'.')
	    else:
	        texte_int1+=(str(name+1)+', ')
	        texte_int2+=(name_int[name]+', ')
	        
        if nb_input==1:
            texte_int=('The integrase '+texte_int1+' corresponds to '+texte_int2)
        else:
            texte_int=('The integrases '+texte_int1+' correspond respectively to '+texte_int2) 
            
        for name in range(0, max_output):
            if name==max_output-1:
	        texte_out1+=(str(name+1))
	        texte_out2+=(name_output[name]+' ')
	    else:
	        texte_out1+=(str(name+1)+', ')
	        texte_out2+=(name_output[name]+', ')
	        
        if nb_input==1:
            texte_out=('For the output '+texte_out1+', the '+texte_out2+'gene is used')
        else:
            texte_out=('For the outputs '+texte_out1+', the '+texte_out2+'genes are used.') 
        
        of.write(texte_int+'\n')    
        of.write(texte_out+'\n\n')
            
    else:
        of.write('Generation of DNA sequence file of computational device is not possible using CALIN website as IDs of output gene are biggest than 4. Please refer to <a href="https://github.com/sguiz/calin" target="_blank" class="nav">github</a> for custom DNA generation.\n\n')



    of.close()
    strain=0
    
    return

"""
TEST MAIN
"""
#f='32001'
#name=random.randint(0,10000)
#main('2', f, 'example_seq')    

    
