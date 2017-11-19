# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:07:29 2017

@author: sarahguiziou

To find biological implementation of a Boolean function
"""

# module to simplify boolean function, Robert Dick algorithm
import qm2
import graphic_combinatorial_modules as gcm
import design_DNA_sequence_comb_customized as dnsc
import os

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'


"""
FUNCTION SIMPLIFICATION

Give as input the ON state with their corresponding number. 
Ex: [0,1,5]

Using the qm function.

Give as outut:
Ex: ['X01', '00X'] that corresponds to not(B).C+not(A).not(B)
Function is wirtten as the disjonctive normal form, each string correspond to a different conjonction.
X: independant to the variable, 1: presence of the variable (A), 0: absence of the variable (not(A)).
In the order of the inputs: Input1(A) Input2(B) Input3(C).

"""


def function_simplification(on_state, off_state):

    dontcares=[]
    pos_imp=qm2.qm(ones=on_state, zeros=off_state ,dc=dontcares)
      
    return pos_imp
      
    
""" 
Test of function_simplification

"""
#print function_simplification([0,1],[2,3])    


"""
BIOLOGICAL IMPLEMENTATION

Inputs:
1:pos_imp: simplification of the boolean function in the disjonctive normal form in the following format:['1X', 'X1']
2:nb_input: number of inputs

Return:
1:nb_strain: the number of strains needed to impelement the boolean function.
2:list_component: the list of computational devices needed for the implementation of the boolean function
From the computational device in strain 1, to the computation device in strain N.
3:list_int_var: correspondance of input with integrase for each strains 
[list_strain1, list_strain2, ...] with list_strain1 a string


"""

def biological_implementation(pos_imp, nb_input):
        
    nb_strain=len(pos_imp)
    
    list_component=[]
    list_int_var=[]
    
    # loop in pos_imp
    for strain in pos_imp:

        one=strain.count('1')
        x=strain.count('X')
        
        nb_variable=nb_input-x
        
        # component correspond to the number of component
        # It is a string
        # The first number correspond to the number of variables
        # the second number correspond to the number of 1
        component=''.join(str(e) for e in [nb_variable, one])
        # ex: 10 not(A), 21 not(A).B
        # list of component
        list_component.append(component)
        
        zero=nb_variable-one
        n=0
        p=1
        int_variable=[]
        
        for variable in strain:

            if variable=='1':
                int_variable.append(str(zero+1+n))
                n+=1
                
            if variable=='0':
                int_variable.append(str(p))
                p+=1
            
            if variable=='X':
                int_variable.append(str(0))   
        
        list_int_var.append(int_variable)
        
    return [nb_strain, list_component, list_int_var]


""" 
Test function biological_implementation
"""

#print biological_implementation(function_simplification([0,1]), 2)    

""" 
FUNCTION DESIGN

Input: 
1: the truth table output: binary number that correlates to a specific truth table.
2: the number of inputs for this Boolean function (limit to 6 inputs due to the generation of the DNA sequence)
3: mypath: the path for the result directory
4: directory_name: the name of the result directory

Using graphic_computational_module, generation of the image file corresponding
to the biological design in each strains.
Generation of a text file resuming the number of strains needed to implement the boolean function
and the connection between inputs and integrases.
"""

def design(output, nb_input, path, directory_name, DNA_file):
    
    # convert the output integer in a binary number
    #output=bin(output)[2:]
    name_output_file=path+'/'+directory_name+'_output.txt'
    of=open(name_output_file, 'w')
    # if len of the output does not correspond to the number of input, add 0 to the number
    if len(output)<(2**(nb_input)):
        for w in range(2**(nb_input)-len(output)):
            output='0'+output
       
       
    # initialize the on_sate which will be a list of state equal to 1. The states are identifiy by an integer.        
    on_state=[]
    off_state=[]
    # index for input state
    x= 0

    # screen all output state (all bit), from the lowest bit to the highest one
    for c in reversed(output):
            
    # if output is equal to 1
        if c=='1': 
            # add X to the on_state list              
            on_state.append(x)
            
        else:
            off_state.append(x)
        
        x+=1
      
    # obtain the simplification from Cain McKluskey algorithm of the on_state
    # pos_imp correspond to a list of string. Number of string correspond to the nb of sub_function
    # Each string is composed of X for whatever, 1 for A or 0 for not(A) and the position of the term correspond to the variable.
    pos_imp=function_simplification(on_state, off_state)
    
    #print pos_imp
    
    # bioloigcal_implementation is to obtain the explicit biological implement
    # inputs are pos_imp (result from QM algo) and the number of inputs
    # outputs are the number of inputs, the number of strains needed, the list of computation devices
    # and the list of the connection between integrases and variables.
    [nb_strain, list_component, list_int_var]=biological_implementation(pos_imp, nb_input)
    
    print list_int_var
    
    of.write("This "+str(nb_input)+'-input Boolean function '+output+' is implemented using ' +str(nb_strain)+' strain(s).\n\n')
    of.write('DNA sequences are generated using <i>E. coli</i> promoters, terminators and GOI sequences.\n')
    
    name_int=['Bxb1', 'Tp901', 'Int5', 'Int7', 'Int3', 'Int4']
    
    texte_int1=''
    texte_int2=''
    
    for name in range(0, nb_input):
        if name==nb_input-1:
	    texte_int1+=(str(name+1))
	    texte_int2+=(name_int[name]+'.')
	else:
	    texte_int1+=(str(name+1)+', ')
	    texte_int2+=(name_int[name]+', ')
	    
    if nb_input==1:
        texte_int=('Integrase '+texte_int1+' corresponds to '+texte_int2)
    else:
        texte_int=('Integrases '+texte_int1+' correspond respectively to '+texte_int2) 
	
    of.write(texte_int+'\n\n')
    of.close()
    strain=0
    
    # for each strain, loop in the list of component and list of integrase variable connections
    for COMP, X in zip(list_component, list_int_var):
        
        strain+=1
        
        # To represent the computation device with the information of the strain number and of var/int connections
        gcm.computation_device(COMP, str(strain), X, 'Strain'+str(strain), path, directory_name)
        # Generation of the DNA sequence
        dnsc.design_DNAsequence(COMP, str(strain), output, path, directory_name, DNA_file)
        
    gcm.legend(nb_input, path, directory_name)

        


"""
TEST function design
"""

#design('00010011', 3, '../results/334', '334')


""" 
FUNCTION MAIN

Take as inputs:
1: the number of inputs (limt to 6 inputs, due to the generation of the DNA sequence)
2: the truth table output: binary number
3: the name of the directory: directory where the results will be stored

Run the design function with the corresponding inputs.

Biological design corresponding to the desired truth table are saved in the input directory.

"""

def main(nb_input, output, directory_name, type_DNA):
    
    nb_input=int(nb_input)
    
    while len(directory_name)<6:
        directory_name='0'+directory_name
    
    mypath='../results/'+directory_name
    
    DNA_file=type_DNA+'.csv'
    
    # to create the directory mypath
    os.makedirs(mypath)    
    
    design(output, nb_input, mypath, directory_name, DNA_file)


""" 
Test function

"""    
#main('3','11111110', 'test34', 'bsubtilis')
    
