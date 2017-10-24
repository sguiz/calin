# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:07:29 2017

@author: sarahguiziou

To find biological implementation from a Boolean function
"""

import time
# module to simplify boolean function, Robert Dick algorithm
import qm2
from PIL import Image
import graphic_combinatorial_modules
import design_DNA_sequence
import random
import os

""" 

Test of Quine McCluskey algorithm

"""


def test_qm_simplify():
    ones = [0,5]
    dontcares = []
    print(qm2.qm(ones, dontcares))
    
#test_qm_simplify()



"""
Simplification of Boolean function

Give as input the on state with their number 
Ex: [0,3,5]

Give as outut:
Ex: ['1X', 'X1']
X independant to the variable, 1 the variable is equal to 1, 0 the variable is equal to 0

"""


def function_simplification(on_state):

    dontcares=[]
    pos_imp=qm2.qm(on_state, dontcares)
      
    return pos_imp
      
    
""" TEst function

"""
#print function_simplification([0])    

"""
Biological implementation

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
        
        n=0
        p=1
        int_variable=[]
        
        for variable in strain:

            if variable=='1':
                int_variable.append(nb_variable-n)
                n+=1
                
            if variable=='0':
                int_variable.append(p)
                p+=1
            
            if variable=='X':
                int_variable.append(0)
                
            
        int_var=''.join(str(e) for e in int_variable)
        
        list_int_var.append(int_var)
    
    
    return [nb_strain, list_component, list_int_var]


""" Test function

"""

#print biological_implementation(function_simplification([0]))    

""" FUNCTION which take as input the number of a component and as output show the image of this component
if good component number.

For now, working from 10 to 22

Number of the component need to have its second number smaller than it first one.
"""


def open_image_component(nb):
    
    if nb[1]>nb[0]:
        print "Error, not good number of component."
        
    else:
        full_path='/Users/sarahguiziou/Documents/Work/Project/Gates/Formalization of design/Combinatorial logic output gene/Soft/version 2.2.17/computation devices'
        title=full_path+'/'+str(nb)+'.jpg'
        img=Image.open(title)
        img.show()
    

# TEst of the function open image component
#open_image_component("11")

""" FINAL FUNCTION

Input: 
1 - OUTPUT, correspond to an integer which is correlate to a specific truth table.
2 - nb_input, correspond to the number of input needed for this function.

"""

def design(output, nb_input, name_directory, directory_number):
    
    # convert the output integer in a binary number
    #output=bin(output)[2:]
    name_output_file=name_directory+'/'+directory_number+'_output.txt'
    of=open(name_output_file, 'w')
    # if len of the output does not correspond to the number of input, add 0 to the number
    if len(output)<(2**(nb_input)):
        for w in range(2**(nb_input)-len(output)):
            output='0'+output
       
       
    # initialize the on_sate which will be a list of state equal to 1. The states are identifiy by an integer.        
    on_state=[]
    # index for input state
    x= 0

    # screen all output state (all bit), from the lowest bit to the highest one
    for c in reversed(output):
            
    # if output is equal to 1
        if c=='1': 
            # add X to the on_state list              
            on_state.append(x)
        
        x+=1
      
    # obtain the simplification from Cain McKluskey algorithm of the on_state
    # pos_imp correspond to a list of string. Number of string correspond to the nb of sub_function
    # Each string is composed of X for whatever, 1 for A or 0 for not(A) and the position of the term correspond to the variable.
    pos_imp=function_simplification(on_state)
    
    #print pos_imp
    
    # bioloigcal_implementation is to obtain the explicit biological implement
    # inputs are pos_imp (result from QM algo) and the number of inputs
    # outputs are the number of inputs, the number of strains needed, the list of computation devices
    # and the list of the connection between integrases and variables.
    info_imp=biological_implementation(pos_imp, nb_input)
    
    nb_strain=info_imp[0]
    list_component=info_imp[1]
    #print list_component
    list_int_var=info_imp[2]
    #print list_int_var
    
    of.write("To implement this Boolean function: "+output+' composed of '+str(nb_input)+' input(s) ')
    of.write(str(nb_strain)+" strains are needed.\n\n")
    
    strain=0
    
    # for each strain, loop in the list of component and list of integrase variable connections
    for COMP, X in zip(list_component, list_int_var):
        
        strain+=1
        
        #of.write("For strain "+str(strain)+", component CM"+COMP+" is needed.\n")
        of.write('For strain '+str(strain)+': \n')
        
        variable=1
        sub=[]
        
        # for the different variable, get the integrase number
        for integrase in X:
            
            if integrase!='0':
                of.write("Input "+str(variable)+" is connected to integrase "+integrase+'\n')
                
                sub.append('Input'+str(variable)+'-> Integrase'+integrase)
            
            
            variable+=1
        
        of.write('\n')
        # Correspond to the sub title of the graph which explicit variable integrase connections   
        subtitle=', '.join(e for e in sub)   
        # To represent the computation device with the information of the strain number and of var/int connections
        graphic_combinatorial_modules.computation_device(COMP, str(strain), subtitle, 'Strain'+str(strain), name_directory, directory_number)
        graphic_combinatorial_modules.legend(nb_input, name_directory, directory_number)
        # open image from the good path
        #full_path='/Users/sarahguiziou/Desktop/New soft'
        #title=full_path+'/'+COMP+'.png'
        #img=Image.open(title)
        #img.show()
        
    of.close()

""" Test function main
"""
#design('0010110000000000001111101110000100101111001100111010011000001011', 6)

def main(nb_input, output, random_number):
    
    nb_input=int(nb_input)
    
#    text_file=open(name_file, 'r')
#    f=list(text_file)
#    x=f[0]
#    y=x.split(' ')
#    output=y[0]
#    nb_input=int(y[1])
    
    directory_number=random_number
    
    while len(directory_number)<6:
        directory_number='0'+directory_number
    
    mypath='../results/'+directory_number
    os.makedirs(mypath)
    
    design(output, nb_input, mypath, directory_number)
    
#
main('4', '0000000110111100', 'G')
    
