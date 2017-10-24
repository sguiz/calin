# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 12:13:18 2017

@author: sarahguiziou
"""
import math
import itertools
import graphic_sequential_modules
from PIL import Image
import os
import random
import sys

def construction_matrix_equivalence(nb_input):
    
    nb_lineage=math.factorial(nb_input)
    
    # generate a matrix of the size of the equivalent matrix fill with zero.
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
                
                #incremente the lineage indice
                lineage+=1
                
            # increment the state indice
            indice+=1
    
    # return the matrice of interest       
    return matrix


#X=construction_matrix_equivalence(2)
#print X
#print X[1][0]
#print X[0][1]


"""
This function  permits to convert the number of the lineage
to a list composed of the inputs in the order of apperence.
"""

def nblineage_to_orderintegrase(nb_input, nb):
    
    liste=[X for X in range(1, nb_input+1)]
    
    all_lineage_function=list(itertools.permutations(liste))
    
    correspondance=all_lineage_function[nb]

    return correspondance

###################TEST########################    
#print nblineage_to_orderintegrase(4, 5)    

"""
This function permits from a number of inputs and a sequential function
which correspond to a matrice to obtain the lineage needed and the genes that
have to be implemented

Input example: 
f=[[1,1,1,1], [1,1,0,0],[1,1,1,0], [1, 1, 0, 1], [1,0,0,0], [1,0,0,0]]

Output example: 
[[0, [1, 1, 1, 1]], [3, [0, 1, 0, 1]], [2, [0, 0, 1, 0]]]
The first number corresponds to the lineage number and the second one to 
which gene is present from the state 0 to N state...

"""

def extraction_of_constructions(nb_input, function):
    
    imp=[]    
    
    matrix=construction_matrix_equivalence(nb_input)
    
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
                    
                    # value correspond to the value of the state of 
                    #the corresponded lineage
                    value=lineage[id_value]
                    
                    # if the value is equal to 1
                    #if value==1:
                    if value!=0:
                        
                        # take the id corresponding to this state in the matrix
                        id_correspondance=matrix[id_lineage][id_value]
                        
                        # loop into all lineages for this column
                        for X in range(0, len(matrix)):

                            # if the id of this cell correspond to the previous one
                            # we set the value to 0
                            if matrix[X][id_value]==id_correspondance:
                                function[X][id_value]=0
                
                # set the lineage to 0
                function[id_lineage]=[0,0,0,0]
            
            # increment the lineage
            id_lineage+=1

    return imp
    
###############TEST##################    

#f=[[1,1,1,1], [1,1,0,0],[1,0,0,0], [1, 0, 2, 0], [1,0,0,0], [1,0,0,0]]
#print extraction_of_constructions(3, f)


def switch_STT_bin_matrix(bin_output, nb_input):

    matrix_eq=construction_matrix_equivalence(nb_input)

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

######TEST
        
#print switch_STT_bin_matrix('1231141110000110', 3)


 
def main(nb_input, f, random_number):

    nb_input=int(nb_input)

    function=switch_STT_bin_matrix(f, nb_input)
    directory_number=random_number
    
    while len(directory_number)<6:
        directory_number='0'+directory_number
        
    mypath='../results/'+directory_number
    os.makedirs(mypath)
    
    name_output_file=mypath+'/'+directory_number+'_output.txt'
    of=open(name_output_file, 'w')

    imp=extraction_of_constructions(nb_input, function)

    
    of.write('For the following function:\n'+f+'\n')
    of.write(str(len(imp))+' strains are needed')

    nb=0
    
    for construction in imp:
        nb+=1
        
        list_integrase=nblineage_to_orderintegrase(nb_input, construction[0])

        title='Strain'+str(nb)
        graphic_sequential_modules.main(nb_input, list_integrase, construction[1], title, mypath, directory_number)
        
        # open image from the good path
#        full_path='/Users/sarahguiziou/Desktop/New soft'
#        title=full_path+'/'+title+'.png'
#        img=Image.open(title)
#        img.show()

    return

#f='02210'    
#main('2', f, '33')    

    
