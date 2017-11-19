#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri May  5 09:06:02 2017

@author: sarahguiziou
"""

import sys
import sequential_logic_design
import combinatorial_logic_design

if sys.argv[1]=='seq':
    sequential_logic_design.main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
else:
    combinatorial_logic_design.main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
