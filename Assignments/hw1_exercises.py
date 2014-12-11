# -*- coding: utf-8 -*-
"""
Created on Tue Sep 09 18:33:02 2014

@author: Gary
"""

import numpy
import math

line = numpy.linspace(4, 23, 43)
print(line[5])

ones_array = numpy.ones((5,17))
zeros_array = numpy.zeros(ones_array.shape)

zeroShape = numpy.shape(zeros_array)
print(zeroShape)
#print(zeros_array.size())

p = 7
r = numpy.array([11.2, 4.7, 6.6])
p_over_r = p / r
sin_res = numpy.sin(p_over_r)
pow_res = numpy.power(sin_res, 3)
result = pow_res[1]
print(result)

dt_values = numpy.array([0.04, 0.02, 0.01])
u_values = numpy.array([1.600, 1.500, 1.475])

order_of_convergence = math.log((u_values[0] - u_values[1]) / (u_values[1] - u_values[2])) / math.log(2)
print(order_of_convergence)
