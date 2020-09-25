#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 16:03:50 2020

@author: leila
3 entrega DNL: Formas normales

http://materias.df.uba.ar/dnla2020c1/files/2020/06/DNL_formasnormales.pdf

"""

#%%

import numpy as np
import os 
import sys
from sympy import symbols, lambdify
from sympy import solve, Matrix, Eq,Poly,degree
from sympy import diff
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#%%

path =  os.getcwd()

guardar_graficos = 1
solo_graphs_zoom = 0 #solo guardar los graficos cerca de los punto fijo 
guardar_output = 0
if guardar_graficos == 1:
    path_g = path + '/' + 'graficos'
    if not os.path.exists(path_g):
        print('Crear carpeta para guardar graficos')
        os.mkdir(path_g)
if guardar_output ==1:
    f = open("README", 'w')
    sys.stdout = f    
os.chdir(path)
    
#Definir parametros para graficos
tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%% Definir funciones del problema

x = symbols('x')
y = symbols('y')
a = symbols('a') 

xpunto = -x + x*x - y*y
ypunto = a*y - y*y*y + x*y

print('3 entrega DNL: Formas normales')
print('xpunto del problema:',xpunto)
print('ypunto del problema:',ypunto)

#%%

print('')
item = 'a)'
print(item)

"""
¿Cuánto debe valer a para que la forma normal a orden 2 pueda escribirse de la
siguiente manera?
"""

A = Matrix([[a,1],[0,-1]])

eig_vec = A.eigenvects()
v1,v2 = eig_vec[0][2][0],eig_vec[1][2][0] 
v1 = [v1[0],v1[1]]
v2 = [v2[0],v2[1]]
lambda1, lambda2 = A.eigenvals()

print(lambda1,lambda2)
print(v1,v2)

















