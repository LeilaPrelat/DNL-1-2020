#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  14 15:00:05 2020

@author: leila

2 entrega DNL: Variedad Central

Este sistema presenta 2 bifurcaciones

http://materias.df.uba.ar/dnla2020c1/files/2020/06/DNL_variedadcentral.pdf
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

print('2 entrega DNL: Variedad Central')
print('xpunto del problema:',xpunto)
print('ypunto del problema:',ypunto)
print('')

#%%

def puntos_fijos(xpunto,ypunto,shh): #shh: printear o no printear
    try:
        nulclina1 = solve(xpunto, y)
        nulclina1[0]
    except IndexError:
        #el len(nulclina)=0
        nulclina1 = solve(xpunto, x)
    
    try:
        nulclina2 = solve(ypunto, y) 
        nulclina2[0]
    except IndexError:
        #el len(nulclina)=0
        nulclina2 =  solve(ypunto, x)  
    
    fnul1_list = []
    for fnul1 in nulclina1:
        fnul1_list.append(fnul1)
        if shh!=1:
            print('nulclina 1: xpunto = 0 sii y =',fnul1)
    
    extra1 = []
    try:
        float(xpunto.subs(x,0))
        if float(xpunto.subs(x,0))==0:
            extra1.append(x)
    except TypeError:
        pass
    
    fnul2_list = []
    for fnul2 in nulclina2:
        fnul2_list.append(fnul2)
        if shh!=1:
            print('nulclina 2: ypunto = 0 sii y =',fnul2)
    
    extra2 = []
    try:
        float(ypunto.subs(x,0))
        if float(ypunto.subs(x,0))==0:
            extra2.append(x)
    except TypeError:
        pass
        
    xpf_list = []
    ypf_list = []

    eq1 = Eq(xpunto,0)
    eq2 = Eq(ypunto,0)
    result = solve([eq1,eq2],(x,y))
    for res in result:
        xpf,ypf = res
        xpf_list.append(xpf)
        ypf_list.append(ypf)

    return [xpf_list,ypf_list,fnul1_list,fnul2_list,extra1,extra2]

#%%

def clasif_pf(xpunto,ypunto,a_ev,shh): #shh: printear o no printear
    a_ev = float(a_ev)
    if a_ev==0: #alguno es cero
        xpunto = xpunto.subs(a,0)
        ypunto = ypunto.subs(a,0)
        [listxpf,listypf,fnul1_list,fnul2_list,extra1,extra2]=puntos_fijos(xpunto,ypunto,1)  
    else:
        [listxpf,listypf,fnul1_list,fnul2_list,extra1,extra2]=puntos_fijos(xpunto,ypunto,1)  
    
    listxpf_ev = []
    listypf_ev = []
    k = 0
    for xpf in listxpf:
        xpf = xpf.subs(a,a_ev)
        listypf[k] = listypf[k].subs(a,a_ev)
        if complex(xpf).imag==0 and complex(listypf[k]).imag==0: #filtrar los valores complejos
            listxpf_ev.append(float(xpf))
            listypf_ev.append(float(listypf[k]))
        k = k + 1
        
    grad_xpunto = [diff(xpunto, x),diff(xpunto, y)]
    grad_ypunto = [diff(ypunto, x),diff(ypunto, y)]
    DF = Matrix([grad_xpunto,grad_ypunto])
    lambda1, lambda2 = DF.eigenvals()
    tau = DF[0] + DF[3]
    delta = DF[0]*DF[3] - DF[1]*DF[2]
    tau = lambdify([x,y,a], tau)
    delta = lambdify([x,y,a], delta)

    def eval3(x,y,a,expr):
        return expr(x,y,a)

    list_lmbd1_ev = []
    list_lmbd2_ev = []    
    list_tau_ev = []
    list_delta_ev = []
    types_pf = []
    if len(listxpf_ev)==0:
        print('No hay puntos fijos para los parametros elegidos')
    
    #valores no repetidos
    list_pf_tot = []
    for k in range(len(listxpf_ev)):
        xpf, ypf = listxpf_ev[k],listypf_ev[k]
        list_pf_tot.append([xpf,ypf])
    
    list_pf_tot2 = np.unique(np.array(list_pf_tot), axis=0)
    listxpf_ev = []
    listypf_ev = []
    for k in range(len(list_pf_tot2)):
        xpf, ypf = list_pf_tot2[k]
        listxpf_ev.append(xpf)
        listypf_ev.append(ypf)
    
    for j in range(len(listxpf_ev)):
        xpf = listxpf_ev[j]
        ypf = listypf_ev[j]
        delta_ev = eval3(xpf,ypf,a_ev,delta)
        tau_ev = eval3(xpf,ypf,a_ev,tau)
        lambda1, lambda2 = DF.eigenvals()
        lmbd1 = lambda1.subs(x,xpf).subs(y,ypf)
        lmbd2 = lambda2.subs(x,xpf).subs(y,ypf)
        lmbd1_ev = complex(lmbd1.subs(a,a_ev))
        lmbd2_ev = complex(lmbd2.subs(a,a_ev))
        deg = tau_ev**2-4*delta_ev
        list_lmbd1_ev.append(lmbd1_ev)
        list_lmbd2_ev.append(lmbd2_ev)
        list_delta_ev.append(delta_ev)
        list_tau_ev.append(tau_ev)
        if delta_ev<0:
            type_pf = 'Punto saddle'
            types_pf.append(type_pf)
            if shh!=1:
                print(type_pf,xpf,ypf)
        elif delta_ev>0:
            if tau_ev <0:
                if deg>0:
                    type_pf = 'Nodo estable'
                    types_pf.append(type_pf)
                    if shh!=1:
                        print(type_pf,xpf,ypf)            
                elif deg<0:
                    type_pf = 'Espiral estable'
                    types_pf.append(type_pf)
                    if shh!=1:
                        print(type_pf,xpf,ypf)
                if deg ==0:
                    type_pf = 'Nodo estrella'
                    types_pf.append(type_pf)                   
                    if shh!=1:
                        print(type_pf,xpf,ypf)
            elif tau_ev>0:
                 if deg>0:
                    type_pf = 'Nodo inestable'
                    types_pf.append(type_pf)     
                    if shh!=1:
                        print(type_pf,xpf,ypf)            
                 elif deg<0:
                    type_pf = 'Espiral inestable'
                    types_pf.append(type_pf)   
                    if shh!=1:
                        print(type_pf,xpf,ypf)
                 if deg ==0:
                    type_pf = 'Nodo estrella'
                    types_pf.append(type_pf)                   
                    if shh!=1:
                        print(type_pf,xpf,ypf)
            elif tau_ev==0:
                type_pf = 'Hopf'
                types_pf.append(type_pf) 
                if shh!=1:
                    print(type_pf,xpf,ypf)
        elif delta_ev==0:
            type_pf = 'Non-lineal'
            types_pf.append(type_pf) 
            if shh!=1:
                print(type_pf,xpf,ypf)
    return [listxpf_ev,listypf_ev,list_delta_ev,list_tau_ev,list_lmbd1_ev,list_lmbd2_ev,types_pf]

#%%

def retrato_fases(xpunto,ypunto,a_ev,item):
    """
    retrato_fases devuelve un grafico del retrato de fases
    parametros: (solo para dimension 2D y para 2 parametros no nulos)
        -xpunto,ypunto del problema
        -a_ev es el valor del parametro a
        -item para guardar el grafico con el nombre del item
        -minx,maxx = -3,3 por default
        (se cambian automaticamente para cubrir los puntos fijos)
    """
    a_ev = float(a_ev)
    if a_ev==0:
        xpunto = xpunto.subs(a,0)
        ypunto = ypunto.subs(a,0)
        [listxpf,listypf,fnul1_list,fnul2_list,extra1,extra2]= puntos_fijos(xpunto,ypunto,1)     
    else:
        [listxpf,listypf,fnul1_list,fnul2_list,extra1,extra2]= puntos_fijos(xpunto,ypunto,1)  
        
    def eval2(x,a,expr):
        return expr(x,a)
    def eval3(x,y,a,expr):
        return expr(x,y,a)
    
    [listxpf_ev,listypf_ev,delta_ev,tau_ev,lmbd1_ev,lmbd2_ev,types_pf] = clasif_pf(xpunto,ypunto,a_ev,0)
    
    xdot = lambdify([x,y,a], xpunto)
    ydot = lambdify([x,y,a], ypunto)
                
    n = int(1e3)
    minx,maxx = -3,3
    if len(listxpf_ev)>0: #vamos a hacer que el eje x siempre cubra los puntos fijos
        if minx > np.min(listxpf_ev):
            minx = np.min(listxpf_ev)-1
        if maxx < np.max(listxpf_ev):
            maxx = np.max(listxpf_ev)+1
    
    info = '$\dot{x}$=' + str(xpunto) + ', $\dot{y}$=' + str(ypunto)
    title = '%s Retrato de fases con a = %.2f'%(item,a_ev)
    print('Graficar retrato de fase') 
    plt.figure(figsize=tamfig)
    plt.title(info+'\n'+title,fontsize=tamtitle)
    ejex = np.linspace(minx,maxx,n)     
    maxisy = []
    minisy = []
    for fnul1s in fnul1_list:
        ejey1 = []
        ejex1 = []
        fnul1 = lambdify([x,a], fnul1s)
        for valuex in ejex:
            try:
                valuey1 = eval2(valuex,a_ev,fnul1)
                valuey1 = complex(valuey1)
                if valuey1.imag==0:
                    ejey1.append(valuey1.real)    
                    ejex1.append(valuex)
            except RuntimeWarning as error:
                print(error)
        plt.plot(ejex1,ejey1,'o',ms=6,label='$\dot{x}$=0: y =' + str(fnul1s))
        cleaned_ejey1 = [x for x in ejey1 if str(x) != 'nan']
        maxisy.append(np.max(cleaned_ejey1))
        minisy.append(np.min(cleaned_ejey1))

    for fnul1s in extra1: 
        plt.plot(np.ones(n)*0,ejey1,'o',ms=6,label='$\dot{x}$=0: x = 0')        
    for fnul2s in fnul2_list:
        ejey2 = []
        ejex2 = []
        fnul2 = lambdify([x,a], fnul2s)
        for valuex in ejex:
            try:
                valuey2 = eval2(valuex,a_ev,fnul2)
                valuey2 = complex(valuey2)
                if valuey2.imag==0:
                    ejey2.append(valuey2.real)    
                    ejex2.append(valuex)
            except RuntimeWarning as error:
                print(error)
        plt.plot(ejex2,ejey2,'*',ms=6,label= '$\dot{y}$=0: y =' + str(fnul2s))
        cleaned_ejey2 = [x for x in ejey2 if str(x) != 'nan']
        try: 
            maxisy.append(np.max(cleaned_ejey2))
            minisy.append(np.min(cleaned_ejey2))
        except ValueError:
            print('empty list:',len(cleaned_ejey2))
    for fnul2s in extra2: 
        plt.plot(np.ones(n)*0,ejey2,'*',ms=6,label='$\dot{y}$=0: x = 0')   
    
    if len(listxpf_ev)>0:
        plt.plot(listxpf_ev,listypf_ev,'k.',ms=12,label = 'Puntos fijos')
    else:
        plt.plot([],[],'w.',label = 'No hay puntos fijos')
    
    narrows = 60
    miny,maxy = np.min(minisy),np.max(maxisy)
    miny,maxy = float(miny),float(maxy)
    xs_list = np.linspace(minx, maxx, narrows)
    ys_list = np.linspace(miny, maxy, narrows)
    X,Y = np.meshgrid(xs_list,ys_list)
    U = eval3(X,Y,a_ev,xdot)
    V = eval3(X,Y,a_ev,ydot)
    plt.streamplot(X, Y, U, V,color='green') 
    if (len(fnul1_list)+len(fnul2_list))>=3:
        plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    else:
        plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.ylabel('y',fontsize=tamletra)
    plt.xlabel('x',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    if solo_graphs_zoom==0:
        if guardar_graficos==1:
            os.chdir(path_g)
            plt.savefig('item_%s' %(item))
            os.chdir(path)
    else:
        plt.close()
     
    # uniq_listypf_ev = set(listypf_ev) #valores no repetidos
    # uniq_listxpf_ev = set(listxpf_ev)
    if 0<len(listypf_ev): #flechas cerca de los puntos fijos
        plt.figure(figsize=tamfig)
        plt.title(info+'\n'+title,fontsize=tamtitle)
        plt.ylabel('y',fontsize=tamletra)
        plt.xlabel('x',fontsize=tamletra)
        min_ypf = np.min(listypf_ev)
        max_ypf = np.max(listypf_ev)
        deltay = np.abs(max_ypf-min_ypf)*0.05
        deltax = np.abs(maxx-minx)*0.05
        for ind in range(len(listypf_ev)):    
            ypf = listypf_ev[ind]
            xpf = listxpf_ev[ind]
           
            # plt.ylim([min_ypf-2,max_ypf+2])
            # plt.xlim([minx,maxx])
            xpf2,ypf2 = np.abs(xpf),np.abs(ypf)
            plt.text(xpf - deltax,ypf - deltay,types_pf[ind],fontsize=10)
            info1 = '$\Delta$ = %.2f,\u03C4 = %.2f' %(delta_ev[ind],tau_ev[ind])
            if lmbd1_ev[ind].imag==0 and lmbd2_ev[ind].imag==0:
                lmbd1,lmbd2 = lmbd1_ev[ind].real,lmbd2_ev[ind].real
                info2 = '$\lambda_1$ = %.2f,$\lambda_2$ = %.2f' %(lmbd1,lmbd2)
            else:
                lmbd11,lmbd12 = lmbd1_ev[ind].real,lmbd1_ev[ind].imag
                lmbd21,lmbd22 = lmbd2_ev[ind].real,lmbd2_ev[ind].imag
                if lmbd12>=0:
                    info_a = '$\lambda_1$ = %.2f + i%.2f' %(lmbd11,lmbd12)               
                else:
                    info_a = '$\lambda_1$ = %.2f - i%.2f' %(lmbd11,-lmbd12)               
                if lmbd22>=0:
                    info_b = '$\lambda_2$ = %.2f + i%.2f' %(lmbd21,lmbd22)
                else:
                    info_b = '$\lambda_2$ = %.2f - i%.2f' %(lmbd21,-lmbd22)
                info2 = info_a+','+info_b
            plt.plot(xpf,ypf,'.',ms=18,label = info1+','+info2)
            if (len(fnul1_list)+len(fnul2_list))>=3:
                plt.legend(loc='upper left',markerscale=1,fontsize=int(tamlegend*0.6))
            else:
                plt.legend(loc='upper left',markerscale=2,fontsize=tamlegend)
        
        nzoom = 40     
        epsig = 0.5
        ys_list = np.linspace(min_ypf-epsig,max_ypf+epsig,nzoom)
        xs_list = np.linspace(minx-epsig,maxx+epsig,nzoom)
        X,Y = np.meshgrid(xs_list,ys_list)
        U = eval3(X,Y,a_ev,xdot)
        V = eval3(X,Y,a_ev,ydot)   
        plt.streamplot(X, Y, U, V, density=[2, 1],color='green')         
        vert = np.linspace(np.min(ys_list),np.max(ys_list),10) 
        hori = np.linspace(np.min(xs_list),np.max(xs_list),10) 
        plt.plot(np.ones(10)*0,vert,'-k')
        plt.plot(hori,np.ones(10)*0,'-k')     
        plt.tick_params(labelsize = tamnum)    
        
        if guardar_graficos==1:
            os.chdir(path_g)
            print('Guardar grafico')
            plt.savefig('item_%s_zoom' %(item))
            os.chdir(path)
            
    return [listxpf_ev,listypf_ev]
    
#%% 

print('')
item = 'a)'
print(item)

"""
Calcule los puntos fijos (si no encuentra una expresión explícita, muestre su existencia
gráficamente)
"""

list_xpf = puntos_fijos(xpunto,ypunto,1)[0]
list_ypf = puntos_fijos(xpunto,ypunto,1)[1]
print('%s Puntos fijos'%(item))
for j in range(len(list_xpf)):
    xpf,ypf = list_xpf[j],list_ypf[j]
    print('x*,y*:', xpf,',',ypf)

barrido_a = [-1.5,-1.05,-1,-0.95,-0.5,-0.2,0,0.2,0.5,1]
j = 0
for value_a in barrido_a:
    print('Retrato de fases con a =%.2f'%(value_a))
    ittem = item + '_%i' %(j+1)
    retrato_fases(xpunto,ypunto,value_a,ittem)
    j = j + 1
    print('')

#%%

print('')
item = 'b)'
print(item)

"""
Muestre que uno de los autovalores del Jacobiano del sistema se anula en cada
bifurcación. Nombre cada bifurcación e identifique el punto fijo que está bifurcando en cada
caso.

bifurcaciones: a = 0, a = -1
"""

barrido_a = [-1.5,-1,-0.5,0,0.5,1]
a_no_hiper = []
y_pf_no_hiper = []
x_pf_no_hiper = []
for value_a in barrido_a:     
    [xpf_ev,ypf_ev,delta_ev,tau_ev,lmbd1_ev,lmbd2_ev,types_pf] = clasif_pf(xpunto,ypunto,value_a,1)
    for j in range(len(types_pf)):
        if types_pf[j]=='Non-lineal':
            print('a =', value_a ,'contine un punto no hiperbolico')
            a_no_hiper.append(value_a)
            y_pf_no_hiper.append(ypf_ev[j])
            x_pf_no_hiper.append(xpf_ev[j])
       
for j in range(len(a_no_hiper)):
    grad_xpunto = [diff(xpunto, x),diff(xpunto, y)]
    grad_ypunto = [diff(ypunto, x),diff(ypunto, y)]
    DF = Matrix([grad_xpunto,grad_ypunto])
    a_ev = a_no_hiper[j]
    xpf = x_pf_no_hiper[j]
    ypf = y_pf_no_hiper[j]
    DF = DF.subs(a,a_ev).subs(x,xpf).subs(y,ypf)
    eig_vec = DF.eigenvects()
    v1,v2 = eig_vec[1][2][0],eig_vec[0][2][0]
    v11 = v1[0],v1[1]
    v22 = v2[0],v2[1]
    print(v11)
    print(v22)
    
    lambda1, lambda2 = DF.eigenvals()
    print(a_ev,':',xpf,ypf)
    print(lambda1,lambda2)
    print('')
    
#cambio de base ---> llevar v1 al ejex y v2 al ejey (rotacion)

v1 = [1,0]
v2 = [0,1]

#%%


def K_var_no_central(list_var_central,orden):
    """
    devuelve la expansion polinomica de K de la 
    variable no central (pol mayor a orden 2
                         y menor a la var orden)
    """
    n = len(list_var_central)
    coef = symbols('a1:%d' %(int(orden*n+2)), commutative=True)
    coef = np.array(coef)
    poly = []
    a = 0
    for k in range(len(list_var_central)-1):    
        for j in range(0,orden+1):
            for m in range(0,orden+1):
                if 2<=j+m<=orden: #queremos que sea mayor a orden2
                    coeff = coef[a]
                    term = (list_var_central[k]**m)*(list_var_central[k+1]**j)
                    poly.append(coeff*term)
                    a = a + 1
    poly2 = 0
    for termino in poly: 
        poly2 =  termino + poly2
    return poly2,coef

#%%

def filtro_f(f,list_variables,order):
    """
    Parameters
    ----------
    f : function
    list_variables: ingresar variables de f (symbols)
    order : numero entero

    Returns
    devuelve los terminos de la funcion
    que tienen grado <= al order
    (considerando todas las variables de la funcion)

    """
    # list_variables = [x,y,a]
    f_2 = f.expand()
    f_3 = f_2.as_expr().as_coefficients_dict()
    monomios = f_3.keys()
    coef1 = f_3.values()
    coef1_2 = []
    for coeff in coef1:
        coef1_2.append(coeff)
    
    f_4 = 0
    k = 0
    for monos in monomios: 
        list_grados = []
        for var in list_variables:
            grado = degree(monos,gen = var)               
            list_grados.append(grado)
        if np.sum(list_grados)<=order:
            f_4 = f_4 + monos*coef1_2[k]
        k = k + 1
    return f_4

#%%

order = 4
K_x,coef = K_var_no_central([y,a],order)
eq1 = Poly(xpunto,x,y).subs(x,K_x)
eq1_filtrada = filtro_f(eq1,[x,y,a],order)

eq2 = Poly(diff(K_x,y)*ypunto,x,y,a)
eq2 = eq2.subs(x,K_x)
eq2_filtrada = filtro_f(eq2,[x,y,a],order)

eq_final = Poly(eq1_filtrada - eq2_filtrada,x,y,a)
eq_final_info =  eq_final.as_expr().as_coefficients_dict()
monomios = eq_final_info.keys()
monomios2 = []
for monos in monomios: 
    monomios2.append(monos)

list_res = {}
for j in range(len(monomios2)):
    monos = monomios2[j]
    for k in range(len(coef)):
        coeff = coef[k]
        res = solve(monos,coeff)
        if len(res)>0:
            if coeff not in list_res:
                list_res[coeff] = res[0]
                
m = len(monomios2)
n = len(list_res)
if m>n:
    for l in range(n-1,m):
        monos = monomios2[l]
        for k in range(len(coef)):
            coeff = coef[k]
            res = solve(monos,coeff)
            print(monos,res)
            
#%%

list_res_keys = []
for coeff in list_res.keys():
    list_res_keys.append(coeff)

list_res_values = []
for val in list_res.values():
    list_res_values.append(val)             

#reordenar con los coeficientes en orden creciente:
list_res_values2 = []
indices = []
for k in range(len(list_res_keys)):
    for j in range(len(coef)):
        if coef[j]==list_res_keys[k]:
            list_res_values2.append(list_res_values[j])
            indices.append(j+1)

#no me sale 

#%%

print('')
item = 'c)'
print(item)

"""
Considere un entorno del origen. En términos del método de la variedad central:
¿para qué valores del parámetro espera poder reducir la dinámica a una descripción
unidimensional? ¿Por qué?
"""

v1 = [1,0]
v2 = [0,1]

def retrato_fases_var_central(xpunto,ypunto,a_ev,xpf_no_hiper,ypf_no_hiper,epsilon,item):
    """
    retrato_fases devuelve un grafico del retrato de fases
    parametros: (solo para dimension 2D y para 2 parametros no nulos)
        -xpunto,ypunto del problema
        -a_ev es el valor del parametro a
        -K es la funcion de la variedad no central
        -item para guardar el grafico con el nombre del item
        -minx,maxx = -3,3 por default
        (se cambian automaticamente para cubrir los puntos fijos)
    """

    a_ev = float(a_ev)
    xpf = xpf_no_hiper
    ypf = ypf_no_hiper
    
    print(a_ev,xpf,ypf)
    print('')
    
    K_x_manual = (2*a-1)*(y-ypf)**2 + xpf #K_x calculado en el overleaf           
    
    epsi = epsilon
    minix,maxx = xpf-epsi,xpf+epsi
    miniy,maxy = ypf-epsi,ypf+epsi
 
    def eval3(x,y,a,expr):
        return expr(x,y,a)
    
    xdot = lambdify([x,y,a], xpunto)
    ydot = lambdify([x,y,a], ypunto)
    
    K_x_manual3 = K_x_manual.subs(a,a_ev)
    
    n = int(1e3)
    ejey = np.linspace(miniy,maxy,n)
    var1 = []
    for valuey in ejey: 
        var1.append(K_x_manual3.subs(y,valuey))
                      
    info = '$\dot{x}$=' + str(xpunto) + ', $\dot{y}$=' + str(ypunto)
    title = '%s Retrato de fases con a = %.4f'%(item,a_ev)
    print('Graficar retrato de fase') 
    plt.figure(figsize=tamfig)
    plt.title(info+'\n'+title,fontsize=tamtitle)

    nzoom = 900   
    ys_list = np.linspace(miniy,maxy,nzoom)
    xs_list = np.linspace(minix,maxx,nzoom)
    X,Y = np.meshgrid(xs_list,ys_list)
    U = eval3(X,Y,a_ev,xdot)
    V = eval3(X,Y,a_ev,ydot)   
    
    plt.streamplot(X, Y, U, V, density=[3, 2], minlength=.1,color='green') 
    plt.plot(np.linspace(minix,maxx,nzoom),np.ones(nzoom)*0,'k-',label = 'v1:(1,0)')
    plt.plot(np.ones(nzoom)*0,np.linspace(miniy,maxy,nzoom),'k-',label = 'v2:(0,1)')
    plt.plot(var1,ejey,'m.',ms=10,label = 'var central:'+str(K_x_manual))
    plt.plot(xpf,ypf,'.',ms=10,label = 'punto fijo no hiperbolico')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
    plt.ylabel('y',fontsize=tamletra)
    plt.xlabel('x',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.xlim([minix,maxx])
    plt.ylim([miniy,maxy])
    if guardar_graficos==1:
        os.chdir(path_g)
        print('Guardar grafico')
        plt.savefig('item_%s_zoom' %(item))
        os.chdir(path)
            
    return [np.mean(U),np.mean(V),np.mean(var1),np.mean(ejey)]

diff_u = []
diff_v = []
entorno_a1 = [-1e-1,-1e-2,-1e-3,0,1e-3,1e-2,1e-1]
entorno_a2 = np.array(entorno_a1) -1
entorno_a2 = [-1-1e-1,-1-1e-2,-1-1e-3,-1,-1+1e-3,-1+1e-2,-1+1e-1]

# Na = int(5*1e2)
# entorno_a = np.linspace(-5,5,Na)
j=0 
for a_ev in entorno_a1:
    if 0 in x_pf_no_hiper and 0 in y_pf_no_hiper and 0 in a_no_hiper:
        ittem = item + '_%i' %(j+1)

        epsi = 0.5

        xpfnohiper,ypnohiper = x_pf_no_hiper[1],y_pf_no_hiper[1]
        [u1,v1,u2,v2] = retrato_fases_var_central(xpunto,ypunto,a_ev,xpfnohiper,ypnohiper,epsi,ittem)
        print(a_ev)
        diff_u_value,diff_v_value = np.abs(u1-u2),np.abs(v1-v2)
        diff_u.append(diff_u_value)
        diff_v.append(diff_v_value)
        print(diff_u_value)
        print(diff_v_value)
        print('')
        j = j + 1 
        
        
for a_ev in entorno_a2:
    if 0 in x_pf_no_hiper and 0 in y_pf_no_hiper and 0 in a_no_hiper:
        ittem = item + '_%i' %(j+1)

        epsi = 5*1e-10
            
        xpfnohiper,ypnohiper = x_pf_no_hiper[0],y_pf_no_hiper[0]
        [u1,v1,u2,v2] = retrato_fases_var_central(xpunto,ypunto,a_ev,xpfnohiper,ypnohiper,epsi,ittem)
        print(a_ev)
        diff_u_value,diff_v_value = np.abs(u1-u2),np.abs(v1-v2)
        diff_u.append(diff_u_value)
        diff_v.append(diff_v_value)
        print(diff_u_value)
        print(diff_v_value)
        print('')
        j = j + 1 

#%%

# diff_tot = []
# for j in range(len(diff_v)):
#     diff_tot.append(diff_v[j]+diff_u[j])
    
# info = '$\dot{x}$=' + str(xpunto) + ', $\dot{y}$=' + str(ypunto)
# title = '%s Retrato de fases con a = %.2e'%(item,a_ev)
# plt.figure(figsize=tamfig)
# plt.title(info+'\n'+title,fontsize=tamtitle)
# # plt.plot(entorno_a,diff_u,'.m',ms=10,label = '|mean(U)-mean(K(y,a))|')
# # plt.plot(entorno_a,diff_v,'.b',label = '|mean(V)-mean(y)|')
# plt.plot(entorno_a,diff_tot,'.g',label = '|mean(V)-mean(y)|+|mean(U)-mean(K(y,a))|')
# plt.ylabel('Comparacion',fontsize=tamletra)
# plt.xlabel('a',fontsize=tamletra)
# plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
# plt.tick_params(labelsize = tamnum)
# plt.grid(1)
# if guardar_graficos==1:
#         os.chdir(path_g)
#         print('Guardar grafico')
#         plt.savefig('item_%s_compar' %(item))
#         os.chdir(path)

#%%

print('')
item = 'e)'
print(item)

"""
Dibuje retratos de fase compatibles con la información obtenida para las distintas
condiciones según el parámetro a cerca de la bifurcación y [opcional]: compárelos con
retratos de fase obtenidos por integración numérica
"""
a_ev = -1

if a_ev==0:
    xo,yo = 0,0
    xpf_nohiper,ypf_nohiper = 0,0
elif a_ev==-1:
    xo,yo = 1,0
    xpf_nohiper,ypf_nohiper = 1,0
else:
    raise TypeError('invalid value for a')
    
def system(array_y,t):    
    xi,yi = array_y
    dxdt = -xi + xi**2 - yi**2
    dydt = a_ev*yi - yi**3 + xi*yi
    return [dxdt,dydt]

if a_ev==0:
    list_cond_initial = [[-2,-2,2,2],[2,-2,-2,2],[-2,2,2,-2],[2,2,-2,-2]]
    epsi = 0.05
elif a_ev==-1:
    list_cond_initial1 = [[-0.8,-0.2,0.8,0.2],[0.8,-0.2,-0.8,0.2],[-0.8,0.2,0.8,-0.2],[0.8,0.2,-0.8,-0.2]] 
    list_cond_initial2 = [[1.001,0.001,0.999,-0.001],[1.001,0.001,0.999,-0.001],[1.001,0.001,0.999,-0.001],[1.001,0.001,0.999,-0.001]]
    list_cond_initial3 = [[1.2,0.001,0.999,-0.001]]
    list_cond_initial = list_cond_initial1
    if list_cond_initial==list_cond_initial1:
        print('Por las condiciones iniciales dadas, convergemos al (0,0)')
        epsi = 0
    if list_cond_initial==list_cond_initial2:
        print('Por las condiciones iniciales dadas, convergemos al (1,0)')
        epsi = 0
    
for j in range(len(list_cond_initial)):
    ittem = item + '_%i' %(j+1)
    x01,y01,x02,y02 = list_cond_initial[j]
    cond_initial1 = [x01,y01]
    cond_initial2 = [x02,y02]
    
    dt = 0.001
    tf = 50
    t1 = np.arange(-tf, 0, dt)
    t2 = np.arange(0, tf, dt)
    
    sol1 = odeint(system, cond_initial1, t1)
    xt1 = sol1[:, 0]
    yt1 = sol1[:, 1]
    minx1,miny1 = np.min(xt1),np.min(yt1)
    maxx1,maxy1 = np.max(xt1),np.max(yt1)
    
    sol2 = odeint(system, cond_initial2, t2)
    xt2 = sol2[:, 0]
    yt2 = sol2[:, 1]
    minx2,miny2 = np.min(xt2),np.min(yt2)
    maxx2,maxy2 = np.max(xt2),np.max(yt2)
    
    minx = np.min([minx1,minx2])
    miny = np.min([miny1,miny2])
    maxx = np.max([maxx1,maxx2])
    maxy = np.max([maxy1,maxy2])
    
    n = int(3*1e2)    
    ys_list = np.linspace(miny,maxy,n)
    #maxx = 2
    xs_list = np.linspace(minx,maxx,n)
    X,Y = np.meshgrid(xs_list,ys_list)
    info = '$\dot{x}$=' + str(xpunto) + ', $\dot{y}$=' + str(ypunto)
    title = '%s Retrato de fases con a = %i'%(item,a_ev)
    plt.figure(figsize=tamfig)
    plt.title(info+'\n'+title,fontsize=tamtitle)
    DX, DY = system([X, Y], 0)
    plt.streamplot(X, Y, DX, DY, density=[2,2], minlength=.1,color='green') 
    if a_ev ==0:
        plt.plot(xt1, yt1,'m.',ms = 8,label = 'int t<0 con x0,y0 = %i,%i'%(x01,y01))
        plt.plot(xt2, yt2,'r.',ms = 8,label = 'int t>0 con x0,y0 = %i,%i'%(x02,y02))
    elif a_ev==-1:
        plt.plot(xt1, yt1,'m.',ms = 8,label = 'int t<0 con x0,y0 = %.1f,%.1f'%(x01,y01))
        plt.plot(xt2, yt2,'r.',ms = 8,label = 'int t>0 con x0,y0 = %.1f,%.1f'%(x02,y02))        
    
    plt.plot(xpf_nohiper,ypf_nohiper,'.k',ms=15)
    plt.ylabel('y',fontsize=tamletra)
    plt.xlabel('x',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)   
    plt.xlim([minx-epsi,maxx+epsi])
    plt.ylim([miny-epsi,maxy+epsi])
    plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
    if guardar_graficos==1:
            os.chdir(path_g)
            print('Guardar grafico')
            if a_ev== -1:
                if list_cond_initial == list_cond_initial1:
                    plt.savefig('item_%s_a%i_cond1' %(ittem,a_ev))
                elif list_cond_initial == list_cond_initial2:
                    plt.savefig('item_%s_a%i_cond2' %(ittem,a_ev))      
            elif a_ev ==0:
                plt.savefig('item_%s_a%i' %(ittem,a_ev))    
            os.chdir(path)

#%%

print('')
item = 'e)extra'
print(item)

a_ev = -1
xo,yo = 1,0
xpf_nohiper,ypf_nohiper = 0,0
    
xpunto2 = -(x-xo) + (x-xo)*(x-xo) - (y-yo)*(y-yo)
ypunto2 = a*(y-yo) - (y-yo)*(y-yo)*(y-yo) + (x-xo)*(y-yo)

barrido_a = [-1.5,-1,-0.5,0,0.5,1]
j = 0
for value_a in barrido_a:
    ittem = item + '_%i' %(j+1)
    retrato_fases(xpunto2,ypunto2,value_a,ittem)
    j = j + 1

#%%

def system(array_y,t):    
    xi,yi = array_y
    dxdt = -xi + xi**2 - yi**2
    dydt = a_ev*yi - yi**3 + xi*yi
    return [dxdt-xo,dydt-yo]

if a_ev==0:
    list_cond_initial = [[-2,-2,2,2],[2,-2,-2,2],[-2,2,2,-2],[2,2,-2,-2]]
    epsi = 0.05
elif a_ev==-1:
    list_cond_initial = [[-0.8,-0.2,0.8,0.2],[0.8,-0.2,-0.8,0.2],[-0.8,0.2,0.8,-0.2],[0.8,0.2,-0.8,-0.2]]
    list_cond_initial = [[1.1,0.1,0.9,-0.1],[1.1,0.1,0.9,-0.1],[1.1,0.1,0.9,-0.1],[1.1,0.1,0.9,-0.1]]
    epsi = 0.01
    
for j in range(len(list_cond_initial)):
    ittem = item + '_%i' %(j+1)
    x01,y01,x02,y02 = list_cond_initial[j]
    cond_initial1 = [x01,y01]
    cond_initial2 = [x02,y02]
    
    dt = 0.001
    tf = 50
    t1 = np.arange(-tf, 0, dt)
    t2 = np.arange(0, tf, dt)
    
    sol1 = odeint(system, cond_initial1, t1)
    xt1 = sol1[:, 0]
    yt1 = sol1[:, 1]
    minx1,miny1 = np.min(xt1),np.min(yt1)
    maxx1,maxy1 = np.max(xt1),np.max(yt1)
    
    sol2 = odeint(system, cond_initial2, t2)
    xt2 = sol2[:, 0]
    yt2 = sol2[:, 1]
    minx2,miny2 = np.min(xt2),np.min(yt2)
    maxx2,maxy2 = np.max(xt2),np.max(yt2)
    
    minx = np.min([minx1,minx2])
    miny = np.min([miny1,miny2])
    maxx = np.max([maxx1,maxx2])
    maxy = np.max([maxy1,maxy2])
    
    n = int(3*1e2)
    ys_list = np.linspace(miny,maxy,n)
    maxx = 2
    xs_list = np.linspace(minx,maxx,n)
    X,Y = np.meshgrid(xs_list,ys_list)
    info = '$\dot{x}$=' + str(xpunto2) + ', $\dot{y}$=' + str(ypunto2)
    title = '%s Retrato de fases con a = %i'%(item,a_ev)
    plt.figure(figsize=tamfig)
    plt.title(info+'\n'+title,fontsize=tamtitle)
    DX, DY = system([X, Y], 0)
    plt.streamplot(X, Y, DX, DY, density=[2,2], minlength=.1,color='green') 
    if a_ev ==0:
        plt.plot(xt1, yt1,'m.',ms = 8,label = 'int t<0 con x0,y0 = %i,%i'%(x01,y01))
        plt.plot(xt2, yt2,'r.',ms = 8,label = 'int t>0 con x0,y0 = %i,%i'%(x02,y02))
    elif a_ev==-1:
        plt.plot(xt1, yt1,'m.',ms = 8,label = 'int t<0 con x0,y0 = %.1f,%.1f'%(x01,y01))
        plt.plot(xt2, yt2,'r.',ms = 8,label = 'int t>0 con x0,y0 = %.1f,%.1f'%(x02,y02))        
    
    plt.plot(xpf_nohiper,ypf_nohiper,'.k',ms=15)
    plt.plot(1,0,'.k',ms=15)
    plt.ylabel('y',fontsize=tamletra)
    plt.xlabel('x',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.xlim([minx-epsi,maxx+epsi])
    plt.ylim([miny-epsi,maxy+epsi])
    plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
    if guardar_graficos==1:
            os.chdir(path_g)
            print('Guardar grafico')
            plt.savefig('item_%s_a%i' %(ittem,a_ev))
            os.chdir(path)            

#%%

