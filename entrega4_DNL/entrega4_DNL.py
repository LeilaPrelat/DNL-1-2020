#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 18:06:23 2020

@author: leila

4 entrega DNL: Mapas

http://materias.df.uba.ar/dnla2020c1/files/2020/06/DNL_mapas.pdf

"""
import matplotlib.pyplot as plt
import numpy as np
import os 
from scipy.optimize import fsolve
# import math
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib as mpl
# from scipy import signal
# from scipy.integrate import odeint

#%%

tamfig = (11,9)
tamlegend = 16
tamletra = 18
tamtitle = 18
tamnum = 16

path =  os.getcwd()

guardar_graficos = 1

if guardar_graficos == 1:
    path_g = path + '/' + 'graficos'
    if not os.path.exists(path_g):
        print('Crear carpeta para guardar graficos')
        os.mkdir(path_g)

#%%

# print('intento de construir el mapa de periodo 1')

# phi_o = 0 #condicion inicial 
# epsilon = 0.05
# T_To = 10

# omega0 = 1
# T = 10
# nmax = 50

# def mapa_f(ti,tf):
#     phif = phi_o
#     delta_t = tf-ti
#     intervalo = int(delta_t/5)
#     list_t = np.linspace(ti,tf,intervalo+1)
#     list_x = []
#     list_y = []
#     list_z = []
#     for n in list_t:
#         n = int(n)
#         phif = phif + epsilon*np.cos(phif) + 2*np.pi*T_To
#         x,y = np.cos(phif),np.sin(phif)
#         z = np.abs(x+1j*y)
#         list_x.append(x)
#         list_y.append(y)
#         list_z.append(z)
#     return list_x,list_y,list_t


# ti = 500
# tf = 400000
# fig = plt.figure(figsize=tamfig)
# ax = fig.gca(projection='3d')
# mpl.rcParams['legend.fontsize'] = tamlegend
# x,y,z = mapa_f(ti,tf)
# ax.plot(x,y,z,'.', label='parametric curve')
# ax.legend(loc='best',markerscale=2)

# plt.figure()
# plt.plot(x,y,'.')

# #%% 

# print('intento de construir el mapa de periodo 1')

# omega0 = 2
# epsilon = 1
# T = 1
# nmax = 100

# def system(coordenadas,t):    
#     x,y = coordenadas
#     delta_tot = 0
#     list_n = np.linspace(-nmax,nmax,2*nmax+1)
    
#     def delta_t(n):
#         if t == n*T:
#             print('golpe')
#             return 1
#         else:
#             return 0
    
#     for n in list_n:
#         n = int(n)
#         delta_tot = delta_tot - delta_t(n)
        
#     z = x + 1j*y
#     z2 = np.abs(z)**2
#     eq1 = x - omega0*y - x*z2
#     eq2 = omega0*x + y - y*z2 + epsilon*delta_tot

#     return [eq1,eq2]

# ti = 0
# tf = 100
# dt = 0.5
# list_t = np.arange(ti, tf, dt)
# cond_initial = [mapa_f(ti,tf)[0][0],mapa_f(ti,tf)[1][0]]
# sol = odeint(system, cond_initial, list_t)
# # sol = np.mod(sol+np.pi,2*np.pi)-np.pi
# xt = sol[:, 0]
# yt = sol[:, 1]

# plt.figure()
# plt.plot(xt,yt)

# fig = plt.figure(figsize=tamfig)
# ax = fig.gca(projection='3d')
# mpl.rcParams['legend.fontsize'] = tamlegend
# x,y,z = mapa_f(ti,tf)
# ax.plot(xt,yt,list_t,'.', label='parametric curve')
# ax.legend(loc='best',markerscale=2)


#%%

item = 'g)'
print('hallemos los puntos fijos')

phi_o = 0 #condicion inicial 
epsilon = 9
T_To = 3

def phi_next(phi_before):
    phif = phi_before + epsilon*np.cos(phi_before) + 2*np.pi*T_To
    return phif

N = int(1e3)
phi_before = np.linspace(-2*np.pi,2*np.pi,N)
phi_next2 = phi_next(phi_before)

plt.figure(figsize=tamfig)
plt.title('g) Puntos fijos con $\epsilon$ = %.2f, $T/T_o$ = %i'%(epsilon,T_To),fontsize=tamtitle)
plt.plot(phi_before,phi_before + 4*np.pi,'-',color = 'purple',lw = 3,label = '$\phi_{n+1} = \phi_n + 4\pi$' )
plt.plot(phi_before,phi_before + 2*np.pi,'-',color = 'orange',lw = 3,label = '$\phi_{n+1} = \phi_n + 2\pi$' )
# plt.plot(phi_before,phi_before,'-',color = 'red',lw = 3,label = '$\phi_{n+1} = \phi_n$' )
plt.plot(phi_before,phi_next2,'.', color = 'green',ms = 8)
plt.xlabel('$\phi_n$',fontsize=tamletra)
plt.ylabel('$\phi_{n+1}$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s' %(item))
    os.chdir(path)


#%% item f)

"""
En un diagrama $(T/T_o,\varepsilon)$ grafique las condiciones halladas en (f) y (h). 
¿Qué significado físico tienen estas soluciones?
"""

item = 'i)f)'

print(item)
print('Orbita de periodo 1')
print('Formula de phi_n y condiciones del item f)')

def phi_f1(T_To,epsilon):
    cte = 2*np.pi/epsilon
    cte2 = 1 - T_To
    
    # cota_inf = 1 - epsilon/(2*np.pi)
    # cota_sup = 1 + epsilon/(2*np.pi)
    
    rta = np.arccos(cte*cte2)/np.pi
    return rta

def phi_f2(T_To,epsilon):
    cte = 2*np.pi/epsilon
    cte2 = 1 - T_To
    
    # cota_inf = 1 - epsilon/(2*np.pi)
    # cota_sup = 1 + epsilon/(2*np.pi)
    
    rta = -np.arccos(cte*cte2)/np.pi
    return rta

N = int(1e3)
maxi = 6
list_T_To = np.linspace(-maxi,maxi,N)
list_epsilon = np.linspace(-maxi,maxi,N)

X,Y = np.meshgrid(list_T_To,list_epsilon)
f1 = np.vectorize(phi_f1)
Z1 = f1(X, Y)

f2 = np.vectorize(phi_f2)
Z2 = f2(X, Y)


# plt.clim(0,2)

def cotas_f(epsilon):
    cota_inf = 1 - epsilon/(2*np.pi) #T_Tocrit
    cota_sup = 1 + epsilon/(2*np.pi) #T_Tocrit
    return cota_inf,cota_sup

list_cota_inf = []
list_cota_sup = []
for epsi in list_epsilon:
    list_cota_inf.append(cotas_f(epsi)[0])
    list_cota_sup.append(cotas_f(epsi)[1])

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del f) Orbitas de periodo 1',fontsize=int(tamtitle*0.9))
im1 = plt.imshow(Z1, extent = limits,interpolation='bilinear')
cbar1 = plt.colorbar(im1)
cbar1.set_label('$\phi_{n,+}/\pi$',size=tamlegend)
plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    # plt.tight_layout(1)
    plt.savefig('item_%s+' %(item))
    os.chdir(path)

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del f) Orbitas de periodo 1',fontsize=int(tamtitle*0.9))
im2 = plt.imshow(Z2, extent = limits,interpolation='bilinear')
cbar2 = plt.colorbar(im2)
cbar2.set_label('$\phi_{n,-}/\pi$',size=tamlegend)
plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    # plt.tight_layout(1)
    plt.savefig('item_%s-' %(item))
    os.chdir(path)
    
print('')

#%% item g)

item = 'i)g)'
print('Orbita de periodo 1 + cotas del item g) (bif. de dup. de periodo)')

print(item)

N = int(1e3)
maxi = 6
list_T_To = np.linspace(-maxi,maxi,N)
list_epsilon = np.linspace(-maxi,maxi,N)

X,Y = np.meshgrid(list_T_To,list_epsilon)
f1 = np.vectorize(phi_f1)
Z1 = f1(X, Y)

f2 = np.vectorize(phi_f2)
Z2 = f2(X, Y)


# plt.clim(0,2)

def cotas_g(T_To):
    cte1 = (2*np.pi)**2
    cte2 = (1-T_To)**2
    cota_inf = -np.sqrt(4+cte1*cte2) #epsicrit
    cota_sup = np.sqrt(4+cte1*cte2) #epsicrit
    return cota_inf,cota_sup

list_cota_inf = []
list_cota_sup = []
for epsi in list_epsilon:
    list_cota_inf.append(cotas_f(epsi)[0])
    list_cota_sup.append(cotas_f(epsi)[1])
    
list_cota_inf2 = []
list_cota_sup2 = []
for T_To in list_T_To:
    list_cota_inf2.append(cotas_g(T_To)[0])
    list_cota_sup2.append(cotas_g(T_To)[1])
    

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del f) y g)',fontsize=int(tamtitle*0.9))
im1 = plt.imshow(Z1, extent = limits,interpolation='bilinear')
cbar1 = plt.colorbar(im1)
cbar1.set_label('$\phi_{n,+}/\pi$',size=tamlegend)    
plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_T_To,list_cota_inf2,color = 'orange',lw = 3,label = 'curva - g)')
plt.plot(list_T_To,list_cota_sup2,color = 'magenta',lw = 3, label = 'curva + g)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s+' %(item))
    os.chdir(path)

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del f) y g)',fontsize=int(tamtitle*0.9))
im2 = plt.imshow(Z2, extent = limits,interpolation='bilinear')
cbar2 = plt.colorbar(im2)
cbar2.set_label('$\phi_{n,-}/\pi$',size=tamlegend)
plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_T_To,list_cota_inf2,color = 'orange',lw = 3,label = 'curva - g)')
plt.plot(list_T_To,list_cota_sup2,color = 'magenta',lw = 3, label = 'curva + g)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s-' %(item))
    os.chdir(path)
    
print('')

#%%  item h)

item = 'i)h)1'
print('Orbita de periodo 2')
print('Formula de phi_n del item h) despreciando epsilon^2 y condiciones del item h)')

print(item)

def cotas_h(epsilon):
    cota_inf = 0.5 - epsilon/(2*np.pi) #T_Tocrit
    cota_sup = 0.5 + epsilon/(2*np.pi) #T_Tocrit
    return cota_inf,cota_sup

def phi_h1(T_To,epsilon):
    cte1 = np.pi - 2*np.pi*T_To
    cte3 = np.pi*T_To
    cte2 = epsilon*np.cos(cte3)
    
    rta = np.arccos(cte1/cte2) - cte3
    rta2 = rta/np.pi
    return rta2

def phi_h2(T_To,epsilon):
    cte1 = np.pi - 2*np.pi*T_To
    cte3 = np.pi*T_To
    cte2 = epsilon*np.cos(cte3)
    
    rta = -np.arccos(cte1/cte2) - cte3
    rta2 = rta/np.pi
    return rta2

N = int(2*1e3)
maxi = 5
list_T_To = np.linspace(-maxi,maxi,N)
list_epsilon = np.linspace(-maxi,maxi,N) #epsilon tiene que ser chico

list_cota_inf = []
list_cota_sup = []
for epsi in list_epsilon:
    list_cota_inf.append(cotas_f(epsi)[0])
    list_cota_sup.append(cotas_f(epsi)[1])

list_cota_inf2 = []
list_cota_sup2 = []
for epsi in list_epsilon:
    list_cota_inf2.append(cotas_h(epsi)[0])
    list_cota_sup2.append(cotas_h(epsi)[1])

X,Y = np.meshgrid(list_T_To,list_epsilon)
f1 = np.vectorize(phi_h1)
Z1 = f1(X, Y)

f2 = np.vectorize(phi_h2)
Z2 = f2(X, Y)

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del h) Orbitas de periodo 2 con aprox',fontsize=int(tamtitle*0.9))
im1 = plt.imshow(Z1, extent = limits,interpolation='bilinear')
cbar1 = plt.colorbar(im1)
cbar1.set_label('$\phi_{n,+}/\pi$',size=tamlegend)
# plt.clim(0,2)

plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_cota_inf2,list_epsilon, color = 'orange',lw = 3,label = 'cota inferior h)')
plt.plot(list_cota_sup2,list_epsilon, color = 'steelblue',lw = 3, label = 'cota superior h)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s+' %(item))
    os.chdir(path)
    

plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del h) Orbitas de periodo 2 con aprox',fontsize=int(tamtitle*0.9))
im2 = plt.imshow(Z2, extent = limits,interpolation='bilinear')
cbar2 = plt.colorbar(im2)
cbar2.set_label('$\phi_{n,-}/\pi$',size=tamlegend)
# plt.clim(0,2)

plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_cota_inf2,list_epsilon, color = 'orange',lw = 3,label = 'cota inferior h)')
plt.plot(list_cota_sup2,list_epsilon, color = 'steelblue',lw = 3, label = 'cota superior h)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='lower right',markerscale=2,fontsize=int(tamlegend*0.9))
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s-' %(item))
    os.chdir(path)    

print('')

#%%

item = 'i)h)2'
print('Formula de phi_n del item h) sin despreciar epsilon^2 y condiciones del item g)')

print(item)

def phi_h_sin_aprox1(T_To,epsilon):
    cte1 = 2*np.pi - 4*np.pi*T_To
    def equation(x):
        f = epsilon*np.cos(x) + epsilon*np.cos(x + epsilon*np.cos(x) + 2*np.pi*T_To) - cte1
        return f

    cond_inicial = phi_h1(T_To,epsilon)*np.pi
    # if math.isnan(cond_inicial) == False: 
    #     cond_inicial = cond_inicial 
    # else:
    #     cond_inicial = 0.5
    sol = fsolve(equation, cond_inicial)
    return sol

def phi_h_sin_aprox2(T_To,epsilon):
    cte1 = 2*np.pi - 4*np.pi*T_To
    def equation(x):
        f = epsilon*np.cos(x) + epsilon*np.cos(x + epsilon*np.cos(x) + 2*np.pi*T_To) - cte1
        return f

    cond_inicial = phi_h2(T_To,epsilon)*np.pi
    # if math.isnan(cond_inicial) == False: 
    #     cond_inicial = cond_inicial 
    # else:
    #     cond_inicial = 0.5
    sol = fsolve(equation, cond_inicial)
    return sol

N = int(1e3)
maxi = 5
X,Y = np.meshgrid(list_T_To,list_epsilon, sparse=True)
list_T_To = np.linspace(-maxi,maxi,N)
list_epsilon = np.linspace(-maxi,maxi,N) #epsilon tiene que ser chico

list_cota_inf = []
list_cota_sup = []
for epsi in list_epsilon:
    list_cota_inf.append(cotas_f(epsi)[0])
    list_cota_sup.append(cotas_f(epsi)[1])

list_cota_inf2 = []
list_cota_sup2 = []
for epsi in list_epsilon:
    list_cota_inf2.append(cotas_h(epsi)[0])
    list_cota_sup2.append(cotas_h(epsi)[1]) 

f1 = np.vectorize(phi_h_sin_aprox1)
Z1 = f1(X, Y)

f2 = np.vectorize(phi_h_sin_aprox2)
Z2 = f2(X, Y)


plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del h) Orbitas de periodo 2 sin aprox',fontsize=int(tamtitle*0.9))
im1 = plt.imshow(Z1, extent = limits,interpolation='bilinear')
cbar1 = plt.colorbar(im1)
cbar1.set_label('$\phi_{n,+}/\pi$',size=tamlegend)

plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_cota_inf2,list_epsilon, color = 'orange',lw = 3,label = 'cota inferior h)')
plt.plot(list_cota_sup2,list_epsilon, color = 'steelblue',lw = 3, label = 'cota superior h)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)

if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s+' %(item))
    os.chdir(path)


plt.figure(figsize=tamfig)
limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
plt.xlabel('$T/T_o$',fontsize=tamletra)
plt.ylabel('$\epsilon$',fontsize=tamletra)
plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del h) Orbitas de periodo 2 sin aprox',fontsize=int(tamtitle*0.9))
im2 = plt.imshow(Z2, extent = limits,interpolation='bilinear')
cbar2 = plt.colorbar(im2)
cbar2.set_label('$\phi_{n,-}/\pi$',size=tamlegend)

plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior f)')
plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior f)')
plt.plot(list_cota_inf2,list_epsilon, color = 'orange',lw = 3,label = 'cota inferior h)')
plt.plot(list_cota_sup2,list_epsilon, color = 'steelblue',lw = 3, label = 'cota superior h)')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
plt.ylim([-maxi,maxi])
plt.xlim([-0.5,maxi])
# plt.grid(1)
if guardar_graficos==1:
    os.chdir(path_g)
    print('Guardar grafico')
    plt.tight_layout(1)
    plt.savefig('item_%s-' %(item))
    os.chdir(path)
    
#%%

# def cond_no_periodo1(T_To,epsilon):
#     phi = phi_h(T_To,epsilon)*np.pi
#     T_Tocrit = 1 - np.cos(phi)*epsilon/(2*np.pi)
#     return T_Tocrit

# f2 = np.vectorize(cond_no_periodo1)
# Z2 = f(X, Y)

# plt.figure(figsize=tamfig)
# limits = [min(list_T_To) , max(list_T_To), min(list_epsilon) , max(list_epsilon)]
# plt.xlabel('$T/T_o$',fontsize=tamletra)
# plt.ylabel('$\epsilon$',fontsize=tamletra)
# plt.title('Item i) Diagrama ($T/T_o,\epsilon$) del h)',fontsize=int(tamtitle*0.9))
# im = plt.imshow(Z2, extent = limits,  cmap='RdBu',interpolation='bilinear')
# cbar = plt.colorbar(im)
# cbar.set_label('$\phi_n/\pi$',size=tamlegend)
# plt.clim(0,2)

# plt.plot(list_cota_inf,list_epsilon, color = 'darkred',lw = 3,label = 'cota inferior')
# plt.plot(list_cota_sup,list_epsilon, color = 'green',lw = 3, label = 'cota superior')
# plt.tick_params(labelsize = tamnum)
# plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
# ejex = np.linspace(np.min(list_T_To),np.max(list_T_To),5)
# ejey = np.linspace(np.min(list_epsilon),np.max(list_epsilon),5)
# plt.plot(ejex,np.ones(5)*0,'k-',lw = 2)
# plt.plot(np.ones(5)*0,ejey,'k-',lw = 2)
# # plt.grid(1)
# if guardar_graficos==1:
#     os.chdir(path_g)
#     print('Guardar grafico')
#     plt.savefig('item_%s4' %(item))
#     os.chdir(path)

#%%

