# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 21:38:20 2020

@author: Mathis Granger
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial
from scipy.special import lpmv

#Angle polaire et azimuthale

theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)

#Création d'un repère shérique en 2D dans le plan theta phi.

theta, phi = np.meshgrid(theta, phi)

#Coordonnées cartésiennes dans cette grilles (transformations du cartésiens en sphériques)

xyz = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])



def SimuSphHarm (ax, l, m, ri):

    
            
    
    def shr (m, l, phi, theta):
 #Définition des HS réelles        
        A = (2*l) / (4*np.pi)
        b = l - np.abs(m)
        B = factorial(b)
        c = l + np.abs(m)
        C = factorial(c)
        D = lpmv(m, l, np.cos(theta))
        E = np.cos( m * phi )
        H = (-1)**m
        G = H*(np.sqrt(2) * np.sqrt( A * (B/C)) * D * E)
        
        return G
#Définition des HS imaginaire   
    def shc (m, l, phi, theta):
         
        A = (2*l) / (4*np.pi)
        b = l - np.abs(m)
        B = factorial(b)
        c = l + np.abs(m)
        C = factorial(c)
        D = lpmv(m, l, np.cos(theta))
        E = np.sin( m * phi )
        H = (-1)**m
        G =  H*(np.sqrt(2) * np.sqrt( A * (B/C)) * D * E)
        return G
 #Switch pour obtenir soit la partie réelle ou imaginaire 1 = réel 0  = imaginaire  
    def realimg(ri):
        if ri == 1:
            R = shr (m,l,phi,theta)
        elif ri == 0:
            R = shc (m,l,phi,theta)
        return R
            
    SH = realimg(ri)    
    
    HSx, HSy, HSz = np.abs(SH)*xyz
    
#On colore la surface en fonction du signe de HS, vert = positif et rouge = négatif
    
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('RdYlGn'))
    cmap.set_clim(-0.35, 0.35)
    
    ax.plot_surface(HSx, HSy, HSz,
                    facecolors=cmap.to_rgba(SH.real),
                    rstride=1, cstride=1)
    
#On trace les axes x, y et z     
    ax_lim = 0.5
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    
#Affichafe du titre et suppression du cadre
    
    ax.set_title("$d_{x^{2}-y^2}$")
    ax_lim = 0.5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')

fig = plt.figure(figsize=plt.figaspect(1))
ax = fig.add_subplot(projection='3d')
l, m , ri = 2 , 2 , 0

SimuSphHarm(ax, l, m, ri)

plt.show()
        