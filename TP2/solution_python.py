#!/bin/python
import numpy as np
import netCDF4
from netCDF4 import Dataset
import xarray
from scipy import signal
from scipy.ndimage.filters import convolve
from scipy.misc import imread, imshow
import matplotlib.pyplot as plt

'''
Parametros
'''
w = np.array([[-1,-1,-1],[-1,8,-1],[-1,-1,-1]]) #MATRIZ CON LA QUE SE REALIZA LA CONVOLUCION
print w

dataDIR = "/home/anij/facu/2019/1erSemestre/SO2/S02/TP2/data.nc"

'''
Abrir el dataset como una matriz XARRAY y guardar la matriz CMI
'''
DS = Dataset(dataDIR)
f = DS.variables['CMI']
print f

h = f[0:12000,0:12000]
print h[10000:10005,10000:10005]
# k = convolve(h,w)
# print h.shape
# print k[10000:10005,10000:10005]
# plt.imshow(h,cmap="RdBu")

# Colormap blue is not recognized. Possible values are: Accent, Accent_r, Blues,
# Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, 
# Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges,
# Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r,
# PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples,
# Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, 
# RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, 
# Spectral_r, Vega10, Vega10_r, Vega20, Vega20_r, Vega20b, Vega20b_r, Vega20c, 
# Vega20c_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, 
# YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, 
# bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, 
# copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, 
# gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, 
# gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, 
# gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, 
# inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, 
# ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, 
# rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, 
# summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, 
# terrain, terrain_r, viridis, viridis_r, winter, winter_r

plt.imshow(h[0:12000,0:12000],cmap=plt.cm.get_cmap('Spectral_r',2))
plt.savefig("pipa.png", format="png",dpi=900)
'''
kleiner = np.array([[10,2,2,2,2],[2,10,2,2,2],[2,2,10,2,2],[2,2,2,10,2],[2,2,2,2,10]])
print (kleiner)

Convolucion

g = signal.convolve2d(kleiner,w,boundary='fill',mode='same')
h = convolve(kleiner,w)
y = signal.convolve2d(kleiner,w,boundary='symm',mode='valid')
print(g)
print g.shape
print(h)
print h.shape
print(y)
print y.shape

y = g[10000:10200,10000:10200]
print y.shape
plt.imshow(y)
plt.show()'''
