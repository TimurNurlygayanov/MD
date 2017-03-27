#
# Copyright 2006-2017 Vladimir Veshnev, Aleksey Geraskin, Marina Nurlygayanova,
# Timur Nurlygayanov
#
# This file is part of MD program.
#
# MD is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MD. If not, see <http://www.gnu.org/licenses/>
#

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange
from pylab import *
import os
import matplotlib.mlab as mlab


colors = {0: 'r', 1: 'c', 2: 'g', 3: 'k', 4: 'y', 5: 'b', 6: 'm'}
selected_particles1 = []
selected_particles2 = []
selected_particles3 = []

def read_data(file_name):
    global im1, im2, im3
    global y1, y2, y3
    global z1, z2, z3
    global selected_particles1, selected_particles2, selected_particles3
    
    im, im1, im2, im3, y, z, y1, y2, y3, z1, z2, z3 = [], [], [], [], [], [], [], [], [], [], [], []
    with open(file_name) as f:
        for line in f:
            if 'SELECTED0' in line:
                print line.split(' ')
                selected_particles1 = []
                selected_particles1 = [int(q) for q in line.split(' ') if (('S' not in q) and (q != '\n'))]
            elif 'SELECTED1' in line:
                print line.split(' ')
                selected_particles2 = []
                selected_particles2 = [int(q) for q in line.split(' ') if (('S' not in q) and (q != '\n'))]
                y1 = y[:]
                y = []
                z1 = z[:]
                z = []
                im1 = im[:]
                im = []
            elif 'SELECTED2' in line:
                selected_particles3 = []
                selected_particles3 = [int(q) for q in line.split(' ') if (('S' not in q) and (q != '\n'))]
                y2 = y[:]
                y = []
                z2 = z[:]
                z = []
                im2 = im[:]
                im = []
            else:
                i_m, yy, zz = line.split()
                im.append(int(i_m))
                y.append(float(yy))
                z.append(float(zz))
    z3 = z
    y3 = y
    im3 = im
    
                
def show_profile(im, y, z, selected_particles, file_name):
    name = ' Profile '
    
    z100, y100 = [], []
    for i, k in enumerate(im):
        if k in selected_particles:
            z100.append(z[i])
            y100.append(y[i])

    my_dpi = 200
    fig = plt.figure(file_name + name, figsize=(5400/my_dpi, 1600/my_dpi),
                     dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title(file_name + '\n' + name)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.grid(True)
    #ax.set_xticks([i for i in xrange(0, 20)], minor=True)
    #ax.set_yticks([i for i in xrange(0, 20)], minor=True)
    ax.grid(which='minor', alpha=0.3)
    #plt.axis([0.0, 0.0, 20.0, 20.0])

    plt.plot(z, y, '{0},'.format('b'))
    plt.plot(z100, y100, '{0},'.format('r'))
    savefig(file_name + ' - ' + name + '.pdf')
    close(fig) 


for file in reversed(os.listdir(".")):
    if '.txt' in file:
        read_data(file)
        short_name = file[:-14]
        
        show_profile(im1, y1, z1, selected_particles1, "1" + file.replace('txt', ''))
        show_profile(im2, y2, z2, selected_particles1, "2" + file.replace('txt', ''))
        show_profile(im3, y3, z3, selected_particles1, "3" + file.replace('txt', ''))
