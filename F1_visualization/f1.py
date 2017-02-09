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


def read_data(file_name):
    global x
    global y
    global z
    
    x, y, z = [], [], []
    with open(file_name) as f:
        for line in f:
            i_m, x_m, y_m, z_m = line.split()
            x.append(float(x_m))
            y.append(float(y_m))
            z.append(float(z_m))


def show_3d_particles():
    fig = plt.figure("3D")

    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.scatter(x, y, z)
    ax.set_xlabel(u'X')
    ax.set_xlim([-1, 1])
    ax.set_xscale('linear')
    ax.set_ylabel(u'Y')
    ax.set_ylim([-1, 1])
    ax.set_yscale('linear')
    ax.set_zlabel(u'Z')
    ax.set_zlim([-1, 1])


def show_2d_xy_particles(file_name, z_min, z_max, x_offset=None):
    x_offset = int(x_offset or 0)
    name = 'F1 XY, Z=[{0};{1}]'.format(z_min, z_max)
    x2, y2 = [], []
    for i, z2 in enumerate(z):
        if (z2 >= z_min and z2 <= z_max):
            x2.append(x[i]+x_offset)
            y2.append(y[i])

    my_dpi = 200
    fig = plt.figure(file_name + name, figsize=(5400/my_dpi, 1600/my_dpi),
                     dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(file_name + '\n' + name)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.grid(True)
    ax.set_xticks([i/10.0 for i in xrange(-10, 10)], minor=True)
    ax.set_yticks([i/10.0 for i in xrange(-10, 10)], minor=True)
    ax.grid(which='minor', alpha=0.3)
    plt.axis([-1.0, x_offset + 1.0, -1.0, 1.0])

    plt.plot(x2, y2, '{0},'.format(colors[x_offset % 7]))
    savefig(file_name + ' - ' + name + '.pdf')
    if x_offset is None:
        close(fig)    


def show_2d_yz_particles(file_name, x_min, x_max):
    name = 'F1 YZ, X=[{0};{1}]'.format(x_min, x_max)
    z2, y2 = [], []
    for i, x2 in enumerate(x):
        if (x2 >= x_min and x2 <= x_max):
            z2.append(z[i])
            y2.append(y[i])

    my_dpi = 200
    fig = plt.figure(file_name + name, figsize=(5400/my_dpi, 1600/my_dpi),
                     dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title(file_name + '\n' + name)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.grid(True)
    ax.set_xticks([i/10.0 for i in range(-10, 10)], minor=True)
    ax.set_yticks([i/10.0 for i in range(-10, 10)], minor=True)
    ax.grid(which='minor', alpha=0.3)
    plt.axis([-1.0, 1.0, -1.0, 1.0])

    plt.plot(z2, y2, 'r,', alpha=0.1)
    savefig(file_name + ' - ' + name + '.pdf')
    close(fig)

def show_projection_by_slice(file_name, pr, arg, arg_min, arg_max,
                             fig_name=None):
    r, f = [], []
    if fig_name:
        name = fig_name
    else:
        name = 'F1 {0}_profile, {1}=[{2};{3}]'.format(pr, arg, arg_min,
                                                      arg_max)

    my_dpi = 200
    fig = plt.figure(file_name + name, figsize=(5400/my_dpi, 1600/my_dpi),
                     dpi=my_dpi)
    ax = fig.add_subplot(1, 2, 1)

    for i in xrange(-100, 100):
        r.append(0)
        f.append(float(i) / 100.0)

    if (pr == 'X'):
        pr_list = x
    elif (pr == 'Y'):
        pr_list = y
    else:
        pr_list = z
    if (arg == 'X'):
        args_list = x
    elif (arg == 'Y'):
        args_list = y
    else:
        args_list = z

    fr = []
    for i, arg2 in enumerate(pr_list):
        if (arg2 > -1.0 and arg2 < 1.0):
            if (args_list[i] >= arg_min and args_list[i] <= arg_max):
                r[100 + int(arg2*100000.0)/1000] += 1
                fr.append(arg2)

    ax.set_xticks([i/20.0 for i in range(-20, 20)], minor=True)
    ax.grid(which='minor', alpha=0.3)
    plt.title(file_name + '\n' + name)
    plt.xlabel(pr)
    
    mu = np.mean(fr)
    sigma = np.std(fr)
    plt.plot(f, r, label='F1 {0}_profile, $\mu$={1}, $\sigma$={2}'.format(pr, mu, sigma))

    """ Add gausian normal distribution graph
        to compare actual distribution with normal distribution
        Important: sigma != sqrt(1.0) - we need to discuss it
    """
    mean = np.mean(fr)
    sigma = np.std(fr)
    xf = np.linspace(-1, 1, len(f))
    yf = mlab.normpdf(xf, mean, sigma)
    yf_max = max(yf)
    for i, _ in enumerate(yf):
        yf[i] *= (max(r)/yf_max)
    label = 'normal distribution, $\mu$={0}, $\sigma$={1}'.format(mean, sigma)
    plt.plot(xf, yf, label=label)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    savefig(file_name + ' - ' + name + '.pdf')
    close(fig)


x, y, z = [], [], []
k = 0
for file in reversed(os.listdir(".")):
    if '.txt' in file:
        read_data(file)
        k += 1
        short_name = file[:-14]
        dw = 0.1
        for w in arange(-dw/2.0, dw/2.0, dw):
            show_2d_xy_particles(short_name, w, w+dw)
            show_2d_xy_particles("test", w, w+dw, k)
            show_2d_yz_particles(short_name, w, w+dw)
            
            show_projection_by_slice(short_name, 'X', 'Y', w, w+dw)
            show_projection_by_slice(short_name, 'X', 'Z', w, w+dw)
            show_projection_by_slice(short_name, 'Y', 'X', w, w+dw)
            show_projection_by_slice(short_name, 'Y', 'Z', w, w+dw)
            show_projection_by_slice(short_name, 'Z', 'X', w, w+dw)
            show_projection_by_slice(short_name, 'Z', 'Y', w, w+dw)
            
            # show_projection_by_slice(short_name, 'X', 'Y', w, w+dw, 'test')
            # show_projection_by_slice(short_name, 'Y', 'X', w, w+dw, 'test')

            show_3d_particles()
            # plt.show()
