from __future__ import division
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate
import os


# helpful functions
def Data_pre(t, x, s):
    col = 0
    dat = np.zeros((t, x))
    for i in range(t):
        for j in range(x):
            dat[i, j] = s[j + col]
            # if (j + column) > 20498:
            #    print j + column
        col += x
    return dat

def Open_file(file):

    with open(file) as fn:
        Input_S = fn.readlines()
    Stim = np.asarray(Input_S)

    return Stim

def Gaussian(arr, amp, var, halfsize):
    fullx = 2*halfsize
    for i in range(0, fullx + 0):
        arr[i] = 10 + amp*np.exp(-(i - halfsize)*(i - halfsize)/var**2)

    return arr
# get names of files in the directory
Stimulus = []
Zcells = []
Xcells = []
SXcells = []
Ycells = []
for c in range(2):

    if c == 1:
        cond = 'Attention'
    else:
        cond = 'Neutral'

    isi = 3
    if isi == 1:
        ms = 'isi1_/' #'isi1_/'
    elif isi == 2:
        ms = 'isi2/' #'isi2_/'
    else:
        ms = 'isi3/' #'isi3_/'


    #root = "/home/coglab/Documents/Results/Retina_p3/TemporalExp/" + ms + cond
    root = "/home/coglab/Documents/Results/Retina_p3/TemporalExp3/" + ms + cond
    #root = "/home/coglab/Documents/Results/Retina_p3/TemporalExp4/" + ms + cond


    Stimulus_file = root + "/Simulus.dat"
    print(Stimulus_file)

    Zon_file = root + "/Zcells.dat"
    X_file = root + "/Xcells.dat"
    SX_file = root + "/SXcells.dat"
    Y_file = root + "/Ycells.dat"

    # Open the files
    Stim = Open_file(Stimulus_file)
    Zon = Open_file(Zon_file)
    Xon = Open_file(X_file)
    SXon = Open_file(SX_file)
    Yon = Open_file(Y_file)

    ncells = 580
    time_iter = 1000
    Stimulus.append(Data_pre(time_iter, ncells, Stim))
    Zcells.append(Data_pre(time_iter, ncells, Zon))
    Xcells.append(Data_pre(time_iter, ncells, Xon))
    SXcells.append(Data_pre(time_iter, ncells, SXon))
    Ycells.append(Data_pre(time_iter, ncells, Yon))
    #print Xcells[10, 10], Ycells[500,223]

    # Rectifying
    #Xcells[Xcells < 231.64] = 0.0
    print(np.sum(SXcells[c][100, :]))
    print(np.sum(Xcells[c][100, :]))
    # Save data for further processing
    #path = '/home/coglab/Documents/Results/Retina_p3/TemporalExp/' + ms + cond
    path = '/home/coglab/Documents/Results/Retina_p3/TemporalExp3/' + ms + cond
    #path = '/home/coglab/Documents/Results/Retina_p3/TemporalExp4/' + ms + cond


    np.savetxt(path + '/Mcell.dat', Ycells[c])
    np.savetxt(path + '/Pcell.dat', SXcells[c])


#print type(Ycells)
print(Ycells[0][150, 257])

ncells = 580
time_iter = 1000
# 3D plot
t = np.linspace(0, 250, time_iter)  # np.linspace(0, 3, 2000)   np.linspace(0, 100, 400)
x = np.linspace(0, ncells - 1, ncells)

# Rectify the outputs
#SXcells[SXcells < 231.64] = 0.0
varc = (6*0.03*60.0*60.0)/23.0
halfx = int(ncells/2)
xspace = np.zeros(ncells)
GausRF = Gaussian(xspace, 1, varc, halfx)

plt.figure(1)
plt.plot(GausRF)
# Plots
plt.figure(2)
plt.plot(t, Stimulus[0][:, 290])
plt.title('Input Stimulus')
plt.xlabel('Time')
plt.ylabel('Activity')
plt.figure(3)
plt.plot(Stimulus[0][100, :])
plt.plot(GausRF)
plt.title('Input Stimulus')
plt.xlabel('Space')
plt.ylabel('Activity')
plt.figure(4)
plt.plot(t, Zcells[0][:, 290])
plt.figure(5)
plt.plot(t, Xcells[0][:, 290])
plt.plot(t, Xcells[1][:, 290])
plt.xlabel('Time')
plt.figure(6)
plt.plot(SXcells[0][100, :])
plt.xlabel('Space')
plt.figure(7)
dif = abs(SXcells[0][0, 290] - SXcells[1][0, 290])
plt.plot(t, SXcells[0][:, 290] - 0*dif)
plt.plot(t, 1*SXcells[1][:, 290] - 1*dif, 'r')
plt.title('Sustained response')
plt.xlabel('Time')
plt.ylabel('Activity')
y_values = np.interp(50, t, SXcells[0][:, 290]) # get the value at time 50
y_max = np.interp(29, t, SXcells[0][:, 290]) # get the value at time 50
y_diff = y_max - y_values
print 'y_value at t=50: ', y_values
print 'y_value at t=29: ', y_max
print 'y diff: ', y_diff

plt.figure(8)
plt.plot(SXcells[1][100, :])
plt.xlabel('Space')
plt.figure(9)
plt.plot(t, Ycells[1][:, 249])
#plt.plot(t, Ycells[0][:, 249])
plt.title('Transient response')
plt.xlabel('Time')
plt.ylabel('Activity')
plt.figure(10)
plt.plot(Ycells[1][4, :])
plt.xlabel('Space')
#plt.show()

# 3D plot
X = range(0, ncells)
Y = range(0, time_iter)
X, Y = np.meshgrid(X, Y)
fig = plt.figure(11)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Ycells[1], rstride=5, cstride=5, cmap='cool', linewidth=0, antialiased=False)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Activity')
#ax.set_zlim(-6, 3)
ax.grid(False)
fig.colorbar(surf, shrink=0.5, aspect=5)
fig = plt.figure(12)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Xcells[1], rstride=5, cstride=5, cmap='cool', linewidth=0, antialiased=False)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Activity')
ax.grid(False)
#ax.set_zlim(-1.01, 7.01)
fig.colorbar(surf, shrink=0.5, aspect=5)
fig = plt.figure(13)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, SXcells[1], rstride=5, cstride=5, cmap='cool', linewidth=0, antialiased=False)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Activity')
ax.grid(False)
#ax.set_zlim(-1.01, 7.01)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

