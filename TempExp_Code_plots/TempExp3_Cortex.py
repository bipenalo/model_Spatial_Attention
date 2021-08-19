from __future__ import division
import matplotlib.pyplot as plt
import matplotlib
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os

def SSqs(f, obs):
    average = np.mean(obs)
    n = len(obs)
    avg_arr = np.ones(n)*average
    SS_tot = np.sum((obs - avg_arr)**2)
    SS_res = sum((obs - f)**2)
    return SS_tot, SS_res

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

def getconditions():
    # get names of files in the directory
    inputS = []
    Pcells = []
    Mcells = []
    for c in range(2):
        #c = 1
        if c == 1:
            cond = 'Attention'
        else:
            cond = 'Neutral'

        for g in range(3): #3
            if g == 0:
                ms = 'isi2/' #isi1_/
            elif g == 1:
                ms = 'isi2/'
            else:
                ms = 'isi3/' #isi3_/



            root = "/home/coglab/Documents/Results/Retina_p3/TemporalExp3/" + ms + cond #change to TemporalExp!
            inp_data = root + "/Stimuli.dat"
            SX_file = root + "/cortexP.dat"
            Y_file = root + "/cortexM.dat"

            # Open the files
            Stimulus = Open_file(inp_data)
            Pcortex = Open_file(SX_file)
            Mcortex = Open_file(Y_file)



            ncells = 580
            time_iter_Stimulus = int(Stimulus.shape[0]/ncells)
            time_iter_Pcortex = int(Pcortex.shape[0]/ncells)
            time_iter_Mcortex = int(Mcortex.shape[0]/ncells)


            inputS.append(Data_pre(time_iter_Stimulus, ncells, Stimulus))
            Pcells.append(Data_pre(time_iter_Pcortex, ncells, Pcortex))
            Mcells.append(Data_pre(time_iter_Mcortex, ncells, Mcortex))

    return inputS, Pcells, Mcells, ncells

inputS, Pcells, Mcells, ncells = getconditions()



# Plots



plt.figure(1)
plt.plot(np.linspace(0, 250, Pcells[0].shape[0]), Pcells[0][:, 245], label='isi1 neutral P')
plt.plot(np.linspace(0, 250, Mcells[0].shape[0]), Mcells[0][:, 234], label='isi1 neutral M')

plt.plot(np.linspace(0, 250, Pcells[3].shape[0]), Pcells[3][:, 245], label='isi1 cue P')
plt.plot(np.linspace(0, 250, Mcells[3].shape[0]), Mcells[3][:, 234], label='isi1 cue M')
plt.legend(loc='best')
plt.title('ISI1')
plt.xlabel('Time')
plt.ylabel('Activity')

plt.figure(2)
plt.plot(np.linspace(0, 250, Pcells[1].shape[0]), Pcells[1][:, 267], label='isi2 neutral P')
plt.plot(np.linspace(0, 250, Mcells[1].shape[0]), Mcells[1][:, 234], label='isi2 neutral M') #250

plt.plot(np.linspace(0, 250, Pcells[4].shape[0]), Pcells[4][:, 267], label='isi2 cue P')
plt.plot(np.linspace(0, 250, Mcells[4].shape[0]), Mcells[4][:, 234], label='isi2 cue M') #250
plt.legend(loc='best')
plt.title('ISI2')
plt.xlabel('Time')
plt.ylabel('Activity')

plt.figure(3)
plt.plot(np.linspace(0, 250, Pcells[2].shape[0]), Pcells[2][:, 245], label='isi3 neutral P')
plt.plot(np.linspace(0, 250, Mcells[2].shape[0]), Mcells[2][:, 234], label='isi3 neutral M')

plt.plot(np.linspace(0, 250, Pcells[5].shape[0]), Pcells[5][:, 245], label='isi3 cue P')
plt.plot(np.linspace(0, 250, Mcells[5].shape[0]), Mcells[5][:, 234], label='isi3 cue M')
plt.legend(loc='best')
plt.title('ISI3')
plt.xlabel('Time')
plt.ylabel('Activity')


# import Yeshurun's Data
loc1 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/temporal_cued.csv'
loc2 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/temporal_neutral.csv'
err1 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/temporal_cued_std.csv'
err2 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/temporal_neutral_std.csv'

isis = [11, 23, 35]
temporal_cued = np.genfromtxt(loc1, delimiter=',')
temporal_neutral = np.genfromtxt(loc2, delimiter=',')
temporal_cued_err = np.genfromtxt(err1, delimiter=',')
temporal_neutral_err = np.genfromtxt(err2, delimiter=',')

temp_cue_err = temporal_cued_err[:, 1] - temporal_cued[:, 1]
temp_neutral_err = temporal_neutral_err[:, 1] - temporal_neutral[:, 1]

# get the maximum value
'''isi1n = np.max(np.array(Mcells[0][120:600, 250]), axis=0)
isi2n = np.amax(Mcells[1][120:600, 250], axis=0)
isi3n = np.amax(Mcells[2][120:600, 250], axis=0)
isi1c = np.amax(Mcells[3][120:600, 250], axis=0)
isi2c = np.amax(Mcells[4][120:600, 250], axis=0)
isi3c = np.amax(Mcells[5][120:600, 250], axis=0)'''
# sum of area in time
'''isi1n = np.sum(Mcells[0][120:600, 250])
isi2n = np.sum(Mcells[1][120:600, 250])
isi3n = np.sum(Mcells[2][120:600, 250])
isi1c = np.sum(Mcells[3][120:600, 250])
isi2c = np.sum(Mcells[4][120:600, 250])
isi3c = np.sum(Mcells[5][120:600, 250])'''
for i in range(6):
    valmax = np.amax(Mcells[i][200:600, 234], axis=0)
    s = np.where(Mcells[i][200:600, 234] == valmax)
    print s
# sum of area in space
isi1n = np.sum(Mcells[0][232, 220:280])
isi2n = np.sum(Mcells[1][275, 200:280])
isi3n = np.sum(Mcells[2][325, 200:300])
isi1c = np.sum(Mcells[3][232, 220:280])
isi2c = np.sum(Mcells[4][275, 200:280])
isi3c = np.sum(Mcells[5][323, 200:300])
#plt.plot(Mcells[1][:, 234])


plt.figure(9)
scale = 0.085
isis1 = [23, 34]
plt.plot(isis1, [1.55, 2.26], label='Cue--data', color='lightsalmon', linewidth=2)
cue_err = np.array([1.67, 2.39]) - np.array([1.55, 2.26])
plt.fill_between(isis1, np.array([1.55, 2.26]) - cue_err, np.array([1.55, 2.26]) + cue_err, alpha=0.5, color='lightsalmon')
#plt.plot(isis1, scale*np.array([0.7*69.44*0.37, 0.7*79.63*0.51]), 'k', label='Cue--Sim', linewidth=2)
plt.plot(isis1, [1.81, 2.88], label='Neutral--data', color='silver', linewidth=2)
neutral_err = np.array([1.92, 2.99]) - np.array([1.81, 2.88])
plt.fill_between(isis1, np.array([1.81, 2.88]) - neutral_err, np.array([1.81, 2.88]) + neutral_err, alpha=0.5, color='silver')
#plt.plot(isis1, scale*np.array([isi2n, isi3n]), 'k--', label='Neutral--Sim', linewidth=2)
#plt.plot(isis1, scale*np.array([0.7*69.44*0.44, 0.7*80.89*0.58]), 'k--', label='Neutral--Sim', linewidth=2)
plt.legend(loc='best', fontsize=15)
#plt.title('Temporal Resolution Task')
plt.xlabel('ISI (ms)', fontsize=18)
plt.xticks([23,29,35], size=17)
plt.ylabel("d'", fontsize=18)
plt.yticks([1.5, 2.0, 2.5,3.0], size=17)

# R2 operation
R2n = SSqs(scale*np.array([0.7*69.44*0.44, 0.7*80.89*0.58]), np.array([1.81, 2.88]))
R2c = SSqs( scale*np.array([0.7*69.44*0.37, 0.7*79.63*0.515]), np.array([1.55, 2.26]))
R2 = 1 - (R2n[1] + R2c[1])/(R2n[0] + R2c[0])
print R2
# 3D plot
X = range(0, 500)
Y = range(0, Pcells[2].shape[0])
X, Y = np.meshgrid(X, Y)
fig = plt.figure(10)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Pcells[2][:, 0:500], rstride=5, cstride=5, cmap='cool', linewidth=0, antialiased=False)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Activity')
#ax.set_zlim(-1.01, 3.01)
ax.grid(False)
fig.colorbar(surf, shrink=0.5, aspect=5)

X = range(0, ncells)
Y = range(0, Mcells[2].shape[0])
X, Y = np.meshgrid(X, Y)
fig = plt.figure(11)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Mcells[2], rstride=5, cstride=5, cmap='cool', linewidth=0, antialiased=False)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Activity')
ax.set_zlim(0.0, 0.9)
ax.grid(False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()