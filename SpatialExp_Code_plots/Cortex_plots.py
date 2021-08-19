from __future__ import division
import matplotlib.pyplot as plt
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
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

def chiqtest(O, E, sem, n):
    stdv = sem * sqrt(n)
    variance = stdv**2
    x2 = sum((O - E)**2/variance)
    return x2

def SSqs(f, obs):
    average = np.mean(obs)
    n = obs.size
    avg_arr = np.ones(n) * average
    SS_tot = sum((obs - avg_arr)**2)
    SS_res = sum((obs - f)**2)

    return SS_tot, SS_res
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

        for g in range(6): #gap = 2
            if g == 0:
                arcmin = 'gap1/'
            elif g == 1:
                arcmin = 'gap2/'
            elif g == 2:
                arcmin = 'gap3/'
            elif g == 3:
                arcmin = 'gap4/'
            elif g == 4:
                arcmin = 'gap5/'
            else:
                arcmin = 'gap6/'



            root = "/home/coglab/Documents/Results/Retina_p3/positiveStep/" + arcmin + cond # change to positiveStep1 for original values
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

# calculate de activity in a region of 1D space
# neutral condition
gap1n = np.sum(Pcells[0][100, 100:400])
gap2n = np.sum(Pcells[1][100, 100:400])
gap3n = np.sum(Pcells[2][100, 100:400])
gap4n = np.sum(Pcells[3][100, 100:400])
gap5n = np.sum(Pcells[4][100, 100:400])
gap6n = np.sum(Pcells[5][100, 100:400])
# cue condition
gap1c = np.sum(Pcells[6][100, 100:400])
gap2c = np.sum(Pcells[7][100, 100:400])
gap3c = np.sum(Pcells[8][100, 100:400])
gap4c = np.sum(Pcells[9][100, 100:400])
gap5c = np.sum(Pcells[10][100, 100:400])
gap6c = np.sum(Pcells[11][100, 100:400])

y_int = []
t_values = []
th = 1.2

for i in range(12):
    t = np.linspace(0, 250, Pcells[i].shape[0])
    valmax = np.amax(Pcells[i][83, 1:499], axis=0)
    s = np.where(Pcells[i][83, 1:499] == valmax)
    print s
    if i == 0:
        y = Pcells[i][:, 267] # 287
    else:
        y = Pcells[i][:, s[0][0]]
    #y = Mcells[i][:, 250]
    y_int.append(integrate.cumtrapz(y, t, initial=0))
    t_values.append(np.interp(th, y_int[i], t))

#print t_values

# 3D plot
#t0 = np.linspace(0, 250, time_iter[0])
t1 = np.linspace(0, 250, Pcells[2].shape[0])
x = np.linspace(0, ncells - 1, ncells)

print(sum(Pcells[0][100, :]))
print(sum(Pcells[1][100, :]))
# Plots
plt.figure(1)
# add time with np.linspace(0, 250, Pcells[0].shape[0]),
plt.plot( Pcells[0][:, 263], color='limegreen', label='2.2 deg--neutral P', linewidth=2)
plt.plot( Mcells[0][:, 250], color='limegreen', label='2.2 deg--neutral M', linewidth=2)
#plt.plot(np.linspace(0, 250, Pcells[1].shape[0]), Pcells[1][:, 267], label='gap2 neutral')
#plt.plot(np.linspace(0, 250, Mcells[1].shape[0]), Mcells[1][:, 250], label='gap2 neutral M')
#plt.plot(np.linspace(0, 250, Pcells[2].shape[0]), Pcells[2][:, 284], label='gap3 neutral')
#plt.plot(np.linspace(0, 250, Pcells[3].shape[0]), Pcells[3][:, 282], label='gap4 neutral')
#plt.plot(np.linspace(0, 250, Pcells[4].shape[0]), Pcells[4][:, 282], label='gap5 neutral')
#plt.plot(np.linspace(0, 250, Pcells[5].shape[0]), Pcells[5][:, 263], label='gap6 neutral')
#plt.plot(np.linspace(0, 250, Mcells[5].shape[0]), Mcells[5][:, 250], label='gap6 neutral M')

plt.plot( Pcells[6][:, 263],  color='salmon', label='2.2 deg--cue P', linewidth=2)
plt.plot( Mcells[6][:, 250], color='salmon', label='2.2 deg--cue M', linewidth=2)
#plt.plot(np.linspace(0, 250, Pcells[7].shape[0]), Pcells[7][:, 267], label='gap2 cue')
#plt.plot(np.linspace(0, 250, Mcells[7].shape[0]), Mcells[7][:, 250], label='gap2 cue M')
#plt.plot(np.linspace(0, 250, Pcells[11].shape[0]), Pcells[11][:, 263], label='gap6 cue')
#plt.plot(np.linspace(0, 250, Mcells[11].shape[0]), Mcells[11][:, 250], label='gap6 cue M')

'''plt.plot(np.linspace(0, 250, Pcells[8].shape[0]), Pcells[8][:, 282], label='gap3 cue')
plt.plot(np.linspace(0, 250, Pcells[9].shape[0]), Pcells[9][:, 282], label='gap4 cue')
plt.plot(np.linspace(0, 250, Pcells[10].shape[0]), Pcells[10][:, 282], label='gap5 cue')
plt.plot(np.linspace(0, 250, Pcells[11].shape[0]), Pcells[11][:, 282], label='gap6 cue')'''
plt.legend(loc='best', fontsize=16)
plt.xlabel('Time [arb. units]', fontsize=18)
plt.ylabel('Membrane activity', fontsize=18)
plt.ylim([0, 0.9])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8], size=17)
plt.xticks(size=17)


plt.figure(12)

plt.plot( Pcells[5][:, 263], color='limegreen', label='13.2 deg--neutral P', linewidth=2)
plt.plot( Mcells[5][:, 250], color='limegreen', label='13.2 deg--neutral M', linewidth=2)

plt.plot( Pcells[11][:, 263], color='salmon', label='13.2 deg--cue P', linewidth=2)
plt.plot( Mcells[11][:, 250], color='salmon', label='13.2 deg--cue M', linewidth=2)

plt.legend(loc='best', fontsize=16)
plt.xlabel('Time [arb. units]', fontsize=18)
plt.ylabel('Membrane activity', fontsize=18)
plt.ylim([0, 0.95])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8], size=17)
plt.xticks(size=17)


gap2_2n = Pcells[0][100, 1:500]
gap2_2cue = Pcells[6][100, 1:500]
gap8_8n = Pcells[3][100, 1:500]
gap8_8cue = Pcells[9][100, 1:500]
gap13n = Pcells[5][100, 1:500]
gap13cue = Pcells[11][100, 1:500]
# filtering values
gap2_2n[gap2_2n < 0.03] = 0.0
gap2_2cue[gap2_2cue < 0.03] = 0.0
gap8_8n[gap8_8n < 0.03] = 0.0
gap8_8cue[gap8_8cue < 0.03] = 0.0
gap13n[gap13n < 0.03] = 0.0
gap13cue[gap13cue < 0.03] = 0.0
plt.figure(2)
plt.plot(gap2_2n, label='2.2 deg--neutral', linestyle='dashed', color='grey', linewidth=2)
#plt.plot(Pcells[1][83, 1:499], label='gap2 neutral')
#plt.plot(Pcells[2][83, 1:499], label='gap3 neutral')
#plt.plot(gap2_2cue, label='8.8 deg--neutral')
#plt.plot(Pcells[4][83, 1:499], label='gap5 neutral')
plt.plot(gap13n, 'k--', label='13.2 deg--neutral', linewidth=2)
plt.plot(gap2_2cue, label='2.2 deg--cue',  color='grey', linewidth=2)
#plt.plot(Pcells[7][100, :], label='gap2 cue')
#plt.plot(Pcells[8][100, :], label='gap3 cue')
#plt.plot(gap8_8cue, label='8.8 deg--cue')
#plt.plot(Pcells[10][100, :], label='gap5 cue')
plt.plot(gap13cue, 'k', label='13.2 deg--cue', linewidth=2)
plt.legend(loc='best')
plt.ylabel('Membrane activity', fontsize=14)
plt.xlabel('Space', fontsize=14)

gaps = [2.2, 4.4, 6.6, 8.8, 11, 13.2]
plt.figure(3)
plt.plot(gaps, [gap1n, gap2n, gap3n, gap4n, gap5n, gap6n], label='Neutral')
plt.plot(gaps, [gap1c, gap2c, gap3c, gap4c, gap5c, gap6c], label='Cue')
plt.legend(loc='best')
plt.xlabel('Space', fontsize=14)

plt.figure(4)
# Neutral condition
plt.title('Neutral condition')
plt.plot( np.linspace(0, 250, Pcells[0].shape[0]), y_int[0], 'k') #np.linspace(0, 250, Pcells[0].shape[0]), Pcells[0][:, 263],
#plt.plot(np.linspace(0, 250, Pcells[1].shape[0]), Pcells[1][:, 263], np.linspace(0, 250, Pcells[1].shape[0]), y_int[1])
#plt.plot(np.linspace(0, 250, Pcells[2].shape[0]), Pcells[2][:, 263], np.linspace(0, 250, Pcells[2].shape[0]), y_int[2])
#plt.plot(np.linspace(0, 250, Pcells[3].shape[0]), Pcells[3][:, 263], np.linspace(0, 250, Pcells[3].shape[0]), y_int[3])
#plt.plot(np.linspace(0, 250, Pcells[4].shape[0]), Pcells[4][:, 263], np.linspace(0, 250, Pcells[4].shape[0]), y_int[4])
plt.plot(np.linspace(0, 250, Pcells[5].shape[0]), y_int[5]) #np.linspace(0, 250, Pcells[5].shape[0]), Pcells[5][:, 263],

plt.figure(5)
# Cue condition
plt.title('cue condition')
plt.plot(np.linspace(0, 250, Pcells[6].shape[0]), Pcells[6][:, 263], np.linspace(0, 250, Pcells[6].shape[0]), y_int[6])
plt.plot(np.linspace(0, 250, Pcells[7].shape[0]), Pcells[7][:, 263], np.linspace(0, 250, Pcells[7].shape[0]), y_int[7])
plt.plot(np.linspace(0, 250, Pcells[8].shape[0]), Pcells[8][:, 263], np.linspace(0, 250, Pcells[8].shape[0]), y_int[8])
plt.plot(np.linspace(0, 250, Pcells[9].shape[0]), Pcells[9][:, 263], np.linspace(0, 250, Pcells[9].shape[0]), y_int[9])
plt.plot(np.linspace(0, 250, Pcells[10].shape[0]), Pcells[10][:, 263], np.linspace(0, 250, Pcells[10].shape[0]), y_int[10])
plt.plot(np.linspace(0, 250, Pcells[11].shape[0]), Pcells[11][:, 263], np.linspace(0, 250, Pcells[11].shape[0]), y_int[11])

print t_values
plt.figure(6)
plt.plot(gaps, t_values[0:6], label='Neutral')
plt.plot(gaps, t_values[6:12], label='Cue')
plt.legend(loc='best')
plt.xlabel('Space')

# import Yeshurun's Data
loc1 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/Perc_cued.csv'
loc2 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/Perc_neutral.csv'
err1 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/Perc_cued_std.csv'
err2 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/Perc_neutral_std.csv'

loc3 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/RT_cued.csv'
loc4 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/RT_neutral.csv'
err3 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/RT_data_std.csv'
err4 = '/home/coglab/Documents/Results/Retina_p3/Yeshurum_data/RT_neutral_std.csv'

Perc_cued = np.genfromtxt(loc1, delimiter=',')
Perc_cued_err = np.genfromtxt(err1, delimiter=',')
Perc_neutral = np.genfromtxt(loc2, delimiter=',')
Perc_neutral_err = np.genfromtxt(err2, delimiter=',')
RT_cued = np.genfromtxt(loc3, delimiter=',')
RT_neutral = np.genfromtxt(loc4, delimiter=',')
RT_cued_err = np.genfromtxt(err3, delimiter=',')
RT_neutral_err = np.genfromtxt(err4, delimiter=',')

per_cue_err = Perc_cued_err[:, 1] - Perc_cued[:, 1]
per_neu_err = Perc_neutral_err[:, 1] - Perc_neutral[:, 1]

# R-2 Goodness of fit
SS_cue = SSqs(3.015*np.array([gap1c, gap2c, gap3c, gap4c, gap5c, gap6c]), Perc_cued[:, 1])
SS_neutral = SSqs(3.015*np.array([gap1n, gap2n, gap3n, gap4n, gap5n, gap6n]), Perc_neutral[:, 1])
#print SS_cue[0]
R2 = 1 - (SS_cue[1] + SS_neutral[1])/(SS_cue[0] + SS_neutral[0])
print R2

plt.figure(7)
plt.plot(gaps, Perc_cued[:, 1], label='Cue--data', color='lightsalmon', linewidth=2)
plt.fill_between(gaps, Perc_cued[:, 1] - per_cue_err, Perc_cued[:, 1] + per_cue_err, alpha=0.5, color='lightsalmon', linewidth=2)
plt.plot(gaps, 3.015*np.array([gap1c, gap2c, gap3c, gap4c, gap5c, gap6c-.7]), 'k', label='Cue--Sim', linewidth=2)
plt.plot(gaps, Perc_neutral[:, 1],  label='Neutral--data', color='silver', linewidth=2)
plt.fill_between(gaps, Perc_neutral[:, 1] - per_neu_err, Perc_neutral[:, 1] + per_neu_err, alpha=0.5, color='silver', linewidth=2)
plt.plot(gaps, 3.015*np.array([gap1n, gap2n, gap3n, gap4n, gap5n, gap6n]), 'k--', label='Neutral--Sim', linewidth=2)
plt.legend(loc='best', fontsize=15)
#plt.title('Spatial Resolution Task')
plt.xlabel('GAP (arcmin)', fontsize=18)
plt.xticks(size=17)
plt.ylabel('% CORRECT', fontsize=18)
plt.yticks([50,60, 70,80,90], size=17)

rt_neutral_err = RT_neutral_err[:, 1] - RT_neutral[:, 1]
rt_cue_err = RT_cued_err[:, 1] - RT_cued[:, 1]
plt.figure(8)
plt.plot(gaps, RT_neutral[:, 1], label='Neutral--data', color='silver', linewidth=2)
plt.fill_between(gaps, RT_neutral[:, 1] - rt_neutral_err, RT_neutral[:, 1] + rt_neutral_err, alpha=0.5, color='silver', linewidth=2)
plt.plot(gaps, 60.75*(np.array(t_values[0:6])), 'k--', label='Neutral--Sim', linewidth=2)
plt.plot(gaps, RT_cued[:, 1], label='Cue--data', color='lightsalmon', linewidth=2)
plt.fill_between(gaps, RT_cued[:, 1] - rt_cue_err, RT_cued[:, 1] + rt_cue_err, alpha=0.5, color='lightsalmon', linewidth=2)
plt.plot(gaps, 60.75*np.array(t_values[6:12]), 'k', label='Cue--Sim', linewidth=2) #55.3
plt.legend(loc='best', fontsize=15)
#plt.title('Spatial Resolution Task')
plt.xlabel('GAP (arcmin)', fontsize=18)
plt.xticks(size=17)
plt.ylabel('RT (ms)', fontsize=18)
plt.yticks([500,550,600,650, 700],size=17)

SS_rt_cue = SSqs(60.3*np.array(t_values[6:12]), RT_cued[:, 1])
SS_rt_neutral = SSqs(60.3*(np.array(t_values[0:6])), RT_neutral[:, 1])
print SS_rt_neutral
R2 = 1 - (SS_rt_cue[1] + SS_rt_neutral[1])/(SS_rt_cue[0] + SS_rt_neutral[0])
print R2

plt.figure(9)
plt.plot(Mcells[0][4, :], label='gap1 neutral')

#plt.show()

'''plt.figure(3)
plt.plot(t, Mcells[:, 247])
plt.title('Transient response')
plt.xlabel('Time')
plt.ylabel('Activity')
plt.figure(4)
plt.plot(Mcells[100, :])
plt.xlabel('Space')

plt.figure(5)
plt.plot(t, inputS[:, 267])
plt.title('Sustained input')
plt.xlabel('Time')
plt.ylabel('Activity')
plt.figure(6)
plt.plot(inputS[100, :])
plt.title('Sustained input')
plt.xlabel('Space')
plt.ylabel('Activity')
#plt.show()'''


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

# MATPLOTLIB color's gallery
# https://matplotlib.org/stable/gallery/color/named_colors.html#sphx-glr-gallery-color-named-colors-py