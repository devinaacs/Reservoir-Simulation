import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from numpy.random import random

def solver(LHS,RHS):
    n = LHS.shape[1]
    # Forward Elimination (Gauss)
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            factor = LHS[i, k] / LHS[k, k]
            for j in range(k+1, n):
                LHS[i, j] = LHS[i, j] - factor * LHS[k, j]
            RHS[i] = RHS[i] - factor * RHS[k]

    # Back Substitution
    x = np.zeros((n, 1))
    x[n-1] = RHS[n-1] / LHS[n-1, n-1]
    for i in range(n - 2, -1, -1):
        suma = RHS[i]
        for j in range(i + 1, n):
            suma = suma - LHS[i, j] * x[j]
        x[i] = suma / LHS[i, i]
    return x

#======================SIMULATION PROPERTIES==================================================
por = 0.1
mu = 80
comp = 0.000001
perm = 1

#======================LENGTH PROPERTIES======================================================
Length = 10000
dx = 50

#======================TIME PROPERTIES========================================================
t_final = 100
dt = 1

#======================INITIAL AND BOUNDARY PRESSURE==========================================
P_init = 200
PL = 1000
PR = 14.7

#======================TRANSMISSIBILITY CALCULATION===========================================
alpha        = dx**2/dt*por*mu*comp/perm             #Transmmissibility at non-Boundary blocks
beta         = alpha * 3/4                           #Transmissibility at Boundary


#This tells the position in which the program will define the x position
#This also tells how many positions are being measured using the len(x) function
x = np.arange(dx/2,Length+dx/2,dx)
t = np.arange(0,t_final+dt,dt)
xnum = len(x)
tnum = len(t)

LHS = np.zeros((xnum,xnum))
RHS = np.ones(xnum)

#This construct the Transmissbility Value in each pressure
for i in range(0,xnum):
    if i == 0:
        LHS[i,i]        = (3+beta)/beta
        LHS[i,i+1]      = -1/beta
    else:
        if i == xnum-1:
            LHS[i,i-1]  = -1/beta
            LHS[i,i]    = (3+beta)/beta
        else:
            LHS[i,i-1]  = -1/alpha
            LHS[i,i]    = (2+alpha)/alpha
            LHS[i,i+1]  = -1/alpha
print(LHS)

#This sets the initial pressure
RHS = RHS*P_init

#Alternative way to set the initial pressure
# for i in range(0,len(RHS)):
    # RHS[i] = RHS[i]*P_init*random()


fig, ax = plt.subplots()                        # Defines an empty plot window
plt.subplots_adjust(bottom=0.25)                # Adjusts the position of the plot area relative to the window edge


P2 = np.zeros((tnum,xnum+2))
# The P2 will be used to store all arrays that are coming out from the loop function below.

for f in range (0,tnum):
    RHS[0]          = RHS[0] + 2*PL/beta
    RHS[xnum-1]     = RHS[xnum-1] + 2*PR/beta
    # RHS             = solver(LHS,RHS)
    RHS             = np.linalg.solve(LHS,RHS) 
    # RHS             = cp.linalg.solve(LHS,RHS) #GPU-Based Library CuPy
    #  This currently cheats my way through as the current naive gaussian doesn't perform partial pivoting to prevent numerical errors.
    #The above code tells the program to apply the boundary condition first on the first block and the last block by applying:
    #   P1 = P1 + 2*PL/T        where T = Transmissibility value in the boundary
    #   PN = PN + 2*PR/T

    P               = np.ones((xnum+2))
    L               = np.ones((xnum+2))
    P[0]            = PL
    L[0]            = 0

    for g in range(1,xnum+1):
        P[g]            = RHS[g-1]
        L[g]            = x[g-1]
    P[len(P)-1]     = PR
    L[len(L)-1]     = Length
    
    for i in range(0,len(P)):
        P2[f,i]         = P[i]

# print(P2)



l, = plt.plot(L,P2[0])                  # Defines the data to be plotted
ax.set_title('Time (sec): %.3f' %(0))   # Set Title Name. The %(var) refers to the variable to be printed in the title. %.2f defines the print decimal limit
ax.set_xlim([0,Length])                 # This sets the limit of where the axis extends to in the X direction                              
ax.set_ylim([0,PL])                     # This sets the limit of where the axis extends to in the Y direction
ax.set_xlabel('Distance (m)')           # Set the label of the X axis
ax.set_ylabel('Pressure (psia)')        # Set the label of the Y axis

axtime = plt.axes([0.2, 0.1, 0.65, 0.03])                                           # Set the [Left Indent, Bottom Indent, Length, Thickness]
time_slider = Slider(axtime, 'Time', 0, t_final, valinit=0, valstep=dt)             # Defines the slider and what variable is being plotted.

# Button to automatically play the plot
resetax = plt.axes([0.75, 0.15, 0.1, 0.04])
button = Button(resetax, 'Play', color='lightgoldenrodyellow',hovercolor='0.975')

def update(val):
    new_time = int(time_slider.val/dt)
    # new_alpha = alpha_slider.val
    ax.set_title('Time (sec): %.3f' %(time_slider.val))
    A = P2[new_time]
    l.set_ydata(A)
    fig.canvas.draw_idle()

# This defines how the plot play overtime when the button play is pressed
# During this time, no time value will update when you click it until the function ends.
def auto(val):
    for j in range(0, len(t)):
        ax.set_title('Time (sec): %.2f' %(t[j]))
        time_slider.set_val(t[j])
        A = P2[j]
        l.set_ydata(A)
        fig.canvas.draw()

time_slider.on_changed(update)
button.on_clicked(auto)
plt.show()