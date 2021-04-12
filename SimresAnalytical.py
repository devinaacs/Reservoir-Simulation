import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

#Define Length and Boundary Conditions=============================================================================
L = 1000
LeftPressure = 1000
RightPressure = 100

#Define Transmissibility
alpha = 30 #Transmissibility

#Define The Analytical Equation====================================================================================
def Equation(PL,PR,x,Length,IterLimit,Alpha,time):
    '''P(x,t) = PL + (PR-PL)[x/L + 2/pi*SUM(1/2*exp(-(n^2 * pi^2 / L^2 * alpha * t))sin(n*pi*x/L))]'''
    A = (PR-PL)
    B = x/Length
    C = 0
    D = 0
    D = C + D
    for n in range (1, IterLimit): # This Limits the Iteration
        C = np.sin(n*np.pi*x/Length)*np.exp(-(n**2*np.pi**2)/(Length**2)*Alpha*time)/n
        D = C + D
    return PL + (A)*(B + 2*D/np.pi)

r       = 21    # Define The Number of Positions to be Measured. The higher the number, the smoother the lines.
t_final = 5000  # final time
dt      = 50    # sec, time interval.
redo    = 300   # number of n to be iterated in Equation function

#CREATING DATA STRUCTURES==========================================================================================
P=np.zeros(r)   # This creates an array filled with zeros with [r] number of points into the memory.

x = np.linspace(0,L,r)          #Creates an array of evenly-spaced distance with [r] number of distances
t = np.arange(0,t_final+dt,dt)  #Creates an array of forward-spacing value from 0 until before t_final+dt with step dt

# Differences between np.linspace and np.arange:
#
#       np.linspace creates a determined amount of elements
#
#               np.linspace(0,50,6) = [0, 10, 20, 30, 40, 50]                       <<< 6 elements
#               np.linspace(0,50,5) = [0, 12.5, 25, 37.5, 50]                       <<< 5 elements
#               np.linspace(0,50,7) = [0, 8.333 16.667, 25, 33.333, 41.667, 50]     <<< 7 elements
#
#
#       np.arange keeps adding the step value until it reaches the final value.
#       note that the final value is not written because it's how the function works.
#
#               np.arange(0,50,5)   = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]        <<< 10 elements
#
#       np.arange function write everything within timestep under 50 (means the final writeable value is 49.99999)
#       to add the final value you needed, you can add its own step variable to increase its write range
#
#               dt = 5
#               np.arange(0,50+dt,dt) = [0, 50, 10, 15, 20, 25, 30, 35, 40, 45, 50] <<< 11 elements

#WRITING PRESSURE VALUE FROM CALCULATION TO EACH X POSITION========================================================
curtime = dt*0
for i in range(0,r):
    P[i] = Equation(LeftPressure, RightPressure, x[i], L, redo, alpha, curtime)

#The Code Below defines the plot of the simulation.================================================================
fig, ax = plt.subplots()                        # Defines an empty plot window
plt.subplots_adjust(bottom=0.25)                # Adjusts the position of the plot area relative to the window edge
l, = plt.plot(x,P)                              # Defines the data to be plotted
ax.set_title('Time (sec): %.2f' %(curtime))     # Set Title Name. The %(var) refers to the variable to be printed in the title. %.2f defines the print decimal limit
ax.set_xlim([0,L])                              # This sets the limit of where the axis extends to in the X direction                              
ax.set_ylim([RightPressure,LeftPressure])       # This sets the limit of where the axis extends to in the Y direction
ax.set_xlabel('Distance (m)')                   # Set the label of the X axis
ax.set_ylabel('Pressure (psia)')                # Set the label of the Y axis

axtime = plt.axes([0.2, 0.1, 0.65, 0.03])                                           # Set the [Left Indent, Bottom Indent, Length, Thickness]
time_slider = Slider(axtime, 'Time', 0, t_final, valinit=curtime, valstep=dt)       # Defines the slider and what variable is being plotted.

axalpha = plt.axes([0.2, 0.05, 0.65, 0.03])                                           # Set the [Left Indent, Bottom Indent, Length, Thickness]
alpha_slider = Slider(axalpha, 'Transmissibility', 0, 100, valinit=alpha, valstep=1)       # Defines the slider and what variable is being plotted.

# This update the value change when the slider changes it value====================================================
def update(val):
    new_time = time_slider.val
    new_alpha = alpha_slider.val
    ax.set_title('Time (sec): %.2f' %(new_time))
    for i in range(0,r):
        P[i] = Equation(LeftPressure, RightPressure, x[i], L, redo, new_alpha, new_time)
    l.set_ydata(P)
    fig.canvas.draw_idle()


# Button to automatically play the plot============================================================================
resetax = plt.axes([0.75, 0.15, 0.1, 0.04])
button = Button(resetax, 'Play', color='lightgoldenrodyellow',hovercolor='0.975')

# This defines how the plot play overtime when the button play is pressed.=========================================
# During this time, no alpha and time value will update when you click it until the function ends.=================
def auto(val):
    new_alpha = alpha_slider.val
    for j in range(0, len(t)):
        ax.set_title('Time (sec): %.2f' %(t[j]))
        time_slider.set_val(t[j])
        for i in range(0,r):
            P[i] = Equation(LeftPressure, RightPressure, x[i], L, redo, new_alpha, t[j])
        l.set_ydata(P)
        fig.canvas.draw()

time_slider.on_changed(update)
alpha_slider.on_changed(update)
button.on_clicked(auto)
plt.show()

# ============================TEST SITE. DOESN'T MATTER ANYMORE========================================

# plt.ion()
# fig, ax = plt.subplots()
# for j in range(0, len(t)):
#     for i in range(0, r):
#         P[i]=Equation(LeftPressure, RightPressure, x[i], L, 600, alpha, t[j])
#     # print(P)
#     plt.cla()
#     ax.plot(x,P)
#     ax.set_title('Time (sec): %.2f' %(t[j]))
#     ax.set_xlim([0, L])
#     ax.set_ylim([0,LeftPressure])
#     ax.set_xlabel('Distance (m)')
#     ax.set_ylabel('Pressure (psia)')
#     fig.canvas.draw()

# print(x)

# print(dx)
# print(Equation(LeftPressure, RightPressure, x[1], L, 100, alpha, 10000))