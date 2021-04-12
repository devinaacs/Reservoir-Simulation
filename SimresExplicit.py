import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

#Define Length and Boundary Conditions
L = 1000
LeftPressure = 1000
RightPressure = 0
Pinit = 200

#Define Transmissibility
por = 0.1
visc = 0.8
compress = 0.000005
perm = 0.05
alpha = por*visc*compress/perm #Transmissibility
print(alpha)

#Define The Explicit Equation
def Explicit(P2,P1,P0,Coeff):
    '''P(t+dt) = P(x) + dt/dx^2 * alpha * (P(x+dx) - 2P(x) + P(x-dx))'''
    return P1+Coeff*(P2-2*P1+P0)

dx      = 100    # Define The Number of Positions to be Measured
t_final = 1      # final timestep
dt      = 0.01   # sec, time interval

Coeff = dt/dx**2/alpha

# This checks the stability.
stability = 4*dt/alpha/dx**2
print(stability)

#Creating X centroid position
x = []
for i in range (0,int(L/dx+2)):
    if i == 0:
        x.append(0)
    else:
        if i == int(L/dx+1):
            x.append(L)
        else:
            pos_init = dx/2 + dx*(i-1)
            x.append(pos_init)

n = len(x)-1
t = np.arange(0,t_final+dt*2,dt)      #

print(x)

P = np.zeros((len(t),len(x)))
# This creates an empty array that have the size of [r] into the memory.

# This sets the initial condition of the initial pressure
P[0,0] = LeftPressure
P[0,n] = RightPressure
for i in range(1,n):
    P[0,i] = Pinit
print(P[0,:])

# This is the explicit equation
for j in range(1, len(t)):
    P[j,0] = LeftPressure
    P[j,n] = RightPressure
    P[j,1] = P[j-1,1] + Coeff/0.75 * (P[j-1,2]-3*P[j-1,1]+2*P[j-1,0])
    for i in range(2,n-1):
        P[j,i] = Explicit(P[j-1,i+1],P[j-1,i],P[j-1,i-1],Coeff)
    P[j,n-1] = P[j-1,n-1] + Coeff/0.75 * (2*P[j-1,n]-3*P[j-1,n-1]+P[j-1,n-2])
    # P[j,n-1] = Explicit(P[j-1,n],P[j-1,n-1],P[j-1,n-2],Coeff/0.75)
    # print(P[j,:])

#The Code Below defines the plot of the simulation.
fig, ax = plt.subplots()                        # Defines an empty plot window
plt.subplots_adjust(bottom=0.25)                # Adjusts the position of the plot area relative to the window edge
l, = plt.plot(x,P[0,:])                         # Defines the data to be plotted
ax.set_title('Time (sec): %.3f' %(0))           # Set Title Name. The %(var) refers to the variable to be printed in the title. 
                                                # %.2f defines the print decimal limit
ax.set_xlim([0,L])                              # This sets the limit of where the axis extends to in the X direction                              
ax.set_ylim([0,LeftPressure])                   # This sets the limit of where the axis extends to in the Y direction
ax.set_xlabel('Distance (m)')                   # Set the label of the X axis
ax.set_ylabel('Pressure (psia)')                # Set the label of the Y axis

axtime = plt.axes([0.2, 0.1, 0.65, 0.03])                                           # Set the [Left Indent, Bottom Indent, Length, Thickness]
time_slider = Slider(axtime, 'Time', 0, t_final, valinit=0, valstep=dt)             # Defines the slider and what variable is being plotted.

# axalpha = plt.axes([0.2, 0.05, 0.65, 0.03])                                               # Set the [Left Indent, Bottom Indent, Length, Thickness]
# alpha_slider = Slider(axalpha, 'Transmissibility', 0, 100, valinit=alpha, valstep=1)      # Defines the slider and what variable is being plotted.

# This update the value change when the slider changes it value
def update(val):
    new_time = int(time_slider.val/dt)
    # new_alpha = alpha_slider.val
    ax.set_title('Time (sec): %.3f' %(time_slider.val))
    A = P[new_time,:]
    l.set_ydata(A)
    fig.canvas.draw_idle()

print(P[1,:])

# Button to automatically play the plot
resetax = plt.axes([0.75, 0.15, 0.1, 0.04])
button = Button(resetax, 'Play', color='lightgoldenrodyellow',hovercolor='0.975')

# # This defines how the plot play overtime when the button play is pressed
# # During this time, no time value will update when you click it until the function ends.
def auto(val):
    for j in range(0, len(t)-1):
        ax.set_title('Time (sec): %.2f' %(t[j]))
        time_slider.set_val(t[j])
        A = P[j,:]
        l.set_ydata(A)
        fig.canvas.draw()

time_slider.on_changed(update)
# alpha_slider.on_changed(update)
button.on_clicked(auto)
plt.show()