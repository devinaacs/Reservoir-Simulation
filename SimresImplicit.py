import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from numpy.random import random

LengthX = 2000                               #ft
LengthY = 2000                               #ft
dLx = 250                                    #ft
dLy = 250                                    #ft
Lx = np.arange(dLx/2,LengthX+dLx,dLx)       
Ly = np.arange(dLy/2,LengthY+dLy,dLy)       

t_init = 0
t_final = 100
dt = 1
tplot = np.arange(t_init,t_final+dt,dt)

#======================SIMULATION PROPERTIES==================================================
por = 0.1
comp = 0.000001
visc = 600
perm = 10
fvf = 1.2

nrows = len(Ly)
ncols = len(Lx)
nnum = nrows*ncols
tnum = len(tplot)

O = np.ones((nrows,ncols))
T = np.ones((len(tplot),len(Lx)+2))

for j in range(0,nrows):
    for i in range(0,ncols):
        O[j,i] = por*visc*comp/perm

RHS = np.ones(nnum)
LHS = np.zeros((nnum,nnum))

#======================GROUP CALCULATION===========================================
alpha = 1/dLx**2*perm/visc/fvf
beta = 1/dLy**2*perm/visc/fvf
gamma = por/dt*comp/fvf

#======================ASSIGNING TRANSMISSIBILITY===========================================
for i in range (0,nnum):
    if i in np.arange(0,nnum,ncols):
        if i == 0:                                      # 1,1 boundary
            LHS[i,i]        = (alpha + beta + gamma)/gamma
            LHS[i,i+1]      = -alpha/gamma
            LHS[i,i+ncols]  = -beta/gamma
        elif i == nnum-ncols:                           # 1,N boundary
            LHS[i,i-ncols]  = -beta/gamma
            LHS[i,i]        = (alpha + beta + gamma)/gamma
            LHS[i,i+1]      = -alpha/gamma
        else:                                           # 1,j boundary
            LHS[i,i-ncols]  = -beta/gamma
            LHS[i,i]        = (alpha + 2*beta + gamma)/gamma
            LHS[i,i+1]      = -alpha/gamma
            LHS[i,i+ncols]  = -beta/gamma
    elif i in np.arange(ncols-1,nnum,ncols):
        if i == ncols-1:                                # N,1 boundary
            LHS[i,i-1]      = -alpha/gamma
            LHS[i,i]        = (alpha + beta + gamma)/gamma
            LHS[i,i+ncols]  = -beta/gamma
        elif i == nnum-1:                               # N,N boundary
            LHS[i,i-ncols]  = -beta/gamma
            LHS[i,i-1]      = -alpha/gamma
            LHS[i,i]        = (alpha + beta + gamma)/gamma
        else:                                           # N,j boundary
            LHS[i,i-ncols]  = -beta/gamma
            LHS[i,i-1]      = -alpha/gamma
            LHS[i,i]        = (2*alpha + beta + gamma)/gamma
            LHS[i,i+ncols]  = -beta/gamma
    elif i < ncols:                                     # i,1 boundary
        LHS[i,i-1]          = -alpha/gamma
        LHS[i,i]            = (2*alpha + beta + gamma)/gamma
        LHS[i,i+1]          = -alpha/gamma
        LHS[i,i+ncols]      = -beta/gamma
    elif i > nnum-ncols:                                # i,N boundary
        LHS[i,i-ncols]      = -beta/gamma
        LHS[i,i-1]          = -alpha/gamma
        LHS[i,i]            = (2*alpha + beta + gamma)/gamma
        LHS[i,i+1]          = -alpha/gamma
    else:                                               # i,j grid
        LHS[i,i-ncols]      = -beta/gamma
        LHS[i,i-1]          = -alpha/gamma
        LHS[i,i]            = (2*alpha + 2*beta + gamma)/gamma
        LHS[i,i+1]          = -alpha/gamma
        LHS[i,i+ncols]      = -beta/gamma

# with np.printoptions(precision=1, suppress=True):
#     print(LHS/1000)


#======================DEFINING WELL===============================================

bhp = 0                           #psia
interface_area = 1000
net_pay = 1                          #ft
well_radius = 3                      #inch
well_grid = 40


well_radius_ft  = well_radius/12
drain_radius = 0.2*np.sqrt(dLx*dLy)
mobility_term = 1/visc/fvf

well_index = 2*np.pi*perm*net_pay/np.log(drain_radius/well_radius_ft)

#==========Calculating Well in Grid
# tau = well_index*mobility_term*(dLx+dLy)/grid_area/(dLx*dLy)
tau = well_index*mobility_term/interface_area/(dLx)
LHS[well_grid-1,well_grid-1] = (LHS[well_grid-1,well_grid-1]*gamma + tau)/gamma

print(LHS)

#======================SETS INITIAL PRESSURE========================================
RHS[:] = 3000

# with np.printoptions(precision=1, suppress=True):
#     print(LHS/1000)

P = np.ones((tnum,nnum))
P2 = np.ones((tnum,nrows,ncols))

for i in range(0,tnum):
    RHS = RHS + bhp*tau
    RHS = np.linalg.solve(LHS,RHS)
    P[i,:] = RHS
    print(RHS)

for i in range(0,tnum):
    for j in range(0,nrows):
        for k in range(0,ncols):
            P2[i][j][k] = P[i][j*ncols+k]

# ======================PLOTTING SYSTEM==========================================
ploty, plotx = np.mgrid[slice(0,LengthY+dLy*2,dLy),
                        slice(0,LengthX+dLx*2,dLx)]

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)                                                    # Adjusts the position of the plot area relative to the window edge

axtime = plt.axes([0.125, 0.1, 0.685, 0.03])                                        # Set the [Left Indent, Bottom Indent, Length, Thickness]
time_slider = Slider(axtime, 'Time', 0, t_final, valinit=t_init, valstep=dt)        # Defines the slider and what variable is being plotted.

ax.set_title('Time (sec): %.3f' %(0))                                               # Set Title Name. The %(var) refers to the variable to be printed in the title. %.2f defines the print decimal limit
ax.set_xlabel('x (ft)')                                                             # Set the label of the X axis
ax.set_ylabel('y (ft)')                                                             # Set the label of the Y axis
cmap = plt.get_cmap('summer')
im = ax.pcolormesh(plotx, ploty, P2[time_slider.val,:,:],vmin=1000,vmax=3200,cmap=cmap)
fig.colorbar(im, ax=ax)

def update(val):
    new_time = int(time_slider.val/dt)
    ax.set_title('Time (sec): %.3f' %(time_slider.val))
    A = P2[new_time]
    im.set_array(A.ravel())
    fig.canvas.draw_idle()

time_slider.on_changed(update)
plt.show()




#===================================ARCHIVE=====================================
# for i in range (0,nnum):
#     if i in np.arange(0,nnum,ncols):
#         if i == 0:
#             LHS[i,i]        = 3
#             LHS[i,i+1]      = 4
#             LHS[i,i+ncols]  = 5
#         elif i == nnum-ncols:
#             LHS[i,i-ncols]  = 1
#             LHS[i,i]        = 3
#             LHS[i,i+1]      = 4
#         else:
#             LHS[i,i-ncols]  = 1
#             LHS[i,i]        = 3
#             LHS[i,i+1]      = 4
#             LHS[i,i+ncols]  = 5
#     elif i in np.arange(ncols-1,nnum,ncols):
#         if i == ncols-1:
#             LHS[i,i-1]      = 2
#             LHS[i,i]        = 3
#             LHS[i,i+ncols]  = 5
#         elif i == nnum-1:
#             LHS[i,i-ncols]  = 1
#             LHS[i,i-1]      = 2
#             LHS[i,i]        = 3
#         else:
#             LHS[i,i-ncols]  = 1
#             LHS[i,i-1]      = 2
#             LHS[i,i]        = 3
#             LHS[i,i+ncols]  = 5
#     elif i < ncols:
#         LHS[i,i-1]          = 2
#         LHS[i,i]            = 3
#         LHS[i,i+1]          = 4
#         LHS[i,i+ncols]      = 5
#     elif i > nnum-ncols:
#         LHS[i,i-ncols]      = 1
#         LHS[i,i-1]          = 2
#         LHS[i,i]            = 3
#         LHS[i,i+1]          = 4
#     else:
#         LHS[i,i-ncols]      = 1
#         LHS[i,i-1]          = 2
#         LHS[i,i]            = 3
#         LHS[i,i+1]          = 4
#         LHS[i,i+ncols]      = 5