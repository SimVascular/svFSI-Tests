import numpy as np
import matplotlib.pyplot as plt

fname = "./u_y0.5_Re1000.dat"
temp = np.loadtxt(fname,skiprows=1)
y = temp[:,0]
u = temp[:,1]

fname = "../../Ref/data_Ghia/ghia_u_y0.5_Re1000.dat"
temp = np.loadtxt(fname,skiprows=1)
yg = temp[:,0]
ug = temp[:,1]


plt.figure(figsize=(8,6))
plt.title("x-velocity @ x=0.5",fontsize=25)
plt.plot(y, u,linewidth=3,label='Current')
plt.plot(yg, ug,'o',markersize=10,label='Ghia, JCP, 1982')
plt.xlabel('y',fontsize=30)
plt.ylabel('u',fontsize=30)
ax1 = plt.gca()
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2)
plt.yticks(fontsize=25)
plt.xticks(fontsize=25)
plt.legend(fontsize=22)
plt.tight_layout()
plt.show()

#######################################################

fname = "./v_x0.5_Re1000.dat"
temp = np.loadtxt(fname,skiprows=1)
x = temp[:,0]
v = temp[:,1]

fname = "../../Ref/data_Ghia/ghia_v_x0.5_Re1000.dat"
temp = np.loadtxt(fname,skiprows=1)
xg = temp[:,0]
vg = temp[:,1]


plt.figure(figsize=(8,6))
plt.title("y-velocity @ y=0.5",fontsize=25)
plt.plot(x, v,linewidth=3,label='Current')
plt.plot(xg, vg,'o',markersize=10,label='Ghia, JCP, 1982')
plt.xlabel('x',fontsize=30)
plt.ylabel('v',fontsize=30)
ax1 = plt.gca()
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2)
plt.yticks(fontsize=25)
plt.xticks(fontsize=25)
plt.legend(fontsize=22)
plt.tight_layout()
plt.show()


