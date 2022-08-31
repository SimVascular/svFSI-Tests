import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = "Times New Roman"

time = np.linspace(1,1600,1600)*0.005

## Velocity data
temp = np.loadtxt("./16-procs/B_NS_Pressure_average.txt",skiprows=3)
po = temp[:,1]

plt.figure(figsize=(8,6))
plt.plot(time, po/1334.16,linewidth=3)
plt.xlabel('t',fontsize=24)
plt.ylabel('Outlet Pressure [mmHg]',fontsize=24)
ax1 = plt.gca()
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.tight_layout()
plt.show()


