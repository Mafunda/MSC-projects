import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline

time = np.array([660.1,1320.2,1980.3,2640.1,3300.5,3960.6,4620.7,5280.8,5940.9,6601])
a = np.array([1.34e-2,2.00e-3,1.25e-2,3.23e-3,1.52e-3,4.86e-2,1.99e-2,2.25e-1,7.59e-1,9.64e-2])
b = np.array([7.20e-4,8.89e-4,9.53e-4,1.09e-3,1.15e-3,1.18e-3,5.83e-3,1.70e-1,1.71e-1,1.71e-1])

time_smooth = np.linspace(time.min(), time.max(), 100000)


a_smooth = spline(time,a,time_smooth)
b_smooth = spline(time,b,time_smooth)

plt.plot(time_smooth,a_smooth,'r', label = "Inflow")
plt.plot(time_smooth,b_smooth,'k', label = "Cumm. Bottom flux")
plt.legend(loc="best")
plt.xlabel("Time(hr)")
plt.ylabel("Volume(mm/hr)")
plt.show()
