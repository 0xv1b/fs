import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl


mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['legend.title_fontsize'] = 'xx-large'
#mpl.style.use('./assets/thesis_template.mplstyle')
#mpl.style.use('./assets/basic.mplstyle')
#mpl.style.use('./assets/paper.mplstyle')
#mpl.style.use(['./assets/basic.mplstyle','./assets/thesis_template.mplstyle'])

def func_current(magnitude,p,phase):
    return(magnitude*np.cos(np.pi/180*(p-phase)))

def CEP_fit(x,a,b,c,D,E,F):
    return (a*np.exp(-(x-b)**2/c**2)    *     np.sin(D*(x-E))     + F)
# a = j_0
# b = D_0
# c = 2*sig
# D = 
# E = 
# F = c

#k_CEP = D(x-E)

def Theory_fit(D_wedge,j_0,D_0,sig,k_CEP,theta, c):
    return (j_0*np.exp(-(D_wedge-D_0)**2/2*sig**2) * np.cos(2*np.pi*k_CEP*D_wedge + theta)  + c)

def sin_fit(x,a,b,c,D,E,F):
    return (a* np.sin(D*(x-E))     + F)

def trans2ins(t):
    return(2*np.tan(2*np.pi/180)/np.sin(38*np.pi/180)*t)

def power2Efield(P,t,d): # P=average laser power, t=pulse lenght (FWHM), d=beam diameter
    peakpower = P/80e6/t
    S = peakpower/(d**2*np.pi/4)
    E = np.sqrt(377*S)
    return(E)

def get_CEP_fit(wedge_positions, currents_detrend, pinit, maxfev, x_min = 0.45, x_max = 1.0):
    fit_x = np.linspace(x_min,x_max,1000)
    fit_y = {}
    
    for key in np.sort([*currents_detrend]):
        popt,pcov = curve_fit(CEP_fit,trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,p0=pinit, maxfev=maxfev)
        fit_y[key] = CEP_fit(fit_x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

    return (fit_x, fit_y)

def get_sin_fit(wedge_positions, currents_detrend, pinit, maxfev, x_min = 0.45, x_max = 1.0):
    sin_x = np.linspace(x_min,x_max,1000)
    sin_y = {}
    
    for key in np.sort([*currents_detrend]):
        sopt,scov = curve_fit(sin_fit, trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,p0=pinit, maxfev=maxfev)
        sin_y[key] = sin_fit(sin_x,sopt[0],sopt[1],sopt[2],sopt[3],sopt[4],sopt[5])
    
    return (sin_x, sin_y)