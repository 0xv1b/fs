import numpy as np

def sin(x,A,B,C,D):
    return(A*np.sin(B*(x-C))+D)

def parabula(x,A,B,C):
    return(A*(x-B)**2+C)

def sec2(x,A,B,C,D):
    return(A/np.cosh(B*(x-C))**2+D)

def gaus(x,A,B,C,D):
    return(A*np.exp(-0.5*((x-B)/C)**2)+D)

def gaus2(x,A,B,C):
    return(A*np.exp(-0.5*((x-B)/C)**2))

def linear(x,A,B):
    return(A*x+B)