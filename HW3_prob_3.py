import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.stats import poisson
from scipy.interpolate import interp1d
print "\n"

def int_func(x):
    return (np.sin(x))*(np.sin(x))

def trapz(n_evals):
    x = np.linspace(0.0,4.0*np.pi,num=n_evals)
    x_sub = np.linspace(0.0,4.0*np.pi,num=10000)
    y = int_func(x)
    y_sub = int_func(x_sub)
    plt.plot(x_sub,y_sub,linewidth=2.0)
    plt.plot(x,y,ls="none",marker="o")
    val = np.sum(0.5*(np.roll(x,-1)[:-1]-x[:-1])*(y[:-1]+np.roll(y,-1)[:-1]))
    plt.title("N= "+str(n_evals)+"   I= "+str(val))
    plt.xlabel("x [rad]",fontsize=16)
    plt.ylabel("f(x)",fontsize=16)
    #plt.legend()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
    exit()
    return np.sum(0.5*(np.roll(x,-1)[:-1]-x[:-1])*(y[:-1]+np.roll(y,-1)[:-1]))
trapz(32)
exit()
exact = 2.0*np.pi

err=[]
R_n_1=[0.0]
n_evals_arr = ((np.ones(7)*2.0)**np.arange(1,8)).astype(int)
#n_evals_arr = (np.ones(9)*3.0)**np.arange(1,10)
#n_evals_arr = np.arange(2,6)
for n_evals in n_evals_arr:
    #trapz(n_evals)
    err.append(abs(trapz(n_evals)-exact)/exact)
    R_n_1.append(trapz(n_evals))
#exit()
err = np.array(err)
R_n_1 = np.array(R_n_1)
#print exact
#print R_n_1

R = np.zeros((len(R_n_1),len(R_n_1)))
R[:,1] = R_n_1

for n in range(1,len(R_n_1)-1):
    for m in range(1,len(R_n_1)-1):
        if m>n:
            continue
        R[n+1,m+1] = R[n+1,m] + (1.0/((4**m)-1.0))*(R[n+1,m] - R[n,m])

R_err = np.copy(R)
R_err[R>0.0] = np.absolute(R[R>0.0] - exact)/exact

print "\n R:",R[1:,1:] #avoid zero-padded values!

print "\nR_err",R_err[1:,1:]

#e_crit = 10.0**(-5.0)
