import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.stats import poisson
from scipy.interpolate import interp1d
print "\n"

def decay_func(x):
    return np.exp(-1.0*x)

def trapz(n_evals):
    x = np.linspace(0.0,1.0,num=n_evals)
    y = decay_func(x)
    #x_sub = np.linspace(0.0,1.0,num=10000)
    #y_sub = decay_func(x_sub)

    #plt.plot(x_sub,y_sub,linewidth=2.0)
    
    #plt.xlabel("x",fontsize=16)
    #plt.ylabel("f(x)",fontsize=16)
    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14)
    #plt.show()
    #exit()
    return np.sum(0.5*(np.roll(x,-1)[:-1]-x[:-1])*(y[:-1]+np.roll(y,-1)[:-1]))

exact = 1.0 - np.exp(-1.0) 

n_evals = 2

err=[]
R_n_1=[0.0]
n_evals_arr = ((np.ones(9)*2.0)**np.arange(1,10)).astype(int)

for n_evals in n_evals_arr:
    err.append(abs(trapz(n_evals)-exact)/exact)
    R_n_1.append(trapz(n_evals))
    #print int(n_evals), np.round(1.0/n_evals,3), trapz(n_evals),abs(trapz(n_evals)-exact)/exact  
err = np.array(err)
R_n_1 = np.array(R_n_1)

plt.loglog(1.0/n_evals_arr,err,marker="o",ls="--")
plt.xlabel(r'h  [$n_\mathrm{evals}^{-1}$]',fontsize=16)
plt.ylabel(r'$\epsilon_\mathrm{r}$',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.show()

R = np.zeros((len(R_n_1),len(R_n_1)))
R[:,1] = R_n_1

for n in range(1,len(R_n_1)-1):
    for m in range(1,len(R_n_1)-1):
        if m>n:
            continue
        R[n+1,m+1] = R[n+1,m] + (1.0/((4**m)-1.0))*(R[n+1,m] - R[n,m])

R_err = np.copy(R)
R_err[R>0.0] = np.absolute(R[R>0.0] - exact)/exact

#print "\n R:",R[1:,1:] #avoid zero-padded values!

#print "\nR_err",R_err[1:,1:]

# Part c: hit or miss
def hit_or_miss_mc(n_evals):
    x = np.random.random(n_evals)
    y = np.random.random(n_evals)
    fn = decay_func(x)

    misses = y[y>fn]
    hits = y[y<fn]
    A=1.0 #entire area
    
    est_integral = A*len(hits)/(1.0*(len(misses)+len(hits)))
    err = np.absolute(est_integral - exact)/exact

    return (est_integral, err)

n_evals_arr = np.array([10,20,50,100,200,500,1000,2000,5000,10000]).astype(int)
I_n=[]
I_n_err=[]
for n_evals in n_evals_arr:
    est_integral, err = hit_or_miss_mc(n_evals)
    I_n.append(est_integral)
    I_n_err.append(err)
    print n_evals, est_integral, err
exit()
I_n=np.array(I_n)
I_n_err=np.array(I_n_err)

plt.loglog(n_evals_arr,I_n_err,marker="o",ls="--",color="orange")
plt.xlabel(r'n  [scales as $h^{-1}$]',fontsize=16)
plt.ylabel(r'$\epsilon_{r}$',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.show()

print "\ndone"
