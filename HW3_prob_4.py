import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.stats import poisson
from scipy.interpolate import interp1d
print "\n"

err=[]

n_evals_arr = ((np.ones(17)*2.0)**np.arange(1,18)).astype(int)

#n_evals = 10000
for n_evals in n_evals_arr:
    
    mc_data = np.random.random((n_evals,10))

    val = np.mean(np.sum(mc_data,axis=1)**2.0)
    exact= 155.0/6.0
    print n_evals, val, np.absolute(val-exact)/exact
    err.append(np.absolute(val-exact)/exact)
err = np.array(err)*100.0
exit()
plt.semilogy(np.sqrt(n_evals_arr),err,marker="o",ls="--")
plt.xlabel(r'$\sqrt{n_{evals}}$',fontsize=16)
plt.ylabel(r'$\epsilon_{r}$ [%]',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.show()

print "\ndone"
