import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.stats import poisson
from scipy.interpolate import interp1d
print "\n"

k = np.arange(20)
k_prime_arr = np.arange(20).astype(int)

greater_16=[]
lesser_16=[]
n_o = 5
for k_prime in k_prime_arr:
    cum_prob_high = np.sum(poisson.pmf(k[k>n_o], k_prime)) #is > 16% ?
    cum_prob_low = np.sum(poisson.pmf(k[k<n_o], k_prime)) #is > 16% ?
    greater_16.append(cum_prob_high)
    if k_prime<5:
        lesser_16.append(0.0)
    else:
        lesser_16.append(cum_prob_low)
greater_16=np.array(greater_16)
lesser_16=np.array(lesser_16)

greater_ind= list(greater_16).index(np.max(greater_16[greater_16 <= 0.16]))
lesser_ind= list(lesser_16).index(np.max(lesser_16[lesser_16 <= 0.16]))

interp_greater = interp1d(greater_16, k, kind='linear')
interp_lesser = interp1d(lesser_16, k, kind='linear')
print interp_greater(0.16)
print "classical CI:",interp_greater(0.16),interp_lesser(0.16)
print "classical length:",abs(interp_lesser(0.16) - interp_greater(0.16))
#exit()
plt.plot(k, poisson.pmf(k, n_o),marker="o",ls="--",label=r"$n_{o}$=5")
#plt.plot(k, poisson.pmf(k, k_prime_arr[greater_ind]),marker="o",ls="--",label=r"$\mu_{1}$="+str(k_prime_arr[greater_ind]))
plt.plot(k, poisson.pmf(k, interp_greater(0.16)),marker="o",ls="--",label=r"$\mu_{1}$="+str(np.round(interp_greater(0.16),2)))
plt.plot(k, poisson.pmf(k, interp_lesser(0.16)),marker="o",ls="--",label=r'$\mu_{2}$='+str(np.round(interp_lesser(0.16),2)))
#plt.plot(k, poisson.pmf(k, k_prime_arr[lesser_ind]),marker="o",ls="--",label=r"$\mu_{2}$="+str(k_prime_arr[lesser_ind]))
plt.xlabel("k successes",fontsize=18)
plt.ylabel("P( k successes | N trials )",fontsize=18)
plt.legend(frameon=False,fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
#exit()

# Part b
def poisson_likelihood(n_o,mu_t):
    return ((mu_t**n_o)*np.exp(-1.0*mu_t))/np.math.factorial(n_o)

def uni_prior(mu_t):
    return np.ones(len(mu_t)).astype(float)

def inverse_prior(mu_t):
    return 1.0/mu_t

def bayes_prob_poisson(n_o, mu_t, uniform):
    if uniform==True:
        non_norm_prob = poisson_likelihood(n_o,mu_t)*uni_prior(mu_t)
    else:
        non_norm_prob = poisson_likelihood(n_o,mu_t)*inverse_prior(mu_t)

    norm_prob = non_norm_prob/np.nansum(non_norm_prob)
    norm_prob = np.nan_to_num(norm_prob)
    
    most_prob = mu_t[list(norm_prob).index(np.max(norm_prob))]
    most_prob_index = list(norm_prob).index(np.max(norm_prob))

    # prob value
    lower = np.array([np.sum(norm_prob[:x]) for x in np.arange(len(mu_t))])
    upper = np.array([np.sum(norm_prob[x:]) for x in np.arange(len(mu_t))])
    lower_CI = (lower[lower<0.16])[-1]
    upper_CI = (upper[upper<0.16])[0]
    
    # mu_t associated with prob value:
    lower_mu = mu_t[list(lower).index(lower_CI)]
    upper_mu = mu_t[list(upper).index(upper_CI)]

    return (lower_mu, most_prob, upper_mu, norm_prob)

n_o = 5
mu_t = np.linspace(0,20,num=1000)

lower_mu, most_prob, upper_mu, norm_prob = \
    bayes_prob_poisson(n_o, mu_t, uniform=True)
plt.plot(mu_t, norm_prob,linewidth=1.5,color="purple")
plt.fill_between(mu_t, 0.0, norm_prob, where=mu_t<3.62,color="purple")
plt.fill_between(mu_t, 0.0, norm_prob, where=mu_t>8.39,color="purple")
plt.xlabel(r'$\mu_{t}$',fontsize=18)
plt.ylabel(r'P( $\mu_{t}$ )',fontsize=18)
plt.legend(frameon=False,fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()
print "Bayes uni:",lower_mu, most_prob, upper_mu
print "len Bayes uni:",upper_mu-lower_mu,"\n"

# Bayesian inverse prior:
lower_mu, most_prob, upper_mu, norm_prob = \
    bayes_prob_poisson(n_o, mu_t, uniform=False)

plt.plot(mu_t, norm_prob,linewidth=1.5,color="orange")
plt.fill_between(mu_t, 0.0, norm_prob, where=mu_t<lower_mu,color="orange")
plt.fill_between(mu_t, 0.0, norm_prob, where=mu_t>upper_mu,color="orange")
plt.xlabel(r'$\mu_{t}$',fontsize=18)
plt.ylabel(r'P( $\mu_{t}$ )',fontsize=18)
plt.legend(frameon=False,fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()

print "Bayes inv:",lower_mu, most_prob, upper_mu
print "len Bayes inv:",upper_mu-lower_mu
print ""
print "RMS:",n_o-np.sqrt(n_o), n_o, n_o+np.sqrt(n_o)
print "len RMS:",n_o+np.sqrt(n_o) - (n_o-np.sqrt(n_o))
