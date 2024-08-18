import numpy as np
import matplotlib.pyplot as plt

#material properies

pc = 10e+07
pm = 5e+07
h = 1

Ec = 151e+09
Em = 70e+09

nu = 0.3

def powerlaw(z, n, pc, pm, h):
    return (pc-pm)*(0.5+z/h)**n + pm

z_values = np.linspace(-h/2, h/2, 500)

n_values = [1, 0.2, 7]

plt.figure(figsize=(10, 6))

#Power Law Gragh
for n in n_values:
    p_values = powerlaw(z_values, n, pc, pm, h)
    plt.plot(z_values, p_values, label=f'n={n}')

plt.title('Power Law Distribution in FGM')
plt.xlabel('z')
plt.ylabel('p(z)')
plt.legend()
plt.grid(True)
#plt.show()

#Module of Elasticiy
E = powerlaw(0.5, 0, Ec, Em, h)

#Reduced Stiffness Matrix
Q11 = Q22 = E/(1-(nu**2))
Q12 = Q21 = (nu*E)/(1-(nu**2))
Q44 = Q55 = Q66 = E/(2*(1+nu))

Q = np.zeros([5, 5])

Q[0, 0] = Q11
Q[1, 1] = Q22
Q[0, 1] = Q12 = Q[1, 0]
Q[2, 2] = Q44
Q[3, 3] = Q55
Q[4, 4] = Q66

#print(Q)

#e (diffrence between Mid-plane and physical neutral plane)

e = 