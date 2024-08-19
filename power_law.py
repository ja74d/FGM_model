import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, integrate

#material properies
pc = 10e+07
pm = 5e+07
h = 1

Ec = 151e+09
Em = 70e+09

nu = 0.3

#Tempreture change trugh the thikness
T = 1

#Tempreture Dependent Matrial Properties
def TDP(P0, P_1, P1, P2, P3):
    P = P0*( P_1*(T**-1) +1 + P1*T + P2*(T**2) +P3*(T**3) )
    return P
#print(TDP(244.27e+09, 0, -1.371e-03, 1.214e-06, -3.681e-10))

#Simple Power law
def powerlaw(z, n, pc, pm, h):
    #in case of no change in matrial properties with T change
    return (pc-pm)*(0.5+z/h)**n + pm
    #in case of material properties dependent on temprucure change
    #return (TDP(0, 0, 0, 0, 0) - TDP(0, 0, 0, 0, 0))*(0.5+z/h)**n + TDP(0, 0, 0, 0, 0)


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
z, n = symbols('z n')
E = (Ec-Em)*((z/h)+0.5)**n + Em

e = ( integrate(z*E, (z, -h/2, h/2)) ) / (integrate(E, (z, -h/2, h/2)))
#print(e)
