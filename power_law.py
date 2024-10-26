import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sympy import symbols, integrate

h = 1
z_values = np.linspace(-h/2, h/2, 500)
n_values = [0.1, 0.2, 0.5, 1, 2, 5, 10]

nu = 0.3

#Tempreture change trugh the thikness in kelvin
T0 = 300 

#no teperature change in thickness direction
Tf = 600
deltaT = Tf - T0
T = T0 + deltaT

#temperature change linearly
Tl = 100
Tu = 600

def LT(Tl, Tu, z, T0):
    Tz = Tl + (0.5+z/h)*(Tu-Tl)
    deltaT = Tz - Tl
    T = T0 + deltaT
    return T
T = LT(Tl, Tu, 0.5, T0)

#Tempreture Dependent Matrial Properties
def TDP(P0, P_1, P1, P2, P3, T):
    P = P0*( P_1*(T**-1) +1 + P1*T + P2*(T**2) +P3*(T**3) )
    return P
#print(TDP(244.27e+09, 0, -1.371e-03, 1.214e-06, -3.681e-10, LT(Tl, Tu, 0.5, T0)))

Ec = TDP( 244.27e+09, 0, -1.371e-03, 1.214e-06, -3.681e-10, LT(Tl, Tu, 0.5, T0) )
Em = TDP( 122.56e+09, 0, -4.589e-4, 0, 0, LT(Tl, Tu, 0.5, T0) )

#Simple Power law
def powerlaw(z, n, pc, pm, h):
    #in case of no change in matrial properties with T change
    return (pc-pm)*(0.5+z/h)**n + pm
    #in case of material properties dependent on temprucure change
    #return (TDP(0, 0, 0, 0, 0) - TDP(0, 0, 0, 0, 0))*(0.5+z/h)**n + TDP(0, 0, 0, 0, 0)

def V_power_law(z, n, h):
    return (0.5+z/h)**n

plt.figure(figsize=(10, 6))

#Power Law Gragh
for n in n_values:
    #p_values = powerlaw(z_values, n, Ec, Em, h)
    v_values = V_power_law(z_values, n, h)
    #plt.plot(z_values, p_values, label=f'n={n}')
    plt.plot(z_values, v_values, label=f'n={n}')

# Use MaxNLocator to automatically set a higher number of ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=25))  # 10 ticks on x-axis
plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=25))  # 10 ticks on y-axis
plt.title('Power Law Distribution in FGM')
plt.xlabel('z/h')
plt.ylabel('Volume fraction index')
plt.legend()
plt.grid(True)
plt.show()

#Module of Elasticiy
E = powerlaw(0.5, 0, Ec, Em, h)
print("E: ", E)

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

#e = ( integrate(z*E, (z, -h/2, h/2)) ) / (integrate(E, (z, -h/2, h/2)))

#print("e: ", e)
