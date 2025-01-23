"""
Desenvolvido por Enzo Petersen

Código de extrapolação e predição da viscosidade cinemática de óleos em diferentes temperaturas.
"""

"""
Equação de Puttagunta et al. (1992) and Aboul-Seoud and Moharam (1999) através da equação Walther (1931):
ln[ln(v+0.8)]=a1+a2.ln(T)


Equação derivada de MacCoull na ASTM-D341:
loglog(Z)=A-B.log(T)
Z=v+0.7+exp(-1.47-1.84.v-0.51.v^2)
v=[Z-0.7]-exp(-0.7487-3.295[Z-0.7]+0.6119[z-0.7]^2-0.3193[Z-0.7]^3)


Cálculo do índice da viscosidade pela ISO-2909:
VI=[(L-U)/(L-H)] x 100
Se a viscosidade cinemática do óleo em 100°C for acima de 70 mm2/s:
L=0.835v^2+14.67v-216
H=0.168v^2+11.85v-97
Senão, L e H são tabelados.
"""

import math
import numpy as np
def viscosity(v1, T1, v2, T2):
    # Cálculo da Viscosidade Cinemática através Puttagunta et al. (1992) and Aboul-Seoud and Moharam (1999): Walther (1931)
    A1 = np.array([ [1, math.log(T1+273.15)], [1, math.log(T2+273.15)] ])
    b1 = np.array([math.log(math.log(v1+0.8)), math.log(math.log(v2+0.8))])
    z1 = np.linalg.solve(A1,b1)
    T=np.arange(15, 45, 0.1)
    v = np.zeros(len(T))
    for i in range(len(T)):
        v[i] = math.exp(math.exp(z1[0] + z1[1] * math.log(T[i]+273.15))) - 0.8
        #print(v[i])
    # Cálculo da Viscosidade Cinemática através da ASTM D341
    A2 = np.array([[1,-math.log10(T1+273.15)], [1,-math.log10(T2+273.15)]])
    Z21 = v1+0.7+math.exp(-1.47-1.84*v1-0.51*(v1**2))
    Z22 = v2+0.7+math.exp(-1.47-1.84*v2-0.51*(v2**2))
    b2 = np.array([math.log10(math.log10(Z21)), math.log10(math.log10(Z22))])
    z2 = np.linalg.solve(A2,b2)
    Z = np.zeros(len(T))
    vASTM = np.zeros(len(T))
    for i in range(len(T)):
        Z[i] = 10**(10**(z2[0] - z2[1]* math.log10(T[i]+273.15)))
        vASTM[i] = Z[i] - 0.7 - math.exp(-0.7487 - 3.295*(Z[i] - 0.7) + 0.6119*((Z[i] - 0.7)**2) - 0.3193*((Z[i] - 0.7)**3))
        #print(vASTM[i])
    # Cálculo do índice de viscosidade através ISO 2909
    L = 12.447  
    H = 9.8279
    VI = ((L-10)/(L-H))*100
    return (v, vASTM, VI)

viscosity(10,40,2.63,100)
