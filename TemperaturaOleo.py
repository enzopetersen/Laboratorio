
"""
Desenvolvido por Enzo Petersen.

Código de modelagem da temperatura do óleo em escoamento através da medição da temperatura externa da tubulação por meio de termopar.
"""

import math

def temperatura_oleo(t_amb, t_par, vazao, rho_oil, mu_oil):
    de = 0.0334  # Diâmetro externo [m]
    di = 0.0207  # Diâmetro interno [m]
    L = 40  # Comprimento da tubulação [m]
    k_steel = 16  # Condutividade térmica do aço [W/m.K]
    k_oil = 0.14  #  Condutividade térmica do óleo [W/m.K]
    vazao = (vazao/1000)/3600  # Vazão do óleo [m^3/s]
    U_oil = vazao/(math.pi*(di**2)/4)  # Velocidade do óleo [m/s]
    cp_oil = 2000  # Capacidade térmica específica do óleo [J/kg.K]
    Pr_oil = mu_oil*cp_oil/k_oil  # Número de Prandtl do óleo [-]
    Re_oil = rho_oil*U_oil*di/mu_oil  # Número de Reynolds do óleo [-]
    if Re_oil < 2300:
        Nu_oil = 3.66  # Número de Nusselt do óleo [-]
    else:
        if Re_oil < 4000:
            Nu_oil = 3.66+((0.068*(Re_oil*Pr_oil*di/L))/(1+0.04*((Re_oil*Pr_oil*di/L)**(2.0/3.0))))
        else:
            Nu_oil = 0.023*(Re_oil**(3.0/4.0))*(Pr_oil**(0.3))
    h_oil = k_oil*Nu_oil/di  # Coeficiente de convecção do óleo [W/m^2.K]
    k_air = 0.026  # Condutividade térmica do ar [W/m.K]
    alpha_air = 2.2e-5  # Difusividade térmica do ar [m^2/s]
    nu_air = 1.5e-5  # Viscosidade cinemática do ar [m^2/s]
    Pr_air = nu_air/alpha_air  # Número de Prandtl do ar [-]
    beta_air = 1/(((t_par+t_amb)/2)+273.15)  # Coeficiente de expansão térmica do ar [-]
    Ra_air = (9.81*beta_air*(t_par-t_amb)*(de**3))/(nu_air*alpha_air)  # Número de Rayleith do ar [-]
    Nu_ar = (0.6+((0.387*(Ra_air**(1/6)))/(((1+(0.559/Pr_air)**(9/16)))**(8/27))))**2  # Número de Nussel do ar [-]
    h_air = k_air*Nu_ar/de  # Coeficiente de convecção do ar [W/m^2.K]
    Q = (t_par-t_amb)/(1/(math.pi*de*h_air))  # Fluxo de calor por unidade de comprimento [W/m]
    t_supi = Q*(math.log(de/di)/(2*math.pi*k_steel)) + t_par  # Temperatura da superfície interna do tubo [°C]
    t_oil = Q*(1/(math.pi*di*h_oil)) + t_supi  # Temperatura do óleo [°C]
    
    return t_oil

print(temperatura_oleo(25, 30, 750, 850, 0.0003))