import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

def colebrook(re, f, k, dh):
    """
    Calcula o fator de fricção usando a equação de Colebrook-White.
    
    re : número de Reynolds
    k  : rugosidade da tubulação (em metros)
    dh : diâmetro hidráulico (em metros) - para tubo circular, é o mesmo que o diâmetro do tubo
    
    Retorna o valor de λ (fator de fricção)
    """
    # Estimativa inicial de f (aproximação de Swamee-Jain pode ser útil aqui)
    
    tolerance = 1e-6  # Maior precisão
    max_iter = 100  # Número máximo de iterações
    
    for i in range(max_iter):
        # Calcular o lado direito da equação de Colebrook
        f_right = -2 * math.log10((k / (dh * 3.7)) + (2.51 / (re * math.sqrt(f))))
        f_new = 1/f_right**2
        # Atualizar o valor de f com a diferença entre f e f_right
        delta = f_new - f
        
        # Verificar se a diferença é suficientemente pequena (convergência)
        if abs(delta) < tolerance:
            return f
        
        # Atualizar o valor de f
        f = f_new
    
    # Se não convergiu dentro do número máximo de iterações
    raise ValueError("Não convergiu após o número máximo de iterações")


L = 10 # Comprimento da distância entre as tomadas de pressão [m]
D = 0.0207 # Diâmetro da tubulação [m]
e = 0.000045 # Rugosidade do tubo
rho = np.array([847.50, 847.50, 847.50, 847.50, 847.50, 847.50, 847.50]) # Densidade do óleo [kg/m^3]        
v = np.array([1.45e-5, 1.45e-5, 1.45e-5, 1.45e-5, 1.45e-5, 1.45e-5 ,1.45e-5]) # Viscosidade Cinemática do óleo [m^2/s]     

Q = np.array([887, 721, 719, 719, 719, 720, 721]) # Vazão do óleo [L/h]    
Q = Q/(1000*3600)
deltap = np.array([0.0897, 0.0741, 0.0744, 0.0749, 0.0754, 0.0758, 0.0765]) # Diferencial de pressão [bar]     
deltap = deltap*100000
A = math.pi*(D**2)/4
U = Q/A
Re_exp = np.zeros(len(U))
for i in range(len(U)):
    Re_exp[i] = U[i]*D/v[i]
Re_laminar = np.arange(100, 2300, 10)
Re_turb = np.arange(1000, 1000000, 1000)
# Experimental
f_exp = np.zeros(len(Re_exp))
for i in range(len(f_exp)):
    f_exp[i] = (2*D*deltap[i])/(rho[i]*L*(U[i]**2))

# Laminar
f_laminar = 64/Re_laminar

# Blasius 
f_blasius = 0.32/(Re_turb**(0.25))

# Colebrook
f_colebrook = np.zeros(len(Re_turb))
for i in range(len(Re_turb)):
    f_colebrook[i] = colebrook(Re_turb[i], 0.01, e, D)

plt.plot(Re_laminar, f_laminar, label='Laminar', color="blue", linestyle='-')  # Curva 1
plt.plot(Re_turb, f_blasius, label='Blasius', color="red", linestyle='-')  # Curva 2
plt.plot(Re_turb, f_colebrook, label='Colebrook', color="green", linestyle='-') 
plt.plot(Re_exp, f_exp, label='Experimental', color="orange", linestyle='', marker='x') 

# Adicionando título e rótulos
plt.title("Diagrama de Moody")
plt.xlabel("Reynolds")
plt.ylabel("Fator de atrito")

plt.xscale('log')
plt.yscale('log')

plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=10))
plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=10))

plt.xticks([10**2, 10**3, 10**4, 10**5, 10**6], ['10²', '10³', '10⁴', '10⁵', '10⁶'])
plt.yticks([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3], ['0.01', '0.02', '0.03', '0.04', '0.05', '0.06', '0.07', '', '', '0.1', '0.2', '0.3'])


plt.grid(True)

# Adicionando a legenda
plt.legend()

# Exibindo o gráfico
plt.show()
