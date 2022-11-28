# -*- coding: latin-1 -*-
import math
import numpy as np
#import pylab as plote
TABELA1 = open ('configura.txt', 'w') #Tabela 1 contendo colunas: Elemento qui, posições e indice
TABELA2 = open ('proprieda.txt', 'w') #Tabela 2 contendo colunas: energia potêncial, temperatura e energia cinética
# retornar tabela 1 para .xyz e tabela 2 para .dat
TABELA2.writelines(['#ENER. POTENCIAL(eV)\t' , 'ENER. CINÉTICA(eV)\t' , 'TEMPERATURA(K)\n'])



def forcas_energia_pot(r0,aresta,corte,num):
    E = np.float(0.997e0)                       # E da equacao de Lennard-Jones em Kj/mol
    s = np.float(3.4e0)                         # sigma da equação de Lennard-Jones
    U = np.float(0.0e0)                         #energia potencial 
    F = np.zeros((num,3),float)                 #???
    unitario = np.zeros((3),float)
    for i in range(num-1): 
        for j in range (i+1,num):
            rx = r0[j][0] - r0[i][0]                  #Calcula a norma
            if np.abs(rx) > (aresta/2.0e0):   
               rx = np.abs(rx) - aresta       
            ry = r0[j][1] - r0[i][1]
            if np.abs(ry) > (aresta/2.0e0):   
               ry = np.abs(ry) - aresta
            rz = r0[j][2] - r0[i][2]
            if np.abs(rz) > (aresta/2.0e0):
               rz = np.abs(rz) - aresta
            r = math.sqrt(rx*rx + ry*ry + rz*rz) # distancia entre 2 atomos, dada em A[Angtrom] = e-10 m
            unitario[0] = rx / r
            unitario[1] = ry / r
            unitario[2] = rz / r
            if r <= corte:
               f1 = s / r                                                                             # 3.4/[A ou e-10m]
               f2 = math.pow(f1,12) - math.pow(f1,6)                                                  #Lennard-Jones
               f3 = np.float(48.0e0) * (E / s) * (math.pow(f1,13) - np.float(0.5e0) * math.pow(f1,7)) #Derivada da equação da energia Potencial
               U = U + np.float(4.0e0) * E * f2 * 1.04e-2                                             #Energia Potencial total do processo,e a soma de todas as energias potenciais (em eV)
               F[i] = F[i] + f3 * unitario * 1.0e-10 / 6.022                                          #Calculos de forças
               F[j] = F[j] - f3 * unitario * 1.0e-10 / 6.022
              
    return F, U
    
def integracao_das_equacoes_de_movimento(F,r0,v0,massa,dt,mdt2):
    v = np.zeros((num,3),float)
    r = np.zeros((num,3),float)
    
    for j in range(num):
        r[j] = r0[j] + (v0[j] * dt + F[j] * mdt2 / massa) * np.float(1.0e10)
        v[j] = v0[j] + F[j] * dt / massa

    return r, v

def posicoes(num,aresta_cubo):
    l = aresta_cubo / math.pow(num,(1.0/3.0))
    nterco = np.int(math.pow(num,(1.0/3.0)))
    if num != np.int(math.pow(nterco,3)):
       nterco = np.int(nterco) + np.int(1)
    r0 = np.zeros((num,3),dtype=float)
    npart = 0
    for n in range(nterco):
        for m in range(nterco):
            for w in range(nterco):
                r0[npart][0] = np.float(n) * l
                r0[npart][1] = np.float(m) * l
                r0[npart][2] = np.float(w) * l
                npart = npart + 1
    return r0

def centro_de_massa_em_repolso(v, num , massa):
    pcm = np.zeros((3),float) #momento linear do centro de massa
    
    for j in range(num):
        pcm = pcm + v[j] * massa

    vcm = pcm / (num * massa)

    for j in range(num):
        v[j] = v[j] - vcm
    
    v0 = v        

    return v0

def distribuicao_aleatoria(massa,num,K):
    #Rb = np.float(1.38064852e-23)                                    #constante de Boltzmann = 1,38064852 Ã— 10-23 m2 kg s-2 K-1
    V0 = math.sqrt(2 * K / (num * massa))                            #Calcula a velocidade inicial das molÃ©culas
    v0 = np.random.random((num,3)) - 0.5e0
    #normaliza vetores
    for i in range(num):
        modulo = np.sqrt(v0[i][0]*v0[i][0] + v0[i][1]*v0[i][1] + v0[i][2]*v0[i][2])
        if modulo != 0.0e0:
           v0[i] = v0[i] / modulo
    v0 = v0 * V0
    return v0
    

def recalcular_temperatura(v0,num,massa):
    K = np.float(0.0e0)
    Rb = np.float(1.38064852e-23)                                                            #constante de Boltzmann = 1,38064852x10-23 (m2 kg s-2 K-1)
    
    for i in range(num):
        V2 = math.sqrt (v0[i][0] * v0[i][0] + v0[i][1] * v0[i][1] + v0[i][2] * v0[i][2] )    #**velocidade em modulo de cada átomo
        K = K + 1.0e0/2.0e0 * massa * V2                                                     #energia cinética de cada partí­cula
    T = 2 * K /(3 * num * Rb)                                                                #recalcula a temperatura
    return T, K

def reescala_veloci(v0,Talvo,T):
    v0 = v0 * np.sqrt(Talvo/T)
    return v0

def correcao_posicao_particulas(num, r0, aresta_cubo):
    for j in range (num):
        if r0[j][0] >= aresta_cubo:
            r0[j][0] = r0[j][0] - aresta_cubo
        else:
            if r0[j][0] < 0.0e0:
               r0[j][0] = r0[j][0] + aresta_cubo
            
        if r0[j][1] >= aresta_cubo: 
            r0[j][1] = r0[j][1] - aresta_cubo
        else:
            if r0[j][1] < 0.0e0:
               r0[j][1] = r0[j][1] + aresta_cubo
        
        if r0[j][2] >= aresta_cubo:
            r0[j][2] = r0[j][2] - aresta_cubo
        else:
            if r0[j][2] < 0.0e0:
               r0[j][2] = r0[j][2] + aresta_cubo
    return r0
    
print("PROGRAMA DE DINÃMICA ÃTOMICA\n")
print("VALORES DEFINIDOS:\n")
print("E=0.997 kJ/mol , s=3.4 A, V0 definida pela temperatura inicial\n")
print("Distancia de equili­brio: 3.81 A \n")


num = int(input("Informe o numero de particulas do sistema:\n"))
T = float(input("Informe a temperatura em Kelvin:\n"))
tmax = float(input("Informe o tempo maximo de simulacao(tmax):\n")) #colocar se é em min ou s
dt = float(input("Informe o passo (dt) :\n"))                       #colocar se é em min ou s
L = float(input("Informe o comprimento da caixa:\n"))               #colocar se é em m

Talvo = T                               # temperatura ideal
r = np.float(3.81e0)                    #distancia de equilibrio
rc = r * np.float(2.5e0)                #Raio de corte #AUMENTAR RAIO DE CORTE PARA TESTES
aresta_cubo = L * rc

#Rmax = float(input("Entre com Rmax:\n")) 
#dr = float(input("Entre com dr:\n"))

r0 = posicoes(num,aresta_cubo)          #Subrotína de distribuição dos átomos na caixa

aresta_cubo2 = aresta_cubo / np.float(2.0e0)

"""
print("X_0\n",X0)
print("Y_0\n",Y0)
print("Z_0\n",Z0)
"""

TABELA1.writelines([str(num), '\n', '\n'])
for i in range (num):
    TABELA1.writelines(['Ar  \t' , str(r0[i][0]), '\t' , str(r0[i][1]), '\t' , str(r0[i][2]), '\n' ])

t=0

uma = 1.66e-27      # unidade de massa atomica em kilogramas
massa = 39.94 * uma # massa do átomo de argonio

K = 3 / 2 * num / 6.022e23 * 8.31 * T # num/6.022e23 = numero de moles
#8.31 J/mol*K e T = K

v0 = distribuicao_aleatoria(massa,num,K)

"""
print("v0x_0 \n", v0x)
print("v0y_0 \n", v0y)
print("v0z_0 \n", v0z)
"""



#chama rotina para calcular forcas e energia potencial do sistema
F, U = forcas_energia_pot(r0,aresta_cubo,rc,num)
"""
print("Fx_0 \n", Fx)
print("Fy_0 \n", Fy)
print("Fz_0 \n", Fz)
print("U_0 \n", U )
"""
#determina o numero de iteracoeses totais, tendo conhecimento do intervalo de tempo e do tempo maximo
itmax = int (tmax/dt)
dt = dt * np.float(1.0e-15) #em fs ou e^(-15)s
dt2 = dt * dt
mdt2 = 0.5 * dt2
numesc = 10

#laÃ§o da dinÃ¢mica molecular        
for i in range(itmax):
    #integra as equacoes de movimento
    r0, v = integracao_das_equacoes_de_movimento(F,r0,v0,massa,dt,mdt2)
    
    #agora vamos fazer o centro de massa do sistema ficar em repouso
    v0 = centro_de_massa_em_repolso(v,num,massa)
    
    #Recalcular a temperatura e energia cinética
    T, C = recalcular_temperatura(v0,num,massa)
    C = C * 6.242e+18

    #verifica se a particula sai do cubo e faz correcoes
    r0 = correcao_posicao_particulas(num,r0,aresta_cubo)    

    #chama rotina para calcular forcas e energia potencial do sistema
    F, U = forcas_energia_pot(r0,aresta_cubo,rc,num)

    #reescala velocidades caso temperatura atual seja distinta da temperatura alvo
    v0 = reescala_veloci(v0,Talvo,T)
    
    esc = int(i / numesc)
    if ((i - esc*numesc) == 0):
       TABELA2.writelines([str(i), '\t', str(U) , '\t' , str(C) , '\t\t' , str(T) , '\n'])
       TABELA1.writelines([str(num), '\n', '\n'])
       for j in range(num):
           TABELA1.writelines(['Ar  \t' , str(r0[j][0]) , '\t' , str(r0[j][1]) , '\t' , str(r0[j][2]) , '\n' ])
    
TABELA1.close()
TABELA2.close()        
        
        
  
"""   
print (U)
print (Fx, Fy, Fz)
print (ax,ay, az)
uma = 1.66e-27 # em kilogramas
massa = 39.948 * uma



resp=int(Rmax/dr)

x = np.zeros(resp)
y = np.zeros(resp)
z = np.zeros(resp)
factor1 = np.zeros(resp)
factor2 = np.zeros(resp)
factor3 = np.zeros(resp)

for i in range(0,resp):   
    factor1 = s/x[i]
    factor2 = math.pow(factor1,12) - math.pow(factor1,6)
    factor3 = 48.0e0 * (E/ x[i]) * (math.pow(factor1,12) - 0.5e0 * math.pow(factor1,6))
    x[i] = 0.1e0 + np.float(i) * dr
    y[i]= 4.0e0 * E * factor2
    z[i]= 48.0e0 * (E / x[i]) * factor2
    
    

plote.figure(1,figsize=(5,5))
plote.title('Grafico Potencial')
#plote.axis([min(x),max(x),min(y),max(y)])
plote.axis([2,max(x),min(y),10])
plote.plot(x,y,'g',linewidth=3)




plote.figure(2,figsize=(5,5))
plote.title('Grafico Forca')
#plote.axis([min(x),max(x),min(y),max(y)])
plote.axis([0,max(x),min(z),2])
plote.plot(x,z,'r',linewidth=3)


"""

