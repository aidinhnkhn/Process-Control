import scipy.optimize as op
import numpy as np
from math import exp
import scipy.integrate
import matplotlib.pyplot as plt




def proportional(PV,SP,Kp=22654):
    return Kp*(SP-PV)

def PID(Kp, Ki, Kd, MV_bar=0):
    
    e_prev = 0
    t_prev = -100
    I = 0
    
    MV = MV_bar
    
    while True:
        # yield MV, wait for new t, PV, SP
        t, PV, SP = yield MV
        
        # PID calculations
        e = SP - PV
        
        P = Kp*e
        I = I + Ki*e*(t - t_prev)
        D = 0
       
        MV = MV_bar + P + I + D
        # update stored data for next iteration
        e_prev = e
        t_prev = t


# model ode's are in this function
def model (t,z,K,PT):

    CA1, T1, CA2,T2 = z
    u= [100,1,350,98.99]
    q ,CAf,Tf,qc= u
    
    # steady state vars which were obtained from previous files
    Tcf=350
    V1=100
    V2=100
    hA1= 1.67e5
    hA2 = 1.67e5
    Cp=0.239
    Cpc=0.239
    deltaH = -4.78e4
    rho = 1000
    rhoc = 1000
    EdivR = 1e4
    k0 = 7.2e10

    CA2setPoint = 0.005
    qc_ss = 98.99

    if t>=PT:
        CA2setPoint*=K
    
    
    qc=controller.send([t,CA2,CA2setPoint])
    if (qc>500): print("ERR")
    
    dCA1 = (q/V1)*(CAf-CA1)-k0*CA1*exp(-EdivR/T1)

    dT1=(q/V1)*(Tf-T1)-(k0*deltaH*CA1)/(rho*Cp)*exp(-EdivR/T1)+(rhoc*Cpc)/(rho*Cp*V1)*qc*(1-exp(-hA1/(rhoc*Cpc*qc)))*(Tcf-T1)

    dCA2=(q/V1)*(CA1-CA2)- k0*CA2*exp(-EdivR/T2)

    dT2 = (q/V2)*(T1-T2)-(k0*deltaH*CA2)/(rho*Cp)*exp(-EdivR/T2)+(rhoc*Cpc)/(rho*Cp*V2)*qc*(1-exp(-hA2/(rhoc*Cpc*qc)))*(T1-T2+exp(-hA1/(rhoc*Cpc*qc))*(Tcf-T1))

    return [dCA1,dT1,dCA2,dT2]

x0 = np.array([0.08506,442,0.005,449.9116]) #initial values
K=1.1  # disturbance
PT = 5 # u(t-PT)

#solves the 4 ode
k=[7022.74,9571.68]## Tyreus-Luyben
#k=[16000,0] #to find Kcu
kp,ki=k
controller = PID(kp,ki,0,98.99) #define controller
controller.send(None)

ans = scipy.integrate.solve_ivp(model,[0,20],x0,args=(K,PT),dense_output=True,method='LSODA',rtol=1e-7,atol=1e-7)

t = np.linspace(0, 20, 3000)

z= ans.sol(t)

#plotting the graphs
plt.plot(t,z[2])
plt.xlabel('Time (min)')
plt.ylabel('CA2 (mol/lit)')
plt.legend(['CA2'])
plt.title('+10% CA2 closed loop')
plt.grid('True')

plt.figure(2)
plt.plot(t,z[0])
plt.xlabel('Time (min)')
plt.ylabel('CA1 (mol/lit)')
plt.legend(['CA1'])
plt.title('+10% CA2 closed loop')
plt.grid('True')

plt.figure(3)
plt.plot(t,z[1])
plt.xlabel('Time (min)')
plt.ylabel('T1 (K)')
plt.legend(['T1'])
plt.title('+10% CA2 closed loop')
plt.grid('True')

plt.figure(4)
plt.plot(t,z[3])
plt.xlabel('Time (min)')
plt.ylabel('T2 (K)')
plt.legend(['T2'])
plt.title('+10% CA2 closed loop')
plt.grid('True')

plt.show()
