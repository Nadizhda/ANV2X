# Objective function of the sum rate that targets the frequency elemnts only,and the gaurd interval (Adaptive numerology scheme) 
# The minimization methods of SciPy work with real arguments only.
#  But minimization on the complex space Cn amounts to minimization on R2n, the algebra of complex numbers never enters the consideration.
#  Thus, adding two wrappers for conversion from Cn to R2n and back, you can optimize over complex numbers.
#  I added two functions one to conver he complex values to real C^n space-----> R^(2n) space
# then we use the absolute value to construct the objective function 
# def main()
import scipy
import argparse
from scipy.optimize import minimize
import math
import numpy as np
import pandas as pd
import cmath
import traceback
import sys
import os

def Add_zeros(vars, mod_symb):
    N_IFFT = vars
    L= len(mod_symb)
    if int(N_IFFT) >L:
        zeros= list(np.zeros(int(N_IFFT) -L))
        mod_symb= mod_symb+ zeros
        N_symb =1
    elif int(N_IFFT) <L:
        
        N_symb =math.ceil( L/N_IFFT)  
        #print(N_symb,L)
        zeros = list(np.zeros(int(N_IFFT) * N_symb -L))
        
        mod_symb= mod_symb+zeros       
    elif  int(N_IFFT) ==L: 
        mod_symb= mod_symb
        N_symb =1
    return mod_symb,N_symb

# FFT formulation
def Fourier_Transf(vars,mod_symb):
    [x,alpha,beta,N_IFFT,y,M] = vars
    L= len(mod_symb)

    [zerpadding_modsymb,N_symbls] = Add_zeros(N_IFFT, mod_symb)
    
    if N_symbls==1:
        Y_k=[]
        for k in range(int(N_IFFT)):
            # print(N_IFFT)
            # print(len(zerpadding_modsymb), int(N_IFFT))
            Y_k.append([sum(zerpadding_modsymb[n]*cmath.exp(-1j*2*math.pi*n*k/int(N_IFFT)) for n in range (int(N_IFFT)))])
    else:
        Y_k=[]
        zerpadding_modsymb = np.array(zerpadding_modsymb)
        zerpadding_modsymb = zerpadding_modsymb.reshape([-1,int(int(N_IFFT))])
        [sz1,sz2] = zerpadding_modsymb.shape
        # print(sz1,sz2)
        
        for sym in range(sz1):
            Y_kd=[]
            for k in range(int(N_IFFT)):
                Y_kd.append([sum(zerpadding_modsymb[sym][n]*cmath.exp(-1j*2*math.pi*n*k/int(N_IFFT)) for n in range(int(N_IFFT)))])
            Y_k.append(Y_kd)
        Y_k = np.array(Y_k).reshape([1,-1])
    return Y_k

def SINR_calc(vars,symb, H_x,f_d, numb):
    [x,alpha,beta,N_IFFT,y,M] = vars
    S=[]
    smbls = np.array(symb)
    H_x = np.array(H_x)
    epsilon = f_d/x
    #print(H_x)
    #print(symb)
    Rxp= [ (np.sin(cmath.pi*epsilon)/(int(N_IFFT)*np.sin(cmath.pi*epsilon/int(N_IFFT)) ))*\
        (np.real(smbls[i]*H_x[i])*np.real(smbls[i]*H_x[i])) * (np.exp(1j*cmath.pi*epsilon*(int(N_IFFT) -1)/int(N_IFFT))) for i in range(len(smbls))]

    for i in range(len(smbls)):
        Int= (np.exp(1j*cmath.pi*epsilon*(int(np.ceil(N_IFFT) -1)/int(N_IFFT)))) * sum([(np.real(smbls[ii]*H_x[ii])*np.real(smbls[ii]*H_x[ii]))*\
        (np.exp(1j*cmath.pi*(ii-i)*(int(N_IFFT) -1)/int(N_IFFT))) for ii in range(len(smbls))]) - Rxp[i]

        S.append(Int)

    if numb==0:
        return Rxp
    elif numb==1:
        return S
    

def objective_function1(vars,mod_sym=None,N_noise=None,H_x=None,f_d=None,Bw=None,t=None,):

    [x,alpha,beta,N_IFFT,y,M] = vars
   
    # x= x[0] if isinstance(x, (list, np.ndarray)) else x
    # N_IFFT= x[4] if isinstance(x, (list, np.ndarray)) else None

    # print(len(mod_sym),len(H_x),len(N_noise))
    Sympls= Fourier_Transf(vars,mod_sym,)
    L=len(Sympls)
    H_xd= Fourier_Transf(vars,H_x,)
    noise=Fourier_Transf(vars,N_noise)

    I =  SINR_calc(vars,Sympls, H_xd,f_d,1)
   
    # I= np.array(I)
    Rxp =  SINR_calc(vars,Sympls, H_xd,f_d,0)
    #print(Rxp)
    #print(I)
    Rxp= np.array(Rxp)

    SINR = [Rxp[0][i]/ (noise[0][i] + I[0][i]) for i in range(int(len(Rxp[0])/2))]
    #print(SINR)
    # SINR = [  (Rxp[i])* cmath.exp(1j*2*cmath.pi*(f_d/ x)*(N_IFFT-1)/N_IFFT)/(N_noise[i] + I[i])*\
    
    #                                         cmath.exp(1j*2*cmath.pi*(f_d/ x)*(N_IFFT-1)/N_IFFT) for i in range(len(Rxp))]#for a in Sympls for b in N_noise for c in H_x for d in mod_sym]
    Rate =  np.concatenate((np.real(SINR), np.imag(SINR)))
    return Rate,SINR

def Data_rate(SINR):
    #print(SINR)
    CQI= range(0,16) 
    SNR= [-1.889, -0.817, 0.954, 2.984,4.889, 7.39, 8.898, 11.02, 13.32, 14.68, 16.62, 18.91,21.58, 24.88, 29.32]
    Rate = [0,0.1532,0.2344, 0.3770, 0.6016, 0.8770, 1.1758, 1.4766,1.9141, 2.4036, 2.7305, 3.3223, 3.9023, 4.5234,5.1152,5.5547]
    SNR2 = [10**(a/10) for a in SNR]
    SE=[]
    CQI_u=[]
    for b in SINR:
        #print(SINR)
        a=np.abs(b)
       # print(a)
        for i in range(len(SNR2)-1):
            if a <SNR2[0]:
                se = Rate[0]
                cqi=CQI[0]
            if SNR2[i]<=a<=SNR2[i+1]:
                se = Rate[i]
                cqi=CQI[i]
            elif a>SNR2[-1]:
                se =Rate[-1]
                cqi=CQI[-1]
        SE.append(se)
        CQI_u.append(cqi)
    return (SE,CQI_u )

def objective_function2(vars, mod_sym=None,N_noise=None,H_x=None,f_d=None,Bw=None,t=None,):

    [x,alpha,beta,N_IFFT,y,M] = vars

    # convert from real to complex and apply the rate 

    Rate,SINR = objective_function1(vars, mod_sym,N_noise,H_x,f_d,Bw,t)

    Result = np.sum(Rate22)/(Bw*((1/x)+y)) # SE optimization 
    
    return   Result 

# Printing the values of the onjective functio 

# def print_componet(N_noise=None,H_x=None,modu_symbol=None,f_d=None,=None,t=None,):
#     # print ('At x =', x,':')
#     objective_value = objective_function2(modu_symbol,N_noise,H_x,f_d,Bw,t,)
#     # print (objective_value)
#     print("The value of the objective function is", objective_value )
# Constraint equations


# # Initial parameters for testing the code.
# alpha_0 = 2.5
# beta_0 = 5
# t = 2e-6
def constraint_equations(vars,Bw):
    x,alpha,beta,N_IFFT,y,M = vars

    for k in range(7,10+1):

        return  np.ceil(N_IFFT - pow(2,k)), #np.ceil(N_IFFT-Bw/x)

#     return  x - 1/z # * y # return [y + z - 1 / x, y - alpha * t, z - beta * y]

def constraint_equations2(vars, t):

    x,alpha,beta,N_IFFT,y,M= vars

    return  y-alpha*t , (1/(4*x))-y

#def constraint_equations3(vars,t):

   # x,alpha,beta,N_IFFT,y,M= vars

   # return   (1/(beta* y+1e-9)) -x   #  1/(beta* y) -x 

# % gaurd band based on the 3GPP definetion  
def constraint_equations4(vars, Bc,Tc ):
    x,alpha,beta,N_IFFT,y,M= vars
    
    return Bc-x, x- (1/Tc)

def N_RBconstraint_equations(vars,Bw ):
    x,alpha,beta,N_IFFT,y,M= vars
    return ( (Bw/(x*12)) -M)
 
def Read_Data (file_path):
    Data=[]
    with open(file_path) as f:
        # Data = np.array([float(line.strip()) for line in f])
        lines = f.readlines()       
        lines = ([e.strip().replace('\n', '') for e in lines if e != ''])
        Data = [float(item) for subitem in lines for item in subitem.split()]
    return Data

def find_SCS(vale):
    scs = [15*1e3*pow(2,i) for i in range (3)]
    diff= [abs(scs[i]*1e3 - vale) for i in range(len(scs))]
    min= np.min(diff)
    inex = diff.index(min)
    scs_value = scs[inex]
    return scs_value

def main(file_path):
    try:

        input_data = Read_Data (file_path)
        
        #print(input_data[0:17],input_data[-1])

        # y0 = 4.69e-6 # 
        Bw = int(input_data[0])

        inital_xvalue=input_data[1:5]


        [x0,alpha0,beta0,N_IFFT0] = inital_xvalue

        # N = int(N_IFFT0)
        M=int(Bw/ (x0*12))

        # print(M)

        f_d= input_data[13]
        
        if f_d==0:
           f_d= 1e-9
        Tc = 1/f_d

        t= input_data[14]
        
        
        N = int(input_data[15])

        y0 = input_data[16] #Considered as the symbol time of 15 KHz and the maximum guard band is Ts/4


        bound = [(int(input_data[5]),int(input_data[6])),(int(input_data[7]),int(input_data[8])),\
                  (int(input_data[9]),int(input_data[10])),(int(input_data[11]),int(input_data[12])),(0,1/(2*int(input_data[6]))), (10,300)] # change the maximum guard time from 
        
        #  Arrangement of data in text
        #  write_text(Bw,scs1,2,5,1024,SCS_min,SCS_max,2,5,5,15,128,1024,fD(i),tau_max, Noise,h_x,rceiv_symbol);
        b1 =np.array([ arg1 for arg1 in (input_data[17:N+17]) ])
        b2 =np.array([ arg2 for arg2 in (input_data[N+17:2*N+17]) ])
        Sympls = [b1[i] +1j*b2[i] for i in range(len (b1)) ] #OFDM modulated symbols for further usage
        #print(Sympls)
        b1 =np.array([ arg1 for arg1 in (input_data[2*N+17:3*N+17]) ])
        b2 =np.array([ arg2 for arg2 in (input_data[3*N+17:4*N+17]) ])
        N_noise= [b1[i] +1j*b2[i] for i in range(len (b1)) ]

        b1 =np.array([ arg1 for arg1 in (input_data[4*N+17:5*N+17]) ])
        b2 =np.array([ arg2 for arg2 in (input_data[5*N+17:6*N+17]) ])   
        H_x= [b1[i] +1j*b2[i] for i in range(len (b1)) ]
        #print(H_x)
        
        # modu_symbol=np.array([ arg1[i]1j*b2[i] for i  in range([6*N+1:6*N+1+(len(input_data[6*N+1:-1])/2)]) ]) 

        end=None

        L = int(len(input_data[6*N+17:end]))

        b1 =np.array([ arg1 for arg1 in (input_data[6*N+17:int(6*N+(L/2)+17)]) ])

        b2 =np.array([ arg2 for arg2 in (input_data[int(6*N+(L/2)+17):6*N+L+17]) ])

        modu_symbol= [b1[i] +1j*b2[i] for i in range(len (b1)) ]

        rms_d = input_data[-1]

        #print(modu_symbol)

        Bc = 1/(5*rms_d) #  Captuuring the effect of the delay through RMS delay spread. 
        
       # print(Bc)
       # print(len(modu_symbol),len(Sympls),len(H_x),len(N_noise))
       

        con= {'type': 'eq','fun': constraint_equations,'args':(Bw,)}

        # con11= {'type': 'eq','fun': constraint_equations11,'args':()}

        #con1= {'type': 'ineq','fun': constraint_equations1,'args':(t)}

        con2= {'type': 'ineq','fun': constraint_equations2,'args':(t,)}

        #con3= {'type': 'ineq','fun': constraint_equations3,'args':(t,)}

        con4= {'type': 'ineq','fun': constraint_equations4,'args':(Bc,Tc,)}

        con5= {'type': 'ineq','fun': N_RBconstraint_equations,'args':(Bw,)}

        constraint = [con,con2,con4,con5]

        # print(inital_xvalue) method="Nelder-Mead",
            #    options={'disp':True, 'fatol':1e-04}

        result = minimize(objective_function2, (np.array([x0,alpha0,beta0,N_IFFT0,y0,M])),\

                        bounds=bound, args=(Sympls,N_noise,H_x,f_d,Bw,t,), constraints=constraint,method='SLSQP', options={'maxiter':10 })
        scs_value = find_SCS(result.x[0]) 
        y_optimal=result.x[1]*t
        z_optimal= 1/result.x[0]#result.x[2] * y_optimal
        M_optimal=result.x[5]
        #if y_optimal > z_optimal/4:
        #   result.fun=0
        
        #print ('optimal value:',result.fun)#/(Bw * ((1/scs_value) + result.x[4])))
        #print ('optimal point', result.x)

        # Print optimal solution
        #print("Optimal value of x:", result.x[0])
        #print("Optimal value of alpha:", result.x[1])
        #print("Optimal value of beta:", result.x[2])
        #print("Optimal value of FFT size:", result.x[3])
        #print("Optimal value of y:", y_optimal)
        #print("Optimal value of z:", z_optimal)
        #print("Optimal value of M:", M_optimal)

        f_result= [result.fun, result.x[0], result.x[1], result.x[2],result.x[3],y_optimal,z_optimal,M]

        print(f_result)

        return f_result

    except Exception as e:
        # Print the full traceback
        traceback.print_exc()
        # Print the error message
        print(f"Error: {e}")

if __name__ == "__main__":
  # We  add it as an arguments, the arguments are saved in a data file (text file),
  #that contains the digitally modulated symbols (in case it is adopted) or the input binary data after transmission through the channel to further set 
  # Alpha and Beta for a given numerology and the proper Tg.

    file_path = 'path to the data file'  
    f_result= main (file_path)

    globals()['computed_result'] = f_result
