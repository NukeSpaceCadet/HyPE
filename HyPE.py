# -*- coding: utf-8 -*-
"""
Hybrid Performance Estimator (HyPE)
Python Ver. 1, 3/31/2018
Written by Lucas Beveridge

This code simlulates the burn of a hybrid fuel rocket motor, and
generates approximate performance data and plots
"""
import csv
from itertools import zip_longest
import numpy as np
import numpy
import argparse
import configparser
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib as mpl

logo="""
                                        ___
 _     _          ______    _______    |   |
| |   | |        (_____ \  (_______)   |___|
| |__ | | _   _   _____) )  _____       _|_ 
|  __)| || | | | |  ____/  |  ___)     |   |
| |   | || |_| | | |       | |_____    |   |
|_|   |_| \__  | |_|       |_______)   |   |
         (____/                        |___| 
                                        \ /
                                        /_\\
Hybrid Performance Estimator            |||
                                       (|||)
Version 0.0.1                         ((|||)) 
Written by Lucas Beveridge, April 2018
Original code written 2010          
"""
print(logo)


""" Parse command line arguments """
parser = argparse.ArgumentParser()
parser.add_argument("input_file_name", help="input file name with complete path and file extension")
parser.parse_args()
args = parser.parse_args()

# Open input file formatted as a config file (similar to *.ini files)
config = configparser.ConfigParser()
config.sections()
data_in=config.read(args.input_file_name)   
print("input file loaded...")

"""
Parse input file
"""

""" Simulation parameters """
out_name=config['Simulation Parameters']['name'] # file name for output
dt=float(config['Simulation Parameters']['dt']) # time step length, s
tb=float(config['Simulation Parameters']['tb']) # total burn time, s
data_name=config['Simulation Parameters']['data'] # file name

""" Engine Parameters """
#ambient pressure (Pa)
Pam=float(config['Engine Parameters']['Pam'])
#initial tank feed pressure (Pa)
Pt=float(config['Engine Parameters']['Pt'])
#Injector discharge coefficient
Cd=float(config['Engine Parameters']['Cd'])
#total injector nozzle orifice area (all orifces) (m2)
Ai=float(config['Engine Parameters']['Ai'])
#physical state of oxidizer, 1=gas, 2=liquid
OXmodl=int(config['Engine Parameters']['OXmodl'])
#regression rate model, 1-6
rmodl=int(config['Engine Parameters']['rmodl'])
#number of burning ports
Np=float(config['Engine Parameters']['Np'])
#port diameter (m)
D=float(config['Engine Parameters']['D'])
#Port (fuel grain) length (m)
L=float(config['Engine Parameters']['L'])
#nozzle throat area (m2)
At=float(config['Engine Parameters']['At'])
#nozzle expansion ratio
epsilon=float(config['Engine Parameters']['epsilon'])
#nozzle expansion correction
lambdaa=float(config['Engine Parameters']['lambda'])

""" Chemical Properties """
#mass of oxidizer that boils off (thus contributing to feed pressure, increasing ullage) (kg)
mboil=float(config['Chemical Properties']['mboil'])
#oxidizer molar mass (kg/mol)
molO=float(config['Chemical Properties']['molO'])
#oxidizer gamma
gammaO=float(config['Chemical Properties']['gammaO'])
#density of liquid oxidizer (kg/m3)
rhoO=float(config['Chemical Properties']['rhoO'])
#oxidizer temperature in tank (K)
T=float(config['Chemical Properties']['T'])
#fuel density (kg/m3)
rhoF=float(config['Chemical Properties']['rhoF'])
#regression rate coefficient
a=float(config['Chemical Properties']['a'])
#flux regression rate exponent
n=float(config['Chemical Properties']['n'])
#length regression rate exponent
me=float(config['Chemical Properties']['me'])

""" Initial Values """
#Initial oxidizer mass (kg)
mOX=float(config['Initial Values']['mOX'])
#initial ullage volume (m3)
Vol=float(config['Initial Values']['Vol'])
#initial oxidizer mass flow rate (kg/s)
mdotO=float(config['Initial Values']['mdotO'])
#initial O/F
OF=float(config['Initial Values']['OF'])
#initial regression rate (m/s)
rdot=float(config['Initial Values']['rdot'])

""" parameters for other regression rate models """
phi=float(config['Alt Regress Params']['phi'])
d=float(config['Alt Regress Params']['d'])
#multiplier for high-regression rate fuels (i.e. paraffin)
mult=float(config['Alt Regress Params']['mult'])

""" declare other variables """
prop=data_name; # Name of propellant combo lookup table with thermo data
Psea=101325 #sea-level pressure (Pa)
nm=mOX/molO #initial number of moles of oxidizer
Rbar=8.314 #universl gas constant (J/mol/kg)
Sp=numpy.pi*L*D #initial port surface area (m2)
Ap=numpy.pi*(D/2)**2 #initial port cross section area (m2)
Ae=epsilon*At #nozzle exit area (m2)
g=9.806 #gravitational acceleration at sea-level (m/s2)

#load thermo data into array 8 cloumns wide
thermodat = np.loadtxt(prop, skiprows=20)
print("thermochemical data loaded from table: "+data_name)

""" other variable intializations """
count=0 #loop counter
F=[] #thrust array
t=0 #actual time
tm=[] #time index array
mdotF=0 #fuel mdot array
ISP=[] #Isp array
Pc=Pam
Pch=[Pam] 
Vexit=[] 
massF=[]
OFr=[OF]
Pexit=[Pam]
Diam=[D]
massO=[mdotO]
mass=[mdotO+mdotF]
dat_len=len(thermodat[:,0]) #find number of rows in thermo data table
print("all variables initialized...")
print("""
beginning simulation...""")

"""
Parse lookup table for thermochemical data
column 0 (parameter): O/F ratio
column 1 (parameter): chamber pressure (psi)
Data column 2: Chamber Pressure (x10^5 Pa)
Data column 3: Chamber Temperature (K)
Data column 4: Chamber M (g/mol)
Data column 5: Chamber Gamma
Data column 6: Exit C* (m/s)
Data column 7: Exit Isp (s)
"""

"""
========================
   BEGIN CALCULATIONS 
========================
"""

for idx in range(int(math.ceil(tb/dt))): #iterate until until tb has been reached in steps of dt
    
    #error check, if chamber pressure exceeds feed pressure, stop program    
    if Pt<Pc:
        print("ERROR: Chamber pressure is greater than feed pressure. Consider a design change. Burn Time="+str(t)+" s.")
        
    count=count+1
    Gox=(mdotO/Np)/Ap #oxidizer mass flux
    Gf=Gox/OF #fuel mass flux
    G=Gf+Gox #total mass flux
    
    #select and calculate regression rate for this time step
    if rmodl==1:
        rdot=mult*a*(G**n)*L**me
    elif rmodl==2:
        rdot=mult*a*(Gox**n)*L**me
    elif rmodl==3:
        rdot=mult*a*(Gox**n)*(L**me)*(1+(2*a*n*rhoF*L**(1+me))/(D*Gox**(1+n)))
    elif rmodl==4:
        rdot=mult*a*(Gox**n)*(L**me)*(1-numpy.exp((-D)/1.4))*(1-numpy.exp((-Pc)/625))
    elif rmodl==5:
        rdot=mult*a*(Gox**n)*(L**me)*(1-numpy.exp((-D)/1.06))
    elif rmodl==6:
        rdot=mult*a*(Gox**n)*(L**me)*(D**d)*Pc**phi
    else:
        print("ERROR: Incorrect regression rate model specification. Must be an integer 1-6")
    
    mdotF=rdot*rhoF*Sp*Np #calculate fuel mass flow rate
    massF.append(mdotF) #place fuel mdot into list
    
    mdot=mdotF+mdotO #calculate total mass flow rate
    mass.append(mdot) #place mass flow rate in list
    
    OF=mdotO/mdotF #calculate new O/F ratio
    OFr.append(OF) #place OF into list
    
    """ extract data from thermo-chemical table based on Pc and O/F """
    # find row where OF range begins
    for i in range(dat_len):
        if thermodat[i,2]>=Pc/(1e5):
            row1=i
            break
    
    #Find sub-row with where Pc matches calculated
    for j in range(dat_len-row1):
        if thermodat[j+row1,0]>=OF:
            row=row1+j-1 # select the next highest row
            break
    
    T=thermodat[row,3] #extract combustion temperature, K
    M=thermodat[row,4] #extract molecular weight of combustion products  
    gamm=thermodat[row,5] #extract gamma for combustion products
    cstar=thermodat[row,6] #extract c* from table
    
    """ calculate values dependant on table data """
    Pc=mdot*cstar/At #calculate new chamber pressure
    Pch.append(Pc) #place new chamber pressure into list
    
    Ve=numpy.sqrt(((2*gamm)/(gamm-1))*((Rbar/(M/1000)*T))*(1-(Pam/Pc)**((gamm-1)/gamm))) #calculate exit velocity
    Vexit.append(Ve) #place Ve into list
    
    """====================
    the following code numerically solves for the exit mach number """
    func = lambda m : (1/m) * numpy.sqrt(((2/(gamm+1))*(1+(gamm-1)/2 * m**2))**((gamm+1)/(gamm-1)))-epsilon
    mach=fsolve(func,3)
    """===================="""
    
    Pe=Pc*((1+(gamm-1)/2)*mach**2)**(-gamm/(gamm-1)) #calculate exit pressure
    Pexit.append(Pe) #place Pe into list
    
    F.append(lambdaa*(mdot*Ve)) #calculate thrust and place into list
    
    isp=Ve/g #calculate actual Isp
    ISP.append(isp) #put Isp into a list
    
    #update tank/feed pressure and oxidizer volume
    r=Pc/Pt
    Pt=(Pt*Vol)/((mdotO*dt/rhoO)+Vol)
    Vol=Vol+((mdotO)/rhoO)*dt
    
    #update oxidizer flow rate based on compressibility
    if OXmodl==1: #for gaseious oxidizer
        Y = (r**(2 / gammaO) * (gammaO / (gammaO - 1)) * ((1 - r**((gammaO - 1) / gammaO)) / (1 - r)))**.5
        rhoO = (mOX)/(Vol)
        mOX=mOX+(mboil-mdotO)*dt
    else:
        Y=1
        
    mdotO = Cd * Y * Ai * (2 * rhoO * numpy.abs((Pt - Pc)))**0.5 #update oxidizer mass flow rate
    massO.append(mdotO) #place new mdot into list
    
    D=D+2*rdot*dt #update port diameter
    Diam.append(D)
    
    Ap=numpy.pi*(D/2)**2 #update port surface area
    Sp=numpy.pi*L*D #update port burning surface area
    
    t=t+dt #update time
    tm.append(t)

"""
========================
 FINISHED CALCULATIONS 
========================
""" 
print("Simulation completed!")

"""
========================
 GENERATE OUTPUT FILES 
========================
"""
Thrust=np.array(F) #convert F into numpy array for performance calculations
Fuelmassf=np.array(massF) #convert massF into numpy array for performance calculations
Oxmassf=np.array(massO) #convert massO into numpy array for performance calculations
ISPnum=np.array(ISP) #convert ISP into numpy array for performance calculations 

total_impulse=str(numpy.sum(Thrust*dt))
total_fuel=str(numpy.sum(Fuelmassf*dt))
total_ox=str(numpy.sum(Oxmassf*dt))
ave_isp=str(numpy.mean(ISPnum))
 
# Generate header
header=["time (s)","Oxidizer mass flow rate (kg/s)", "Fuel mass flow rate (kg/s)", "Port Diameter (m)","O/F ratio","Exhaust velocity (m/s)","Specific impulse (s)","Thrust (N)"]
output=[tm,massO,massF,Diam,OFr,Vexit,ISP,F] #combine relevant data into a single array
export_data = zip_longest(*output, fillvalue='') #rotate data into columns for spreadsheet

with open(out_name+".csv", 'w', newline='') as csvfile: #open csv file for outputting data
    datawriter = csv.writer(csvfile, delimiter=',',
                            quotechar=' ', quoting=csv.QUOTE_MINIMAL)                      
    datawriter.writerow(["total impulse (N*s):",total_impulse])
    datawriter.writerow(["total fuel mass (kg):",total_fuel])
    datawriter.writerow(["total oxidizer mass (kg):",total_ox])
    datawriter.writerow(["average Isp (s.):",ave_isp])
    datawriter.writerow(header) #save header to first row
    datawriter.writerows(export_data) #save data into subsequent rows
csvfile.close()

print("data saved to output file: "+out_name+".csv")

"""
========================
     GENERATE PLOTS 
========================
"""

# Make a square figure and axes

plt.subplot(2, 1, 1)
lines=plt.plot(tm, F, '-')
plt.title('Thrust & Isp')
plt.ylabel('Thrust (N)')
plt.setp(lines, color='r', linewidth=1.0)

plt.subplot(2, 1, 2)
plt.plot(tm, ISP, '-')
plt.xlabel('Time (s)')
plt.ylabel('Isp (s)')

mpl.rcParams['savefig.dpi'] = 100

plt.savefig(out_name+"_plots.png")

print("plots saved to: "+out_name+"_plots.png")

print("""
Basic performance data:""")
print("total impulse (N*s): "+total_impulse)
print("total fuel mass (kg): "+total_fuel)
print("total oxidizer mass (kg): "+total_ox)
print("average Isp (s.): "+ave_isp)