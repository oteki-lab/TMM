
import os
import csv
import datetime
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
from scipy import interpolate

del matplotlib.font_manager.weight_dict['roman']
matplotlib.font_manager._rebuild()
plt.rcParams["font.family"]         = "Times New Roman"     #全体のフォントを設定
plt.rcParams["mathtext.fontset"]    = "stix"                #数式のフォントを設定 
plt.rcParams["font.size"]           = 12                    #フォントの大きさ
plt.rcParams["xtick.direction"]     = "in"                  #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"]     = "in"                  #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True                  #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True                  #y軸補助目盛りの追加
plt.rcParams["xtick.major.width"]   = 1.5                   #x軸主目盛り線の線幅
plt.rcParams["ytick.major.width"]   = 1.5                   #y軸主目盛り線の線幅
plt.rcParams["xtick.minor.width"]   = 1.0                   #x軸補助目盛り線の線幅
plt.rcParams["ytick.minor.width"]   = 1.0                   #y軸補助目盛り線の線幅
plt.rcParams["xtick.major.size"]    = 4                     #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"]    = 4                     #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"]    = 2                     #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"]    = 2                     #y軸補助目盛り線の長さ
plt.rcParams["axes.linewidth"]      = 1.5                   #囲みの太さ

def graph_setting(figsize, lims, labels, label_positions, major_ticks, minor_ticks, invert_axis):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(lims['x'])
    ax.set_ylim(lims['y'])
    ax.set_xlabel(labels['x'])
    ax.set_ylabel(labels['y'])
    ax.xaxis.set_label_position(label_positions['x'])
    ax.yaxis.set_label_position(label_positions['y'])
    ax.xaxis.set_ticks_position(label_positions['x'])
    ax.yaxis.set_ticks_position(label_positions['y'])
    ax.tick_params(axis='both', which='major', **major_ticks)#direction='in', bottom=True, top=True, left=True, right=True)#, length=6, width=2, colors='r', grid_color='r', grid_alpha=0.5)
    ax.tick_params(axis='both', which='minor', **minor_ticks)#direction='in', bottom=True, top=True, left=True, right=True)#, length=6, width=2, colors='r', grid_color='r', grid_alpha=0.5)
    if invert_axis['x']: ax.invert_xaxis()
    if invert_axis['y']: ax.invert_yaxis()
    return fig, ax
    

# Define parameters
h=6.626e-34    # Planck's constant
c=2.998e8      # Speed of light
q=1.602e-19    # Elementary charge

wl_min = 0.5                                # [um]
wl_max = 1.8                                # [um]
nb_lambda = int((wl_max-wl_min)*1e3)+1      # Number of wavelength
wl = np.linspace(0.5,1.8,nb_lambda)*1e-6    # Wavelength vector
length_step = 1e-9                          # discretization step of the stack, for intensity calculation

cal_rta = 1         # calculate the absorption
cal_abs = 1         # calculate the R, T, A
cal_field = 1       # calculate the field intensity in each layer
transform_E = 0     # Transform the data into a function of energy

## List of layers
n_0 = 1*np.ones(nb_lambda)                  # Incident medium
#n_end=retindice_chen(wl*1e6,1.7)   
mirror = pd.read_csv('mirror.csv')
nmirror = interpolate.interp1d(mirror['wl'], mirror['n1'] + 1j*mirror['n2'], kind='cubic')
n_end = nmirror(wl*1e6)                     # Last medium

## list of indices of the stack
#n[0]=retindice_semicond(wl*1e6,4.702)
n = []
GaAs = pd.read_csv('GaAs.csv')
nGaAs = interpolate.interp1d(GaAs['wl'], GaAs['n1'] + 1j*GaAs['n2'], kind='cubic')
n.append(nGaAs(wl*1e6))  #n = pd.DataFrame({ 0 :GaAs['n1'] + 1j*GaAs['n2']})

## list of thicknesses of the stack
d=[]
d.append(1782e-9)

## Initialization
nb_layers=len(n)
length_stack=np.sum(d)
nb_steps=[np.round(x/length_step) for x in d]
stack=np.arange(0,np.sum(nb_steps))#,length_stack)
n_tot = np.vstack((n_0, n))
n_tot = np.vstack((n_tot, n_end))
Omega=np.empty((nb_lambda,2,2), dtype=np.complex128)
r=np.zeros(nb_lambda, dtype=np.complex128)
t=np.zeros(nb_lambda, dtype=np.complex128)
thetaz=np.zeros((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
Delta=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
Upsilon=np.empty((nb_layers,nb_lambda,2,2), dtype=np.complex128)
I_z=np.empty((nb_layers,int(np.sum(nb_steps)),nb_lambda))
A_z=np.empty((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
LDOS_z=np.empty((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
I_poynting=np.zeros((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
I_poynting_z=np.empty((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
A_poynting_z=np.empty((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
thetaz_poynting=np.zeros((nb_layers,int(np.sum(nb_steps)),nb_lambda), dtype=np.complex128)
Omega_m_prime=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
Omega_m=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)

## Calculation of R and T
if len(n)==0:
    r=(n_0-n_end)/(n_0+n_end)
    t=2*n_0/(n_0+n_end)
else:
    theta = []  #Phase shift in a layer
    for x,nx in enumerate(n):
        theta.append([2*np.pi*nx[i]*d[x]/wl[i] for i,ni in enumerate(nx)])
    
    for j in range(0,nb_lambda-1,1):
        for i in range(nb_layers-1,-1,-1):
            Delta[i+1][j]=(1/(2*n_tot[i+1][j]))*(n_tot[i+1][j]+n_tot[i+2][j]*np.array([[1,-1],[-1,1]]))

            if i==nb_layers-1:
                Omega_m[-1][j]=(Delta[-1][j]) # transfer matrix between layers 1 and i
            else:
                Omega_m[i+1][j]=Delta[i+1][j]*Omega_m_prime[i+1][j]
            
            Upsilon[i][j]=[[np.exp(-1j*theta[i][j]), 0], [0, np.exp(1j*theta[i][j])]]    
            Omega_m_prime[i][j]=np.dot(Upsilon[i][j], Omega_m[i+1][j])

        Delta[0][j]=(1/(2*n_tot[0][j]))*(n_tot[0][j]+n_tot[1][j]*np.array([[1,-1],[-1,1]]))
        Omega_m[0][j]=np.dot(Delta[0][j], Omega_m_prime[0][j])
        Omega_m_prime[-1][j]=[[1, 0], [0, 0]]
        Omega[j]=Omega_m[0][j] # Transfer matrix of the system
        r[j]=Omega[j][1][0]/Omega[j][0][0] # Reflection coefficient
        t[j]=1/Omega[j][0][0] # Transmission coefficient

R=abs(r)**2                         # Reflectance
T=(n_end.real/n_0.real)*abs(t)**2   # Transmittance
A=1-R-T                             # Absorption of the whole stack
                                                                                                                                                                                                      

## Calculation of absorptivity for each layer
if cal_abs==1:
    for j in range(0,nb_lambda-1,1):
        for i in range(0,nb_layers,1):
            I_poynting[i][j]=abs(t[j])**2/(n_0[j]).real*(n_tot[i+1][j]*np.conj(Omega_m_prime[i][j][0][0]+Omega_m_prime[i][j][1][0])*(Omega_m_prime[i][j][0][0]-Omega_m_prime[i][j][1][0])).real
    Abs_layer=I_poynting[0:nb_layers-1,:]-I_poynting[1:nb_layers,:]

## Calculation of intensity of electric field with depth
if cal_field==1:
    for j in range(0,nb_lambda-1,1):
        for i in range(0,nb_layers,1):
            for k in range(0,int(nb_steps[i]),1):
                thetaz[i][k][j]=(k-1/2)*length_step*2*np.pi*n[i][j]/wl[j]
                I_z[i][k][j]=(n[i][j]).real/(n_0[j]).real*abs(t[j]*(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz[i][k][j])+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz[i][k][j])))**2

                #A_z[i][k][j]=(4*np.pi*(n[i][j]/wl[j]).imag)*I_z[i][k][j]
                #LDOS_z[i][k][j]=abs(1+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz[i][k][j])/(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz[i][k][j])))**2

            #for k in range(0,int(nb_steps[i]),1):
            #    thetaz_poynting[i][k][j]=(k-1)*length_step*2*np.pi*n[i][j]/wl[j]
            #    I_poynting_z[i][k][j]=abs(t[j])**2/(n_0[j]).real#*(n_tot[i+1][j]*np.conj(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz_poynting[i][k][j])+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz_poynting[i][k][j]))*(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz_poynting[i][k][j])-Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz_poynting[i][k][j]))).real

            #for k in range(0,int(nb_steps[i]-1),1):
            #    A_poynting_z[i][k][j]=(I_poynting_z[i][k][j]-I_poynting_z[i][k+1][j])/length_step
    
    I=I_z[0]
    #A_local=A_z[0]
    #LDOS=LDOS_z[0]
    #A_poynting_local=A_poynting_z[0]
    #I_poynting_local=I_poynting_z[0]
    for i in range(1,nb_layers,1):
        I=np.concatenate([I,I_z[i]],axis=0)
        #A_local=np.concatenate([A_local,A_z[i]],axis=0)
        #LDOS=np.concatenate([LDOS,LDOS_z[i]],axis=0)
        #A_poynting_local=np.concatenate([A_poynting_local,A_poynting_z[i]],axis=0)
        #I_poynting_local=np.concatenate([I_poynting_local,I_poynting_z[i][2:-1,:]],axis=0)   
    #I_poynting_local=(I_poynting_local[2:-1,:]+I_poynting_local[1:-2,:])/2


## Plots
# Total, QDs, Mirror
if cal_abs==1:
    fig, ax = graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 1]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': 'R, T, A (%)'},
        'label_positions': {'x':'bottom', 'y': 'left'},
        'major_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'minor_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'invert_axis': {'x': False, 'y': False}
    })
    ax.plot(wl*1e6, (A+T)*100, label="Total", color='gold',  linestyle='solid', linewidth = 1.0)
    ax.plot(wl*1e6, A*100, label="QDs", color='red',  linestyle='solid', linewidth = 1.0)
    ax.plot(wl*1e6, T*100, label="Ag", color='black',  linestyle='solid', linewidth = 1.0)
    plt.legend()
    plt.tight_layout()
    plt.show()

# R, T, A
if cal_rta==1:
    fig, ax = graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 1]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': 'R, T, A (%)'},
        'label_positions': {'x':'bottom', 'y': 'left'},
        'major_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'minor_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'invert_axis': {'x': False, 'y': False}
    })
    ax.plot(wl*1e6, R*100, label="R", color='gold',  linestyle='solid', linewidth = 1.0)
    ax.plot(wl*1e6, T*100, label="T", color='aqua',  linestyle='solid', linewidth = 1.0)
    ax.plot(wl*1e6, A*100, label="A", color='red',  linestyle='solid', linewidth = 1.0)
    plt.legend()
    plt.tight_layout()
    plt.show()

if cal_field==1:
    fig, ax = graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0.0, max(stack)/1e3]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': r'Depth $\, \mathrm{(\mu m)}$'},
        'label_positions': {'x':'top', 'y': 'left'},
        'major_ticks': {'direction': 'out', 'bottom': False, 'top': True, 'left': True, 'right': False},
        'minor_ticks': {'direction': 'out', 'bottom': False, 'top': False, 'left': False, 'right': False},
        'invert_axis': {'x': False, 'y': True}
    })
    img = ax.contourf(wl*1e6,stack/1e3, I, levels=np.arange(round(np.min(I)),round(np.max(I))+0.1,0.1), cmap='plasma')
    cbar = plt.colorbar(img, ticks=np.arange(round(np.min(I)),round(np.max(I))+2,2))
    cbar.ax.tick_params(axis='y', direction='out')
    cbar.ax.minorticks_off()
    cbar.set_label(r'Normalized Field Intensity $\, \|E\|^2/\|E_0\|^2$')
    plt.show()

"""
    for i in range(0,nb_layers-1):
        plot([wl[0], wl[-1]]*1e6, [sum(d[0:i]), sum(d[0:i])]*1e6, 'linewidth',1, 'color','w')
"""