import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate

del matplotlib.font_manager.weight_dict['roman']
matplotlib.font_manager._rebuild()
plt.rcParams["font.family"]         = "Times New Roman"     #全体のフォントを設定
plt.rcParams["mathtext.fontset"]    = "stix"                #数式のフォントを設定
plt.rcParams["font.size"]           = 12                    #フォントの大きさ
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
    _fig, _ax = plt.subplots(figsize=figsize)
    _ax.set_xlim(lims['x'])
    _ax.set_ylim(lims['y'])
    _ax.set_xlabel(labels['x'])
    _ax.set_ylabel(labels['y'])
    _ax.xaxis.set_label_position(label_positions['x'])
    _ax.yaxis.set_label_position(label_positions['y'])
    _ax.xaxis.set_ticks_position(label_positions['x'])
    _ax.yaxis.set_ticks_position(label_positions['y'])
    _ax.tick_params(axis='both', which='major', **major_ticks)
    _ax.tick_params(axis='both', which='minor', **minor_ticks)
    if invert_axis['x']:
        _ax.invert_xaxis()
    if invert_axis['y']:
        _ax.invert_yaxis()
    return _fig, _ax


# Define parameters
h=6.626e-34    # Planck's constant
c=2.998e8      # Speed of light
q=1.602e-19    # Elementary charge

wl_min = 0.92                                   # [um]
wl_max = 1.34                                   # [um]
nb_lambda = int((wl_max-wl_min)*1e3)+1          # Number of wavelength
wl = np.linspace(wl_min,wl_max,nb_lambda)*1e-6  # Wavelength vector
length_step = 1e-9                              # discretization step of the stack, for intensity calculation

cal_rta = 1         # calculate the absorption
cal_abs = 1         # calculate the R, T, A
cal_field = 1       # calculate the field intensity in each layer
transform_E = 0     # Transform the data into a function of energy


## list of indices
indice = {}
for m in ['GaAs', 'QD', 'mirror']:
    ti = pd.read_csv(m+'.csv')
    indice[m] = interpolate.interp1d(ti['wl'], ti['n1'] + 1j*ti['n2'], kind='cubic')

## List of layers (material, indices, thickness)
n_0 = np.ones(nb_lambda)                  # Incident medium
layers = pd.DataFrame([
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     171e-9],
    ['QD',      indice['QD'](wl*1e6),       1e-9],
    ['GaAs',    indice['GaAs'](wl*1e6),     62e-9]
], columns=['m', 'n', 'd'])
n_end = indice['mirror'](wl*1e6)

m = layers.m.values.tolist()
n = layers.n.values.tolist()
d = layers.d.values.tolist()

## Initialization
nb_layers=len(n)
length_stack=np.sum(d)
nb_steps=[np.round(x/length_step) for x in d]
stack=np.arange(0,np.sum(nb_steps)+1)#,length_stack)
n_tot = np.vstack((n_0, n, n_end))
Omega=np.empty((nb_lambda,2,2), dtype=np.complex128)
r=np.zeros(nb_lambda, dtype=np.complex128)
t=np.zeros(nb_lambda, dtype=np.complex128)
thetaz=np.zeros((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
Delta=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
Upsilon=np.empty((nb_layers,nb_lambda,2,2), dtype=np.complex128)
I_z=np.empty((nb_layers,int(np.sum(nb_steps)+1),nb_lambda))
A_z=np.empty((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
LDOS_z=np.empty((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
I_poynting=np.zeros((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
I_poynting_z=np.empty((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
A_poynting_z=np.empty((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
thetaz_poynting=np.zeros((nb_layers,int(np.sum(nb_steps)+1),nb_lambda), dtype=np.complex128)
Omega_m_prime=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)
Omega_m=np.empty((nb_layers+1,nb_lambda,2,2), dtype=np.complex128)

material=[m[0]]
for i in range(1,nb_layers+1,1):
    for j in range(int(np.sum(d[:i-1])*1e9), int(np.sum(d[:i])*1e9)):
        material.append(m[i-1])


## Calculation of R and T
if len(n)==0:
    r=(n_0-n_end)/(n_0+n_end)
    t=2*n_0/(n_0+n_end)
else:
    theta = []  #Phase shift in a layer
    for x,nx in enumerate(n):
        theta.append([2*np.pi*nx[i]*d[x]/wl[i] for i,ni in enumerate(nx)])

    for j in range(0,nb_lambda,1):
        for i in range(nb_layers-1,-1,-1):
            Delta[i+1][j]=(1/(2*n_tot[i+1][j]))*(n_tot[i+1][j]+n_tot[i+2][j]*np.array([[1,-1],[-1,1]]))

            if i==nb_layers-1:
                Omega_m[-1][j]=Delta[-1][j] # transfer matrix between layers 1 and i
            else:
                Omega_m[i+1][j]=np.dot(Delta[i+1][j], Omega_m_prime[i+1][j])

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
    for j in range(0,nb_lambda,1):
        for i in range(0,nb_layers,1):
            I_poynting[i][j]=abs(t[j])**2/(n_0[j]).real*(n_tot[i+1][j]*np.conj(Omega_m_prime[i][j][0][0]+Omega_m_prime[i][j][1][0])*(Omega_m_prime[i][j][0][0]-Omega_m_prime[i][j][1][0])).real
    Abs_layer=I_poynting[0:nb_layers-1,:]-I_poynting[1:nb_layers,:]

## Calculation of intensity of electric field with depth
if cal_field==1:
    for j in range(0,nb_lambda,1):
        for i in range(0,nb_layers,1):
            for k in range(0,int(nb_steps[i])+1,1):
                thetaz[i][k][j]=(k-1/2)*length_step*2*np.pi*n[i][j]/wl[j]
                I_z[i][k][j]=(n[i][j]).real/(n_0[j]).real*abs(t[j]*(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz[i][k][j])+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz[i][k][j])))**2
                A_z[i][k][j]=(4*np.pi*(n[i][j]/wl[j]).imag)*I_z[i][k][j]
                #LDOS_z[i][k][j]=abs(1+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz[i][k][j])/(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz[i][k][j])))**2

            #for k in range(0,int(nb_steps[i])+2,1):
            #    thetaz_poynting[i][k][j]=(k-1)*length_step*2*np.pi*n[i][j]/wl[j]
            #    I_poynting_z[i][k][j]=abs(t[j])**2/(n_0[j]).real#*(n_tot[i+1][j]*np.conj(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz_poynting[i][k][j])+Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz_poynting[i][k][j]))*(Omega_m_prime[i][j][0][0]*np.exp(1j*thetaz_poynting[i][k][j])-Omega_m_prime[i][j][1][0]*np.exp(-1j*thetaz_poynting[i][k][j]))).real

            #for k in range(0,int(nb_steps[i]+1),1):
            #    A_poynting_z[i][k][j]=(I_poynting_z[i][k][j]-I_poynting_z[i][k+1][j])/length_step

    I=I_z[0][:int(nb_steps[0])+1]
    A_local=A_z[0][:int(nb_steps[0])+1]
    #LDOS=LDOS_z[0]
    #A_poynting_local=A_poynting_z[0]
    #I_poynting_local=I_poynting_z[0]
    for i in range(1,nb_layers,1):
        I=np.concatenate([I,I_z[i][:int(nb_steps[i])]],axis=0)
        A_local=np.concatenate([A_local,A_z[i][:int(nb_steps[i])]],axis=0)
        #LDOS=np.concatenate([LDOS,LDOS_z[i]],axis=0)
        #A_poynting_local=np.concatenate([A_poynting_local,A_poynting_z[i]],axis=0)
        #I_poynting_local=np.concatenate([I_poynting_local,I_poynting_z[i][2:-1,:]],axis=0)
    #I_poynting_local=(I_poynting_local[2:-1,:]+I_poynting_local[1:-2,:])/2

df_m = pd.DataFrame(material, columns=['material'])
df_d = pd.DataFrame(stack, columns=['depth'])
df_I = pd.DataFrame(I, columns=(wl*1e9).astype(np.int64))


## Plots
# Total, QDs, Mirror
if cal_abs==1:
    fig, ax = graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 100]},
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
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 100]},
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
    for i in range(0,nb_layers-1,1):
        ax.plot([wl[0]*1e6, wl[-1]*1e6], [sum(d[0:i])*1e6, sum(d[0:i])*1e6], color='white',  linestyle='solid', linewidth = 1.0)
    plt.show()

if cal_field==1:
    #df_mdI = pd.concat([df_m, df_d, df_I], axis=1)
    df_mI = pd.concat([df_m, df_I], axis=1)
    fig, ax = graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 100]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': 'R, T, A (%)'},
        'label_positions': {'x':'bottom', 'y': 'left'},
        'major_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'minor_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'invert_axis': {'x': False, 'y': False}
    })
    ax.plot(wl*1e6, (df_mI[df_mI['material'] == 'QD'].mean()), label="Total", color='gold',  linestyle='solid', linewidth = 1.0)
    plt.legend()
    plt.tight_layout()
    plt.show()