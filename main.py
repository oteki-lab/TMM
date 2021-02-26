import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import TMM as tmm
import drawGraph as dg
dg.graph_manager()

cal_rta = 0         # calculate the absorption
cal_abs = 0         # calculate the R, T, A
cal_field = 1       # calculate the field intensity in each layer
cal_QD = 0       # calculate the field intensity in each layer

wl_min = 0.92                                   # [um]
wl_max = 1.34                                   # [um]

## list of indices
indice = {}
for m in ['GaAs', 'QD', 'mirror']:
    ti = pd.read_csv(m+'.csv')
    indice[m] = interpolate.interp1d(ti['wl'], ti['n1'] + 1j*ti['n2'], kind='cubic')

## List of layers (material, indices, thickness)
layers = pd.DataFrame([
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    171e-9],
    ['QD',      1e-9],
    ['GaAs',    62e-9]
], columns=['m', 'd'])
m_end = 'mirror'

wl, R, T, A, df_d, df_m, df_I = tmm.cal_TMM(indice, wl_min, wl_max, layers, m_end)


df_mI = pd.concat([df_m, df_I], axis=1)
df_mdI = pd.concat([df_m, df_d, df_I], axis=1)

## Plots
# Total, QDs, Mirror
if cal_abs==1:
    fig, ax = dg.graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 100]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': 'Absorption (%)'},
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
    fig, ax = dg.graph_setting(**{
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

if cal_QD==1:
    fig, ax = dg.graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 100]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': r'Normalized Field Intensity $\, \|E\|^2/\|E_0\|^2$'},
        'label_positions': {'x':'bottom', 'y': 'left'},
        'major_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'minor_ticks': {'direction': 'in', 'bottom': True, 'top': True, 'left': True, 'right': True},
        'invert_axis': {'x': False, 'y': False}
    })
    ax.plot(wl*1e6, (df_mI[df_mI['material'] == 'QD'].mean()), label="Total", color='gold',  linestyle='solid', linewidth = 1.0)
    plt.legend()
    plt.tight_layout()
    plt.show()

if cal_field==1:
    fig, ax = dg.graph_setting(**{
        'figsize': (8,6),
        'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [df_d.melt().value.min()/1e3, df_d.melt().value.max()/1e3]},
        'labels': {'x':r'Wavelength $\, \mathrm{(\mu m)}$', 'y': r'Depth $\, \mathrm{(\mu m)}$'},
        'label_positions': {'x':'top', 'y': 'left'},
        'major_ticks': {'direction': 'out', 'bottom': False, 'top': True, 'left': True, 'right': False},
        'minor_ticks': {'direction': 'out', 'bottom': False, 'top': False, 'left': False, 'right': False},
        'invert_axis': {'x': False, 'y': True}
    })
    img = ax.contourf(wl*1e6,df_d.melt().value/1e3, df_I, levels=np.arange(round(df_I.min().min()),round(df_I.max().max())+0.1,0.1), cmap='plasma')
    
    cbar = plt.colorbar(img, ticks=np.arange(round(df_I.min().min()),round(df_I.max().max())+2,2))
    cbar.ax.tick_params(axis='y', direction='out')
    cbar.ax.minorticks_off()
    cbar.set_label(r'Normalized Field Intensity $\, \|E\|^2/\|E_0\|^2$')
    if cal_QD==1:
        for d in df_mdI[df_mdI['material'] == 'QD'].depth:
            ax.plot([wl[0]*1e6, wl[-1]*1e6], [d/1e3, d/1e3], color='white',  linestyle='solid', linewidth = 1.0)
    
    func = interpolate.interp2d(wl*1e6,df_d.melt().value/1e3, df_I)
    plt.gca().format_coord = (lambda x,y: 'Wavelength={x:.3f}um, Depth={y:.4f}um, I={z:.1f}'.format(x=x, y=y, z=np.take(func(x, y), 0)))
    plt.show()
