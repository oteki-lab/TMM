import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import TMM as tmm

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
    fig, ax = graph_setting(**{
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

if cal_QD==1:
    fig, ax = graph_setting(**{
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
    fig, ax = graph_setting(**{
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
