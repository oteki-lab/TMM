import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import TMM as tmm
import drawGraph as dg
dg.graph_manager()

def calc(height, stack):
    ## wavelength range
    wl_min = 0.92                                   # [um]
    wl_max = 1.34                                   # [um]

    ## flags
    cal_rta = 0         # calculate the absorption
    cal_abs = 0         # calculate the R, T, A
    cal_field = 0       # calculate the field intensity in each layer
    cal_QD = 0       # calculate the field intensity in each layer

    # structure parameters
    back_b = 50
    surf_b = 25
    spacer = 25
    intensity = 1
    height_QD = 1
    max_N = stack

    ## list of indices
    indice = {}
    for m in ['GaAs', 'QD', 'mirror']:
        ti = pd.read_csv(m+'.csv')
        indice[m] = interpolate.interp1d(ti['wl'], ti['n1'] + 1j*ti['n2'], kind='cubic')

    ## without QD layers
    ## List of layers (material, indices, thickness)
    layers = pd.DataFrame([
        ['GaAs',    height*1e-9]
    ], columns=['m', 'd'])
    m_end = 'mirror'

    # calculate normalized field intensity
    wl, R, T, A, df_d, df_m, df_I = tmm.cal_TMM(indice, wl_min, wl_max, layers, m_end)

    df_mdI = pd.concat([df_m, df_d, df_I], axis=1)

    # I for evaluation
    target = 1192
    I_mean = df_mdI.loc[:,['material', 'depth', target]]
    I_mean = I_mean.rename(columns={target: 'intensity'})
    d_max = I_mean.iloc[-1]['depth']

    # sort in descending order of I
    I_table = I_mean[(I_mean['depth']>=back_b) & (I_mean['depth']<=(d_max-surf_b))]
    I_sorted = I_table.sort_values('intensity', ascending=False)

    # set QD layers
    df_dicts = I_sorted.to_dict(orient='records')
    qd_table = [df_dicts[0]]
    n_stack=max_N
    for i, df_dict in enumerate(df_dicts):
        if all(abs(df_dict['depth']-np.array([d.get('depth') for d in qd_table]))>spacer) and df_dict['intensity']>=intensity and n_stack>1:
            qd_table.append(df_dict)
            n_stack=n_stack-1
    qd_table = pd.DataFrame.from_dict(qd_table).sort_values('depth').reset_index(drop=True)

    del layers, wl, R, T, A, df_d, df_m, df_I, df_mdI, target, I_mean, I_table, I_sorted, df_dicts, n_stack
    ## with QD layers
    ## List of layers (material, indices, thickness)
    d_temp=0
    layer_list = []
    qd_dicts = qd_table.to_dict(orient='records')
    for qd_dict in qd_dicts:
        if d_temp<qd_dict.get('depth'):
            layer_list.append(['GaAs', (qd_dict.get('depth')-d_temp-height_QD)*1e-9])
            layer_list.append(['QD', height_QD*1e-9])
            d_temp = qd_dict.get('depth')
    layer_list.append(['GaAs', (d_max-d_temp)*1e-9])
    layers = pd.DataFrame(layer_list, columns=['m', 'd'])
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
            'lims': {'x':[wl[0]*1e6, wl[-1]*1e6], 'y': [0, 14]},
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


    return df_mI[df_mI['material'] == 'QD'].mean()[1192]