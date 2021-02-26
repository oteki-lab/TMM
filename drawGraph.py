import matplotlib
import matplotlib.pyplot as plt

def graph_manager():
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