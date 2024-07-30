import pandas as pd
import matplotlib.colors as colors
import matplotlib.cm as cm

def get_hex_color(val_lst,cmap = 'RdBu_r'):
    """
    Converts a z-value to a hex color using a given colormap.
    
    Parameters:
    z_value (float): The z-value to convert.
    cmap (str): The name of the colormap to use (default: 'RdBu_r').
    
    Returns:
    str: The hex color string.
    """
    # Normalize the z-value to [0, 1]
    hex_lst = []
    z_min = min(val_lst)
    z_max = max(val_lst)

    norm = colors.Normalize(vmin=z_min, vmax=z_max)

    for val in val_lst:
        norm_value = norm(val)
        
        # Get the RGBA color tuple from the colormap
        cmap = cm.get_cmap(cmap)
        rgba_color = cmap(norm_value)
        
        # Convert the RGBA color tuple to a hex color string
        hex_color = colors.rgb2hex(rgba_color)
        hex_lst.append(hex_color)

    return hex_lst

# get zval and corresponding label and hex color
df_wm_zval = pd.read_csv('../output/tables/enigmav_fa.csv')
hex_vals = get_hex_color(df_wm_zval.zval.values)
df_wm_zval['hex_code'] = hex_vals
df_wm_labels = pd.read_csv('../output/tables/enigmav_fa_labels.csv')
df_fwm = pd.DataFrame()
df_fwm['label'] = df_wm_labels.label[df_wm_labels.roi.isin(df_wm_zval.roi)]
df_fwm['label_code'] = df_wm_labels.label_code[df_wm_labels.roi.isin(df_wm_zval.roi)]
df_fwm['hemi'] = df_wm_labels.hemi[df_wm_labels.roi.isin(df_wm_zval.roi)]
df_fwm['hex_code'] = df_wm_zval.hex_code[df_wm_zval.roi.isin(df_wm_labels.roi)]
df_fwm['zval'] = df_wm_zval.zval[df_wm_zval.roi.isin(df_wm_labels.roi)]
df_wm_comb = df_fwm.dropna(subset=['label']) 
lblc = []
lbl = []
hm = []
zval = []
hexc = []
for col in range(len(df_wm_comb)):
    col_lbl = df_wm_comb['label'].iloc[col].split(',')
    lblcode = df_wm_comb['label_code'].iloc[col].split(',')
    if len(col_lbl) > 1:
        lbl.append(col_lbl[0])
        lbl.append(col_lbl[1])
        lblc.append(lblcode[0])
        lblc.append(lblcode[1])
        hm.append(df_wm_comb['hemi'].iloc[col])
        hm.append(df_wm_comb['hemi'].iloc[col])
        zval.append(df_wm_comb['zval'].iloc[col])
        zval.append(df_wm_comb['zval'].iloc[col])
        hexc.append(df_wm_comb['hex_code'].iloc[col])
        hexc.append(df_wm_comb['hex_code'].iloc[col])
    else:
        lbl.append(col_lbl[0])
        label_code = df_wm_comb['label_code'].iloc[col].zfill(2)
        lblc.append(label_code)
        hm.append(df_wm_comb['hemi'].iloc[col])
        zval.append(df_wm_comb['zval'].iloc[col])
        hexc.append(df_wm_comb['hex_code'].iloc[col])

data = {'label_code': lblc,
        'label': lbl,
        'zval': zval,
        'hex_code': hexc,
        'hemi':hm}

df_wm_20 = pd.DataFrame(data)
df_wm_20 = df_wm_20.sort_values('label_code')
df_wm_20.to_csv('../output/tables/enigmav_zval_wm_20_3d.csv')