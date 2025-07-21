import xarray as xr
import xgcm #type:ignore

def gridmetrics(data):
    at, au = data['e1t'] * data['e2t'], data['e1u'] * data['e2u']
    av, af = data['e1v'] * data['e2v'], data['e1f'] * data['e2f']
    vt, vu, vv, vw = data['e3t'] * at, data['e3v'] * au, data['e3v'] * av, data['e3w'] * at
    data = data.update({'at': at, 'au': au, 'av': av, 'af': af, 'vt': vt, 'vu': vu, 'vv': vv, 'vw': vw})
    data = data.set_coords(['at', 'au', 'av', 'af', 'vt', 'vu', 'vv', 'vw'])
    metrics = {
        ('X',): ['e1t', 'e1u', 'e1v', 'e1f'], # X distances
        ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'], # Y distances
        ('Z',): ['e3t', 'e3u', 'e3v', 'e3w'], # Z distances
        ('X', 'Y'): ['at', 'au', 'av', 'af'], # Areas
        ('X', 'Y', 'Z'): ['vt', 'vu', 'vv', 'vw'] # Volumes
        }
    grid = xgcm.Grid(data, metrics=metrics)
    #data = data.update({'llon_rc': (['y_r', 'x_c'], data.llon_cc.data), 'llat_rc': (['y_r', 'x_c'], data.llat_cc.data)})
    #data = data.set_coords(['llon_rc', 'llat_rc'])
    bathy_i = grid.interp(data.bathy, {'X', 'Y'})
    data['bathymetry'] = (('y_r', 'x_r'), bathy_i.data)
    return data, grid, metrics

def gridmetrics2D(data):
    at, au = data['e1t'] * data['e2t'], data['e1u'] * data['e2u']
    av, af = data['e1v'] * data['e2v'], data['e1f'] * data['e2f']
    data = data.update({'at': at, 'au': au, 'av': av, 'af': af})
    data = data.set_coords(['at', 'au', 'av', 'af'])
    metrics = {
        ('X',): ['e1t', 'e1u', 'e1v', 'e1f'], # X distances
        ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'], # Y distances
        ('X', 'Y'): ['at', 'au', 'av', 'af'], # Areas
        }
    grid = xgcm.Grid(data, metrics=metrics)
    bathy_i = grid.interp(data.bathy, {'X', 'Y'})
    data['bathymetry'] = (('y_r', 'x_r'), bathy_i.data)
    return data, grid, metrics

def bathy(zgr_path, data):
    ds_zgr = xr.open_dataset(zgr_path)
    try:
        bathy = ds_zgr.bathy_meter.isel(t = 0)
    except:
        bathy = ds_zgr.hdept.isel(t = 0)
    data = data.update({'bathy': (['y_c', 'x_c'], bathy.data)})        
    return data

def do_all(data, zgr_path):
    #data = meshmask(data, zgr_path)
    data = bathy(zgr_path, data)
    data, grid, metrics = gridmetrics2D(data)
    
    return data, grid, metrics