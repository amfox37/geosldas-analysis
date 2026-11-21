# FOV > 0 on CF0360: today's rule vs proposed (box-overlap table + box-scaled distance, fallback nearest center
# within max_dist), for SMAP-like (M36 land-cell centers, FOV 20 km) and H SAF-like (12.5 km spaced points, FOV 12.5 km).
import sys, numpy as np
from scipy import ndimage
sys.path.insert(0, '/Users/amfox/Desktop/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSldas_App/util/shared/python')
from read_GEOSldas import read_tilecoord
from EASEv2 import EASEv2_ind2latlon, EASEv2_latlon2ind
R = 6371.0; DTR = np.pi/180
V = '/Users/amfox/Desktop/geosldas-analysis/data/viirs_obs_study/'
cs = read_tilecoord(V+'validation_cf0360/tiles/M21C_test_CF0360_jpss_5d4d636_A2/output/CF0360x6C_GLOBAL/rc_out/M21C_test_CF0360_jpss_5d4d636_A2.ldas_tilecoord.bin')
L = np.asarray(cs['typ']) == 100; NL = int(L.sum())
c = {k: np.asarray(cs[k])[L].astype(np.float64) for k in ['com_lon','com_lat','min_lon','max_lon','min_lat','max_lat']}
d = np.load(V+'validation_cf0360/raster/cf0360_cmg_dominant_tile.npz'); dom = d['dom'][::-1].astype(np.int64)
landc = dom <= NL
NX, NY = 1440, 1080; dlon, dlat = 0.25, 180./NY
def hij(lon, lat):
    return np.clip(np.floor((lon+180)/dlon).astype(int),0,NX-1), np.clip(np.floor((lat+90)/dlat).astype(int),0,NY-1)
def table(keys_per_tile):           # keys_per_tile: (tile index array, cell key array) -> CSR
    t, k = keys_per_tile; o = np.lexsort((t, k)); ks = k[o]
    return t[o], np.searchsorted(ks, np.arange(NX*NY)), np.searchsorted(ks, np.arange(NX*NY), side='right')
ti, tj = hij(c['com_lon'], c['com_lat'])
T_center = table((np.arange(NL), tj*NX+ti))
i0, j0 = hij(c['min_lon'], c['min_lat']); i1, j1 = hij(c['max_lon'], c['max_lat'])
ni, nj = i1-i0+1, j1-j0+1; rep = ni*nj; tt = np.repeat(np.arange(NL), rep)
off = np.arange(rep.sum()) - np.repeat(np.cumsum(rep)-rep, rep)
T_box = table((tt, (j0[tt]+off//ni[tt])*NX + i0[tt]+off%ni[tt]))
hx, hy = np.maximum((c['max_lon']-c['min_lon'])/2,1e-6), np.maximum((c['max_lat']-c['min_lat'])/2,1e-6)
def search(lon, lat, T, metric, need_box):
    tl, st, sp = T; oi, oj = hij(lon, lat); oc = oj*NX+oi; n = sp[oc]-st[oc]
    out = np.full(lon.size, -1); best = np.full(lon.size, 1e30)
    for k in range(n.max()):
        m = np.flatnonzero(n > k); t = tl[st[oc[m]]+k]
        dx, dy = lon[m]-c['com_lon'][t], lat[m]-c['com_lat'][t]
        s = np.abs(dx)+np.abs(dy) if metric == 'mink' else np.hypot(dx/hx[t], dy/hy[t])
        if need_box:
            inb = (lon[m]>=c['min_lon'][t])&(lon[m]<=c['max_lon'][t])&(lat[m]>=c['min_lat'][t])&(lat[m]<=c['max_lat'][t])
            s = np.where(inb, s, 1e30)
        b = s < best[m]; best[m[b]] = s[b]; out[m[b]] = t[b]
    return out
def keep_test(lon, lat, t, fov):    # GEOSldas: drop if outside bbox AND too far (either component > max_dist)
    ok = t >= 0; tt = np.where(ok, t, 0)
    inb = (lon>=c['min_lon'][tt])&(lon<=c['max_lon'][tt])&(lat>=c['min_lat'][tt])&(lat<=c['max_lat'][tt])
    ydeg = fov/DTR/R; xdeg = ydeg/np.cos(lat*DTR)
    near = (np.abs(lon-c['com_lon'][tt])<=xdeg)&(np.abs(lat-c['com_lat'][tt])<=ydeg)
    return np.where(ok & (inb | near), t, -1)
def today(lon, lat, fov): return keep_test(lon, lat, search(lon, lat, T_center, 'mink', False), fov)
def proposed(lon, lat, fov):
    t = search(lon, lat, T_box, 'ell', True)
    fb = t < 0
    t[fb] = keep_test(lon[fb], lat[fb], search(lon[fb], lat[fb], T_box, 'mink', False), fov)
    return t
def truth_cell(lon, lat):
    j = np.clip(((90-lat)/0.05).astype(int),0,3599); i = np.clip(((lon+180)/0.05).astype(int),0,7199); return j, i
def report(name, lon, lat, fov):
    j, i = truth_cell(lon, lat); lat = 89.975-0.05*j; lon = -179.975+0.05*i   # snap to 0.05 deg cell center (truth = its main tile)
    j, i = truth_cell(lon, lat); tr = dom[j, i]; isl = tr <= NL
    near = ndimage.binary_dilation(landc, structure=np.ones((5,5),bool))[j, i] & ~isl   # water point within ~2 CMG cells of land
    print(f'\n{name}: {lon.size} points, FOV {fov} km; land points {isl.sum()}, coastal water points {near.sum()}')
    bands = [(-90,90),(60,75),(75,90)]
    for rn, f in [('today', today), ('proposed', proposed)]:
        t = f(lon, lat, fov); kept = t >= 0; right = kept & (t+1 == tr)
        s = '  '.join(f'{a}..{b}: kept {100*kept[isl&(lat>=a)&(lat<b)].mean():5.1f}% right {100*right[isl&(lat>=a)&(lat<b)].mean():5.1f}%' for a,b in bands)
        print(f'  {rn:9s} land: {s} | coastal water kept {100*kept[near].mean():5.1f}%  | all points kept {kept.sum()}')
# SMAP: M36 land cell centers
m36 = read_tilecoord('/Users/amfox/Desktop/LS_DAv8_M36_as_test3/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_DAv8_M36_as_test3.ldas_tilecoord.bin')
Lm = np.asarray(m36['typ']) == 100
r, cc = EASEv2_latlon2ind(np.asarray(m36['com_lat'])[Lm], np.asarray(m36['com_lon'])[Lm], 'M36')
la, lo = EASEv2_ind2latlon(r, cc, 'M36'); la = np.asarray(la, float); lo = np.asarray(lo, float)
for fov in [20., 5., 0.]: report('SMAP-like (M36 land-cell centers)', lo, la, fov)
# H SAF-like: points ~12.5 km apart (0.1125 deg in lat, 0.1125/cos(lat) in lon), kept if land or coastal water
pts = []
for latv in np.arange(-59.95, 89.95, 0.1125):
    dl = 0.1125/np.cos(latv*DTR); lons = np.arange(-180+dl/2, 180, dl); pts.append(np.c_[lons, np.full(lons.size, latv)])
P = np.vstack(pts); j, i = truth_cell(P[:,0], P[:,1])
sel = landc[j, i] | ndimage.binary_dilation(landc, structure=np.ones((5,5),bool))[j, i]
P = P[sel]
for fov in [12.5, 0.]: report('H SAF-like (12.5 km points)', P[:,0], P[:,1], fov)
