"""
Score obs->tile rules on CF0360 (FOV = 0, obs at every 0.05 deg CMG cell center) against the BCS tile raster.

Truth: the dominant tile of each CMG cell (36 raster pixels), from cf0360_cmg_dominant_tile.npz made on Discover
(GEOSldas root discover_brief_cf0360_raster.md). Today's GEOSldas rule is the cached obs2tile map (checked against
ObsFcstAna). The other rules choose among the land tiles whose bounding box contains the obs.
Run: python -I score_cf0360_tile_rules.py   (~1 min, ~6 GB memory)
"""
import sys, numpy as np
sys.path.insert(0,'/Users/amfox/Desktop/geosldas-analysis/projects/snow_da')
from lib.obs2tile import read_tilecoord, land_tiles, _cmg_lonlat32
V='/Users/amfox/Desktop/geosldas-analysis/data/viirs_obs_study/'
tc=land_tiles(read_tilecoord(V+'validation_cf0360/tiles/M21C_test_CF0360_jpss_5d4d636_A2/output/CF0360x6C_GLOBAL/rc_out/M21C_test_CF0360_jpss_5d4d636_A2.ldas_tilecoord.bin'))
NL=tc['N_tile']; d=np.load(V+'validation_cf0360/raster/cf0360_cmg_dominant_tile.npz')
dom=d['dom'][::-1].astype(np.int64); domn=d['domn'][::-1]; ndist=d['ndist'][::-1]
st=np.load(V+'validation_cf0360/raster/cf0360_raster_tile_stats.npz')
lon32,lat32=_cmg_lonlat32(); NY,NX=lat32.size,lon32.size
mnlo,mxlo,mnla,mxla=[tc[k].astype(np.float64) for k in ('min_lon','max_lon','min_lat','max_lat')]
clo,cla=tc['com_lon'].astype(np.float64),tc['com_lat'].astype(np.float64)
i0=np.searchsorted(lon32,tc['min_lon'],'left'); i1=np.searchsorted(lon32,tc['max_lon'],'right')
la=lat32[::-1]; j0=NY-np.searchsorted(la,tc['max_lat'],'right'); j1=NY-np.searchsorted(la,tc['min_lat'],'left')
nc=np.maximum(i1-i0,0); nr=np.maximum(j1-j0,0)
# all (cell,tile) pairs
rep=nc*nr; T=np.repeat(np.arange(NL),rep); off=np.arange(rep.sum())-np.repeat(np.cumsum(rep)-rep,rep)
J=j0[T]+off//nc[T]; I=i0[T]+off%nc[T]; K=J*NX+I
x=lon32[I].astype(np.float64); y=lat32[J].astype(np.float64)
hwx=np.maximum((mxlo-mnlo)/2,1e-6)[T]; hwy=np.maximum((mxla-mnla)/2,1e-6)[T]
bx=((mnlo+mxlo)/2)[T]; by=((mnla+mxla)/2)[T]
area=((mxlo-mnlo)*(mxla-mnla))[T]
fill=st['cnt'][T+1]/np.maximum(area*14400,1)          # fraction of bbox covered by the tile
cosl=np.cos(np.deg2rad(y))
rules={
 'A nearest center |dlon|+|dlat|':np.abs(x-clo[T])+np.abs(y-cla[T]),
 'B center, dlon*cos(lat), Euclid':np.hypot((x-clo[T])*cosl,y-cla[T]),
 'C center dist / box half-width (sum)':np.abs(x-clo[T])/hwx+np.abs(y-cla[T])/hwy,
 'D center dist / box half-width (ellipse)':np.hypot((x-clo[T])/hwx,(y-cla[T])/hwy),
 'E smallest bbox area':area,
 'F highest bbox fill (pixels/box)':-fill,
 'G ellipse norm / fill':np.hypot((x-clo[T])/hwx,(y-cla[T])/hwy)/np.sqrt(fill),
}
landc=dom<=NL; pure=landc&(domn==36); mix=landc&(ndist>1)
latc=np.repeat(lat32,NX).reshape(NY,NX)
old=np.load(V+"tilemap/cmg2tile_v4_CF0360_M21C_test_CF0360_jpss_5d4d636_A2.ldas_tilecoord.npy").astype(np.int64); old=np.where(old>=0,old+1,0)
for nm,mm in [("today (pert-cell nearest center + bbox)",old)]:
    c=(mm==dom); print(nm, "mapped%", round(100*(mm[dom<=NL]>0).mean(),1), "correct%", round(100*c[dom<=NL].mean(),1))
print(f'{"rule":42s} correct%: all  pure  mixed  60-75N  >75N')
for name,s in rules.items():
    o=np.lexsort((T,s,K)); k=K[o]; f=np.r_[True,k[1:]!=k[:-1]]
    m=np.zeros(NY*NX,np.int64); m[k[f]]=T[o][f]+1; m=m.reshape(NY,NX)
    c=(m==dom)
    def p(msk): return 100*c[msk].mean()
    print(f'{name:42s} {p(landc):9.1f} {p(pure):5.1f} {p(mix):6.1f} {p(landc&(latc>=60)&(latc<75)):7.1f} {p(landc&(latc>=75)):5.1f}')

