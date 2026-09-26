"""Close-pair innovation correlation for CYGNSS L1 vs distance, footprint overlap, same-track (fixed OL, 3 yr)."""
import glob, sys, numpy as np, netCDF4 as nc, pandas as pd, datetime as dt
sys.path.insert(0,'/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/obs_scaling_params')
from obs_scaling.tile_io import read_tilecoord
P='/discover/nobackup/projects/land_da/cygl1_operator_test'
OL=P+'/OLv8_M36_AZ_fixedop/output/SMAP_EASEv2_M36_GLOBAL'
tc=read_tilecoord(glob.glob(OL+'/rc_out/*tilecoord.bin')[0])
c=nc.Dataset(P+'/scaling_params/cygnss_l1_z_score_clim/AZ_CYGNSS_L1_zscore_all_pentads.nc4'); CL={v:np.ma.filled(c[v][:].astype('f8'),np.nan) for v in ['o_mean','o_std','m_mean','m_std','n_data','m_min','m_max']}; c.close()
_pre={}
def pre(day):
    if day not in _pre:
        fs=glob.glob(f'{P}/CYGNSS_L1/Y{day[:4]}/M{day[4:6]}/*_{day}_all_cyg.nc4')
        if not fs: _pre[day]=None; return None
        d=nc.Dataset(fs[0]); g=lambda v: d[v][:]
        ts,tcnt=g('tile_start'),g('tile_count'); ig,jg,w=g('tile_ig'),g('tile_jg'),g('coefficient_weight')
        fp=[]
        for s,n in zip(ts,tcnt):
            ww=np.asarray(w[s:s+n],'f8'); tot=ww.sum()
            fp.append(dict(zip((np.asarray(ig[s:s+n])*1000+np.asarray(jg[s:s+n])).tolist(),(ww/tot if tot>0 else ww).tolist())))
        _pre[day]=dict(lon=np.asarray(g('sp_lon'),'f8'),lat=np.asarray(g('sp_lat'),'f8'),y=np.asarray(g('observed_y_db'),'f8'),
                       sc=np.asarray(g('sc_num')),ch=np.asarray(g('ch_id')),t=np.asarray(g('ddm_timestamp_utc_sec'),'f8'),inc=np.asarray(g('sp_inc_angle'),'f8'),fp=fp)
        d.close()
        for k in [k for k in _pre if k<day and k!=day][:-2]: del _pre[k]   # keep memory small
    return _pre[day]
def match(day,lon,lat,y):
    out=[]
    for dd in (day,(dt.datetime.strptime(day,'%Y%m%d')-dt.timedelta(days=1)).strftime('%Y%m%d'),(dt.datetime.strptime(day,'%Y%m%d')+dt.timedelta(days=1)).strftime('%Y%m%d')):
        p=pre(dd)
        if p is None: continue
        k=np.where((np.abs(p['lon']-lon)<2e-3)&(np.abs(p['lat']-lat)<2e-3)&(np.abs(p['y']-y)<1e-3))[0]
        if len(k): return p,k[0]
    return None,None
rows=[]; nobs=nmatch=0
files=sorted(glob.glob(OL+'/ana/ens_avg/Y20*/M*/*.nc4'))
for fi,f in enumerate(files):
    t=f.split('.')[-2]
    if t[:4]=='2023': continue
    d=nc.Dataset(f)
    if d.dimensions['n_obs'].size==0: d.close(); continue
    m=d['species'][:]==13
    if m.sum()<2: d.close(); continue
    lon,lat,o,fc,tl=[np.asarray(d[v][:][m],'f8') for v in ('lon','lat','obs','fcst','tilenum')]; d.close()
    k=tl.astype(int)-1; i_,j_=tc.i_indg[k],tc.j_indg[k]
    doy=dt.datetime.strptime(t[:8],'%Y%m%d').timetuple().tm_yday; p_=min((doy-1)//5,72)
    om,os_,mm,ms,nd=[CL[v][p_,j_,i_] for v in ('o_mean','o_std','m_mean','m_std','n_data')]
    sc=np.clip(mm+ms/os_*(o-om),CL['m_min'][j_,i_],CL['m_max'][j_,i_]); ok=np.isfinite(sc)&(nd>=20)
    dinn=sc-fc
    idx=np.where(ok)[0]
    if len(idx)<2: continue
    L=np.c_[lon[idx],lat[idx]]; D=np.hypot(L[:,0,None]-L[None,:,0],L[:,1,None]-L[None,:,1])
    ii,jj=np.where(np.triu(D<0.6,1))
    if len(ii)==0: continue
    info={}
    for q in set(ii.tolist())|set(jj.tolist()):
        info[q]=match(t[:8],lon[idx[q]],lat[idx[q]],o[idx[q]]); nobs+=1; nmatch+=info[q][0] is not None
    for a,b in zip(ii,jj):
        pa,ka=info[a]; pb,kb=info[b]
        if pa is None or pb is None: continue
        fa,fb=pa['fp'][ka],pb['fp'][kb]
        ov=sum(min(v,fb.get(key,0)) for key,v in fa.items())      # footprint overlap fraction (sum of min weights)
        same_track=(pa['sc'][ka]==pb['sc'][kb]) and (pa['ch'][ka]==pb['ch'][kb]) and abs(pa['t'][ka]-pb['t'][kb])<30
        same_sc=pa['sc'][ka]==pb['sc'][kb]
        rows.append((t,D[a,b],ov,same_track,same_sc,abs(pa['t'][ka]-pb['t'][kb]),abs(pa['inc'][ka]-pb['inc'][kb]),dinn[idx[a]],dinn[idx[b]],tl[idx[a]]!=tl[idx[b]]))
    if fi%1000==0: print(f'{fi}/{len(files)} files, {len(rows)} pairs, match rate {nmatch/max(nobs,1):.3f}',flush=True)
df=pd.DataFrame(rows,columns=['t','dist','overlap','same_track','same_sc','dt_s','dinc','d1','d2','diff_tile'])
df.to_parquet('/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/obs_error_correlation/l1_pairs_3yr.parquet')
print('pairs',len(df),'match rate',nmatch/max(nobs,1))
