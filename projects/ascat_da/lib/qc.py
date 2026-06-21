"""
QC configuration and filter functions for Legacy BUFR and H121 CDR obs.

QC dicts are plain Python dicts — scripts pass {**QC_DEFAULT_H121, 'tc_max': 5}
to override individual thresholds without touching the rest.

Setting any threshold to None disables that filter.
"""

import numpy as np

# ── Default QC configs ────────────────────────────────────────────────────────

QC_DEFAULT_BUFR = dict(
    ssm_min=0,
    ssm_max=100,
    smpf_ok={0},           # soilMoistureProcessingFlag values to accept
    smcf_ok={0, 4},        # soilMoistureCorrectionFlag values to accept
    tpcx_max=10,           # topographicComplexity ≤ N%
    iwfr_max=10,           # inundationAndWetlandFraction ≤ N%
)

QC_DEFAULT_H121 = dict(
    ssm_min=0,
    ssm_max=100,
    sf_water_bit=1,        # reject if surface_flag open-water bit is set
    pf_bad_bits=3,         # reject if processing_flag bits 0x01 or 0x02 are set
    pf_ok=None,            # optional exact processing_flag values to accept
    cf_ok=None,            # GEOS HSAF reader does not screen correction_flag
    tc_max=10,             # reject topographic_complexity >= N%
    wf_max=10,             # reject wetland_fraction >= N%
    subsfc_max=5,          # reject subsurface_scattering_probability >= N% (Hahn et al. 2026, Sect. 3.5)
    sens_min=1.0,          # reject surface_soil_moisture_sensitivity <= N dB (Hahn et al. 2026, Sect. 4.1.1)
    bsflag_bad_bits=4,     # reject if backscatter40_flag bit 4 (noise_out_of_limits) is set.
                           # GEOSldas assigns every accepted obs the same fixed
                           # this_obs_param%errstd regardless of flag status (no
                           # per-obs noise down-weighting exists), so this is the
                           # only available lever against noisy backscatter.
)


# ── Filter functions ──────────────────────────────────────────────────────────

def apply_bufr_qc(arrays, qc=None):
    """Return boolean mask of obs passing QC.

    arrays: dict with keys matching BUFR field names (ssm, smpf, smcf, tpcx, iwfr).
    qc:     QC config dict; defaults to QC_DEFAULT_BUFR.
    """
    if qc is None:
        qc = QC_DEFAULT_BUFR

    ssm  = np.asarray(arrays['ssm'],  float)
    smpf = np.asarray(arrays['smpf'], float)
    smcf = np.asarray(arrays['smcf'], float)
    tpcx = np.asarray(arrays['tpcx'], float)
    iwfr = np.asarray(arrays['iwfr'], float)

    mask = np.ones(len(ssm), dtype=bool)

    if qc.get('ssm_min') is not None:
        mask &= ssm >= qc['ssm_min']
    if qc.get('ssm_max') is not None:
        mask &= ssm <= qc['ssm_max']
    if qc.get('smpf_ok') is not None:
        mask &= np.isin(smpf, list(qc['smpf_ok']))
    if qc.get('smcf_ok') is not None:
        mask &= np.isin(smcf, list(qc['smcf_ok']))
    if qc.get('tpcx_max') is not None:
        mask &= (tpcx >= 0) & (tpcx <= qc['tpcx_max'])
    if qc.get('iwfr_max') is not None:
        mask &= (iwfr >= 0) & (iwfr <= qc['iwfr_max'])

    return mask


def apply_h121_qc(arrays, qc=None):
    """Return boolean mask of obs passing QC.

    arrays: dict with keys matching H121 variable names
            (ssm, sf, pf, cf, tc, wf, subsfc).
    qc:     QC config dict; defaults to QC_DEFAULT_H121.
    """
    if qc is None:
        qc = QC_DEFAULT_H121

    ssm   = np.asarray(arrays['ssm'],   float)
    sf    = np.asarray(arrays.get('sf', np.zeros(len(ssm))), float)
    pf    = np.asarray(arrays['pf'],    float)
    cf    = np.asarray(arrays.get('cf', np.zeros(len(ssm))), float)
    tc    = np.asarray(arrays['tc'],    float)
    wf    = np.asarray(arrays['wf'],    float)
    subsfc = np.asarray(arrays['subsfc'], float)
    sens  = np.asarray(arrays.get('sens', np.full(len(ssm), np.inf)), float)
    bsflag = np.asarray(arrays.get('bsflag', np.zeros(len(ssm))), float)

    mask = np.ones(len(ssm), dtype=bool)

    if qc.get('ssm_min') is not None:
        mask &= ssm >= qc['ssm_min']
    if qc.get('ssm_max') is not None:
        mask &= ssm <= qc['ssm_max']
    if qc.get('sf_water_bit') is not None:
        mask &= (sf.astype(np.int64) & int(qc['sf_water_bit'])) == 0
    if qc.get('pf_bad_bits') is not None:
        mask &= (pf.astype(np.int64) & int(qc['pf_bad_bits'])) == 0
    elif qc.get('pf_ok') is not None:
        mask &= np.isin(pf, list(qc['pf_ok']))
    if qc.get('cf_ok') is not None:
        mask &= np.isin(cf, list(qc['cf_ok']))
    if qc.get('tc_max') is not None:
        mask &= (tc >= 0) & (tc < qc['tc_max'])
    if qc.get('wf_max') is not None:
        mask &= (wf >= 0) & (wf < qc['wf_max'])
    if qc.get('subsfc_max') is not None:
        mask &= subsfc < qc['subsfc_max']
    if qc.get('sens_min') is not None:
        mask &= sens > qc['sens_min']
    if qc.get('bsflag_bad_bits') is not None:
        mask &= (bsflag.astype(np.int64) & int(qc['bsflag_bad_bits'])) == 0

    return mask
