""" Wilson-Daubechies-Mayer wavelet transform.
"""
import numpy as np
from scipy.special import betainc
from math import ceil

# Constants.
PI = 3.141592653589793238462643383279502884
SQRT2 = 1.414213562373095048801688724209698079

def _rolloff_idx(flat_width, trans_width):
    """Computer start and end index of Meyer roll-off edge.
    
    Parameters
    ----------
    flat_width : float
        Width of half flat interval.
    trans_width : float
        Width if transition interval.

    Returns
    -------
    (sidx, eidx) : tuple
        Start and end indices of transition interval.
    """
    sidx = ceil(flat_width)
    eidx = ceil(flat_width + trans_width)
    return sidx, eidx


def _gen_rolloff(flat_width, trans_width, steepness):
    """Compute Meyer roll-off edge.
    """
    sidx, eidx = _rolloff_idx(flat_width, trans_width)
    idxs = np.arange(sidx, eidx, dtype=int)

    edge = np.cos(PI / 2 *
            betainc(steepness,
                    steepness,
                    (idxs - flat_width) / trans_width
            )
        )

    return (sidx, eidx), edge

def rwdm(x, nt, nf, trans_width, steepness, axis=-1):   #TODO: delete nt or nf
    """Compute WDM transform for real input.
    
    Parameters
    ----------
    x : array_like
        Input array.
    nt : int
        N_t.
    nf : int
        N_f.
    trans_width : float
        Width of transition interval.
    steepness : float
        Steepness of transition interval.
    axis : int
        Axis to be transformed.

    Returns
    -------
    w : np.ndarray
        WDM coefficients, stored in the last two axes.
    """
    if nf & 1:
        raise ValueError("nf must be even.")
    if trans_width < 1.:
        trans_width = 1.
    flat_width = (nt - trans_width) / 2

    x = np.asarray(x)
    if axis != -1:
        x = np.moveaxis(x, axis, -1)
    
    (trans_sidx, trans_eidx), edge = _gen_rolloff(flat_width, trans_width, steepness)

    # Transform.
    trans_sidx_r = trans_sidx + nt
    trans_eidx_r = trans_eidx + nt
    trans_sidx_l = -trans_eidx + nt + 1
    trans_eidx_l = -trans_sidx + nt + 1

    xf = np.fft.rfft(x, norm="ortho")
    ans_shape = xf.shape[:-1] + (nf//2, 2*nt)
    ans = np.empty(ans_shape, dtype=float)

    ## m == 0
    xf_ = xf[..., nt::-1].conj()
    xf_[..., trans_sidx_l:trans_eidx_l] *= edge[::-1]
    xf_[..., :trans_sidx_l] = 0.

    ans0 = np.fft.irfft(xf_, norm="ortho")
    ans[..., 0, ::2] = ans0[..., ::2] * SQRT2

    ## m in [1, nf//2)
    for m in range(1, nf//2):  
        k_sidx = (m - 1) * nt
        k_eidx = k_sidx + 2 * nt   
        xf_ = xf[..., k_sidx:k_eidx].copy()
        xf_[..., trans_sidx_r:trans_eidx_r] *= edge
        xf_[..., trans_sidx_l:trans_eidx_l] *= edge[::-1]   # need edge_rev cache?
        xf_[..., trans_eidx_r:] = 0.
        xf_[..., :trans_sidx_l] = 0.

        ans_m = np.fft.ifft(xf_, norm="ortho") 
        neo = 1 ^ (m & 1)
        ans_m[..., neo::2].real = ans_m[..., neo::2].imag 
        ans_m = 2 * ans_m.real
        if neo == 0:
            ans_m[..., ::2] *= -1

        ans[..., m, :] = ans_m

    ## m == nf//2
    xf_ = xf[..., -(nt+1):].copy()
    xf_[..., trans_sidx_l:trans_eidx_l] *= edge[::-1]     
    xf_[..., :trans_sidx_l] = 0.

    ans1 = np.fft.irfft(xf_, norm="ortho")
    eo = (nf//2) & 1
    ans[..., 0, 1::2] = ans1[..., eo::2] * SQRT2

    return ans  #TODO: restore moved axis

def irwdm(w, trans_width, steepness):  #TODO: tran_width check, axis. 
    """Compute the inverse of real WDM transform ``rwdm``.
    
    Parameters
    ----------
    w : array_like
        Input WDM coefficients, stored in the last two axis.
    trans_width : float
        Width of transition interval.
    steepness : float
        Steepness of transition interval.
    axis : int

    Returns
    -------
    np.ndarray
        Inverse of w.
    """
    w = np.asarray(w)
    w_shape = w.shape
    nt = (w_shape[-1] + 1) // 2
    nf = w_shape[-2] * 2
    flat_width = (nt - trans_width) / 2
    n = nf * nt
    xf_shape = w_shape[:-2] + (n//2+1,)

    xf = np.zeros(xf_shape, dtype=complex)
    (trans_sidx, trans_eidx), edge = _gen_rolloff(flat_width, trans_width, steepness)
    trans_sidx_r = trans_sidx + nt
    trans_eidx_r = trans_eidx + nt
    trans_sidx_l = -trans_eidx + nt + 1
    trans_eidx_l = -trans_sidx + nt + 1

    # m == 0
    w_ = np.fft.fft(w[..., 0, ::2], norm="ortho")
    xf[..., :trans_sidx] = w_[..., :trans_sidx]
    xf[..., trans_sidx:trans_eidx] += w_[..., trans_sidx:trans_eidx] * edge

    # m in [1, nf//2)
    for m in range(1, nf//2):
        k_shift = (m - 1) * nt
        k_sidx_r = trans_sidx_r + k_shift
        k_eidx_r = trans_eidx_r + k_shift
        k_sidx_l = trans_sidx_l + k_shift
        k_eidx_l = trans_eidx_l + k_shift 

        neo = 1 ^ (m & 1)
        w_ = np.array(w[..., m, :], dtype=complex)
        w_[..., neo::2] *= 1j
        if neo == 0:
            w_[..., ::2] *= -1
        w_ = np.fft.fft(w_, norm="ortho")
        xf[..., k_sidx_l:k_eidx_l] += w_[..., trans_sidx_l:trans_eidx_l] * edge[::-1]
        xf[..., k_eidx_l:k_sidx_r]  = w_[..., trans_eidx_l:trans_sidx_r]
        xf[..., k_sidx_r:k_eidx_r] += w_[..., trans_sidx_r:trans_eidx_r] * edge
    
    # m = nf//2
    m = nf // 2
    k_shift = (m - 1) * nt
    k_sidx_l = trans_sidx_l + k_shift
    k_eidx_l = trans_eidx_l + k_shift

    neo = 1 ^ (m & 1)
    w_ = np.fft.fft(w[..., 0, 1::2], 
                    norm="ortho") * np.exp(-1j * PI * np.arange(nt) / nt)
    xf[..., k_sidx_l:k_eidx_l] += w_[..., trans_sidx_l:trans_eidx_l] * edge[::-1]
    xf[..., k_eidx_l:-1] = w_[..., trans_eidx_l:]
    xf[..., -1] = -w_[..., 0]

    ans = np.fft.irfft(xf, norm="ortho")
    # ans = xf
    return ans
    
# Debug
# nt = 6
# nf = 6
# n = nt * nf
# x = np.arange(2*n).reshape(2,-1)
# trans_width = 2

# w = rwdm(x, nt, nf, trans_width, 2)
# xf_ = np.fft.rfft(x, norm="ortho")
# xf = irwdm(w, trans_width, 2)


# pass
