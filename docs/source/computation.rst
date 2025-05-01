Computation
===========

Here are some numerical computing strategies:

- The algorithm only does conversion between wavelet coefficients and sqeuence in frequency domain. 
  (I)FFT are used to transform between time and frequency domain.
- Outside the transition band, value of Meyer kernel is either :math:`1` or :math:`0`. 
  So, calculate the transition band only.

Transform
---------

.. math::
    w_{mn}=A\frac{(-1)^{mn}}2\sum_{k=0}^{N-1}\left(\tilde{x}[k-mN_t]\pm\tilde{x}[k+mN_t]\right)
    \tilde{\phi}[k]\mathrm{e}^{\mathrm{i}\pi nk/N_t}.

Summation regarding :math:`k` contains :math:`2N_t` items actually and can be calculated with IFFT.
For real sequence :math:`x[j]`, the above formula is reduced to

.. math::
    \begin{gather}
    w_{mn} = A(-1)^{mn} \mathcal{O}_{mn}\sum_{k}\tilde{x}[k-mN_t]\tilde{\phi}[k]\mathrm{e}^{\mathrm{i}\pi nk/N_t}. \\
    O_{mn} = \mathcal{R},\qquad\text{if}\; m+n\;\text{even}\quad \text{else}\; \mathcal{I}.
    \end{gather}



Inverse transform
-----------------

.. math::
    \tilde{x}[k] = \frac12 \sum_{m}\left(\tilde{\phi}[k-mN_t]\pm\tilde{\phi}[k+mN_t]\right)
    \sum_{n}A_{mn}w_{mn}\mathrm{e}^{-\mathrm{i}\pi nk/N_t}.

Summation regarding :math:`n` can be calculated with FFT. For real sequence :math:`x[j]`, 
items with negative :math:`k` are redudant. Therefore, :math:`\tilde{\phi}[k+mN_t]` is ignored.

.. attention::
    Special treatment for :math:`m=0` and :math:`N_f/2`.
