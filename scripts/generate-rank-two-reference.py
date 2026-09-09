"""Independent angular quadrature at 80 decimal digits (mpmath).

Factor out the peak before quadrature so its tolerance also resolves tiny tails.
"""
import mpmath as mp
mp.mp.dps = 80
for q in (8, 400, 2000):
    q = mp.mpf(q)
    integral = mp.quad(lambda a: mp.exp(-q/(2*(mp.cos(a)**2+2*mp.sin(a)**2))+q/4),
                       [0, mp.pi/4, mp.pi/2])
    print(q, mp.exp(-q/4)*2*integral/mp.pi)
