# Copyright 2026 Roby Joehanes; GNU GPL version 3.
"""Synthetic independent references: chi-square(1) + epsilon*chi-square(rank-1).

Requires mpmath. Uses positive conditional convolution, not characteristic
inversion or the production truncation bound. Run from the repository root.
"""
import csv
from mpmath import mp


def survival(q, epsilon, rank, digits):
    with mp.workdps(digits):
        q, epsilon = mp.mpf(q), mp.mpf(epsilon)
        shape = mp.mpf(rank - 1) / 2
        normalizer = mp.power(2, shape) * mp.gamma(shape)
        leading = mp.erfc(mp.sqrt(q / 2))

        def integrand(x):
            if x == 0:
                return mp.mpf(0)
            tail = mp.erfc(mp.sqrt((q - epsilon*x)/2))
            density = mp.exp((shape-1)*mp.log(x)-x/2)/normalizer
            return tail/leading*density

        limit = q/epsilon
        knots = [mp.mpf(x) for x in [0, 40, 100, rank-1, 250, 400, 800, 1600] if x < limit]
        knots.append(limit)
        value = leading*mp.quad(integrand, sorted(set(knots)))
        value += mp.gammainc(shape, limit/2, mp.inf)/mp.gamma(shape)
        return +value


if __name__ == '__main__':
    rows = []
    for q, epsilon, rank in [('500', '0.0001', 160), ('700', '0.0001', 160),
                             ('500', '0.000001', 160), ('1400', '0.0001', 160)]:
        first = survival(q, epsilon, rank, 70)
        second = survival(q, epsilon, rank, 110)
        with mp.workdps(70):
            assert abs(first-second)/second < mp.mpf('1e-50')
            row = [q, epsilon, rank, mp.nstr(second, 60)]
        rows.append(row)
        print(','.join(map(str, row)), flush=True)
    with open('test/resources/tail-remainder-reference.csv', 'w', newline='', encoding='utf8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['q', 'epsilon', 'rank', 'survival'])
        writer.writerows(rows)
