"""Independent radial integral and analytic limiting cases in free space."""

import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree
from scipy.stats import ncx2

from mobile_water import RT, Site, penalty

sigma = 0.47
center = np.zeros(3)
ligand = cKDTree([[2.35, 0, 0]])
b_iso = sigma**2 * 8 * np.pi**2
fixed, collapsed = Site.prepare(center, 0, None).evaluate(ligand)
assert np.isclose(fixed, collapsed, atol=1e-12)
away = cKDTree([[100, 0, 0]])
assert Site.prepare(center, b_iso, None).evaluate(away) == (0.0, 0.0)
assert np.isclose(penalty(0), -RT * np.log(0.25))
noncentrality = (2.35 / sigma) ** 2


def cdf(r):
    return ncx2.cdf((r / sigma) ** 2, 3, noncentrality)


# Integrate by parts: E[A(R)] = 1 - integral CDF(r) A'(r) dr.
def integrand(r):
    t = (r - 2.0) / 0.6
    return cdf(r) * 6 * t * (1 - t) / 0.6


expected = penalty(1 - quad(integrand, 2.0, 2.6, epsabs=1e-12)[0])
for order, tolerance in [(11, 0.006), (17, 0.004)]:
    _, estimate = Site.prepare(center, b_iso, None, order).evaluate(ligand)
    print(f"order {order}: {estimate:.8f}, radial integral {expected:.8f}")
    assert abs(estimate - expected) < tolerance
print("collapsed, saturation, far-field and independent radial integral checks passed")
