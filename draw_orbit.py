from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

def get_ellipse(v0, ra, mu):
	"""v0 is starting velocity, ra is apogee distance, mu = GM"""
	a_inv = (2 / ra) - (v0 * v0 / mu)# vis viva
	a = 1 / a_inv
	T = np.sqrt(4 * np.pi * np.pi * a * a * a / mu)
	e = (ra / a) - 1
	return a, abs(e), T

def iterate_E(_M):
	lambda _M: newton(lambda E: E - (e*np.sin(E)) - _M, _M, lambda E: 1 - (e * np.cos(E)))

def plot_orbit(t_array, mu, v0, pos0, pos1):
	"""posn are arrays of (x, y)"""
	r0 = np.linalg.norm(pos1 - pos0)
	a, e, T = get_ellipse(v0, r0, mu)
	M = 2 * np.pi * t_array / T
	iterate =  np.vectorize()
	E = iterate(M)
	r_array = a * (1 - (e * np.cos(E)))
	y = a * np.sqrt(1 - (e * e)) * np.sin(E) # rsin(theta) = y
	x = np.sqrt(np.square(r_array) - np.square(y))
	return x, y

if __name__ == '__main__':
	r = 1.
	v = 0.5
	print get_ellipse(v, r, 1)