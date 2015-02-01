from __future__ import division
import numpy as np
from itertools import combinations, tee

def forward_euler(R, V, M, A, dt, static):
	pairsA, pairsVR = tee(combinations(xrange(len(R)), 2))
	for i, j in pairsA:
		vect_dist = R[j] - R[i]
		r_norm = np.linalg.norm(vect_dist)
		dr = vect_dist / r_norm / r_norm / r_norm
		if not static[i]:
			A[i] = M[j] * dr
		if not static[j]:
			A[j] = -1 * M[i] * dr
		
	for i, j in pairsVR:
		if not static[i]:
			R[i] += V[i] * dt
			V[i] += A[i] * dt
		if not static[j]:
			R[j] += V[j] * dt
			V[j] += A[j] * dt

def main():
	r = np.array([[0,0,0], [1,1,0]], dtype='f')
	v = np.zeros_like(r)
	m = np.array([1, 1])
	static = [True, False]
	a = np.zeros_like(r)
	forward_euler(r, v, m, a, 1, static)
	forward_euler(r, v, m, a, 1, static)
	print r

if __name__ == '__main__':
	main()
