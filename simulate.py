import numpy as np
import os
import datetime
from itertools import combinations
G = 6.67e-11

class Simulation(object):
	"""Contains the state of a simulation with methods for evolving it"""
	def __init__(self, evolver, dt, T, N=None, r=None, v=None, m=None, a=None, static=None, dim=3):
		super(Simulation, self).__init__()
		self.dt = dt
		self.T = T
		self._evolver = evolver
		self.dim = dim
		if any([i is not None for i in (r,v,m,a,N,static)]):
			initialise_particles(N, r, v, a, m, N, static)

	def initialise_particles(self, N=None, r=None, v=None, a=None, m=None, static=None):
		"""r,v,m are 3xN numpy arrays"""
		if N is not None:
			shape = (N, self.dim)
			self.r = np.zeros(shape)
			self.v = self.r.copy()
			self.a = self.r.copy()
			self.m = np.ones(N)
			self.N = N
			self.static = [False]*N
		else:
			if not (r.shape == v.shape == m.shape) or r.shape[0] != self.dim:
				raise TypeError('Input particle properties must be of the same length and of {}-dimensions'.format(self.dim))
			self.r = r
			self.v = v
			self.m = m
			self.a = a
			self.static = static
			self.N = r.shape[0]

	def evolve(self):
		self._evolver(self.r, self.v, self.m, self.a, self.dt, self.static)

	def start(self, record_int=None):
		for t in xrange(0, int(self.T/self.dt)):
			self.evolve()
			if record_int is not None:
				if record_int == 0:
					self.record()
				elif not t % record_int:
					self.record()


	def append_file(self, f, x):
		np.savetxt(f, x)
	
	def open_files(self, d):
		date = datetime.datetime.now().strftime('%y-%m-%d-%H-%M-%S')
		base = '{}-{}-{}-{}'.format(self.dt, self.T, self.N, date)
		rf = os.path.join(d, base+'-r.dat')
		vf = os.path.join(d, base+'-v.dat')
		self.record_r, self.record_v = open(rf, 'a'), open(vf, 'a')
		self.rf, self.vf = rf, vf
		return rf, vf

	def close_files(self):
		self.record_r.close()
		self.record_v.close()

	def record(self):
		self.append_file(self.record_r, self.r)
		self.append_file(self.record_v, self.v)

def read_output(x, N):
	_T = np.loadtxt(x)
	return _T.reshape([_T.shape[0]/N,N,_T.shape[1]]) # time, particle, vector

def get_kinetic(v_array, m_list):
	return (np.apply_along_axis(lambda x: (x**2).sum(), 2, v_array) * m_list).sum(axis=1) / 2

def get_potential(r_array, m_list):
	u_array = np.zeros(r_array.shape[0])
	for i in xrange(r_array.shape[0]):
		pairs = combinations(xrange(r_array.shape[1]), 2)
		U = 0
		for p in pairs:
			dist = r_array[i, p[1], :] - r_array[i, p[0], :]
			U += -1 * m_list[p[0]] * m_list[p[1]] / np.linalg.norm(dist)
		u_array[i] = U
	return u_array
		
def plot_xyprojection(sim_obj):
	r = read_output(sim_obj.rf, sim_obj.N)
	v = read_output(sim_obj.vf, sim_obj.N)
	K = get_kinetic(v, sim_obj.m)
	U = get_potential(r, sim_obj.m)
	E = K + U

	f = plt.figure()
	f.set_facecolor('w')
	ax = f.add_subplot(211)
	ax.plot((E - E[0]) / E[0], 'b-')
	ax.set_ylabel(r'Relative Energy Error ($\frac{E - E_0}{E_0}$)')
	ax.set_xlabel(r'Time $t$')
	
	if sim_obj.N == 2:
		_r = r[:, 1, :]
		ax2 = ax.twinx()
		dist = np.linalg.norm(_r, axis=1)
		ax2.plot(dist, 'r-')
		ax2.set_ylabel(r'Separation $r$')

	ax3 = f.add_subplot(212)
	colors = list('bgrcmykw')
	print r.shape
	for p in xrange(r.shape[1]):
		x,y = r[:, p, 0], r[:, p, 1]
		ax3.scatter(x, y, marker='.', alpha=0.5, s=S.m[p]*100, color=colors[p])
	

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	from algorithms import forward_euler
	S = Simulation(forward_euler, dt=0.0001, T=2)
	S.initialise_particles(N=2)
	S.r[1] = [1, 0, 0]
	S.v[1] = [0, 0.5, 0]

	S.static[0] = True

	rs, vs = S.open_files('output')
	S.start(2)
	S.close_files()
	plot_xyprojection(S)
	plt.show()