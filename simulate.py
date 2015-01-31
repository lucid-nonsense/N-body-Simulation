import numpy as np
import os
import datetime
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

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	from algorithms import forward_euler
	N = 2
	S = Simulation(forward_euler, dt=0.01, T=10)
	S.initialise_particles(N)
	S.r[1] = [1, 0, 0]
	S.v[1] = [0, 0.5, 0]
	S.static[0] = True

	rs, vs = S.open_files('')
	S.start(0)
	S.close_files()
	r = read_output(rs, N)[:, 1, :2]
	x,y = r[:, 0], r[:, 1]
	print x.shape, y.shape
	plt.scatter(x, y, marker='.')
	plt.show()