import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math


def spherical_asperities(L, N, R, surface_outfile=''):
	'''
	Creates gaussian surfaces with spherical asperities of given radius
	:param L: Length of the square surface
	:param R: Radius of the asperities
	:param N: number of intervals
	:return: writes an point cloud data file and returns surface data as an numpy array
	'''
	delta = L / N
	surface_data = np.zeros((N + 1, N + 1))
	height = lambda x, y, xc, yc: np.sqrt(R ** 2 - (x - xc) ** 2 - (y - yc) ** 2) \
		if (R ** 2 - (x - xc) ** 2 - (y - yc) ** 2) > 0 else 0
	lines = []
	for i in range(surface_data.shape[0]):
		for j in range(surface_data.shape[1]):
			x = delta * i
			y = delta * j
			xc = 2 * R * (math.ceil((delta * i) / (2 * R)) - 1) + R if L >= R else L/2.0
			yc = 2 * R * (math.ceil((delta * j) / (2 * R)) - 1) + R if L >= R else L/2.0
			surface_data[i, j] = height(x, y, xc, yc)
			# if math.isnan(surface_data[i, j]) or surface_data[i, j] > 5:
			# print("X = " + str(x) + ", Y = " + str(y) + " ,XC = " + str(xc) + " ,YC = " + str(yc) + ",Z = " + str(surface_data[i, j]))
			lines.append(str(delta * i) + '\t' + str(delta * j) + '\t' + str(surface_data[i, j]) + '\n')

	if surface_outfile:
		with open(surface_outfile.as_posix(), 'w') as file:
			file.writelines(lines)

	return surface_data


def planar(L, N, surface_outfile=''):
	'''
	creates a planar surface of length L and N number of divisions
	:param L: Length of the square surface
	:param N: number of intervals
	:return: write a asc file at surface_outfile string path
	'''

	delta = L / N
	lines = []
	for i in range(N + 1):
		for j in range(N + 1):
			x = delta * i
			y = delta * j
			lines.append(str(delta * i) + '\t' + str(delta * j) + '\t0.0\n')

	if surface_outfile:
		with open(surface_outfile.as_posix(), 'w') as file:
			file.writelines(lines)


def plot3D(surface_file):
	with open(surface_file, 'r') as file:
		lines = file.readlines()

	lines = [[float(x) for x in line.split()] for line in lines]
	X = list(set([line[0] for line in lines]))
	Y = list(set([line[1] for line in lines]))

	Z = np.zeros((len(X), len(Y)))
	for i, x in enumerate(X):
		for j, y in enumerate(Y):
			for line in lines:
				if x == line[0] and y == line[1]:
					Z[i, j] = line[2]
					break

	X, Y = np.meshgrid(X, Y)
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.plot_surface(X, Y, Z, cmap=cm.viridis, linewidth=0, antialiased=False)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	plt.show()

def paper_plot(surface_file):

	with open(surface_file, 'r') as file:
		lines = file.readlines()

	lines = [[float(x) for x in line.split()] for line in lines]
	X = list(set([line[0] for line in lines]))
	Y = list(set([line[1] for line in lines]))

	Z = np.zeros((len(X), len(Y)))
	for i, x in enumerate(X):
		for j, y in enumerate(Y):
			for line in lines:
				if x == line[0] and y == line[1]:
					Z[i, j] = line[2]
					break

	X, Y = np.meshgrid(X, Y)
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.scatter(X, Y, Z)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	plt.show()