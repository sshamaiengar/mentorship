# -*- coding: utf-8 -*-

from numpy import deg2rad, linalg, tan, sin, cos, random, array, std, mean, where
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys
import csv
from rosenbluth import rosenbluth, partition

#x = epsilon, y = reduced, fit_function = lambda x: a*x + b
def chi_square(x, y, error, fit_function):
	vals = []
	for i in range(len(x)):
		vals.append((y[i] - fit_function(x[i]))**2/error[i])
	return sum(vals)

if len(sys.argv) > 1:
	q_squared = []
	energies = []
	thetas = []
	cross_sections = []
	uncertainties = []

	file = sys.argv[1]
	try:
		with open(file, 'rb') as f:
			data = csv.reader(f, delimiter=" ")
			# skip header row
			next(data)
			for row in data:
				q_squared.append(float(row[0]))
				energies.append(float(row[1]))
				thetas.append(float(row[2]))
				cross_sections.append(float(row[3]))
				uncertainties.append(float(row[4]))
	except IOError as e:
		print(e)
else:
	q_squared = [4.0,4.0,4.0,4.0,4.0]
	energies = [3.4, 3.956, 4.507, 5.507, 9.8]
	thetas = [57.572, 43.707, 35.592, 26.823,13.248]
	cross_sections = [1.297e-2, 2.77e-2, 4.929e-2, 1.023e-1, 6.18e-1]
	uncertainties = [2.243e-4,4.407e-4,7.853e-4,1.370e-3,8.073e-3]

q_squared, energies, thetas, cross_sections, total_errors = partition(q_squared, energies, thetas, cross_sections, uncertainties)
eps = []
red = []
chi_squares = []
samples = 10000

#just do one Q^2 for now
for i in range(len(q_squared[0])):
	q2 = q_squared[0][i]
	energy = energies[0][i]
	theta = thetas[0][i]
	cross_section = cross_sections[0][i]
	error = total_errors[0][i]

	epsilon, tau, reduced = rosenbluth(q2, energy, theta, cross_section)

	eps.append([epsilon]*samples)
	red.append(random.normal(reduced, 0.02, samples))

eps = array(eps)
red = array(red)

for i in range(samples):

	eps_samples = eps[:,i]
	red_samples = red[:,i]
	print(eps_samples)
	reg = stats.linregress(eps_samples, red_samples)
	chi_squares.append(chi_square(eps_samples, red_samples, total_errors[0], lambda x: reg[0]*x+reg[1]))

#plot chi square distribution
plt.delaxes()
plt.figure(figsize=(12,9))
ax = plt.subplot(111)
ax.set_xlabel(r'$\chi^2$', fontsize=30)
ax.set_ylabel(r'Frequency', fontsize=30)
plt.hist(chi_squares, bins=100, normed=1, facecolor='gray', alpha=0.5, linewidth=2, histtype="stepfilled")

plt.show()

