# -*- coding: utf-8 -*-

# test working of chi square distribution with ~30 points along a straight line

from numpy import deg2rad, linalg, tan, sin, cos, random, array, std, mean, where, linspace
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
		vals.append((fit_function(x[i]) - y[i])**2/error[i]**2)
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
eps_original = []
red_original = []
eps = []
red = []
chi_squares = []
samples = 10000

# TEST: make 30 points along a line, then randomize
if '-t' in sys.argv:
	test = True
else:
	test = False
if (test):
	a = random.random()
	b = random.random()
	x = linspace(1,10,num=30)
	y = a*x + b
	total_errors = [[0.002]*samples]*30
	for i in range(len(total_errors)):
		eps.append([x[i]]*samples)
		red.append(random.normal(y[i], total_errors[i][0], samples))
	eps=array(eps)
	red=array(red)

else:
	#just do one Q^2 for now
	for i in range(len(q_squared[0])):
		q2 = q_squared[0][i]
		energy = energies[0][i]
		theta = thetas[0][i]
		cross_section = cross_sections[0][i]
		error = total_errors[0][i]
		epsilon, tau, reduced, error = rosenbluth(q2, energy, theta, cross_section, error)
		eps_original.append(epsilon)
		red_original.append(reduced)
		eps.append([epsilon]*samples)
		red.append(random.normal(reduced, error, samples))

	#try randomizing from first fit function (2 points at first epsilon, etc.)
	reg = stats.linregress(eps_original, red_original)
	eps,red=[],[]
	for i in range(len(q_squared[0])):
		epsilon = eps_original[i]
		print(epsilon)
		eps.append([epsilon]*samples)
		red.append(random.normal(reg[0]*epsilon+reg[1], total_errors[0][i], samples))


eps = array(eps)
red = array(red)


for i in range(samples):

	eps_samples = eps[:,i]
	red_samples = red[:,i]
	reg = stats.linregress(eps_samples, red_samples)
	chi_squares.append(chi_square(eps_samples, red_samples, total_errors[0], lambda x: reg[0]*x+reg[1]))

#plot chi square distribution
plt.delaxes()
plt.figure(figsize=(12,9))
ax = plt.subplot(111)
ax.set_xlabel(r'$\chi^2$', fontsize=30)
ax.set_ylabel(r'Frequency', fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if test:
	plt.hist(chi_squares, bins=100, normed=1, facecolor='gray', alpha=0.5, linewidth=2, histtype="stepfilled")
else:
	plt.hist(chi_squares, bins=100, normed=1, facecolor='gray', alpha=0.5, linewidth=2, histtype="stepfilled")

plt.show()

