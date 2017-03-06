# -*- coding: utf-8 -*-

# test working of chi square distribution with ~30 points along a straight line

from numpy import random, array, linspace, diff, ones_like
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys
import csv
from rosenbluth import rosenbluth, partition

#x = epsilon, y = reduced, fit_function = lambda x: a*x + b
def chi_square(x, y, error, fit_function):

	""" Calculate the chi-squared statistic for a set of reduced cross section and epsilon values.

	Arguments:
		x 				(list) - Epsilon values
		y 				(list) - Reduced cross section values
		error 			(list) - Errors associated with reduced cross sections
		fit_function	(lambda) - Function to calculate expected reduced cross section value given epsilon

	Returns:
		sum				(float) - The calculated chi-square statistic

	"""

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
		with open(file, 'r') as f:
			data = csv.reader(f, delimiter=" ")
			# skip header row
			next(data)
			for row in data:
				# allow comments
				if not row or "".join(row)[0] == "#":
					# only normalize 1.6 GeV data or 8 GeV where angle ~90
					if len(row) > 1 and ("1.6" in row[1] or "norm" in row[1]):
						next_normalized= True
					else:
						next_normalized = False
					continue
				

				q_squared.append(float(row[0]))
				energies.append(float(row[1]))
				thetas.append(float(row[2]))

				# specify a normalization factor as a command line argument
				if next_normalized and len(sys.argv) > 2:
					norm = float(sys.argv[2]) # 0.958, etc.
				else:
					# don't normalize
					norm = 1

				cross_sections.append(float(row[3])*norm)
				uncertainties.append(float(row[4])*norm)
				if len(row) == 6:
					energies_prime.append(float(row[5]))
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
actual_chi_squared = 0.0

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

	#looks redundant with next for loop but it has to go through once to get all the epsilons and reduced for the fit
	for j in range(len(q_squared)):
		eps_original = []
		red_original = []
		eps = []
		red = []
		chi_squares = []
		samples = 10000
		actual_chi_squared = 0.0

		for i in range(len(q_squared[j])):
			q2 = q_squared[j][i]
			energy = energies[j][i]
			theta = thetas[j][i]
			cross_section = cross_sections[j][i]
			error = total_errors[j][i]
			epsilon, tau, reduced, error = rosenbluth(q2, energy, theta, cross_section, error)
			total_errors[j][i] = error
			error = total_errors[j][i]
			eps_original.append(epsilon)
			red_original.append(reduced)
			# print("(" + str(epsilon) + ", " + str(reduced) + ")" + str(error))

		#try randomizing from first fit function (2 points at first epsilon, etc.)
		reg = stats.linregress(eps_original, red_original)
		actual_chi_squared = chi_square(eps_original, red_original, total_errors[j], lambda x: reg[0]*x+reg[1])
		# print(eps_original, red_original)

		#distribution using each of equal-epsilon points, based on fit of only three points
		#remove one at a time the equal-angle points, then uncomment below and run

		# eps_original.insert(0,eps_original[0])
		# red_original.insert(0,red_original[0])
		# total_errors[0].insert(0, total_errors[0][0])
		# q_squared[0].append(1.75)

		#-----------------------------------------------------

		eps,red=[],[]
		for i in range(len(q_squared[j])):
			epsilon = eps_original[i]
			# print(epsilon, red_original[i], total_errors[j][i])
			eps.append([epsilon]*samples)
			red.append(random.normal(reg[0]*epsilon+reg[1], total_errors[j][i], samples))

		eps = array(eps)
		red = array(red)

		#reset monte carlo points log file
		with open('monte_carlo_points.csv','w') as out:
				writer = csv.writer(out, delimiter=' ')

		for i in range(samples):

			eps_samples = eps[:,i]
			red_samples = red[:,i]
			# for logging the points 
			
			if i < 100:
				with open('monte_carlo_points.csv', 'a') as out:
						writer = csv.writer(out, delimiter=' ')
						writer.writerow(["({},{})".format(eps_samples[c], red_samples[c]) for c in range(len(eps_samples))])
			reg = stats.linregress(eps_samples, red_samples)
			chi_squares.append(chi_square(eps_samples, red_samples, total_errors[j], lambda x: reg[0]*x+reg[1]))

		#plot chi square distribution
		plt.delaxes()
		fig = plt.figure(figsize=(12,9))
		rc('font',**{'family':'serif'})
		ax = plt.subplot(111)
		ax.set_xlabel(r'$\chi^2$', fontsize=30)
		ax.set_ylabel(r'Probability', fontsize=30)
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		ax.annotate('$Q^2 = %.3f$ \n $df = %d$ \n $\chi^2 = %.3f$' % (q2, len(q_squared[j])-1, actual_chi_squared), xy=(0, 0), xycoords='data',
		xytext=(0.6, 0.7), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=30)

		area = 0.0
		if test:
			vals, bins, _ = plt.hist(chi_squares, bins=100, facecolor='blue', alpha=1.0, linewidth=2, histtype="stepfilled")

		else:
			vals, bins, _ = plt.hist(chi_squares, bins=100, facecolor='gray', alpha=1.0, linewidth=2, histtype="stepfilled", normed=True)
			
			# total area of histogram
			area = sum(diff(bins)*vals)
			# print(area)
			
			# area of histogram to actual chi squared
			lastBin = 0
			while bins[lastBin] <= actual_chi_squared:
				lastBin+=1
				if lastBin == len(bins):
					lastBin-=1
					break

			# print percentage from partial/total area
			partialArea = sum(diff(bins[:lastBin])*vals[:lastBin-1])
			print(actual_chi_squared)
			print("{:.4f} %".format(100*partialArea/area))


		fig.tight_layout()
		plt.show()



