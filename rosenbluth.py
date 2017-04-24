# -*- coding: utf-8 -*-

from numpy import deg2rad, linalg, tan, sin, cos, random, array, std, mean, where
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys
import csv


# constants
fsc = 0.00729735256
mass = 0.938 		# proton mass in GeV
pmm = 2.79 			# proton magnetic moment


def rosenbluth(q_squared, energy, theta, cross_section, error, energy_prime=None):

	""" Takes the constants used in the Rosenbluth formula and returns the computed variables

	Args:
		q_squared 		(float)	- The value of Q^2
		energy 			(float)	- The beam energy in GeV
		theta 			(float)	- The spectrometer angle in degrees
		cross_section 	(float)	- The differential cross section in (nb/sr)
		error 			(float)	- The total error
		energy_prime 	(float)	- The (final) beam energy in GeV

	Returns:
		(tuple) - 
			epsilon 		(float)	- The computed value of epsilon
			tau 			(float)	- The computed value of tau
			reduced 		(float)	- The computed value of the reduced cross section
			error 			(float) - The computed error associated with the reduced cross section

	"""
	
	tau = q_squared / 4 / mass ** 2
	epsilon = (1 + 2 * (1 + tau) * tan(deg2rad(theta)/2) ** 2) ** -1
	eta = 1 + (energy / mass) * (1 - cos(deg2rad(theta)))

	mott_cross_section = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * 0.389 * 1e6 * eta ** -1
	reduced = cross_section/mott_cross_section * epsilon * (1 + tau)
	error = error/mott_cross_section * epsilon * (1 + tau)
	return (epsilon, tau, reduced, error)

def form_factors(epsilon, reduced):

	""" Performs a linear regression of (epsilon, reduced cross section) data 
		to get the proton electric and magnetic form factor values squared (Rosenbluth separation)

	Args:
		epsilon 		(list) - The values of epsilon
		reduced			(list) - The values of the reduced cross sections

	Returns:
		(tuple) - 
			(float) - The proton electric form factor squared (G_E^2)
			(float) - The proton magnetic form factor squared times tau (G_M^2)

	"""

	regression = stats.linregress(epsilon, reduced)
	return (regression[0], regression[1])

def partition(q2, e, theta, cross_section, error, e_prime=None):

	""" Takes lists of cross section data (typically read from datafile) and 
		splits them into sublists based on values of Q^2

	Args:
		q2 				(list) - The list of all Q^2 values
		e 				(list) - The list of all beam energies
		theta 			(list) - The list of all spectrometer angles
		cross_section 	(list) - The list of all differential cross sections
		error 			(list) - The list of all total errors associated with differential cross sections
		e_prime 		(list) - The list of all final beam energies

	Returns:
		q2_partitions 				(list) - The list of lists of homogeneous values of Q^2
		e_partitions 				(list) - The list of lists of beam energies corresponding to each Q^2
		theta_partitions 			(list) - The list of lists of spectrometer angles corresponding to each Q^2
		cross_section_partitions 	(list) - The list of lists of differential cross sections corresponding to each Q^2
		error_partitions 			(list) - The list of lists of total errors of differential cross sections corresponding to each Q^2

	"""

	q2_partitions = []
	e_partitions = []
	theta_partitions = []
	cross_section_partitions = []
	error_partitions = []

	last_same = -1
	for i in range(len(q2)-1):
		if q2[i+1] != q2[i]:
			q2_partitions.append(q2[last_same+1:i+1])
			e_partitions.append(e[last_same+1:i+1])
			theta_partitions.append(theta[last_same+1:i+1])
			cross_section_partitions.append(cross_section[last_same+1:i+1])
			error_partitions.append(error[last_same+1:i+1])
			last_same = i
	q2_partitions.append(q2[last_same+1:])
	e_partitions.append(e[last_same+1:])
	theta_partitions.append(theta[last_same+1:])
	cross_section_partitions.append(cross_section[last_same+1:])
	error_partitions.append(error[last_same+1:])

	return q2_partitions, e_partitions, theta_partitions, cross_section_partitions, error_partitions


def plot_form_factors(ge2_vals, gm2_vals, q2, save_path=None):

	""" Plots the distributions of proton form factors for a given Q^2 side-by-side, based on Monte Carlo simulations

	Args:
		ge2_vals			(list) - The values of the electric form factor squared
		gm2_vals 			(list) - The values of the magnetic form factor squared

	Returns:
		(tuple) - 
			(float) - The uncertainty/standard deviation of the electric form factor measurements
			(float) - The uncertainty/standard deviation of the magnetic form factor measurements

	"""

	plt.delaxes()
	plt.figure(figsize=(24,10))
	ax1 = plt.subplot(121)
	mu1, sigma1 = stats.norm.fit(ge2_vals)

	n, bins, patches = plt.hist(ge2_vals, bins=50, normed=1, facecolor='red', linewidth=3, histtype="stepfilled")
	y = mlab.normpdf(bins, mu1, sigma1)
	plt.plot(bins, y, "-", color="black", linewidth=4)
	plt.axvspan(mu1-sigma1, mu1+sigma1, color="white", alpha=0.5)
	plt.xticks(fontsize=20, rotation=30)
	plt.yticks(fontsize=20)

	ax2 = plt.subplot(122)
	mu2, sigma2 = stats.norm.fit(gm2_vals)
	n, bins, patches = plt.hist(gm2_vals, bins=50, normed=1, facecolor='blue', linewidth=3, histtype="stepfilled")
	y = mlab.normpdf(bins, mu2, sigma2)
	plt.plot(bins, y, "-", color="black", linewidth=4)
	plt.axvspan(mu2-sigma2, mu2+sigma2, color="white", alpha=0.5)
	plt.xticks(fontsize=20, rotation=30)
	plt.yticks(fontsize=20)

	ax1.set_xlabel(r'$G_E^2$', fontsize=30, labelpad=20)
	ax1.set_ylabel(r'Frequency', fontsize=30, labelpad=20)
	ax1.annotate('$Q^2 = %.3f$ \n $\mu = %f$ \n $\sigma = %f$' % (q2, mu1, sigma1), xy=(mu1, 0), xycoords='data',
		xytext=(0.6, 0.7), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=30)
	ax1.text(0.5, 1.05, 'a', transform=ax1.transAxes, 
            size=40, weight='bold')
	ax2.set_xlabel(r'$G_M^2$', fontsize=30, labelpad=20)
	ax2.set_ylabel(r'Frequency', fontsize=30, labelpad=20)
	ax2.annotate('$Q^2 = %.3f$ \n $\mu = %f$ \n $\sigma = %f$' % (q2, mu2, sigma2), xy=(mu2, 0), xycoords='data',
		xytext=(0.65, 0.7), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=30)
	ax2.text(0.5, 1.05, 'b', transform=ax2.transAxes, 
            size=40, weight='bold')
	plt.subplots_adjust(bottom=0.15)
	if not save_path:
		plt.show()
	else:
		plt.savefig(save_path+"/ff_q2_"+str(q2)+".png", bbox_inches="tight")
	return sigma1, sigma2

def dipole_form_factor(q2):

	""" Returns the value of the dipole form factor for a given Q^2 """

	return (1 + q2 / 0.71) ** -2


# won't run if importing above functions
if __name__ == '__main__':

	# to use a data file, specify with a flag
	if len(sys.argv) > 1:
		q_squared = []
		energies = []
		thetas = []
		cross_sections = []
		uncertainties = []

		# keeps track of which cross sections in the data to normalize based on annotations
		next_normalized = False

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
		except IOError as e:
			print(e)
	else:
		# default/test values if no datafile specified
		q_squared = [4.0,4.0,4.0,4.0,4.0]
		energies = [3.4, 3.956, 4.507, 5.507, 9.8]
		thetas = [57.572, 43.707, 35.592, 26.823,13.248]
		cross_sections = [1.297e-2, 2.77e-2, 4.929e-2, 1.023e-1, 6.18e-1]
		uncertainties = [2.243e-4,4.407e-4,7.853e-4,1.370e-3,8.073e-3]

	# reset output file with header row
	with open('out.csv','w') as out:
		writer = csv.writer(out, delimiter=' ')
		writer.writerow(['Q^2', 'G_E^2', 'σ_E', 'G_M^2', 'σ_M', 'G_E/G_D', 'G_M/(μG_D)'])

	# solve for ge^2 and gm^2 with least squares regression
	a = []
	b = []

	# 
	gm2_points = []
	ge2_points = []
	q2_vals = []

	# get lists of lists of values for separate Q^2
	q_squared, energies, thetas, cross_sections, total_errors = partition(q_squared, energies, thetas, cross_sections, uncertainties)
	
	for i in range(len(q_squared)):
		if len(q_squared[i]) <= 1:
			break
		else:
			# for each Q^2, calculate form factors and use Monte Carlo to calculate errors (uncertainties)
			for j in range(len(q_squared[i])):
				q2 = q_squared[i][j]
				energy = energies[i][j]
				theta = thetas[i][j]
				cross_section = cross_sections[i][j]
				error = total_errors[i][j]

				result = rosenbluth(q2, energy, theta, cross_section, error)
				# result:
				# [0] -> epsilon
				# [1] -> tau
				# [2] -> reduced
				# [3] -> reduced error

				# error adjusted with reduced cross section
				total_errors[i][j] = result[3]
				tau = result[1]
				epsilon = result[0]

				a.append(result[0])
				b.append(result[2])

			regression = form_factors(a, b)
			ge_squared = regression[0]
			gm_squared = regression[1]/tau
			print("Q^2 = " + str(q2))
			
			print("G_E^2 = " + str(ge_squared))
			print("G_M^2 = " + str(gm_squared))
			print("G_E / G_D = " + str(ge_squared**0.5/dipole_form_factor(q2)))
			print("G_M / mu G_D = " + str(gm_squared**0.5/pmm/dipole_form_factor(q2)))
				
			# Monte Carlo simulation
			# sample from normal distribution around reduced cross section and epsilon
			# use standard deviation of distribution of computed form factors as error associated with form factor measurements
			eps = []
			red = []
			ge2 = []
			gm2 = []
			samples = 10000
			for j in range(len(a)):
				# put values of epsilon into rows (not randomized)
				s = [a[j]]*samples
				eps.append(s)
				# put randomized values of reduced cross sections into rows
				s = random.normal(b[j], total_errors[i][j], samples)
				red.append(s)
			eps = array(eps)
			red = array(red)


			for j in range(samples):
				# accessing i-th column of matrix (i-th set of epsilons and randomized reduced cross sections)
				eps_samples = eps[:,j]
				red_samples = red[:,j]

				# recalculate form factors
				res = form_factors(eps_samples, red_samples)

				ge2.append(res[0])
				gm2.append(res[1]/tau)

			# clear the lists from Monte Carlo	
			a = []
			b = []
			
			# log and plot results
			sigma1, sigma2 = plot_form_factors(ge2, gm2, q2, "Figures")
			with open('out.csv', 'a') as out:
				writer = csv.writer(out, delimiter=' ')
				writer.writerow([q2, ge_squared, sigma1, gm_squared, sigma2, ge_squared**0.5/dipole_form_factor(q2), gm_squared**0.5/dipole_form_factor(q2)/pmm])
				# writer.writerow([q2, ge_squared, sigma1, gm_squared, sigma2])

