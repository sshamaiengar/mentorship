# -*- coding: utf-8 -*-


# problem could be 8 gev spec at 90 degs
# try leaving those two runs out
# compare answers
# ALSO
# use asymmetry data to get formula for G_E, then substitute into reduced cross section and do one-parameter fit


from numpy import deg2rad, linalg, tan, sin, cos, random, array, std, mean, where
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys
import csv


# constants
fsc = 0.00729735256
# mass = 0.93827 # GeV
mass = 0.938 # GeV
# pmm = 2.793 # proton magnetic moment
pmm = 2.79 # proton magnetic moment


# returns epsilon, tau, and reduced cross section (can be used to calculate form factors with linear fit)
def rosenbluth(q_squared, energy, theta, cross_section, error, energy_prime=None):
	
	tau = q_squared / 4 / mass ** 2
	epsilon = (1 + 2 * (1 + tau) * tan(deg2rad(theta)/2) ** 2) ** -1
	eta = 1 + (energy / mass) * (1 - cos(deg2rad(theta)))

	# 1 GeV^-2 = 0.389 mb
	if energy_prime:
		ideal_scattering = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * energy_prime/energy * 0.389 * 1e6
		reduced = cross_section/ideal_scattering * epsilon * (1 + tau)
		# reduced = cross_section * (1+tau)/5.18 * epsilon/tau * energy**3/energy_prime * sin(deg2rad(theta / 2)) ** 4/cos(deg2rad(theta / 2)) ** 2 * dipole_form_factor(q2) ** -2
		error = error/ideal_scattering * epsilon * (1 + tau)
	else:
		ideal_scattering = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * 0.389 * 1e6 * eta ** -1
		reduced = cross_section/ideal_scattering * epsilon * (1 + tau)
		error = error/ideal_scattering * epsilon * (1 + tau)
	return (epsilon, tau, reduced, error)

# calculates the form factors based on epsilon and reduced cross section
# add true weights from table in notebook
def form_factors(epsilon, reduced):
	regression = stats.linregress(epsilon, reduced)
	# returns slope, intercept
	return (regression[0], regression[1])

# split data based on Q^2 to make separate calculations and plots
def partition(q2, e, theta, cross_section, error, e_prime=None):
	q2_partitions = []
	e_partitions = []
	theta_partitions = []
	cross_section_partitions = []
	error_partitions = []

	# if using data4.csv or dataset with E'
	if e_prime:
		e_prime_partitions = []
	last_same = -1
	for i in range(len(q2)-1):
		if q2[i+1] != q2[i]:
			q2_partitions.append(q2[last_same+1:i+1])
			e_partitions.append(e[last_same+1:i+1])
			theta_partitions.append(theta[last_same+1:i+1])
			cross_section_partitions.append(cross_section[last_same+1:i+1])
			error_partitions.append(error[last_same+1:i+1])
			if e_prime:
				e_prime_partitions.append(e_prime[last_same+1:i+1])
			last_same = i
	q2_partitions.append(q2[last_same+1:])
	e_partitions.append(e[last_same+1:])
	theta_partitions.append(theta[last_same+1:])
	cross_section_partitions.append(cross_section[last_same+1:])
	error_partitions.append(error[last_same+1:])

	# if using data4.csv or dataset with E'
	if e_prime:
		e_prime_partitions.append(e_prime[last_same+1:])
		return q2_partitions, e_partitions, theta_partitions, cross_section_partitions, error_partitions, e_prime_partitions
	return q2_partitions, e_partitions, theta_partitions, cross_section_partitions, error_partitions


# returns uncertainties (sigma on gaussian fit of histogram) for ge2 and gm2 respectively
def plot_form_factors(ge2_vals, gm2_vals, q2):

	#increase font size on y axis label and numbers
	# get rid of bars in between
	# white background, remove title
	# information about mu, sigma, in plot, not title

	# f1, axes = plt.subplots(1, 2, figsize=(20,10))
	plt.delaxes()
	plt.figure(figsize=(24,10))
	# ax1 = f1.add_subplot(121)
	ax1 = plt.subplot(121)
	# ax1.hist(ge2, bins=50, facecolor='red')
	mu1, sigma1 = stats.norm.fit(ge2_vals)

	n, bins, patches = plt.hist(ge2_vals, bins=50, normed=1, facecolor='red', linewidth=3, histtype="stepfilled")
	y = mlab.normpdf(bins, mu1, sigma1)
	plt.plot(bins, y, "-", color="black", linewidth=4)
	plt.axvspan(mu1-sigma1, mu1+sigma1, color="white", alpha=0.5)
	plt.xticks(fontsize=20, rotation=30)
	plt.yticks(fontsize=20)

	# ax2 = f1.add_subplot(122)
	ax2 = plt.subplot(122)
	# ax2.hist(gm2, bins=50, facecolor='blue')
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
	# ax1.annotate('$2\sigma$', xy=(mu1, ax1.get_ylim()[1]/3), xycoords='data', xytext=(mu1, ax1.get_ylim()[1]/3), textcoords='data', horizontalalignment='center', fontsize=30)
	ax2.set_xlabel(r'$G_M^2$', fontsize=30, labelpad=20)
	ax2.set_ylabel(r'Frequency', fontsize=30, labelpad=20)
	ax2.annotate('$Q^2 = %.3f$ \n $\mu = %f$ \n $\sigma = %f$' % (q2, mu2, sigma2), xy=(mu2, 0), xycoords='data',
		xytext=(0.65, 0.7), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=30)
	ax2.text(0.5, 1.05, 'b', transform=ax2.transAxes, 
            size=40, weight='bold')
	# ax2.annotate('$2\sigma$', xy=(mu2, ax2.get_ylim()[1]/3), xycoords='data', xytext=(mu2, ax2.get_ylim()[1]/3), textcoords='data', horizontalalignment='center', fontsize=30)
	plt.subplots_adjust(bottom=0.15)
	# plt.show()
	return sigma1, sigma2

def dipole_form_factor(q2):
	return (1+q2/0.71)**-2


# only run this stuff if specifically executing this file
# won't run if importing above functions
if __name__ == '__main__':


	# sys.argv is a list of the command line flags
	# to use a data file, specify with a flag
	if len(sys.argv) > 1:
		q_squared = []
		energies = []
		thetas = []
		cross_sections = []
		uncertainties = []
		energies_prime = []

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
						# only normalize 1.6 GeV data
						if len(row) > 1 and "1.6" in row[1]:
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
					# print("XS = " + str(float(row[3])*norm))
					uncertainties.append(float(row[4])*norm)
					# print("Err = " + str(float(row[4])*norm))
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

	# reset output file
	with open('out.csv','w') as out:
		writer = csv.writer(out, delimiter=' ')
		writer.writerow(['Q^2', 'G_E^2', 'σ_E', 'G_M^2', 'σ_M', 'G_E/G_D', 'G_M/(μG_D)'])

	# print(q_squared)
	# print(energies)
	# print(thetas)
	# print(cross_sections)
	# print(uncertainties)

	# solve for ge^2 and gm^2 with least squares regression
	a = []
	b = []


	gm2_points = []
	ge2_points = []
	q2_vals = []

	# if not using data4.csv or dataset with E'
	if not energies_prime:
		q_squared, energies, thetas, cross_sections, total_errors = partition(q_squared, energies, thetas, cross_sections, uncertainties)
	else:
		q_squared, energies, thetas, cross_sections, total_errors, energies_prime = partition(q_squared, energies, thetas, cross_sections, uncertainties, energies_prime)
	
	for i in range(len(q_squared)):
		if len(q_squared[i]) <= 1:
			break
		else:
			for j in range(len(q_squared[i])):
				q2 = q_squared[i][j]
				energy = energies[i][j]
				theta = thetas[i][j]
				cross_section = cross_sections[i][j]
				error = total_errors[i][j]

				# if not using data4.csv or dataset with E'
				if not energies_prime:
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
				else:
					energy_prime = energies_prime[i][j]
					result = rosenbluth(q2, energy, theta, cross_section, error, energy_prime)
					total_errors[i][j] = result[3]
					tau = result[1]
					epsilon = result[0]
					a.append(epsilon)
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
			eps = []
			red = []
			ge2 = []
			gm2 = []
			samples = 10000
			for j in range(len(a)):
				#epsilon
				s = [a[j]]*samples
				eps.append(s)
				#reduced
				s = random.normal(b[j], total_errors[i][j], samples)
				red.append(s)
			eps = array(eps)
			red = array(red)


			for j in range(samples):
				# accessing i-th column of matrix
				eps_samples = eps[:,j]
				red_samples = red[:,j]
				res = form_factors(eps_samples, red_samples)
				#negative slope of fit becomes 0

				if res[0] < 0:
					ge2.append(res[0])
					pass
				else:

					#divide out tau
					ge2.append(res[0])
				gm2.append(res[1]/tau)

			a = []
			b = []
			# ge2_points.append(ge2)
			# gm2_points.append(gm2)
			# q2_vals.append(q2)
			sigma1, sigma2 = plot_form_factors(ge2, gm2, q2)
			with open('out.csv', 'a') as out:
				writer = csv.writer(out, delimiter=' ')
				writer.writerow([q2, ge_squared, sigma1, gm_squared, sigma2, ge_squared**0.5/dipole_form_factor(q2), gm_squared**0.5/dipole_form_factor(q2)/pmm])
				# writer.writerow([q2, ge_squared, sigma1, gm_squared, sigma2])

