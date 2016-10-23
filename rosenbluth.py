# -*- coding: utf-8 -*-

from numpy import deg2rad, linalg, tan, sin, cos, random, array, std, mean
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys
import csv

# returns epsilon, tau, and reduced cross section (can be used to calculate form factors with linear fit)
def rosenbluth(q_squared, energy, theta, cross_section):
	
	tau = q_squared / 4 / mass ** 2
	epsilon = (1 + 2 * (1 + tau) * tan(deg2rad(theta)/2) ** 2) ** -1
	eta = 1 + (energy / mass) * (1 - cos(deg2rad(theta)))

	# 1 GeV^-2 = 0.389 mb
	ideal_scattering = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * 0.389 * 1e6 * eta ** -1
	reduced = cross_section/ideal_scattering * epsilon * (1 + tau)
	return (epsilon, tau, reduced)

# calculates the form factors based on epsilon and reduced cross section
# add true weights from table in notebook
def form_factors(epsilon, reduced):
	regression = stats.linregress(epsilon, reduced)
	# returns slope, intercept
	return (regression[0], regression[1])

# split data based on Q^2 to make separate calculations and plots
def partition(q2, e, theta, cross_section, error):
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


# returns uncertainties (sigma on gaussian fit of histogram) for ge2 and gm2 respectively
def plot_form_factors(ge2_vals, gm2_vals, q2):

	#increase font size on y axis label and numbers
	# get rid of bars in between
	# white background, remove title
	# information about mu, sigma, in plot, not title

	# f1, axes = plt.subplots(1, 2, figsize=(20,10))
	plt.delaxes()
	plt.figure(figsize=(24,9))
	# ax1 = f1.add_subplot(121)
	ax1 = plt.subplot(121)
	# ax1.hist(ge2, bins=50, facecolor='red')
	mu1, sigma1 = stats.norm.fit(ge2_vals)

	n, bins, patches = plt.hist(ge2_vals, bins=50, normed=1, facecolor='red', linewidth=3, histtype="stepfilled")
	y = mlab.normpdf(bins, mu1, sigma1)
	plt.plot(bins, y, "-", color="orange", linewidth=4)
	plt.axvspan(mu1, mu1+sigma1, color="white", alpha=0.5)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	# ax2 = f1.add_subplot(122)
	ax2 = plt.subplot(122)
	# ax2.hist(gm2, bins=50, facecolor='blue')
	mu2, sigma2 = stats.norm.fit(gm2_vals)
	n, bins, patches = plt.hist(gm2_vals, bins=50, normed=1, facecolor='blue', linewidth=3, histtype="stepfilled")
	y = mlab.normpdf(bins, mu2, sigma2)
	plt.plot(bins, y, "-", color="orange", linewidth=4)
	plt.axvspan(mu2, mu2+sigma1, color="white", alpha=0.5)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	# ax1.set_title(r'$\mathrm{Histogram\ of\ G_E^2: Q^2 = %.3f},\ \mu = %f,\ \sigma = %f$' % (q2,mu1,sigma1), fontsize=20)
	ax1.set_xlabel(r'$G_E^2$', fontsize=20)
	ax1.set_ylabel(r'Frequency', fontsize=20)
	ax1.annotate('$Q^2 = %f$ \n $\mu = %f$ \n $\sigma = %f$' % (q2, mu1, sigma1), xy=(mu1, 0), xycoords='data',
		xytext=(0.7, 0.8), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=20)
	# ax2.set_title(r'$\mathrm{Histogram\ of\ G_M^2: Q^2 = %.3f},\ \mu = %f,\ \sigma = %f$' % (q2,mu2,sigma2), fontsize=20)
	ax2.set_xlabel(r'$G_E^2$', fontsize=20)
	ax2.set_ylabel(r'Frequency', fontsize=20)
	ax2.annotate('$Q^2 = %f$ \n $\mu = %f$ \n $\sigma = %f$' % (q2, mu2, sigma2), xy=(mu2, 0), xycoords='data',
		xytext=(0.7, 0.8), textcoords="axes fraction", verticalalignment='bottom', horizontalalignment='left', fontsize=20)
	plt.show()
	return sigma1, sigma2

# sys.argv is a list of the command line flags
# to use a data file, specify with a flag

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

# reset output file
with open('out.csv','w') as out:
	writer = csv.writer(out, delimiter=' ')
	writer.writerow(['Q^2', 'G_E^2', 'σ_E', 'G_M^2', 'σ_M'])

# print(q_squared)
# print(energies)
# print(thetas)
# print(cross_sections)
# print(uncertainties)

# constants
fsc = 0.00729735256
mass = 0.938 # GeV

# solve for ge^2 and gm^2 with least squares regression
a = []
b = []


#keep track of previous q^2 to make separate estimates of form factors

gm2_points = []
ge2_points = []
q2_vals = []

q_squared, energies, thetas, cross_sections, total_errors = partition(q_squared, energies, thetas, cross_sections, uncertainties)
for i in range(len(q_squared)):
	if len(q_squared[i]) <= 1:
		break
	else:
		for j in range(len(q_squared[i])):
			q2 = q_squared[i][j]
			energy = energies[i][j]
			theta = thetas[i][j]
			cross_section = cross_sections[i][j]
			# total_errors.append(uncertainties[i][j])

			result = rosenbluth(q2, energy, theta, cross_section)
			tau = result[1]

			a.append(result[0])
			b.append(result[2])

		regression = form_factors(a, b)
		ge_squared = regression[0]
		gm_squared = regression[1]/tau
		print("Q^2 = " + str(q2))
		print("G_e^2 = " + str(ge_squared))
		print("G_m^2 = " + str(gm_squared))
		print
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
			s = random.normal(b[j], total_errors[i][j]/100, samples)
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
			writer.writerow([q2, ge_squared, sigma1, gm_squared, sigma2])

