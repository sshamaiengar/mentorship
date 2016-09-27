from numpy import deg2rad, linalg, tan, sin, cos, random, array
import matplotlib.pyplot as plt
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

def plot_form_factors(ge2_vals, gm2_vals, q2):
	#label the Q^2
	f1 = plt.figure()
	f2 = plt.figure()
	ax1 = f1.add_subplot(111)
	ax1.hist(ge2, bins=50, facecolor='red')
	ax2 = f2.add_subplot(111)
	ax2.hist(gm2, bins=50, facecolor='blue')

	# plt.hist(ge2, bins=50, facecolor='red')
	ax1.set_title(r'$\mathrm{Histogram\ of\ G_E^2: Q^2 = '+str(q2)+'}$', {'fontsize':20})
	ax1.set_xlabel(r'$G_E^2$', {'fontsize':20})

	ax2.set_title(r'$\mathrm{Histogram\ of\ G_M^2: Q^2 = '+str(q2)+'}$', {'fontsize':20})
	ax2.set_xlabel(r'$G_M^2$', {'fontsize':20})

	plt.show()

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

percent_errors = []

for i in range(len(energies)):

	q2 = q_squared[i]
	energy = energies[i]
	theta = thetas[i]
	cross_section = cross_sections[i]
	percent_errors.append(uncertainties[i])

	result = rosenbluth(q2, energy, theta, cross_section)
	tau = result[1]
	# a.append([result[0], result[1]])
	# b.append(result[2])

	a.append(result[0])
	b.append(result[2])

# regression = linalg.lstsq(a, b)[0]
# ge_squared = regression[0]
# gm_squared = regression[1]
	if i < len(energies)-1:
		if q_squared[i+1] != q2:
			regression = form_factors(a, b)
			ge_squared = regression[0]
			gm_squared = regression[1]/tau
			print("Q^2 = " + str(q2))
			print("G_e^2 = " + str(ge_squared))
			print("G_m^2 = " + str(gm_squared))


			# Monte Carlo simulation

			# sample from normal distribution around reduced cross section and epsilon
			eps = []
			red = []

			ge2 = []
			gm2 = [] # actually tau*gm2

			samples = 10000

			for j in range(len(a)):

				# use uncertainties instead of 2%

				#epsilon
				# s = random.normal(a[j][0], 0.000000002, samples)
				s = [a[j]]*samples
				eps.append(s)
				#reduced
				s = random.normal(b[j], percent_errors[j], samples)
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
					# ge2.append(0)
					pass
				else:

					#divide out tau
					ge2.append(res[0])
				gm2.append(res[1]/tau)

			a = []
			b = []
			percent_errors = []

			plot_form_factors(ge2, gm2, q2)
	elif i == len(energies)-1:
		regression = form_factors(a, b)
		ge_squared = regression[0]
		gm_squared = regression[1]/tau

		print("Q^2 = " + str(q2))
		print("G_e^2 = " + str(ge_squared))
		print("G_m^2 = " + str(gm_squared))


		# Monte Carlo simulation

		# sample from normal distribution around reduced cross section and epsilon
		eps = []
		red = []

		ge2 = []
		gm2 = [] # actually tau*gm2

		samples = 10000

		for j in range(len(a)):

			# use uncertainties instead of 2%

			#epsilon
			# s = random.normal(a[j][0], 0.000000002, samples)
			s = [a[j]]*samples
			eps.append(s)
			#reduced
			s = random.normal(b[j], percent_errors[j], samples)
			red.append(s)

		eps = array(eps)
		red = array(red)


		for j in range(samples):
			# accessing j-th column of matrix
			eps_samples = eps[:,j]
			red_samples = red[:,j]
			res = form_factors(eps_samples, red_samples)

			#negative slope of fit becomes 0
			if res[0] < 0:
				# ge2.append(0)
				pass
			else:

				#divide out tau
				ge2.append(res[0])
			gm2.append(res[1]/tau)

		a = []
		b = []
		percent_errors = []

		plot_form_factors(ge2, gm2, q2)


