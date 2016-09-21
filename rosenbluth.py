from numpy import deg2rad, linalg, tan, sin, cos, random, array
import matplotlib.pyplot as plt
import scipy.stats as stats

# constants
fsc = 0.00729735256
mass = 0.938 # GeV

# q_squared = (2 * mass * energy ** 2 * (1 - cos(deg2rad(theta)))) / (mass + energy * (1 - cos(deg2rad(theta))))
q_squared = 4.0

# solve for ge^2 and gm^2 with least squares regression
a = []
b = []

energies = [3.4, 3.956, 4.507, 5.507, 9.8]
thetas = [57.572, 43.707, 35.592, 26.823,13.248]
cross_sections = [1.297e-2, 2.77e-2, 4.929e-2, 1.023e-1, 6.18e-1]
# energies = [3.4, 9.8]
# thetas = [57.572, 13.248]
# cross_sections = [1.297e-2, 6.18e-1]

# returns epsilon, tau, and reduced cross section (can be used to calculate form factors with linear fit)
def rosenbluth(energy, theta, cross_section):
	
	tau = q_squared / 4 / mass ** 2
	epsilon = (1 + 2 * (1 + tau) * tan(deg2rad(theta)/2) ** 2) ** -1
	eta = 1 + (energy / mass) * (1 - cos(deg2rad(theta)))

	# 1 GeV^-2 = 0.389 mb
	ideal_scattering = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * 0.389 * 1e6 * eta ** -1
	reduced = cross_section/ideal_scattering * epsilon * (1 + tau)
	return (epsilon, tau, reduced)

# calculates the form factors based on epsilon and reduced cross section
def form_factors(epsilon, reduced):
	regression = stats.linregress(epsilon, reduced)
	# returns slope, intercept
	return (regression[0], regression[1])


for i in range(len(energies)):
	# energy = float(input("Beam energy (GeV): "));
	# theta = float(input("Scattering angle (deg): "))
	# cross_section = float(input("Cross section (nb/sr): "))

	energy = energies[i]
	theta = thetas[i]
	cross_section = cross_sections[i]

	result = rosenbluth(energy, theta, cross_section)

	a.append([result[0], result[1]])
	b.append(result[2])
	regression = linalg.lstsq(a, b)[0]
	ge_squared = regression[0]
	gm_squared = regression[1]



print("G_e^2 = " + str(ge_squared))
print("G_m^2 = " + str(gm_squared))

# Monte Carlo simulation

# sample from normal distribution around reduced cross section and epsilon
eps = []
red = []

ge2 = []
gm2 = [] # actually tau*gm2

samples = 10000

for i in range(len(a)):

	#epsilon
	s = random.normal(a[i][0], 0.02*a[i][0], samples)
	eps.append(s)
	#reduced
	s = random.normal(b[i], 0.02*b[i], samples)
	red.append(s)

eps = array(eps)
red = array(red)


# number of points in fit
for i in range(samples):
	# accessing i-th column of matrix
	eps_samples = eps[:,i]
	red_samples = red[:,i]
	res = form_factors(eps_samples, red_samples)

	#negative slope of fit becomes 0
	if res[0] < 0:
		ge2.append(0)
		pass
	else:

		#divide out tau
		ge2.append(res[0])
	gm2.append(res[1])


f1 = plt.figure()
f2 = plt.figure()
ax1 = f1.add_subplot(111)
ax1.hist(ge2, bins=50, facecolor='red')
ax2 = f2.add_subplot(111)
ax2.hist(gm2, bins=50, facecolor='blue')

# plt.hist(ge2, bins=50, facecolor='red')
ax1.set_title(r'$\mathrm{Histogram\ of\ G_E^2:}\ \sigma=2\%$', {'fontsize':20})
ax1.set_xlabel(r'$G_E^2$', {'fontsize':20})

ax2.set_title(r'$\mathrm{Histogram\ of\ \tau G_M^2:}\ \sigma=2\%$', {'fontsize':20})
ax2.set_xlabel(r'$\tau G_M^2$', {'fontsize':20})

plt.show()


