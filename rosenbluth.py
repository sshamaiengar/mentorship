from numpy import *

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

for i in range(len(energies)):
	# energy = float(input("Beam energy (GeV): "));
	# theta = float(input("Scattering angle (deg): "))
	# cross_section = float(input("Cross section (nb/sr): "))

	energy = energies[i]
	theta = thetas[i]
	cross_section = cross_sections[i]

	tau = q_squared / 4 / mass ** 2
	# print("tau: " + str(tau))
	epsilon = (1 + 2 * (1 + tau) * tan(deg2rad(theta)/2) ** 2) ** -1
	# print("theta rads: " + str(deg2rad(theta)))
	# print("epsilon: " + str(epsilon))
	eta = 1 + (energy / mass) * (1 - cos(deg2rad(theta)))
	# print("eta: " + str(eta))

	# 1 GeV^-2 = 0.389 mb
	ideal_scattering = (1 ** 2 * fsc ** 2 * cos(deg2rad(theta / 2)) ** 2) / (4 * energy ** 2 * sin(deg2rad(theta / 2)) ** 4) * 0.389 * 1e6 * eta ** -1
	# print("Mott: " + str(ideal_scattering))
	a.append([epsilon, tau])
	b.append(cross_section/ideal_scattering * epsilon * (1 + tau))
	# print("reduced: " + str(cross_section/ideal_scattering * epsilon * (1 + tau)))

regression = linalg.lstsq(a, b)[0]

ge_squared = regression[0]
gm_squared = regression[1]

print("G_e^2 = " + str(ge_squared))
print("G_m^2 = " + str(gm_squared))