from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys, csv
from os import path
from rosenbluth import pmm, dipole_form_factor
from numpy import arange

ge = []
ge_error = []
gm = []
gm_error = []
q2 = []

def marker():

	""" Cycles through markers to use when plotting the form factors

	Yields:
		s		(string) - The matplotlib symbol for the next marker to use

	"""

	syms = ['s', 'o', '^']
	while True:
		for s in syms:
			yield s

def error(ff2, error2):

	""" Calculates the error associated with the form factor given the form factor squared and its error

	Arguments:
		ff2				(float) - The form factor squared
		error2 			(float) - The error associated with the form factor squared

	Returns:
		(float) The error associated with the form factor 

	"""

	upper_bound_squared = ff2 + error2
	upper_bound = upper_bound_squared ** 0.5
	ff = ff2 ** 0.5
	return upper_bound - ff

if len(sys.argv) > 1:

	# read from 1 or more results files within the Figures subdirectory
	files = sys.argv[1:]

	for i in range(len(files)):

		# store the values from each file in lists of lists
		ge.append([])
		gm.append([])
		ge_error.append([])
		gm_error.append([])
		q2.append([])
		file_path = path.relpath("Figures/"+files[i])
		try:

			with open(file_path, 'r') as f:
				data = csv.reader(f, delimiter=" ")
				# skip header row
				next(data)
				for row in data:
					#do mu ge/gm

					q2[i].append(float(row[0]))
					gd = dipole_form_factor(float(row[0]))
					ge[i].append(float(row[1])/gd**2)
					ge_error[i].append(float(row[2])/gd**2)
					gm[i].append(float(row[3])/gd**2)
					gm_error[i].append(float(row[4])/gd**2)


		except IOError as e:
			print(e)

	plt.style.use("classic")
	rc('font',**{'family':'serif'})

	# plot G_E/G_D
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$G_E^2/G_{SD}^2$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		plt.errorbar([j+i/20 for j in q2[i]], ge[i], yerr=ge_error[i], fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)
	
	# plot a fit line where G_E = (1+x/0.662)^{-2} and G_D = (1+x/0.71)^{-2}
	x = arange(0., 8., 0.1)
	y = (1+x/0.662)**-4 / (1+x/0.71)**-4
	plt.plot(x, y)

	# shade unphysical region
	if ax.get_ylim()[0] < 0:
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)
	
	file_names = [i.split("+")[1][:-4] for i in files]
	plt.savefig("Figures/ge-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")
	# plt.show()

	# plot G_M/mu G_D
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$G_M^2/\mu^2 G_{SD}^2$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		plt.errorbar([j+i/20 for j in q2[i]], [gm[i][j]/pmm**2 for j in range(len(gm[i]))], yerr=[gm_error[i][j]/pmm**2 for j in range(len(gm_error[i]))], fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)

	# plot a fit line where G_M = mu
	x = arange(0., 8., 0.1)
	y = [1.0]*len(x)
	plt.plot(x, y)

	# shade unphysical region
	if ax.get_ylim()[0] < 0:
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)

	file_names = [i.split("+")[1][:-4] for i in files]
	plt.savefig("Figures/gm-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")
	# plt.show()

	# plot mu G_E/G_M
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$\mu^2 G_E^2/G_M^2$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		# recalculate error for a product/ratio of variables (see http://www.utm.edu/staff/cerkal/Lect4.html)
		ratio_err = [pmm**2*ge[i][j]/gm[i][j] * ((ge_error[i][j]/ge[i][j])**2 + (gm_error[i][j]/gm[i][j])**2)**0.5 for j in range(len(ge[i]))]

		plt.errorbar([j+i/20 for j in q2[i]], [pmm**2*ge[i][j]/gm[i][j] for j in range(len(ge[i]))], yerr=ratio_err, fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)

	# plot a fit line where mu G_E/G_M = (1-x/8.02)
	x = arange(0., 8., 0.1)
	y = (1-x/8.02)**2
	plt.plot(x, y)

	# shade unphysical region
	if ax.get_ylim()[0] < 0:
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)

	file_names = [i.split("+")[1][:-4] for i in files]
	plt.savefig("Figures/ge-gm-" + "&".join(file_names) + ".png", bbox_inches="tight")
	# plt.show()
	