from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys, csv
from os import path
from rosenbluth import pmm, dipole_form_factor
from numpy import arange
from monte_carlo_ratio_errors import calc_ff_ratio_err
ge = []
ge_error = []
gm = []
gm_error = []
q2 = []
ge_gm_error = []

def marker():

	""" Cycles through markers to use when plotting the form factors

	Yields:
		s		(string) - The matplotlib symbol for the next marker to use

	"""

	syms = ['o', '^', 's']
	while True:
		for s in syms:
			yield s

def color():

	""" Cycles through colors to use when plotting the form factors

	Yields:
		c 		(string) - The matplotlib identifier for the next color to use

	"""
	colors = ['b', 'r', 'g', 'k']
	while True:
		for c in colors:
			yield c

def error(ff2, error2):

	""" Calculates the error associated with the form factor given the form factor squared and its error

	Args:
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
		ge_gm_error.append([])
		if "/" not in files[i]:
			file_path = path.relpath("Figures/"+files[i])
		else:
			file_path = path.relpath(files[i])
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
					ge_gm_error[i].append(float(row[7]))

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
	c = color()
	fill = True
	c1 = next(c)
	m1 = next(mark)
	for i in range(len(q2)):

		if fill:
			plt.errorbar([j+i/20 for j in q2[i]], ge[i], yerr=ge_error[i], fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor=c1, markeredgecolor=c1, markeredgewidth=2)
		else:
			plt.errorbar([j+i/20 for j in q2[i]], ge[i], yerr=ge_error[i], fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor='none', markeredgecolor=c1, markeredgewidth=2)
			c1 = next(c)
			m1 = next(mark)
			
		fill = not fill
		# c1 = next(c)

	# plot a fit line where G_E = (1+x/0.662)^{-2} and G_D = (1+x/0.71)^{-2}
	x = arange(0., 8., 0.1)
	y = (1+x/0.662)**-4 / (1+x/0.71)**-4
	plt.plot(x, y, '-g', label='Carlson-Griffioen')

	# add a legend
	# if len(files) > 1:
	# 	handles, labels = ax.get_legend_handles_labels()
	# 	# handles.append(fit_curve)
	# 	print(handles, labels)
	# 	# handles = [h[0] for h in handles[0]]
	# 	ax.legend(handles, labels, loc='upper left', numpoints=1, fontsize=20)

	# shade unphysical region
	lim = ax.get_ylim()
	if ax.get_ylim()[0] < 0:
		
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)
	ax.set_ylim(lim)
	# file_names = [i.split("+")[1][:-4] for i in files]
	# plt.savefig("Figures/Form_Factors/ge-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")
	plt.show()

	# plot G_M/mu G_D
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$G_M^2/\mu^2 G_{SD}^2$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	c = color()
	fill = True
	c1 = next(c)
	m1 = next(mark)
	for i in range(len(q2)):
		if fill:
			plt.errorbar([j+i/20 for j in q2[i]], [gm[i][j]/pmm**2 for j in range(len(gm[i]))], yerr=[gm_error[i][j]/pmm**2 for j in range(len(gm_error[i]))], fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor=c1, markeredgecolor=c1, markeredgewidth=2)
		else:
			plt.errorbar([j+i/20 for j in q2[i]], [gm[i][j]/pmm**2 for j in range(len(gm[i]))], yerr=[gm_error[i][j]/pmm**2 for j in range(len(gm_error[i]))], fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor='none', markeredgecolor=c1, markeredgewidth=2)
			c1 = next(c)
			m1 = next(mark)
			
		fill = not fill
		# c1 = next(c)

	# plot a fit line where G_M = mu
	x = arange(0., 8., 0.1)
	y = [1.0]*len(x)
	plt.plot(x, y, '-g')

	# shade unphysical region
	if ax.get_ylim()[0] < 0:
		
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)

	# file_names = [i.split("+",maxsplit=1)[1][:-4] for i in files]
	# plt.savefig("Figures/Form_Factors/gm-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")
	plt.show()

	# plot mu G_E/G_M
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$\mu^2 G_E^2/G_M^2$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	c = color()
	c1 = next(c)
	fill = True
	m1 = next(mark)
	for i in range(len(q2)):
		# recalculate error for a product/ratio of variables (see http://www.utm.edu/staff/cerkal/Lect4.html)
		# ratio_err = [pmm**2 * ge[i][j] / gm[i][j] * ( (ge_error[i][j] / ge[i][j]) ** 2 + (gm_error[i][j] / gm[i][j]) ** 2) ** 0.5 for j in range(len(ge[i]))]
		ratio_err = ge_gm_error[i]
		if fill:
			plt.errorbar([j+i/20 for j in q2[i]], [pmm**2*ge[i][j]/gm[i][j] for j in range(len(ge[i]))], yerr=ratio_err, fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor=c1, markeredgecolor=c1, markeredgewidth=2)
		else:
			plt.errorbar([j+i/20 for j in q2[i]], [pmm**2*ge[i][j]/gm[i][j] for j in range(len(ge[i]))], yerr=ratio_err, fmt=m1, elinewidth=2, markersize=10, ecolor=c1, markerfacecolor='none', markeredgecolor=c1, markeredgewidth=2)
			c1 = next(c)
			m1 = next(mark)
			
		fill = not fill
		# c1 = next(c)

	# plot a fit line where mu G_E/G_M = (1-x/8.02)
	x = arange(0., 8., 0.1)
	y = (1-x/8.02)**2
	plt.plot(x, y, '-g')

	# shade unphysical region
	if ax.get_ylim()[0] < 0:
		plt.axhspan(ax.get_ylim()[0], 0, facecolor='gray', alpha=0.25)

	# file_names = [i.split("+")[1][:-4] for i in files]
	# plt.savefig("Figures/Form_Factors/ge-gm-" + "&".join(file_names) + ".png", bbox_inches="tight")
	plt.show()
	