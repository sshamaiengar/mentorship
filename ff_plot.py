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
	syms = ['s', 'o', '^']
	while True:
		for s in syms:
			yield s

def error(ff2, error2):
	upper_bound_squared = ff2 + error2
	upper_bound = upper_bound_squared ** 0.5
	ff = ff2 ** 0.5
	return upper_bound - ff

if len(sys.argv) > 1:
	files = sys.argv[1:]
	for i in range(len(files)):
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
					ge[i].append(float(row[1])**0.5/gd)
					ge_error[i].append(error(float(row[1])/gd**2, float(row[2]))/gd**2)
					gm[i].append(float(row[3])**0.5/gd)
					gm_error[i].append(error(float(row[3])/gd**2, float(row[4]))/gd**2)

					# FIX ERRORS:
					# need to sqrt(ff^2 + ff_err^2) and sqrt(ff^2 - ff_err^2)
					# then get sqrt(ff^2)
					# take differences to get new error

					# divide gm by mu for gm/gsd
		except IOError as e:
			print(e)

	plt.style.use("classic")
	rc('font',**{'family':'serif'})

	# ge/gd
	# NOT mu ge/gd, ge should go to 1
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$G_E/G_{SD}$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		# print(q2[i], ge[i], ge_error[i])
		# plt.plot(q2[i], ge[i], marker = next(mark))
		plt.errorbar([j+i/20 for j in q2[i]], ge[i], yerr=ge_error[i], fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)
	x = arange(0., 8., 0.1)
	y = (1+x/0.662)**-2 / (1+x/0.71)**-2
	plt.plot(x, y)
	
	# plt.show()
	file_names = [i.split("+")[1][:-4] for i in files]
	# print(file_names)
	plt.savefig("Figures/ge-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")

	# gm/gd
	# not gm/gd, gm should go to 2.7 (mu)
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$G_M/\mu G_{SD}$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		# print(q2[i], ge[i], ge_error[i])
		# plt.plot(q2[i], ge[i], marker = next(mark))
		plt.errorbar([j+i/20 for j in q2[i]], [gm[i][j]/pmm for j in range(len(gm[i]))], yerr=[gm_error[i][j]/pmm for j in range(len(gm_error[i]))], fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)
	x = arange(0., 8., 0.1)
	y = [1.0]*len(x)
	plt.plot(x, y)
	# plt.show()
	file_names = [i.split("+")[1][:-4] for i in files]
	# print(file_names)
	plt.savefig("Figures/gm-gd-" + "&".join(file_names) + ".png", bbox_inches="tight")

	# mu ge/gm
	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)
	ax.set_xlabel(r'$Q^2$', fontsize=30)
	ax.set_ylabel(r'$\mu G_E/G_M$', fontsize=30)
	ax.set_xlim([0,8.0])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	mark = marker()
	for i in range(len(q2)):
		# print(q2[i], ge[i], ge_error[i])
		# plt.plot(q2[i], ge[i], marker = next(mark))

		# recalculate error (multiplication is special case, see http://www.utm.edu/staff/cerkal/Lect4.html)
		ratio_err = [pmm*ge[i][j]/gm[i][j] * ((ge_error[i][j]/ge[i][j])**2 + (gm_error[i][j]/gm[i][j])**2)**0.5 for j in range(len(ge[i]))]

		plt.errorbar([j+i/20 for j in q2[i]], [pmm*ge[i][j]/gm[i][j] for j in range(len(ge[i]))], yerr=ratio_err, fmt=next(mark), elinewidth=2, markersize=10, markerfacecolor="None", markeredgewidth=2)
	x = arange(0., 8., 0.1)
	y = (1-x/8.02)
	plt.plot(x, y)
	# plt.show()
	file_names = [i.split("+")[1][:-4] for i in files]
	# print(file_names)
	plt.savefig("Figures/ge-gm-" + "&".join(file_names) + ".png", bbox_inches="tight")
	
	