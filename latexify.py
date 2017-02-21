import sys, csv
from os import path

if len(sys.argv) == 2:
	file = sys.argv[1]
	file_path = path.relpath("Figures/"+file)

	print(file_path)
	try:

		with open(file_path, 'r') as f:
			data = csv.reader(f, delimiter=" ")
			# skip header row
			next(data)
			for row in data:
				print("{:.2f}".format(float(row[0])), end='')
				for i in [row[1], row[2]]:
					string = "{:.2e}".format(float(i)).split("e")
					print(" & \\({}\\times10^{{{:d}}}\\)".format(string[0], int(string[1])), end='')
				print(" & \\({:.3f}\)".format(float(row[5])), end='')
				for i in [row[3], row[4]]:
					string = "{:.2e}".format(float(i)).split("e")
					print(" & \\({}\\times10^{{{:d}}}\\)".format(string[0], int(string[1])), end='')
				print(" & \\({:.3f}\)".format(float(row[6])), end='')
				print("\\\\")
	except IOError as e:
		print(e)