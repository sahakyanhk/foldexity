import os

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
os.chdir(script_directory)


def parse_scop_fasta(scop_fasta, output_file):
	with open(output_file, "w") as f_out:
		#f_out.write("scopid\tfoldid\ttaxid\tseqlen\n")
		writedata = False
		seq = ""

		for oline in open(scop_fasta):
			if oline.startswith(">"):
				if writedata:
					f_out.write(f"{scopid}\t{foldid}\t{taxid}\t{len(seq)}\n")
					seq = ""

				sline1 = oline.split(' ')
				scopid = sline1[0].replace(">", "")
				foldid = sline1[1]
				if "[TaxId: " in oline:
					taxid = oline.split("[TaxId: ")[1].split("]")[0]
				else:
					taxid = "missing"
				writedata = True
			else:
				seq += oline.strip()

		# Write last entry
		f_out.write(f"{scopid}\t{foldid}\t{taxid}\t{len(seq)}\n")


scop_fasta = "scop95.fasta"

parse_scop_fasta(scop_fasta, "scop95.dat")


#mapping
f2fold = {}
sf2fold ={}
scopid2sf = {}
scopid2fold = {}

for line in open("scop95.dat"):
	splited_line = line.split('\t')
	scopid = splited_line[0]
	foldid = splited_line[1]
	splited_foldid = foldid.split('.')
	fold = '.'.join(splited_foldid[0:2])
	sf = '.'.join(splited_foldid[0:3])
	f = '.'.join(splited_foldid[0:4])
	#scopid to superfamily map
	if sf in scopid2sf:
		scopid2sf[sf].append(scopid)
	else:
		scopid2sf[sf] = []
		scopid2sf[sf].append(scopid)
	#scopid to fold map
	if fold in scopid2fold:
		scopid2fold[fold].append(scopid)
	else:
		scopid2fold[fold] = []
		scopid2fold[fold].append(scopid)
	#superfamily to fold map
	if sf in fold:
		sf2fold[fold].append(sf)
	else:
		sf2fold[fold] = []
		sf2fold[fold].append(sf)
	#family to fold map
	if f in fold:
		f2fold[fold].append(f)
	else:
		f2fold[fold] = []
		f2fold[fold].append(f)


print(len(scopid2sf.keys()))
print(len(scopid2fold.keys()))
# print(len(sf2fold.keys()))
# print(len(f2fold.keys()))

os.system("cut -f2 scop95.dat | awk -F. '{print $1,$2,$3,$4}' | uniq -c | sort -nk1 > nr_uniq_families.dat")
os.system("cut -f2 scop95.dat | awk -F. '{print $1,$2,$3}' | uniq -c | sort -nk1 > nr_uniq_superfamilies.dat")
os.system("cut -f2 scop95.dat | awk -F. '{print $1,$2}' | uniq -c | sort -nk1 > nr_uniq_folds.dat")
