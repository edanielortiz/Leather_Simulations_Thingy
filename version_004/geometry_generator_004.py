import string
import random
import os
import shutil as sht
import numpy as np
import replaceline as rl
from distutils.dir_util import copy_tree
import subprocess as sp
import matplotlib.pyplot as plt

print"Enter the length on x of the system (floating number in meters):"
world_x_length = input()
print"Enter the length on y of the system (floating number in meters):"
world_y_length = input()
print"Enter the length on z of the system (floating number in meters):"
world_z_length = input()
print"Enter number of collagen molecules in the system (MUST BE INTEGER):"
n_of_collagen_mol = input()
print"Enter number of elastin molecules in the system (MUST BE INTEGER):"
n_of_elastin_mol = input()
print"Enter number of electrons (MUST BE INTEGER):"
n_e = input()
print"Enter lowest energy value (floating number in KeV):"
low_E = input()
print"Enter highest energy value (floating number in KeV):"
high_E = input()
print"Enter the number of energy steps (MUST BE INTEGER and count the limits):"
stepss = input()
print"Enter the full path to Geant4 installation (PATH_TO_GEANT4/lib[64]/Geant4*):"
g4path = raw_input()
print"Enter the number of cores of your computer for shorter computation times."
print"In case you dont know this number please just enter zero (0)."
print"(number MUST BE INTEGER)"
n_cores = input()

world_x_length = world_x_length * 10000000000.0
world_y_length = world_y_length * 10000000000.0
world_z_length = world_z_length * 10000000000.0

#COLLAGEN

file_collagen = open("1cgd.pdb", "r")
input_text = file_collagen.readlines()

header_collagen = []
first_str_collagen = []
n_atm_collagen = []
second_str_collagen = []
x_input_collagen = []
y_input_collagen = []
z_input_collagen = []
last_str_collagen = []
conect_collagen = []

i = 0

while i <= 459:
    header_collagen.append(input_text[i])
    i += 1

#print i

while i <= 1537:

        #print i

        if i >= 1207:
            conect_collagen.append(input_text[i])
        else:
            if input_text[i][0:3] == "TER":
                first_str_collagen.append(input_text[i])
                n_atm_collagen.append(0)
                second_str_collagen.append(" ")
                x_input_collagen.append(1.0)
                y_input_collagen.append(1.0)
                z_input_collagen.append(1.0)
                last_str_collagen.append(" ")
            else:
                first_str_collagen.append(input_text[i][0:6])
                n_atm_collagen.append(int(input_text[i][6:11]))
                second_str_collagen.append(input_text[i][11:30])
                x_input_collagen.append(float(input_text[i][30:38]))
                y_input_collagen.append(float(input_text[i][38:46]))
                z_input_collagen.append(float(input_text[i][46:54]))
                last_str_collagen.append(string.replace(input_text[i][54::], "\n", "1\n"))

        i += 1

file_collagen.close()

#ELASTIN

file_elastin = open("1q5a.pdb", "r")
input_text = file_elastin.readlines()

first_str_elastin = []
n_atm_elastin = []
second_str_elastin = []
x_input_elastin = []
y_input_elastin = []
z_input_elastin = []
last_str_elastin = []
conect_elastin = []

i = 2271

#print i

while i <= 11098:

        #print i

        if i >= 11099:
            conect_elastin.append(input_text[i])
        else:
            if input_text[i][0:3] == "TER":
                first_str_elastin.append(input_text[i])
                n_atm_elastin.append(0)
                second_str_elastin.append(" ")
                x_input_elastin.append(0.0)
                y_input_elastin.append(0.0)
                z_input_elastin.append(0.0)
                last_str_elastin.append(" ")
            else:
                first_str_elastin.append(input_text[i][0:6])
                n_atm_elastin.append(int(input_text[i][6:11]))
                second_str_elastin.append(input_text[i][11:30])
                x_input_elastin.append(float(input_text[i][30:38]))
                y_input_elastin.append(float(input_text[i][38:46]))
                z_input_elastin.append(float(input_text[i][46:54]))
                last_str_elastin.append(string.replace(input_text[i][54::], "\n", "1\n"))

        i += 1

file_elastin.close()

#print first_str_elastin



file_geometry = open("geometry.pdb", "w+")

i = 0

while i < len(header_collagen):
    file_geometry.write(header_collagen[i])
    i += 1


#Oxigen for geometry fixation

x_positive_o = (world_x_length/2.0) - 1.5
y_positive_o = (world_y_length/2.0) - 1.5
z_positive_o = (world_z_length/2.0) - 1.5
n_atom_positive_o = 1


x_negative_o = ((world_x_length/2.0) - 1.5)*(-1.0)
y_negative_o = ((world_y_length/2.0) - 1.5)*(-1.0)
z_negative_o = ((world_z_length/2.0) - 1.5)*(-1.0)
n_atom_negative_o = 2


#ATOM      4  O   PRO A   1       6.313   2.168   7.380  0.77 18.55           O\n
first_str_o = "ATOM  "
second_str_o = "  O   PRO A   1    "
last_str_o = "  0.77 18.55           O\n"


#positive one
file_geometry.write(first_str_o)
file_geometry.write("\t\t{0:20d}\t\t".format(n_atom_positive_o))
file_geometry.write(second_str_o)
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(x_positive_o))
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(y_positive_o))
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(z_positive_o))
file_geometry.write("\t\t"+last_str_o)

#negative one
file_geometry.write(first_str_o)
file_geometry.write("\t\t{0:20d}\t\t".format(n_atom_negative_o))
file_geometry.write(second_str_o)
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(x_negative_o))
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(y_negative_o))
file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(z_negative_o))
file_geometry.write("\t\t"+last_str_o)


#COLLAGEN

x_randomized = 0.0
y_randomized = 0.0
z_randomized = 0.0



k = 1
atom_counter = 3

while k <= n_of_collagen_mol:

    flag = 0

    while flag == 0:

        x_randomized = (random.uniform(-1.0, 1.0))*world_x_length/2.0
        y_randomized = (random.uniform(-1.0, 1.0))*world_y_length/2.0
        z_randomized = (random.uniform(-1.0, 1.0))*world_z_length/2.0

        for j in range(len(x_input_collagen)):

            if (x_input_collagen[j] + x_randomized - 20.0) < ((-1.)*world_x_length/2.0):
                flag = 0
                break
            else:
                if (x_input_collagen[j] + x_randomized + 20.0) > (world_x_length/2.0):
                    flag = 0
                    break
                else:
                    flag = 1

            if (y_input_collagen[j] + y_randomized - 20.0) < ((-1.)*world_y_length/2.0):
                flag = 0
                break
            else:
                if (y_input_collagen[j] + y_randomized + 20.0) > (world_y_length/2.0):
                    flag = 0
                    break
                else:
                    flag = 1

            if (z_input_collagen[j] + z_randomized - 20.0) < ((-1.)*world_z_length/2.0):
                flag = 0
                break
            else:
                if (z_input_collagen[j] + z_randomized + 20.0) > (world_z_length/2.0):
                    flag = 0
                    break
                else:
                    flag = 1

    #print "collagen point found"
    i = 0

    while i < len(first_str_collagen):

        if first_str_collagen[i][0:3] == "TER":
            file_geometry.write(first_str_collagen[i])

        else:
            file_geometry.write(first_str_collagen[i])
            file_geometry.write("\t\t{0:20d}\t\t".format(atom_counter))
            file_geometry.write(second_str_collagen[i])
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(x_input_collagen[i] + x_randomized))
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(y_input_collagen[i] + y_randomized))
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(z_input_collagen[i] + z_randomized))
            file_geometry.write("\t\t"+last_str_collagen[i])

        i += 1
        atom_counter += 1

    k += 1

#ELASTIN

print "collagen done"

k = 1

#print len(first_str_elastin)

while k <= n_of_elastin_mol:

    flag = 0

    #print"first elastin"

    while flag == 0:
        x_randomized = (random.uniform(-1.0, 1.0))*world_x_length/2.0
        y_randomized = (random.uniform(-1.0, 1.0))*world_y_length/2.0
        z_randomized = (random.uniform(-1.0, 1.0))*world_z_length/2.0
        #print "what"
        #print len(x_input_elastin)
        for j in range(len(x_input_elastin)):

            #print "for elastin"

            if (x_input_elastin[j] + x_randomized - 20.0) < ((-1.)*world_x_length/2.0):
                flag = 0
                #print "break 1"
                break
            else:
                if (x_input_elastin[j] + x_randomized + 20.0) > (world_x_length/2.0):
                    flag = 0
                    #print "break 2"
                    break
                else:
                    flag = 1

            if (y_input_elastin[j] + y_randomized - 20.0) < ((-1.)*world_y_length/2.0):
                flag = 0
                #print "break 3"
                break
            else:
                if (y_input_elastin[j] + y_randomized + 20.0) > (world_y_length/2.0):
                    flag = 0
                    #print "break 4"
                    break
                else:
                    flag = 1

            if (z_input_elastin[j] + z_randomized - 20.0) < ((-1.)*world_z_length/2.0):
                flag = 0
                #print "break 5"
                break
            else:
                if (z_input_elastin[j] + z_randomized + 20.0) > (world_z_length/2.0):
                    flag = 0
                    #print "break 6"
                    break
                else:
                    flag = 1

    #print "elastin point found"
    i = 0

    while i < len(first_str_elastin):

        if first_str_elastin[i][0:3] == "TER":
            file_geometry.write(first_str_elastin[i])

        else:
            file_geometry.write(first_str_elastin[i])
            file_geometry.write("\t\t{0:20d}\t\t".format(atom_counter))
            file_geometry.write(second_str_elastin[i])
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(x_input_elastin[i] + x_randomized))
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(y_input_elastin[i] + y_randomized))
            file_geometry.write("\t\t\t{0:20.3f}\t\t\t".format(z_input_elastin[i] + z_randomized))
            file_geometry.write("\t\t"+last_str_elastin[i])

        i += 1
        atom_counter += 1

    k += 1

print "elastin done"

i = 0

while i < len(conect_collagen):

    file_geometry.write(conect_collagen[i])
    i += 1

file_geometry.close()

print "file geometry.pdb done"
#os.mkdir()


#GENERATION OF FOLDERS FOR THE SIMS AND COMPILATION OF THE APPLICATION

sdirname1 = "leather_sim_"
sdirname2 = "_keV"

#Creation of the simulation dirs

curr_dir = os.getcwd()
parent_dir = curr_dir+"/simulations"
srcfiles_dir = curr_dir+"/PDB_OF"
os.mkdir(srcfiles_dir)



rl.replaceline_and_save(fname = curr_dir+"/pdb4dna_sourcefiles/src/DetectorConstruction.cc",
        findln = "G4double world_x_size =",
        newline = "G4double world_x_size = "+str(world_x_length)+"*1*angstrom;")

rl.replaceline_and_save(fname = curr_dir+"/pdb4dna_sourcefiles/src/DetectorConstruction.cc",
        findln = "G4double world_y_size =",
        newline = "G4double world_y_size = "+str(world_y_length)+"*1*angstrom;")

rl.replaceline_and_save(fname = curr_dir+"/pdb4dna_sourcefiles/src/DetectorConstruction.cc",
        findln = "G4double world_z_size =",
        newline = "G4double world_z_size = "+str(world_z_length)+"*1*angstrom;")




#Copy of the geometry file and compilation of the Geant4 Application


sht.copyfile(curr_dir+"/geometry.pdb", srcfiles_dir+"/geometry.pdb")
os.chdir(srcfiles_dir)
cmake_g4dir = "-DGeant4_DIR="+g4path
sp.call(["cmake", cmake_g4dir, "../pdb4dna_sourcefiles"])

if n_cores != 0:
    cores_call = "-j"+str(n_cores)
    sp.call(["make", cores_call])
else:
    sp.call(["make"])

os.chdir(curr_dir)
os.mkdir(parent_dir)

E_span = np.linspace(low_E, high_E, num=stepss)  # Energy swipe
#n_e = 10000							  Number of electrons (events)

j = 0

for E in E_span:
	j += 1
	fullsdirname = "%s%.4f%s"%(sdirname1, E, sdirname2)


	#Creation of the Simulation Directory
	print"Creating %s"%(fullsdirname)
	sim_dir = parent_dir+"/"+fullsdirname
	os.mkdir(sim_dir)
	copy_tree(srcfiles_dir, sim_dir)


	#Format of the Simulation
	#	pdb4dna.in
	rl.replaceline_and_save(fname = sim_dir+"/pdb4dna.in",
			findln = "/gun/energy ",
			newline = "/gun/energy "+str(E)+" keV")

	rl.replaceline_and_save(fname = sim_dir+"/pdb4dna.in",
			findln = "/run/beamOn ",
			newline = "/run/beamOn "+str(n_e))


	#Run the just created simulation
	os.chdir(sim_dir)
	sp.call(["./pdb4dna", "-mac", "pdb4dna.in"]) #, "-gui", "qt"])


# +----------------------------------------------------------------------------+
# |							BEGIN DATA ANALYSIS 							   |
# +----------------------------------------------------------------------------+

# ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

# COLLAGEN PLOTTING

os.chdir(curr_dir)
os.mkdir(curr_dir+"/results")

n_of_collagen_broken_oo = np.array([])
n_of_collagen_broken_o = np.array([])
n_of_collagen_non_o_broken = np.array([])
n_of_elastin_broken = np.array([])
n_of_water_inter = np.array([])
n_of_no_interactions = np.array([])
array_of_depths = np.array([])

j = 0
for E in E_span:
    j += 1
    fullsdirname = "%s%.4f%s"%(sdirname1, E, sdirname2)
    sim_dir = parent_dir+"/"+fullsdirname


    energy_relation = np.zeros((1,6))
    input_from_file = np.zeros((1,6))

    #input_from_file = np.loadtxt(sim_dir+"/energy_and_bond_type.txt")
    #input_from_file = np.genfromtxt(sim_dir+"/energy_and_bond_type.txt",delimiter='\t\t\t',dtype='float64')
    f_fix = open(sim_dir+"/energy_and_bond_type.txt", "r")
    in_text = f_fix.readlines()
    f_fix.close()

    conta = 1

    for i in range(len(in_text)):

        if len(np.fromstring(in_text[i], sep="\t\t\t")) != 6:
            continue
        else:

            if conta == 1:

                input_from_file = np.reshape(np.fromstring(in_text[i], sep="\t\t\t"), (1,6))
                conta += 1

            else:

                input_from_file = np.append(input_from_file, np.reshape(np.fromstring(in_text[i], sep="\t\t\t"), (1,6)), 0)
                conta += 1


    n_inter = float(len(np.loadtxt(sim_dir+"/interaction_counter.txt")))

    k = 0

    for i in range(len(input_from_file[:])):		#Reduce the size of the array to only enev

        if input_from_file[i][0] != 0:

            k += 1

            if k == 1:

                energy_relation	 = np.reshape(input_from_file[i][:], (1,6))

            else:

                energy_relation = np.append(energy_relation, np.reshape(input_from_file[i][:], (1,6)), 0)

    bond_info = np.array([0,0,0,0,0,0])
    for i in range(len(energy_relation[:])):

        if energy_relation[i][4] == 1.:

            if energy_relation[i][1] == 1.:

                if energy_relation[i][0] > 6.1984:

                    bond_info[0] += 1.

                elif (energy_relation[i][0] < 6.1984) and (energy_relation[i][0] > 2.4151):

                    bond_info[1] += 1.

                else:

                    if energy_relation[i][0] > 0.:

                        bond_info[2] += 1.

                    else:

                        bond_info[5] += 1.

            elif energy_relation[i][1] == 2.:

                if energy_relation[i][0] > 0.:

                    bond_info[2] += 1.

                else:

                    bond_info[5] += 1.

            else:

                if energy_relation[i][0] == 0.:

                    bond_info[5] += 1.

                else:

                    bond_info[4] += 1.

        elif energy_relation[i][4] == 2.:

            if energy_relation[i][0] > 3.58640:

                bond_info[3] += 1.

            else:

                if energy_relation[i][0] == 0.:

                    bond_info[5] += 1.

                else:

                    bond_info[4] +=1.
        else:

            if energy_relation[i][0] == 0.:

                bond_info[5] += 1.

            else:

                bond_info[4] += 1.


	#bond_info 0 is broken c = o in collagen
    #bond_info 1 is broken c - o in collagen
    #bond_info 2 is other bond broken in collagen
    #bond_info 3 is bond broken in elastin
    #bond_info 4 is interactions in water
    #bond_info 5 is no interactions

    n_of_accounted_interactions = bond_info[0] + bond_info[1] + bond_info[2] + bond_info[3] + bond_info[4] + bond_info[4] + bond_info[5]
    n_of_other_no_interactions = n_e - n_of_accounted_interactions
    bond_info[5] = bond_info[5] + n_of_other_no_interactions

	#bond_info.append(n_of_broken_oo)
	#bond_info.append(n_of_broken_o)
	#bond_info.append(n_of_non_o)
	#bond_info.append(n_of_nonbroken)

	#bond_info[:] = [y / float(n_e) for y in bond_info]
    x_hist = np.array([1,2,3,4,5,6])


    x_labels = ["C=O Collagen","C-O Collagen","Other Collagen", "Bonds in Elastin", "Water", "No interactions"]

    fig1 = plt.figure(figsize = (10,7))
    ax1 = fig1.add_subplot(111)

    rects1 = ax1.bar(x_hist, bond_info/float(n_e))

    ax1.set_xticks(x_hist)
    ax1.set_xticklabels(x_labels)
    ax1.set_xlabel("Types of Interactions")
    ax1.set_ylabel("Frequency")
    ax1.set_ylim(0,1.0)
    ax1.set_title("Interactions with the Collagen Molecule at "+str(E)+" keV")



    fig1.savefig(curr_dir+"/results/bonds_histogram_at_"+str(E)+"keV.png")
    plt.close()

	#Zoomed in histogram

    fig2 = plt.figure(figsize = (10,7))
    ax2 = fig2.add_subplot(111)

    rects2 = ax2.bar(x_hist[0:4], bond_info[0:4]/float(n_e))

    ax2.set_xticks(x_hist[0:4])
    ax2.set_xticklabels(x_labels[0:4])
    ax2.set_xlabel("Types of Interactions")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Interactions with the simulated leather at "+str(E)+" keV")



    fig2.savefig(curr_dir+"/results/zoomed_bonds_histogram_at_"+str(E)+"keV.png")
    plt.close()


	#Save the number of bond interactions for this step for entire energy simulations
    n_of_collagen_broken_oo = np.append(n_of_collagen_broken_oo, bond_info[0])
    n_of_collagen_broken_o = np.append(n_of_collagen_broken_o, bond_info[1])
    n_of_collagen_non_o_broken = np.append(n_of_collagen_non_o_broken, bond_info[2])
    n_of_elastin_broken = np.append(n_of_elastin_broken, bond_info[3])
    n_of_water_inter = np.append(n_of_water_inter, bond_info[4])
    n_of_no_interactions = np.append(n_of_no_interactions, bond_info[5])

	#n_of_collagen_broken_oo.append(bond_info[0]/float(n_e))
	#n_of_collagen_broken_o.append(bond_info[1]/float(n_e))
	#n_of_collagen_non_o_broken.append(bond_info[2]/float(n_e))
	#n_of_elastin_broken.append(bond_info[3]/ float(n_e))
    #n_of_water_inter.append(bond_info[4]/ float(n_e))
    #n_of_no_interactions.append(bond_info[5]/ float(n_e))


    num_bins = 100
    fig4 = plt.figure(figsize = (15,15))
    ax4 = fig4.add_subplot(111)
    n, bins, patches = ax4.hist((energy_relation[:,5] - (world_z_length/2.0))*(-1), num_bins, facecolor='blue', alpha = 0.5)
    ax4.set_xlabel("Depth in Angstrom")
    ax4.set_ylabel("Frequency")
    ax4.set_title("Depth Reach with "+str(E)+" keV electron beam")
    fig4.savefig(curr_dir+"/depth_reach_at_"+str(E)+"keV.png")



# ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#Plotting of Energy vs frequency of interactions

#Number of broken C=0

fig3 = plt.figure(figsize = (10,7))
ax = plt.axes()
ax.plot(E_span, n_of_collagen_broken_oo, label = "Broken Collagen C=O")
ax.plot(E_span, n_of_collagen_broken_o, label = "Broken Collagen C-O")
ax.plot(E_span, n_of_collagen_non_o_broken, label = "Other Collagen Bonds Broken")
ax.plot(E_span, n_of_elastin_broken, label = "Bonds Broken in Elastin")
ax.plot(E_span, n_of_water_inter, label = "Interactions in Water")
#ax.plot(E_span, n_of_no_interactions, label = "No Interactions")
ax.set_ylim(bottom = 0.)
ax.set_xlabel("Energy in keV")
ax.set_ylabel("Frequency")
ax.set_title("Interactions with the simulated leather vs Energy")
ax.legend()
fig3.savefig(curr_dir+"/results/frequency_of_interactions_vs_energy.png")





#plt.show()
