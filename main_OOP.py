# This is a Python implementation of Tau-leap gillespie stochastic simulaitonwith partial deterministic reactions
# Refer to  readme.ipynb for backgrounds

import numpy as np
import method_oop as method			# method_oop has global variable, class and functions defined
import sys


def main():
	# run simulations for num_sim and dosage profile
	# Global parameters in method_oop can also be altered for other scenarios

	num_sim = 100                                   # No. of simulations
	dosage_prof = [0, 1, 2, 3]                      # Dosage profile ID - in parallel with method.dosage_prof_dic

	for i in range(num_sim):
		sim_ID = i
		method.Type_dose = len(dosage_prof)
		File_name = file_open(sim_ID)
		for j in dosage_prof:
			f= File_name[j]
			run_sim_spec(sim_ID, j, f)
			f.close()


def run_sim_spec(sim_ID, j, f, **kwargs):

	# Create instance of gillespie simulation and run the simulation, accordingly with the initial conditions

	if len(kwargs) ==0:
		# Initial condition, if there's no input
		t_curr = 0.
		gill_sim= method.Gill_sim()
		gill_sim.S = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1E12, 1E12, 1E12, 1E12], [3E2, 3E2, 3E2, 3E2]])[:,j]  # state_matrix S: update once every */
		gill_sim.init_F()

	else:
		# If there's an input
		try:
			t_curr= kwargs['t_domain'][0]
			method.t= kwargs['t_domain'][1]
		except KeyError:
			print "no t_domain key"
		try:
			method.Type_dose = 1
			dosage_prof = kwargs['dosage_prof']
		except KeyError:
			print "no dosage profile key"
		try:
			gill_sim.S= kwargs['init_cond']
		except KeyError:
			print "no initial species condition"

	f.write("time\tAg\tB\tIC\tIgM\tIgG\tTfh\n")
	gill_sim.calc_eps()

	# Simulation t<6d has to add Ag dosage everyday
	# for t_end<= 6
	if method.t < method.dose_per:
		i = np.ceil(t_curr)
		gill_sim.update(i, j, [t_curr,method.t], f)

	#t_start < 6d && t_end>6ds
	elif t_curr<method.dose_per and method.t>method.dose_per:
		i = np.ceil(t_curr)
		gill_sim.update(i, j, [t_curr, method.t], f)

	# t_start > 6d
	else:
		i = np.ceil(t_curr)
		gill_sim.update(i, j, [t_curr, method.t], f)

def file_open(sim_ID):
	# Open up dosing schedule file, state files  */
	State_EI = open("State_EI%d.txt" % (sim_ID),'w')
	State_ED = open("State_ED%d.txt" % (sim_ID),'w')
	State_Const = open("State_Const%d.txt" % (sim_ID),'w')
	State_PB = open("State_PB%d.txt" % (sim_ID),'w')
	File_name = [State_EI, State_ED, State_Const, State_PB]
	return File_name

def file_close(file_list):
	for f in file_list:
		f.close()
	return

if __name__ == "__main__":
	import cProfile
	cProfile.run('main()',sort="tottime", filename= "profileStats")
