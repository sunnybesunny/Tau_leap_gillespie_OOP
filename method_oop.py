# This is a Python implementation of Tau-leap gillespie stochastic simulaiton with partial deterministic reactions

import numpy as np
import random
import sys
import method
import os

import math
import matplotlib.pyplot as plt


###### Global parameters ######
t = 6				    # total t span, modified by link */
dose_per = 6 			# dosing period : seven dosage is given on day 0,1,2,3,4,5,6
Type_spe = 6  			# Type of species: C_Ag,C_B*,C_Ic,C_IgG,C_IgM,T_fH */
Type_dose = 1 			# Type of dosing: input one by one
gil_num_RXN = 6 		# # of rxn*/
# a = 2  				# Version 1 a=2*/
Volume = 0.000001 		# volume of the lymph node, 1uL  */
n_A = 6.02e23           # Avogadro's number
dosage_prof_dic = {'EI': 0, 'ED': 1, 'Const': 2, 'Boost': 3} # 4 dosage kinetics tested


######  Tau-leap parameters ######
nc = 10 				# of critical rxns- labeled as critical, if within 10 rxns from depletion*/
eps = 0.05 				# epsilon */

###### deterministic update parameters ######
num_RXN_det = 2         # number of reactios that are being updated deterministically
num_spec_det = 2        # number of species that are being updated deterministically

###### simulation function ######


def Poirand(lamb):
	# Create poission random number.
	# Approximate to normal rv if lambda > 1000

	if lamb>1000:
		k= int(np.random.normal(lamb, np.sqrt(lamb)))
		return k
	else:
		try:
			k= np.random.poisson(lamb)
			return k
		except ValueError:
			k = 0
			produ = random.random()
			while produ > np.exp(-lamb):
				produ = produ * random.random()
				k = k + 1
			return k

def beta(ti):
	# Affinity - receptor binding affinity determining functions; affinity increases upto 100 folds in 21 days
	# In combined simulaiton, this is replaced with affinity determined by affinity maturation process
	if ti<=21:
		return ((8.64E7) / (Volume * n_A)) * (1 + (99 / 21) * ti)
	else:
		return 8.64E9 / (Volume * n_A)

###### Simulaiton class ######

class Gill_sim:
	# Gillespie simulation cycle; variables and functions
	def __init__(self):

		###### Tau-leap variable ######
		self.V = np.array( [[-1, 0, 0, 0, 0, 0], [-1, 0, 1, 0, -1, 0], [-1, 0, 1, -1, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, -1, 0, 0],
					  [0, 0, 0, 0, -1, 0]])                                 # Stochiometry
		self.C = np.zeros((gil_num_RXN,))                                   # rate constant matrix c */
		# self.V_s = np.zeros((gil_num_RXN, Type_spe))
		# self.C_s = np.zeros((gil_num_RXN,))
		self.H_s = np.zeros((gil_num_RXN,))                                 # Permutation matrix H  */
		self.S = np.zeros((Type_spe,), dtype=np.float32)                    # State_matrix S: */
		self.P = np.zeros((gil_num_RXN,))                                   # Propensity */
		self.det_P = np.zeros((num_RXN_det,))                               # Deterministic reaction propensity
		self.F = np.zeros((dose_per, Type_dose))                            # Vaccine dosage
		self.t_range = np.zeros((2,))                                       # time range
		self.det_rate = {'S1': 0.0233, 'S2': 0.0233, 'S3': 0.048, 'S4': 0}  # Reaction rate of deterministic reactions
		self.HOR = np.zeros((Type_spe,))                                    # Highest order of reaction species i is involved  */
		self.eps_array = np.zeros((Type_spe,))                              # epsilon array
		self.I_nc = np.zeros((gil_num_RXN,), dtype=int)                     # Index of non critical reactions */
		self.I_c = np.zeros((gil_num_RXN,), dtype=int)                      # Index of critical reactions*/
		# / SpecInRxn[num_RXN][] */
		self.num_crit = 0                                                   # Number of critical reactions
		self.num_non_crit = 0                                               # Number of non critical reactions
		# num_reactant_init = 0  # # of reactants for t<6 */
		self.num_reactant = Type_spe                                        # Number of reactants
		self.num_gil_cycle = 0                                              # Number of gillespie cycle
		self.flag_early_term = 0                                             # flag for early termination, default = -1, Species in the loop depleted =0, species in the deterministic approximation depleted =1
		# of reactants for t>6
		# n_Ag = S[0]
		# n_PC = S[1]
		# n_IC = S[2]
		# n_IgG = S[3]
		# n_IgM = S[4]
		# n_Tfh = S[5]
		self.n_GCB=0                                                        # Number of Germical center B cell


	class RXN(object):
		# Reaction class includes group of characteristic variables
		def __init__(self, RXN_num, num_specs, SpecsinRXN, crit):
			self.RXN_num = RXN_num
			self.num_specs = num_specs
			self.SpecsinRXN = SpecsinRXN  # np array  of SpecsinRXN[Type_spe]
			self.crit = crit

	RXN_info = [RXN(0, 1, [0], -1), RXN(1, 2, [0, 4, 2], -1), RXN(2, 2, [0, 3, 2], -1),
	            RXN(3, 1, [1,3], -1), RXN(4, 1, [3], -1), RXN(5, 1, [4], -1)]

	class MyError(Exception):
		def __init__(self, value, expr):
			self.value = value
			self.expr = expr

		def __str__(self):
			return repr(self.expr)


	###### Class function ######

	def check(self, bool, msg):

		# if False is passed, error is raised
		if bool:
			raise Gill_sim.MyError(bool, msg)
		else:
			pass

	def calc_C(self,ti):
		# Calculate reaction rates
		# unit: day
		if ti < 6:
			self.C = np.array([4, (8.64E7) / (Volume * n_A), 0, 0, 0.0233,0.0233 ])
		# rate constant matrix c(/day) */
		else:
			self.C = np.array(
				[4, (8.64E7) / (Volume * n_A), beta(ti), 2.5E6, 0.0233, 0.0233])
		self.check(len(np.where(self.C < 0)[0]), "error:negative C matrix")
		# self.check(len(np.where(method.V < -2)[0]), "error: V matrix smaller than -2")

	def calc_eps(self):
		# The highest order of rxns and calculate eps value accordingly
		self.HOR = np.array([2, 1, 1, 2, 2, 1])
		for i in range(Type_spe):
			# print self.HOR[i] < 0
			self.check(self.HOR[i] < 0, "error: HOR[%d]negative" %(i))
			if self.HOR[i] == 1:
				self.eps_array[i] = eps
			elif self.HOR[i] == 2:
				self.eps_array[i] = eps / 2
			else:
				self.eps_array[i] = 0  # when it is non reacting species

	def calc_H(self):
		# Calculate permutation
		# r1: Ag(0) decay: C_Ag[0] -> 0 (k: 1/day)
		# r2: IC(2) formation from IgM(4): C_Ag[0] +C_IgM[4]  -> C_IC[2] (k: 10^3/(M*S)*8.64*10^4(day/sec))
		# r3: IC(2) formation from IgG(3):  C_Ag[0] +C_IgG[3]  -> C_IC[2] (k: 10^5/(M*S)*8.64*10^4(day/sec))
		# r4: IgG(3) prod: C_PC[1] -> C_PC+C_IgG[3] (k:2.5E6 Pcells/day from model ): multiply by vol factor(lymph vol/serum vol= 0.5*10^-3) for comparision with serum measured values
		# r5: IgG[3] decay: IgG[3] -> 0 (k:0.0233/day)
		# r6: IgM[4] decay: IgM --> 0 ( k: 0.02333/day )

		for i in range(Type_spe):
			self.check(self.S[i] < 0, "error: S[%d]negative, %f" % (i, self.S[i]))
		self.H_s[0] = self.S[0]
		self.H_s[1] = self.S[4]*self.S[0]
		self.H_s[2] =np.int(0.025* self.S[3]*self.S[0])
		self.H_s[3] = self.S[1]
		self.H_s[4] = self.S[3]
		self.H_s[5] = self.S[4]
		self.check(len(np.where(self.H_s < 0)[0]), "error:self.H_s negative,%e \n")

	def calc_P(self,ti):
		# Calculate propensity
		self.calc_H()
		self.calc_C(ti)
		self.check(len(np.where(self.H_s < 0)[0]), "error:self.H_s negative,%e \n")
		self.check(len(np.where(self.C < 0)[0]), "error:self.C negative,%e \n")
		self.P[0:gil_num_RXN] = self.H_s[0:gil_num_RXN] * self.C[0:gil_num_RXN]
		self.check(len(np.where(self.P < 0)[0]), "error:self.P negative,%e \n")

	def sum_P(self):
		# Calculate the sum of propensity
		sum_prop = self.P[0:gil_num_RXN].sum()
		return sum_prop

	def sum_det_P(self,t_old, t_before):
		# Make deterministic update every 0.5 day or if it's only possible to carry deterministic update
		# Two conditions: Propensity, residual species
		# First condition: checks if there're enough speces
		# Second condition: checks if it's time for det_update : whether delta_adi time has passed or sum_P() ==0 so that only det_update is possible
		# return True(rxn will not occur) , False(reaction will occur)
		delta_adi = 0.5
		verdict = True # rxn does not occur

		if self.n_GCB !=0:
			if (self.S[2] >= int(self.S[2] * delta_adi * self.det_rate['S1'])) and (
						self.S[5] >= int(delta_adi * self.det_rate['S3']*(1-self.S[5]/self.n_GCB))) :
				if t_old - t_before >= delta_adi or self.sum_P() == 0.:
					verdict = False
		else:
			if (self.S[2] >= int(self.S[2] * delta_adi * self.det_rate['S1'])):
				if t_old - t_before >= delta_adi or self.sum_P() == 0.:
					verdict = False
		return verdict

	def calcNCRxnNum_Update(self,tau):
		# Calculate number of non-critical rxns occuring
		# update the species accordingly
		for j in range(self.num_non_crit):
			index_2 = self.I_nc.astype(int)[j]
			k = Poirand(self.P[index_2] * tau)
			self.S[0:Type_spe] += np.multiply(self.V[index_2][:], k).astype(np.float32)

	def calcTauNonC(self):
		# Compute tau using (24): maximum leap time allowed by the Leap Condition for the non-critical reactions
		tau = sys.float_info.max
		if (self.num_non_crit == 0):
			return tau
		self.check(len(np.where(self.I_nc[range(self.num_non_crit)] < 0)[0]), "error:self.I_nc negative")
		self.check(len(np.where(self.I_nc[range(self.num_non_crit)] > gil_num_RXN)[0]),
		             "error:self.I_nc more than Type_spe")
		I_rs_set = [set(x.SpecsinRXN) for x in [method.RXN_info[z] for z in self.I_nc[range(self.num_non_crit)]] if
		            np.all(
			            self.S[x.SpecsinRXN] > 0) ]
		I_rs = set.union(*I_rs_set)
		index_1 = list(I_rs)
		rxn_ind = self.I_nc[range(self.num_non_crit)]
		nom = np.maximum(self.eps_array[index_1] * self.S[index_1], 1)
		denom = np.fabs(np.einsum('ij,i->j', self.V[np.ix_(rxn_ind, index_1)], self.P[rxn_ind]))
		denom2 = np.einsum('ij,i->j', self.V[np.ix_(rxn_ind, index_1)] * self.V[np.ix_(rxn_ind, index_1)],
		                   self.P[rxn_ind])
		ind = np.where(denom != 0)[0]
		tau = min(np.minimum(np.divide(nom[ind], denom[ind]), np.divide(nom[ind] * nom[ind], denom2[ind])))
		return tau

	def calcTauC(self):
		# Compute tau using (24): maximum leap time allowed by the Leap Condition for the non-critical reactions
		tauc = 0  # next critical rxn time*/
		sum_prop_crit = 0.0  # sum of propencities*/
		if (self.num_crit == 0):
			tauc = sys.float_info.max
		else:
			self.check(np.where(self.I_c[range(self.num_crit)] <= 0)[0], "error:self.I_c negative")
			sum_prop_crit = self.P[self.I_c[range(self.num_crit)]].sum()
			rand = random.random()
			#				printf("%10.8f  ", genrand64_real3());
			if (sum_prop_crit > 0):
				tauc = -np.log(rand) / sum_prop_crit
			elif (sum_prop_crit == 0):
				tauc = sys.float_info.max - 1
			else:
				pass
		return tauc, sum_prop_crit

	def RxnPartition(self):
		# Partition rxns into critical and non critical reactions
		# If all critical/non critical, I_c=[all rxn index]/I_nc=[np.inf]*method.num_RXN
		# Classify as critical any Rj for which aj(x) > 0 and which is within nc firings of exhausting any reactant.
		# Classify all other Rj non-critical.
		self.num_crit = 0
		self.num_non_crit = 0
		self.I_nc[:] = np.inf
		self.I_c[:] = np.inf
		# print np.ravel(np.fabs(method.V)).max()
		for i in range(gil_num_RXN):
			if len(np.where(self.S[method.RXN_info[i].SpecsinRXN[:]] < np.ravel(np.fabs(self.V)).max()*nc)[0]) or self.P[
				method.RXN_info[i].RXN_num] == 0:
				method.RXN_info[i].crit = 1
				self.I_c[self.num_crit] = method.RXN_info[i].RXN_num
				#                    prf("critical rxns %d\n",I_c[num_crit]);
				self.num_crit += 1
			else:
				method.RXN_info[i].crit = 0
				self.I_nc[self.num_non_crit] = method.RXN_info[i].RXN_num
				self.num_non_crit += 1

		self.check(self.num_crit > gil_num_RXN, "error:num_crit pass num_RXN")
		self.check(self.num_non_crit > gil_num_RXN, "error:num_non_crit pass num_RXN")
		self.check(self.num_non_crit + self.num_crit != gil_num_RXN, "error:total sum not equal to num_RXN")


	def select_reaction(self,sum_prop, randN, RXN_list, num_RXN):
		# Choose reaction within critical regime
		reaction = -1
		sp = 0.0
		randN = randN * sum_prop
		for i in range(num_RXN):
			index = RXN_list[i]
			sp += self.P[index]
			if (randN < sp):
				reaction = index
				break
		self.check(reaction < 0, "error:negative reaction selected")
		self.check(reaction > gil_num_RXN, "error:reaction index over number of rxns")
		return reaction

	def init_F(self):
		# Get Ag dosage
		cwd = os.path.dirname(os.path.abspath(__file__))
		F_dosing_schedule = open("%s/F_dosing_schedule.txt" %cwd, "r")
		self.check(F_dosing_schedule == None, "file open error\n")
		self.F = np.loadtxt(F_dosing_schedule)
		F_dosing_schedule.close()

	def det_step( self, t_old, t_before):
		# deterministic update of species. updated every delta_adi, unless gillespie rxns are propensity zero
		#t_old : current time tag / t_before: time at which the last deterministic update has made
		delta_adi = 0.5

		# Case 1: no gilespie update, kept only deterministic updates
		if self.sum_P() == 0:
			t_old += delta_adi
			delT = delta_adi
			t_before = t_old
		# Case 2: there's gilespie update and when it is time to update
		elif t_old - t_before >= delta_adi:
			delT = t_old - t_before
			t_before = t_old
		else:
			delT = 0
		self.S[2] = max(0, (self.S[2] - int(self.S[2] * delT * self.det_rate['S1'])))
		# apply logistic growth
		if self.n_GCB!= 0:
			cap =  self.n_GCB
			self.S[5] = min(cap,(self.S[5] + int(delT * self.det_rate['S3']*(1-self.S[5]/cap))))
		else:
			self.S[5]= 0
		# f.write("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n" % (
		# 	t_old, method.S[0], method.S[1], method.S[2], method.S[3], method.S[4],
		# 	method.S[5]))
		return t_old, t_before

	def gillespie_cycle(self,f):
		# Differentiate the case of simulation and execute
		t_before = self.t_range[0]
		t_old = self.t_range[0]  # initialize time pt*/
		t_end = self.t_range[1]  # end point for gillespie algorithm*/
		self.calc_C(t_old)
		while (t_old <= t_end):
			S_old = self.S[0:Type_spe].copy()
			self.calc_P(t_old)
			verdict = self.sum_det_P(t_old, t_before)
			f.write("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n" % (
				t_old, self.S[0], self.S[1], self.S[2], self.S[3], self.S[4], self.S[5]))
			# f_aff.write("%.5e\t%.5e\n" % (
			# 	t_old,self.C[2]))
			self.num_gil_cycle += 1
			# initialize ( calculated calc_H automatically)

			# Case 1: Propensity ==0 for both Gillespie and deterministic update
			if self.sum_P() == 0 and verdict:  # screen for early termination of the cycle
				self.flag_early_term = -1
				#Var.t_sim = t_old
				break

			# Case 2,3: Non-zero Gilespie RXN propensity with Zero propensity det RXN (case2)/ non-zero propensity det (case3)
			elif self.sum_P() != 0:
				t_old = self.gillespie_step( t_old)
				if self.sum_det_P(t_old, t_before):  # if det_rxn propensity is ze
					pass
				else:
					t_old, t_before = self.det_step( t_old, t_before)
			# Case 4: Zero propensity gilespie RXN && non-zero Prop det RXN
			else:
				t_old, t_before = self.det_step(t_old, t_before)
		if t_old > t_end:
			self.S[0:Type_spe] = S_old #prevent species update made by time step outside of bounds
		# Var.t_sim = t_end

	def gillespie_step(self,t_old):
		# Explicit tau-leap
		self.RxnPartition()  # calculate method.num_crit/ method.num_non_crit/method.I_nc[:]/method.I_c[:]
		taunc = self.calcTauNonC()
		self.check(taunc < 0, "error:self.calcTauNonC negative")
		tauc, sum_prop_crit = self.calcTauC()

		if (tauc < taunc):
			tau = tauc
			r = random.random()
			reaction = self.select_reaction(sum_prop_crit, r, self.I_c, self.num_crit)
			# */
			self.S[0:Type_spe] += self.V[reaction, 0:Type_spe]
			self.check(np.where(self.S < 0)[0], "error:self.S negative")
		else:
			tau = taunc
			self.calcNCRxnNum_Update(taunc)
		t_old += tau

		return t_old

	def add_dose(self, k,i):
		# Add antigen by vaccine dosage
		self.S[0] += self.F[k, i]

	def update(self,i,k,t_domain,f):
		# Give antigen dosage and run gillespie cycle
		t_curr =t_domain[0]
		t_end = t_domain[1]
		while (i <= min(dose_per, np.floor(t_end))):
			# method.check(ac.calc_H(k) == 0, "error: calc_H(%d)" % k)
			self.t_range[0] = t_curr
			self.t_range[1] = i
			self.gillespie_cycle(f)
			# print self.num_gil_cycle
			if self.flag_early_term == -1:
				print "Species depleted"
				return
			t_curr = i
			self.add_dose(k, int(i))
			i += 1
		self.t_range[0] = t_curr
		self.t_range[1] = t_end
		self.gillespie_cycle(f)
		# print method.num_gil_cycle
		if method.flag_early_term == -1:
			print "Species depleted"
			return
