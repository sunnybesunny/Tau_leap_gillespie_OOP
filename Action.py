import numpy as np
import method
import sys
def calc_eps():
	method.HOR = np.array([2,1,1,2,2,1])
	for  i in range(method.Type_spe):
		method.check(method.HOR[i]<0, "error: HOR[%d]negative" %i)
		if method.HOR[i] == 1 :
			method.eps_array[i]= method.eps
		elif  method.HOR[i] == 2:
			method.eps_array[i]= method.eps/2
		else:
			method.eps_array[i]= 0   #when it is non reacting species

def Calc_H(j):
	method.check(len(np.where(method.S[:,j]<0)[0]),"error:S negative")
	method.H_s[0] = method.S[0,j]
	method.H_s[1] =method.S[4,j]*method.S[0,j]
	method.H_s[2] =method.S[3,j]*method.S[0,j]
	method.H_s[3] =method.S[1,j]
	method.H_s[4] =method.S[5,j]
	method.H_s[5] = method.S[4,j]
	method.H_s[6] = method.S[3,j]
	method.H_s[7] = method.S[2,j]
	method.check(len(np.where(method.H_s < 0)[0]), "error:method.H_s negative,%e \n")

def Calc_P(scheme):
	Calc_H(scheme)
	method.check(len(np.where(method.H_s< 0)[0]),"error:method.H_s negative,%e \n" )
	method.check(len(np.where(method.C< 0)[0]),"error:method.C negative,%e \n" )
	method.P[0:method.num_RXN]=  method.H_s[0:method.num_RXN]*method.C[0:method.num_RXN]
	method.check(len(np.where(method.P<0)[0]), "error:method.P negative,%e \n")

def Sum_P():
	sum_prop= method.P.sum()
	return sum_prop

def CalcNCRxnNum(tau,scheme):
	for j in range(method.num_non_crit):
		index_2=method.I_nc.astype(int)[j]
		k= np.random.poisson(method.P[index_2]*tau)
		method.S[0:method.Type_spe,scheme] += np.multiply(method.V[index_2][:],k)

def CalcTauNonC (dose,time ):
	# Compute tau using (24): maximum leap time allowed by the Leap Condition for the non-critical reactions
	# tau = DBL_MAX-1
	tau =sys.float_info.max
	if(method.num_non_crit==0):
		return tau
	method.check(len(np.where(method.I_rs_init<0)[0]), "error:method.I_rs_init negative")
	method.check(len(np.where(method.I_rs_init>method.Type_spe)[0]), "error:method.I_rs_init more than method.Type_spe" )
	method.check(len(np.where(method.I_nc[range(method.num_non_crit)] < 0)[0]),"error:method.I_nc negative" )
	method.check(len(np.where(method.I_nc[range(method.num_non_crit)]  > method.num_RXN)[0]),"error:method.I_nc more than method.Type_spe")

	if time<6:
		index_1=method.I_rs_init
	else:
		index_1 = method.I_rs
	rxn_ind= method.I_nc[range(method.num_non_crit)]
	nom= np.maximum(method.eps_array[index_1]*method.S[index_1,dose],1)
	denom = np.fabs(np.einsum('ij,i->j',method.V[np.ix_(rxn_ind,index_1)],method.P[rxn_ind]))
	denom2 = np.einsum('ij,i->j',method.V[np.ix_(rxn_ind,index_1)]*method.V[np.ix_(rxn_ind,index_1)],method.P[rxn_ind])
	tau= min(np.minimum(nom/denom, nom*nom/denom2))
	return tau

def CalcTauC ():
	# Compute tau using (24): maximum leap time allowed by the Leap Condition for the non-critical reactions
	# tau = DBL_MAX-1
	if (method.num_crit == 0):
		tauc = sys.float_info.max
	else:
		method.check(np.where(method.I_c[range(method.num_crit)] <= 0)[0], "error:method.I_c negative")
		sum_prop_crit = method.P[method.I_c[range(method.num_crit)]].sum()
		random = method.rand_gen()
		#				printf("%10.8f  ", genrand64_real3());
		if (sum_prop_crit > 0):
			tauc = -np.log(random) / sum_prop_crit
		elif (sum_prop_crit == 0):
			tauc = sys.float_info.max - 1

		else:
			pass
	return tauc,sum_prop_crit

def RxnPartition(scheme):
	# Partition rxns into critical and non critical reactions
	# If all critical/non critical, I_c=[all rxn index]/I_nc=[np.inf]*method.num_RXN
	# 2. RxnPartition
	# * In state x at time method.t, evaluate all the aj(x)
	#         calcprop(x)
	# * Classify as critical any Rj for which aj(x) > 0 and which is within method.nc firings of exhausting any reactant.
	#     Classify all other Rj non-critical.
	#         for i in num_spec
	#             if method.S[i]<method.nc
	#                 find RXN_including_species_i : better way of doing?
	#                 method.I_c= include RXN, non duplicate
	#             method.I_nc= RxN-I_C: all the rest
	RXN_info = [method.RXN(0,1,[0],-1),method.RXN(1,2,[0,4],-1),method.RXN(2,2,[0,3],-1),method.RXN(3,1,[1],-1),method.RXN(4,1,[5],-1),method.RXN(5,1,[4],-1),method.RXN(6,1,[3],-1),method.RXN(7,1,[2],-1)]
	method.num_crit=0
	method.num_non_crit=0
	method.I_nc[:]=np.inf
	method.I_c[:]=np.inf
	for i in range(method.num_RXN):
		if len(np.where(method.S[RXN_info[i].SpecsinRXN[:],scheme] <method.nc)[0]):
			RXN_info[i].crit=1
			method.I_c[method.num_crit]=RXN_info[i].RXN_num
			#                    prf("critical rxns %d\n",I_c[num_crit]);
			method.num_crit +=1
		else:
			RXN_info[i].crit=0
			method.I_nc[method.num_non_crit]=RXN_info[i].RXN_num
	#                        prf("non critical rxns %d\n",I_nc[num_non_crit]);
			method.num_non_crit +=1

	method.check(method.num_crit>method.num_RXN,"error:num_crit pass num_RXN")
	method.check(method.num_non_crit> method.num_RXN,"error:num_non_crit pass num_RXN")
	method.check(method.num_non_crit+method.num_crit !=method.num_RXN,"error:total sum not equal to num_RXN")
#        prf("num_crit:%d , num_non_crit:%d", num_crit, num_non_crit);

def select_reaction(sum_prop,randN,RXN_list,num_RXN):
	reaction = -1
	sp = 0.0
	randN = randN  * sum_prop
	for i in range(num_RXN):
		index= RXN_list[i]
		sp += method.P[index]
		if(randN < sp):
			reaction = index
			break

	method.check(reaction<0,"error:negative reaction selected")
	method.check(reaction> method.num_RXN,"error:reaction index over number of rxns")
	return reaction

def init_F():
	F_dosing_schedule=open("F_dosing_schedule.txt","r")
	method.check(F_dosing_schedule == None,"file open error\n")
	method.F = np.loadtxt(F_dosing_schedule)
	F_dosing_schedule.close()

def init_V_C(t):
	if t<6:
		method.V = np.array([[-1,0,0,0,0,0],[-1,0,1,0,-1,0],[-1,0,1,-1,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,-1,0],[0,0,0,-1,0,0],[0,0,-1,0,0,0]])
		#	double C[num_RXN]= {1,(8.64E7)/(Volume*n_A),0,0,0.0233,0,0.0233};    (t<6)   # rate constant matrix c(/day) */
		method.C= np.array([1,(8.64E7)/(method.Volume*method.n_A),0,0,2,0.0233,0,0.0233])
 # rate constant matrix c(/day) */
	else:
		method.V = np.array( [[-1,0,0,0,0,0], [-1,0,1,0,-1,0], [-1,0,1,-1,0,0], [0,0,0,1,0,0], [0,0,0,0,0,1], [0,0,0,0,-1,0], [0,0,0,-1,0,0], [0,0,-1,0,0,0]])
		method.C = np.array([1,(8.64E7)/(method.Volume*method.n_A),(8.64E9)/(method.Volume*method.n_A), 2.5E7,2,0.02333, 0.02333, 0.02333])
	method.check(len(np.where(method.C < 0)[0]), "error:negative C matrix")
	method.check(len(np.where(method.V < -2)[0]), "error: V matrix smaller than -2")

def Adispe(j):
	decay_rate = {'S2': 0.0233, 'S3':0.0233, 'S4':0.0233}
	method.S[2,j ] -= method.S[2,j ] *0.5*decay_rate['S2']
	method.S[3,j] -= method.S[3,j ] *0.5*decay_rate['S3']
	method.S[4,j] -= method.S[4,j ]*0.5*decay_rate['S4']

def Gillespie_tau(f,scheme):
	t_old= method.t_range[0] # initialize time pt*/
	t_end= method.t_range[1] # end point for gillespie algorithm*/
	tau=0 # next rxn time*/
	tauc=0# next critical rxn time*/
	taunc=0 # next non critical rxn time*/
	r=0 #random # for reaction*/
	sum_prop_crit= 0.0	#sum of propencities*/
	reaction=0
	random=0#random # draw from uniform distribution for next rxn time*/
	count=0
	init_V_C(t_old)
	while(t_old<=t_end):
		sum_prop_crit=0
		tau=0
		tauc=0
		taunc=0
		r=0
		reaction=0
		method.num_gil_cycle+=1
		Calc_P(scheme)
		if Sum_P() ==0: #screen for early termination of the cycle
			method.flag_early_term =1
			return
		RxnPartition(scheme)
		taunc= CalcTauNonC(scheme,t_old)
		method.check(taunc<0,"error:CalcTauNonC negative")
		tauc,sum_prop_crit=CalcTauC()
		fmt= "%.5e\t"*method.Type_spe

		f.write("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n" %(t_old,method.S[0,scheme],method.S[1,scheme],method.S[2,scheme],method.S[3,scheme],method.S[4,scheme],method.S[5,scheme]))

		if (tauc < taunc):
			tau= tauc
			r= method.rand_gen()
			reaction = select_reaction(sum_prop_crit,r,method.I_c,method.num_crit)
		# */
			method.S[0:method.Type_spe,scheme] +=  method.V[reaction,0:method.Type_spe]
			method.check(np.where(method.S<0)[0],"error:method.S negative")
		else:
			tau= taunc
			CalcNCRxnNum(taunc, scheme)

		t_old= t_old+tau
		count = count+1
# # temporary place to save it as buffer
#         print ("%f\n", sum_prop) */