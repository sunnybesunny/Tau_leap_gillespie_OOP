import numpy as np

# Reactions
#r1: Ag(0) decay: C_Ag[0] -> 0 (k: 1/day)
#r2: IC(2) formation from IgM(4): C_Ag[0] +C_IgM[4]  -> C_IC[2] (k: 10^3/(M*S)*8.64*10^4(day/sec))
#r3: IC(2) formation from IgG(3):  C_Ag[0] +C_IgG[3]  -> C_IC[2] (k: 10^5/(M*S)*8.64*10^4(day/sec))
#r4: IgG(3) prod: C_PC[1] -> C_PC+C_IgG[3] (k:2.5E6 Pcells/day from model ): multiply by vol factor(lymph vol/serum vol= 0.5*10^-3) for comparision with serum measured values
#r5: Tfh(5) rep */ : T_fh[5] -> a T_fh[5] (k:? )
#r6: IgM[4] decay: IgM[4] -> 0 (k:0.0233/day)
#r7: IgG[3] decay: IgG[3] -> 0 (k:0.0233/day)
#r8: IC[2] decay: IC[2] -> 0  (k:0.0233/day)

#Variables
#	double V[num_RXN][Type_spe]= {{-1,0,0,0,0,0},{-1,0,1,0,-1,0},{-1,0,1,-1,0,0},{0,0,0,1,0,0},{0,0,0,0,0,a-1},{0,0,0,0,-1,0},{0,0,0,-1,0,0},{0,0,-1,0,0,0}}; #stoichiometry matrix */
#	double C[num_RXN]= {1,(8.64E7)/(Volume*n_A),0,0,0.0233,0,0.0233};    (t<6)   # rate constant matrix c(/day) */
#	double C[num_RXN]= {1,(8.64E7)/(Volume*n_A),(8.64E9)/(Volume*n_A),2.5E7,?(Need to know T cell biology),0.02333,0.02333,0.02333};    (t>6)   # rate constant matrix c(/day) */
#	double V_s[num_RXN][Type_spe]={{0.}}; ? what are these doing?
#	double C_s[num_RXN]= {0.};
#	double H_s[num_RXN]= {0.}; #Permutation matrix H  */
#	#Ag,B*,IC,IgG,IgM,Tfh */
#	double S[Type_spe][Type_dose]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1E13,1E13,1E13,1E13},{0,0,0,0}};  # state_matrix S: update once every, initial condition C_IgM(4) 10^-10M/S *(Vol*n_A) */
#	#double S[Type_spe][Type_dose]={{0.}};   state_matrix S: update once every */
#	double P[num_RXN]={0.}; # propensity */
#	double F[dose_per][Type_dose]={{0.}};
#	double	t_range[1]={0};
#	double t_curr= 0.;
#	char * File_name[]= {"State_EI.txt","State_ED.txt","State_Const.txt","State_PB.txt"};
#    double eps_array[Type_spe]= {0};
#    int HOR[Type_spe]= {2,1,1,2,2,1};
#    int I_rs_init[Type_spe-1] ={0,1,2,3,4,5};
#    int num_reactant =6;
#
#    int I_rs[Type_spe] ={0};#reacting species t>6  */
#    int I_nc[num_RXN]={0}; #non critical reactions */
#    int I_c[num_RXN]={0};  #critical reactions*/
#    int num_crit =0;
#    int num_non_crit=0;
 #   int num_reactant_init=0 ; # # of reactants for t<6 */

t = 0.05  # total t span ,should be parameterized after linking   */
dose_per = 1  # dosing period*/
num_RXN_init = 8  # # of rxn before t<6d*/
Type_spe = 6  # # of Type of species: C_Ag,C_B*,C_Ic,C_IgG,C_IgM,T_fH */
Type_dose = 4 #Type of dosing
num_RXN = 8 # # of rxn*/
a = 2  # Version 1 a=2*/
Volume = 0.000001  # volume of the lymph node, 1uL  */
n_A= 6.02e23
nc = 10  ## of critical rxns- labeled as critical, if within 10 rxns from depletion*/
eps = 0.05  # epsilon */

# Variables */
V = np.zeros((num_RXN, Type_spe))  # stoichiometry matrix */
C = np.zeros((num_RXN, ))  # rate constant matrix c */

V_s = np.zeros((num_RXN, Type_spe))
C_s = np.zeros((num_RXN, ))
H_s = np.zeros((num_RXN, ))  # Permutation matrix H  */
S = np.zeros((Type_spe, Type_dose))  # state_matrix S: */
P = np.zeros((num_RXN, ))  # propensity */
F = np.zeros((dose_per, Type_dose)) #dosage function
t_range = np.zeros((2, ))

HOR = np.zeros((Type_spe, ))  # highest order of reaction species i is involved  */
eps_array = np.zeros((Type_spe, ))
I_rs_init = np.zeros((Type_spe - 2, ),dtype=int)  # reacting species t<6  */
I_rs = np.zeros((Type_spe, ),dtype=int)  # reacting species t>6  */
I_nc = np.zeros((num_RXN, ),dtype=int)  # non critical reactions */
I_c = np.zeros((num_RXN, ),dtype=int)  # critical reactions*/
# / SpecInRxn[num_RXN][] */
num_crit = 0
num_non_crit = 0
num_reactant_init = 0  # # of reactants for t<6 */
num_reactant = Type_spe
num_gil_cycle=0
flag_early_term =0 # flag for early termination
# / # of reactants for t>6 */
class MyError(Exception):
	def __init__(self, value,expr):
		self.value = value
		self.expr = expr
	def __str__(self):
		return repr(self.expr)

class RXN(object):
	def __init__(self, RXN_num, num_specs, SpecsinRXN, crit):
		self.RXN_num = RXN_num
		self.num_specs = num_specs
		self.SpecsinRXN = SpecsinRXN  # np array  of SpecsinRXN[Type_spe]
		self.crit = crit

def check(bool,msg): # if False is passed, error is raised
	if bool:
		raise MyError(bool,msg)
	else:
		pass

def Poirand(lamb):
# need to document how this works
	k=0
	r=0
	L=np.exp(-lamb)
	p=1
	condition = True
	while condition:
		k+=1
		r= rand_gen()
		p *= r
		condition == (p > L)
		k -= 1
	return k

def rand_gen():
	floatinfo = np.finfo(float)
	epsilon = floatinfo.eps
	a = np.random.uniform(0+epsilon, 1-epsilon)
	return a
