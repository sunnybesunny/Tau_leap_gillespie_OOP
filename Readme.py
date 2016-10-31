# #  explicit tau leap method from Simulation Methods in Systems Biology (Daniel T. Gillespie) in Formal Methods
# for Computational Systems Biology
# Two parameters : n_c and eps are controlled


#  1.
#     *Choose  ε(0.05) ,nc(10)
#         eps= 0.05; ,nc= 10;, SpecinRxn= [RXN_num][speciesnum] /
#     *For each species Si, determine the functions εi(ε, xi) according to the rules (21).
#         HOR(i) = 1, set εi = ε /
#         HOR(i) = 2, set εi = ε/2,
#         Two Si reactant molecules set εi = ε/(2 + (xi − 1)^-1)
#
#         for i E species,
#         HOR(i)= k;
#         seteps(i)
#         eps= seteps(i) {
#         switch i ,
#         case 1
#             eps=eps;
#         case 2
#             eps = eps/2;
#         case 2.5
#             eps = eps/(2+(xi-1)^-1)
#         }
#
#  * Initialize t <-0 and x <- x0
#         t=0; x= x0;
#
# 2.
# * In state x at time t, evaluate all the aj(x)
#         calcprop(x)
# * Classify as critical any Rj for which aj(x) > 0 and which is within nc firings of exhausting any reactant.
#     Classify all other Rj non-critical.
#         for i in num_spec
#             if S[i]<nc
#                 find RXN_including_species_i : better way of doing?
#                 I_c= include RXN, non duplicate
#             I_nc= RxN-I_C: all the rest
#
# 3.
# *Compute τ' using (24) (except if there are no non-critical reactions take
# τ' = ∞). This is the maximum leap time allowed by the Leap Condition for the non-critical reactions.
#  τ' = min{max {εixi, 1}/|sumj∈Jncr (νij aj(x))|,{εixi, 1}^2/|sumj∈Jncr (νij^2 aj(x))|}
#
#     Tau= CalcTau (S,eps,I_c,C)
#         tau = Inf
#         for i in species
#             nom= max(eps*S(i),1)
#             denom= sum over j in I_nc (C[i][j]*P[j])
#             denom2= sum over j in I_nc (C[i][j]^2*P[j])
#             if tau> min(nom/denom, nom^2/denom2)
#                 tau = min(nom/denom, nom^2/denom2)
#             end
#         end
#
# 4.
#     *Compute τ''using SSA scheme (except if there are no critical reactions take τ'' = ∞)
#     This is the time to the next critical reactions.
#         Already implemented
#
# 5.
#     *Take τ = min(τ', τ''). Compute X(t+τ) using (28).
#         X(t + τ)= x + sumj∈Jncr Pj (aj(x)τ ) νj
#     tau = min(tau', tau'')
#     for j in I_nc
#         k= poisson(aj*tau)
#         S(j)= S+ k*C(j)
#  This gives the new state of the system if no critical reaction fired during the leap.
#
# 6.
#     *If τ'' <= τ', compute jc from (27) and replace X(t + τ) ← X(t + τ) + νjc
#     if (tau'' < tau')
#     select reaction
#     X(t + τ) ← X(t + τ) + νjc
#     end
#
# 7.
#     *Update x ← X(t + τ) and t ← t + τ.
#     X= X(t+tau)
#     t= T+tau
#     if (t< tend)
#         goback to step 2)
#     else stop
#     Return to Step 2, or else stop
# */
# Reactions
# r1: Ag(0) decay: C_Ag[0] -> 0 (k: 1/day)
# r2: IC(2) formation from IgM(4): C_Ag[0] +C_IgM[4]  -> C_IC[2] (k: 10^3/(M*S)*8.64*10^4(day/sec))
# r3: IC(2) formation from IgG(3):  C_Ag[0] +C_IgG[3]  -> C_IC[2] (k: 10^5/(M*S)*8.64*10^4(day/sec))
# r4: IgG(3) prod: C_PC[1] -> C_PC+C_IgG[3] (k:2.5E6 Pcells/day from model ): multiply by vol factor(lymph vol/serum vol= 0.5*10^-3) for comparision with serum measured values
# r5: Tfh(5) rep */ : T_fh[5] -> a T_fh[5] (k:? )
# r6: IgM[4] decay: IgM[4] -> 0 (k:0.0233/day)
# r7: IgG[3] decay: IgG[3] -> 0 (k:0.0233/day)
# r8: IC[2] decay: IC[2] -> 0  (k:0.0233/day)

# Variables
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