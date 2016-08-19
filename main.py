/*  explicit tau leap method
1.
    *Choose  ε(0.05) ,nc(10)
        eps= 0.05; ,nc= 10;, SpecinRxn= [RXN_num][speciesnum] /
    *For each species Si, determine the functions εi(ε, xi) according to the rules (21).
        HOR(i) = 1, set εi = ε /
        HOR(i) = 2, set εi = ε/2,
        Two Si reactant molecules set εi = ε/(2 + (xi − 1)^-1)

        for i E species,
        HOR(i)= k;
        seteps(i)
        eps= seteps(i) {
        switch i ,
        case 1
            eps=eps;
        case 2
            eps = eps/2;
        case 2.5
            eps = eps/(2+(xi-1)^-1)
        }

 * Initialize t <-0 and x <- x0
        t=0; x= x0;

2.
* In state x at time t, evaluate all the aj(x)
        calcprop(x)
* Classify as critical any Rj for which aj(x) > 0 and which is within nc firings of exhausting any reactant.
    Classify all other Rj non-critical.
        for i in num_spec
            if S[i]<nc
                find RXN_including_species_i : better way of doing?
                I_c= include RXN, non duplicate
            I_nc= RxN-I_C: all the rest

3.
*Compute τ' using (24) (except if there are no non-critical reactions take
τ' = ∞). This is the maximum leap time allowed by the Leap Condition for the non-critical reactions.
 τ' = min{max {εixi, 1}/|sumj∈Jncr (νij aj(x))|,{εixi, 1}^2/|sumj∈Jncr (νij^2 aj(x))|}

    Tau= CalcTau (S,eps,I_c,C)
        tau = Inf
        for i in species
            nom= max(eps*S(i),1)
            denom= sum over j in I_c (C[i][j]*P[j])
            denom2= sum over j in I_c (C[i][j]^2*P[j])
            if tau> min(nom/denom, nom^2/denom2)
                tau = min(nom/denom, nom^2/denom2)
            end
        end

4.
    *Compute τ''using SSA scheme (except if there are no critical reactions take τ'' = ∞)
    This is the time to the next critical reaction.
        Already implemented

5.
    *Take τ = min(τ', τ''). Compute X(t+τ) using (28).
        X(t + τ)= x + sumj∈Jncr Pj (aj(x)τ ) νj
    tau = min(tau', tau'')
    for j in I_nc
        k= poisson(aj*tau)
        S(j)= S+ k*C(j)
 This gives the new state of the system if no critical reaction fired during the leap.

6.
    *If τ'' <= τ', compute jc from (27) and replace X(t + τ) ← X(t + τ) + νjc
    if (tau'' < tau')
    select reaction
    X(t + τ) ← X(t + τ) + νjc
    end

7.
    *Update x ← X(t + τ) and t ← t + τ.
    X= X(t+tau)
    t= T+tau
    if (t< tend)
        goback to step 2)
    else stop
    Return to Step 2, or else stop
*/

#include <stdio.h>              /* File I/O */
#include <math.h>               /* for logf*/
#include "header_tau.h"
#include "dbg.h"
#include <limits.h>

FILE *State_EI;
FILE *State_ED;
FILE *State_Const;
FILE *State_PB;

//	double V[num_RXN][Type_spe]= {{-1,0,0,0,0,0},{-1,0,1,0,-1,0},{-1,0,1,-1,0,0},{0,0,0,1,0,0},{0,0,0,0,0,a-1}}; /*stoichiometry matrix */
//	double C[num_RXN]= {1,0,0,0,0};      /* rate constant matrix c */
//    double S[Type_spe][Type_dose]={{0}};  /* state_matrix S: update once every */

	double V[num_RXN][Type_spe]= {{-1,0,0,0,0,0},{-1,0,1,0,-1,0},{-1,0,1,-1,0,0},{0,0,0,1,0,0},{0,0,0,0,0,a-1}}; /*stoichiometry matrix */
	double C[num_RXN]= {1,(8.64E7)/Volume,0,3.11E7,0};      /* rate constant matrix c */
    double S[Type_spe][Type_dose]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1E7,1E7,1E7,1E7},{0,0,0,0}};  /* state_matrix S: update once every */
	double V_s[num_RXN][Type_spe]={{0.}};
	double C_s[num_RXN]= {0.};
	double H_s[num_RXN]= {0.}; /*Permutation matrix H  */
	/*Ag,B*,IC,IgG,IgM,Tfh */
	/*double S[Type_spe][Type_dose]={{0.}};   state_matrix S: update once every */
	double P[num_RXN]={0.}; /* propensity */
	double F[dose_per][Type_dose]={{0.}};
	double	t_range[2]={0.};
	double t_curr= 0.;
	char * File_name[]= {"State_EI.txt","State_ED.txt","State_Const.txt","State_PB.txt"};
    double eps_array[Type_spe]= {0};
    int HOR[Type_spe]= {2,1,0,2,2,1};
    int I_rs_init[Type_spe-1] ={0,1,3,4,5};
    int I_rs[Type_spe] ={0};/*reacting species t<6  */
    int I_nc[num_RXN]={0}; /*non critical reactions */
    int I_c[num_RXN]={0};  /*critical reactions*/
    int num_crit =0;
    int num_non_crit=0;
    int num_reactant_init=5;
    int num_reactant =5;


int main(void){

	int i,j;

//    /*r1: Ag(0) decay, r2: IC(2) formation from IgM(4), r3:IC(2) formation from IgG(3), r4: IgG(3) prod, r5: Tfh(5) rep */
//	double V[num_RXN][Type_spe]= {{-1,0,0,0,0,0},{-1,0,1,0,-1,0},{-1,0,1,-1,0,0},{0,0,0,1,0,0},{0,0,0,0,0,a-1}}; /*stoichiometry matrix */
//	double C[num_RXN]= {1,(8.64E7)/Volume,0,2.5E7,0};      /* rate constant matrix c */
//	double V_s[num_RXN][Type_spe]={{0.}};
//	double C_s[num_RXN]= {0.};
//	double H_s[num_RXN]= {0.}; /*Permutation matrix H  */
//	/*Ag,B*,IC,IgG,IgM,Tfh */
//	double S[Type_spe][Type_dose]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1E13,1E13,1E13,1E13},{0,0,0,0}};  /* state_matrix S: update once every */
//	/*double S[Type_spe][Type_dose]={{0.}};   state_matrix S: update once every */
//	double P[num_RXN]={0.}; /* propensity */
//	double F[dose_per][Type_dose]={{0.}};
//	double	t_range[1]={0};
//	double t_curr= 0.;
//	char * File_name[]= {"State_EI.txt","State_ED.txt","State_Const.txt","State_PB.txt"};
//    double eps_array[Type_spe]= {0};
//    int HOR[Type_spe]= {2,1,0,2,2,1};
//    int I_rs_init[Type_spe-1] ={0,1,3,4,5};
//    int num_reactant =5;
//
//    int I_rs[Type_spe] ={0};/*reacting species t>6  */
//    int I_nc[num_RXN]={0}; /*non critical reactions */
//    int I_c[num_RXN]={0};  /*critical reactions*/
//    int num_crit =0;
//    int num_non_crit=0;
 /*   int num_reactant_init=0 ; // # of reactants for t<6 */


	for (i = 0; i < Type_dose; i++)
	{
		outputFiles[i] = fopen(File_name[i], "a+");
//	    printf("%s\n",File_name[i]);
	}

    check(calc_eps()==0,"Error: calc_eps()");



/*Open up dosing schedule file, state files  */
    State_EI=fopen("State_EI.txt", "w");
    State_ED=fopen("State_ED.txt", "a+");
    State_Const=fopen("State_Const.txt", "a+");
    State_PB=fopen("State_PB.txt", "a+");


	/*Initialize F*/
	init_F();

	/* initialize V_s,C_s,H_s  */
    check(init_V_C()==0,"Error: init_V_C()");

//
//------------------ test set start-------------------------------------------------------------------- //
//
//for (j=0;j<Type_dose;j++) /* Operate Gillespie for 1 day time domain
////		{ */
	i=0;
	j=0;


	fprintf(outputFiles[j],"%s\t %s\t %s\t %s\t %s\t %s\t %s\n","time","Ag","B","IC","IgM","IgG","Tfh");

			/* for t<= 6 */
			while(i<dose_per)
				{
				    check(F[i][j]>=0,"F[%d][%d] negative ",i,j);
					S[0][j]+= F[i][j]; /* Update S[0] by respective dosing by reading F*/
                    check(Calc_H(j)==0,"error: Calc_H(%d)",j);
					t_range[0]=0;
					t_range[1]=t;
					check(t_range[1]>=0,"error:t_range[1] negative");
				//	printf("%f",t_range[1]);
                    int get= Gillespie_tau(outputFiles[j],j);
//                    printf("%d",get);
					check(get==0,"error:gillespie_tau") ; /* UPdate S*/
					i=i+1;
				}
/*			t_curr= dose_per; */

			/* for t> 6 */

//				while(t_curr<= t[1])
//					{
//						Calc_H(*S,*H_s); /* Update H_s from update S
//						Calc_C_s(*C_s);
//						t_range[0]=0; /*should be relevant time point
//						t_range[1]=1; /*should be relevant time point
//						Gillespie_Fast(*t_range,*V_s,*C_s,*H_s,*S,*P,File_name[j],num_RXN);
//						t_curr= t_range[1];
//					}
//		}

//------------------ test set end--------------------------------------------------------------------- //


//	original code
//
//		for (j=0;j<Type_dose;j++) /* Operate Gillespie for 1 day time domain
//			{
//			Calc_H(S,H_s,j);
//
//				/* for t<= 6
//				for(i=0;i<dose_per;i++)
//					{
//						S[0][j]+= F[i][j]; /* Update S[0] by respective dosing by reading F
//						Calc_H(S,H_s,j);
//						t_range[0]=i;
//						t_range[1]=i+1;
//						Gillespie_Fast(*t_range,V_s,C_s,H_s,S,P,File_name[j],num_RXN_init,j); /* UPdate S
//					}
//				t_curr= dose_per;
//
//				/* for t> 6
////				while(t_curr<= t[1])
////					{
////						Calc_H(*S,*H_s); /* Update H_s from update S
////						Calc_C_s(*C_s);
////						t_range[0]=0; /*should be relevant time point
////						t_range[1]=1; /*should be relevant time point
////						Gillespie_Fast(*t_range,*V_s,*C_s,*H_s,*S,*P,File_name[j],num_RXN);
////						t_curr= t_range[1];
////					}
//			} */
	   fclose(State_EI);
	   fclose(State_ED);
	   fclose(State_Const);
	   fclose(State_PB);

	   return (0);

	   error:
	       return (-1);
}
