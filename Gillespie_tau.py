#include <math.h>     /* This program simulates dynamics of affinity maturation in germinal centers */
#include <time.h>	 	        /* for seeding the rd # generator */
#include "header_tau.h"
#include "MT19937-64.h"

int Gillespie_tau(FILE * f,int scheme){

		double t_new=0; /* time marker for gillespie algorithm*/
		double t_old= (double)t_range[0]; /* initialize time pt*/
		double t_end= (double)t_range[1]; /* end point for gillespie algorithm*/
		double tau=0; /* next rxn time*/
		double tauc=0;/* next critical rxn time*/
		double taunc=0; /* next non critical rxn time*/
		double r=0; /*random # for reaction*/
		double sum_prop_crit= 0.0;	/*sum of propencities*/
		int reaction=0;
		double random=0;/*random # draw from uniform distribution for next rxn time*/
		int count=0;
		unsigned long long init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL}, length=4; /*initialize rd generator */
		init_by_array64(init, length);
		int i,j;

		while(t_old<t_end){
                /*printf("in loop"); */
				/*Calc Propensity */
				sum_prop_crit=tau= tauc=taunc=r=0; reaction=0;
            check(Calc_H(scheme)==0,"error: Calc_H(%d)",scheme);
////////////////////////////////////////////////////////////////////////////////////
////					Validation
//					for(i=0;i<num_RXN;i++){
//						printf("Gil_fast.H\n %e \n",H_s[i]);
//					}
/////////////////////////////////////////////////////////////////////////////////////

            check(Calc_P()==0,"error: Calc_P");
/*
2. RxnPartition
* In state x at time t, evaluate all the aj(x)
        calcprop(x)
* Classify as critical any Rj for which aj(x) > 0 and which is within nc firings of exhausting any reactant.
    Classify all other Rj non-critical.
        for i in num_spec
            if S[i]<nc
                find RXN_including_species_i : better way of doing?
                I_c= include RXN, non duplicate
            I_nc= RxN-I_C: all the rest
*/
            check(RxnPartition(scheme)==0,"error: RxnPartition(%d)",scheme);
/*
3.  CalcTauNonC
*Compute τ' using (24): maximum leap time allowed by the Leap Condition for the non-critical reactions
 τ' = min{max {εixi, 1}/|sumj∈Jncr (νij aj(x))|,{εixi, 1}^2/|sumj∈Jncr (νij^2 aj(x))|}
*/
            taunc= CalcTauNonC(scheme);
            check(taunc!=-1,"error:CalcTauNonC");
////						Validation
//					for(i=0;i<num_RXN;i++){
//						printf("Gil_fast.P\n %e \n",P_curr[i]);
//					}
/*
4.  CalcTauC
    *Compute τ''using SSA scheme (except if there are no critical reactions take τ'' = ∞)
    This is the time to the next critical reaction.
        Already implemented
*/
    int index_1=0;
            if(num_crit==0){ tauc = DBL_MAX-1;}
            else {
//                    printf("num_crit: %e",num_crit);
                    for(j=0; j<num_crit;j++){
                    index_1=I_c[j];
                    check(I_c[j]>=0,"error:I_c[%d] negative",j);
            //        printf("%d",index_1);
                    sum_prop_crit += P[index_1];
                        }
           // printf("sumpropcrit %e",sum_prop_crit);
//        Validation
//        printf("Gil_fast.sumprop\n %e\n",sum_prop);
//				Sample Tau: generate random number between (0,1) from MT twister
				random= genrand64_real3();
//				printf("%10.8f  ", genrand64_real3());

				if(sum_prop_crit > 0){
					tauc = -logf(random)/sum_prop_crit;
				}
				else if (sum_prop_crit == 0) {
                    tauc = DBL_MAX-1;
				}
                else { printf("sum_prop_crit<=0\n" );
                        return -1;
				}
        }


/*
5.  Update_S
    *Take τ = min(τ', τ''). Compute X(t+τ) using (28).
        X(t + τ)= x + sumj∈Jncr Pj (aj(x)τ ) νj
    tau = min(tau', tau'')
    for j in I_nc
        k= poisson(aj*tau)
        S(j)= S+ k*C(j)
 This gives the new state of the system if no critical reaction fired during the leap.
 */

                    fprintf(f,"%.17g\t %f\t %f\t %f\t %f\t %f\t %f\n",t_old,S[0][scheme],S[1][scheme],S[2][scheme],S[3][scheme],S[4][scheme],S[5][scheme]);

//
//                if (count % 1000000==0){
//                    fprintf(f,"%.17g\t %f\t %f\t %f\t %f\t %f\t %f\n",t_old,S[0][scheme],S[1][scheme],S[2][scheme],S[3][scheme],S[4][scheme],S[5][scheme]);
//                    }

//				fprintf(f,"%e\t %e\t %e\t %e\t %e\t %e\t %e\n",t_old,S[0][scheme],S[1][scheme],S[2][scheme],S[3][scheme],S[4][scheme],S[5][scheme]);
//                printf("tauc taunc%e %e \n",tauc,taunc);
                check(tauc >0,"error: tauc<0 %e", tauc);
                check(taunc >0,"erro: taunc<0 %e", taunc);
                tau= MIN(tauc,taunc); /* smaller of twos */
                printf("tau %e \n",tau);
                CalcNCRxnNum(tau,scheme);
/*
6.
    *If τ'' <= τ', compute jc from (27) and replace X(t + τ) ← X(t + τ) + νjc
    if (tau'' < tau')
    select reaction
    X(t + τ) ← X(t + τ) + νjc
    end
    */
            if (tau < taunc){
                    printf("tau<taunc");
                    r= genrand64_real3();
        /*						Validation
        				printf("Gil_fast.rand\n %e\n",r);
        */
                    reaction = select_reaction(sum_prop_crit,r,I_c,num_crit);
                    check(reaction != -1,"error:select_reaction");
        /*

        						Validation
        				printf("Gil_fast.reaction\n %e\n",reaction);
        */
                    for(i=0; i<Type_spe; i++){
                        S[i][scheme] +=  V_s[reaction][i];
                        check(S[i][scheme]>=0,"error:S[%d][%d] negative",i,scheme);
                    }
            }
/*
    7.
    *Update x ← X(t + τ) and t ← t + τ.
    X= X(t+tau)
    t= T+tau
    if (t< tend)
        goback to step 2)
    else stop
    Return to Step 2, or else stop
*/
        t_new= t_old+tau;
        t_old = t_new;
        count = count+1;
/* temporary place to save it as buffer
        printf("%f\n", sum_prop); */
	}

	return (0);

    error:
        return -1;

}