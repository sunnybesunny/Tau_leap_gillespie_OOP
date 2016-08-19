#include "../tau_leap/header_tau.h"
#include "../tau_leap/dbg.h"

int calc_eps(){
    int i;

    for (i=0; i< Type_spe; i++ ){

        check(HOR[i]>=0,"error: HOR_array_negative");
       switch ( HOR[i] ) {
            case 1:
              eps_array[i]= eps;// Code
              break;
            case 2:
              eps_array[i]= eps/2;// Code
              break;
            default:
                eps_array[i]= 0; // when it is non reacting species
              break;
            }
        check(eps_array[i]>=0,"error: eps_array_negative");
    }
    return(0);

    error:
        return -1;

}

#include "header_tau.h"
#include "dbg.h"
#include <math.h>

int Calc_H(int j)
	{
        int i;
        for (i=0; i<Type_spe; i++){
            check(S[i][j]>=0,"error:S[%d][%d] negative",i,j);
            check(S[i][j]<DBL_MAX, "error:S[%d][%d] overflow",i,j);
            check(S[4][j]*S[0][j]<DBL_MAX, "error:H_s[1] overflow");
            check(S[3][j]*S[0][j]<DBL_MAX, "error:H_s[2] overflow");
        }

		H_s[0] = S[0][j];
		H_s[1] =S[4][j]*S[0][j] ;
		H_s[2] =S[3][j]*S[0][j] ;
		H_s[3] =S[1][j] ;
		H_s[4] =S[5][j] ;

     //printf("%f\n",H_curr[2]);
/////////////////////////////////////////////////////////////////////
//		int i;
//		for(i=0;i<num_RXN;i++){
//			printf("Calc_H.H\n %e \n",H_s[i]);
//		}
/////////////////////////////////////////////////////////////////////
    return(0);

    error:
        return -1;
	}


#include "header_tau.h"
#include "dbg.h"

int Calc_P(){

	int i;
	for(i=0;i<num_RXN;i++)
		{
        check(H_s[i]>=0,"error:H_s[%d] negative,%e \n",i,H_s[i]);
        check(C_s[i]>=0,"error:C_s[%d] negative,%e \n",i,C_s[i]);

		P[i]=  H_s[i]*C_s[i];
	//	printf("%e \n", P[i]);
        check(P[i]<DBL_MAX,"error:P[%d] overflow,%e",i,P[i]);

//        printf("Calc_P %f\n",P[i]);
		}
    return(0);

	error:
	    return -1;

	}
#include "header_tau.h"
#include "dbg.h"

double Sum_P()
{
	int i;
	double sum_prop=0;

	for (i=0;i<num_RXN;i++)
		{
			sum_prop += P[i];

		}
	return sum_prop;

}

#include <stdio.h>              /* This program simulates dynamics of affinity maturation in germinal centers */
#include <stdlib.h>
#include <math.h>
#include <time.h>	 	        /* for seeding the rd # generator */
#include <limits.h>
#include "header_tau.h"
#include "MT19937-64.h"
#include <float.h>
#include "dbg.h"

void CalcNCRxnNum(double tau,int scheme){

    int index_2=0;
    int i,j,k=0;

    for(j=0; j<num_non_crit;j++){

        index_2=I_nc[j];
//                                printf("lambda %e\n",P[index_2]*tau);
        k= Poirand(P[index_2]*tau);
//                                printf("k %e\n",k);
            for(i=0; i<Type_spe; i++){
                S[i][scheme] += k* V_s[index_2][i];
        //                                    printf("%e\n", S[i][scheme]);
            }
    }


}

#include <math.h>               /* for logf*/
#include <limits.h>
#include "header_tau.h"
#include <float.h>
#include "dbg.h"
/*
 Taunc= CalcTauNonC (S,eps,I_c,C)
        tau = Inf
        for i in species
            nom= max(eps*S(i),1)
            denom= sum over j in I_c (C[i][j]*P[j])
            denom2= sum over j in I_c (C[i][j]^2*P[j])
            if tau> min(nom/denom, nom^2/denom2)
                tau = min(nom/denom, nom^2/denom2)
            end
        end
*/
double CalcTauNonC (int dose){
    int i,j,index_1,index_2; index_1=index_2=0;
    double nom,denom,denom2; nom=denom=denom2=0;
    double tau = DBL_MAX-1;

    if(num_non_crit==0){
            return tau;
    }

        for(i=0;i<num_reactant_init;i++){
                nom=denom=denom2=0;
                index_1=index_2=0;

                check(I_rs_init[i]>=0, "error:I_rs_init[%d] negative",i)
                check(I_rs_init[i]<= Type_spe, "error:I_rs_init[%d] more than Type_spe",i)
                index_1=I_rs_init[i];
//            printf("%e %e ",eps_array[index_1],S[index_1][dose]);
//            printf("%e ",eps_array[index_1]*S[index_1][dose]);
            nom= MAX(eps_array[index_1]*S[index_1][dose],1);
//            printf(" nom %e \n", nom);

            for(j=0; j<num_non_crit; j++ ){
                check(I_nc[j]>=0, "error:I_nc[%d] negative",i)
                check(I_nc[j]<= Type_spe, "error:I_nc[%d] more than Type_spe",j)
                index_2=I_nc[j];
                denom+= fabs((V_s[index_2][index_1]*P[index_2]));
                denom2+= (V_s[index_2][index_1]*V_s[index_2][index_1]*P[index_2]);
            }
//            printf(" denom %e \n", denom);
//            printf(" denom2 %e \n", denom2);
        if (tau> MIN(nom/denom, nom*nom/denom2)){
//            printf(" tau %e \n", tau);
//            printf(" nom/denom %e \n", nom/denom);
//            printf(" nom*nom/denom2 %e \n", nom*nom/denom2);
//            printf(" new time %e \n", MIN(nom/denom, nom*nom/denom2));
            tau = MIN(nom/denom, nom*nom/denom2);


        }

        }
        return tau;

    error:
        return -1;
 }

#include "dbg.h"

int select_reaction(double sum_propencity, double randN,int RXN_list[num_RXN],int RXN_list_length){
	int reaction = -1;
	double sp = 0.0;
	int i;
	int index;
	randN = randN * sum_propencity;
	for(i=0; i<RXN_list_length; i++){
        index= RXN_list[i];
		sp += P[index];
		if(randN < sp){
			reaction = index;
			break;
		}
	}
    check(reaction>=0,"error:negative reaction selected");
    check(reaction<=num_RXN,"error:reaction index over number of rxns");
	return reaction;

	error:
	    return -1;
}