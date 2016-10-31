import numpy as np
import method
import Action as ac
def main():

	method.S=np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1E13,1E13,1E13,1E13],[3E2,3E2,3E2,3E2]])  # state_matrix S: update once every */
	method.HOR= np.array([2,1,1,2,2,1])
	method.I_rs_init=np.array([0,2,3,4]) #only reacting species
	method.num_reactant_init=6
	method.num_reactant =6
	t_curr= 0.

	#Open up dosing schedule file, state files  */
	State_EI=open("State_EI.txt", "w")
	State_ED=open("State_ED.txt", "a+")
	State_Const=open("State_Const.txt", "a+")
	State_PB=open("State_PB.txt", "a+")

	File_name= [State_EI,State_ED,State_Const,State_PB]
	#Initialize F*/
	ac.init_F()
	#------------------ test set start-------------------------------------------------------------------- #
	#
	#for (j=0;j<Type_dose;j++) # Operate Gillespie for 1 day time domain
	#		{ */
	i=0
	j=0

	f = File_name[j]
	f.write("time\tAg\tB\tIC\tIgM\tIgG\tTfh\n")
	ac.calc_eps()
		# for t<= 6 */
	while(i<min(method.dose_per,method.t)):
		method.S[0,j]+= method.F[j,i] # Update S[0] by respective dosing by reading F*/
		method.check(ac.Calc_H(j)==0,"error: Calc_H(%d)" %j)
		method.t_range[0]=i
		method.t_range[1]=min(i+1,method.t)
		method.check(method.t_range[1]<=0,"error:method.t_range[1] negative")
		#	printf("%f",t_range[1]);
		ac.Gillespie_tau(f,j)
		print method.num_gil_cycle
		#                    printf("%d",get);
		# UPdate S*/
		i+=1

	#			t_curr= dose_per; */

		# for t> 6 */

	#				while(t_curr<= t[1])
	#					{
	#						Calc_H(*S,*H_s); # Update H_s from update S
	#						Calc_C_s(*C_s);
	#						t_range[0]=0; #should be relevant time point
	#						t_range[1]=1; #should be relevant time point
	#						Gillespie_Fast(*t_range,*V_s,*C_s,*H_s,*S,*P,File_name[j],num_RXN);
	#						t_curr= t_range[1];
	#					}
	#		}

	#------------------ test set end--------------------------------------------------------------------- #


	#	original code
	#
	#		for (j=0;j<Type_dose;j++) # Operate Gillespie for 1 day time domain
	#			{
	#			Calc_H(S,H_s,j);
	#
	#				# for t<= 6
	#				for(i=0;i<dose_per;i++)
	#					{
	#						S[0][j]+= F[i][j]; # Update S[0] by respective dosing by reading F
	#						Calc_H(S,H_s,j);
	#						t_range[0]=i;
	#						t_range[1]=i+1;
	#						Gillespie_Fast(*t_range,V_s,C_s,H_s,S,P,File_name[j],num_RXN_init,j); # UPdate S
	#					}
	#				t_curr= dose_per;
	#
	#				# for t> 6
	#				while(t_curr<= t[1])
	#					{
	#						Calc_H(*S,*H_s); # Update H_s from update S
	#						Calc_C_s(*C_s);
	#						t_range[0]=0; #should be relevant time point
	#						t_range[1]=1; #should be relevant time point
	#						Gillespie_Fast(*t_range,*V_s,*C_s,*H_s,*S,*P,File_name[j],num_RXN);
	#						t_curr= t_range[1];
	#					}
	#			} */
	State_EI.close()
	State_ED.close()
	State_Const.close()
	State_PB.close()


if __name__ == "__main__":
    import cProfile
    cProfile.run('main()')


