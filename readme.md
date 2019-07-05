
## Simulating Antigen dosage and key immune cell kinetics using modified explicit tau-leap Gillespie simulation 

This is a gillespie simulation that describes the kinetics of vaccine dosage, immune complex formation, antibody production and decay. This simulation is one of three stochastic simulation modules that were developed to study the effect of vaccine kinetics on immune response. 

### Engineering vaccine for better antibody response 

Vaccination is undoubtedly one of the most effective preventive measures against infectious diseases. By exposing our body to pathogentic materials (Antigens) vaccination induces immune response and establishes memory of invasion that natural infeciton would do, but in a minute and maneagable scale. Recent development in bio compatible materials have promising and positive implications for vaccines, that through engineering the propoerty of vessel materials which deliver antigens, immune response - such as antibody quantity can be greatly enhanced. In fact, [Tam et al](https://www.ncbi.nlm.nih.gov/pubmed/27702895), explored vairous antigen delivery kinetics and found that spacing out the dosage into multiple smaller dosages results in greater antibody response. Moreover, they also found that if the kinetics in which these dosages follow exponentially increasing pattern,it can induce upto 20 times more antibody production. While Tam et al paper eludicated the mechanism through which confers more antibody production, namely through antigen availability, further study is warranted so as to understand how temporal vaccine kinetics affects affinity maturation through which affinity and quantity of antibodies are determined. Affinity maturation is a seqeunce of stochastic events which is in essence an accelerated darwinian evolution that leads to stronger affinity antibodies through mutation and seletion cycle. In order to understand how dosage kinetics alter affinity maturaiton and to predict the resulting immune response, I developed a explicit tau-leap gillespie simulation with a intent of combining this module with stochastic affinity maturation simulaiton, previously developed by [ShenShen Wang](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4357364/) 

### Kinetics of interest

The reacitons in which the simulation implement are followings. Note that parenthesis next to species(ex- Ag(0)) denotes the number at which the species are stored in the state array:

Ag[0]:   Antigen
PC[1]:   Plasma Bcell
IC[2]:   Immune complex
IgG[3]:  Immunoglobulin G
IgM[4]:  Immunoglobulin M
Tfh[5]:  Follicular helper T cell

   
    r1 Ag[0] decay:                 Ag[0] -> 0                       (k: 1/day)
    r2 IC[2] formation from IgM(4): Ag[0] +IgM[4] -> IC[2]           (k: 10^3/(M*S)*8.64*10^4(day/sec))
    r3 IC[2] formation from IgG(3): Ag[0] +0.0125* IgG[3] -> IC[2]   (k: 10^5/(M*S)*8.64*10^4(day/sec))
    r4 IgG[3] production:           PC[1] -> PC[1]+IgG[3]            (k: 2.5E6 Pcells/day from model)
    r5 IgG[3] decay:                IgG[3] -> 0                      (k: 0.0233/day)
    r6 IgM[4] decay:                IgM[4] --> 0                     (k: 0.02333/day)
    r7 Tfh[5] growth:               Tfh[5] --> 2Tfh[5]                     (Logistic growth, max rate: 0.048/day,
                                                                     carrying capacity: Germinal center B cell)
    r8 IC[2]                        IC[2]--> 0                       (k: 0.0233/day)  
    r9 Ag[0] supply                 Ag 				(Vaccine dosage: given one dosage each day for seven days)	

Reference: [Tam et al](https://www.ncbi.nlm.nih.gov/pubmed/27702895)   
    
### Tau leap gillespie simulation with partial deterministic approximation 

Tau leap gillespie simulation was chosen for its advantage of speed and less trade off with accuracy. The implementation was constructed based on the theoretical explanation and pseudo code in [Simulation Methods in Systems Biology (Daniel T. Gillespie) in Formal Methods for Computational Systems Biology](https://link.springer.com/chapter/10.1007/978-3-540-68894-5_5).  The initial explicit tau-leap simulation revealed the separation of disparate time scales: namely, r5,r6,r7,r8 from other reactions. However, though r5 and r6 occurs at much less frequency, due to the abundance of involved specieis - IgG and IgM- at different time of simulaiton, the approximation of these reactions were not appropriate, as it would distort the reaction propensity and thus lead to inaccurate kinetics simulation. On the other hand, reaction propensity for r7 and r8 was kept small throughout the simulation and therefore, chosen for deterministic approximation. The approximation further enhanced the speed of the simulaiton without sacrificing the accuracy. 
The speed and accuracy of the simulation can be manipulated by varying two tau-leap parameters: nc and epsilon
Additional model parameters are: number of simulations, dosage profiles and initial conditions of species.

