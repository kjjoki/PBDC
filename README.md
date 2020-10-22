# PBDC
The proximal bundle method for nonsmooth DC optimization

PBDC is a proximal double bundle solver (Fortran 95) for nonsmooth DC programming (difference of two convex functions) by Kaisa Joki. 

The software utilizes code PVMM by Prof. Ladislav Luksan that is licensed by the GNU Lesser General Public License (LGPL). In addition, code PLQDF1 by Prof. Ladislav Luksan is used to solve quadratic direction finding problem.

The software utilizes OpenMP at each round of 'main iteration' to calculate subproblems in parallel. To turn down OpenMP, see instructions in tpbdc.f95. In addition, there is a possibility to use simple stepsize determination after each 'main iteration'. 

The software is free for academic teaching and research purposes but I ask you to refer the reference given below if you use it. To use the software modify tpbdc.f95 and functions.f95 as needed. If you have any questions conserning the software, please contact directly the author Kaisa Joki (email: kjjoki@utu.fi).


# Codes include:                                                                     
         
   tpbdc.f95          - Main program for PBDC 
   
   constants.f95      - Double precision (also some parameters) 
   
   bundle1.f95        - Bundle of DC component f_1
   
   bundle2.f95        - Bundle of DC component f_2                                    
        
   functions.f95      - User-specified DC components f_1 and f_2 together with subgradients of DC components. Contains also user-specified initial values for parameters        
   
   norm_min.f95 	    - Solver for the norm minimization problem.
   
   fun.f95 	          - Defines objective funtion and gradient of the norm minimization problem.
   
   pbdc.f95 	        - PBDC method.
	
   plqdf1.f 	        - Quadratic solver by L. Luksan.
  
   pvmm.f               - Variable metric method by L. Luksan.
   
   mqsubs.f 	        - Basic modules for PVMM by L. Luksan.
   
   pqsubs.f 	        - Matrix modules for PVMM by L. Luksan.
	
   Makefile 	        - Makefile.
	
   testproblems.f95 	- Test problems used in [1].
   
   
# References:                                                                        
                                                                                              
[1] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. Mäkelä: "A proximal bundle method for nonsmooth DC optimization utilizing nonconvex cutting planes". J. Glob. Optim. 68 (2017), pp. 501-535, https://doi.org/10.1007/s10898-016-0488-3                                       
  
   
