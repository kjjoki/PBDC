        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |         PBDC - THE PROXIMAL BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION          | | 
        !| |                                 (version 2)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                       by Kaisa Joki (last modified August 2015)                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |      NEW features :                                                              | |
        !| |                                                                                  | |
        !| |           * Possibility to use simple stepsize determination after               | |
        !| |             each 'main iteration'.                                               | |
        !| |                                                                                  | |
        !| |           * During each round of 'main iteration' utilizes OpenMP                | |
        !| |             to calculate subproblems in parallel. However if you DO NOT          | |
        !| |             WANT to use this feature then                                        | |
        !| |                1) in Makefile DELETE '-fopenmp' from line 5                      | |
        !| |                2) in pbdc.f95 COMMENT lines 412-415                              | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |     The software is free for academic teaching and research purposes but I       | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    Utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver.          |
        !|                                                                                      |
        !|    Utilizes PVMM by Ladislav Luksan as a norm minimization solver. This subroutine   | 
        !|    uses PQSUBS and MQSUBS by Ladislav Luksan. PVMM is a VARIABLE METRIC ALGORITHM    |
        !|    for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                          |
        !|                                                                                      |
        !|    The subroutine PVMM together with PQSUBS and MQSUBS is licensed by                |
        !|    the GNU Lesser General Public License (LGPL).                                     |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !|                                                                                      |
        !|                                                                                      |
        !|   Codes include:                                                                     |
        !|                                                                                      |
        !|   tpbdc.f95          - Main program for PBDC (this file)                             |
        !|   constants.f95      - Double precision (also some parameters)                       |
        !|   bundle1.f95        - Bundle of DC component f_1                                    |
        !|   bundle2.f95        - Bundle of DC component f_2                                    |
        !|   functions.f95      - User-specified DC components f_1 and f_2 together with        |
        !|                        subgradients of DC components. Contains also user-specified   |
        !|                        initial values for parameters                                 |
        !|   norm_min.f95       - Solver for the norm minimization problem                      |
        !|   fun.f95            - Defines objective funtion and gradient of the norm            |
        !|                        minimization problem                                          |
        !|   pbdc.f95           - PBDC method                                                   |
        !|                                                                                      |
        !|   plqdf1.f           - Quadratic solver by Ladislav Luksan                           |
        !|   pvmm.f             - Variable metric method by Ladislav Luksan                     |
        !|   mqsubs.f           - Basic modules for PVMM (by Ladislav Luksan)                   |
        !|   pqsubs.f           - Matrix modules for PVMM (by Ladislav Luksan)                  |
        !|                                                                                      |
        !|   Makefile           - Makefile                                                      |
        !|                                                                                      |
        !|   testproblems.f95   - Contains test problems used in [1]                            |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software MODIFY   tpbdc.f95   and   functions.f95   as needed           |
        !|                                                                                      |
        !|                                                                                      |
        !|   References:                                                                        |
        !|                                                                                      |
        !|   [1] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. MÃ¤kelÃ¤:               |
        !|       "New Proximal Bundle Method for Nonsmooth DC Optimization."                    |
        !|       TUCS Technical Report No. 1130, Turku Centre for Computer Science,             | 
        !|       Turku, 2015.                                                                   |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       


      PROGRAM tpbdc
      
         USE constants, ONLY : dp   ! double precision (i.e. accuracy)
         USE functions              ! INFORMATION from the USER
         USE bundle1                ! The BUNDLE of the DC component f_1
         USE bundle2                ! The BUNDLE of the DC component f_2
         USE pbdc                   ! PBDC method
         
        IMPLICIT NONE   
        
        ! 'user_n' is the number of variables in the problem (USER specifies this in MODULE functions.f95)
        REAL (KIND=dp), DIMENSION(user_n) :: x_0           ! The starting point in the PBDC method (specified by USER)
        REAL (KIND=dp), DIMENSION(user_n) :: x_solution    ! The solution obtained to the problem
        
        REAL(KIND=dp) :: f_solution         ! The objective function value at the solution 'x_solution'
        
        INTEGER, DIMENSION(6) :: counter    ! Contains the values of different counteres:
                                            !   counter(1) = iter_counter:         the number of 'main iterations' executed
                                            !   counter(2) = subprob_counter:      the number of subproblems solved
                                            !   counter(3) = f_counter:            the number of function values evaluated for DC component
                                            !   counter(4) = subgrad1_counter:     the number of subgradients calculated for f_1 
                                            !   counter(5) = subgrad2_counter:     the number of subgradients calculated for f_2 
                                            !   counter(6) = stop_cond_counter:    the number of times approximate stopping condition was tested during the algorithm 
 
    
        INTEGER :: iprint                   ! Variable that specifies print option (specified by USER): 
                                            !   iprint = 0 : print is suppressed
                                            !   iprint = 1 : basic print of final result 
                                            !   iprint = -1: basic print of final result (without the solution vector)
                                            !   iprint = 2 : extended print of final result 
                                            !   iprint = -2: extended print of final result (without the solution vector)
                                            !   iprint = 3 : basic print of intermediate results and extended print of final results
                                            !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                            !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                            !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                            !   iprint = 5 : prints each step of the PBDC algorithm (i.e. everything is printed step by step) (NOT recommended)
                                            !
                                            ! If 'iprint' <= -5 .OR. 'iprint' >= 6 then DEFAULT value 'iprint'=1 is used    
                                            
        
        INTEGER :: mrounds                  ! The maximum number of 'main iterations' (specified by USER).
                                            ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used

        INTEGER :: mit                      ! The maximum number of rounds during one 'main iteration' (specified by USER).
                                            ! If 'mit' <=0 then DEFAULT value 'mit'=1000 is used
        
        INTEGER :: termination              ! The reason for termination in PBDC:
                                            !   1 - the stopping condition is satisfied (i.e. criticality)
                                            !   2 - the approximate stopping condition is satisfied (i.e. eps-criticality)
                                            !   3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                            !   4 - the maximum number of 'main iterations' is executed 

        LOGICAL :: agg_used                 ! If .TRUE. then aggregation is used in PBDC (specified by USER).
        LOGICAL :: stepsize_used            ! If .TRUE. then simple stepsize determination is done in PBDC after each 'main iteration' (specified by USER).
        
        INTEGER :: i

          mrounds = 100000      ! maximum number of 'main iterations'
          mit = 100000          ! maximum number of rounds during one 'main iteration'
          iprint = -3           ! basic print of intermediate results and extended print of final results (without solution vector)
          
          agg_used = .TRUE.           ! Aggregation is used
          stepsize_used = .FALSE.     ! Simple stepsize determination is not used


!--------------------------------------------------------------------------------------
          
! The starting points for the test problems presented in MODULE testproblems.f95 and in the article [1]  
          
! Problem 6.1.
!        x_0 = (/ 2.0_dp, 2.0_dp /)
 
! Problem 6.2. (L1 version of Rosenbrock) 
!         x_0 = (/ -1.20_dp, 1.0_dp /)

! Problem 6.3.
!         x_0 = (/ 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp /)

! Problem 6.4  (for any user_n)
          DO i = 1, user_n/2
             x_0(i)  = 1.0_dp * i
          END DO
          DO i = (user_n/2 + 1), user_n
             x_0(i)  = -1.0_dp * i
          END DO             

! Problem 6.5.
!         x_0 = 0.0_dp

! Problem 6.6.
!        x_0 = (/ 10.0_dp, 1.0_dp /)

! Problem 6.7.
!         x_0 = (/ -2.0_dp, 1.0_dp /)
          
! Problem 6.8.        
!         x_0 = (/ 0.5_dp, 0.5_dp, 0.5_dp /)

! Problem 6.9.
!         x_0= (/ 4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp/)   

! Problem 6.10.
!          DO i = 1, user_n
!             x_0(i)  = 0.1_dp * i
!         END DO


!--------------------------------------------------------------------------------------


         WRITE(*,*) '------------------------------------------------------------------'
         WRITE(*,*) '** START ** START ** START ** START ** START ** START ** START **'  
         WRITE(*,*) '------------------------------------------------------------------'         
         WRITE(*,*) ' '

         CALL bundle_algorithm( x_0, x_solution, f_solution, mit, &
                            & mrounds, termination, counter, agg_used, stepsize_used, iprint )



         WRITE(*,*) ' '
         WRITE(*,*) '------------------------------------------------------------------'
         WRITE(*,*) '** END ** END ** END ** END ** END ** END ** END ** END ** END **'  
         WRITE(*,*) '------------------------------------------------------------------'
         WRITE(*,*) ' '      


      END PROGRAM tpbdc






















