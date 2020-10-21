      MODULE functions  
      
        USE constants, ONLY : dp   ! double precision (i.e. accuracy)    
        IMPLICIT NONE
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                 INFORMATION SUPPLIED BY THE USER:                        | | 
        !| |                                                                          | |
        !| |                                                                          | |
        !| |    * Different PARAMETERS:                                               | |
        !| |        - the Lipschitz constants:      'user_L1' and 'user_L2'           | |       
        !| |        - the number of variables:      'user_n'                          | |       
        !| |        - the size of bundle B_1:       'user_size_b1'                    | |       
        !| |        - the size of bundle B_2:       'user_size_b2'                    | |       
        !| |        - the descent parameter:        'user_m'                          | |       
        !| |        - the decrease parameter:       'user_r_dec'                      | |       
        !| |        - the increase parameter:       'user_r_inc'                      | |       
        !| |        - the criticality tolerance:    'user_crit_tol'                   | |       
        !| |        - the proximity measure:        'user_eps'                        | |       
        !| |                                                                          | |       
        !| |    * Computation of the value of the DC functions f_1 and f_2:           | |
        !| |        - f1(y)      the value of DC component f_1 at a point y           | |
        !| |        - f2(y)      the value of DC component f_2 at a point y           | |               
        !| |                                                                          | |               
        !| |                                                                          | |               
        !| |    * Computation of the subgradient of the DC functions f_1 and f_2:     | |
        !| |        - subgradient_f1(y)   the subgradient of DC component f_1 at y    | |
        !| |        - subgradient_f2(y)   the subgradient of DC component f_2 at y    | |       
        !| |                                                                          | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  INFORMATION ABOUT PARAMETERS:                            | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------       
        
    
        !****************** GLOBAL PARAMETERS *******************************************
      

        INTEGER, PARAMETER :: user_n = 200                       ! The number of variables in the problem
        
        INTEGER, PARAMETER :: user_size_b1 = MIN(user_n+5,500)   ! The biggest possible size of the bundle B_1 
                                                                 ! If user_size_b1 <= 0 then DEFAULT value MIN(user_n+5,1000) is used 
                                                                   
        INTEGER, PARAMETER :: user_size_b2 = 3                   ! The biggest possible size of the bundle B_2 
                                                                 ! If user_size_b2 <= 0 then DEFAULT value 3 is used  
              
        REAL(KIND=dp), PARAMETER :: user_m = 0.2_dp              ! The descent parameter  If user_m <= 0.0_dp .OR. user_m >= 1.0_dp 
                                                                 !                        then DEFAULT value 0.2_dp is used
                                                                
        REAL(KIND=dp), PARAMETER :: user_r_dec = 0.97_dp         ! The decrease parameter    
        
        !If user_r_dec <= 0.0_dp .OR. user_r_dec >= 1.0_dp then DEFAULT value is used.
        !                               
        !   DEFAULT value:                          
        !     If user_n < 10:           user_r_dec = 0.75_dp    
        !     If 10 <= user_n < 300:    user_r_dec = the first two decimals of n/(n+5)
        !     If user_n >= 300:         user_r_dec = 0.99_dp
        !
        !   Some examples of the DEFAULT value of the parameter 'user_r_dec':
        !     If user_n=10:     user_r_dec = 0.66_dp                          
        !     If user_n=20:     user_r_dec = 0.80_dp                         
        !     If user_n=50:     user_r_dec = 0.90_dp                         
        !     If user_n=100:    user_r_dec = 0.95_dp                         
        !     If user_n=150:    user_r_dec = 0.96_dp     
        !     If user_n=200:    user_r_dec = 0.97_dp                      
        !     If user_n=250:    user_r_dec = 0.98_dp    
        !
        
        REAL(KIND=dp), PARAMETER :: user_r_inc = 10000000.0_dp  ! The increase parameter     If user_r_inc <= 1.0_dp 
                                                                !                            then DEFAULT value 10000000.0_dp is used
        
    
        !****************** PARAMETRES NEEDED IN STOPPING CONDITIONS ********************
              
        ! ** The criticality tolerance **      
        REAL(KIND=dp), PARAMETER :: user_crit_tol = user_n*0.015_dp         ! If n <= 0.0_dp then DEFAULT value is used
        
        !   DEFAULT value:  
        !     If user_n < 150:           user_crit_tol = user_n*0.005_dp      
        !     If 150 <= user_n <= 200:   user_crit_tol = user_n*0.015_dp      
        !     If user_n > 200:           user_crit_tol = user_n*0.05_dp       
 
        
        ! ** The proximity measure **
        REAL(KIND=dp), PARAMETER :: user_eps = 0.1_dp           ! If user_eps <= 0.0_dp then DEFAULT value 0.1_dp is used
        

        !****************** LIPSCHITZ CONSTANTS *****************************************
        
                 ! set F_0 = { x\in\R^n | f(x) <= f(x_0)},  where x_0 is a starting point used in the algorithm 
                 ! 'user_eps' is the proximity measure selected by the user
             
        ! The (overestimated) Lipschitz constant of the DC component f_1 on the set {x\in\R^n | d(x, F_0) <= user_eps }      
        REAL(KIND=dp), PARAMETER :: user_L1 = 1000.0_dp 
        
        ! The (overestimated) Lipschitz constant of the DC component f_2 on the set {x\in\R^n | d(x, F_0) <= user_eps }     
        REAL(KIND=dp), PARAMETER :: user_L2 = 1000.0_dp       
         
        
        ! NOTICE: If user_L1 <= 0.0_dp or user_L2 <= 0.0_dp then DEFAULT value 1000.0_dp is used        
    
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |     To EXECUTE the bundle algorithm PBDC the USER needs to DETERMINE:    | | 
        !| |                                                                          | |       
        !| |        f1(y)               - value of DC component f_1 at a point y      | |
        !| |        f2(y)               - value of DC component f_2 at a point y      | |
        !| |                                                                          | |
        !| |        subgradient_f1(y)   - subgradient of DC component f_1 at y        | |
        !| |        subgradient_f2(y)   - subgradient of DC component f_2 at y        | |
        !| |                                                                          | |
        !| |                                                                          | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
    
        
        CONTAINS
        

        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                    TEST PROBLEM 6.4. in [1]                               | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------       

        
        !********************************************************************************
        !                                                                               |
        !              FUNCTION VALUES OF THE DC COMPONENTS f_1 AND f_2                 |
        !                                                                               |
        !********************************************************************************

           FUNCTION f1(y) RESULT(f)        
                !
                ! Calculates the function value of the DC component f_1 at a point 'y'
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'
                REAL(KIND=dp) :: apu                            ! help varaible
                INTEGER :: i                                    ! help variable
                
                apu = ABS(y(1))
                
                DO i = 2, user_n
                    IF (apu < ABS(y(i)) ) THEN 
                         apu = ABS(y(i))
                    END IF
                END DO
                
                f = user_n * apu
                
           END FUNCTION f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f2(y) RESULT(f)         
                !
                ! Calculates the function value of DC component f_2 at a point 'y'
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_2 at a point 'y'
                INTEGER :: i                                    ! help variable
                
                f = 0.0_dp
                
                DO i = 1, user_n
                   f =  f + ABS(y(i))                                   
                END DO              


           END FUNCTION f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           

        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS f_1 AND f_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_f1(y) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_1 is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'
                REAL(KIND=dp) :: apu                            ! help varaible
                INTEGER :: i, ind                                   ! help variable
                
                grad = 0.0_dp 
                
                apu = ABS(y(1))
                ind = 1
                
                DO i = 2, user_n
                    IF (apu < ABS(y(i))) THEN 
                         apu = ABS(y(i))
                         ind = i
                    END IF
                END DO
                
                IF (y(ind) <= 0.0_dp) THEN 
                    grad(ind) = - user_n
                ELSE 
                    grad(ind) = user_n
                END IF
                
           END FUNCTION subgradient_f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION subgradient_f2(y) RESULT(grad)              
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_2 is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'
                REAL(KIND=dp) :: apu                            ! help varaible
                INTEGER :: i                                    ! help variable
                
                DO i = 1, user_n
                    IF (y(i) <= 0.0_dp ) THEN 
                         grad(i) = -1.0_dp
                    ELSE
                         grad(i) = 1.0_dp
                    END IF
                END DO

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
           



           
      END MODULE functions     