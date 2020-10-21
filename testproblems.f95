   
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |             TEST PROBLEMS FOR THE PROXIMAL BUNDLE METHOD PBDC                    | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|  References:                                                                         |
        !|                                                                                      |
        !|  [1] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. MÃ¤kelÃ¤:                |
        !|      "New Proximal Bundle Method for Nonsmooth DC Optimization."                     |
        !|      TUCS Technical Report No 1130, Turku Centre for Computer Science, Turku, 2015.  | 
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
   
   
   
   
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                   THE PROBLEM 6.1. in [1]                                 | |
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
                REAL(KIND=dp) :: a1, a2, a3, a4, a5, a6          
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'
                
                a1 = y(1)**4 + y(2)**2
                a2 = (2.0_dp - y(1))**2 + (2.0_dp - y(2))**2 
                a3 = 2.0_dp * EXP(-y(1)+y(2))

                a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                a5 = a5 + 4.0_dp
                a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp                
                
                f = MAX(a1,a2,a3) + a4 + a5 + a6
    
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
                REAL(KIND=dp) :: a4, a5, a6                     
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_2 at a point 'y'
                
                a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                a5 = a5 + 4.0_dp
                a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp
                
                f = MAX( (a4 + a5) , (a4 + a6) , (a5 + a6) )
                
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
                REAL(KIND=dp) :: a1, a2, a3 
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'

                a1 = y(1)**4 + y(2)**2
                a2 = (2.0_dp - y(1))**2 + (2.0_dp - y(2))**2 
                a3 = 2.0_dp * EXP(-y(1)+y(2))               
                
                IF (a1 >= a2) THEN
                    IF (a1 >= a3) THEN
                       grad(1) = 4.0_dp * y(1)**3 
                       grad(2) = 2.0_dp * y(2)                  
                    ELSE
                       grad(1) = -2.0_dp * EXP(-y(1)+y(2))
                       grad(2) = 2.0_dp * EXP(-y(1)+y(2))                   
                    END IF              
                ELSE
                    IF (a2 >= a3) THEN
                       grad(1) = 2.0_dp * y(1) - 4.0_dp
                       grad(2) = 2.0_dp * y(2) - 4.0_dp                 
                    ELSE
                       grad(1) = -2.0_dp * EXP(-y(1)+y(2))
                       grad(2) = 2.0_dp * EXP(-y(1)+y(2))                   
                    END IF                  
                END IF
                
                grad(1) = grad(1) + 2.0_dp * y(1) - 2.0_dp
                grad(1) = grad(1) + 4.0_dp * y(1) - 5.0_dp
                grad(1) = grad(1) + 2.0_dp * y(1)
  
                grad(2) = grad(2) + 2.0_dp * y(2) - 4.0_dp
                grad(2) = grad(2) + 2.0_dp * y(2) - 2.0_dp
                grad(2) = grad(2) + 4.0_dp * y(2) - 4.0_dp
                
                
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
                REAL(KIND=dp) :: a4, a5, a6 
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'

                grad = 0.0_dp
                
                a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                a5 = a5 + 4.0_dp
                a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp                
                
                IF ( (a4 + a5) >= (a4 + a6) ) THEN
                    IF ( (a4 + a5) >= (a5 + a6)) THEN
                       grad(1) = 6.0_dp * y(1) - 7.0_dp
                       grad(2) = 4.0_dp * y(2) - 6.0_dp
                    ELSE
                       grad(1) = 6.0_dp * y(1) - 5.0_dp
                       grad(2) = 6.0_dp * y(2) - 6.0_dp
                    END IF              
                ELSE
                    IF ( (a4 + a6) >= (a5 + a6) ) THEN
                       grad(1) = 4.0_dp * y(1) - 2.0_dp
                       grad(2) = 6.0_dp * y(2) - 8.0_dp
                    ELSE
                       grad(1) = 6.0_dp * y(1) - 5.0_dp
                       grad(2) = 6.0_dp * y(2) - 6.0_dp 
                    END IF                  
                END IF

           END FUNCTION subgradient_f2
                
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           

        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  L1 VERSION OF ROSENBROCK FUNCTION                        | |
        ! |                   THE PROBLEM 6.2. in [1]                                 | |
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
                
                IF ((ABS(y(1)) - y(2) ) > 0.0_dp) THEN
                     f = ABS(y(1) - 1.0_dp) + 200.0_dp * (ABS(y(1)) - y(2) )
                ELSE
                     f = ABS(y(1) - 1.0_dp)
                END IF
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
                
                f = 100.0_dp * ( ABS(y(1)) - y(2) )

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
    
                IF ( y(1) <= 1.0_dp ) THEN
                     grad(1) = -1.0_dp
                ELSE
                     grad(1) = 1.0_dp
                END IF
                
                IF( ( ABS(y(1))-y(2) ) > 0.0_dp  ) THEN
                     IF(y(1) <= 0.0_dp) THEN 
                          grad(1) = grad(1) - 200.0_dp
                     ELSE
                          grad(1) = grad(1) + 200.0_dp  
                     END IF                  
                     grad(2) = -200.0_dp
                ELSE         
                     grad(1) = grad(1) 
                     grad(2) = 0.0_dp
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
                
                IF(y(1) <= 0) THEN
                    grad(1) = -100.0_dp 
                ELSE
                    grad(1) = 100.0_dp
                END IF
                
                grad(2) = -100.0_dp

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
       

        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                        TEST PROBLEM 6.3. in [1]                           | |
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
                ! NOTICE: The dimension of 'y' has to be 'user_n' ('user_n' needs to be 4).
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'
                
                f = ABS(y(1) - 1.0_dp)
                f = f + ABS(y(3) - 1.0_dp )                 
                f = f + 4.95_dp * ( ABS(y(2) + y(4) - 2.0_dp ) )
                f = f + 10.1_dp * ( ABS(y(2)-1.0_dp) + ABS(y(4)-1.0_dp) )       
                
                IF ((ABS(y(1)) - y(2) ) > 0.0_dp) THEN
                     f = f + 200.0_dp * (ABS(y(1)) - y(2) )
                END IF 
                
                IF ((ABS(y(3)) - y(4) ) > 0.0_dp) THEN
                     f = f + 180.0_dp * (ABS(y(3)) - y(4) )
                END IF          



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
                
                f = 4.95_dp * ( ABS(y(2) - y(4)) )
                f = f + 90.0_dp * ( ABS(y(3)) - y(4) )
                f = f + 100.0_dp * ( ABS(y(1)) - y(2) ) 
            
                

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
                
                grad(1) = 0.0_dp
                grad(2) = 0.0_dp
                grad(3) = 0.0_dp
                grad(4) = 0.0_dp
                
                IF ( y(1) <= 1.0_dp ) THEN
                     grad(1) = -1.0_dp
                ELSE
                     grad(1) = 1.0_dp
                END IF
                
                IF( ( ABS(y(1))-y(2) ) > 0.0_dp  ) THEN
                     IF(y(1) <= 0.0_dp) THEN 
                          grad(1) = grad(1) - 200.0_dp
                     ELSE
                          grad(1) = grad(1) + 200.0_dp  
                     END IF                  
                     grad(2) = -200.0_dp
                ELSE         
                     grad(1) = grad(1) 
                     grad(2) = 0.0_dp
                END IF

                IF( ( ABS(y(3))-y(4) ) > 0.0_dp  ) THEN
                     IF(y(3) <= 0.0_dp) THEN 
                          grad(3) =  - 180.0_dp
                     ELSE
                          grad(3) =  180.0_dp  
                     END IF                  
                     grad(4) = -180.0_dp
                ELSE         
                     grad(3) = 0.0_dp 
                     grad(4) = 0.0_dp
                END IF
                
                IF ( y(3) <= 1.0_dp ) THEN
                     grad(3) = grad(3) - 1.0_dp
                ELSE
                     grad(3) = grad(3) + 1.0_dp
                END IF          
                
                IF ( y(2) <= 1.0_dp ) THEN
                     grad(2) = grad(2) - 10.1_dp  
                ELSE
                     grad(2) = grad(2) + 10.1_dp
                END IF          
                
                IF ( y(4) <= 1.0_dp ) THEN
                     grad(4) = grad(4) - 10.1_dp  
                ELSE
                     grad(4) = grad(4) + 10.1_dp
                END IF  
                
                IF ( (y(2) + y(4) - 2.0_dp) <= 0.0_dp ) THEN
                     grad(2) = grad(2) - 4.95_dp
                     grad(4) = grad(4) - 4.95_dp  
                ELSE
                     grad(2) = grad(2) + 4.95_dp
                     grad(4) = grad(4) + 4.95_dp
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
                
                IF(y(1) <= 0.0_dp) THEN
                    grad(1) = -100.0_dp 
                ELSE
                    grad(1) = 100.0_dp
                END IF
                
                grad(2) = -100.0_dp
                
                IF(y(3) <= 0.0_dp) THEN
                    grad(3) = -90.0_dp 
                ELSE
                    grad(3) = 90.0_dp
                END IF
                
                grad(4) = -90.0_dp              
                
                IF( (y(2) - y(4) ) <= 0.0_dp) THEN
                    grad(2) = grad(2) - 4.95_dp 
                    grad(4) = grad(4) + 4.95_dp
                ELSE
                    grad(2) = grad(2) + 4.95_dp 
                    grad(4) = grad(4) - 4.95_dp                
                END IF
                
           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
       

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

           
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.5. in [1]                       | |
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
                INTEGER :: j                                    ! help variable

                apu = ABS(summa(y,1)) 

                DO j = 2, 20
                    IF (apu <= ABS(summa(y,j)) ) THEN 
                         apu = ABS(summa(y,j))
                    END IF
                END DO
                
                f = 20.0_dp * apu
                
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
                INTEGER :: j                                ! help variable
                
                f = 0.0_dp
                
                DO j = 1, 20
                    f = f + ABS( summa(y,j) )
                END DO              


           END FUNCTION f2

           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION summa(y, j) RESULT(f)           
                !
                ! Calculates the sum used in f_1 and f_2 for parameter t_j
                !
                ! NOTICE: The dimension of 'y' has to be 'n'. 'j' needs to be an integer from interval [1,20]
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) ::  j    ! used to determines the parameter t_j
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f, t, apu                          ! the function value of the sum used in f_1 and parameter t_j
                INTEGER :: i                                    ! help variable
                
                f = 0.0_dp
                apu = 1.0_dp/ user_n
                
                DO i = 1, user_n 
                    t = (0.05_dp * j )**(i-1)
                    f = f + (y(i)-apu)*t            
                END DO 
            
           END FUNCTION summa
           
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
                INTEGER :: i, j, ind                            ! help variable
                
                grad = 0.0_dp 

                apu = ABS(summa(y,1)) 
                ind = 1

                DO j = 2, 20
                    IF (apu <= ABS(summa(y,j)) ) THEN 
                         apu = ABS(summa(y,j))
                         ind = j
                    END IF
                END DO

                DO i = 1, user_n
                     grad(i) = (0.05_dp * ind)**(i-1)
                END DO              
                
                IF ( summa(y,ind) <= 0.0_dp ) THEN 
                         grad = -20.0_dp * grad
                ELSE
                         grad = 20.0_dp * grad
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
                INTEGER :: i, j                                 ! help variable
                
                grad = 0.0_dp 

                DO j = 1, 20
                    IF (summa(y,j) <= 0.0_dp) THEN 
                        DO i = 1, user_n
                            grad(i) = grad(i) - (0.05_dp * j)**(i-1)
                        END DO                  
                    ELSE 
                        DO i = 1, user_n
                            grad(i) = grad(i) + (0.05_dp * j)**(i-1)
                        END DO                   
                    END IF
                END DO

           END FUNCTION subgradient_f2
            
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -              
        
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.6. in [1]                       | |
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

                f = y(2) + 0.1_dp * (y(1)**2 + Y(2)**2)  

                IF (-y(2) <= 0.0_dp ) THEN 
                         f = f
                    ELSE
                         f = f + 10.0_dp * (-y(2))
                END IF

                
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
                
                IF (y(1) <= 0.0_dp ) THEN 
                         f = f - y(1)
                    ELSE
                         f = f + y(1) 
                END IF
                
                IF (y(2) <= 0.0_dp ) THEN 
                         f = f - y(2)
                    ELSE
                         f = f + y(2)
                END IF              


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
                
                grad(1) = 0.2_dp * y(1)
                
                grad(2) = 1.0_dp + 0.2_dp * y(2) 

                IF ( -y(2) <= 0.0_dp ) THEN 
                         grad(2) = grad(2)
                    ELSE
                         grad(2) = grad(2) - 10.0_dp
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


        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.7. in [1]                       | |
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
                REAL(KIND=dp) :: f, a1, a2, a3, a4              ! the function value of the DC component f_1 at a point 'y'
                REAL(KIND=dp) :: apu                            ! help varaible
                INTEGER :: i                                    ! help variable

                f = ABS(y(1)-1.0_dp) + 200.0_dp * MAX(0.0_dp, ABS(y(1))-y(2) )                  
                
                a1 = y(1)**2 + y(2)**2 + ABS(y(2))
                a2 = y(1) + y(1)**2 + y(2)**2 + ABS(y(2)) - 0.5_dp
                a3 = ABS( y(1) - y(2) ) + ABS(y(2)) - 1.0_dp
                a4 = y(1) + y(1)**2 + y(2)**2
                
                f = f + 10.0_dp * MAX(a1,a2,a3,a4)
                
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
                
                f = 10.0_dp * ( y(1)**2 + y(2)**2 + ABS(y(2)) )
                f = f + 100.0_dp * ( ABS(y(1)) - y(2) )

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
                REAL(KIND=dp), DIMENSION(4) :: a        
                REAL(KIND=dp) :: apu                            ! help varaible
                INTEGER :: i, ind                               ! help variable
                
                grad = 0.0_dp 
                
                IF ( (y(1)-1.0_dp) <= 0.0_dp ) THEN 
                     grad(1) = -1.0_dp
                ELSE
                     grad(1) = 1.0_dp
                END IF 

                a(1) = y(1)**2 + y(2)**2 + ABS(y(2))
                a(2) = y(1) + y(1)**2 + y(2)**2 + ABS(y(2)) - 0.5_dp
                a(3) = ABS( y(1) - y(2) ) + ABS(y(2)) - 1.0_dp
                a(4) = y(1) + y(1)**2 + y(2)**2             
                
                ind = 1
                DO i = 2, 4
                   IF ( a(ind) < a(i) ) THEN 
                       ind = i
                   END IF 
                END DO
                
                IF (ind == 1) THEN
                    grad(1) = grad(1) + 20.0_dp * y(1)
                    grad(2) = grad(2) + 20.0_dp * y(2)
                    IF (y(2) <= 0.0_dp) THEN
                        grad(2) = grad(2) - 10.0_dp
                    ELSE
                        grad(2) = grad(2) + 10.0_dp
                    END IF
                END IF
                
                IF (ind == 2) THEN
                    grad(1) = grad(1) + 10.0_dp + 20.0_dp * y(1)
                    grad(2) = grad(2) + 20.0_dp * y(2)
                    IF (y(2) <= 0.0_dp) THEN
                        grad(2) = grad(2) - 10.0_dp
                    ELSE
                        grad(2) = grad(2) + 10.0_dp
                    END IF
                END IF  
                
                IF (ind == 3) THEN
                    IF ( ( y(1) - y(2) ) <= 0.0_dp ) THEN
                        grad(1) = grad(1) - 10.0_dp 
                        grad(2) = grad(2) + 10.0_dp
                    ELSE
                        grad(1) = grad(1) + 10.0_dp
                        grad(2) = grad(2) - 10.0_dp                 
                    END IF
                    
                    IF (y(2) <= 0.0_dp) THEN
                        grad(2) = grad(2) - 10.0_dp                 
                    ELSE
                        grad(2) = grad(2) + 10.0_dp                 
                    END IF
                END IF  
                
                IF (ind == 4) THEN
                    grad(1) = grad(1) + 10.0_dp + 20.0_dp * y(1)
                    grad(2) = grad(2) + 20.0_dp * y(2)
                END IF
                
                IF ( (ABS(y(1)) - y(2) ) >= 0.0_dp ) THEN
                    grad(2) = grad(2) - 200.0_dp
                    IF ( y(1) <= 0.0_dp ) THEN
                        grad(1) = grad(1) - 200.0_dp
                    ELSE
                        grad(1) = grad(1) + 200.0_dp
                    END IF 
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
                
                grad = 0.0_dp
                
                IF (y(2) <= 0.0_dp ) THEN 
                         grad(2) = 20.0_dp * y(2) - 10.0_dp 
                    ELSE
                         grad(2) = 20.0_dp * y(2) + 10.0_dp 
                END IF
                
                grad(2) = grad(2) - 100.0_dp
                grad(1) = 20.0_dp * y(1)

                IF (y(1) <= 0.0_dp ) THEN 
                    grad(1) = grad(1) - 100.0_dp
                ELSE
                    grad(1) = grad(1) + 100.0_dp
                END IF
            
           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   


        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.8. in [1]                       | |
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
                REAL(KIND=dp), DIMENSION(4) :: a                ! help variable
                INTEGER :: i, ind                               ! help variable

                f = 9.0_dp - 8.0_dp * y(1) -6.0_dp * y(2) - 4.0_dp * y(3)  
                f = f + 2.0_dp * ABS(y(1)) + 2.0_dp * ABS(y(2))+ 2.0_dp * ABS(y(3)) 
                f = f + 4.0_dp * y(1)**2 + 2.0_dp * y(2)**2 + 2.0_dp * y(3)**2 

                a(1) = y(1) + y(2) + 2.0_dp * y(3) - 3.0_dp
                a(2) = -y(1)
                a(3) = -y(2)
                a(4) = -y(3)
                
                f = f + 10.0_dp * MAX(0.0_dp, a(1), a(2), a(3), a(4))
                
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
                
                f = ABS( y(1) -y(2) ) + ABS( y(1) -y(3) )

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
                REAL(KIND=dp), DIMENSION(5) :: a                ! help varaible
                INTEGER :: i, ind                               ! help variable
                
                grad(1) = -8.0_dp + 8.0_dp * y(1)
                grad(2) = -6.0_dp + 4.0_dp * y(2)
                grad(3) = -4.0_dp + 4.0_dp * y(3)

                IF ( y(1) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 2.0_dp
                    ELSE
                         grad(1) = grad(1) + 2.0_dp
                END IF                  
                
                IF ( y(2) <= 0.0_dp ) THEN 
                         grad(2) = grad(2) - 2.0_dp
                    ELSE
                         grad(2) = grad(2) + 2.0_dp
                END IF  

                IF ( y(3) <= 0.0_dp ) THEN 
                         grad(3) = grad(3) - 2.0_dp
                    ELSE
                         grad(3) = grad(3) + 2.0_dp
                END IF                  
                
                a(1) = y(1) + y(2) + 2.0_dp * y(3) - 3.0_dp
                a(2) = -y(1)
                a(3) = -y(2)
                a(4) = -y(3)
                a(5) = 0.0_dp
                
                ind = 5
                
                DO i = 1, 4
                   IF ( a(ind) < a(i) ) THEN 
                       ind = i
                   END IF 
                END DO
                
                IF (ind == 1) THEN
                    grad(1) = grad(1) + 10.0_dp 
                    grad(2) = grad(2) + 10.0_dp 
                    grad(3) = grad(3) + 20.0_dp 
                END IF
                
                IF (ind == 2) THEN
                    grad(1) = grad(1) - 10.0_dp 
                END IF  
                
                IF (ind == 3) THEN
                    grad(2) = grad(2) - 10.0_dp 
                END IF  
                
                IF (ind == 4) THEN
                    grad(3) = grad(3) - 10.0_dp 
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
                
                grad = 0.0_dp
                
                IF ( (y(1)-y(2)) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 1.0_dp 
                         grad(2) = grad(2) + 1.0_dp
                    ELSE
                         grad(1) = grad(1) + 1.0_dp
                         grad(2) = grad(2) - 1.0_dp 
                END IF
                
                IF ( (y(1)-y(3)) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 1.0_dp 
                         grad(3) = grad(3) + 1.0_dp
                    ELSE
                         grad(1) = grad(1) + 1.0_dp
                         grad(3) = grad(3) - 1.0_dp     
                END IF

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   


        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.9. in [1]                       | |
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
                REAL(KIND=dp) :: ti , f_help                    ! help variables
                INTEGER :: i, ind                               ! help variables

                f = y(1)**2 + (y(1)-1.0_dp)**2 + 2.0_dp*(y(1)-2.0_dp)**2 
                f = f + (y(1)-3.0_dp)**2 + 2.0_dp * y(2)**2 + (y(2)-1.0_dp)**2 
                f = f + 2.0_dp*(y(2)-2.0_dp)**2 + y(3)**2 + (y(3)-1.0_dp)**2 
                f = f + 2.0_dp*(y(3)-2.0_dp)**2 + (y(3)-3.0_dp)**2 
                f = f + 2.0_dp*y(4)**2 + (y(4)-1.0_dp)**2 + 2.0_dp*(y(4)-2.0_dp)**2 

                
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
                REAL(KIND=dp) :: ti , f_help, a, b              ! help variables
                INTEGER :: i, ind                               ! help variables

                f = 0.0_dp 

                a = (y(1) -2.0_dp)**2 + y(2)**2
                b = (y(3) -2.0_dp)**2 + y(4)**2
                
                IF ( a >= b) THEN 
                    f = f + a
                ELSE 
                    f = f + b
                END IF

                a = (y(1) -2.0_dp)**2 + (y(2)-1.0_dp)**2
                b = (y(3) -2.0_dp)**2 + (y(4)-1.0_dp)**2
                
                IF ( a >= b) THEN 
                    f = f + a
                ELSE 
                    f = f + b
                END IF

                a = (y(1) -3.0_dp)**2 + y(2)**2
                b = (y(3) -3.0_dp)**2 + y(4)**2
                
                IF ( a >= b) THEN 
                    f = f + a
                ELSE 
                    f = f + b
                END IF

                a = (y(1))**2 + (y(2)-2.0_dp)**2
                b = (y(3))**2 + (y(4)-2.0_dp)**2
                
                IF ( a >= b) THEN 
                    f = f + a
                ELSE 
                    f = f + b
                END IF

                a = (y(1)-1.0_dp)**2 + (y(2)-2.0_dp)**2
                b = (y(3)-1.0_dp)**2 + (y(4)-2.0_dp)**2
                
                IF ( a >= b) THEN 
                    f = f + a
                ELSE 
                    f = f + b
                END IF              

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
                REAL(KIND=dp) :: ti , f_help, f                 ! help variables
                INTEGER :: i, ind                               ! help variables

                grad(1) = 10.0_dp * y(1) - 16.0_dp
                grad(2) = 10.0_dp * y(2) - 10.0_dp 
                grad(3) = 10.0_dp * y(3) - 16.0_dp
                grad(4) = 10.0_dp * y(4) - 10.0_dp 
    
                
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
                REAL(KIND=dp) :: ti , f_help, f , a,b           ! help variables
                INTEGER :: i, ind                               ! help variables

                grad = 0.0_dp 

                a = (y(1) -2.0_dp)**2 + y(2)**2
                b = (y(3) -2.0_dp)**2 + y(4)**2
                
                IF ( a >= b) THEN 
                    grad(1) = grad(1) + 2.0_dp * (y(1)-2.0_dp)
                    grad(2) = grad(2) + 2.0_dp * y(2)
                ELSE 
                    grad(3) = grad(3) + 2.0_dp * (y(3)-2.0_dp)
                    grad(4) = grad(4) + 2.0_dp * y(4)
                END IF

                a = (y(1) -2.0_dp)**2 + (y(2)-1.0_dp)**2
                b = (y(3) -2.0_dp)**2 + (y(4)-1.0_dp)**2
                
                IF ( a >= b) THEN 
                    grad(1) = grad(1) + 2.0_dp * (y(1)-2.0_dp)
                    grad(2) = grad(2) + 2.0_dp * (y(2)-1.0_dp)
                ELSE 
                    grad(3) = grad(3) + 2.0_dp * (y(3)-2.0_dp)
                    grad(4) = grad(4) + 2.0_dp * (y(4)-1.0_dp)
                END IF

                a = (y(1) -3.0_dp)**2 + y(2)**2
                b = (y(3) -3.0_dp)**2 + y(4)**2
                
                IF ( a >= b) THEN 
                    grad(1) = grad(1) + 2.0_dp * (y(1)-3.0_dp)
                    grad(2) = grad(2) + 2.0_dp * (y(2))
                ELSE 
                    grad(3) = grad(3) + 2.0_dp * (y(3)-3.0_dp)
                    grad(4) = grad(4) + 2.0_dp * (y(4)) 
                END IF

                a = (y(1))**2 + (y(2)-2.0_dp)**2
                b = (y(3))**2 + (y(4)-2.0_dp)**2
                
                IF ( a >= b) THEN 
                    grad(1) = grad(1) + 2.0_dp * (y(1))
                    grad(2) = grad(2) + 2.0_dp * (y(2)-2.0_dp)
                ELSE 
                    grad(3) = grad(3) + 2.0_dp * (y(3))
                    grad(4) = grad(4) + 2.0_dp * (y(4)-2.0_dp)
                END IF

                a = (y(1)-1.0_dp)**2 + (y(2)-2.0_dp)**2
                b = (y(3)-1.0_dp)**2 + (y(4)-2.0_dp)**2
                
                IF ( a >= b) THEN 
                    grad(1) = grad(1) + 2.0_dp * (y(1)-1.0_dp)
                    grad(2) = grad(2) + 2.0_dp * (y(2)-2.0_dp)
                ELSE 
                    grad(3) = grad(3) + 2.0_dp * (y(3)-1.0_dp)
                    grad(4) = grad(4) + 2.0_dp * (y(4)-2.0_dp)
                END IF          

           END FUNCTION subgradient_f2  
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                            TEST PROBLEM 6.10. in [1]                      | |
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
                INTEGER :: i, ind                               ! help variable

                f = 0.0_dp

                DO i = 1, user_n
                   f = f + y(i)**2
                END DO
                
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
                
                DO i = 2, user_n
                    f = f + ABS( y(i) - y(i-1) )
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
                INTEGER :: i, ind                               ! help variable
                
                grad = 0.0_dp
                
                DO i = 1, user_n
                     grad(i) = grad(i) + 2.0_dp * y(i)
                END DO
    
                
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
                
                grad = 0.0_dp
                
                DO i = 2, user_n
                    IF ( y(i) - y(i-1) <= 0.0_dp) THEN
                        grad(i-1) = grad(i-1) + 1.0_dp 
                        grad(i) = grad(i) - 1.0_dp 
                    ELSE
                        grad(i-1) = grad(i-1) - 1.0_dp  
                        grad(i) = grad(i) + 1.0_dp      
                    END IF
                END DO
                
           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
           
               