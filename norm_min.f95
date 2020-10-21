      MODULE norm_min

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                  THE SOLVER FOR THE NORM MINIMIZATION PROBLEM                    | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        !|                                                                                      |
        !|    Utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver.          |
        !|                                                                                      |
        !|    Utilizes PVMM by Ladislav Luksan as a norm minimization solver. This subroutine   | 
        !|    uses PQSUBS and MQSUBS by Ladislav Luksan. PVMM is a VARIABLE METRIC ALGORITHM    |
        !|    for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                          |
        !|                                                                                      |
        !|    The subroutine PVMM together with PQSUBS and MQSUBS are licensed by               |
        !|    the GNU Lesser General Public License (LGPL).                                     |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        USE constants   ! Double precision (i.e. accuracy) and varibles: NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH
        USE bundle1     ! The BUNDLE of the DC component f_1
        USE bundle2     ! The BUNDLE of the DC component f_2
        
        USE functions, ONLY : user_size_b1, user_size_b2, user_n  ! Contain INFORMATION from the USER:
                                                                  ! The biggest possible size of both bundles and also the number of variables.
                
        IMPLICIT NONE 
        
        EXTERNAL PVMM       ! Contains the Variable metric method PVVML by Ladislav Luksan and it is needed to solve 
                            ! the norm minimization problem at Step 3 of the 'main iteration' algorithm.
        
        EXTERNAL PQSUBS     ! Basic modules for PVVM (by Ladislav Luksan)
        EXTERNAL MQSUBS     ! Matrix modules for PVVM (by Ladislav Luksan) 
        
        INTEGER, PARAMETER :: m_size = user_size_b1 + user_size_b2  ! The biggest possible size of the matrix (and the vector) used to define the subgradient
                                                                    ! (and the objective function) of the norm minimization problem.  (NOTICE: aggregated elements of B_1 and B_2 are NOT used) 
        INTEGER, PARAMETER :: lkm = user_n                          ! The parameter 'lkm' is the number of variables in the original minimization problem
                                                                    ! and it is needed/used in the EXTERNAL subroutine FUNDER.

        ! Matrices which are used to define the objective function of the norm minimization problem and its gradient. 
        ! Both needed/used in the EXTERNAL subroutine FUNDER.
        REAL(KIND=dp), DIMENSION(lkm, m_size), SAVE :: f_matrix      ! used to define the objective function in the norm minimization problem.
        REAL(KIND=dp), DIMENSION(m_size, m_size), SAVE :: g_matrix   ! used to define gradient of the objective function in the norm minimization problem.
    
        
        CONTAINS
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |    SOLVER FOR NORM PROBLEM      : norm_solver(bxi, obj, NF, n1, n2, B1, B2)      | | 
        !| |                                                                                  | |
        !| |    FORMS MATRICES USED IN NORM PROBLEM : norm_matrix(dimensio, n1, n2, B1, B2)   | |
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*



        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                   SOLVER FOR THE NORM MINIMIZATION PROBLEM                     |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
        
            SUBROUTINE norm_solver(bxi, obj, NF, n1, n2, B1, B2)
                !
                ! Solves the norm minimization problem at Step 3 of the 'main iteration' algorithm.
                !
                ! INPUT:  * 'B1' and 'B2':  The bundles B_1 and B_2 of the DC component f_1 and f_2 
                !         * 'n1' and 'n2':  The size of the bundles B_1 and B_2 (also the 'current bundle element' is (i.e. has to be) taken into account in both bundles)
                !         * 'NF'         :  The number of variables in the objective function of the norm minimization problem (NOTICE: NF = n1 + n2)
                !
                ! OUTPUT: * 'bxi':  The solution vector bxi = bxi_1^* - bxi_2^* = f_matrix * lambda ('lambda' is approximation to the minimum)
                !         * 'obj':  The value of the objective funtion at the solution 'lambda' (i.e. '|| bxi ||')
                !
            !***********************************************************************************
                IMPLICIT NONE
                
            !**************************** NEEDED FROM USER *************************************    
                TYPE(kimppu1), INTENT(IN) :: B1  ! the bundle B_1 for the DC component f_1
                TYPE(kimppu2), INTENT(IN) :: B2  ! the bundle B_2 for the DC component f_2
                                
                REAL(KIND=dp), DIMENSION(give_n_b1(B1)), INTENT(OUT) :: bxi     ! Output: bxi = bxi_1^* - bxi_2^* = f_matrix * lambda
                REAL(KIND=dp), INTENT(OUT) :: obj                               ! Output: the value of the objective function at the solution 'lambda' 
                
                INTEGER, INTENT(IN) :: NF        ! the number of variables in objective function of the norm min problem (NOTICE: NF = n1 + n2)
                INTEGER, INTENT(IN) :: n1, n2    ! the size of the bundle B_1 and B_2  (also the 'current element' is (i.e. has to be) taken into account in both bundles)
                
            
            !***************************** LOCAL VARIABLES WHEN CALLING PVMML ******************


                ! .. Help variables ..                        
                INTEGER :: dimensio      ! the length of subgradients   
                INTEGER :: i,j, ind
                
                ! .. Other variables ..
                
                INTEGER, PARAMETER :: NB = 3  ! the simple bounds are accepted because NB > 0 
                INTEGER, PARAMETER :: NC = 2  ! the number of linear constraints (2 equality constraints)
                
                INTEGER :: IPRNT               ! Specifies print in PVMML:  IPRNT = 0  print is  suppressed
                                               !                            IPRNT = 1  basic print of final results
                                               !                            IPRNT = -1 extended print of final results 
                                               !                            IPRNT = 2  basic print of intermediate and final results
                                               !                            IPRNT = -2  extended print of intermediate and final results
                INTEGER :: ITERM               ! OUTPUT from PVMML: indicates the cause of termination              
                
                REAL(KIND=dp) :: GMAX, F                 ! OUTPUT from PVMML : F    - the value of the objective at the solution
                                                         !                     GMAX - maximum absolute value of a partial derivative of the objective function              
                
                INTEGER, DIMENSION(7) :: IPAR  ! INTEGER parameters (If IPAR = 0, then default values are used in PVMML)                
                
                REAL(KIND=dp), DIMENSION(NF) :: lambda  ! ON INPUT:  a vector with the initial estimate to the solution
                                                        ! ON OUTPUT: the approximation to the minimum               

                REAL(KIND=dp), DIMENSION(7) :: RPAR      ! REAL parameters (IF RPAR = 0.0_dp, then default values are used in PVMML)


                ! .. Array Arguments in PVMML..
                INTEGER, DIMENSION(NF)       :: IX     ! the vector containing types of simple bounds (IX=1 because XL(I) <= X(I) for each variable)
                REAL(KIND=dp), DIMENSION(NF) :: XL     ! the vector containing lower bounds (XL(I) = 0.0_dp for each variable)
                REAL(KIND=dp), DIMENSION(NF) :: XU     ! vector containing upper bounds     (This is not needed)
                
                REAL(KIND=dp), DIMENSION(NC) :: CF     ! optional/auxiliary: the vector containing values of constraint functions
                INTEGER, DIMENSION(NC)       :: IC     ! the vector containing constraint types: 5 - the equality constraint (NOTICE: CL(I) = C(I) = CU(I) )
                REAL(KIND=dp), DIMENSION(NC) :: CL     ! the vector containing lower bounds for constraints 
                REAL(KIND=dp), DIMENSION(NC) :: CU     ! the vector containing upper bounds for constraints
                REAL(KIND=dp), DIMENSION(NF*NC) :: CG  ! the vector whose columns are normals of linear contsraints     

                
                ! .. Initialization ..
                IPAR = 0         ! default values are used
                RPAR = 0.0_dp    ! default values are used
                IPRNT = 0        ! print (because IPRNT=0, print is suppressed)
                
                IX = 1            ! the simple bounds are lower bounds
                XL = 0.0_dp       ! values of the lower bounds
                XU = 10.0_dp      ! values of the upper bounds (NOT needed here!)
                
                !CF = (/  /)                 ! optional         
                IC = (/ 5, 5 /)              ! constraints are equality constraints
                CL = (/ 1.0_dp, 1.0_dp /)    ! lower bounds for constraints
                CU = (/ 1.0_dp, 1.0_dp /)    ! upper bounds for constraints
                
                CG = 0.0_dp         ! Initialization of CG          
                DO i = 1, n1        ! the first constraint
                    CG(i) = 1.0_dp              
                END DO
                
                ind = 2*n1 + n2 + 1  
                DO i = ind, (ind + n2)  ! the second constraint
                    CG(i) = 1.0_dp
                END DO
                
                ! initial estimate to the solution
                lambda = 0.0_dp         ! initialization of lambda
                lambda(1) = 1.0_dp
                lambda(n1+1) = 1.0_dp
                
                ! the length of subgradients (same as the number of variables in the original minimization problem).
                dimensio = give_n_b1(B1) 

                CALL norm_matrix(dimensio, n1, n2, B1, B2) ! Calculates the matrices used to define the objective function and 
                                                           ! its gradient in the norm minimization problem


                ! Execution of the solver PVMML by Ladislav Luksan
                CALL PVMML(NF,NB,NC,lambda,IX,XL,XU,CF,IC,CL,CU,CG,IPAR,&     ! objective F = || f_matrix * lambda ||^2
                            & RPAR, F, GMAX, IPRNT, ITERM )
          

                IF ( F < 0.0_dp) THEN
                    obj = 0.0_dp
                ELSE
                    obj = SQRT(F)    ! the original objective function 'obj' = || f_matrix * lambda || 
                END IF

                bxi = 0.0_dp         ! initialization
                
                DO j = 1, NF
                    DO i = 1, dimensio
                       bxi(i) = bxi(i) + f_matrix(i,j)*lambda(j)
                    END DO
                END DO 
                            
         
            END SUBROUTINE norm_solver
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
            
            
            
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |            CALCULATE MATRICES f_matrix and g_matrix USED IN FUNDER             |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
            SUBROUTINE norm_matrix(dimensio, n1, n2, B1, B2) 
             ! Calculates the matrices used to define the objective of norm minimization problem
             ! and its gradient. These matrices are saved to the variables 'f_matrix' and 'g_matrix'
             !
             ! NOTICE: * 'dimensio'     : the length of the subgradient
             !         * 'n1' and 'n2'  : the size of the bundles B_1 and B_2 (the 'current bundle element' is also taken into account in both bundles) (aggregation is NOT used)
             !         * 'B1' and 'B2*  : the bundles B_1 and B_2 of the DC components f_1 and f_2
             !
                 INTEGER, INTENT(IN) :: dimensio, n1, n2  ! dimensio  - the length of the subgradients
                                                          ! n1 and n2 - the size of the bundles B_1 and B_2 (the 'current bundle element' is also taken into account in both bundles)
                                                                            
                 INTEGER :: i, j, k, ind
     
                 ! .. Help variables ..
                 REAL(KIND=dp), DIMENSION(dimensio) :: grad_b1, grad_b2                  
                 REAL(KIND=dp), DIMENSION(n1+n2,n1+n2) :: M             ! M = grad^T * grad  
                 REAL(KIND=dp), DIMENSION(dimensio, n1+n2) :: grad      ! the matrix whose columns are subgradients
                                                                        ! contains all the subgradients (both from B_1 and B_2)     (aggregated element is NOT included)
                                                                        ! first there are the ones from B_1 and then the ones from B_2

                 ! .. Bundles ..                                                        
                 TYPE(kimppu1), INTENT(IN) :: B1          ! the bundle of the DC component f_1
                 TYPE(kimppu2), INTENT(IN) :: B2          ! the bundle of the DC component f_2                                                                      

                 
                 ! .. Matrix 'grad' is formed ..
                 DO i = 0, give_size_b1(B1)             ! all subgradients from the bundle B_1 are looked through (except for the aggregated element)
                    grad_b1 = give_subgrad_b1(B1,i)     ! the subgradient of the bundle element i from the bundle B_1
                    DO j = 1, dimensio                  ! inserts the current grad_b1 into the matrix 'grad' (the starting position is (i+1):th column)
                       grad(j, i+1) = grad_b1(j)
                    END DO                    
                 END DO 
               
                 ind = give_size_b1(B1) + 2    ! the starting position for subgradients from the bundle B_2
               
                 DO i = 0, give_size_b2(B2)             ! all subgradients from the bundle B_2 are looked through (except for the aggregated element)
                    grad_b2 = give_subgrad_b2(B2,i)     ! the subgradient of the bundle element i from the bundle B_2
                    DO j = 1, dimensio                  ! inserts the current grad_b2 into the matrix 'grad' (NOTICE: the negative sign)
                       grad(j, ind + i) = - grad_b2(j)  ! (the position is (ind+i):th column in the matrix)
                    END DO                
                 END DO         
                 
                 ! .. Matrix 'grad' is saved to 'f_matrix' ..
                 DO j = 1, n1+n2
                    DO i = 1, dimensio
                      f_matrix(i,j) = grad(i,j)
                    END DO
                 END DO              
                 
                 M =  MATMUL( TRANSPOSE(grad), grad )      ! M = grad^T * grad
                 
                 ! .. Matrix 'M' is saved to 'g_matrix' ..
                 DO j = 1, n1+n2
                    DO i = 1, n1+n2
                      g_matrix(i,j) = M(i,j)
                    END DO
                 END DO
                                     
             
            END SUBROUTINE norm_matrix      
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
        

      END MODULE norm_min        