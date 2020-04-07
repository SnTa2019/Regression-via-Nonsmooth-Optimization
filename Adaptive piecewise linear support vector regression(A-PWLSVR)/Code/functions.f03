!*************************************************************************
!*                                                                       *
!*     Computation of DC components f_1 and f_2 and their subgradients   *
!*     for the PWLSVR problem (The original code is part of DBDC method  *
!*     by Kaisa Joki, last modified 08.01.2020 by Napsu Karmitsa).       *
!*                                                                       *
!*************************************************************************

        MODULE functions    
      
        USE constants, ONLY : dp  ! double precision (i.e. accuracy)   
        USE initspr, ONLY : &
          ai => a4, &
          m => ntrain, &
          nft, &
          user_n => maxdim, & 
          k => maxmin2, &         ! k=1,...,maxmin
          maxlin, &
          q, &
          mk, &
          eps, &
          CC, &
          Cnull 
        USE initdbdc, ONLY : & 
          problem1, &
          problem2, &
          user_m, &
          user_c, & 
          user_r_dec, &       
          user_r_inc, &
          user_eps_1, &
          user_size_b1, &
          user_size_b2, &  
          user_crit_tol, &
          user_m_clarke, &
          user_eps, &
          user_size
                   
        IMPLICIT NONE
        
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                                                                          | |
        !| |    * PROBLEM specification:                                              | |
        !| |                                                                          | |       
        !| |        - The DC componen f1:                         'problem1'          | |
        !| |        - The DC componen f2:                         'problem2'          | |
        !| |        - the number of variables:                    'user_n'            | |               
        !| |        - the number of data points:                  'm'                 | |               
        !| |        - the number of features in data:             'nft'               | |               
        !| |        - the number of minimum functions:            'k'                 | |               
        !| |        - the number of lines in minimum functions:   'mk'                | |               
        !| |        - the total number of lines:                  'q'                 | |               
        !| |        - the epsilon used in SVM:                    'eps'               | |               
        !| |        - the coefficient in front of maximum terms:  'CC'                | |               
        !| |        - the coefficient in front of of the quadratic terms:  'Cnull'    | |               
        !| |        - the data matrix:                            'ai'                | |               
        !| |                                                                          | |
        !| |    * Different PARAMETERS:                                               | |
        !| |                                                                          | |       
        !| |        - the stopping tolerance:      'user_crit_tol'                    | |
        !| |                                                                          | |       
        !| |        MAIN ITERATION:                                                   | |       
        !| |        - the size of bundle B_1:      'user_size_b1'                     | |       
        !| |        - the size of bundle B_2:      'user_size_b2'                     | |       
        !| |        - the descent parameter:       'user_m'                           | |       
        !| |        - the decrease parameter:      'user_c'                           | |       
        !| |        - the decrease parameter:      'user_r_dec'                       | |       
        !| |        - the increase parameter:      'user_r_inc'                       | |           
        !| |        - the enlargement parameter:   'user_eps_1'                       | |       
        !| |                                                                          | |       
        !| |        CLARKE STATIONARY ALGORITHM:                                      | |       
        !| |        - the size of bundle B:        'user_size'                        | | 
        !| |        - the proximity measure:       'user_eps'                         | |               
        !| |        - the descent parameter:       'user_m_clarke'                    | |       
        !| |                                                                          | |       
        !| |                                                                          | |       
        !| |    * Computation of the value of the DC functions f_1 and f_2:           | |
        !| |        - f1(y, problem1)   the value of DC component f_1 at a point y    | |
        !| |        - f2(y, problem2)   the value of DC component f_2 at a point y    | |           
        !| |                                                                          | |               
        !| |                                                                          | |               
        !| |    * Computation of the subgradient of the DC components f_1 and f_2:    | |
        !| |        - subgradient_f1(y, problem1)    the subgradient of f_1 at y      | |
        !| |        - subgradient_f2(y, problem2)    the subgradient of f_2 at y      | |       
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  INFORMATION ABOUT PARAMETERS:                            | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------       
	
        !      OBJECTIVE FUNCTION:	
		! 
		!   0.5 * Cnull * sum_{l=1}^k      max      || x^[lj} ||^2  +  CC * sum_{i=1}^{m}  max (  0 ,  |    max           min       (  (x^{lj})^T a^i  + y_l  ) - b_i  | - eps  )
		!                              j=1,...,mk(l)                                                     l=1,...,k    j=1,...,mk(l)
        !
	        
        
         
        CONTAINS

        
        !********************************************************************************
        !                                                                               |
        !              FUNCTION VALUES OF THE DC COMPONENTS f_1 AND f_2                 |
        !                                                                               |
        !********************************************************************************

           FUNCTION f1(y, problem1) RESULT(f)        
                !
                ! Calculates the function value of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                 ! the objective function f_1 for which the value is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'
                REAL(KIND=dp) :: quadratic_sum                  
                REAL(KIND=dp) :: other_sum                  
                REAL(KIND=dp) :: largest_term, term                 
                INTEGER :: ind_term, largest_ind                         
                INTEGER :: term_num, term_start                           
                INTEGER :: i, j, l                            
                
                   
                SELECT CASE(problem1)

                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1)  
                       
                      f = 0.0_dp
					  
					  term_num = 0		
					  
					  ! the sum of the maximums from the quadratic terms
					  quadratic_sum = 0.0_dp
					  DO i = 1, k
					       
						   ! the first term in maximum
						   term = 0.0_dp
						   term_start = term_num * nft
						   DO j = 1, nft-1
						       term = term + y(term_start+j)*y(term_start+j)
						   END DO
						   largest_term = term
						   largest_ind = term_num
						   
						   DO l = 2, mk(i)
						      term_num = term_num + 1
						      term = 0.0_dp
						      term_start = term_num * nft
						      DO j = 1, nft-1
						          term = term + y(term_start+j)*y(term_start+j)
						      END DO
							  
							  IF (largest_term<term) THEN
							       largest_term = term
								   largest_ind = term_num
							  END IF
							  
                           END DO							  
						 
						   quadratic_sum = quadratic_sum + largest_term
					  
					  END DO

					  other_sum = 0.0_dp			  
					  DO i = 1, m
					     other_sum = other_sum + f1_part_i(y,i)
					  END DO
					  					 
					 f = Cnull * quadratic_sum + CC * other_sum
					 
                   !-------------------------------------  
                   
             

                
                END SELECT
                
           END FUNCTION f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f2(y, problem2) RESULT(f)           
                !
                ! Calculates the function value of DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                 ! the objective function f_2 for which the value is calculated               
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                          ! the function value of the DC component f_2 at a point 'y'
                REAL(KIND=dp) :: other_sum                  ! help variable
                INTEGER :: i                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                    CASE(1) 
					
					  other_sum = 0.0_dp			  
					  DO i = 1, m
					     other_sum = other_sum + f2_part_i(y,i)
					  END DO
					 
					 f =  CC * other_sum
					 !PRINT*,'hihu2',other_sum
					 
                   !-------------------------------------  
                   
                       
                
                END SELECT              
            


           END FUNCTION f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		   
		   
           
        !********************************************************************************
        !                                                                               |
        !         FUNCTION VALUES OF THE DC COMPONENTS f1_part_i AND f2_part_i          |
        !                                                                               |
        !********************************************************************************

           FUNCTION f1_part_i(y, ind) RESULT(f)        
                !
                ! Calculates the function value of the DC component f1_part_i at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f1_part_i is calculated
                INTEGER, INTENT(IN) :: ind                      ! the objective function f1_part_i for which the value is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f1_part_i at a point 'y'
                REAL(KIND=dp) :: part1, osa1, osa2             
                REAL(KIND=dp) :: osa1a, osa1b             
                REAL(KIND=dp) :: part2               
                REAL(KIND=dp) :: term_sum               
                REAL(KIND=dp) :: smallest_term                  
                REAL(KIND=dp) :: term                 
                INTEGER :: ind_term, smallest_ind                           
                INTEGER :: i, j, l                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1) 
				      smallest_term = convex_piece(y,ind,1)
                      smallest_ind = 1					  
							 		 
					  DO i = 2, k
					     term = convex_piece(y,ind,i)
						 IF (smallest_term > term) THEN
						    smallest_term = term
							smallest_ind = i
						 END IF
					  END DO
					  
					  ! PART 1
					  
					  osa1 = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind ) THEN
						    osa1 = osa1 
						 ELSE
						    osa1 = osa1 + 2.0_dp * convex_piece(y,ind,i)					  
					     END IF
					  END DO
					  osa1 = osa1 - ai(ind,nft)
					  
					  osa2 = 0.0_dp
					  DO i =1, k
					      osa2 = osa2 + 2.0_dp * convex_piece(y,ind,i)					  
					  END DO
					  osa2 = osa2 + ai(ind,nft)					  
					  
					  IF (osa1 > osa2) THEN
					     part1 = osa1
					  ELSE
					     part1 = osa2
					  END IF
					  
					  
					  ! PART 2 
					  part2 = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind) THEN 
						    part2 = part2 + convex_piece(y,ind,i)
						 ELSE
						    part2 = part2 + 2.0_dp * convex_piece(y,ind,i)						 
						 END IF
					  END DO
					 
					 part2 =  part2 + eps
					 
					 
					 ! f1_part_i
					 
					 IF (part1 > part2) THEN
					     f = part1
					 ELSE
					     f = part2
					 END IF
					 
                   !-------------------------------------  
                   
             

                
                END SELECT
                
           END FUNCTION f1_part_i      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f2_part_i(y, ind) RESULT(f)           
                !
                ! Calculates the function value of DC component f2_part_i at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f2_part_i is calculated
                INTEGER, INTENT(IN) :: ind                      ! the objective function f2_part_i for which the value is calculated               
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                          ! the function value of the DC component f2_part_i at a point 'y'
                REAL(KIND=dp) :: term_sum               
                REAL(KIND=dp) :: smallest_term                  
                REAL(KIND=dp) :: term                 
                INTEGER :: ind_term, smallest_ind                           
                INTEGER :: i, j, l                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1) 
				   
					  smallest_term = convex_piece(y,ind,1)
                      smallest_ind = 1					  
								 		 
					  DO i = 2, k
					     term = convex_piece(y,ind,i)
						 IF (smallest_term > term) THEN
						    smallest_term = term
							smallest_ind = i
						 END IF
					  END DO
					  
					  term_sum = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind) THEN 
						    term_sum = term_sum + convex_piece(y,ind,i)
						 ELSE
						    term_sum = term_sum + 2.0_dp * convex_piece(y,ind,i)						 
						 END IF
					  END DO
					 
					 f =  term_sum + eps
					 
                   !-------------------------------------  
                   
                       
                
                END SELECT              
            


           END FUNCTION f2_part_i
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		   
		   
		   
		   
        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS f_1 AND f_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_f1(y, problem1) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                 ! the objective function f_1 for which the subgradient is calculated                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: quad_grad                  
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: other_grad                  
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: l_grad                  
                REAL(KIND=dp) :: largest_term, term                 
                INTEGER :: ind_term, largest_ind                           
                INTEGER :: term_num, term_start                           
                INTEGER :: i, j, l                            
                
                SELECT CASE(problem1)

                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1)  
                       
                      grad = 0.0_dp
					  
					  term_num = 0		
					  
					  ! the sum of the maximums from the quadratic terms
					  quad_grad = 0.0_dp
					  DO i = 1, k
					       
						   ! the first term in maximum
						   term = 0.0_dp
						   term_start = term_num * nft
						   DO j = 1, nft-1
						       term = term + y(term_start+j)*y(term_start+j)
						   END DO
						   largest_term = term
						   largest_ind = term_num
						   
						   DO l = 2, mk(i)
						      term_num = term_num + 1
						      term = 0.0_dp
						      term_start = term_num * nft
						      DO j = 1, nft-1
						          term = term + y(term_start+j)*y(term_start+j)
						      END DO
							  
							  IF (largest_term<term) THEN
							       largest_term = term
								   largest_ind = term_num
							  END IF
							  
                           END DO				

						   term_start = largest_ind * nft
						   l_grad = 0.0_dp
                           DO l = 1, nft-1
                               l_grad(term_start+l) = l_grad(term_start+l)+ 2.0_dp*y(term_start+l)
                           END DO						   
						 
						   quad_grad = quad_grad + l_grad
					  
					  END DO

					  other_grad = 0.0_dp			  
					  DO i = 1, m
					     other_grad = other_grad + subgrad_f1_part_i(y,i)
					  END DO
					  					 
					 
					 grad = Cnull * quad_grad + CC * other_grad
					 
                   !-------------------------------------       
                   
 
                    
 
                END SELECT                      
                
           END FUNCTION subgradient_f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION subgradient_f2(y, problem2) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                 ! the objective function f_2 for which the subgradient is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: sum_grad   ! help variable
                INTEGER :: i                                    ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                    CASE(1) 
					
					  sum_grad = 0.0_dp			  
					  DO i = 1, m
					     sum_grad = sum_grad + subgrad_f2_part_i(y,i)
					  END DO
					 
					 grad =  CC * sum_grad
					 
                   !-------------------------------------  
                             
                   
                
                END SELECT  

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    


		   
        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS f_1 AND f_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgrad_f1_part_i(x, ind) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f1_part_i at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x    ! a point where the subgradient of the DC component f1_part_i is calculated
                INTEGER, INTENT(IN) :: ind                      ! the objective function f1_part_i for which the subgradient is calculated                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(x)) :: grad       ! the subgradient of the DC component f1_part_i at a point 'y'
                REAL(KIND=dp), DIMENSION(SIZE(x)) :: grad_sum   ! the subgradient of the DC component f1_part_i at a point 'y'
                REAL(KIND=dp) :: part1, osa1, osa2             
                REAL(KIND=dp) :: osa1a, osa1b             
                REAL(KIND=dp) :: part2               
                REAL(KIND=dp) :: term_sum               
                REAL(KIND=dp) :: smallest_term                  
                REAL(KIND=dp) :: term                 
                INTEGER :: ind_term, smallest_ind                           
                INTEGER :: i, j, l                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1) 
				   
					  smallest_term = convex_piece(x,ind,1)
                      smallest_ind = 1					  
								 		 
					  DO i = 2, k
					     term = convex_piece(x,ind,i)
						 IF (smallest_term > term) THEN
						    smallest_term = term
							smallest_ind = i
						 END IF
					  END DO
					  
					  ! PART 1
					  
					  osa1 = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind ) THEN
						    osa1 = osa1 
						 ELSE
						    osa1 = osa1 + 2.0_dp * convex_piece(x,ind,i)					  
					     END IF
					  END DO
					  osa1 = osa1 - ai(ind,nft)
					  
					  osa2 = 0.0_dp
					  DO i =1, k
					      osa2 = osa2 + 2.0_dp * convex_piece(x,ind,i)					  
					  END DO
					  osa2 = osa2 + ai(ind,nft)					  
					  
					  IF (osa1 > osa2) THEN
					     part1 = osa1
					  ELSE
					     part1 = osa2
					  END IF
					  
					  
					  ! PART 2 
					  part2 = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind) THEN 
						    part2 = part2 + convex_piece(x,ind,i)
						 ELSE
						    part2 = part2 + 2.0_dp * convex_piece(x,ind,i)						 
						 END IF
					  END DO
					 
					 part2 =  part2 + eps
					 
					 
					 ! f1_part_i
					 
					 grad_sum = 0.0_dp
					 IF (part1 > part2) THEN
					 
					     DO i = 1, k
						    grad_sum = grad_sum + 2.0_dp*convex_piece_grad(x,ind,i)
						 END DO
						 IF (osa1 > osa2) THEN
						     l = smallest_ind
						    grad_sum = grad_sum - 2.0_dp*convex_piece_grad(x,ind,l)
						 END IF
						 
					 ELSE
					 
					     DO i =1, k
					        IF (i == smallest_ind) THEN 
						       grad_sum = grad_sum + convex_piece_grad(x,ind,i)
						    ELSE
						       grad_sum = grad_sum + 2.0_dp*convex_piece_grad(x,ind,i)						 
						    END IF
					      END DO					     
					     
					 END IF
					 
					 
					 grad = grad_sum
					 
                   !-------------------------------------    
                   
 
                    
 
                END SELECT                      
                
           END FUNCTION subgrad_f1_part_i      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION subgrad_f2_part_i(x, ind) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f2_part_i at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x    ! a point where the subgradient of the DC component f2_part_i is calculated
                INTEGER, INTENT(IN) :: ind                      ! the objective function f2_part_i for which the subgradient is calculated
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(x)) :: grad       ! the subgradient of the DC component f2_part_i at a point 'y'
                REAL(KIND=dp), DIMENSION(SIZE(x)) :: grad_sum               
                REAL(KIND=dp) :: smallest_term                  
                REAL(KIND=dp) :: term                 
                INTEGER :: ind_term, smallest_ind                           
                INTEGER :: i, j, l                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1) 
				   
					  smallest_term = convex_piece(x,ind,1)
                      smallest_ind = 1					  
								 		 
					  DO i = 2, k
					     term = convex_piece(x,ind,i)
						 IF (smallest_term > term) THEN
						    smallest_term = term
							smallest_ind = i
						 END IF
					  END DO
					  
					  grad_sum = 0.0_dp
					  DO i =1, k
					     IF (i == smallest_ind) THEN 
						    grad_sum = grad_sum + convex_piece_grad(x,ind,i)
						 ELSE
						    grad_sum = grad_sum + 2.0_dp * convex_piece_grad(x,ind,i)						 
						 END IF
					  END DO
					 
					 grad = grad_sum
					 
                   !------------------------------------- 
                              
                   
                
                END SELECT  

           END FUNCTION subgrad_f2_part_i
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
		   
           
           
           FUNCTION convex_piece(x, ind, k_ind) RESULT(f)           
                !
                ! Calculates the convex function max_{j=1,...,M_k} { -(x^{kj})^T a^{ind} - y_{kj} } where k = k_ind
                !
                ! NOTICE: The dimension of 'x' has to be 'nft*q'. 
				!         'ind' needs to be an integer from interval [1,m]
				!         'k_ind' needs to be an integer from interval [1,k]				
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) ::  ind                     ! the data point 
                INTEGER, INTENT(IN) ::  k_ind                   ! the maximum function
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f 	                            ! the function value of the convex function
				
				REAL(KIND=dp) :: biggest                       ! the value of the line having the largest value
				REAL(KIND=dp) :: line_value                    ! the value of the line 
				INTEGER :: ind_biggest                         ! the index of line having the largest value
				INTEGER :: lines                               ! number of lines in the convex function
				INTEGER :: line_num                            ! line number
				INTEGER :: line_start                          ! a place where the line start 
                INTEGER :: j, i                                ! help variable
                
				!PRINT*,'kind',k_ind
                f = 0.0_dp				
				
				lines = mk(k_ind)           ! the number of lines in maximum function
				
				!the number of lines before the first line in the maximum function
				line_num = 0
				DO i =1, k_ind-1
				   line_num = line_num + mk(i)
				END DO 
				
				! the value of the first line
				line_value = 0.0_dp
				line_start = line_num * nft 
				DO j = 1, nft-1
                   line_value = line_value - x(line_start+j)*ai(ind,j)
				END DO
				line_value = line_value - x(line_start+nft)
				
				biggest = line_value 
				ind_biggest = line_num
                
				! the values of the other lines and the maximum value
				DO i = 2, lines
				   line_value = 0.0_dp
				   line_num = line_num + 1
				   line_start = line_num * nft				   
				   DO j = 1, nft-1
                       line_value = line_value - x(line_start+j)*ai(ind,j)
				   END DO
				   line_value = line_value - x(line_start+nft)		

				   IF (biggest < line_value) THEN
				      biggest = line_value
				      ind_biggest = line_num 
				   END IF
				   
				END DO
				
                f = biggest
            
           END FUNCTION convex_piece


           FUNCTION convex_piece_grad(x, ind, k_ind) RESULT(grad)           
                !
                ! Calculates a subgradient for the convex function max_{i=1,...,M_k} { -(x^{kj})^T a^{ind} - y_{kj} } where k = k_ind
                !
                ! NOTICE: The dimension of 'x' has to be 'nft*q'. 
				!         'ind' needs to be a n integer from interval [1,m]
				!         'k' needs to be a n integer from interval [1,k]
               !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x    ! a point where the function value of the convex piece is calculated
                INTEGER, INTENT(IN) ::  ind                       ! the data point 
                INTEGER, INTENT(IN) ::  k_ind                   ! the maximum function
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(x)) :: grad        ! the subgradient of the convex part
				REAL(KIND=dp) :: biggest                       ! the value of the line having the largest value
				REAL(KIND=dp) :: line_value                    ! the value of the line 
				INTEGER :: ind_biggest                         ! the index of line having the largest value
				INTEGER :: lines                               ! number of lines in the convex function
				INTEGER :: line_num                            ! the number of lines in the whole before the current line
				INTEGER :: line_start                          ! a place where the line start 
                INTEGER :: j, i                                ! help variable
                		
				lines = mk(k_ind)
				
				line_num = 0
				DO i =1, k_ind-1
				   line_num = line_num + mk(i)
				END DO 
				
				! the value of the first line
				line_value = 0.0_dp                                      
				line_start = line_num * nft 
				DO j = 1, nft-1
                   line_value = line_value - x(line_start+j)*ai(ind,j)
				END DO
				line_value = line_value - x(line_start+nft)
				
				biggest = line_value 
				ind_biggest = line_num
                
				!the values of the other line and the maximum value
				DO i = 2, lines
				   line_value = 0.0_dp
				   line_num = line_num + 1
				   line_start = line_num * nft				   
				   DO j = 1, nft-1
                       line_value = line_value - x(line_start+j)*ai(ind,j)
				   END DO
				   line_value = line_value - x(line_start+nft)		

				   IF (biggest < line_value) THEN
				      biggest = line_value
				      ind_biggest = line_num 
				   END IF
				   
                 END DO
				 
				 grad = 0.0_dp
				 
				 line_start = ind_biggest * nft
				 DO j = 1,nft-1
				     grad(line_start+j) = grad(line_start+j) - ai(ind,j)
				 END DO
				 grad(line_start+nft) = grad(line_start+nft) - 1.0_dp
				 				 
            
           END FUNCTION convex_piece_grad		   
		   
     
           
      END MODULE functions     
