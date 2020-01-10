!*************************************************************************
!*                                                                       *
!*     Bundle of DC component f_1. The code is part of DBDC method       *
!*     by Kaisa Joki.                                                    *
!*                                                                       *
!*************************************************************************

      MODULE bundle2
        USE constants, ONLY : dp   ! double precision (i.e. accuracy)
        IMPLICIT NONE 

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |           THE BUNDLE ELEMENT AND THE BUNDLE OF THE DC COMPONENT F_2              | |  
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        TYPE bundle_element2 ! bundle element of F_2
           PRIVATE
           REAL (KIND=dp), DIMENSION(:), POINTER  :: subgrad   ! subgradient
           REAL (KIND=dp), DIMENSION(:), POINTER  :: direction ! search direction for the subproblem           
           REAL (KIND=dp) :: lin_error     ! linearization error
           REAL (KIND=dp) :: change        ! value of the predicted decrease  
           REAL (KIND=dp) :: subprob_value ! value of the subproblem objective
        END TYPE bundle_element2        

        TYPE kimppu2 ! bundle of F_2
           PRIVATE
           TYPE(bundle_element2), DIMENSION(:), POINTER :: b_elements ! bundle elements
           TYPE(bundle_element2) :: current_element ! bundle element at the current iteration point ('current element')
           
           INTEGER :: n         ! number of variables in vector x (also the length of subgradients)
           INTEGER :: b_maxsize ! 'maximum size of the bundle' - 1 ,(i.e. b_maxsize=size(b_elements) NOTICE: the current_element is stored separately)   
           INTEGER :: b_size    ! the current size of the bundle without the 'current element' (the actual size of the bundle is 'b_size+1')
           INTEGER :: glob_ind  ! the position of the bundle element giving the global solution    
           INTEGER :: indeksi   ! the place where the next element is tried to be added in the bundle element table 'b_elements'       
           
           LOGICAL :: full     ! tells whether this bundle is full or not          
        END TYPE kimppu2


        CONTAINS
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAIN SUBROUTINES:                                  | | 
        !| |                                                                                  | |
        !| |    INITIALIZATION               : init_bundle_b2(set, set_size, grad_length)     | |
        !| |    ADD ELEMENT                  : add_element_b2(set, grad, alpha)               | |
        !| |    FIRST CURRENT ELEMENT        : add_first_element_b2(set, grad)                | |
        !| |    UPDATE BUNDLE                : update_b2(set, new_grad, d, value_change)      | |
        !| |    SOLUTION FOR SUBPROBLEM i    : add_solution(set, i , d, delta, obj)           | |   
        !| |    ADD INDEX OF GLOBAL SOLUTION : add_glob_index(set)                            | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |              CONTAIN FUNCTIONS GIVING DIFFERENT VALUES:                          | |   
        !| |                                                                                  | |
        !| |    GLOBAL SOLUTION                   : give_solution(set)                        | |
        !| |    SOLUTION OF SUBPROBLEM i          : give_subprob_solution(set,i)              | |
        !| |    PREDICTED DECREASE                : give_decrease(set)                        | |
        !| |    PREDICT. DECRAESE OF SUBPROBLEM i : give_subprob_decrease(set,i)              | |
        !| |    INDEX OF SOLUTION                 : give_solution_ind(set)                    | |
        !| |    INDEX OF LAST ELEMENT             : give_last_element_ind_b2(set)             | |
        !| |    BUNDLE SIZE                       : give_size_b2(set)                         | |
        !| |    MAX BUNDLE SIZE                   : give_max_size_b2(set)                     | |
        !| |    NUMBER OF VARIABLES               : give_n_b2(set)                            | |
        !| |    IS BUNDLE FULL?                   : is_full_b2(set)                           | |
        !| |    SUBGRADIENT OF ELEMENT i          : give_subgrad_b2(set, i)                   | | 
        !| |    LIN. ERROR OF ELEMENT i           : give_linerr_b2(set, i)                    | | 
        !| |    MAXIMUM NORM VALUE                : max_norm_value(set)                       | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        
        
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                SUBROUTINES                                       | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*               
        

        !**************************************************************************************
        !                                                                                     |
        !                               INITIALIZATION                                        |
        !                                                                                     |
        !**************************************************************************************
           
           SUBROUTINE init_bundle_b2(set, set_size, grad_length) 
               !
               ! Initializes the bundle 'set'. Now the size of the bundle is 'set_size' and the length of subgradients is 'grad_length'.
               !
               ! NOTICE: * 'grad_length' >= 1
               !         * IF (set_size < 2 ) THEN the size of the bundle is set to be 1 and only the 'current element' is stored 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set          ! bundle
               INTEGER, INTENT(IN):: set_size, grad_length  ! the size of the bundle and the length of subgradients
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, allocstat
            
               IF (set_size < 2) THEN
                    set%b_maxsize = 0               ! only the 'current element' is stored 
                    set%full = .TRUE. 
               ELSE     
                    set%b_maxsize = set_size - 1    ! the biggest possible size of the bundle without the current element
                    set%indeksi = 1
                    set%full = .FALSE.
               END IF

               set%b_size = 0           ! the number of stored bundle elements in the table 'b_elements' ( ! without the current element ! )  
               set%n = grad_length      ! the number of variables (this is also the length of subgradients)
            
               ALLOCATE(set%b_elements(set%b_maxsize), STAT=allocstat)   ! initializes the maximum size of the bundle element table 'b_elements'

               DO i=1, set%b_maxsize
                    ALLOCATE(set%b_elements(i)%subgrad(grad_length), &   ! initializes the length of subgradients in the table 'b_elements'
                        & STAT=allocstat) 
                    ALLOCATE(set%b_elements(i)%direction(grad_length), & ! initializes the length of searh directions in the table 'b_elements'
                        & STAT=allocstat)       
               END DO  
               
               ALLOCATE(set%current_element%subgrad(grad_length), &      ! initialize the length of the subgradient in the 'current element'
                        & STAT=allocstat) 
               ALLOCATE(set%current_element%direction(grad_length), &    ! initialize the length of the searh direction in the 'current element'
                        & STAT=allocstat) 
                        
           END SUBROUTINE init_bundle_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                     ADD ELEMENT INTO TO THE BUNDLE                                  | 
        !                                                                                     |
        !**************************************************************************************
        
           SUBROUTINE add_element_b2(set, grad, alpha)
               !
               ! Adds the new element (grad, alpha) into the bundle 'set' (i.e. into the bundle element table 'b_elements').
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.
               !         * IF the size of the whole bundle is 1, THEN nothing is added to the bundle element table 'b_elements'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad ! the subgradient                       
               REAL(KIND=dp), INTENT(IN) :: alpha                 ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i
               
               IF (set%b_maxsize > 0 ) THEN ! executed if bundle size is larger than 0 (i.e. something can be stored into the table 'b_elements')
                   IF ( set%indeksi > set%b_maxsize ) THEN
                        set%indeksi = 1
                   END IF

                   i = set%indeksi  
                   
                   ! In the algorithm we use the case where the 'bundle element' yielding the previous global solution of the search direction problem cannot be replaced
                   
                   IF( set%full .AND. (i == set%glob_ind ) ) THEN   ! the bundle element which yields the previous global solution of 
                       i = i + 1                                    ! the search direction problem cannot be replaced
                       IF ( i > set%b_maxsize ) THEN                ! we make sure that the updated index is not outside the bundle element table 'b_elements'
                            i = 1
                       END IF 
                   END IF              
                   
                   set%b_elements(i)%lin_error = alpha   ! adds a new linearization error into the position i                      
                   set%b_elements(i)%subgrad = grad      ! adds a new subgradient into the position i

                   set%indeksi = i + 1                   ! the position where the next element is tried to be added
                   
                   IF ( .NOT. set%full ) THEN            ! if the bundle wasn't full during the previous round then the size of the bundle is increased with 1
                      set%b_size = set%b_size + 1
                   END IF
                               
                   IF(set%b_size == set%b_maxsize) THEN  ! we test: Is the bundle full ?
                       set%full = .TRUE.
                   ELSE
                       set%full = .FALSE.
                   END IF 
                   
               END IF
           END SUBROUTINE add_element_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                  INITIALIZE/ADD THE FIRST CURRENT ELEMENT                           |
        !                                                                                     |     
        !**************************************************************************************        
           
           SUBROUTINE add_first_element_b2(set, grad)
               !
               ! Adds the element '(grad, 0)' calculated at the first iteration point x_0 into the bundle 'set'.
               !
               ! NOTICE: * the dimension of the 'grad' has to be 'set%n'.
               !         * the linearization error of the first current element is always zero.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad ! subgradient at the first iteration point x_0                         
               
               set%current_element%subgrad = grad
               set%current_element%lin_error = 0.0_dp ! linearization error is zero at the iteration point x_0      
           
           END SUBROUTINE add_first_element_b2      
           
           
  
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************
        
           SUBROUTINE update_b2(set, new_grad, d, value_change)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this also the linearization errors are updated in the bundle
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * f2(x_{k+1}) - f2(x_k) is the 'value_change'
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: new_grad   ! subgradient calculated at the new iteration point                   
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=dp), INTENT(IN) :: value_change                ! f2(x_{k+1}) - f2(x_k), i.e. the value change in the objective function 
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i            
               
               ! the old 'current element' is added to the bundle set 'b_elements' and after that 
               ! the 'current_element' can be updated with the new element 
               ! (the linearization error of the current element is always zero, thus it is not changed)           
               CALL add_element_b2(set, set%current_element%subgrad, 0.0_dp)
               set%current_element%subgrad = new_grad
               
           
               !linearization error update on the bundle set 'b_elements'
           
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error + &
                          & value_change - DOT_PRODUCT(set%b_elements(i)%subgrad, d)
               END DO
           
           END SUBROUTINE update_b2
           
           
        !**************************************************************************************
        !                                                                                     |
        !            ADD DIRECTION, PREDICTED DECREASE AND SUBPROBLEM OBJECTIVE VALUE         |
        !                                FOR THE SUBPROBLEM i                                 |
        !                                                                                     |
        !************************************************************************************** 
        
            SUBROUTINE add_solution(set, i , d, delta, obj )
               !
               ! Adds the search direction 'd', the predicted decrease 'delta' (delta1+delta2)
               ! and the objective value 'obj' related to the subproblem 'i'
               !
               ! NOTICE: * 0 <= i <= set%b_size  (other indices are outside the current bundle)
               !         * If i=0 values are added to the current element
               !         * the dimension of 'd' has to be 'set%n'
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set              ! bundle
               INTEGER :: i                                     ! the index of the subproblem
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN) :: d ! the search direction
               REAL(KIND=dp), INTENT(IN) :: delta, obj          ! the predicted decraese and the objective value of the subproblem
               
               IF ( (i > set%b_size) .OR. (i < 0) ) THEN        ! subproblem index is outside (i.e. does not belong to) the bundle  
                   WRITE(*,*) 'CANNOT ADD SOLUTION! index ', i ,&
                         & 'outside the bundle (size:', set%b_size, ')'
               ELSE IF (i == 0) THEN                            ! values are added to the current element
                   set%current_element%direction = d
                   set%current_element%change = delta
                   set%current_element%subprob_value = obj                 
               ELSE                                             ! values are added to the bundle element table 'b_elements' (position is 'i')
                   set%b_elements(i)%direction = d
                   set%b_elements(i)%change = delta
                   set%b_elements(i)%subprob_value = obj
               END IF             
            END SUBROUTINE add_solution
            
            
           
        !**************************************************************************************
        !                                                                                     |
        !               ADD INDEX OF SUBPROBLEM YIELDING THE GLOBAL SOLUTION                  |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE add_glob_index(set) 
               !
               ! Calculates the index of the subproblem which 
               ! gives the global solution of the 'search direction' problem
               !
               ! NOTICE: * 'glob_index' is from the interval from '0' to 'b_size'.
               !         * IF 'glob_index' = 0, THEN the 'current element' gives the solution.
               !         * OTHERWISE the global solution is from the bundle element table 'b_elements'
               !           and 'glob_index' is the position of 'b_elements' which yields the global solution.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set   ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: ind, i                    
               
               IF (set%b_size == 0 ) THEN   ! true if only the 'current element' is in the bundle (i.e. we have only one element in the bundle).
                    set%glob_ind = 0        ! In this case the subproblem related to the 'current element' yields the global solution.
               ELSE            
                    ind = 1 
                    DO i = 2, set%b_size    ! calculation of the index of the element yielding the minimum value of the objective in the table 'b_elements'
                        IF ( set%b_elements(ind)%subprob_value &
                             & >  set%b_elements(i)%subprob_value ) THEN
                              ind = i
                        END IF
                    END DO

                    IF( set%b_elements(ind)%subprob_value &                ! the minimum value of the objective from the table 'b_elements' is compared
                             & > set%current_element%subprob_value ) THEN  ! with the 'current element'
                        ind = 0
                    END IF  
                    
                    set%glob_ind = ind
               END IF                   
           END SUBROUTINE add_glob_index
           
              

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                        FUNCTIONS GIVING DIFFERENT VALUES                         | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*                 

        
        !**************************************************************************************
        !                                                                                     |
        !                            GLOBAL SOLUTION VECTOR                                   |
        !                                                                                     |
        !**************************************************************************************     
        
           PURE FUNCTION give_solution(set) RESULT(solution)
               !
               ! Gives the search direction 'solution' which is the global solution.
               !
               ! NOTICE: * 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !           since otherwise the index of global solution is not right.
               !         * the dimension of 'solution' is 'set%n'
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set             ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp), DIMENSION(set%n) :: solution  ! global solution vector
               
               IF (set%glob_ind > 0) THEN    ! global solution is from the bundle element table 'b_elements'
                    solution = set%b_elements(set%glob_ind)%direction
               ELSE                          ! global solution is from the 'current_element'
                    solution = set%current_element%direction
               END IF
           END FUNCTION give_solution
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                           SOLUTION VECTOR OF SUBPROBLEM i                           |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_subprob_solution(set,i) RESULT(solution)
               !
               ! Gives the search direction 'solution' of subproblem 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= set%b_size (other indices are outside the current bundle B_2).
               !         * IF i=0 THEN gives the solution of the subproblem which is get by using the 'current element' of B_2.
               !         * The dimension of 'solution' is 'set%n'.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER, INTENT(IN) :: i         ! index of the subproblem
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp), DIMENSION(set%n) :: solution ! solution vector of the subproblem 'i'
               
               IF ( (i > 0) .AND. (i <= set%b_size) ) THEN   ! solution vector is from the bundle element table 'b_elements'  
                    solution = set%b_elements(i)%direction
               ELSE IF ( i == 0) THEN                        ! solution vector is from the 'current element'
                    solution = set%current_element%direction
               END IF
           END FUNCTION give_subprob_solution   
           


        !**************************************************************************************
        !                                                                                     |
        !                              PREDICTED DECREASE                                     |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_decrease(set) RESULT(dec)
               !
               ! Gives the value of predicted decrease (i.e. delta_1 + delta_2) at the global solution.
               !
               ! NOTICE: 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !         since otherwise the index of global solution is not right.            
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: dec             ! value of the predicted decrease at the global solution 

               IF (set%glob_ind > 0) THEN                       ! predicted decrease is from the bundle element table 'b_elements'
                    dec = set%b_elements(set%glob_ind)%change   
               ELSE                                             ! predicted decrease is from the 'current element'
                    dec = set%current_element%change          
               END IF
           END FUNCTION give_decrease       
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                         PREDICTED DECREASE OF SUBPROBLEM i                          |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_subprob_decrease(set,i) RESULT(dec)
               !
               ! Gives the value of predicted decrease 'dec' in the subproblem 'i' (i.e. delta_1 + delta_2 in the subproblem 'i').
               ! 
               ! NOTICE: * 0 <= 'i' <= set%b_size (other indices are outside the current bundle B_2).
               !         * IF i=0 THEN gives the predicted decrease of the subproblem which is get by using the 'current element' of B_2.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER, INTENT(IN) :: i         ! index of the subproblem   
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: dec             ! value of the predicted decrease of subproblem 'i'
               
               IF ( (i > 0) .AND. (i <= set%b_size) ) THEN  ! index is from the bundle element table 'b_elements'
                    dec = set%b_elements(i)%change 
               ELSE IF (i==0) THEN                          ! index is from the 'current element'
                    dec = set%current_element%change
               END IF
           END FUNCTION give_subprob_decrease              
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              INDEX OF SOLUTION                                      |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_solution_ind(set) RESULT(ind)
               !
               ! Gives the index 'ind' of the subproblem which yields the global solution.
               !
               ! NOTICE: * 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !           since otherwise the index of global solution is not right.                       
               !         * 'ind' is zero if the current element gives the solution.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: ind                   ! index of the global solution 
               
               ind = set%glob_ind
           END FUNCTION give_solution_ind
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                             INDEX OF THE LAST ELEMENT                               |
        !                                                                                     |
        !**************************************************************************************            
         
           PURE FUNCTION give_last_element_ind_b2(set) RESULT(ind)
               !
               ! Gives the index 'ind' of the place where the last bundle element was added in the bundle element table 'b_elements'.
               !
               ! NOTICE: The index 'ind' is zero if there is nothing in the bundle element table 'b_elements'.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: ind                   ! index of the place where the last element was added
               
               IF (set%b_size /= 0) THEN        ! there is something in the bundle element table 'b_elements'
                   ind = set%indeksi - 1 
               ELSE                             ! there in nothing in the bundle element table 'b_elements'
                   ind = 0
               END IF 
           END FUNCTION give_last_element_ind_b2    



        !**************************************************************************************
        !                                                                                     |
        !                                BUNDLE SIZE                                          |
        !                                                                                     |
        !************************************************************************************** 
       
           PURE FUNCTION give_size_b2(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element ! ).
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: bundle_size            ! size of the bundle
               
               bundle_size = set%b_size
           END FUNCTION give_size_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                               MAX BUNDLE SIZE                                       |
        !                                                                                     |
        !************************************************************************************** 
       
           PURE FUNCTION give_max_size_b2(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element ! ).
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: bundle_size            ! size of the bundle
               
               bundle_size = set%b_maxsize
           END FUNCTION give_max_size_b2           
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              NUMBER OF VARIABLES                                    |
        !                                                                                     |
        !**************************************************************************************            
       
           PURE FUNCTION give_n_b2(set) RESULT(variable_number)
               !
               ! Gives the number of varibles in the minimization problem (this is also the length of subgradients).
               ! Number of variables is (i.e. has to be) same as in kimppu1 when used in algorithm.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: variable_number        ! number of variables
               
               variable_number = set%n
           END FUNCTION give_n_b2          
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              IS BUNDLE FULL?                                        |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION is_full_b2(set) RESULT(isfull)
               !
               ! Returns .TRUE. if bundle is full otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************            
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isfull                ! tells whether the bundle is full or not
 
               isfull = set%full 
            END FUNCTION is_full_b2     
            
            
            
        !**************************************************************************************
        !                                                                                     |
        !                           SUBGRADIENT OF SUBPROBLEM i                               |
        !                                                                                     |
        !**************************************************************************************             
            
           FUNCTION give_subgrad_b2(set, i) RESULT(grad)
               ! 
               ! Gives the subgradient of the bundle element at position 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= 'set%b_size' (otherwise we are outside the bundle).
               !         * If 'i'=0 then gives the subgradient of the 'current element'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set         ! bundle
               INTEGER, INTENT(IN) :: i                 ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp), DIMENSION(set%n) :: grad  ! subgradient at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! subgradient is from the bundle element table 'b_elements'
                   grad = set%b_elements(i)%subgrad
               ELSE IF (i == 0) THEN                        ! subgradient is from the 'current element'
                   grad = set%current_element%subgrad
               ELSE                                         ! otherwise we are outside the bundle         
                   WRITE(*,*) 'CANNOT RETURN SUBGRADIENT! index ' &
                        & , i , 'outside the bundle (size without current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_subgrad_b2
                   
           
           
        !**************************************************************************************
        !                                                                                     |
        !                       LINEARIZATION ERROR OF SUBPROBLEM i                           |
        !                                                                                     |
        !**************************************************************************************            
 
           FUNCTION give_linerr_b2(set, i) RESULT(error)
               !
               ! Gives the linearization error of the bundle element at position 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= 'set%b_size' (otherwise we are outside the bundle).
               !         * If 'i'=0 then gives the linearization error of the 'current element'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER, INTENT(IN) :: i         ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: error           ! the linearization error at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! the linearization error is from the bundle element table 'b_elements'
                   error = set%b_elements(i)%lin_error
               ELSE IF (i==0) THEN                          ! the linearization error is from the 'current element'
                   error = set%current_element%lin_error
               ELSE                                         ! otherwise we are outside the bundle 
                   WRITE(*,*) 'CANNOT RETURN LINEARIZATION ERROR! index '&
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_linerr_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              MAXIMUM NORM VALUE                                     |
        !                                                                                     |
        !**************************************************************************************            
        
           PURE FUNCTION max_norm_value(set) RESULT(max_norm)
               !
               ! Gives the value of the maximum subgradient norm.
               !
               IMPLICIT NONE    
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: max_norm         ! value of the maximum subgradient norm
               REAL(KIND=dp) :: norm             ! 'help variable'
               INTEGER :: i
               
               max_norm = DOT_PRODUCT(set%current_element%subgrad, &   ! the square of the norm ||\bxi_2(x_k)|| 
                                         & set%current_element%subgrad) 
                                     
               DO i = 1, set%b_size
                   norm = DOT_PRODUCT(set%b_elements(i)%subgrad, &     ! the square of the norm ||\bxi_2(y_i)||
                                  & set%b_elements(i)%subgrad)
                   IF (max_norm < norm) THEN  
                        max_norm = norm
                   END IF
               END DO
               
               max_norm = SQRT(max_norm)  ! the maximum norm ||\bxi_{2,max}|| (WITHOUT square ! )
           END FUNCTION max_norm_value
           





      END MODULE bundle2