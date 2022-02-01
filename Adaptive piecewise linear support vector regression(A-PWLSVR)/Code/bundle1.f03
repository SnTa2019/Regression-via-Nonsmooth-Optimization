!*************************************************************************
!*                                                                       *
!*     Bundle of DC component f_1. The code is part of DBDC method       *
!*     by Kaisa Joki.                                                    *
!*                                                                       *
!*************************************************************************

      MODULE bundle1  
      
        USE constants, ONLY : dp   ! double precision (i.e. accuracy)
        !USE omp_lib
        IMPLICIT NONE 
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |           THE BUNDLE ELEMENT AND THE BUNDLE OF THE DC COMPONENT F_1              | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   

        TYPE bundle_element1 ! bundle element of F_1
           PRIVATE
           REAL (KIND=dp), DIMENSION(:), POINTER  :: subgrad    ! subgradient of the bundle element
           REAL (KIND=dp) :: lin_error                          ! linearization error of the bundle element
        END TYPE bundle_element1        

        TYPE kimppu1 ! bundle of F_1
           PRIVATE
           TYPE(bundle_element1), DIMENSION(:), POINTER :: b_elements   ! bundle elements (does NOT contain the 'current element' and 'agg_element')
           TYPE(bundle_element1) :: current_element ! bundle element calculated at the current iteration point ('current element') 
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           TYPE(bundle_element1) :: agg_element ! the aggregated bundle element ('agg_element')
           
           INTEGER :: n         ! number of variables (also the length of subgradients)
           INTEGER :: b_maxsize ! 'maximum size of the bundle' - 1, (i.e. b_maxsize=size(b_elements) NOTICE: the 'current element' and 'agg_element' are stored separately)        
           INTEGER :: b_size    ! the current size of the bundle without the 'current element' and 'agg_element' (the actual size of the bundle is 'b_size+1' and 'agg_element is NOT taken into account in this value)
           INTEGER :: indeksi   ! the place where the next bundle element is tried to be added in the bundle element table 'b_elements'  
          
           LOGICAL :: full      ! tells whether the bundle is full or not             
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           LOGICAL :: agg       ! tells whether the aggregated bundle element was inserted into the bundle during the previous round           
        END TYPE kimppu1 


        CONTAINS
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                           CONTAINS SUBROUTINES:                                  | | 
        !| |                                                                                  | |
        !| |    INITIALIZATION         : init_bundle_b1(set, set_size, grad_length)           | |
        !| |    ADD ELEMENT            : add_element_b1(set, grad, alpha)                     | |
        !| |    ADD AGGREGATED ELEMENT : add_agg_element_b1(set, grad, alpha)                 | |
        !| |    ADD 1. CURRENT ELEMENT : add_first_element_b1(set, grad)                      | |
        !| |    UPDATE BUNDLE          : update_b1(set, new_grad, d, value_change)            | |
        !| |    RESET AGGREGATION IN BUNDLE : reset_agg_b1(set)                               | |		
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |              CONTAINS FUNCTIONS GIVING DIFFERENT VALUES:                         | |   
        !| |                                                                                  | |
        !| |    MATRIX OF SUBGRADIENTS      : grad_matrix(set)                                | |
        !| |    MATRIX OF LIN. ERRORS       : lin_error_matrix(set)                           | |
        !| |    MATRIX OF SUBGRADIENTS +AGG : grad_matrix_agg(set)                            | |
        !| |    MATRIX OF LIN. ERRORS  +AGG : lin_error_matrix_agg(set)                       | |       
        !| |    BUNDLE SIZE                 : give_size_b1(set)                               | |
        !| |    NUMBER OF VARIABLES         : give_n_b1(set)                                  | |
        !| |    IS BUNDLE FULL?             : is_full_b1(set)                                 | |
        !| |    IS AGGREGATION USED?        : is_agg_used(set)                                | |
        !| |    SUBGRADIENT OF ELEMENT i    : give_subgrad_b1(set, i)                         | | 
        !| |    LIN. ERROR OF ELEMENT i     : give_linerr_b1(set, i)                          | | 
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
           
           SUBROUTINE init_bundle_b1(set, set_size, grad_length) 
               !
               ! Initializes the bundle 'set'. Now the size of the bundle is 'set_size' and the length of subgradients is 'grad_size'.
               ! 
               ! 
               ! NOTICE: * 'grad_length' >= 1
               !         * IF (set_size < 2 ) THEN the size of the bundle is set to be 1 and only the 'current element' is stored (If aggregation is used, then also the aggregated element 'agg_element' is stored)               
               !         * 'set_size' does NOT include the 'aggregated element'. So if aggregation is used, then the actual size of the bundle is 'set_size+1'. 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set          ! bundle
               INTEGER, INTENT(IN):: set_size, grad_length  ! bundle size and the length of subgardients
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i, allocstat
               
                           
               IF (set_size < 2) THEN
                    set%b_maxsize = 0                       ! only the 'current element' and the 'agg_element' (if used) are stored 
                    set%full = .TRUE.
                
               ELSE     
                    set%b_maxsize = set_size - 1            ! the biggest possible size of the bundle without the 'current element' (the 'agg_element' is not taken into account here)
                    set%indeksi = 1
                    set%full = .FALSE.
               END IF
               
               set%b_size = 0                    ! the number of stored bundle elements in the table 'b_elements' ( ! without the 'current element' and 'agg_element' ! )  
               set%n = grad_length               ! the number of variables (this is also the length of subgradients)
               set%agg = .FALSE.
               
               ALLOCATE(set%b_elements(set%b_maxsize), STAT=allocstat)  ! initializes the maximum size of the bundle table 'b_elements' 
               ALLOCATE(set%current_element%subgrad(grad_length), &     ! initializes the length of the subgradient in the 'current element'
                         & STAT=allocstat)             
               ALLOCATE(set%agg_element%subgrad(grad_length), &         ! initializes the length of the subgradient in the 'aggregated element'
                         & STAT=allocstat)  

               DO i=1, set%b_maxsize
                    ALLOCATE(set%b_elements(i)%subgrad(grad_length), &  ! initializes the length of subgradients in the table 'b_elements'
                        & STAT=allocstat)     
               END DO 
               
           END SUBROUTINE init_bundle_b1
           
           

        !**************************************************************************************
        !                                                                                     |
        !                     ADD ELEMENT INTO TO THE BUNDLE                                  | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_element_b1(set, grad, alpha)
               !
               ! Adds the element '(grad, alpha)' into the bundle 'set' (i.e. into the bundle element table 'b_elements').
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.
               !         * IF the size of the bundle is 1, THEN nothing is added to the bundle element table 'b_elements'.         
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad   ! the subgradient                     
               REAL(KIND=dp), INTENT(IN) :: alpha                   ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i
               
               IF (set%b_maxsize > 0 ) THEN ! executed if bundle is larger than 1 (i.e. something can be stored into the table 'b_elements')
                   IF ( set%indeksi > set%b_maxsize ) THEN
                       set%indeksi = 1         
                   END IF
               
                   i = set%indeksi
                   set%b_elements(i)%subgrad = grad     ! adds the new subgradient into position i
                   set%b_elements(i)%lin_error = alpha  ! adds the new linearization error into position i
                   set%indeksi = i + 1                  ! the position where the next element is tried to be added
                   
                   IF ( .NOT. set%full ) THEN           ! if the bundle was not full during the previous round, then the size of the bundle is increased with 1
                       set%b_size = set%b_size + 1
                   END IF                  
                   
                   IF(set%b_size == set%b_maxsize) THEN  ! we test: Is the bundle full ?
                       set%full = .TRUE.
                   ELSE
                       set%full = .FALSE.
                   END IF 
                   
               END IF
               
           END SUBROUTINE add_element_b1
                 
           
           
        !**************************************************************************************
        !                                                                                     |
        !                     ADD AGGREGATED ELEMENT INTO TO THE BUNDLE                       | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_agg_element_b1(set, grad, alpha)
               !
               ! Adds the aggregated element '(grad, alpha)' into the bundle 'set'.
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.           
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad   ! the aggregated subgradient                      
               REAL(KIND=dp), INTENT(IN) :: alpha                   ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************
                           
               set%agg_element%subgrad = grad
               set%agg_element%lin_error = alpha 
               set%agg = .TRUE.
               
           END SUBROUTINE add_agg_element_b1           
       
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                  INITIALIZE/ADD THE FIRST CURRENT ELEMENT                           |
        !                                                                                     |     
        !**************************************************************************************            
           
           SUBROUTINE add_first_element_b1(set, grad)
               !
               ! Adds the element '(grad, 0)' calculated at the first iteration point x_0 into the bundle 'set'.
               !
               ! NOTICE: * the dimension of the 'grad' has to be 'set%n'.
               !         * the linearization error of the first current element is always zero.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                   ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad    ! the subgradient    
               !**************************** OTHER VARIABLES **************************************            
               
               set%current_element%subgrad = grad
               set%current_element%lin_error = 0.0_dp ! the linearization error is zero at the iteration point x_0
               
           END SUBROUTINE add_first_element_b1  

           
           
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE update_b1(set, new_grad, d, value_change)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this updates also all the other linearization errors in the bundle 'set' 
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * f1(x_{k+1}) - f1(x_k) is the 'value_change'             
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: new_grad   ! the subgradient calculated at the new iteration point               
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=dp), INTENT(IN) :: value_change                ! f1(x_{k+1}) - f1(x_k), i.e. the value change in the objective function
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i            
               
               ! The old 'current element' is added into the bundle set 'b_elements' and after that 
               ! the 'current element' can be updated with the new element
               ! (the linearization error of the current element is always zero, thus it is not changed).                  
               CALL add_element_b1(set, set%current_element%subgrad, 0.0_dp)
               set%current_element%subgrad = new_grad
           
               ! Linearization error update in the bundle set 'b_elements'
                
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error + &
                          & value_change - DOT_PRODUCT(set%b_elements(i)%subgrad, d)
               END DO
               
               IF (set%agg) THEN
                   set%agg_element%lin_error = set%agg_element%lin_error + &
                          & value_change - DOT_PRODUCT(set%agg_element%subgrad, d)           ! update of aggregated element  
               END IF
               
           
           END SUBROUTINE update_b1 
           
          
        !**************************************************************************************
        !                                                                                     |
        !                   RESET AGGREGATED ELEMENT IN THE BUNDLE                            |
        !                                                                                     |
        !************************************************************************************** 
        
           SUBROUTINE reset_agg_b1(set)
               !
               ! deletes the aggregated element from the bundle 'set' 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set ! bundle

               set%agg = .FALSE.            ! the aggregated element is deleted
               
           END SUBROUTINE reset_agg_b1   		  
		  
           
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                        FUNCTIONS GIVING DIFFERENT VALUES                         | |   
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   


        !**************************************************************************************
        !                                                                                     |
        !                             MATRIX OF SUBGRADIENTS                                  |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+1)'.
               !         * The subgradient corresponding to the current element is in the last position.
               !         * The aggregated element is NOT taken into account.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=dp), DIMENSION(set%n*(set%b_size+1)) :: m  ! the subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, j, length, start
               
               length = set%n               ! the length of subgradients
                       
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! the place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO
           END FUNCTION grad_matrix
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                           VECTOR OF LINEARIZATION ERRORS                            |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+1'.
               !         * The linearization error corresponding to the current element is in the last position.    
               !         * The aggregated element is NOT taken into account.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=dp), DIMENSION(set%b_size+1) :: m  ! the linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the last position

           END FUNCTION lin_error_matrix

        !**************************************************************************************
        !                                                                                     |
        !                         MATRIX OF SUBGRADIENTS WITH AGGREGATION                     |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix_agg(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' (with the aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+2)'.
               !         * The subgradient corresponding to the current element is in the 'last but one' position.
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=dp), DIMENSION(set%n*(set%b_size+2)) :: m  ! subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, j, length, start
               
               length = set%n               ! the length of subgradients
               
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO

               start = (set%b_size + 1) * length  ! the place where the 'aggregated element' is copied             
               DO j = 1, length 
                   m(start + j) = set%agg_element%subgrad(j)
               END DO
               
           END FUNCTION grad_matrix_agg
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                  VECTOR OF LINEARIZATION ERRORS WITH AGGREGATION                    |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix_agg(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' (with aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction. 
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+2'.
               !         * The linearization error corresponding to the current element is in the 'last but one' position.  
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=dp), DIMENSION(set%b_size+2) :: m  ! linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the 'last but one' position
               m(set%b_size+2) = set%agg_element%lin_error      ! the lin. error corresponding to the 'aggregated element' is copied into the last position

           END FUNCTION lin_error_matrix_agg           
        


           
        !**************************************************************************************
        !                                                                                     |
        !                                    BUNDLE SIZE                                      |
        !                                                                                     |
        !**************************************************************************************
       
           PURE FUNCTION give_size_b1(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element and the aggregated element ! )
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: bundle_size           ! size of the bundle
               
               bundle_size = set%b_size
               
           END FUNCTION give_size_b1
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              NUMBER OF VARIABLES                                    |
        !                                                                                     |
        !**************************************************************************************        
               
           PURE FUNCTION give_n_b1(set) RESULT(variable_number)
               !
               ! Gives the number of varibles in the minimization problem (this is also the length of subgradients)
               ! Number of varibles is (i.e. has to be) same as in kimppu2 when used in the algorithm.
               !
               IMPLICIT NONE 
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               INTEGER :: variable_number        ! the number of variables
               
               variable_number = set%n
               
           END FUNCTION give_n_b1          
           
           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS BUNDLE FULL?                                     |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_full_b1(set) RESULT(isfull)
               !
               ! Returns .TRUE. if bundle 'set' (i.e. the bundle element table 'b_element') is full otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isfull                 ! tells whether the bundle is full or not
 
               isfull = set%full 
               
            END FUNCTION is_full_b1         

           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS AGGREGATION USED?                                |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_agg_used(set) RESULT(isUsed)
               !
               ! Returns .TRUE. if aggregation is used in the bundle 'set' otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isUsed                 ! tells whether the aggregation is used or not
 
               isUsed = set%agg
               
            END FUNCTION is_agg_used            
        
        
        !**************************************************************************************
        !                                                                                     |
        !                            SUBGRADIENT OF ELEMENT i                                 |
        !                                                                                     |
        !**************************************************************************************
        
           FUNCTION give_subgrad_b1(set, i) RESULT(grad)
               !
               ! Gives the subgradient of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the subgradient of the 'current element'.    
               !         * If 'i=-1' then gives the subgradient of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)

               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set         ! bundle
               INTEGER, INTENT(IN) :: i                 ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp), DIMENSION(set%n) :: grad  ! subgradient at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! subgradient is from the bundle element table 'b_elements'
                   grad = set%b_elements(i)%subgrad
               ELSE IF (i == 0) THEN                        ! subgradient is from the 'current element'
                   grad = set%current_element%subgrad
               ELSE IF (i == -1) THEN                       ! subgradient is from the 'aggregated element'
                   grad = set%agg_element%subgrad                  
               ELSE                                         ! otherwise we are outside the bundle     
                   WRITE(*,*) 'CANNOT RETURN SUBGRADIENT! index ' &
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_subgrad_b1

           
           
        !**************************************************************************************
        !                                                                                     |
        !                          LINEARIZATION ERROR OF ELEMENT i                           |
        !                                                                                     |
        !**************************************************************************************        
    
           FUNCTION give_linerr_b1(set, i) RESULT(error)
               !
               ! Gives the linearization error of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the linearization error of the 'current element'.    
               !         * If 'i=-1' then gives the linearization error of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               INTEGER, INTENT(IN) :: i         ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: error           ! the linearization error at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! the linearization error is from the bundle element table 'b_elements'
                   error = set%b_elements(i)%lin_error
               ELSE IF (i==0) THEN                          ! the linearization error is from the 'current element'
                   error = set%current_element%lin_error
               ELSE IF (i==-1) THEN                         ! the linearization error is from the 'aggregated element'      
                   error = set%agg_element%lin_error               
               ELSE                                         ! otherwise we are outside the bundle
                   WRITE(*,*) 'CANNOT RETURN LINEARIZATION ERROR! index '&
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_linerr_b1         

      END MODULE bundle1
