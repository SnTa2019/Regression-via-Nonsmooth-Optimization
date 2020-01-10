      MODULE constants
        IMPLICIT NONE

        ! ** Double precision (i.e accuracy) **
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  
        
        ! ** Needed in the solver of the norm minimization problem ** 
        INTEGER :: next
        INTEGER, SAVE :: NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH
        
  
      END MODULE constants