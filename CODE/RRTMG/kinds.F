
      MODULE kinds

!***********************************************************************
!
!     This module defines variable precision for all common data
!     types.
!
!-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

!-----------------------------------------------------------------------

      INTEGER, PARAMETER :: &
          char_len  = 80, &
          log_kind  = KIND(.TRUE.), &
          int_kind  = SELECTED_INT_KIND  (09), &
          real_kind = SELECTED_REAL_KIND (06), &
          dbl_kind  = SELECTED_REAL_KIND (13)

      END MODULE kinds

