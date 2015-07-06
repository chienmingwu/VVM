      MODULE timer  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This module contains utilities for timing individual blocks of code.
!
! The module is comprised of three subroutines.  These are discussed below:
!
! 1) timer_start.  This subroutine starts a timer.
!
!    input:   event_name   CHARACTER (LEN=72)  name of event
!             clear_timer  LOGICAL,OPTIONAL    clear accumulated times
!
!      The block of code being timed is associated with an event.  
!      The subroutine argument is a character string called event_name.  
!      The character string provides a name for the event (or block of code)
!      to be timed.   A timer event has a logical attribute called "running" 
!      to indicate if the timer is "on" or "off".  This is analogous to 
!      a stopwatch being "on" or "off".
!
!      If the event does not already exist, its start time is initialized 
!      to the current time, and the timer is turned on.  On the other hand, 
!      if the event already exists and the timer is off, then the start 
!      time is re-initialized to the current time, and the timer is turned 
!      on.  If it already exists and is on, nothing happens.
!
!      The optional argument, clear_timer, allows the accumulated times 
!      associated with the event to be cleared.
!
! 2) timer_stop.   this subroutine stops a timer.
!
!    input:   event_name   CHARACTER (LEN=72) name of event
!
!      This subroutine records the current time and turns the timer off.
!      It subtracts the start time of the input event from the current time 
!      and records accumulated time for the event.
!
! 3) timer_write.  this subroutine writes a list of the timings.
!
!    input:   event_name   CHARACTER (LEN=72) name of event
!
!      This subroutine writes a list of times for each event.
!     
!      If an event called "total time" has been initiated, then
!      the subroutine writes the fractional time of the total time 
!      for each event.
!
! For example...
!
!     CALL timer_start ("total time")
!     DO i = 1,100
!     
!        CALL timer_start ("event 1")
!            < code block 1 >
!        CALL timer_stop  ("event 1")
!
!        CALL timer_start ("event 2")
!            < code block 2 >
!        CALL timer_stop  ("event 2")
!     ENDDO
!
!     CALL timer_stop ("total_time")
!     CALL timer_write ()
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE kinds
!      USE component_intra ! will need this for parallel
      IMPLICIT NONE
      SAVE
      PRIVATE

!#if mpi == 1
!        include 'mpif.h'
!#endif

      TYPE timer_node
         CHARACTER (LEN=72) :: event_name
         LOGICAL (KIND=log_kind) :: running
         REAL (KIND= dbl_kind) ::  cpu_start, cpu_stop, cpu_total
         REAL (KIND= dbl_kind) :: wall_start,wall_stop,wall_total
         TYPE (timer_node),POINTER :: next
      END TYPE timer_node

      TYPE (timer_node),POINTER :: timer_events

! public member functions
      PUBLIC ::                                                      &
         timer_start,                                                &
         timer_stop,                                                 &
         timer_write
      CONTAINS
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc
      SUBROUTINE timer_start (event_name,clear_timer)
! 19 SEPTEMBER 2002   
      CHARACTER (LEN=*),INTENT (IN) :: event_name
      LOGICAL (KIND=log_kind),OPTIONAL,INTENT (IN) :: clear_timer
      LOGICAL (KIND=log_kind) :: event_added,event_found,string_equal
      INTEGER (KIND=int_kind) :: clock,hz
      TYPE (timer_node),POINTER :: current,temp

      event_added = .FALSE.
      event_found = .FALSE.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! THE FIRST TIME timer_start IS CALLED, INITIALIZE LIST
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (.NOT.ASSOCIATED (timer_events)) THEN
         event_added = .TRUE.
         ALLOCATE (timer_events); ALLOCATE (timer_events%next)
         timer_events%event_name = ' '
         current => timer_events%next
         NULLIFY (current%next)
      ELSE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! LOOK THROUGH THE LIST OF EVENTS TO FIND EVENTS WHICH ALREADY EXIST
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         current => timer_events
         DO WHILE ((.NOT.event_found).AND.(ASSOCIATED (current)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF AN EVENTS NAME IS ALREADY ON THE LIST, THEN IT ALREADY EXISTS.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            string_equal = (TRIM (current%event_name)==TRIM (event_name))
            IF (string_equal) THEN
               event_found = .TRUE.
            ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF NOT FOUND, THEN KEEP LOOKING
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (.NOT.event_found) THEN
               current => current%next
            ENDIF
         ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF NOT FOUND, THEN LOOK THROUGH LIST TO FIND POSITION TO ADD NEW EVENT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF (.NOT.event_found) THEN
            current => timer_events
            DO WHILE ((.NOT.event_added).AND.(ASSOCIATED (current%next)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ADD A NEW NODE IN THE MIDDLE OF THE LIST
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               IF (LLT (TRIM (current%event_name),TRIM (event_name)).AND.    &
                   LLT (TRIM (event_name),TRIM (current%next%event_name))) THEN
                  event_added = .TRUE.
                  temp => current%next
                  NULLIFY (current%next); ALLOCATE (current%next)
                  current => current%next
                  current%next => temp
               ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF NOT ADDED, THEN KEEP LOOKING
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               IF (.NOT.event_added) THEN
                  current => current%next
               ENDIF
            ENDDO
         ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF THE END OF THE LIST IS REACHED AND NOT ADDED 
! AND NOT DOES NOT ALREADY EXIST THEN ADD TO END OF LIST
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF ((.NOT.event_added).AND.(.NOT.event_found)) THEN
            event_added = .TRUE.
            ALLOCATE (current%next)
            current => current%next
            NULLIFY (current%next)
         ENDIF
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF NEW EVENT ADDED, THEN INITIALIZE STUFF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (event_added) THEN
         current%event_name = event_name
         current%running    = .FALSE.
         current% cpu_total = 0.0
         current%wall_total = 0.0
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! If (NEW EVENT ADDED) OR (THE EVENT WAS FOUND ON THE LIST BUT IS NOT
! CURRENTLY running) THEN TURN TIMER ON AND GET THE BEGINNING TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF ((event_added).OR.((event_found).AND.(.NOT.current%running))) THEN
         current%running   = .TRUE.
         current% cpu_start = 0.0 ! figure this out later

         CALL system_clock (count=clock)
         CALL system_clock (count_rate=hz)
         current%wall_start = REAL (clock)/REAL (hz)
      ENDIF

      IF (PRESENT (clear_timer)) THEN
         IF (clear_timer) THEN
            current% cpu_start = 0.0
            current% cpu_stop  = 0.0
            current% cpu_total = 0.0
            current%wall_start = 0.0
            current%wall_stop  = 0.0
            current%wall_total = 0.0
         ENDIF
      ENDIF

      END SUBROUTINE timer_start
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc
      SUBROUTINE timer_stop (event_name)
! 19 SEPTEMBER 2002    
      LOGICAL (KIND=log_kind) :: event_found,string_equal
      INTEGER (KIND=int_kind) :: clock,hz
      CHARACTER (LEN=*), intent(in) :: event_name
      TYPE (timer_node),POINTER :: current

      event_found = .FALSE.

      IF (.NOT.ASSOCIATED (timer_events)) THEN
         PRINT *,' timer_stop :: timer_stop called with no events initiated '
         STOP
      ELSE
         current => timer_events
         DO WHILE ((.NOT.event_found).AND.(ASSOCIATED (current)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! THE EVENT IS FOUND
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            string_equal = (TRIM (current%event_name)==TRIM (event_name))
            IF (string_equal) THEN
               event_found = .TRUE.
            ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF NOT FOUND THEN KEEP LOOKING
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (.NOT.event_found) THEN
               current => current%next
            ENDIF
         ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF THE END OF THE LIST IS REACHED AND EVENT NOT FOUND THEN ERROR
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF (.NOT.event_found) THEN
            PRINT *,' timer_stop :: timer_stop called with event_name = "',  &
                                                           event_name,'"'
            PRINT *,'               this event has not been initiated  '
            STOP
         ENDIF
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF THE TIMER IS CURRENTLY running, THEN TURN THE TIMER OFF, 
! GET THE STOPPING TIME AND TOTAL ACCUMULATED TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (current%running) THEN
         current%running   = .FALSE.
      current% cpu_stop = 0.0 ! figure this out later

      CALL system_clock (count=clock)
      CALL system_clock (count_rate=hz)
      current%wall_stop = REAL (clock)/REAL (hz)

      current% cpu_total = current% cpu_total +                              &
                                  (current% cpu_stop-current% cpu_start)
      current%wall_total = current%wall_total +                              &
                                  (current%wall_stop-current%wall_start)
      ENDIF

      END SUBROUTINE timer_stop
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc
      SUBROUTINE timer_write ( component_name )
! 19 SEPTEMBER 2002   
      CHARACTER(len=*), intent(in) :: component_name

      LOGICAL (KIND=log_kind) :: total_found,string_equal
      INTEGER (KIND=int_kind) :: task, my_task, npe
!#if mpi == 1
!      INTEGER (KIND=int_kind) :: ierr
!#endif
      INTEGER (KIND=int_kind) :: clock,hz
      REAL (KIND=dbl_kind) ::  cpu_now, cpu_save, cpu_temp
      REAL (KIND=dbl_kind) :: wall_now,wall_save,wall_temp
      REAL (KIND=dbl_kind) :: frac

      TYPE (timer_node),POINTER :: current
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GET THE CURRENT TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      cpu_now = 0.0 ! figure this out later

      CALL system_clock (count=clock)
      CALL system_clock (count_rate=hz)
      wall_now = REAL (clock)/REAL (hz)


!#if mpi==1
!      CALL MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
      IF (.NOT.ASSOCIATED (timer_events)) THEN
         PRINT *,' timer_write :: timer_write called with no events initiated '
         STOP
      ELSE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! LOOK FOR AN EVENT CALLED 'total time'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         total_found = .FALSE.
         current => timer_events%next
         DO WHILE ((.NOT.total_found).AND.(ASSOCIATED (current)))
            string_equal = (TRIM (current%event_name)=='total time')
            IF (string_equal) THEN
               total_found = .TRUE.
               IF (current%running) THEN
                   cpu_save = current% cpu_total +                           &
                                           ( cpu_now-current% cpu_start)
                  wall_save = current%wall_total +                           &
                                           (wall_now-current%wall_start)
               ELSE
                   cpu_save = current% cpu_total
                  wall_save = current%wall_total
               ENDIF
            ENDIF
            IF (.NOT.total_found) THEN
               current => current%next
            ENDIF
         ENDDO

!         my_task = get_my_task(component_name)
!         npe = get_my_npe(component_name)
         my_task = 0
         npe = 1
         DO task = 0,npe-1
            IF (task==my_task) THEN
               PRINT  *,' '
               IF (total_found) THEN
                  IF (my_task==0) THEN
                     PRINT  *,' TIMINGS (process:event,running,',            &
                                       'cpu,wall,100*(wall/total wall))'
                  ENDIF
                  PRINT 20,task,current%event_name,current%running,          &
                                                      cpu_save,wall_save
                  PRINT  *,' '
               ELSE
                  IF (my_task==0) THEN
                     PRINT  *,' TIMINGS (process:event,running,',            &
                                                             'cpu,wall)'
                  ENDIF
               ENDIF

               current => timer_events%next
               DO WHILE (ASSOCIATED (current))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF THE TIMER IS RUNNING, COMPUTE THE ACCUMULATED TIME USING THE CURRENT TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  IF (current%running) THEN
                      cpu_temp = current% cpu_total +                        &
                                           ( cpu_now-current% cpu_start)
                     wall_temp = current%wall_total +                        &
                                           (wall_now-current%wall_start)
                  ELSE
                      cpu_temp = current% cpu_total
                     wall_temp = current%wall_total
                  ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IF 'total time' EVENT EXISTS, THEN COMPUTE WALL TIME FRACTION
! USING TOTAL TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  IF (total_found) THEN 
                     string_equal = (TRIM (current%event_name)=='total time')
                     IF (.NOT.string_equal) THEN
                        frac = 100.*wall_temp/(wall_save+1.E-10_dbl_kind)
                        PRINT 10, task,current%event_name,current%running,   &
                                                 cpu_temp,wall_temp,frac
                     ENDIF
                  ELSE
                     PRINT 20, task,current%event_name,current%running,      &
                                                 cpu_temp,wall_temp
                  ENDIF
                  current => current%next
               ENDDO
            ENDIF
!#if mpi==1
!            CALL MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
         ENDDO
      ENDIF
!#if mpi==1
!      CALL MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif

   10 FORMAT (1x,i5,' : ',a20,l1,2f15.5,f8.2)
   20 FORMAT (1x,i5,' : ',a20,l1,2f15.5     )

      END SUBROUTINE timer_write
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc
      CHARACTER*03 FUNCTION integer_to_string (i) RESULT (string)
      INTEGER :: i,n
      REAL*8 :: xi
      CHARACTER*01 :: char_i(3)
      CHARACTER*01 :: integers(12)

      integers = (/'0','1','2','3','4','5','6','7','8','9','0','1'/)

      xi = FLOAT (i)/1000.0
      DO n = 1,3
         char_i(n) = integers(INT (10.*xi+0.00001)+1)
         xi = 10.*xi - FLOAT (INT (10.*xi)) + 0.00001
      ENDDO
      string = char_i(1)//char_i(2)//char_i(3)

      END FUNCTION integer_to_string
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc 
      END MODULE timer
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc
