
! --------------------------------------------------------
! MODULE NANE: 
!     time manager model
!
! PURPOSE :
!     To provide some basic operations for time stamp
!
! Initial author: Hua Yuan, /04/2014/
! --------------------------------------------------------

MODULE timemanager

   use precision
   use MOD_TimeInvariants, only: lons
   implicit none

   type :: timestamp
      integer :: year, day, sec
   end type timestamp
  
   INTERFACE ASSIGNMENT (=) 
      module procedure assignidate
      module procedure assigntstamp
   END INTERFACE

   INTERFACE OPERATOR (+)
      module procedure addsec
   END INTERFACE

   INTERFACE OPERATOR (-)
      module procedure subtstamp
   END INTERFACE

   INTERFACE OPERATOR (<=)
      module procedure lessequal
   END INTERFACE
   
   INTERFACE OPERATOR (<)
      module procedure lessthan
   END INTERFACE

   INTERFACE OPERATOR (==)
      module procedure isnull
   END INTERFACE

   INTERFACE calendarday
      module procedure calendarday_date
      module procedure calendarday_stamp
   END INTERFACE
   
   logical, save :: isgreenwich
 
CONTAINS
   
   SUBROUTINE initimetype(greenwich)
      
      implicit none
      logical, intent(in) :: greenwich

      isgreenwich = greenwich
      
   END SUBROUTINE initimetype
 
   SUBROUTINE assignidate(tstamp, idate)
      
      implicit none
      type(timestamp), intent(inout) :: tstamp
      integer,         intent(in)    :: idate(3)
      
      tstamp%year = idate(1)
      tstamp%day  = idate(2)
      tstamp%sec  = idate(3)

   END SUBROUTINE assignidate
 
   SUBROUTINE assigntstamp(tstamp1, tstamp2)
      
      implicit none
      type(timestamp), intent(out) :: tstamp1
      type(timestamp), intent(in)  :: tstamp2
      
      tstamp1%year = tstamp2%year
      tstamp1%day  = tstamp2%day
      tstamp1%sec  = tstamp2%sec

   END SUBROUTINE assigntstamp
   
   FUNCTION addsec(tstamp, sec)

      implicit none
      type(timestamp), intent(in) :: tstamp
      integer,         intent(in) :: sec
      type(timestamp) :: addsec
      integer         :: maxday  

      addsec = tstamp
      addsec%sec = addsec%sec + sec
      if (addsec%sec > 86400) then
         addsec%sec = addsec%sec - 86400
         if( isleapyear(addsec%year) ) then
            maxday = 366
         else
            maxday = 365
         endif
         addsec%day = addsec%day + 1
         if(addsec%day > maxday) then
            addsec%year = addsec%year + 1
            addsec%day = 1
         endif
      endif
      return

   END FUNCTION addsec

   FUNCTION subtstamp(tstamp1, tstamp2)

      implicit none
      type(timestamp), intent(in) :: tstamp1
      type(timestamp), intent(in) :: tstamp2
      integer :: subtstamp
     
      subtstamp = tstamp1%sec - tstamp2%sec
      if (subtstamp < 0) then
         subtstamp = subtstamp + 86400
      end if
      return

   END FUNCTION subtstamp
  
   LOGICAL FUNCTION lessequal(tstamp1, tstamp2)

      implicit none
      type(timestamp), intent(in) :: tstamp1
      type(timestamp), intent(in) :: tstamp2

      integer(kind=4) :: ts1, ts2

      ts1 = tstamp1%year*1000 + tstamp1%day
      ts2 = tstamp2%year*1000 + tstamp2%day

      lessequal = .false.

      if (ts1 < ts2) lessequal = .true.
      
      if (ts1==ts2 .AND. tstamp1%sec<=tstamp2%sec) then
         lessequal = .true.
      end if
 
      return

   END FUNCTION lessequal
 
   LOGICAL FUNCTION lessthan(tstamp1, tstamp2)

      implicit none
      type(timestamp), intent(in) :: tstamp1
      type(timestamp), intent(in) :: tstamp2

      integer(kind=4) :: ts1, ts2
      
      ts1 = tstamp1%year*1000 + tstamp1%day
      ts2 = tstamp2%year*1000 + tstamp2%day

      lessthan = .false.

      if (ts1 < ts2) lessthan = .true.

      if (ts1==ts2 .AND. tstamp1%sec<tstamp2%sec) then
         lessthan = .true.
      end if

      return

   END FUNCTION lessthan
   
   LOGICAL FUNCTION isnull(tstamp, nullstr)

      implicit none
      type(timestamp), intent(in) :: tstamp
      character(4),    intent(in) :: nullstr

      if (tstamp%year < 0 .OR. tstamp%day < 0 .OR. tstamp%sec < 0) then
         isnull = .true.
      else
         isnull = .false.
      end if
      return
      
   END FUNCTION isnull

   LOGICAL FUNCTION isleapyear(year)

      implicit none
      integer, intent(in) :: year

      if( (mod(year,4)==0 .AND. mod(year,100)/=0) .OR. &
         mod(year,400)==0 ) then
         isleapyear = .true.
      else
         isleapyear = .false.
      end if
      return
   END FUNCTION isleapyear

   SUBROUTINE julian2monthday(year, day, month, mday)

      implicit none
      integer, intent(in)  :: year, day
      integer, intent(out) :: month, mday
      
      integer :: i, months(0:12)
      
      if ( isleapyear(year) ) then
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      end if

    ! calculate month and day values
      do i = 1, 12
         if (day .LE. months(i)) then
            month = i; exit
         end if
      end do
      mday = day - months(i-1)

   END SUBROUTINE julian2monthday
   
   LOGICAL FUNCTION isendofday(idate, sec)

      implicit none
      integer, intent(in) :: idate(3)
      real(r8),intent(in) :: sec

      type(timestamp) :: tstamp1
      type(timestamp) :: tstamp2

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      if (tstamp2%day /= tstamp1%day) then
         isendofday = .true.
      else
         isendofday = .false.
      end if
      return 

   END FUNCTION isendofday
   
   LOGICAL FUNCTION isendofmonth(idate, sec)

      implicit none
      integer, intent(in) :: idate(3)
      real(r8),intent(in) :: sec

      type(timestamp) :: tstamp1
      type(timestamp) :: tstamp2
      integer :: month1, month2, day

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      call julian2monthday(tstamp1%year, tstamp1%day, month1, day)
      call julian2monthday(tstamp2%year, tstamp2%day, month2, day)

      if (month1 /= month2) then
         isendofmonth = .true.
      else
         isendofmonth = .false.
      end if
      return 

   END FUNCTION isendofmonth
 
   LOGICAL FUNCTION isendofyear(idate, sec)

      implicit none
      integer, intent(in) :: idate(3)
      real(r8),intent(in) :: sec

      type(timestamp) :: tstamp1
      type(timestamp) :: tstamp2
      integer :: month1, month2, day

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      if (tstamp1%year /= tstamp2%year) then
         isendofyear = .true.
      else
         isendofyear = .false.
      end if
      return 

   END FUNCTION isendofyear
   
   SUBROUTINE adj2begin(idate)
      
      implicit none
      integer, intent(inout) :: idate(3)
      
      if (idate(3) == 86400) then
         idate(3) = 0
         idate(2) = idate(2) + 1
         if (isleapyear(idate(1)) .AND. idate(2)==367) then
            idate(1) = idate(1) + 1; idate(2) = 1
         end if
         if ( .NOT. isleapyear(idate(1)) .AND. idate(2)==366) then
            idate(1) = idate(1) + 1; idate(2) = 1
         end if
      end if

   END SUBROUTINE adj2begin
  
   SUBROUTINE adj2end(idate)
      
      implicit none
      integer, intent(inout) :: idate(3)
      
      if (idate(3) == 0) then
         idate(3) = 86400
         idate(2) = idate(2) - 1
         if (idate(2) == 0) then
            idate(1) = idate(1) - 1
            if ( isleapyear(idate(1)) ) then
               idate(2) = 366
            else
               idate(2) = 365
            end if
         end if
      end if

   END SUBROUTINE adj2end

   SUBROUTINE localtime2gmt(idate, long)

      implicit none
      integer, intent(inout) :: idate(3)
      real(r8),intent(in)    :: long
      
      integer  maxday                   
      real(r8) tdiff

      tdiff = long/15.*3600.
      idate(3) = idate(3) - int(tdiff)
      
      if (idate(3) < 0) then
         
         idate(3) = 86400 + idate(3)
         idate(2) = idate(2) - 1

         if (idate(2) < 1) then
            idate(1) = idate(1) - 1
            if ( isleapyear(idate(1)) ) then
               idate(2) = 366
            else
               idate(2) = 365
            endif
         endif
      endif
      
      if (idate(3) > 86400) then

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1
         
         if ( isleapyear(idate(1)) ) then
            maxday = 366
         else
            maxday = 365
         endif

         if(idate(2) > maxday) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif
      
   END SUBROUTINE localtime2gmt

   SUBROUTINE ticktime(deltim, idate)

      implicit none

      real(r8),INTENT(in)    :: deltim   
      integer, INTENT(inout) :: idate(3)
      integer maxday                   

      idate(3) = idate(3) + nint(deltim)
      if (idate(3) > 86400) then

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1
         
         if ( isleapyear(idate(1)) ) then
            maxday = 366
         else
            maxday = 365
         endif

         if(idate(2) > maxday) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif

   END SUBROUTINE ticktime
   
   real(r8) FUNCTION calendarday_date(date, long)

      implicit none
      integer, intent(in) :: date(3)
      real(r8),optional   :: long
      
      integer idate(3)
      real(r8) longitude

      idate(:) = date(:)
      
      if (.NOT. present(long)) then
         longitude = lons(1)
      else 
         longitude = long
      end if

      if ( .not. isgreenwich ) then
         call localtime2gmt(idate, longitude)
      end if
      
      calendarday_date = float(idate(2)) + float(idate(3))/86400.
      return 

   END FUNCTION calendarday_date
   
   real(r8) FUNCTION calendarday_stamp(stamp, long)

      implicit none
      type(timestamp), intent(in) :: stamp
      real(r8),        optional   :: long

      integer idate(3)
      real(r8) longitude
      
      idate(1) = stamp%year 
      idate(2) = stamp%day
      idate(3) = stamp%sec
      
      if (.NOT. present(long)) then
         longitude = lons(1)
      else 
         longitude = long
      end if
      
      if ( .not. isgreenwich ) then
         call localtime2gmt(idate, longitude)
      end if
      
      calendarday_stamp = float(idate(2)) + float(idate(3))/86400.
      return 

   END FUNCTION calendarday_stamp

END MODULE timemanager
