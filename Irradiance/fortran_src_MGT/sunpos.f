        SUBROUTINE sunpos(jday, hour, ycoor, xcoor, szendeg, sazideg)



*   PROGRAM SUNPOS calculates the position of the sun as azimuth, al-  *
*   titude, declination. It also lets you calculate the "mass" of the  *
*   atmosphere for a given solar altitud angle, and enables to calcu-  *
*   late the attenuation of the solar beam due to the specific solar   *
*   altitude for a given atmospheric transmittivity.                   *
*** currently there is not every feature implemented (no transmittance etc...***


*
* common header
*
      REAL*4        lsn         ! local solar noon (standard time)
      REAL*4        xcoor       ! X-coordinate in real numbers
      REAL*4        ycoor       ! Y-coordinate in real numbers
      INTEGER*2     elev        ! station's elevation (in meters)
      INTEGER*2     ltm         ! local time meridian (SET to 0 by SBE for GMT TIME)
      REAL*4        hour        ! hour in decimal, added by sbe 
*
c      COMMON /header/lsn,xcoor,ycoor,ltm,elev,hour
*
* common time
*
      INTEGER*2     mth         ! month
      INTEGER*2     day         ! day
      INTEGER*2     jday        ! Julian Day
      INTEGER*2     hr          ! hour (local standard time)
      INTEGER*2     min         ! minute (local standard time)
      INTEGER*2     hangle      ! hourangle (local standard time)   
*
*
c      COMMON /time/mth,day,jday,hr,hangle
*
* common solar
*
      REAL*4        pi          ! PI
      REAL*4        d2r         ! degree to radians conversion
      REAL*4        r2d         ! radians to degree conversion
      REAL*4        mass        ! mass of the atmosphere for given solaralt
      REAL*4        decrad      ! declination of the sun
      REAL*4        latrad      ! latitude of site in radians
      REAL*4        saltrad     ! solar altitude in radians
      REAL*4        saltdeg     ! solar altitude in degrees
      REAL*4        szenrad     ! solar zenith angle in radians
      REAL*4        szendeg     ! solar zenith angle in degrees
      REAL*4        sazirad     ! solar azimuth angle in radians
      REAL*4        sazideg     ! solar azimuth angle in degrees
      REAL*4        tbbase      ! base function of beam related transmittance
      REAL*4        tbfac       ! factor to multiply beam transmittance function
      REAL*4        basemax     ! maximum base function for vertical sun pos.
*
c      COMMON /solar/pi,d2r,r2d,mass,decrad,latrad,saltrad,saltdeg,
c     +              szenrad,szendeg,sazirad,sazideg,tbbase,
c     +              tbfac,basemax
*
* local variables
*
      REAL*4        ha          ! hour (floating from hr & min)
      REAL*4        harad       ! hourangle in radians
      REAL*4        decdeg      ! declination in degrees

* -------------------------------------------------------------------- *
*
      ltm = 0
      hr = int(hour)
      min = int( (hour - real(hr)) * 60.)
*
      pi      = 3.14159265358979
      d2r     = pi/180.0
      r2d     = 1/d2r
      lsn     = 12.0+((FLOAT(ltm)-xcoor)/15.0)
      latrad  = ycoor*d2r
      decrad  = 23.45*d2r*sin(d2r*360.*(284.+FLOAT(jday))/365.)
      decdeg  = decrad*r2d
*
      ha      = FLOAT(hr)+(FLOAT(min)/60.0) !convert to floating
      hangle  = (lsn-ha)*60.0               !solrad is given for the hour preceeding the time given
      harad   = hangle*0.0043633            !convert hangle (in minutes) into radians
*                                           !This is the same as multiplying the #hrs by 15 (deg./hr),
*                                           !and convert d2r (*0.017453292)
*
      saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad)
     +               *cos(harad)))
      saltdeg = saltrad * r2d
      sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
      sazideg = sazirad * r2d
      IF (saltdeg.LT.0.0 .OR. saltrad.GT.180.0) THEN    ! sun is below horizon
         saltdeg = 0.0
         saltrad = 0.0
         szendeg = 90.0
         szenrad = 90.0*d2r
         mass    = 1229**.5             ! if solaralt=0 -> sin(0)=0
      ELSE
         szendeg = 90-saltdeg
         szenrad = szendeg*d2r
         mass = (1229.0+(614.0*sin(saltrad))**2)**.5-(614*sin(saltrad))
      ENDIF
*
      tbbase     = exp(-.65*mass)+exp(-.09*mass)

        return

      END


