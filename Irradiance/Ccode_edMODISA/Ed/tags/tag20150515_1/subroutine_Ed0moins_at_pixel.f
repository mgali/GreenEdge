       SUBROUTINE ed0moins(jday,rtime,lat,lon,o3,tcl,cf,ED_inst,ED_LUT
     &,thetas)
!    PROGRAM get_ed0moins_at_pixel

   
!   Define parameters
      parameter (nwl=83, nthetas=19,   nozone=8, ntaucl=8)

!    Declare variables:
      integer jday   
      real*4 lat, lon, rtime, o3, tcl, cf
      character clat*10, clon*10, co3*10, ctcl*10,ccf*10
      character cjday*3, ctime*10
      real*4  ed_cloud(nwl), ed_clear(nwl), thetas

!   Declare input Look-up-tables
      real*4 ED_LUT(nwl, nthetas, nozone, ntaucl) ! LUT for Ed0-(wl, thetas, o3, taucl)
!    Outputs:
      real*4 ED_inst(nwl), PAR, UV 

      character LUTPATH*50, fn*150
!      character LUTPATH*500, fn*150

!   Read user inputs
!    CALL GETARG(1,cjday)   
!    CALL GETARG(2,ctime)
!    CALL GETARG(3,clat)
!    CALL GETARG(4,clon)
!    CALL GETARG(5,co3)
!    CALL GETARG(6,ctcl)
!    CALL GETARG(7,ccf)

!   Convert Character to Integer or Real
!    READ(cjday,FMT='(I10)') jday
!    READ(ctime,FMT='(F10.2)') rtime
!    READ(clat,FMT='(F10.2)') lat
!    READ(clon,FMT='(F10.2)') lon
!    READ(co3,FMT='(F10.2)') o3
!    READ(ctcl,FMT='(F10.2)') tcl
!    READ(ccf,FMT='(F10.2)') cf

!    LUTPATH='/Data/PP_model/data/SBDART/'

   !write(*,*) 'Get Ed0- spectrum'
!    Spectral irradiance Look-up-table
!    CALL read_ed0moins_LUT(LUTPATH, ED_LUT)

!   Calculate the Sun Zenith Angle for the given location and date/time

      CALL sunpos(jday, rtime,
     &                  lat, lon, thetas, azim)

      if (thetas.lt.90.) then
! Get the Ed spectra for the given cloud optical thickness
      CALL interpol_ed0moins(ED_LUT,
     &                  thetas, o3, tcl, ed_cloud)
! Get the Ed spectra for clear sky
      CALL interpol_ed0moins(ED_LUT,
     &                  thetas, o3, 0.0, ed_clear)

!   Calcualte the Ed spectra for the given cloud fraction
      do i = 1,nwl
          ED_inst(i) = ed_cloud(i)*cf + ed_clear(i)*(1.-cf)
!        the Basic unit of Ed is microEinstein/s/m2/nm
      enddo
      else
        do i = 1,nwl
          ED_inst(i)=0.0
        enddo
      endif

! integrate the PAR and UV
      PAR = 0
      UV = 0
      do i=23,nwl-1
          h = (ED_inst(i)+ED_inst(i+1))/2.   
          PAR = PAR + (h)*5
      enddo
      do i=1,22
          h = (ED_inst(i)+ED_inst(i+1))/2.  
        UV = UV + (h)*5
      enddo

 
!    write(*,FMT='(89F)') thetas,o3, tcl, cf, PAR, UV, 
!      &(ED_inst(k), k=1,nwl) 

        return
      END
   
C -----------------------------------------------------------------
C Lecture de la LUT
C
      SUBROUTINE read_ed0moins_LUT(LUTPATH,LUT)
      
      parameter (nwl=83, nthetas=19,   nozone=8, ntaucl=8)

      real*4 LUT(nwl, nthetas, nozone, ntaucl)
      character LUTPATH*500, fn*100
      integer i,j,k,l

!       fn=trim(LUTPATH)//'Ed0moins_LUT.dat'

      OPEN(10, file=LUTPATH,status="old",
     &form="formatted",access="sequential")
c     &       form="unformatted",access="direct")
        do i=1,nthetas
          do j=1,nozone
            do k=1,ntaucl
              read(10, *) (LUT(l,i,j,k), l=1,nwl)
            enddo
          enddo
        enddo
        CLOSE(10) 
      END


C -----------------------------------------------------------------
C Interpolation de Ed
C
        SUBROUTINE interpol_ed0moins(LUT, thetas, ozone, taucl, ed)

        parameter (nwl=83, nthetas=19,   nozone=8, ntaucl=8)

        real*4  LUT(nwl, nthetas, nozone, ntaucl)
        real*4  thetas, ozone, taucl, ed(nwl)
        real*4   ed_tmp3(nwl,2,2), ed_tmp2(nwl,2) 
        real*4  xthetas(nthetas), xozone(nozone), xtaucl(ntaucl) 
        integer ithetas, iozone, itaucl
        real*4 rthetas, rozone, rtaucl
        integer i,j,k,l
        integer zthetas, zozone, ztaucl


         data xthetas /0.0, 5.0, 10.0, 15.0,
     &              20.0, 25.0, 30.0, 35.0,
     &          40.0, 45.0, 50.0, 55.0,
     &          60.0, 65.0, 70.0, 75.0,
     &          80.0, 85.0, 90.0/
         data xozone /200.0, 250.0, 300.0, 350.0,
     &              400.0, 450.0, 500.0, 550.0/
        data xtaucl /0.0, 1.0, 2.0, 4.0,
     &              8.0, 16.0, 32.0, 64.0/

!   If input values overflow the maximum value used to generate the LUT then reset it to its maximum values
      if (thetas.GE.90.) then
          thetas = 89.99
      endif
      if (ozone.GE.550.0) then
          ozone = 549.99
      endif
      if (taucl.GE.64.0) then
          taucl = 63.99
      endif

!    Get the indices for thetaS
         do i=1,nthetas-1
           if ((thetas.GE.xthetas(i)).AND.(thetas.LT.xthetas(i+1))) then
               ithetas=i
           endif
         enddo
         if (thetas.lt.xthetas(1)) then 
            ithetas=1
            rthetas=0.
         else
        rthetas=(thetas-xthetas(ithetas))/
     &           (xthetas(ithetas+1)-xthetas(ithetas))
         endif

!    Get the indices for Ozone
         do i=1,nozone-1
         if ((ozone.GE.xozone(i)).AND.(ozone.LT.xozone(i+1))) then
            iozone=i
         endif
         enddo
         if (ozone.lt.xozone(1)) then 
         iozone=1
         rozone=0.
         else
         rozone=(ozone-xozone(iozone))/(xozone(iozone+1)-xozone(iozone))
         endif
!    Get the indices for Cloud Optical Thickness
         do i=1,ntaucl-1
         if ((taucl.GE.xtaucl(i)).AND.(taucl.LT.xtaucl(i+1))) then
            itaucl=i
         endif
         enddo
         if (taucl.lt.xtaucl(1)) then 
         itaucl=1
         rtaucl=0.
         else
         rtaucl=(taucl-xtaucl(itaucl))/(xtaucl(itaucl+1)-xtaucl(itaucl))
         endif
! Start interpolation
         do i=1,2
      zthetas=ithetas+(i-1)
            do j=1,2
         zozone=iozone+(j-1)
         do l=1,nwl
                        ed_tmp3(l,i,j)=(1.-rtaucl)
     &                                *LUT(l,zthetas, zozone, itaucl)
     &                       +rtaucl*LUT(l,zthetas, zozone, itaucl+1)
         enddo
            enddo
         enddo
         do i=1,2
          do l=1,nwl
                ed_tmp2(l,i)=(1.-rozone)*ed_tmp3(l,i,1)
     &                          +rozone*ed_tmp3(l,i,2)
          enddo
          enddo
          do l=1,nwl
          ed(l)=(1.-rthetas)*ed_tmp2(l,1)
     &             +rthetas*ed_tmp2(l,2)
      if (ed(l).gt.10000) then ! data overflow!!!
        ed(l)=0.0
      endif
      enddo
      return
      END


!   PROGRAM SUNPOS calculates the position of the sun as azimuth, al-  *
!*   titude, declination. It also lets you calculate the "mass" of the  *
!*   atmosphere for a given solar altitud angle, and enables to calcu-  *
!*   late the attenuation of the solar beam due to the specific solar   *
!*   altitude for a given atmospheric transmittivity.                   *
!*** currently there is not every feature implemented (no transmittance etc...***

      SUBROUTINE sunpos(jday, hour, ycoor, xcoor, szendeg, sazideg)

*
* common header
*
      REAL*4        lsn         ! local solar noon (standard time)
      REAL*4        xcoor       ! X-coordinate in real numbers
      REAL*4        ycoor       ! Y-coordinate in real numbers
      INTEGER*2     elev        ! station's elevation (in meters)
      INTEGER*2     ltm         ! local time meridian (SET to 0 by SBE for GMT TIME)
      REAL*4       hour   ! hour in decimal, added by sbe 
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

