!       /opt/intel/bin/ifort get_ed0.f interpol_ed0LUT_5nm_v2.f 
!       read_ed0moins_LUT_5nm_v2.f read_ed0plus_LUT_5nm_v2.f sunpos.f -o
!       get_ed0


!       This programm returns the spectral irradiance a 5nm that are
!       stored in a LUT generated using SBDART
!
!       It takes as input the following arguments:
!          #1 = 1 for Ed0 above surface 0 for Ed below sea surface
!          #2 = Day of Year
!          #3 = time in UTC from 0 to 24.0 (need decimal!!)
!          #4 = latitude in decimal
!          #5 = longitude in decimal
!          #6 = cloud optical thickness (0 to 64.0)
!          #7 = ozone concentration (100 to 550 DU)
!          #8 = cloud fraction (0 to 1.0)
!          #9 = surface albedo (0 to 1.0)

!	OUTPUT:
!		The input parameters follow by the integrated PAR and UV and 
!		then the Ed0 spectra from 290 nm to 700 nm at 5 nm resolution
!

!
!       AUTHOR: Simon Bélanger ,UQAR
!
!       Date of creation: 21 April 2017
!

        PROGRAM get_ed0

!       Define parameters
        parameter (nwl=83, nthetas=19,  nozone=10, ntaucl=8, nalb=7)

!       Declare user's input variables
        character clat*10, clon*10, cjday*3,ctime*10,above*1
        character ctcl*10, co3*10, ccf*10, csurfalb*10

!       Declare internal variables
        real*4 lat, lon, rtime
        real*4 thetas, tcl, o3, cf, surfalb
        real*4 wl(nwl), ed(nwl)
        real*4 ed_cloud(nwl), ed_clear(nwl)
        character pathSBDART*50

!       Declare input Look-up-tables
        real*4 ED_LUT(nwl, nthetas, nozone, ntaucl, nalb) ! LUT for Ed0+, surfalib(wl, thetas, o3, taucl, surfalb)

        CALL GETARG(1,above)
        CALL GETARG(2,cjday)
        CALL GETARG(3,ctime)
        CALL GETARG(4,clat)
        CALL GETARG(5,clon)
        CALL GETARG(6,ctcl)
        CALL GETARG(7,co3)
        CALL GETARG(8,ccf)
        CALL GETARG(9,csurfalb)

        if (csurfalb  == '' ) then
            print *, 'get_ed0 needs 9 arguments'
            print *, '#1 = 1 for Ed0 above surface 
     &0 for Ed below sea surface'
            print *, '#2 = Day of Year'
            print *, '#3 = time in UTC from 0 to 24.0 (need decimal!!)'
            print *, '#4 = latitude in decimal'
            print *, '#5 = longitude in decimal'
            print *, '#6 = cloud optical thickness (0 to 64.0)'
            print *, '#7 = ozone concentration (100 to 550 DU)'
            print *, '#8 = cloud fraction (0 to 1.0)'
            print *, '#9 = surface albedo (0 to 1.0)'
            print *, 'IMPORTANT TO PUT THE DECIMALS!!!!'
            print *, 'Example: ./get_ed0 1 172 18. 70. -140. 
     &3. 330. 1. 0.5'
            stop
        endif


        READ(cjday,FMT='(I10)') jday
        READ(ctime,FMT='(F10.2)') rtime
        READ(clat,FMT='(F10.2)') lat
        READ(clon,FMT='(F10.2)') lon
        READ(ctcl,FMT='(F10.2)') tcl
        READ(co3,FMT='(F10.2)') o3
        READ(ccf,FMT='(F10.2)') cf
        READ(csurfalb,FMT='(F10.2)') surfalb
        pathSBDART='/Users/martigalitapias/Documents/SBDART/'

!       Spectral irradiance Look-up-table
        if (above.eq.'1') then
                CALL read_ed0plus_LUT_5nm_v2(pathSBDART, ED_LUT)
        else
                CALL read_ed0moins_LUT_5nm_v2(pathSBDART, ED_LUT)
        endif

!       Define the wavelenght vector       
        do i=1,nwl
                wl(i) = (i-1)*5 + 290.
        enddo


!       Compute the sun position
        CALL sunpos(jday, rtime, lat, lon, thetas, azim)


        CALL interpol_ed0LUT_5nm_v2(ED_LUT,
     &       thetas, o3, tcl, surfalb, ed_cloud)
        CALL interpol_ed0LUT_5nm_v2(ED_LUT,
     &       thetas, o3, 0.0, surfalb, ed_clear)

         do i = 1,nwl
                ed(i) = ed_cloud(i)*cf + ed_clear(i)*(1.-cf)
!        the  unit of Ed is microEinstein/s/m2/nm
         enddo

!              Integrate Ed over UV and PAR
        PAR = 0
        UV = 0
        do i=21,nwl-1
                h = (ed(i)+ed(i+1))/2.
                PAR = PAR + (h)*5
        enddo
        do i=1,20
                h = (ed(i)+ed(i+1))/2.
                UV = UV + (h)*5
        enddo

!       Write the results
        write(*,FMT='(1A, 1I3.0, 92F10.4)')
     &above,jday,rtime,lat, lon, o3, cf, tcl, 
     &surfalb, PAR, UV,
     &(ed(k), k=1,nwl)


        END

