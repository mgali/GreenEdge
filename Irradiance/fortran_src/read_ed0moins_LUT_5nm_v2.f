
!	This routine reads returns Ed0- LUT which has 5 dimensions for 
!		1. Wavelength = 290 : 700 : 5 
!		2. ThetaS = 0 : 90 : 5
!		3. Ozone = 100 : 550 : 50
!		4. Cloud optical Thickness = 0 to 64 = c(0,1,2,4,8,16,32,64)
!               5. Surface Albedo = 0.05 : 0.9 : 0.15

	SUBROUTINE read_ed0moins_LUT_5nm_v2(PATH,LUT)
	
	parameter (nwl=83, nthetas=19,	nozone=10, ntaucl=8, nalb=7)

	real*4 LUT(nwl, nthetas, nozone, ntaucl, nalb)
	character PATH*50, fn*100
	integer i,j,k,l,m

	fn=trim(PATH)//'Ed0moins_LUT_5nm_v2.dat'

	OPEN(10, file=TRIM(fn),status="old",
     &       form="formatted",access="sequential")
	do i=1,nthetas
		do j=1,nozone
			do k=1,ntaucl
                           do m=1,nalb
				read(10, *) (LUT(l,i,j,k,m), l=1,nwl)
                           enddo
			enddo
		enddo
	enddo
	CLOSE(10) 
	END
