	subroutine getinp(P,CN,A,jstype,D,pL,Y,ek,cfact,pfact,soilty,
     1                  ieroty,dp,om)
C---------------------------------------------------------------
C 	Input parameters
C---------------------------------------------------------------
	implicit double precision (a-h, o-z)
	common/rain/rfix,rti(5000),rfi(5000),rcum(5000,2),ref(5000,2),ncum
	character*20 soilty

c----------------------------------
c read rainfall amount, watershed characteristics
c   P = rainfall amount in mm
c   CN = SCS curve number
c   A = watershed area (ha)
c   jstype = storm type (1=I, 2=II, 3=III, 4=Ia, 5=user curve)
c   D = storm duration (h)
c   pL = maximum flow path length (m)
c   Y = slope (%)
c----------------------------------
	read(1,*)P,CN,A,jstype,D,pL,Y
	read(1,'(A)')
c---------------------------------
c read inputs for soil erosion calculations
c   soilty = soil type (Character)
c   ek = soil erodibility
c   cfact= C factor
c   pfact= P factor
c------------------------------------
	read(1,'(A)')soilty
	read(1,*) ek, cfact, pfact, dp
c      convert dp to um
	 dp=dp*10000.d0
c----------------------------------
c-  ieroty = select method to estimate storm erosion
c-           Selections:
c-         0 or not present = Foster's method for R-factor
c-         1 = Using Williams R-factor
c-         2 = Using R-factor from GLEAMS with daily rainfall
      read(1,*,END=22) ieroty
	if ((ieroty.lt.3).and.(ieroty.ge.0)) go to 24
   22  ieroty=1
   24  continue
c------------------------------------------------

c--------------------------------------
c om = % organic matter, read if ek <0
	   om = 2.0d0
       if (ek.lt.0.d0) then
	   read(1,*,end=32)om
	 endif
c--------------------------------------

  32  continue
      if(jstype.gt.4) then
c-- For user defined case, read "tmid" first and then 24-h P/P24 curve
c-- "tmid" is stored in first position of the rcum(i,j) array
         read(1,*,end=40)rcum(1,1)
		 print*,rcum(1,1)
		 if (rcum(1,1).eq.0) then 
             print*,'ERROR: the first value must be tmid(h), followed'
     &      ,'in the next lines by the cumulative rainfall  that must'	 
     &      ,'begin with (0,0) and end with (24,1)'
             STOP
		 endif			 
         rcum(1,2)=0.5d0
         do 35 i=2,5000
              read(1,*,end=40)(rcum(i,j),j=1,2)
			  print*,i,(rcum(i,j),j=1,2)
			  if(rcum(i,1).eq.24) goto 40
35       continue
40       ncum=i
         if ((rcum(ncum,1).ne.24).or.(rcum(ncum,2).ne.1)) then
          print*,'ERROR: the cumulative rainfall must begin with ', 
     &      '(0,0) and end with (24,1)',i
          STOP
         endif
       endif
	return
	end

	subroutine results(P,CN,Q,A,tc,xIa,jstype,D,pL,Y,qp,tp,qdepth,
     1                   ieroty)
C---------------------------------------------------------------
C	Output summary of hydrology results nicely
C---------------------------------------------------------------
	implicit double precision (a-h, o-z)
	character*4 stype(5)
	data stype/'I','IA','II','III','user'/

	write(2,*)' '
	write(2,*)' HYDROGRAPH CALCULATION FOR WATERSHED-SCS METHOD'
	write(2,*)' '
	write(2,*)'INPUTS'
	write(2,*)'------'
	write(2,100)P
	write(2,200)stype(jstype)
	write(2,250)D
	write(2,300)CN
	write(2,400)A
	write(2,500)pL
	write(2,600)Y*100.d0
	write(2,610)ieroty
	write(2,*)' '
	write(2,*)'OUTPUTS'
	write(2,*)'-------'
	write(2,700)Q, Q*A*10.d0
	write(2,800)xIa
	write(2,900)tc,tc*60.d0
	write(2,*)' '

100	format('Storm Rainfall=',f8.2,' mm')
200	format('SCS storm type= ',a4)
250	format('Storm duration=',f6.1,' h')
300	format('SCS Curve number=',f6.1)
400	format('Watershed area=',f8.2,' ha')
500	format('Maximum flow path length=',f8.2,' m')
600	format('Average slope of flow path=',f8.2,' %')
610   format('MUSLE type=',i3,' where:',/,
     1    2x,'1=Williams, 2=GLEAMS, 3=Foster  (See Manual)')
700	format('Runoff volume=',f8.2,' mm=',f8.2,' m3')
800	format('Initial Abstraction=',f8.2,' mm')
900	format('Concentration time=',f8.2,' h= ',f8.2,' min')

	return
	end

      subroutine vfsout(dp,ieroty,sconc,sconc1,sconc2,A,pL,qp,tp,
     1                  tc,D,ti,nhyet,nhyd)
c--------------------------------
c	Output for VFSMOD input files
c--------------------------------
	implicit double precision (a-h,o-z)
	common/hydgph/u(5000,2),qh(5000,2)
	common/rain/rfix,rti(5000),rfi(5000),rcum(5000,2),ref(5000,2),ncum

c--------------------------------
c Output of VFSMOD input file: *.isd
c--------------------------------
	npart=7
	coarse=1.0d0
	if (ieroty.eq.1) then
	   ci= sconc1/1000.d0
	elseif (ieroty.eq.2) then
	   ci= sconc2/1000.d0
	elseif (ieroty.eq.3) then
	   ci= sconc/1000.d0
	else
	   ci= sconc1/1000.d0
	endif
c--rmc  05/08/03 when runoff is small, sediment concentration by sediment
c----- yields methods that do not consider runoff in calculation (Foster's
c----- ,CREAMS) can be very large. Override user selection of the method
c----- and slect Williams' that considers runoff and typically avoids this
c----- problem. Issue warning.
	if(ci.ge.0.25d0) then
         ci=sconc1/1000.d0
         write(*,160)
         write(*,*)'WARNING: small runoff in this case produces large',
     1    ' sediment concentration with the sediment yield method #', 
     2    ieroty,'selected. Using Williams method instead--see manual'  
         write(*,160)
         write(10,*)
         write(10,160)
         write(10,*)'WARNING: small runoff in this case produces large',
     1    ' sediment concentration with the sediment yield method #', 
     2    ieroty,'selected. Using Williams method instead--see manual'  
         write(10,160)
      endif
	por=0.434d0
	write (15,101) Npart,coarse,ci,por
	dpp=dp/10000.d0
	sg=2.65d0
	write (15,102)dpp,sg

c--------------------------------
c Output of VFSMOD runoff hydrograph: *.iro
c--------------------------------
	swidth=A*10000.d0/pL
	slength=pL
	write (12,103) swidth,slength
	nbcroff=nhyd
	bcropeak=qp
      nstep1=100
      if(nhyd.le.nstep1)then
            write (12,104)nbcroff+1,bcropeak
            write(2,*)' '
            write(2,250)ti
            do 20 ii=1, nbcroff-1
               tt=qh(ii,1)*3600.d0
               if (ii.eq.1) then
                  write(12,105)tt,qh(ii,2)
                 else 
                  write(12,106)tt,qh(ii,2)
               endif
20          continue
         else
            nwrite1=nhyd/nstep1+1
            write(12,104)nhyd/nwrite1+2,bcropeak
            do 29 ii=1,nhyd-1
               tt=qh(ii,1)*3600.d0
               if (ii.eq.1) then
                  write(12,105)tt,qh(ii,2)
                else
                  do 25 k=1,nstep1
                        if(ii.eq.k*nwrite1) then
                              write(12,106)tt,qh(ii,2)
c                              write(*,'(2i4,2e12.5)')tt,qh(ii,2)
                        endif
25                continue
               endif
29          continue            
      endif
c---write 0 entry after last step
      tend1=qh(nhyd,1)*3600.d0
      write(12,106)tend1,qh(nhyd,2)      
      write(12,107)tend1+300.d0,0.d0

c----------------------------------------
c Output VFSMOD rainfall hyetograph: *.irn
c----------------------------------------
      nstep2=100
      if(nhyet.le.nstep2)then
            write(14,201)nhyet+1,rfix
            do 31 ii=1,nhyet-1
               tt=rti(ii)*3600.d0
               if (ii.eq.1) then
                  write(14,203) tt,rfi(ii)
                else
                  write(14,204) tt,rfi(ii)
               endif
31          continue
         else
            nwrite2=nhyet/nstep2+1
            write(14,201)nhyet/nwrite2+2,rfix
            do 33 ii=1,nhyet-1
               tt=rti(ii)*3600.d0
               if (ii.eq.1) then
                  write(14,203) tt,rfi(ii)
                else
                  do 32 k=1,nstep2
                        if(ii.eq.k*nwrite2) then
                              write(14,204) tt,rfi(ii)
c                              write(*,'(2i4,2e12.5)')ii,k,tt,rfi(ii)
                        endif
32                continue
               endif
33          continue            
      endif
c---write 0 entry after last step
      tend2=rti(nhyet)*3600.d0
      write(14,204)tend2,rfi(nhyet)      
      write(14,205)tend2+300.d0,0.d0

c----------
c Output message at end of program -----------------
c----------
	WRITE(*,*)
	WRITE(*,*)'...FINISHED...','UH v3.0.1 2/2012'
	WRITE(*,*)

c-------------------
c  format statements
c-------------------
101	format (2x,i4,2x,f8.1,2x,f11.4,2x,f7.4,8x,
     1   'Npart, Coarse, Ci(g/cm3), Por')
102	format (2x,f10.7,2x,f7.1,21x,'Dp(cm), SG(g/cm3)')
103	format (2x,f7.1,2x,f7.1,21x,'Swidth(m), Slength(m)')
104	format (2x,i4,2x,e12.5,19x,'nbcroff, bcropeak (m3/s)')
105	format (2x,e12.5,2x,e12.5,10x,' time(s), ro(m3/s)')
106	format (2x,e12.5,2x,e12.5)
107	format (2x,e12.5,2x,e12.5,/,30('-'))
160	FORMAT(72('-'))
201	format(i4,2x,e12.5,20x,' NRAIN, RPEAK(m/s)')
203	format(2x,e12.5,3x,e12.5,10x,'time(s), rainfall rate (m/s)')
204	format(2x,e12.5,3x,e12.5)
205	format(2x,e12.5,3x,e12.5,/,30('-'))
250	format('Time to ponding=',f8.3,' h')
260	format('Duration of rainfall excess=',f8.3,' h')
270	format('Time to peak after shifting=',f8.3,' h')
280	format('Time correction to match hyetograph=',f8.3,' h')

	return
	end
