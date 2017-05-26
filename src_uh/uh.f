	program uh
C---------------------------------------------------------------
c      version 3.0.2, Last Modified: See Modifications below
C      WRITTEN FOR: ASAE'99 Toronto paper, March 8, 2002 
C      Written by: R. Munoz-Carpena (rmc)   &   J. E. Parsons, BAE (jep)
C                  University of Florida        BAE, NC State University
C                  Gainesville, FL 32611        Raleigh, NC 27695-7625(USA)
C                  e-mail: carpena@ufl.edu      
C---------------------------------------------------------------
C Program to create input files for VFSmod, based on NRCS-TR55 and Haan
C et al, 1996, with additional work done on coefficients for unit peak
C flow calculation.
c
c    Date      Modification                                Initials
c   -------    ------------------------------                ------ 
c   2/17/99    Check for 0.1<Ia/P<0.5                           rmc
c   2/18/99    Added hyetograph output for 6 h storm            jep
c   2/18/99    Modify File Inputs for Erosion                   jep
c   2/20/99    Roughed in MUSLE                                 jep
c   3/01/99    Checked erosion parameters and units             rmc
c   3/02/99    Additional work on Musle - units close           jep
c   3/03/99    Added hyetographs for storm types I & IA         rmc
c   3/05/99    output irs file for VFSMOD                       jep
c   3/06/99    Input/Output files as in VFSMOD                  jep
c   3/10/99    Checked Input/Output files as in VFSMOD          rmc
c   3/10/99    Cleanup - created hydrograph.f for 
c              hydrograph subroutines, created io.f for
c              input and output related processing              jep
c   3/28/99    Erosion part: fixes in I30 calculation 
c              after Chow and checked for consistency in
c              units, clean up; Hydro: added delay time         rmc
c   8/27/99    Added option to select different methods 
c              for applying MUSLE, default is Foster, 
c              2=Williams, 3=GLEAMS
c   10/01/99   Fixed array so that storm duration (D)        
c              can now be up to 24h                             rmc
c   10/26/99   implemented the project file concept as in vfsm  jep
c    3/09/00   Version changed to 0.9, general program cleanup  rmc
c   16/06/00   Version changed to 1.0, erosion output organized rmc
c   16/03/02   Version changed to 1.06 to couple with VFSMOD, 
c              author affiliation changed                       rmc
c    4/18/03   Fixed K - computed if we enter -1, other use     jep
c              entered value, also fixed dp output format
c    4/19/03   dp now being read in                             jep
c    4/20/03   Runoff calculation for low CN revised            rmc
c    5/01/03   Added chacked for small runoff case to switch    rmc
c              to Williams sediment calculation that includes
c              runoff.
c   11/10/03   Reordered Erosion ieroty 1=Williams, 2=Gleams
c                3=Foster to coincide with changes in Shell     jep
c   11/13/03   Fixed coef. on Type Ia - did not add new
c                hyet curves
c   01/10/05   Added changes suggested by U. of Guelph group    rmc 
c              v2.4.1
c   09/15/11   Rewritten hydrograph calculation using convolution
c              of excess rain steps, v3.0.0                     rmc
c   02/15/12   Added user table for 24-h hyetograph, v3.0.1     rmc
c   02/15/14   Fixed user hyetograph option (jstype=5)v3.0.2    rmc
c
c---------------------------------------------------------------
c    Compiling for Win32 and Unix environments:
c       1. The i/o for these operating systems is different.
c       2. Change the comments in the finput.f program to reflect
c          your operating system.   3/9/00
c---------------------------------------------------------------
c common/hydgph:
c       rot(208), runoff time (units)
c       roq(208), runoff rate (m3/s)
c       u(208,2), unit hydrograph
c common/rain/:
c       rfix, maximum rain intensity (mm/h)
c       rti(200), rainfall time (hrs)
c       rfi(200), rainfall intensity (mm/h)
c       rcum(100,2), cumm rainfall (mm)
c       ref(100), excess rainfall intensity (mm/h)
c       ncum: number of steps if user hyetograph is read
c other:
c       nref=number of excess hyetograph steps
c       mref=number of unit hydrograph steps
c       nhyet=number of hyetograph steps
c       vol(m3), volro (mm)=runoff volume 
C---------------------------------------------------------------
	implicit double precision (a-h, o-z)
	character*20 soilty
	CHARACTER*75 LISFIL(6)

	common/hydgph/u(5000,2),qh(5000,2)
	common/rain/rfix,rti(5000),rfi(5000),rcum(5000,2),ref(5000,2),ncum

c------------------------------------------
c  get inputs and open files
c------------------------------------------
	call  FINPUT(LISFIL)
	call getinp(P,CN,A,jstype,D,pL,Y,ek,cfact,pfact,soilty,
     1            ieroty,dp,om)
c--------------------------
c Calculate runoff volume by SCS method
c--------------------------
	call runoff(P,CN,xIa,Q)
	volro=Q
c--------------------------
c Calculate concentration time by SCS method
c--------------------------
	call calctc(pL,CN,Y,tc)
c--------------------------
c Calculate peak flow and time by SCS-TR55 method
c--------------------------
	call q_peak(A,Q,xIa,P,tc,jstype,qp,tp)
c--------------------------
c Output hydrology results
c--------------------------
	call results(P,CN,Q,A,tc,xIa,jstype,D,pL,Y,qp,tp,qdepth,
     1             ieroty)
c--------------------------
c Calculate SCS-unit hydrograph
c--------------------------
	call unit_hyd(Q,A,qp,tp,D,tc,qp5,tp5,mref)
	vol=volro*A*10000.d0/1000.d0

c--------------------------
c Calculate storm hyetograph from SCS storm type
c--------------------------
 	call hyetgh(jstype,P,D,volro,qdepth,vol,qp,A,xIa,rtpeak,er,er1,
     1            erCoolm,ti,nref,tc,a1,b1,bigE,raimax30,nhyet)

c--------------------------
c Calculate storm hydrograph
c--------------------------
      call tab_hyd(Q,A,mref,nref,qp,nhyd)

c--------------------------
c do the modified usle to get erosion stuff
c--------------------------
      call musle(er,er1,erCoolm,ek,Y,pl,cfact,pfact,A,vol,tc,P,D,soilty,
     1           dp,sconc,sconc1,sconc2,om,a1,b1,bigE,raimax30,qp)

c------------------------------------
c write vfsmod compatible input files
c------------------------------------
      call vfsout(dp,ieroty,sconc,sconc1,sconc2,A,pL,qp,tp,tc,D,ti,
     1           nhyet,nhyd)

	close(1)
	close(2)

	stop
	end

