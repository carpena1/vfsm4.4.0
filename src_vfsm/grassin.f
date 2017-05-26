        SUBROUTINE GRASSIN(ICOARSE,COARSE,LISFIL,INARGS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     Get filter properties to be used in later from filename.igr             C
C     NOTE: all units in CGS system (cm,g,s), including modified Manning's n  C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
 	DIMENSION PARTC(6,3)
	CHARACTER*19 PARTLBL(6),PARTLB
	CHARACTER*75 LISFIL(13)
	DATA (PARTC(I,1), I=1,6),(PARTC(I,2), I=1,6),(PARTC(I,3), I=1,6)
c-----particle diameters (cm)---------------
     &	/.0002D0, .001D0, .003D0, .03D0, .02D0, .0029D0,
c-----fall velocities (cm/s)----------------
     & .0004D0, .0094D0, .0408D0 ,3.0625D0, 3.7431D0, .076D0,
c-----particle densities (g/cm3)------------
     &  2.6D0, 2.65D0, 1.80D0, 1.60D0, 2.65D0, 2.65D0/
	DATA (PARTLBL(I), I=1,6)/'    Clay','    Silt',' Sm.Aggregate',
     &  ' Lg.Aggregate','    Sand',' Silt (USDA)'/	

c---------Read inputs from grass.in file ----------------------

	IF(INARGS.EQ.1) THEN
	  WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(7)
	ENDIF
	READ(12,*)SS, VN, H, VN2,ICO

c---------Read inputs from sediment.in file ----------------------

	IF(INARGS.EQ.1) THEN
	  WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(11)
	ENDIF
	READ(16,*)NPART,COARSE,CI,POR
	
C--if the particle class is other than those defined above, calculate values--

	IF(NPART.LE.6)THEN
		DO 10 I=1,3
			PART(I)=PARTC(NPART,I)
10		CONTINUE
		PARTLB=PARTLBL(NPART)		
	  ELSE
	    READ(16,*)DP,SG
		PART(1)=DP
		PART(3)=SG
c----Barfield et al. (1981) to calculate fall velocity in GRASSF-----------
c----note: dp in mm--------------------------------------------------------
C		VC1=-0.34246727D0
C		VC2=0.98912185D0
C		VC3=1.1461280D0
C		CDP=DLOG(DP/10.D0)
C		PART(2)=10.D0**(VC1*CDP*CDP+VC2*CDP+VC3)
c----Fair and Geyer method (1954), based on Stokes-----------------------
c----note: all units in SI------------------------------------------------
		AGRAV=9.8066352D0
		VISCOS=0.000001003352832D0
		MAXITER=100
		DIA=PART(1)/100.D0
		FV=AGRAV*(SG-1.D0)*DIA**2.D0/(18.D0*VISCOS)
		C1=DIA*VISCOS
		REYN=FV*C1
		IF(REYN.LE.0.1D0) GO TO 70
		C2=SQRT(4.D0*AGRAV*(SG-1.D0)*DIA/3.D0)
		ITER=0
		EPS=1.D0
		FVOLD=10000000.D0
		DO 69 WHILE(EPS.GT.0.00001D0.AND.ITER.LT.MAXITER)
		    ITER=ITER+1
		    CD=24.D0/REYN+3.D0/SQRT(REYN)+0.34D0
		    FV=C2/SQRT(CD)
		    REYN=FV*C1
		    EPS=DABS(FV-FVOLD)
		    FVOLD=FV
69		CONTINUE
70		PART(2)=FV*100
		PARTLB='    (user) coarse'
		IF(DP.LT.0.0037D0)PARTLB='    (user) fine'
	ENDIF		

c-------If particle is fine(<37 microns) don't run the wedge part (icoarse=0)-
 
	ICOARSE=1
	IF(COARSE.LE.0.D0)ICOARSE=0	

c-------Output some input values, the rest are given in INPUTS subroutine ---

	COARSEP= COARSE*100.D0
	WRITE(13,*)'Sediment parameters'
	WRITE(13,*)'-----------------------------------------'
      	WRITE(13,200)'Sediment inflow concentration(Ci)=',CI,'g/cm3'
	WRITE(13,99)PARTLB
	WRITE(13,200)'     Particle size (diameter, dp)=',PART(1),'cm'
	WRITE(13,200)'      Particle fall velocity (Vf)=',PART(2),'cm/s'
	WRITE(13,200)'Particle weight density (gamma-s)=',PART(3),'g/cm3'
	WRITE(13,200)'% of particles with dp>0.0037 cm =',COARSEP,'%'
      	WRITE(13,200)'Porosity of deposited sediment(P)=',POR*100.D0,'%'
	WRITE(13,*)
	WRITE(14,101)
	WRITE(14,102)
	WRITE(14,107)

99	FORMAT('                   Particle class =',a19)
101   	FORMAT('   Time      qin      q1       Rs1       Vm1        ', 
     & 'df1       q2       Rs2       Vm2       df2      q3        Rs3',
     & '        Vm3       df3      qout')
102	FORMAT('    (s)  (cm3/s/cm)(cm3/s/cm)  (cm)     (cm/s)     (cm)',
     & '    (cm3/s/cm)  (cm)     (cm/s)     (cm)   (cm3/s/cm)  (cm)',
     & '      (cm/s)      (cm) (cm3/s/cm)')
107	FORMAT(147('-'))
200   	FORMAT(A35,F12.6,A11)

	RETURN
	END
	
