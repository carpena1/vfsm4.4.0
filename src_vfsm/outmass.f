        SUBROUTINE OUTMASS(VL,FWIDTH,SWIDTH,SLENGTH,TRAI,LISFIL,ISCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      This subroutine processes the output hydrograph and find components    C
C       of the water balance and hydrograph. The results go in "filename.osm" C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,OS,RAINL,ZW,TW,RVH
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/GRASSD3/SUSMASS,WEDGEMASS,NFUP
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD
      COMMON/WQ1/VKD,CCP
      COMMON/IWQ2/NDGDAY,IDG,IWQ
      COMMON/WQ3/DGKREF,FC,DGPIN,DGML,DGT(366),DGTHETA(366)
      CHARACTER*75 DUMMY
      CHARACTER*75 LISFIL(13)

C-------------- Summarize results from filename.ohy file -------------

      CLOSE(11)
      CLOSE(13)
      OPEN(11,FILE=LISFIL(6),STATUS='OLD')
      OPEN(13,FILE=LISFIL(8),STATUS='OLD')

      WRITE(10,*)'INPUTS'
      WRITE(10,*)'------'      
      SUM0=0.D0
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      TIME0=0.D0
5      READ(11,'(A)')DUMMY
      IDX=INDEX(DUMMY,'OUTFLOW')
      IF(IDX.NE.0) GOTO 6
      WRITE(10,*)DUMMY
      GOTO 5
6     READ(11,'(A)')DUMMY
      READ(11,'(A)')DUMMY
      BIG1=0.D0
      NZERO=1
      INI=0
      RAIN_E0=0.D0
      OUTF0=0.D0
      UPIN0=0.D0
      VF0=0.D0
      TINI=0.D0
      TEND=0.D0
      DO 10 I=1,101
            READ(11,*,END=16)TIME1,OUTF1,CUMFLOW,RAIN_E1,UPIN1,CUMIF,
     &             VF1
            TIMEINCR=TIME1-TIME0
            AREA0=TIMEINCR*(UPIN1+UPIN0)/2.D0
            SUM0=SUM0 + AREA0
            AREA1=TIMEINCR*(OUTF1+OUTF0)/2.D0
            SUM1=SUM1 + AREA1
            AREA2=TIMEINCR*(RAIN_E1+RAIN_E0)/2.D0
            SUM2=SUM2 + AREA2
            AREA3=TIMEINCR*(VF1+VF0)/2.D0
            SUM3=SUM3 + AREA3
            UPIN0=UPIN1
            OUTF0=OUTF1
            RAIN_E0=RAIN_E1
            VF0=VF1
            TIME0=TIME1
            C1= OUTF1
            BIG1=DMAX1(BIG1,C1)
            IF(C1.EQ.BIG1)THEN
                  TBIG=TIME1
                  QBIG=C1
            ENDIF
            IF(NZERO.EQ.1.AND.C1.GT.0.D0.AND.INI.EQ.0) THEN
                  TINI=TIME1
                  NZERO=0
                  INI=1
               ELSEIF(NZERO.EQ.0.AND.C1.EQ.0.D0) THEN
                  NZERO=1
                  TEND=TIME1
            ENDIF
            IF(C1.GT.0.D0)TEND=TIME1
10      CONTINUE
12      READ(13,'(A)')DUMMY
      IDX=INDEX(DUMMY,'gsI')
      IF(IDX.NE.0) GOTO 16
      WRITE(10,*)DUMMY
      GOTO 12
16    CONTINUE
      WRITE(10,*)'OUTPUTS'
      WRITE(10,*)'-------'      
      WRITE(10,*)'Water balance'
      WRITE(10,*)'-------------'

C---------------Calculate total rainfall for event in m3 ----------

      TOTRAIN=TRAI*VL*FWIDTH
      WRITE(10,600)TOTRAIN      
      Vout=SUM1
      Vie=SUM2*VL*FWIDTH
      Vin=SUM0
      VFE= SUM3*VL*FWIDTH
      
C---------------Water and sediment balance for the event----------

      WAT_IN=Vin+TOTRAIN
      VF=WAT_IN - Vout
      if(AGA.LE.0.D0)THEN
          VF=0.d0
        ELSEIF(VF.LT.0.d0) then
          VF=VFE
      ENDIF           
      WAT_OUT= Vout+VF
      WAT_BAL=WAT_IN-WAT_OUT
      IF (WAT_IN.GT.0.d0) THEN
		  WAT_ERR=WAT_BAL/WAT_IN*100
	    ELSE
		  WAT_ERR=0.0d0
	  ENDIF		 

      WRITE(10,300)Vin
      WRITE(10,150)Vout
      WRITE(10,400)VF
      WRITE(10,*)' '
      WRITE(10,*)'Hydrology'
      WRITE(10,*)'-----------'      
      WRITE(10,700)TINI
      WRITE(10,800)TBIG,QBIG
      WRITE(10,900)TEND
      WRITE(10,*)' '
      WRITE(10,*)'Sediment'
      WRITE(10,*)'--------'      
      TTE=1.D0
      IF(GSIMASS.GT.0.D0) TTE=(GSIMASS-GSOMASS)/GSIMASS      
      FWID=FWIDTH*100.D0
      WRITE(10,1100)GSIMASS,GSIMASS*FWID
      WRITE(10,1200)GSOMASS,GSOMASS*FWID
      WRITE(10,1010)TTE*100

c-------Sediment wedge final shape ------------------------------

      WRITE(10,*)'                    ----------------'
      WRITE(10,1250)
      VLT=VLCM-XTOLD
      X1=YTOLD/SC
      WRITE(10,1300)YTOLD
      WRITE(10,1400)XTOLD
      WRITE(10,1500)VLT
      WRITE(10,1600)X1
      WRITE(10,1700)DEP

c--------Mass balance -------------------------------------------

      GAMMASB=(1.D0-POR)*PART(3)
      IF(NFUP.EQ.0.AND.YTOLD.LT.H) THEN
c------------Triangular wedge
            BMASS=(DEP*VLT+(XTOLD+X1)*YTOLD*0.5D0)*GAMMASB
         ELSEIF(NFUP.EQ.0.AND.YTOLD.EQ.H) THEN
c------------Trapezoidal wedge
            X3=YTOLD/Se
            X2=XTOLD-X3
            TRBOT=X1+X2+X3
            TRTOP=X2
            BMASS=.5D0*(TRTOP+TRBOT)*YTOLD*gammasb
c------------if strip is filled calculate as rectangle (H*VLCM)
         ELSE
            BMASS=(H*VLCM+(H/SC)*H*0.5D0)*GAMMASB
      ENDIF
      IF (GSIMASS.GT.0.d0) THEN
		  SED_ERR=(BMASS-(GSIMASS-GSOMASS))/(GSIMASS-GSOMASS)*100.D0
	    ELSE
			SED_ERR=0.d0
	  ENDIF		 
	              
      WRITE(10,1800) SED_ERR  
      WRITE(10,1825)WAT_ERR
      IF(NFUP.EQ.1) THEN
            WRITE(10,1850)
        ELSEIF(YTOLD.EQ.H) THEN
            WRITE(10,1875)
      ENDIF            

c-(02/1999)--Filter performance summary output file  -------------

      WRITE(15,*)' '
      WRITE(15,*)'      Summary of Buffer Performance Indicators:'
      WRITE(15,*)' '
      WRITE(15,1900)SWIDTH*SLENGTH
      WRITE(15,2000)SLENGTH
      WRITE(15,2100)SWIDTH
      WRITE(15,2200)VL
      WRITE(15,2225)FWIDTH
      WRITE(15,2250)VN1
      WRITE(15,*)' '
      WRITE(15,2300)VL/SLENGTH*100.d0
      WRITE(15,2400)TOTRAIN/(FWIDTH*VL)*1000.d0
      WRITE(15,2450)TOTRAIN
      WRITE(15,2500)Vin/(SWIDTH*SLENGTH)*1000.d0
      WRITE(15,2550)Vin
      WRITE(15,2600)Vout/(SWIDTH*SLENGTH+FWIDTH*VL)*1000.d0
      WRITE(15,2650)Vout
      WRITE(15,2675)VF
c--- rmc 04/20/03 --Fix for Vin(Q)=0 or Vout-0
      SMIN=GSIMASS*FWID/1000.d0
      SMOUT=GSOMASS*FWID/1000.d0
      if(Vin.le.0) then
          WRITE(15,2685)0.d0
          WRITE(15,*)' '
          WRITE(15,2700)SMIN
          WRITE(15,2800)0.d0          
       else
          WRITE(15,2685)Vout/Vin
          WRITE(15,*)' '
          WRITE(15,2700)SMIN
          WRITE(15,2800)GSIMASS*FWID/(Vin*1000.d0)        
      endif
      WRITE(15,2900)SMOUT
      if(Vout.le.0) then
            WRITE(15,3000)0.d0
            WRITE(15,3050)(GSIMASS-GSOMASS)*FWID/1000.d0
            WRITE(15,3075)0.d0            
       else
            WRITE(15,3000)GSOMASS*FWID/(Vout*1000.d0)
            WRITE(15,3050)(GSIMASS-GSOMASS)*FWID/1000.d0
            WRITE(15,3075)SMOUT/SMIN 
      endif
c--- rmc 04/20/03 --end of fix
      
      WRITE(15,*)' '
      WRITE(15,3100)VLT/100.d0
      WRITE(15,3200)XTOLD/100.d0
      IF(NFUP.EQ.1) THEN
            WRITE(15,1850)
        ELSEIF(YTOLD.EQ.H) THEN
            WRITE(15,1875)
      ENDIF            

c-(08/2008)-Water quality summary output file  -------------
      
      IF(IWQ.GT.0) THEN
            IF (SMIN.EQ.0.D0) THEN
                  FPH=Vin*1000.D0/(VKD*0.001D0)
              ELSE
                  FPH=Vin*1000.D0/(VKD*SMIN)
            ENDIF
            IF(WAT_IN.gt.0.d0)THEN
				 PDQ=VF/WAT_IN*100.D0
			 ELSE
			     PDQ=100.D0
			ENDIF
            PDSED=TTE*100.D0
            DELTAP=24.79D0+.54D0*PDQ+.52D0*PDSED-2.42D0*
     &       DLOG(FPH+1.D0)-.89D0*CCP
            IF (DELTAP.GT.100.D0.OR.PDQ.GE.100.d0) DELTAP=100.D0
            IF (PDQ.EQ.0.d0.AND.PDSED.EQ.0.d0) DELTAP=0.D0            
            IF (DELTAP.LT.0.d0) DELTAP=0.D0            
            WRITE(18,*)
            WRITE(18,*)'Outputs for Water Quality'
            WRITE(18,3500)      
            WRITE(18,4000)VIN
			IF(SMIN.gt.999999.d0) THEN
                WRITE(18,4105)SMIN
			  ELSE
                WRITE(18,4100)SMIN
			ENDIF
			IF(dabs(FPH).le.9999.d0) THEN
                WRITE(18,4200)FPH
			  ELSE
                WRITE(18,4205)FPH
			ENDIF
            WRITE(18,4300)PDQ
            WRITE(18,4400)PDSED
			IF(Vin.gt.0.d0) THEN
				RRED=100.D0*(1-Vout/Vin)
			  ELSE
				RRED=100.d0
			ENDIF
			IF(dabs(RRED).gt.200.d0) THEN
                WRITE(18,4455)RRED
			  ELSE
                WRITE(18,4450)RRED
			ENDIF
            WRITE(18,*)
            WRITE(18,4500)DELTAP
            IF(IDG.GT.0.AND.IDG.LE.4) THEN
                 SAREA=SLENGTH*SWIDTH
                 DGMI=DGPIN*SAREA
                 DGMO=DGMI*(1.d0-DELTAP/100.d0)
                 IF(VIN.le.0.d0.and.SMIN.le.0.d0) THEN
                     DGSI=0.d0
                   ELSE
                     DGSI=DGMI*VKD/(VIN*1000.d0+SMIN*VKD)
                 ENDIF
                 DGMF=DGMI*DELTAP/100.d0
                 DGMFSED=DGSI*SMIN*PDSED/100.d0
c-(02.18.14)---rmc----- for high sorption pesticides check if all pesticide is sediment-bonded
				 IF(DGMFSED.GT.DGMF) DGMFSED=DGMF
                 DGMFF=DGMF-DGMFSED
                 DGCF=DGMFF/(VF*1000.d0)
                 DGROB=(1.d0-OS)*2.65d0
                 DGMML=(OS+VKD*DGROB)*DGCF*DGML/100.d0*VL*FWIDTH
                 IF(ISNAN(DGMML)) DGMML=0.d0                 
				 IF((DGMFSED+DGMML).GT.DGMF) DGMML= DGMF-DGMFSED
                 DGMRES=DGMFSED+DGMML
                 DGBETA=-0.7d0
                 DO 40 I=1,NDGDAY
c------------------ Scheme for 4 degradation rates (IDG=1: EU-FOCUS; 2: US-EPA K=Kref; 
c-------------------3: k=k(T); 4:k=k(theta))
					 IF(IDG.EQ.1) THEN
		                 DGC1=65.4d0/(8.314d0/1000.d0)
		                 DGC2=1.d0/293.15d0
	                     DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(I)+273.15d0))))
	                     DGKTHETA=(DGTHETA(I)/FC)**DGBETA						 
  					   ELSEIF(IDG.EQ.2) THEN
						 DGKTEMP=1.D0
						 DGKTHETA=1.D0
					   ELSEIF(IDG.EQ.3) THEN
		                 DGC1=49.5d0/(8.314d0/1000.d0)
		                 DGC2=1.d0/298.15d0
	                     DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(I)+273.15d0))))
  						 DGKTHETA=1.D0
  					   ELSEIF(IDG.EQ.4) THEN
  						 DGKTEMP=1.D0
	                     DGKTHETA=(DGTHETA(I)/FC)**DGBETA
					 ENDIF			 	 					 
                     DGKI=DGKREF*DGKTEMP*DGKTHETA
                     DGMRES=DGMRES*DEXP(-DGKI)
                     IF(ISNAN(DGMRES)) DGMRES=0.d0
40               CONTINUE
                 WRITE(18,*)
                 WRITE(18,3525)'Pesticide mass balance/degradation(IDG='
     &           ,IDG
                 WRITE(18,3550)      
                 WRITE(18,*)'Total mass in filter:'
                 WRITE(18,*)
                 write(18,4600)DGMI
                 write(18,4700)DGMO
                 write(18,4800)DGMFSED
                 write(18,4900)DGMML
                 write(18,5000)DGMML+DGMFSED
                 write(18,5100)DGMRES,NDGDAY
                 WRITE(18,*)
                 WRITE(18,*)'Normalized values by source area:'
                 WRITE(18,*)
                 WRITE(18,5150)SAREA
                 write(18,5200)DGPIN
                 write(18,5300)DGMO/SAREA
                 write(18,5400)DGMFSED/SAREA
                 write(18,5500)DGMML/SAREA
                 write(18,5600)(DGMML+DGMFSED)/SAREA
                 write(18,5700)DGMRES/SAREA,NDGDAY
                 WRITE(18,*)
                 WRITE(18,*)'Pesticide mass partition for outflow'
                 WRITE(18,*)
c----------------Solid/liquid phase partitioning assuming linear equilibrium based on Kd (VKD)
c----------------between outgoing sediment (SMOUT) and water (VOUT).
				 IF(VOUT.GT.0.D0) THEN
					 DGMOP=(DGMO*SMOUT*VKD)/(VOUT*1000.D0+SMOUT*VKD)
				     DGMOD=DGMO-DGMOP
				   ELSE
                      DGMOP=0.D0
					  DGMOD=0.D0
				 ENDIF
                 write(18,5800) DGMOP
                 write(18,5900)DGMOD
           ENDIF
      ENDIF 

c-------Output message at end of program -----------------

      IF(ISCR.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*)'...FINISHED...','VFSMOD v4.3.2 01/2016'
        WRITE(*,*)
      ENDIF            
      
150   FORMAT('Volume from outflow = ', E14.4,' m3')
200   FORMAT('Volume from i_e     = ', E14.4,' m3',F14.2,'%')
300   FORMAT('Volume from up-field= ', E14.4,' m3')
400   FORMAT('Volume infiltrated  = ', E14.4,' m3')
500   FORMAT(F8.2,3E12.4,I6)
600   FORMAT('Volume from rainfall= ', E14.4,' m3')
700   FORMAT('Time to beginning   = ', E14.4,' s')
800   FORMAT('Time and q at peak  = ', E14.4,' s ',E14.4,' m3/s')
900   FORMAT('Time to end runoff  = ', E14.4,' s')
1010  FORMAT('Trapping efficiency =      ', F5.1,' %')
1100  FORMAT('Sediment inflow     = ', E14.4,' g/cm',E14.4,' g')
1200  FORMAT('Sediment outflow    = ', E14.4,' g/cm',E14.4,' g')
1250  FORMAT('Sediment deposition :')
1300  FORMAT(8x,'- Sediment wedge depth,   Y(t)   =',F7.2,'cm')  
1400  FORMAT(8x,'- Sediment wedge length,  X2(t)  =',F7.2,'cm')  
1500  FORMAT(8x,'- Effective filter length,L(t)   =',F7.2,'cm')  
1600  FORMAT(8x,'- Sediment tail at field, X1(t)  =',F7.2,'cm')  
1700  FORMAT(8x,'- Sediment depth in low section  =',F7.2,'cm')
1800  FORMAT(8x,'- Rough sediment balance error   =',F7.2,'%')
1825  FORMAT(8x,'- Rough water balance error      =',F7.2,'%')
1850  format(/,66('*'),/,'*  WARNING: Strip filled up!',37x,'*',
     &     /,66('*'))
1875  format(/,67('*'),/,'*  WARNING: Top of vegetation reached',
     &     ' - trapezoidal wedge started *',/,67('*'))
1900  FORMAT(1x,F10.1,' m^2 = Source Area (input)')
2000  FORMAT(1x,F10.2,' m   = Source Flow Length (input)')
2100  FORMAT(1x,F10.2,' m   = Source Area Width (input)')
2200  FORMAT(1x,F10.2,' m   = Filter Strip Length (input)')
2225  FORMAT(1x,F10.2,' m   = Filter Strip Width (input)')
2250  FORMAT(1x,F10.3,'     = Mean Filter Mannings Roughness (input)')
2300  FORMAT(1x,F10.2,' %   = ',
     1        'Ratio of Filter Length to Source Flow Length')
2400  FORMAT(1x,F10.3,' mm  = Total Rainfall')
2450  FORMAT(1x,F10.3,' m3  = Total Rainfall on Filter')
2500  FORMAT(1x,F10.3,' mm  = Total Runoff from Source (mm depth over',
     &    ' Source Area)')
2550  FORMAT(1x,F10.3,' m3  = Total Runoff from Source')
2600  FORMAT(1x,F10.3,' mm  = Total Runoff out from Filter (mm depth',
     &    ' over Source+Filter)')
2650  FORMAT(1x,F10.3,' m3  = Total Runoff out from Filter')
2675  FORMAT(1x,F10.3,' m3  = Total Infiltration in Filter')
2685  FORMAT(1x,F10.3,'     = Runoff Delivery Ratio')
2700  FORMAT(F12.3,' kg = Mass Sediment Input to Filter')
2800  FORMAT(F12.3,' g/L= Concentration Sediment in Runoff from',
     &    ' source Area')
2900  FORMAT(F12.3,' kg = Mass Sediment Output from Filter')
3000  FORMAT(F12.3,' g/L= Concentration Sediment in Runoff exiting',
     &   ' the Filter')
3050  FORMAT(F12.3,' kg = Mass Sediment retained in Filter')
3075  FORMAT(2x,F10.3,'    = Sediment Delivery Ratio')
3100  FORMAT(1x,F10.2,' m   = Effective Filter Length')
3200  FORMAT(1x,F10.2,' m   = Wedge Distance')
3500  FORMAT(26('-'))
3525  FORMAT(A40,I1,')') 
3550  FORMAT(50('-'))
4000  FORMAT(F10.3,' m3 = Runoff inflow')
4100  FORMAT(F10.3,' Kg = Sediment inflow')
4105  FORMAT(E10.3,' Kg = Sediment inflow')
4200  FORMAT(F10.3,'    = Phase distribution, Fph')
4205  FORMAT(E10.3,'    = Phase distribution, Fph')
4300  FORMAT(F10.3,' %  = Infiltration (dQ)')
4400  FORMAT(F10.3,' %  = Sediment reduction (dE)')
4450  FORMAT(F10.3,' %  = Runoff inflow reduction')
4455  FORMAT(E10.3,' %  = Runoff inflow reduction')
4500  FORMAT(F10.3,' %  = Pesticide reduction (dP)')
4600  FORMAT(E15.6,' mg = Pesticide input (mi)')
4700  FORMAT(E15.6,' mg = Pesticide output (mo)')
4800  FORMAT(E15.6,' mg = Pesticide retained in filter, sediment-',
     &       'bonded (mf,sed)')
4900  FORMAT(E15.6,' mg = Pesticide retained in filter, mixing ',
     &       'layer (mml)')
5000  FORMAT(E15.6,' mg = Pesticide surface residue at the end of',
     &       ' this event (mres)')
5100  FORMAT(E15.6,' mg = Pesticide surface residue at next event',
     &       ' (after degradation,',I3,' days)') 
5150  FORMAT(F15.2,' m^2  = Source Area (input)')
5200  FORMAT(E15.6,' mg/m2= Pesticide input (mi)')
5300  FORMAT(E15.6,' mg/m2= Pesticide output (mo)')
5400  FORMAT(E15.6,' mg/m2= Pesticide retained in filter, sediment-',
     &       'bonded (mf,sed)')
5500  FORMAT(E15.6,' mg/m2= Pesticide retained in filter, mixing ',
     &       'layer (mml)')
5600  FORMAT(E15.6,' mg/m2= Pesticide surface residue at the end ',
     &       'of the event (mres)')
5700  FORMAT(E15.6,' mg/m2= Pesticide surface residue at next event',
     &       ' (after degradation,',I3,' days)') 
5800  FORMAT(E15.6,' mg = Pesticide output in solid phase (mop)')
5900  FORMAT(E15.6,' mg = Pesticide output in liquid phase (mod)')

      RETURN
      END


