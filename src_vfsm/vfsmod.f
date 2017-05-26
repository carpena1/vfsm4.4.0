      PROGRAM VFSMOD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     WRITTEN FOR: Ph.D.dissertation and later modified for distribution     C
C     RE-WRITTEN : combined version, July 1993                               C
C     Last Updated: 01/2016 v4.3.2                                           C
C     Written by: Rafael Munoz-Carpena         John E. Parsons               C
C                 ABE-University of Florida    BAE, NC State University      C
C                 Gainesville, FL 32611        Raleigh, NC 27695-7625 (USA)  C
C                 e-mail: carpena@ufl.edu      john_parsons@ncsu.edu         C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     PROGRAM CALCULATE OVERLAND FLOW, SEDIMENT AND CONTAMINANT FILTRATION   C
C     THROUGH A VEGETATIVE FILTER STRIP OF AN INFLOW HYDROGRAPH FROM AN      C
C     ADJACENT FIELD DURING A STORM EVENT. THE PROGRAM HANDLES THE CASE OF   C
C     VARYING ROUGHNESS COEFFICIENT (Manning's n) AND SLOPE AT THE NODES,    C
C     AND TIME DEPENDENT INFILTRATION FOR THE DOMAIN.                        C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     This program solves the kinetic wave aproximation of the Saint-        C
C     Vennant's (1881) equations for overland flow (KW) for the 1-D case     C
C     as  presented by Lighthill and Whitham (1955) such as:                 C
C                                                                            C
C            dh     dq                                                       C
C           ---- + ---- = ie(t)            (Continuity equation)             C
C            dt     dx                                                       C
C                                                                            C
C           So = Sf                        (Momentum equation)               C
C                                                                            C
C     Then the momentum equation is used as a link between the two           C
C     variables since then we have (Manning's):                              C
C                         1/2    5/3                                         C
C           q = q(h) = (So  /n) h                                            C
C                                                                            C
C     Where h is depth of overland flow [L], q is the flow per unit width    C
C     of the plane [L^2/T], So is the slope of the plane, Sf is the          C
C     hydraulic or friccion slope, and n is Manning's roughness cofficient   C
C     [L^1/6]. The initial and boundaty conditions can be summarized as:     C
C                                                                            C
C              h=0 ;  0 <=x<= L ; t <=0                                      C
C              h=ho;   x = 0    ; t > 0                                      C
C                                                                            C
C     where ho can be 0, a constant or a time dependent funtion.             C
C                                                                            C
C     The numerical method is based on a N+2 upwinding Petrov-Galerkin       C
C     finite element method approximation  for the spacial derivatives       C
C     and a time weighting finite difference approximation for the time      C
C     derivatives.                                                           C
C                                                                            C
C     The non-linearity of the equation {q=q(h)} is taken care of using      C
C     the Picard iterative scheme inside every time step lagging 2/3 of      C
C     the power of h in q ,[5/3 = 2/3(m)+1 (m+1)] for the iteration level    C
C      m, such as:                                                           C
C                                  m+1     m                                 C
C                            [A] {h} = {b(h)}                                C
C                                                                            C
C     In this program the core of the time step solution is taken care       C
C     of following this steps:                                               C
C                                                                            C
C     1- Form the system matrix [A] of constant coefficients                 C
C     2- Perfom LUD decomposition over this matrix [A]                       C
C     3- Form the system matrix [BM] of constant coefficients                C
C     4- Form r.h.s of equation (vector {b}=[BM]{xo} for each time step)     C
C     5- solve for [A],{b} to get a {x} for that time step                   C
C     6- Repeat 4 & 5 until convergence of that time step                    C
C     7- Repeat 3 & 6 until completion of desired number of time steps       C
C                                                                            C
C     The ie(t) term of the continuity equation is calculated for each       C
C     time step by the extention of the Green-Ampt method as proposed by     C
C     Mein& Larson, (1967) and Chu (1975),                                   C
C                                                                            C
C     The overland flow solution is linked to a submodel to calculate        C
C     sediment transport  on a grass filter. The information from the        C
C     submodel is used to assemble the vector {b} during the procedure       C
C     described above.                                                       C
C     The sediment filtration is based on the method proposed by:            C
C                                                                            C
C  1. Tollner et al. (1976). "Suspended sediment filtration capacity of      C
C     simulated vegetation". Trans. ASAE. 19(4):698-682.                     C
C  2. Tollner et al. (1977). "Sediment deposition patterns in simulated      C
C     grass filters". Trans. ASAE. 20(5):940-944.                            C
C  3. Barfield et. al (1979)"Filtration of sediment by simulated vegetation  C
C     I. Trans. ASAE, 22(3):540-548.                                         C
C  4. Hayes et. al (1979)"Filtration of sediment by simulated vegetationII"  C
C     Trans. ASAE, 22(5):1063-1067                                           C
C  5. Hayes et. al (1984)"Performance of grass filters under laboratory and  C
C     field Conditions".Trans. ASAE, 27(5):1321-1331, to account for         C
C     triangular upslope deposition and particle and size distribution       C
C  6. Wilson et al (1981)"A Hydrology and sedimentology model: Part I.       C
C     Modeling techniques. U.of Kentucky. Lexington. This is a major         C
C     rewrite of the prodedures involved.                                    C
C  7. Haan et al (1994)"Design Hydrology and Sedimentology for Small         C
C     Catchments". Prentice-Hall. Chapter 9C contains updated and clearly    C
C     presented procedures for sediment trapping and wedge formation         C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C                 SUBROUTINES BY ORDER OF APPAREANCE                         C
C                                                                            C
C     INI, INPUTS,QUAD, FORMA, SHAPEF, ASSM, BCA, FORMB, MODIFY, FACTOR,     C
C     SOLVE, FLOW, UPDATE, CONVER, KWWRITE, GASUB, OUTMASS, GRASSED,         C
C     GRASSIN, OCF, EINSTEIN, STEP3, POINTS                                  C
C                                                                            C
C          DEFINITION OF GLOBAL VARIABLES FOR OVERLAND FLOW SOLUTION         C
C                                                                            C
C      A(I,J)= SYSTEM MATRIX, SQUARE OF DIMENSIONS NxN, ie. [A]              C
C      B(I)= RIGHT HAND SIDE VECTOR OF DIMENSIONS 1xN , ie  {b}              C
C      DPSI(L) = DERIVATIVE OF BASIS FUNCTIONS                               C
C      DR= DURATION OF THE RAINFALL (s)                                      C
C      DT = INCREMENT OF TIME (s)                                            C
C      DX= SPACE STEP (m)                                                    C
C      MAXITER= MAXIMUM NUMBER OF ITERATIONS ALOWED                          C
C      MFLAG= CONVERGENCE FLAG (0, NO CONVERGENCE; 1, CONVERGENCE)           C
C      N=ACTUAL NUMBER OF NODES IN THE DOMAIN                                C
C      NDT = NUMBER OF TIME STEPS                                            C
C      NELEM= ACTUAL NUMBER OF ELEMENTS IN THE DOMAIN                        C
C      NL= ORDER OF THE INTEGRATION RULE OVER EACH ELEMENT                   C
C      NMAX= MAXIMUM NUMBER OF EQUATIONS AND VARIABLES THAT CAN BE SOLVED    C
C      NPOL=NUMBER OF NODAL POINTS OVER EACH ELEMENT (POLYNOMIAL DEG +1)     C
C      OUT= 0, print values at the downstream end of the plane (hydrograph)  C
C           1, Print values for all nodes at each time step                  C
C      QK(MAXEQN)= NODAL ALPHA IN MANNING'S UNIFORM FLOW EQUATION            C
C      R= Lateral inflow (m/s)                                               C
C      RN = MANNING'S ROUGHNESS COEFFICIENT                                  C
C      SO = SLOPE OF THE ELEMENT                                             C
C      SR= DURATION OF THE SIMULATION (s)                                    C
C      THETAW= TIME-WEIGHT FACTOR                                            C
C      VL= LENGTH OF THE VFS, DIMENSION ALONG FLOW PATH (m)                  C
C      W(L) = GAUSS QUADRATURE WEIGHTS                                       C
C      X(I)= SOLUTION VECTOR (H[m]), DIMENSION 1xN, AT TIME STEP L+1         C
C      XI(L)= GAUSS QUADRATURE POINT                                         C
C      XM(I)= SOLUTION VECTOR, DIMENSION 1xN, AT ITERATION M, t STEP L+1     C
C      X0(I)= SOLUTION VECTOR, DIMENSION 1xN, AT  TIME STEP L                C
C                                                                            C
C                                                                            C
C         DEFINITION OF GLOBAL VARIABLES FOR INFILTRATION SOLUTION           C
C                                                                            C
C     AGA = Green-Ampt's "A", saturated hydraulic conductivity, Ks (m/s)     C
C     BCROFF(200,2)= Boundary condition at the upstream node (inflow from    C
C          adjacent field.                                                   C
C     BGA = Green-Ampt's "B" = Ks*Sav*M   (m2/s)                             C
C     CU= Chu's surface ponding indicator at end rain period (<0, ponded)    C
C     F= Cumulative infiltration (m)                                         C
C     FPI= Instantaneous infiltration rate (m/s)                             C
C     L= rainfall period                                                     C
C     LO = index to show if time step is in the same rainfall period (LO=L)  C
C     NPOND= surface ponding indicator at beginning rain period (1, ponded)  C
C     NSTART= Indicates that the start of water on filter surface            C
C     NEND= Indicates that the end of runoff is reached                      C
C     PS= Cumulative recipitation in m.                                      C
C     PSOLD= Cumulative recipitation in m for last rainfall period.          C
C     PST= Total cumulative recipitation in m.                               C
C     PSI(L) = BASIS FUNCTIONS (QUADRATIC LAGRANGIAN POLYNOMIALS)            C
C     RAIN(200,2)= Times (s) and rainfall rates (m/s) over the VFS.          C
C     RO= Cumulative runoff rate at the node (without considering BCRO)      C
C     SM= Maximum surface storage (m)                                        C
C     STO= Cumulative surface storage (m)                                    C
C     TP, TPP= Chu's (1978) tp and tp' coefficients                          C
C     TRAI= Total cumulative rainfal (m)                                     C
C                                                                            C
C           DEFINITION OF VARIABLES FOR SEDIMENT SOLUTION                    C
C                                                                            C
C     Ss= spacing of the filter media elments (cm)                           C
C     Sc= filter main slope                                                  C
C     n= Manning's n= 0.0072 for cilindrical media (s/cm^1/3)                C
C     q= overland flow (cm2/s)                                               C
C     df= depth of flow at D(t) (cm)                                         C
C     Vm= depth averaged velocity at D(t)(cm/s)                              C
C     Rs= hydraulic radius of the filter (cm)                                C
C     dp= particle size, diameter (cm)                                       C
C     gamma, gammas=  water and sediment weight density (g/cm3)              C
C     gs2=gsd: sediment load entering downstream section (g/s/cm)            C
C     Rss= hydraulic radius of the filter at B(t) (cm)                       C
C     dfs= depth of flow at B(t) (cm)                                        C
C     Vms= depth averaged velocity at B(t) (cm/s)                            C
C     Se= equilibrium slope at B(t)                                          C
C     f= fraccion trapped in the depodition wedge                            C
C     ico= flag to select feedback to overland flow solution of new slopes   C
C           and roughness (0=no, 1=yes)                                      C
C     coarse= % of particles from incoming sediment with diameter > 0.0037   C
C            cm (coarse fraction that will be routed through wedge).         C
C                                                                            C
C     NOTE: units in sediment transport calculations are in CGS system       C
C           (cm,g,s), including Manning's n                                  C
C                                                                            c
c     Change Log: See CHANGES file in source code directory                  c
C                                                                            c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/PAR/QK(200),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/CINT/XI(20,20),W(20,20)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z
      COMMON/GA2/LO,NPOND,NPFORCE,NPFORCE0
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,OS,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC   
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/GRASSD3/SUSMASS,WEDGEMASS,NFUP
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD
      COMMON/KWW/TIMELAST,OFLOWLAST,OSUMFLOW,QFIELDLAST,FSUMFLOW
      COMMON/WQ1/VKD,CCP
      COMMON/IWQ2/NDGDAY,IDG,IWQ
      COMMON/WQ3/DGKREF,FC,DGPIN,DGML,DGT(366),DGTHETA(366)
      DIMENSION A(MAXEQN,MAXBND),B(MAXEQN),B0(MAXEQN)
      DIMENSION X(MAXEQN),X0(MAXEQN),XM(MAXEQN),Q0(MAXEQN),QM(MAXEQN)
      DIMENSION PGPAR(4)
      DIMENSION BCROFF(200,2),RAIN(200,2),NODEX(4)
      CHARACTER*75 LISFIL(13)

C------Print banner, get I/O filenames and open them -----------

      CALL FINPUT(LISFIL,INARGS,ISCR)

c--for debug--gasubwt-----
c       open(19,file='wt_debug.txt', status='unknown')
c--for debug--gasubwt-----

C-------Initialize matrices--------------------------

      CALL INI(A,B,X,XM,X0,Q0,QM,SSE,NODEX)

C------------Read inputs for sediment problem---------------------------

      CALL GRASSIN(ICOARSE,COARSE,LISFIL,INARGS)

C-----Get inputs and parameters for hydrology problem------

      CALL INPUTS(N,NBAND,NRAIN,RAIN,NBCROFF,BCROFF,TE,QMAX,VL,FWIDTH,
     &  SWIDTH,SLENGTH,PGPAR,NCHK,LISFIL,INARGS,ISCR)


C------Get the Gauss quadrature parameters-------

      CALL QUAD

C-------Assemble the system matrix A -------------

      CALL FORMA(A,NBAND,PGPAR)

C-------Perform LU decomposition over A -----------

      CALL FACTOR(A,N,NBAND)

C-------Numerical time dependent solution-----------

      MAXIT=100
      TIME=0.D0
      RO=0.D0
      ROLD=0.D0
      PS=0.D0
      STO= 0.D0
      F=0.D0
      XTOLD=0.D0
      YTOLD=0.D0
      TOTF=0.D0
      TRAI=0.D0
      NPOND=0
      NPFORCE=0
      NPFORCE0=0
      z=0.d0
      LO=0
      PZERO=1.0D-8
      NSTART=0
      NEND=0
      MFLAG=0
      NWRITE=NDT/100
      IF(NDT.LT.100) NWRITE=NDT
      VLCM=VL*100.D0
      QTEMP=0.d0
      NFUP=0
      TW=0.d0
      DO 5 I=1,4
            QSED(I)=0.D0
5     CONTINUE                
      DO 40 LCOUNT=1,NDT
            TIME=DT*LCOUNT
C------- For each time step select the rainfall intensity period (L)
            R=0.D0
            DO 10 I=1,NRAIN-1
                  IF(TIME.GT.RAIN(I,1).AND.TIME.LE.RAIN(I+1,1)) L=I
10          CONTINUE
C------- Interpolate BC (in depth,m) at the first node of system
C------- (from incoming hydrograph)  
            BCRO=0.D0
            BCROQ=0.D0			
            DO 15 I=1,NBCROFF-1
            IF(TIME.GT.BCROFF(I,1).AND.TIME.LE.BCROFF(I+1,1)) THEN
                  BCROQ=(TIME-BCROFF(I,1))/(BCROFF(I+1,1)-BCROFF(I,1))*
     &                 (BCROFF(I+1,2)-BCROFF(I,2))+BCROFF(I,2)
                  BCRO=((BCROQ/FWIDTH)/QK(1))**(3.D0/5.D0)
            ENDIF
15          CONTINUE

C------- Get effective rainfall and control execution of overland flow for ---
C------- an infiltrating surface by calling Green-Ampt model. If a certain ---
C------- node (nchk, selected by the user) is flooded, i.e. x(nchk)>0 the ----
C------- surface is assumed flooded (NPFORCE=1,rmc03/2011), and the max. GA --
C--------infiltration capacity is selected overriding other rainfall excess --
C--------calculations 
            NEND=0
            IF(X(NCHK).GT.PZERO) NPFORCE=1
            IF(BCRO.EQ.0.D0.AND.X(NCHK).EQ.0.D0.AND.NSTART.EQ.1)NEND=1
C-----------rmc-05/2003 consider infiltration for infiltrating plane 
c-----------(mod for thetai>porosity) ---
            IF(AGA.GT.0.D0.AND.WTD.eq.999.d0)THEN
c---------------No shallow water table present
                CALL GASUB(TIME,DT,L,R,RAIN,NEND,TRAI)
             ELSEIF(AGA.GT.0.D0.AND.WTD.ne.999.d0)THEN
c---------------Shallow water table depth present at WTD
                CALL GASUBWT(TIME,DT,L,R,RAIN,NEND,TRAI)
             ELSE
                R=RAIN(L,2)
                TRAI=TRAI+DT*(R+ROLD)*0.5D0
                ROLD=R
            ENDIF
            NSTART=1
            IF(R.LE.0.D0.AND.BCRO.EQ.0.D0.AND.X(NCHK).EQ.0.D0)NSTART=0
c            IF(R.GT.0.D0.OR.BCRO.GT.0.D0)NSTART=1
c            print*,NPFORCE,NEND,NSTART,MFLAG

C-----------Form of r.h.s vector for that time step -------

            CALL FORMB(B0,X0,Q0,N,BCRO,PGPAR)

C---------------Start Picard iteration-----------------

            M=0
            MFLAG=0
C-----------if no excess rain or field inflow skip overland flow (NSTART=0)- 
            IF(NSTART.EQ.0)MFLAG=1
            DO 20 WHILE (M.LT.MAXITER.AND.MFLAG.EQ.0)
                  M= M+1

C----------------------Update {b} = {bm} -------------------

                  CALL UPDATE(N,B0,B)
                  CALL MODIFY(QM,B,BCRO,PGPAR)

C----------------------Feed the vector to the solver-----------

                  CALL SOLVE(A,B,X,N,NBAND)

C----------------------Check for convergence-------------------

                  CALL CONVER(N,X,XM,MFLAG)

C----------------------Update Xm = X m+1 ----------------------

                  CALL UPDATE(N,X,XM)

C---------------Find flow component at iteration step---------

                  CALL FLOW(N,X,QM)

20          CONTINUE

C--------------Update h and q for next time level----------------------

            CALL UPDATE(N,X,X0)
            CALL FLOW(N,X,Q0)

C--------------Calculate sediment & polutants for 100 time steps along the---
C--------------simulation, each time using the average flow of the last------
C--------------NWRITE values in between, following this steps:---------------
C--------------a) Call sediment subroutine if there is inflow (change units--            
C--------------from q(m2/s)-->qsed(cm2/s))-----------------------------------
C--------------b) Call water quality subroutine if user selected (IWQ=1)-----
C--------------c) Write outputs to files ------------------------------------
            
            QTEMP=QTEMP+Q0(N)
            DO 25 J=1,3
                  ND=NODEX(J)
                  QSED(J)=QSED(J)+Q0(ND)
25          CONTINUE
            DO 30 I=1,100
                IF(LCOUNT.EQ.I*NWRITE) THEN
                  QTEMP=QTEMP/NWRITE
                  DO 27 J=1,3
                        QSED(J)=QSED(J)/NWRITE*10000.D0
                        IF(QSED(J).LT.0.D0)QSED(J)=DABS(QSED(J))
27                CONTINUE
                  QSED(4)=QTEMP*10000.D0
                  QIN=QK(1)*BCRO**(5.D0/3.D0)*10000.D0
                  QOUT=Q0(N)
                  CALL GRASSED(TIME,N,QIN,QOUT,NODEX,
     &                        ICOARSE,COARSE,FWIDTH,ISCR)
                  IF((BCRO.NE.0.D0).or.(QOUT.ne.0.d0)) THEN
                     IF(IWQ.GT.1) CALL WQSUB(IWQ,TIME,N)
                  ENDIF
                  TOLD=TIME
                  CALL KWWRITE(N,LCOUNT,M,QTEMP,X,BCROQ,FWIDTH)
                  QTEMP=0.D0
                  DO 28 J=1,4
                        QSED(J)=0.D0
28                CONTINUE                
                ENDIF
30          CONTINUE
40     CONTINUE
      
C--------------Write a summary of results at the end of the run --

      CALL OUTMASS(VL,FWIDTH,SWIDTH,SLENGTH,TRAI,LISFIL,ISCR)

      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)

c--for debug--gasubwt-----
c      close(19)
c--for debug--gasubwt-----

      STOP
      END
