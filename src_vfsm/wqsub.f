              SUBROUTINE WQSUB(IWQ,TIME,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C					                                                          C 
C  Water quality component skeleton                                           C 
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	PARAMETER (MAXEQN=1001,MAXBND=40)	
	IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

	
C-------Read in main parameters of the program--------------

	WRITE(18,*)IWQ,TIME,N
	
	RETURN
	END

