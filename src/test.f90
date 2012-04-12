      PROGRAM TEST
      USE DISLIN
      IMPLICIT NONE

      INTEGER, PARAMETER :: NPT=100,N1=10003,N2=20001
      REAL, DIMENSION (N1) :: XRAY,YRAY,ZRAY
      INTEGER, DIMENSION (N2) :: I1,I2,I3
      REAL :: FPI,STEP,X,Y,ZLEV
      INTEGER :: N,I,J,NTRI

!     Generating 10000 surface points XRAY(N),YRAY(N),ZRAY(N)
      FPI=3.1415927/180.
      STEP=360./(NPT-1)
      N=0
      DO I=1,NPT
        X=(I-1)*STEP
        DO J=1,NPT
          N=N+1
          Y=(J-1)*STEP
          XRAY(N)=X
          YRAY(N)=Y
          ZRAY(N)=2*SIN(X*FPI)*SIN(Y*FPI)
        END DO
      END DO

!     Calculating triangulation (I1,I2,I3) from X,Y and Z
!     NTRI is the returned number of calculated triangles
      CALL TRIANG(XRAY,YRAY,N,I1,I2,I3,2*N+1,NTRI)

!     Plotting contours of the triangulation on the screen
      CALL METAFL ('CONS')
      CALL DISINI ()

      CALL GRAF(0.,360.,0.,90.,0.,360.,0.,90.)

      DO I=0,8
        ZLEV=-2.+I*0.5
        CALL CONTRI (XRAY,YRAY,ZRAY,N,I1,I2,I3,NTRI,ZLEV)
      END DO
      CALL DISFIN
      END 
