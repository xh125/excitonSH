
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% THIS PROGRAM IS USED TO SIMULATE CHARGE TRANSPORT WITH THE SELF-CONSISTENT SURFACE HOPPING METHOD %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% ONE-DIMENSIONAL MOLECULAR STACK WITH LOCAL ELECTRON-PHONON COUPLINGS AND SYSTEM-BATH INTERACTIONS %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% BY THE WANG GROUP AT DEPARTMENT OF CHEMISTRY, ZHEJIANG UNIVERSITY; 2017/03/09; LJWANG@ZJU.EDU.CN  %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%!
!% MODULE %!
!%%%%%%%%%%!

MODULE PARAS
  IMPLICIT NONE

  INTEGER NSITE,NAVER,NSNAP,NSTEP,IAVER,ISNAP,ISTEP,ISURFACE,ICENTER
  DOUBLE PRECISION AU2CM,AU2EV,AU2J,AU2FS,AU2PS,AU2AMU,AU2ANG,KB,SQRT3,SQRT5,SQRT7
  DOUBLE PRECISION MASS,TEMP,GAMMA,DT,ALPHA,K,TAU
  DOUBLE PRECISION,ALLOCATABLE :: H0(:,:),HH(:,:)
  DOUBLE PRECISION,ALLOCATABLE :: X(:),V(:),E(:),P(:,:),D(:,:,:),G(:)
  DOUBLE PRECISION,ALLOCATABLE :: X0(:),V0(:),E0(:),P0(:,:),D0(:,:,:),G1(:)
  DOUBLE PRECISION,ALLOCATABLE :: PES(:,:,:),INF(:,:,:),CSIT(:,:),WSIT(:,:),PSIT(:,:),XSIT(:,:),KSIT(:,:),MSD(:),IPR(:),MSDS(:,:)
  COMPLEX*16 EYE
  COMPLEX*16,ALLOCATABLE :: C(:),W(:),W0(:)
ENDMODULE

!%%%%%%%%%%%!
!% PROGRAM %!
!%%%%%%%%%%%!

PROGRAM SFSH_LOCAL
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE
  DOUBLE PRECISION T0,T1
  DOUBLE PRECISION SUMG0,SUMG1,MINDE,FLAGD

  !%%%%%%%%%%%%%%%!
  !% PREPARATION %!
  !%%%%%%%%%%%%%%%!

  CALL CPU_TIME(T0)
  CALL SET_CONSTANTS()
  CALL READ_PARAMETERS()
  CALL TREAT_PARAMETERS()
  CALL INIT_RANDOM_SEED()

  ALLOCATE(HH(1:NSITE,1:NSITE))
  ALLOCATE(X(1:NSITE),V(1:NSITE),E(1:NSITE),P(1:NSITE,1:NSITE),D(1:NSITE,1:NSITE,1:NSITE),G(1:NSITE))
  ALLOCATE(X0(1:NSITE),V0(1:NSITE),E0(1:NSITE),P0(1:NSITE,1:NSITE),D0(1:NSITE,1:NSITE,1:NSITE),G1(1:NSITE))
  ALLOCATE(PES(0:NSITE,1:NSNAP,1:NAVER),INF(1:3,1:NSNAP,1:NAVER),CSIT(1:NSITE,1:NSNAP),&
            WSIT(1:NSITE,1:NSNAP),PSIT(1:NSITE,1:NSNAP),XSIT(1:NSITE,1:NSNAP),&
            KSIT(1:NSITE,1:NSNAP),MSD(1:NSNAP))
  ALLOCATE(IPR(1:NSNAP),MSDS(1:NSNAP,1:NAVER))
  ALLOCATE(C(1:NSITE),W(1:NSITE),W0(1:NSITE))

  !%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% LOOP OVER REALIZATIONS %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%!

  PES=0.0D0
  INF=0.0D0
  CSIT=0.0D0
  WSIT=0.0D0
  PSIT=0.0D0
  XSIT=0.0D0
  KSIT=0.0D0
  MSD=0.0D0
  IPR=0.0D0
  MSDS=0.0D0

  DO IAVER=1,NAVER
    WRITE(6,'(A,I4.4,A)') '###### IAVER=',IAVER,' ######'

    !%%%%%%%%%%%%%%%%%%!
    !% INITIALIZATION %!
    !%%%%%%%%%%%%%%%%%%!

    CALL INIT_COORDINATE_VELOCITY(X,V)
    CALL INIT_DYNAMICAL_VARIABLE(X,C,E,P,W)
    CALL CALCULATE_NONADIABATIC_COUPLING(E,P,D)
    X0=X; V0=V; E0=E; P0=P; D0=D; W0=W

    !%%%%%%%%%%%%%%%%%%%%%%%!
    !% LOOP OVER SNAPSHOTS %!
    !%%%%%%%%%%%%%%%%%%%%%%%!

    DO ISNAP=1,NSNAP
      DO ISTEP=1,NSTEP

        !%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% UPDATE X,V,C,E,P,D,W,G %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%!

        CALL RK4_NUCLEI(P0,X,V,DT)
        CALL RK4_ELECTRON_DIABATIC(X0,C,DT)
        CALL CALCULATE_EIGEN_ENERGY_STATE(X,E,P)
        CALL CALCULATE_NONADIABATIC_COUPLING(E,P,D)
        CALL CONVERT_DIABATIC_ADIABATIC(P,C,W)
        CALL CALCULATE_HOPPING_PROBABILITY(W0,V0,D0,DT,G,G1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% CALCULATE SUMG0,SUMG1,MINDE %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        SUMG0=(ABS(W0(ISURFACE))**2-ABS(W(ISURFACE))**2)/ABS(W0(ISURFACE))**2
        SUMG1=SUM(G1)
        IF(ISURFACE.EQ.1) THEN
          MINDE=(E0(ISURFACE+1)-E0(ISURFACE))*AU2EV
        ELSEIF(ISURFACE.EQ.NSITE) THEN
          MINDE=(E0(ISURFACE)-E0(ISURFACE-1))*AU2EV
        ELSEIF((E0(ISURFACE+1)-E0(ISURFACE)).LT.(E0(ISURFACE)-E0(ISURFACE-1))) THEN
          MINDE=(E0(ISURFACE+1)-E0(ISURFACE))*AU2EV
        ELSE
          MINDE=(E0(ISURFACE)-E0(ISURFACE-1))*AU2EV
        ENDIF

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% CHANGE POTENTIAL ENERGY SURFACE %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        IF(ISURFACE.EQ.1) THEN
          ISITE=ISURFACE+1
        ELSEIF(ISURFACE.EQ.NSITE) THEN
          ISITE=ISURFACE-1
        ELSEIF((E0(ISURFACE+1)-E0(ISURFACE)).LT.(E0(ISURFACE)-E0(ISURFACE-1))) THEN
          ISITE=ISURFACE+1
        ELSE
          ISITE=ISURFACE-1
        ENDIF
        G(ISITE)=SUMG0-(SUM(G1)-G1(ISITE))
        IF(G(ISITE).LT.0.0D0) G(ISITE)=0.0D0
        IF(SUM(G).GE.1.0D0) G=G/SUM(G)

        CALL NONADIABATIC_TRANSITION(E0,P0,D0,G,W,V,C)

        !%%%%%%%%%%%%%%%%%%%!
        !% ADD BATH EFFECT %!
        !%%%%%%%%%%%%%%%%%%%!

        CALL ADD_BATH_EFFECT(D0,P0,DT,X,V)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% RESET DYNAMICAL VARIABLE %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        X0=X; V0=V; E0=E; P0=P; D0=D; W0=W
      ENDDO

      !%%%%%%%%%%%%%%%%%%%%%!
      !% STORE INFORMATION %!
      !%%%%%%%%%%%%%%%%%%%%%!

      PES(0,ISNAP,IAVER)=E(ISURFACE)
      INF(1,ISNAP,IAVER)=SUMG0
      INF(2,ISNAP,IAVER)=SUMG1
      INF(3,ISNAP,IAVER)=MINDE
      FLAGD=0.0D0
      DO ISITE=1,NSITE
        CSIT(ISITE,ISNAP)=CSIT(ISITE,ISNAP)+ABS(C(ISITE))**2
        WSIT(ISITE,ISNAP)=WSIT(ISITE,ISNAP)+ABS(W(ISITE))**2
        PSIT(ISITE,ISNAP)=PSIT(ISITE,ISNAP)+P(ISITE,ISURFACE)**2
        XSIT(ISITE,ISNAP)=XSIT(ISITE,ISNAP)+X(ISITE)
        KSIT(ISITE,ISNAP)=KSIT(ISITE,ISNAP)+0.5D0*MASS*V(ISITE)**2
        PES(ISITE,ISNAP,IAVER)=E(ISITE)
        MSDS(ISNAP,IAVER)=MSDS(ISNAP,IAVER)+P(ISITE,ISURFACE)**2*(ISITE-ICENTER)**2
        FLAGD=FLAGD+P(ISITE,ISURFACE)**4
      ENDDO
      IPR(ISNAP)=IPR(ISNAP)+1/FLAGD
      MSD(ISNAP)=MSD(ISNAP)+MSDS(ISNAP,IAVER)
    ENDDO
  ENDDO
  CSIT=CSIT/NAVER
  WSIT=WSIT/NAVER
  PSIT=PSIT/NAVER
  XSIT=XSIT/NAVER*AU2ANG
  KSIT=KSIT/NAVER*AU2EV
  MSD=MSD/NAVER
  IPR=IPR/NAVER

  !%%%%%%%%%%%%%%%%%%%%!
  !% SAVE INFORMATION %!
  !%%%%%%%%%%%%%%%%%%%%!

  OPEN(UNIT=99,FILE='PES')
  DO IAVER=1,1
    DO ISNAP=1,NSNAP
      WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(PES(ISITE,ISNAP,IAVER),ISITE=0,NSITE)
    ENDDO
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='INF')
  DO IAVER=1,1
    DO ISNAP=1,NSNAP
      WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(INF(ISITE,ISNAP,IAVER),ISITE=1,3)
    ENDDO
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='CSIT')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(CSIT(ISITE,ISNAP),ISITE=1,NSITE)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='WSIT')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(WSIT(ISITE,ISNAP),ISITE=1,NSITE)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='PSIT')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(PSIT(ISITE,ISNAP),ISITE=1,NSITE)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='XSIT')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(XSIT(ISITE,ISNAP),ISITE=1,NSITE)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='KSIT')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(KSIT(ISITE,ISNAP),ISITE=1,NSITE),KB*TEMP/2*AU2EV
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='MSD')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,MSD(ISNAP)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='IPR')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,IPR(ISNAP)
  ENDDO
  CLOSE(99)

  OPEN(UNIT=99,FILE='MSDS')
  DO ISNAP=1,NSNAP
    WRITE(99,'(9999E12.5)') DT*NSTEP*ISNAP*AU2FS,(MSDS(ISNAP,IAVER),IAVER=1,NAVER)
  ENDDO
  CLOSE(99)

  CALL CPU_TIME(T1)
  WRITE(6,'(A,F10.2,A)') 'TOTAL TIME IS',(T1-T0)/3600,'HOURS'
ENDPROGRAM





!###############################################################!
!# BBBB    A    SSSS IIIII  CCCC  CCCC  OOO  DDDD  EEEEE  SSSS #!
!# B   B  A A  S       I   C     C     O   O D   D E     S     #!
!# BBBB   AAA   SSS    I   C     C     O   O D   D EEEEE  SSS  #!
!# B   B A   A     S   I   C     C     O   O D   D E         S #!
!# BBBB  A   A SSSS  IIIII  CCCC  CCCC  OOO  DDDD  EEEEE SSSS  #!
!###############################################################!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% SET PHYSICAL CONSTANTS                              %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_UNIT       %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_MASS_UNIT  %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/PHYSICAL_CONSTANT %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE SET_CONSTANTS()
  USE PARAS
  IMPLICIT NONE

  AU2CM=2.194887656D5
  AU2EV=2.7211D1
  AU2J=4.35974417D-18
  AU2FS=2.418884326505D-2
  AU2PS=2.418884326505D-5
  AU2AMU=5.4858D-4
  AU2ANG=5.291772108D-1
  KB=1.3806504D-23
  EYE=(0.0D0,1.0D0)
  SQRT3=DSQRT(3.0D0)
  SQRT5=DSQRT(5.0D0)
  SQRT7=DSQRT(7.0D0)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% READ PARAMETERS FROM INPUT %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE READ_PARAMETERS()
  USE PARAS
  IMPLICIT NONE

  OPEN(UNIT=99,FILE='INPUT-SSH')
  READ(99,*) NSITE
  READ(99,*) TAU
  READ(99,*) ALPHA
  READ(99,*) MASS
  READ(99,*) K
  READ(99,*) TEMP
  READ(99,*) GAMMA
  READ(99,*) DT
  READ(99,*) NSTEP
  READ(99,*) NSNAP
  READ(99,*) NAVER
  CLOSE(99)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%!
!% TREAT PARAMETERS %!
!%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE TREAT_PARAMETERS()
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE

  ALLOCATE(H0(1:NSITE,1:NSITE))
  H0=0.0D0
  DO ISITE=1,NSITE-1
    H0(ISITE,ISITE+1)=TAU
    H0(ISITE+1,ISITE)=TAU
  ENDDO

  TAU=TAU/AU2CM
  H0=H0/AU2CM
  ALPHA=ALPHA/AU2CM*AU2ANG
  MASS=MASS/AU2AMU
  GAMMA=GAMMA*AU2PS
  K=K/AU2AMU*AU2PS**2
  DT=DT/AU2FS
  KB=KB/AU2J
  ICENTER=ANINT(NSITE/2.0D0)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% INIT COORDINATE AND VELOCITIE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE INIT_COORDINATE_VELOCITY(XX,VV)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE
  DOUBLE PRECISION XX(1:NSITE),VV(1:NSITE),GAUSSIAN_RANDOM_NUMBER
  EXTERNAL GAUSSIAN_RANDOM_NUMBER

  DO ISITE=1,NSITE
    XX(ISITE)=GAUSSIAN_RANDOM_NUMBER(0.0D0,DSQRT(KB*TEMP/K))
    VV(ISITE)=GAUSSIAN_RANDOM_NUMBER(0.0D0,DSQRT(KB*TEMP/MASS))
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE EIGENENERGY AND EIGENSTATE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE CALCULATE_EIGEN_ENERGY_STATE(XX,EE,PP)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,IERR
  DOUBLE PRECISION XX(1:NSITE),EE(1:NSITE),PP(1:NSITE,1:NSITE),FV1(1:NSITE),FV2(1:NSITE)

  HH=H0
  DO ISITE=1,NSITE
    HH(ISITE,ISITE)=HH(ISITE,ISITE)+ALPHA*XX(ISITE)
  ENDDO
  CALL RS(NSITE,NSITE,HH,EE,1,PP,FV1,FV2,IERR)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CONVERT WAVEFUNCTION FROM DIABATIX TO ADIABATIC BASIS %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE CONVERT_DIABATIC_ADIABATIC(PP,CC,WW)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE
  DOUBLE PRECISION PP(1:NSITE,1:NSITE)
  COMPLEX*16 CC(1:NSITE),WW(1:NSITE)

  WW=0.0D0
  DO ISITE=1,NSITE
    DO JSITE=1,NSITE
      WW(ISITE)=WW(ISITE)+PP(JSITE,ISITE)*CC(JSITE)
    ENDDO
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% INIT DYNAMICAL VARIABLE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE INIT_DYNAMICAL_VARIABLE(XX,CC,EE,PP,WW)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE
  DOUBLE PRECISION XX(1:NSITE),EE(1:NSITE),PP(1:NSITE,1:NSITE),FLAGR,FLAGD
  COMPLEX*16 CC(1:NSITE),WW(1:NSITE)
  
  CC=0.0D0; CC(ICENTER)=1.0D0
  CALL CALCULATE_EIGEN_ENERGY_STATE(XX,EE,PP)
  CALL CONVERT_DIABATIC_ADIABATIC(PP,CC,WW)

  CALL RANDOM_NUMBER(FLAGR)
  FLAGD=0.0D0
  DO ISITE=1,NSITE
    FLAGD=FLAGD+PP(ICENTER,ISITE)**2
    IF(FLAGR.LE.FLAGD) THEN
      ISURFACE=ISITE
      EXIT
    ENDIF
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE NONADIABATIC COUPLING %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 630          %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE CALCULATE_NONADIABATIC_COUPLING(EE,PP,DD)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE,KSITE
  DOUBLE PRECISION EE(1:NSITE),PP(1:NSITE,1:NSITE),DD(1:NSITE,1:NSITE,1:NSITE)

  DD=0.0D0
  DO ISITE=1,NSITE
  DO JSITE=1,NSITE
  IF(JSITE.NE.ISITE) THEN
    DO KSITE=1,NSITE
      DD(ISITE,JSITE,KSITE)=ALPHA*PP(KSITE,ISITE)*PP(KSITE,JSITE)/(EE(JSITE)-EE(ISITE))
    ENDDO
  ENDIF
  ENDDO
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE DERIVATIVE OF COORDINATE AND VELOCITIE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 630                           %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE DERIVS_NUCLEI(PP,XX,VV,DX,DV)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE
  DOUBLE PRECISION PP(1:NSITE,1:NSITE),XX(1:NSITE),VV(1:NSITE),DX(1:NSITE),DV(1:NSITE)

  DO ISITE=1,NSITE
    DV(ISITE)=(-K*XX(ISITE)-ALPHA*PP(ISITE,ISURFACE)**2)/MASS-GAMMA*VV(ISITE)
    DX(ISITE)=VV(ISITE)
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% RK4 METHOD TO OBTAIN COORDINATE AND VELOCITIE AFTER A TIME INTERVAL %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/RUNGE_KUTTA_METHODS               %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE RK4_NUCLEI(PP,XX,VV,TT)
  USE PARAS
  IMPLICIT NONE

  DOUBLE PRECISION PP(1:NSITE,1:NSITE),TT,TT2,TT6
  DOUBLE PRECISION XX(1:NSITE),XX0(1:NSITE),DX1(1:NSITE),DX2(1:NSITE),DX3(1:NSITE),DX4(1:NSITE)
  DOUBLE PRECISION VV(1:NSITE),VV0(1:NSITE),DV1(1:NSITE),DV2(1:NSITE),DV3(1:NSITE),DV4(1:NSITE)

  TT2=TT/2.0D0; TT6=TT/6.0D0
  CALL DERIVS_NUCLEI(PP,XX,VV,DX1,DV1)
  XX0=XX+TT2*DX1; VV0=VV+TT2*DV1
  CALL DERIVS_NUCLEI(PP,XX0,VV0,DX2,DV2)
  XX0=XX+TT2*DX2; VV0=VV+TT2*DV2
  CALL DERIVS_NUCLEI(PP,XX0,VV0,DX3,DV3)
  XX0=XX+TT*DX3; VV0=VV+TT*DV3
  CALL DERIVS_NUCLEI(PP,XX0,VV0,DX4,DV4)
  XX=XX+TT6*(DX1+2.0D0*DX2+2.0D0*DX3+DX4)
  VV=VV+TT6*(DV1+2.0D0*DV2+2.0D0*DV3+DV4)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% ADD BATH EFFECT TO COORDINATE AND VELOCITIE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 462 AND 638              %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE ADD_BATH_EFFECT(DD,PP,TT,XX,VV)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE
  DOUBLE PRECISION DD(1:NSITE,1:NSITE,1:NSITE),PP(1:NSITE,1:NSITE),TT,XX(1:NSITE),VV(1:NSITE)
  DOUBLE PRECISION SIGMAR,KK,R1,R2,R3,R4,Z1,Z2,Z3,Z4,GAUSSIAN_RANDOM_NUMBER_FAST
  EXTERNAL GAUSSIAN_RANDOM_NUMBER_FAST

  SIGMAR=DSQRT(2.0D0*KB*TEMP*GAMMA*TT/MASS)
  DO ISITE=1,NSITE
    KK=K
    DO JSITE=1,NSITE
      IF(JSITE.NE.ISURFACE) KK=KK+2.0D0*ALPHA*DD(JSITE,ISURFACE,ISITE)*PP(ISITE,ISURFACE)*PP(ISITE,JSITE)
    ENDDO

    R1=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
    R2=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
    R3=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
    R4=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
    Z1=R1
    Z2=TT*(R1/2.0D0+R2/SQRT3/2.0D0)
    Z3=TT**2*(R1/6.0D0+R2*SQRT3/12.0D0+R3/SQRT5/12.0D0)
    Z4=TT**3*(R1/24.0D0+R2*SQRT3/40.0D0+R3/SQRT5/24.0D0+R4/SQRT7/120.0D0)
    XX(ISITE)=XX(ISITE)+(Z2-GAMMA*Z3+(-KK/MASS+GAMMA**2)*Z4)
    VV(ISITE)=VV(ISITE)+(Z1-GAMMA*Z2+(-KK/MASS+GAMMA**2)*Z3+(2.0D0*GAMMA*KK/MASS-GAMMA**3)*Z4)
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE DERIVATIVE OF WAVEFUNCTION %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 631               %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE DERIVS_ELECTRON_DIABATIC(XX,CC,DC)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE
  DOUBLE PRECISION XX(1:NSITE)
  COMPLEX*16 CC(1:NSITE),DC(1:NSITE)

  DO ISITE=1,NSITE
    DC(ISITE)=ALPHA*XX(ISITE)*CC(ISITE)
    DO JSITE=1,NSITE
      IF(ABS(JSITE-ISITE).EQ.1) DC(ISITE)=DC(ISITE)+TAU*CC(JSITE)
    ENDDO
    DC(ISITE)=DC(ISITE)*(-EYE)
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% RK4 METHOD TO OBTAIN WAVEFUNCTION AFTER A TIME INTERVAL %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/RUNGE_KUTTA_METHODS   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE RK4_ELECTRON_DIABATIC(XX,CC,TT)
  USE PARAS
  IMPLICIT NONE

  DOUBLE PRECISION TT,TT2,TT6
  DOUBLE PRECISION XX(1:NSITE)
  COMPLEX*16 CC(1:NSITE),CC0(1:NSITE),DC1(1:NSITE),DC2(1:NSITE),DC3(1:NSITE),DC4(1:NSITE)

  TT2=TT/2.0D0; TT6=TT/6.0D0
  CALL DERIVS_ELECTRON_DIABATIC(XX,CC,DC1)
  CC0=CC+TT2*DC1
  CALL DERIVS_ELECTRON_DIABATIC(XX,CC0,DC2)
  CC0=CC+TT2*DC2
  CALL DERIVS_ELECTRON_DIABATIC(XX,CC0,DC3)
  CC0=CC+TT*DC3
  CALL DERIVS_ELECTRON_DIABATIC(XX,CC0,DC4)
  CC=CC+TT6*(DC1+2.0D0*DC2+2.0D0*DC3+DC4)
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE HOPPING PROBABILITY %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 631        %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE CALCULATE_HOPPING_PROBABILITY(WW,VV,DD,TT,GG,GG1)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE
  DOUBLE PRECISION VV(1:NSITE),DD(1:NSITE,1:NSITE,1:NSITE),GG(1:NSITE),GG1(1:NSITE),TT,SUMVD
  COMPLEX*16 WW(1:NSITE)

  GG=0.0D0
  GG1=0.0D0
  DO ISITE=1,NSITE
  IF(ISITE.NE.ISURFACE) THEN
    SUMVD=0.0D0
    DO JSITE=1,NSITE
      SUMVD=SUMVD+VV(JSITE)*DD(ISURFACE,ISITE,JSITE)
    ENDDO
    GG(ISITE)=2.0D0*TT*REAL(CONJG(WW(ISURFACE))*WW(ISITE))*SUMVD/REAL(CONJG(WW(ISURFACE))*WW(ISURFACE))
    GG1(ISITE)=GG(ISITE)
    IF(GG(ISITE).LT.0.0D0) GG(ISITE)=0.0D0
  ENDIF
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% NONADIABATIC TRANSITION %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: NOTEBOOK PAGE 635  %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE NONADIABATIC_TRANSITION(EE,PP,DD,GG,WW,VV,CC)
  USE PARAS
  IMPLICIT NONE

  INTEGER ISITE,JSITE
  DOUBLE PRECISION EE(1:NSITE),PP(1:NSITE,1:NSITE),DD(1:NSITE,1:NSITE,1:NSITE),GG(1:NSITE),VV(1:NSITE),SUMVD,SUMDD,SUMGG,FLAGR,FLAGD
  COMPLEX*16 WW(1:NSITE),CC(1:NSITE)

  CALL MORE_RANDOM()
  CALL RANDOM_NUMBER(FLAGR)
  SUMGG=0.0D0
  DO ISITE=1,NSITE
  IF(ISITE.NE.ISURFACE) THEN
    SUMGG=SUMGG+GG(ISITE)
    IF(FLAGR.LT.SUMGG) THEN
      SUMVD=0.0D0
      SUMDD=0.0D0
      DO JSITE=1,NSITE
        SUMVD=SUMVD+VV(JSITE)*DD(ISURFACE,ISITE,JSITE)
        SUMDD=SUMDD+DD(ISURFACE,ISITE,JSITE)**2
      ENDDO
      FLAGD=1.0D0+2.0D0*(EE(ISURFACE)-EE(ISITE))*SUMDD/MASS/SUMVD**2

      IF(FLAGD.GE.0.0D0) THEN
        FLAGD=SUMVD/SUMDD*(-1.0D0+DSQRT(FLAGD))
        DO JSITE=1,NSITE
          VV(JSITE)=VV(JSITE)+FLAGD*DD(ISURFACE,ISITE,JSITE)
        ENDDO
        ISURFACE=ISITE
      ENDIF

      EXIT
    ENDIF
  ENDIF
  ENDDO
ENDSUBROUTINE





!###########################################################################!
!#  GGG  EEEEE N   N EEEEE RRRR    A   L      CCCC  OOO  DDDD  EEEEE  SSSS #!
!# G     E     NN  N E     R   R  A A  L     C     O   O D   D E     S     #!
!# GGGG  EEEEE N N N EEEEE RRRR   AAA  L     C     O   O D   D EEEEE  SSS  #!
!# G   G E     N  NN E     R  R  A   A L     C     O   O D   D E         S #!
!#  GGG  EEEEE N   N EEEEE R   R A   A LLLLL  CCCC  OOO  DDDD  EEEEE SSSS  #!
!###########################################################################!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% A BETTER VERSION OF RANDOM_SEED %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE INIT_RANDOM_SEED()
  IMPLICIT NONE

  INTEGER II,NN,VALUE(1:8)
  INTEGER,ALLOCATABLE :: SEED(:)
  DOUBLE PRECISION FLAGD

  CALL RANDOM_SEED(SIZE=NN)
  ALLOCATE(SEED(NN))
  CALL DATE_AND_TIME(VALUES=VALUE)
  SEED = VALUE(8)+37*(/(II-1,II=1,NN)/)
  CALL RANDOM_SEED(PUT=SEED)
  DEALLOCATE(SEED)

  DO II=1,VALUE(6)*3600+VALUE(7)*60+VALUE(8)
    CALL RANDOM_NUMBER(FLAGD)
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% MAKE RANDOM NUMBER MORE RANDOM %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE MORE_RANDOM()
  IMPLICIT NONE

  INTEGER II,VALUE(1:8)
  DOUBLE PRECISION FLAGD

  CALL DATE_AND_TIME(VALUES=VALUE)
  DO II=1,VALUE(8)/100
    CALL RANDOM_NUMBER(FLAGD)
  ENDDO
ENDSUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% GAUSSIAN RANDOM NUMBER GENERATOR USING BOX-MULLER METHOD   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/GAUSSIAN_RANDOM_VARIABLE %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

FUNCTION GAUSSIAN_RANDOM_NUMBER(MEAN,SIGMA)
  IMPLICIT NONE

  DOUBLE PRECISION GAUSSIAN_RANDOM_NUMBER,MEAN,SIGMA,PI,R1,R2

  PI=4.0D0*DATAN(1.0D0)
  CALL MORE_RANDOM()
  CALL RANDOM_NUMBER(R1)
  CALL MORE_RANDOM()
  CALL RANDOM_NUMBER(R2)
  GAUSSIAN_RANDOM_NUMBER=MEAN+SIGMA*DSQRT(-2.0D0*DLOG(R1))*DCOS(2.0D0*PI*R2)
ENDFUNCTION

FUNCTION GAUSSIAN_RANDOM_NUMBER_FAST(MEAN,SIGMA)
  IMPLICIT NONE

  DOUBLE PRECISION GAUSSIAN_RANDOM_NUMBER_FAST,MEAN,SIGMA,PI,R1,R2

  PI=4.0D0*DATAN(1.0D0)
  CALL RANDOM_NUMBER(R1)
  CALL RANDOM_NUMBER(R2)
  GAUSSIAN_RANDOM_NUMBER_FAST=MEAN+SIGMA*DSQRT(-2.0D0*DLOG(R1))*DCOS(2.0D0*PI*R2)
ENDFUNCTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% DIAGONALIZE A SYMMETRIC MATRIX USING RS   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% THIS IS A REARRANGED VERSION OF DSUB1.F90 %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)

  IF(N.LE.NM) GOTO 10
  IERR=10*N
  GOTO 50

  10 IF(MATZ.NE.0) GOTO 20
  CALL TRED1(NM,N,A,W,FV1,FV2)
  CALL TQLRAT(N,W,FV2,IERR)
  GOTO 50

  20 CALL TRED2(NM,N,A,W,FV1,Z)
  CALL TQL2(NM,N,W,FV1,Z,IERR)
  50 RETURN
ENDSUBROUTINE

SUBROUTINE TRED2(NM,N,A,D,E,Z)
  IMPLICIT REAL*8 (A-H,O-Z) 
  DIMENSION A(NM,N),D(N),E(N),Z(NM,N)

  DO 100 I=1,N
  DO 100 J=1,I
  Z(I,J)=A(I,J)
  100 CONTINUE
  IF(N.EQ.1) GOTO 320
  DO 300 II=2,N
  I=N+2-II
  L=I-1
  H=0.0D0
  SCALE=0.0D0
  IF(L.LT.2) GOTO 130
  DO 120 K=1,L
  120 SCALE=SCALE+DABS(Z(I,K))
  IF(SCALE.NE.0.0D0) GOTO 1120
  130 E(I)=Z(I,L)
  GOTO 290
  1120 DO 150 K=1,L
  Z(I,K)=Z(I,K)/SCALE
  H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
  F=Z(I,L)
  G=-DSIGN(DSQRT(H),F)
  E(I)=SCALE*G
  H=H-F*G
  Z(I,L)=F-G
  F=0.0D0
  DO 240 J=1,L
  Z(J,I)=Z(I,J)/(SCALE*H)
  G=0.0D0
  DO 180 K=1,J
  180 G=G+Z(J,K)*Z(I,K)
  JP1=J+1
  IF(L.LT.JP1) GOTO 220
  DO 200 K=JP1,L
  200 G=G+Z(K,J)*Z(I,K)
  220 E(J)=G/H
  F=F+E(J)*Z(I,J)
  240 CONTINUE
  HH=F/(H+H)
  DO 260 J=1,L
  F=Z(I,J)
  G=E(J)-HH*F
  E(J)=G
  DO 260 K=1,J
  Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
  DO 280 K=1,L
  280 Z(I,K)=SCALE*Z(I,K)
  290 D(I)=H
  300 CONTINUE
  320 D(1)=0.0D0
  E(1)=0.0D0
  DO 500 I=1,N
  L=I-1
  IF(D(I).EQ.0.0D0) GOTO 380
  DO 360 J=1,L
  G=0.0D0
  DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
  DO 360 K=1,L
  Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
  Z(I,I)=1.0D0
  IF(L.LT.1) GOTO 500
  DO 400 J=1,L
  Z(I,J)=0.0D0
  Z(J,I)=0.0D0
  400 CONTINUE
  500 CONTINUE
  RETURN
ENDSUBROUTINE

SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8 MACHEP
  DIMENSION D(N),E(N),Z(NM,N)

  MACHEP=1.0D-30
  IERR=0

  ONE=1.0D0
  TWO=2.0D0
  ZERO=0.0D0

  IF(N.EQ.1) GOTO 1001
  DO 100 I=2,N
  100 E(I-1)=E(I)
  F=ZERO
  B=ZERO
  E(N)=ZERO
  DO 240 L=1,N
  J=0
  H=MACHEP*(DABS(D(L))+DABS(E(L)))
  IF(B.LT.H)B=H
  DO 110 M=L,N
  IF(DABS(E(M)).LE.B) GOTO 120
  110 CONTINUE
  120 IF(M.EQ.L) GOTO 220
  130 IF(J.EQ.30) GOTO 1000
  J=J+1
  L1=L+1
  G=D(L)
  P=(D(L1)-G)/(TWO*E(L))
  R=DSQRT(P*P+ONE)
  D(L)=E(L)/(P+DSIGN(R,P))
  H=G-D(L)
  DO 1120 I=L1,N
  1120 D(I)=D(I)-H
  F=F+H
  P=D(M)
  C=ONE
  S=ZERO
  MML=M-L
  DO 200 II=1,MML
  I=M-II
  G=C*E(I)
  H=C*P
  IF(DABS(P).LT.DABS(E(I))) GOTO 150
  C=E(I)/P
  R=SQRT(C*C+ONE)
  E(I+1)=S*P*R
  S=C/R
  C=ONE/R
  GO TO 160
  150 C=P/E(I)
  R=DSQRT(C*C+ONE)
  E(I+1)=S*E(I)*R
  S=ONE/R
  C=C*S
  160 P=C*D(I)-S*G
  D(I+1)=H+S*(C*G+S*D(I))
  DO 180 K=1,N
  H=Z(K,I+1)
  Z(K,I+1)=S*Z(K,I)+C*H
  Z(K,I)=C*Z(K,I)-S*H
  180 CONTINUE
  200 CONTINUE
  E(L)=S*P
  D(L)=C*P
  IF(DABS(E(L)).GT.B) GOTO 130
  220 D(L)=D(L)+F
  240 CONTINUE
  DO 300 II=2,N
  I=II-1
  K=I
  P=D(I)
  DO 260 J=II,N
  IF(D(J).GE.P) GOTO 260
  K=J
  P=D(J)
  260 CONTINUE
  IF(K.EQ.I) GOTO 300
  D(K)=D(I)
  D(I)=P
  DO 280 J=1,N
  P=Z(J,I)
  Z(J,I)=Z(J,K)
  Z(J,K)=P
  280 CONTINUE
  300 CONTINUE
  GOTO 1001
  1000 IERR=L
  1001 RETURN
ENDSUBROUTINE

REAL*8 FUNCTION SSUM(N,A,NSTP)
  REAL*8 A
  DIMENSION A(N)

  SSUM=A(1)
  IF(N.LE.1) GOTO 100
  DO I=2,N
  SSUM=SSUM+A(I)
  ENDDO
  100 CONTINUE
ENDFUNCTION

SUBROUTINE TRED1(NM,N,A,D,E,E2)
  IMPLICIT REAL*8 (A-H,O-Z) 
  DIMENSION A(NM,N),D(N),E(N),E2(N)

  ZERO=0.0D0
  ONE=1.0D0
  TWO=2.0D0

  DO 100 I=1,N
  100 D(I)=A(I,I)

  DO 300 II=1,N
  I=N+1-II
  L=I-1
  H=ZERO
  SCALE=ZERO
  IF(L.LT.1) GOTO 130
  DO 120 K=1,L
  120 SCALE=SCALE+DABS(A(I,K))
  IF(SCALE.NE.ZERO) GOTO 140
  130 E(I)=ZERO
  E2(I)=ZERO
  GOTO 290
  140 DO 150 K=1,L
  A(I,K)=A(I,K)/SCALE
  150 H=H+A(I,K)*A(I,K)
  E2(I)=SCALE*SCALE*H
  F=A(I,L)
  G=-DSIGN(DSQRT(H),F)
  E(I)=SCALE*G
  H=H-F*G
  A(I,L)=F-G
  IF(L.EQ.1) GOTO 270
  F=ZERO
  DO 240 J=1,L
  G=ZERO
  DO 180 K=1,J
  180 G=G+A(J,K)*A(I,K)
  JP1=J+1
  IF(L.LT.JP1) GOTO 220
  DO 200 K=JP1,L
  200 G=G+A(K,J)*A(I,K)
  220 E(J)=G/H
  F=F+E(J)*A(I,J)
  240 CONTINUE

  H=F/(H+H)
  DO 260 J=1,L
  F=A(I,J)
  G=E(J)-H*F
  E(J)=G
  DO 260 K=1,J
  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
  260 CONTINUE
  270 DO 280 K=1,L
  280 A(I,K)=SCALE*A(I,K)
  290 H=D(I)
  D(I)=A(I,I)
  A(I,I)=H
  300 CONTINUE
  RETURN
ENDSUBROUTINE

SUBROUTINE TQLRAT(N,D,E2,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION D(N),E2(N)

  ONE=1.0D0
  TWO=2.0D0
  ZERO=0.0D0

  MACHEP=1.0D-30
  IERR=0
  IF(N.EQ.1) GOTO 1001
  DO 100 I=2,N
  100 E2(I-1)=E2(I)
  F=ZERO
  B=ZERO
  E2(N)=ZERO
  DO 290 L=1,N
  J=0
  H=MACHEP*(DABS(D(L))+DSQRT(E2(L)))
  IF(B.GT.H) GOTO 105
  B=H
  C=B*B
  105 DO 110 M=L,N
  IF(E2(M).LE.C) GOTO 120
  110 CONTINUE
  120 IF(M.EQ.L) GOTO 210
  130 IF(J.EQ.30) GOTO 1000
  J=J+1
  L1=L+1
  S=DSQRT(E2(L))
  G=D(L)
  P=(D(L1)-G)/(TWO*S)
  R=SQRT(P*P+ONE)
  D(L)=S/(P+DSIGN(R,P))
  H=G-D(L)
  DO 140 I=L1,N
  140 D(I)=D(I)-H
  F=F+H
  G=D(M)
  IF(G.EQ.ZERO) G=B
  H=G
  S=ZERO
  MML=M-L
  DO 200 II=1,MML
  I=M-II
  P=G*H
  R=P+E2(I)
  E2(I+1)=S*R
  S=E2(I)/R
  D(I+1)=H+S*(H+D(I))
  G=D(I)-E2(I)/G
  IF(G.EQ.ZERO)G=B
  H=G*P/R
  200 CONTINUE
  E2(L)=S*G
  D(L)=H
  IF(H.EQ.ZERO) GOTO 210
  IF(DABS(E2(L)).LE.DABS(C/H)) GOTO 210
  E2(L)=H*E2(L)
  IF(E2(L).NE.ZERO) GOTO 130
  210 P=D(L)+F

  IF(L.EQ.1) GOTO 250
  DO 230 II=2,L
  I=L+2-II
  IF(P.GE.D(I-1)) GOTO 270
  D(I)=D(I-1)
  230 CONTINUE
  250 I=1
  270 D(I)=P
  290 CONTINUE
  GOTO 1001
  1000 IERR=L
  1001 RETURN
ENDSUBROUTINE

SUBROUTINE RSP(NM,N,A,W,FV1,Z,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION A(NM,2),W(NM),FV1(NM),Z(NM,N)

  ZERO=0.0D0
  ONE=1.0D0

  DO 100 I=1,N
  DO 50 J=1,N
  Z(I,J)=ZERO
  50 CONTINUE
  Z(I,I)=ONE
  W(I)=A(I,2)
  FV1(I)=A(I,1)
  100 CONTINUE
  CALL IMTQL2(NM,N,W,FV1,Z,IERR)
ENDSUBROUTINE

SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8 MACHEP
  DIMENSION D(N),E(N),Z(NM,N)

  ZERO=0.0D0
  ONE=1.0D0
  TWO=2.0D0

  MACHEP=DBLE(2.**(-37))
  IERR=0
  IF(N.EQ.1) GOTO 1001
  DO 100 I=2,N
  100 E(I-1)=E(I)
  E(N)=ZERO
  DO 240 L=1,N
  J=0
  105 DO 110 M=L,N
  IF(M.EQ.N) GOTO 120
  IF(DABS(E(M)).LE.MACHEP*(DABS(D(M))+DABS(D(M+1)))) GOTO 120
  110 CONTINUE
  120 P=D(L)
  IF(M.EQ.L) GOTO 240
  IF(J.EQ.30) GOTO 1000
  J=J+1
  G=(D(L+1)-P)/(TWO*E(L))
  R=SQRT(G*G+ONE)
  G=D(M)-P+E(L)/(G+DSIGN(R,G))
  S=ONE
  C=ONE
  P=ZERO
  MML=M-L
  DO 200 II=1,MML
  I=M-II
  F=S*E(I)
  B=C*E(I)
  IF(DABS(F).LT.DABS(G)) GOTO 150
  C=G/F
  R=DSQRT(C*C+ONE)
  E(I+1)=F*R
  S=ONE/R
  C=C*S
  GOTO 160
  150 S=F/G
  R=DSQRT(S*S+ONE)
  E(I+1)=G*R
  C=ONE/R
  S=S*C
  160 G=D(I+1)-P
  R=(D(I)-G)*S+2.0*C*B
  P=S*R
  D(I+1)=G+P
  G=C*R-B
  DO 180 K=1,N
  F=Z(K,I+1)
  Z(K,I+1)=S*Z(K,I)+C*F 
  Z(K,I)=C*Z(K,I)-S*F
  180 CONTINUE
  200 CONTINUE
  D(L)=D(L)-P
  E(L)=G
  E(M)=ZERO
  GOTO 105
  240 CONTINUE
  DO 300 II=2,N
  I=II-1
  K=I
  P=D(I)
  DO 260 J=II,N
  IF(D(J).GE.P) GOTO 260
  K=J
  P=D(J)
  260 CONTINUE
  IF(K.EQ.I) GOTO 300
  D(K)=D(I)
  D(I)=P
  DO 280 J=1,N
  P=Z(J,I)
  Z(J,I)=Z(J,K)
  Z(J,K)=P
  280 CONTINUE
  300 CONTINUE             
  GOTO 1001
  1000 IERR=L
  1001 RETURN
ENDSUBROUTINE

SUBROUTINE RST(NM,N,W,E,Z,IERR)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION W(N),E(N),Z(NM,N)

  ZERO=0.0D0
  ONE=1.0D0
  TWO=2.0D0

  IF(N.LE.NM) GOTO 10
  IERR=10*N
  GOTO 50
  10 DO 40 I=1,N
  DO 30 J=1,N
  Z(J,I)=ZERO
  30 CONTINUE
  Z(I,I)=ONE
  40 CONTINUE
  CALL IMTQL2(NM,N,W,E,Z,IERR)
  50 RETURN
ENDSUBROUTINE

