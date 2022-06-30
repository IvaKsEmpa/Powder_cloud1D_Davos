  module precisions
  implicit none
  INTEGER, parameter :: dp=kind(1.0d0)
  end module precisions
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  MODULE GlobalParam
    USE precisions
    IMPLICIT NONE
    INTEGER::  ImpreE, ImpreF, Nv_Prim, argunit = 6
      INTEGER:: isave, penteX
    INTEGER  ::  Nx, iterfinal, cond_lim
     INTEGER  :: method_source_term
    INTEGER  ::  H_pv, U_pv,  Phi_pv, Hp_pv
    REAL (KIND = DP)  ::  X0, H_0, amplitude,  frottcoeff, disscoeff
    REAL (KIND = DP) ::  CFL, TIMEOUT, period_time, pi, angle, g, phi2
    REAL (KIND = DP), PARAMETER ::  EPS = 1.d-8, rhoAir = 1.3d0, rhoSnow =10.d0
    REAL (KIND = DP), PARAMETER :: injection_speed = 35.d0,  injection_density = 7.d0
    REAL (KIND = DP), PARAMETER :: rhoPowder =10.D0 
  END MODULE GlobalParam
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      

  PROGRAM code2D
  USE precisions
  use GlobalParam
    IMPLICIT NONE  
    INTEGER :: iv, I, IT, ix
    REAL (KIND = DP), ALLOCATABLE :: Ein(:), Sound_ax(:), U1(:), Qaverage(:)
    REAL (KIND = DP), ALLOCATABLE ::  Pression(:)
    REAL (KIND = DP), ALLOCATABLE :: CONS(:,:), FLUX(:,:), Prim(:,:), pente(:,:)
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:):: MaxVP, MinVp, X
    REAL (KIND = DP)  :: T1_CPU, T2_CPU, TIME,  TIME2
    REAL (KIND = DP) :: Lx, DX, DT, dt2
    REAL (KIND = DP) :: Hstar, Pstar, Ustar, Estar, Cxmax
    REAL(KIND=DP) :: xmax, XMIN, XMAX1, UMAX
    
    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:)    :: NamesOfPrimVar
  !----------------------------------------------------------------------------------------------------
    pi = 4.0d0*ATAN(1.0d0)
    Nv_Prim  = 4
    H_pv     = 1
    U_pv     = 2
    Phi_pv   = 3
    Hp_pv = 4
  !---------------------------------------------------------------------------------------------------- 
    CALL CPU_TIME(T1_CPU)
  !----------------------------------------------------------------------------------------------------
    CALL LECTURE_DONNEES(Lx)

  !----------------------------------------------------------------------------------------------------
    DX = Lx/DFLOAT(Nx)
    X0 = 0.D0
    isave = -1
    print*, 'Nx =', Nx,'DX =', DX
   !----------------------------------------------------------------------------------------------------
   ALLOCATE( Ein(0:Nx+1), Sound_ax(0:Nx+1))
   Allocate(Pression(0:Nx+1) , Qaverage(1:Nx))  
   ALLOCATE( Prim(Nv_Prim,0:Nx+1), CONS(Nv_Prim,1:Nx), FLUX(Nv_Prim,0:Nx), pente(Nv_Prim,0:Nx+1))  
   ALLOCATE( NamesOfPrimVar(Nv_Prim), U1(1:Nx))
   ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim), X(1:Nx))
   flux(:,:)=0.d0; cons(:,:) = 0.d0; prim(:,:) = 0.d0; U1(:)=0.d0
   Pression(:) = 0.d0; Qaverage(:) = 0.d0
  !----------------------------------------------------------------------------------------------------
    NamesOfPrimVar(H_pv)         = "Depht Lenght"
    NamesOfPrimVar(U_pv)         = "Velocity (x)"
    NamesOfPrimVar(Phi_pv)       = " enstrophy"
    NamesOfPrimVar(Hp_pv)       = " true hight"

  !----------------------------------------------------------------------------------------------------
    DO iv = 1, Nv_Prim
       WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx
    END DO
    WRITE(6,*) " >>> End of LECTURE_DONNEES"
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
   CALL INITIALISATION(DX,X, Prim, SOUND_ax,Ein,CONS, Pression,  Lx)
 
  !----------------------------------------------------------------------------------------------------
   TIME = 0.D0
   period_time = 1.d0;
   TIME2  = period_time ;
   IT = 1
  CALL PutonScreen()
  !----------------------------------------------------------------------------------------------------
   XMIN  = MINVAL(X(:))
   XMAX1 = MaxVal(X(:))

    
  WRITE(6,'( 2(A10,E15.6))')  " Xmin = ", XMIN, &
       &                     " Xmax = ", XMAX1
   !----------------------------------------------------------------------------------------------------

    CALL Ecriture_donnees(X,Prim, Ein, Pression, time, DX, U1, Lx, Qaverage)

       !OPEN(argunit+1,FILE = './resu/Umax.out')
       OPEN(argunit+2,FILE = './resu/QaverOverT.out')

   !-------------------------do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT)!-------------------
   ! BOUCLE SUR LE TEMPS

    Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))
    !----------------------------------------------------------------------------------------------------
    Umax = MAXVAL(dabs(Prim(U_pv,1:Nx)))
    Cxmax = MAXVAL(dabs(Sound_ax(1:Nx)) ) 
   !---------------------------------------------------------------------------------------------------- 
    DT   = dx/(umax + cxmax)
    DT   = CFL*DT; dt2 = 0.5d0*dt
   !----------------------------------------------------------------------------------------------------
    IF( It == 1 ) WRITE(6,*) " Dt = ", Dt

    TIME = TIME + DT
!----------------------------------------------------------------------------------------------------
if (method_source_term == 1) then 
CALL euler_method(DT2, time, x, CONS, Prim, Ein, SOUND_ax,Pression, it, U1, Qaverage)
else if (method_source_term == 2) then 
CALL rk2(DT2, CONS, Prim, Ein, SOUND_ax,Pression, it, U1)
endif

CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax,CONS, Pression, it)

IF (cond_lim == 1) THEN 
CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax, TIME)
ELSE IF (cond_lim == 6) THEN 
CALL Powder(Prim, Ein, Pression, SOUND_ax)
ENDIF

CALL PENTE1_x( Prim, PENTE)
!PENTE = 0.D0
!----------------------------------------------------------------------------------------------------
call HLLC_x_sub1(prim,flux, pente, cons, it, dt) 

!print*, '2 dt = ', dt
!----------------------------------------------------------------------------------------------------
call godunov_x_sub1(cons,flux,dt,dx)

 !----------------------------------------------------------------------------------------------------
CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax,CONS, Pression, it)
!   !----------------------------------------------------------------------------------------------------

!   !----------------------------------------------------------------------------------------------------

IF (cond_lim == 1) THEN 
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
 CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
 CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,TIME)
 ELSE IF (cond_lim == 6) THEN 
CALL Powder(Prim, Ein, Pression, SOUND_ax)
ENDIF

!----------------------------------------------------------------------------------------------------

if (method_source_term == 1) then 
CALL euler_method(DT2, time, x, CONS, Prim, Ein, SOUND_ax,Pression, it, U1, Qaverage)
else if (method_source_term == 2) then 
CALL rk2(DT2, CONS, Prim, Ein, SOUND_ax,Pression, it, U1)
endif


 !----------------------------------------------------------------------------------------------------
CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax,CONS, Pression, it)
!   !----------------------------------------------------------------------------------------------------

!   !----------------------------------------------------------------------------------------------------

IF (cond_lim == 1) THEN 
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
 CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
 CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,TIME)
 ELSE IF (cond_lim == 6) THEN 
CALL Powder(Prim, Ein, Pression, SOUND_ax)
ENDIF


  !----------------------------------------------------------------------------------------------------
    IT = IT + 1
  !----------------------------------------------------------------------------------------------------
    IF (TIME2.LE.TIME) THEN
      CALL Ecriture_donnees(X, Prim, Ein, Pression, time, DX, U1, LX, Qaverage)
    WRITE(argunit+1,'(2(E20.13,1X))') time,  MAXVAL(Prim(U_pv,1:Nx))
    
    END IF
  !----------------------------------------------------------------------------------------------------
    IF (TIME2.LE.TIME) THEN
      PRINT*, 'EN', IT, 'ITERATIONS, ', ' TIME:', TIME
      TIME2 = TIME + period_time
    END IF
  !----------------------------------------------------------------------------------------------------
    CALL PutonScreen()

    
  !----------------------------------------------------------------------------------------------------
    ENDDO TIME_LOOP

    do ix = 1, Nx
     WRITE(argunit+2,'(2(E20.13,1X))') X(ix),  Qaverage(ix)
    enddo
  !----------------------------------------------------------------------------------------------------
   ! FIN BOUCLE SUR LE TEMPS

   CALL Ecriture_donnees(X,Prim, Ein, Pression, time, DX, U1, Lx, Qaverage)
  !----------------------------------------------------------------------------------------------------
   DEALLOCATE(Prim, Ein, SOUND_ax,X, Qaverage, &
    & Pression, CONS, FLUX, MinVP, MaxVp, NamesOfPrimVar, U1)  
   close(argunit+1)
   close(argunit+2)
  !----------------------------------------------------------------------------------------------------
   CALL CPU_TIME(T2_CPU)
  !----------------------------------------------------------------------------------------------------
    PRINT*, 'L EXECUTION DU PROGRAMME A PRIS', T2_CPU - T1_CPU
    PRINT*, 'EN', IT-1, 'ITERATIONS, ', ' TIME:', TIME
  !----------------------------------------------------------------------------------------------------
    STOP
  !----------------------------------------------------------------------------------------------------
  CONTAINS
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !!!suivant x

  Subroutine HLLC_x_sub1(prim,flux, pente, cons, it, dt)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, it
  real(kind=dp) :: dt, dt2
  REAL (KIND = DP) :: cons(Nv_Prim,1:Nx)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:)::consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED
  REAL (KIND = DP) :: FLUX(Nv_Prim,0:Nx)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), pente(Nv_Prim,0:Nx+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:) :: EinG, EinD, PressionG, PressionD
  real(kind=dp)    :: ul,ur,hl,hr,pr,pl,phil,phir, cl,cr,ml,mr,sl,sr, hpl, hpr
  real(kind=dp)    :: EL, ER, mpl, mpr
  real(kind=dp)    :: pstar,ustar,Estar,hstar,hpstar

  ALLOCATE(consG(Nv_Prim,1:Nx),consD(Nv_Prim,1:Nx), FLUXG(Nv_Prim,0:Nx))
  Allocate(FLUXD(Nv_Prim,0:Nx),PrimG(Nv_Prim,0:Nx+1))
  Allocate(PrimD(Nv_Prim,0:Nx+1), CONS_PRED(Nv_Prim,1:Nx))
  ALLOCATE(EinG(0:Nx+1), EinD(0:Nx+1), PressionG(0:Nx+1), PressionD(0:Nx+1))
  dt2 = 0.5d0*dt; PrimG(:,:) = 0.D0; PrimD(:,:)= 0.D0; 
  FLUXG(:,:) = 0.D0; FLUXD(:,:)= 0.D0; CONSG(:,:) = 0.D0; CONSD(:,:)= 0.D0
  EinG= 0.D0; EinD= 0.D0; PressionG= 0.D0; PressionD= 0.D0; CONS_PRED = 0.D0

do ix = 1, Nx
do iv = 1, Nv_Prim
PrimG(iv, ix) = Prim(iv, ix) - 0.5d0*Pente(iv, ix)
PrimD(iv, ix) = Prim(iv, ix) + 0.5d0*Pente(iv, ix)
enddo
EinG(ix) = 0.5d0*(g*PrimG(H_pv, ix)*dcos(angle) &
  &+ (PrimG(phi_pv, ix)+phi2)*PrimG(H_pv, ix)**2.d0)
EinD(ix) = 0.5d0*(g*PrimD(H_pv, ix)*dcos(angle) + &
  & (PrimD(phi_pv, ix) +phi2)*PrimD(H_pv, ix)**2.d0)

PressionG(ix) = 0.5d0*g*dcos(angle)*PrimG(H_pv, ix)**2.d0 &
& + PrimG(H_pv, ix)**3.d0*(PrimG(phi_pv, ix)+phi2)
PressionD(ix) = 0.5d0*g*dcos(angle)*PrimD(H_pv, ix)**2.d0 &
&+ PrimD(H_pv, ix)**3.d0*(PrimD(phi_pv, ix)+phi2)    
enddo

IF (cond_lim == 1) THEN 
CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax)
CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax)
CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax, TIME)
CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax, TIME)
ELSE IF (cond_lim == 6) THEN 
CALL Powder(PrimG, EinG, PressionG, SOUND_ax)
CALL Powder(PrimD, EinD, PressionD, SOUND_ax)
ENDIF

CALL PRIM_TO_CONS_FLUX_sub1x(PrimG, EinG,PressionG, CONSG, FLUXG)
CALL PRIM_TO_CONS_FLUX_sub1x(PrimD, EinD,PressionD, CONSD, FLUXD)

CALL CONS_PREDICTION_sub1x(DX, DT, CONSG, CONS_PRED, FLUXG, FLUXD)

if (method_source_term == 1) then 
CALL euler_method(DT2, time, x,CONS_PRED, PrimG, EinG, SOUND_ax,PressionG, it, U1, Qaverage)
else if (method_source_term == 2) then 
CALL rk2(DT2, CONS_PRED, PrimG, EinG, SOUND_ax, PressionG, it, U1)
endif

CALL NOUVELL_VARIABLE_PRIM(PrimG, EinG, SOUND_ax,CONS_PRED, PressionG, it)

IF (cond_lim == 1) THEN 
CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax, TIME)
ELSE IF (cond_lim == 6) THEN 
CALL Powder(PrimG, EinG, PressionG, SOUND_ax)
ENDIF

CALL PRIM_TO_CONS_FLUX_sub1x(PrimG, EinG,PressionG, CONSG, FLUXG)

CALL CONS_PREDICTION_sub1x(DX, DT, CONSD, CONS_PRED, FLUXG, FLUXD)

if (method_source_term == 1) then 
CALL euler_method(DT2, time, x, CONS_PRED, PrimD, EinD, SOUND_ax,PressionD, it, U1, Qaverage)
else if (method_source_term == 2) then 
CALL rk2(DT2, CONS_PRED, PrimD, EinD, SOUND_ax, PressionD, it, U1)
endif

CALL NOUVELL_VARIABLE_PRIM(PrimD, EinD, SOUND_ax, CONS_PRED, PressionD, it)

IF (cond_lim == 1) THEN 
CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax)
ELSE IF (cond_lim == 2) THEN 
CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax)
ELSE IF (cond_lim == 3) THEN 
CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax, TIME)
ELSE IF (cond_lim == 6) THEN 
CALL Powder(PrimD, EinD, PressionD, SOUND_ax)
ENDIF

do ix=0,nx

!! etat gauche et droite
ul=primD(U_pv,ix);  ur=primG(U_pv,ix+1)
hl=primD(h_pv,ix);  hr=primG(h_pv,ix+1)
hpl=primD(hp_pv,ix);  hpr=primG(hp_pv,ix+1)

phil=primD(phi_pv,ix);  phir=primG(phi_pv,ix+1)

!! calcul des énergie
El=(ul*ul+g*hl*dcos(angle)+(phil+phi2)*hl**2.d0 )*0.5D0
Er=(ur*ur+g*hr*dcos(angle)+(phir +phi2)*hr**2.d0)*0.5D0
!! calcul des pression
pl=g*dcos(angle)*hl*hl*0.5d0+hl**3.d0*(phil+phi2) 
pr=g*dcos(angle)*hr*hr*0.5d0+hr**3.d0*(phir + phi2)
!! calcul des vitesses du son
cl=dsqrt(g*hl*dcos(angle)+3.d0*(phil+phi2)*hl**2.d0);
cr=dsqrt(g*hr*dcos(angle)+3.d0*(phir+phi2)*hr**2.d0)
! davis
sr=dmax1(ul+cl,ur+cr)
Sl=DMIN1(ul-cl,ur-cr)
! etat star 
ml=hl*(ul-sl)
mr=hr*(ur-sr)


mpl=hpl*(ul-sl)
mpr=hpr*(ur-sr)


ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)



if (ustar.ge.0.d0) then
if (sl.ge.0.d0) THEN
!!etat gauche
flux(1,ix)=hl*ul
flux(4,ix)=hpl*ul
flux(2,ix)=hl*ul*ul+pl
flux(3,ix)=hl*El*ul+pl*ul
ELSE
!! etat star gauche
hstar=ml/(ustar-sl)

hpstar = mpl/(ustar - sl)

Estar=El+(pl*ul-pstar*ustar)/ml
!! remplissage des flux

flux(1,ix)=hstar*ustar!hl*ul + sl*(hstar - hl)
flux(4,ix)=hpstar*ustar
flux(2,ix)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
flux(3,ix)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 

endif
ELSE
if (sr.ge.0.d0) then
!!etat droit etoile
!!etat star droit
hstar=mr/(ustar-sr)

hpstar=mpr/(ustar-sr)

Estar=Er+(pr*ur-pstar*ustar)/mr
!remplissage des flux

flux(1,ix)=hstar*ustar!hl*ul + sl*(hstar - hl)
flux(4,ix)=hpstar*ustar
flux(2,ix)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
flux(3,ix)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 

ELSE
!!etat droit
flux(1,ix)=hr*ur
flux(4,ix)=hpr*ur
flux(2,ix)=hr*ur*ur+pr
flux(3,ix)=hr*Er*ur+pr*ur
end if
end if

end do

DEALLOCATE(consG,consD, FLUXG,FLUXD, PrimG, PrimD, &
  & CONS_PRED, EinG, EinD, PressionG, PressionD )
return
end subroutine HLLC_x_sub1

  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_x_sub1(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(Nv_Prim,1:nx), FLUX(Nv_Prim,0:Nx)
  real(Kind=dp) :: dt, dx
  INTEGER :: k, ix, iy

  do ix=1,nx
!------------------------------------------------------------------------------------------------------------------------
    do k=1,Nv_Prim
    cons(k,ix)=cons(k,ix)+dt/dx*(Flux(k,ix-1)-flux(k,ix))
   end do
  end do

  end subroutine godunov_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE PutonScreen()
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  REAL (KIND = DP)     :: MinVp(Nv_Prim), MaxVp(Nv_Prim)

IF ( MOD(IT,ImpreE) == 0 ) THEN
do iv = 1, Nv_Prim
MinVp(iv)     = MINVAL(Prim(iv,1:Nx))
MaxVp(iv)     = MAXVAL(Prim(iv,1:Nx))
enddo

WRITE(argunit,*)
WRITE(argunit,'(65(":"))') 
WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
WRITE(argunit,'(":: | Time  = ", E10.3, 1x, &
     &   "  ||     kt =  ",I9,3x,"| ",10x,"::")' ) TIME, it
WRITE(argunit,'(":: |    Dt = ", E10.3,1x, &
     &   "  ||    Cfl =  ",E10.3,"  | ",10x,"::")' ) Dt , CFL
WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
WRITE(argunit,'("::",61(" "),"::")') 
WRITE(argunit,'("::   ",12x, "    ",  2(3x,A12), 11x," ::" )') " Minimum ", " Maximum "

DO iv = 1, Nv_Prim
   WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )') &
        & TRIM(NamesOfPrimVar(iv)), MinVp(iv) , MaxVp(iv)
END DO

WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )')  
WRITE(argunit,'(65(":"))') 
END IF

END SUBROUTINE PutonScreen
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE Ecriture_donnees(X,Prim, Ein, Pression, time, DX, U1, Lx , Qaverage)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy , MyUnit = 30, il
  REAL (KIND = DP) :: X(1:Nx), U1(1:Nx), Vcore, Qaverage(1:Nx)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), Ein(0:Nx+1), Pression(0:Nx+1)
  REAL (KIND = DP) :: time, DX, dt, Lx, rhoP, VairVel, VairTurb
  CHARACTER(LEN=3) :: NB, Zero="000"

  ! ECRITURE DES RESULTATS

  OPEN(MyUnit+1,FILE = './resu/powder628.txt')
  OPEN(MyUnit+3,FILE = './resu/qdt.out')
 
    
 ! OPEN(MyUnit+2,FILE = './resu/FuncTime.out')

  Do ix = 1, Nx

if (cond_lim == 6) then 

      if ( x(ix) >= Lx/6.d0 ) then
    U1(ix) = (0.4d0*dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e


    !if ( x(ix) >= Lx/6.d0 ) then
    !U1(ix) = (0.5d0*dsin(0.4d0*pi*time)**2.d0)*&
    !&dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    !& (35.D0)*time)**2.D0/(2.d0*(17.5d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e


    

    ! if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    !U1(ix) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    !&dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    !& (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e
   
   else

     U1(ix) = 0.D0
   endif 


else 

   U1(ix) = 0.D0
endif 

 if (U1(ix) .le. eps) then
   U1(ix) = 0.D0
  endif 


    Vcore = injection_speed*injection_density/rhoPowder 
    !Qaverage1(ix) = Qaverage1(ix)+Vcore*U1(ix)*DT/dx 

rhoP = rhoAir + Prim(h_pv,ix)*RhoPowder/Prim(Hp_pv, ix)
  if (  rhoP .gt.2.d0*rhoAir .and. Prim(U_pv,ix) .gt. 1.d-8 ) then 
     VairVel = (0.1d0*Prim(U_pv,ix)**2.d0)*(rhoP-rhoAir)*rhoAir/rhoPowder ! *(phi2 +Phi)*H**2.d0
     VairTurb = 0.2d0*(phi2 +prim(phi_pv,ix))*Prim(H_pv,ix)**2.d0*(rhoP-rhoAir)*rhoAir/rhoPowder 
    else 
      !Vair = 0.25d0*dabs(u)*(rhoP-rhoAir)*rhoAir/rhoPowder
     VairVel =  0.d0
     VairTurb = 0.d0
    endif 




 
        WRITE(MyUnit+1,'(9(E20.13,1X))') X(ix), Prim(H_pv,ix), Prim(U_pv,ix), &
        & prim(phi_pv,ix), U1(ix), Prim(Hp_pv, ix), rhoP, &
        & VairVel , &
        & VairTurb

        WRITE(MyUnit+3,'(2(E20.13,1X))') X(ix),Qaverage(ix)
     
        ! if ( dabs(x(ix)- 200.d0) .le. dx) then
        !   WRITE(MyUnit+2,'(5(E20.13,1X))') time, Prim(H_pv,ix), Prim(U_pv,ix), &
        ! & prim(phi_pv,ix), Prim(Hp_pv, ix)
        ! endif 
  
END DO
    WRITE(MyUnit+1,*)  
    WRITE(MyUnit+3,*)  

!               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
 !close(MyUnit+1) 
return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE LECTURE_DONNEES(Lx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   REAL (KIND = DP) :: Lx
   OPEN(UNIT=21, FILE = 'data1D.inp', STATUS = 'OLD')  
    READ(21,*) cond_lim, method_source_term ! 1 box/ 2 absorbtion/ 3 batteur/ 4 jump; cond_lim ;   test_shear
    READ(21,*) angle                 ! inclination angle
    READ(21,*) Nx                ! NUMBER OF CELLS
    READ(21,*) Lx                ! DOMAIN LENGTH 
    READ(21,*) TIMEOUT                ! OUTPUT TIME
    READ(21,*) iterfinal              ! Iteration final
    READ(21,*) g, CFL                 ! acceleration due to gravity
    READ(21,*) H_0                    ! stationary unstable solution (H_0, U_0, Phi_0 = 0)
    READ(21,*) phi2                   ! patit enstrophy
    READ(21,*) frottcoeff, disscoeff  ! cf, cr
    READ(21,*) amplitude              ! amplitude des perturbation
    READ(21,*) ImpreE, ImpreF
    close(21)
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE INITIALISATION(DX, X, Prim, SOUND_ax,Ein,CONS, Pression,  Lx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix
  REAL (KIND = DP) :: DX, Lx, U_0
  REAL (KIND = DP) :: X(1:Nx)
  REAL (KIND = DP) :: Prim(Nv_Prim, 0:Nx+1)
  REAL (KIND = DP) :: Pression(0:Nx+1), CONS(Nv_Prim,1:Nx), Ein(0:Nx+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1)
  REAL (KIND = DP) :: SOUND_bx(0:Nx+1), SOUND_ay(0:Nx+1), SOUND_by(0:Nx+1)

  pi = 4.0d0*ATAN(1.0d0)
  U_0 = DSQRT(g*Dtan(angle)*H_0/frottcoeff) 
  
DO ix = 1, Nx
X(ix) = 0.5D0*DX + (ix-1)*DX  

IF (cond_lim == 1 ) THEN ! ! cond_lim: 1 box/ 2 absorption / 3 batteur/ 4 hydraulic jump
 Prim(H_pv, ix) = H_0*(1.d0 + amplitude*dsin(2.d0*Pi*x(ix)/Lx) )
 Prim(U_pv,ix) = U_0
 Prim(Phi_pv,ix) = 0.d0
ELSE IF (cond_lim == 2 ) THEN 
  Prim(H_pv, ix) = H_0
  Prim(U_pv,ix) = U_0
  Prim(Phi_pv,ix) =  0.d0
ELSE IF (cond_lim == 3 ) THEN 
  Prim(H_pv, ix)= H_0
  Prim(U_pv,ix) = U_0
  Prim(Phi_pv,ix) = 0.d0
else if (cond_lim == 6) then 
!  if ( x(ix) .le. 50.d0 ) then 
    Prim(H_pv, ix) = H_0
    Prim(U_pv,ix) = 0.d0
    Prim(Phi_pv,ix) = 0.d0
 ! else 
  !  Prim(H_pv, ix) = H_0
   ! Prim(U_pv,ix) = 0.d0
    !Prim(Phi_pv,ix) = 0.d0 
  !endif
ENDIF

Prim(Hp_pv, ix) =   Prim(H_pv, ix) 
SOUND_ax(ix) = DSQRT( g*dcos(angle)*Prim(H_pv, ix) + &
&  3.d0*(Prim(Phi_pv, ix) +phi2)*Prim(H_pv, ix)**2.d0)
Pression(ix) = g*dcos(angle)*Prim(H_pv, ix)*Prim(H_pv, ix)/2.d0 &
& + (Prim(Phi_pv, ix)+phi2 )*Prim(H_pv, ix)**3.d0
Ein(ix) = ( Prim(H_pv, ix)*g*dcos(angle) +&
& (Prim(Phi_pv, ix)+phi2)*Prim(H_pv, ix)**2.D0)/2.d0

ENDDO


! VARIABLE CONSERVATIVES 
DO ix = 1, Nx
CONS(1,ix) = Prim(H_pv, ix)
CONS(4,ix) = Prim(Hp_pv, ix)
CONS(2,ix) = Prim(H_pv, ix)*Prim(U_pv, ix)
CONS(3,ix) = Prim(H_pv, ix)*(Ein(ix)+ Prim(U_pv, ix)**2.d0/2.d0)
END DO
return
END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax,CONS, Pression, it)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix,it
REAL (KIND = DP) :: phi
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), CONS(Nv_Prim,1:Nx)
REAL (KIND = DP) ::  Ein(0:Nx+1), Pression(0:Nx+1)
REAL (KIND = DP) :: SOUND_ax(0:Nx+1)

! CALCUL NOUVELS VARIABLES PRIMITIVES
DO ix = 1, Nx

Prim(H_pv,ix) = CONS(1,ix)
Prim(Hp_pv,ix)=CONS(4,ix)
Prim(U_pv,ix) = CONS(2,ix)/CONS(1,ix)
!------------------------------------------------------------------------------------------------------------------------

if (dabs(Prim(U_pv,ix)).le.1.d-8) then 
Prim(U_pv,ix)=0.d0; cons(2, ix) = 0.d0
endif
!------------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------------
Ein(ix) = CONS(3,ix)/CONS(1,ix)- 0.5d0*Prim(U_pv,ix)**2.d0 

Prim(Phi_pv,ix) = (2.d0*Ein(ix)  - &
  & g*dcos(angle)*Prim(H_pv,ix))/Prim(h_pv,ix)**2.d0 - phi2  
!if (Prim(P11_pv,ix, iy).le.-1.d-12) print*, 'p11', Prim(P11_pv,ix, iy)
Prim(Phi_pv,ix) = dmax1(Prim(Phi_pv,ix), 1.d-9)

SOUND_ax(ix) = DSQRT( g*Prim(H_pv,ix)*dcos(angle) +&
& 3.d0*(Prim(Phi_pv,ix)  +phi2)*Prim(H_pv,ix)**2.D0)

Pression(ix) = g*dcos(angle)*(Prim(H_pv,ix))**2.d0/2.d0&
& + (Prim(phi_pv,ix)+phi2)*Prim(H_pv,ix)**3.d0
  
END DO
return
END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE Condition_lim_box(Prim, Ein, Pression, SOUND_ax)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iv
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), Ein(0:Nx+1),Pression(0:Nx+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1)


Prim(:,0)    = Prim(:,Nx)

Ein(0)       = Ein(Nx)
SOUND_ax(0)  = SOUND_ax(Nx)
Pression(0)  = Pression(Nx)

!------------------------------------------------------------------------------------------------------------------------

Prim(:,Nx+1)   = Prim(:,1)

Ein(Nx+1)      = Ein(1)
SOUND_ax(Nx+1)  = SOUND_ax(1)

Pression(Nx+1)  = Pression(1)
   
return
END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  
SUBROUTINE CondLimABSORPTION(Prim, Ein, Pression,SOUND_ax)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER ::  ix
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), Ein(0:Nx+1), Pression(0:Nx+1)
REAL (KIND = DP) ::  SOUND_ax(0:Nx+1) 

!------------------------------------------------------------------------------------------------------------------------

Prim(:,0)  = Prim(:,1)

Ein(0)  = Ein(1)
SOUND_ax(0) = SOUND_ax(1)

Pression(0)  = Pression(1)

Prim(:,Nx+1) = Prim(:,Nx)

Ein(Nx+1)  = Ein(Nx)
SOUND_ax(Nx+1) = SOUND_ax(Nx)
Pression(Nx+1)  = Pression(Nx)

return
END SUBROUTINE

 !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  
SUBROUTINE Powder(Prim, Ein, Pression, SOUND_ax)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER ::  ix
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), Ein(0:Nx+1), Pression(0:Nx+1)
REAL (KIND = DP) ::  SOUND_ax(0:Nx+1) 
!------------------------------------------------------------------------------------------------------------------------

Prim(:,0)  = Prim(:,1)

! Prim(H_pv,0)  = H_0 !Prim(H_pv,1,  iy)
! Prim(Hp_pv,0)  = Prim(H_pv,0)
! Prim(U_pv,0)  = 0.d0
! Prim(Phi_pv,0)  = 0.d0 !Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2 + eps 



SOUND_ax(0) = DSQRT( g*dcos(angle)*Prim(H_pv, 0) + &
  &3.d0*(Prim(Phi_pv,0) +phi2)* Prim(H_pv, 0)**2.d0)

Pression(0) = g*dcos(angle)*Prim(H_pv, 0)**2.d0/2.d0 + &
& (Prim(Phi_pv, 0)  +phi2)*Prim(H_pv, 0)**3.d0
Ein(0) = ( Prim(H_pv, 0)*g*dcos(angle) +&
& (phi2+Prim(Phi_pv, 0) )*Prim(H_pv, 0)**2.d0 )/2.d0


Prim(:,Nx+1) = Prim(:,Nx)

Ein(Nx+1)  = Ein(Nx)
SOUND_ax(Nx+1) = SOUND_ax(Nx)

Pression(Nx+1)  = Pression(Nx)

return
END SUBROUTINE
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,TIME)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER ::  ix
  REAL (KIND = DP) :: q_0,u_0, xi, omega, TIME
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), Ein(0:Nx+1), Pression(0:Nx+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1)
    pi = 4.0d0*ATAN(1.0d0)
    omega = 6.73d0
!------------------------------------------------------------------------------------------------------------------------
   if (angle .le. 1.d-8 .and. frottcoeff .le. 1.d-8 .and. disscoeff .le. 1.d-8) then 
      u_0 = 3.d0
   else 
      u_0 = DSQRT(H_0*g*Dtan(angle)/frottcoeff)
   end if 

   q_0 = H_0*u_0  


Prim(H_pv,0) = H_0*(1.d0 + amplitude*DSIN(omega*TIME))


Prim(Hp_pv,0) = Prim(H_pv,0) 
Prim(U_pv,0)    = q_0/Prim(H_pv,0) 


Prim(Phi_pv,0)  = 0.d0


Ein(0)          =  ( Prim(H_pv, 0)*g*dcos(angle) +&
& (phi2+Prim(Phi_pv, 0) )*Prim(H_pv, 0)**2.d0 )/2.d0

SOUND_ax(0)     = DSQRT( g*dcos(angle)*Prim(H_pv, 0) + &
  &3.d0*(Prim(Phi_pv,0) +phi2)* Prim(H_pv, 0)**2.d0) 

Pression(0)     = g*dcos(angle)*Prim(H_pv, 0)**2.d0/2.d0 + &
& (Prim(Phi_pv, 0)  +phi2)*Prim(H_pv, 0)**3.d0

Prim(H_pv,Nx+1)   = Prim(H_pv,Nx)
Prim(Hp_pv,Nx+1)   = Prim(Hp_pv,Nx)
Prim(U_pv,Nx+1)   = Prim(U_pv,Nx)
Prim(Phi_pv,Nx+1) = Prim(Phi_pv,Nx)

Ein(Nx+1)         = Ein(Nx)
SOUND_ax(Nx+1)    = SOUND_ax(Nx)
Pression(Nx+1)    = Pression(Nx)

return
END SUBROUTINE
  !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE euler_method(DT,time, X, CONS, Prim, Ein, &
    & SOUND_ax,Pression, it, U1, Qaverage)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, IT,k, MyUnit1=10
  REAL (KIND = DP) :: DT, traceP, time, H, phi,u, hp
  REAL (KIND = DP) :: rhoP 
  REAL (KIND = DP) :: Vair, Vcore
  REAL (KIND = DP) :: CONS(Nv_Prim,1:Nx), U1(1:Nx)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1), X(1:Nx)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1), Qaverage(1:Nx)
  REAL (KIND = DP) :: Ein(0:Nx+1), Pression(0:Nx+1),fracp11,fracp22,erreur
  REAL (KIND = DP), ALLOCATABLE :: TS(:)

  Allocate(TS(Nv_Prim))
  
 TS(:) = 0.d0

!------------------------------------------------------------------------------------------------------------------------
DO  ix = 1, Nx

!------------------------------------------------------------------------------------------------------------------------ 
  if (cons(1,ix).le.1.d-8) then
  print*, 'pas de l eau', cons(1,ix), ix, it
  stop
  end if
!------------------------------------------------------------------------------------------------------------------------
!TERME SOURCE
H = Prim(H_pv,ix); hP = Prim(hp_pv, ix)
u =  Prim(U_pv,ix)


if (cond_lim == 6) then

      
    if ( x(ix) >= Lx/6.d0 ) then
    U1(ix) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0))

    ! if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    !U1(ix) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    !&dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    !& (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e
   
   else

     U1(ix) = 0.D0
   endif 


 if (U1(ix) .le. eps) then
   U1(ix) = 0.D0
  endif 



rhoP = rhoAir + H*RhoPowder/hP


   if (  rhoP .gt.2.d0*rhoAir .and. u .gt. 1.d-8) then 
     Vair = (0.1d0*U + 0.9d0*dsqrt((phi2 +Phi)*H**2.d0 ) )*H/hP ! *(phi2 +Phi)*H**2.d0
    else 
      !Vair = 0.25d0*dabs(u)*(rhoP-rhoAir)*rhoAir/rhoPowder
     Vair =  0.d0
    endif 


   ! if (dabs(U1(ix)) .le. eps) then 
    !  Vcore = 0.D0
    !else

    Vcore = injection_speed*injection_density/rhoPowder

   !end if 

else
  Vair = 0.D0
  Vcore = 0.d0
  U1(ix) = 0.d0
endif 

H = Prim(H_pv,ix)
Phi = Prim(Phi_pv,ix); 
u =  Prim(U_pv,ix); 

  TS(1) = Vcore*U1(ix)
  TS(4) = Vcore*U1(ix) + Vair
  TS(2) = g*Dsin(angle)*H - frottcoeff*dsqrt( u**2.d0)*u  + 5.d0*U1(ix)**2.d0*Vcore!  
  TS(3) =   g*Dsin(angle)*H*U - (frottcoeff+disscoeff*Phi/(phi2 +Phi))*dabs(u)*u**2.d0&
  & + 5.d0*(U1(ix)**2.d0)*Vcore*U + Vcore*U1(ix)*(g*dcos(angle)*h &
    &-U**2.d0/2.d0 +3.d0/2.d0*(phi2+Phi)*h**2.d0 )  ! 
  Qaverage(ix) = Qaverage(ix)+Vcore*U1(ix)*dt

  
!-----------------------------------------------------------------------------------------------------------------

do k=1,Nv_Prim
CONS(k,ix) = CONS(k,ix) + DT*TS(k) 
end do 
ENDDO
!------------------------------------------------------------------------------------------------------------------------
CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, CONS, Pression, it)
!------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------

DEALLOCATE(TS)
return
  END SUBROUTINE
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 SUBROUTINE rk2(DT, CONS, Prim, Ein, SOUND_ax, Pression, it, U1)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, IT,k
  REAL (KIND = DP) :: DT, dt2, H, phi,  u
  REAL (KIND = DP) :: rhoP
  REAL (KIND = DP) ::  Vair, Vcore , U1(1:Nx)
  REAL (KIND = DP) :: CONS(Nv_Prim,1:Nx)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1)
  REAL (KIND = DP) :: Ein(0:Nx+1), Pression(0:Nx+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:)  :: TS1, TS2, CONS_k1

  ALLOCATE(TS1(Nv_Prim,1:Nx), TS2(Nv_Prim,1:Nx), CONS_k1(Nv_Prim,1:Nx))
 
  cons_k1 = 0.d0; TS1 = 0.d0; TS2 = 0.d0
  dt2 = 0.5d0*dt
!------------------------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------------------------
DO  ix = 1, Nx
!------------------------------------------------------------------------------------------------------------------------ 
if (cons(1,ix).le.1.d-8) then
print*, 'pas de l eau', cons(1,ix), ix, it
stop
end if
!------------------------------------------------------------------------------------------------------------------------
!TERME SOURCE

if (cond_lim == 6) then

  Vair = dabs(u)*rhoAir/rhoPowder !0.2d0*dabs(u)*rhoAir/rhoPowder
       
   !rhoP = rhoAir + H*RhoPowder/hP
   
    !U1(ix,iy) =  0.D0

    if ( U1(ix) .le. eps) then 
      Vcore = 0.D0
    else

    Vcore = injection_speed*injection_density/rhoPowder
   end if 

else
  Vair = 0.D0
  Vcore = 0.d0
  U1(ix) = 0.d0

endif 

!TERME SOURCE
H = Prim(H_pv,ix)
Phi = Prim(Phi_pv,ix); u =  Prim(U_pv,ix);
!print*, 'cr=', disscoeff

TS1(1,ix) = Vcore 
TS1(4,ix) = Vcore + Vair  
TS1(2,ix) = g*Dsin(angle)*H - frottcoeff*dsqrt( u**2.d0)*u  + U1(ix)*Vcore ! g*Dtan(angle)*H  
TS1(3,ix) =  g*Dsin(angle)*H*U - (frottcoeff+disscoeff*Phi/(phi2 +Phi))*dabs(u)*u**2.d0&
  & + (U1(ix))*Vcore*U + Vcore*(g*dcos(angle)*h -U**2.d0/2.d0)  ! g*Dtan(angle)*H*U 
!------------------------------------------------------------------------------------------------------------------------
do k=1,Nv_Prim
CONS_k1(k,ix) = CONS(k,ix) + DT*TS1(k,ix) 
end do 

ENDDO
  
!------------------------------------------------------------------------------------------------------------------------
CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, CONS_k1, Pression, it)
!------------------------------------------------------------------------------------------------------------------------
DO  ix = 1, Nx

H = Prim(H_pv,ix)
Phi= Prim(Phi_pv,ix);
u =  Prim(U_pv,ix); 

TS2(1,ix) = Vcore 
TS2(4,ix) = Vcore + Vair  
TS2(2,ix) = g*Dsin(angle)*H - frottcoeff*dsqrt( u**2.d0)*u  + U1(ix)*Vcore ! g*Dtan(angle)*H  
TS2(3,ix) =  g*Dsin(angle)*H*U - (frottcoeff+disscoeff*Phi/(phi2 +Phi))*dabs(u)*u**2.d0&
&+ U1(ix)*Vcore*U + Vcore*(g*dcos(angle)*h -U**2.d0/2.d0)  ! g*Dtan(angle)*H*U 
!------------------------------------------------------------------------------------------------------------------------

do k=1,Nv_Prim
CONS(k,ix) = CONS(k,ix) + DT2*(TS1(k,ix) + TS2(k,ix) )
end do 


ENDDO
!------------------------------------------------------------------------------------------------------------------------
CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, CONS, Pression, it)
!------------------------------------------------------------------------------------------------------------------------
DEALLOCATE(TS1, TS2, cons_k1)
return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PENTE1_x( Prim, PENTE)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix
REAL (KIND = DP) :: pent, pente11, pente12
REAL (KIND = DP) :: PENTE(Nv_Prim,0:Nx+1), Prim(Nv_Prim,0:Nx+1)
pente = 0.d0
 DO ix= 1,Nx
    DO iv = 1, Nv_Prim
      PENTE11 = (Prim(iv, ix) - Prim(iv, ix- 1)) 
      PENTE12=   (Prim(iv, ix + 1) - Prim(iv, ix)) 
      call minmod(pente11,pente12,pent)    
      PENTE(iv, ix)=pent
  END DO
END DO

return
END SUBROUTINE

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !limiteur de pente:

subroutine vleer(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,prod   
if ((DABS(s1) .lt. 1.d-8) .and. (DABS(s1) .lt. 1.d-8)) then
  slim = 0.d0
  return
endif

prod=s1*s2
if(prod>0.d0.and.dabs(s1+s2)>1.d-6) then
  slim=2.d0*prod/(s1+s2)
else
  slim=0.d0
endif
return
end subroutine vleer

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine minmod(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2

if(s1*s2 > 0.d0) then
  slim=dabs(s1)
  if(dabs(s2)<slim) slim=dabs(s2)
  if(s1 < 0.d0)slim=-slim
else
  slim=0.d0
endif
return
end subroutine minmod
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!    module limiteur de pentes superbee

subroutine superbee(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2      

if(s1*s2.gt.1.d-6) then
  ss1=dabs(s1)
  ss2=dabs(s2)
  slim=dmax1(dmin1(2.d0*ss1,ss2),dmin1(ss1,2.d0*ss2))
  if(s1.lt.0.d0)slim=-slim
else
  slim=0.d0
endif
return
end subroutine superbee
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !    module limiteur de pentes Van Albada

subroutine valbada(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim  

if(s1.ne.0.d0.or.s2.ne.0.d0) then
  slim=s1*s2*(s1+s2)/(s1*s1+s2*s2)
else
  slim=0.d0
endif
return
end subroutine valbada
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !    module limiteur de pentes MC

subroutine vmc(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2   

if(s1*s2.gt.0.d0) then
  ss1=dabs(s1)
  ss2=dabs(s2)
  slim=dmin1(2.d0*ss1,2.d0*ss2,0.5d0*(ss1+ss2))
  if(s1.lt.0.d0) slim=-slim
else
  slim=0.d0
endif
return
end subroutine vmc

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PRIM_TO_CONS_FLUX_sub1x(Prim, Ein,Pression, CONS, FLUX)
USE precisions
  use GlobalParam
IMPLICIT NONE
REAL (KIND = DP):: Ein(0:Nx+1), SOUND_ax(0:Nx+1), Pression(0:Nx+1)
REAL (KIND = DP) :: CONS(Nv_Prim,1:Nx),  FLUX(Nv_Prim,0:Nx), Prim(Nv_Prim,0:Nx+1)
REAL (KIND = DP) :: h, u, phi, hp
INTEGER :: ix

! VARIABLE CONSERVATIVES 
DO ix = 1, Nx
hp = Prim(hp_pv, ix)
h = Prim(H_pv, ix); u = Prim(U_pv, ix); 
phi= Prim(Phi_pv, ix);  

CONS(1,ix) = h
CONS(4,ix) = hp
CONS(2,ix) = h*u
CONS(3,ix) = h*(Ein(ix)+ (u**2.d0 )/2.d0)

flux(1,ix)=h*u
flux(4,ix)=hp*u
flux(2,ix)=h*u*u+Pression(ix)
flux(3,ix)=h*u*(Ein(ix) + 0.5d0*u**2.d0)+Pression(ix)*u

END DO

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE CONS_PREDICTION_sub1x( DX, DT, CONS, CONS_PRED, FLUXGG, FLUXDD)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix
REAL (KIND = DP) :: DX, DT, hu, h,  hE
REAL (KIND = DP) :: CONS(Nv_Prim,1:Nx), FLUXGG(Nv_Prim,0:Nx)
REAL (KIND = DP) :: FLUXDD(Nv_Prim,0:Nx), CONS_PRED(Nv_Prim,1:Nx)

CONS_PRED = 0.d0

DO ix=1,Nx

CONS_PRED(:,ix) = CONS(:,ix) - 0.5D0*DT/DX*(FLUXDD(:,ix) - FLUXGG(:,ix))  

END DO


return
END SUBROUTINE


  END PROGRAM code2D


   
