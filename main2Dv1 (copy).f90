  module precisions
  implicit none
  INTEGER, parameter :: dp=kind(1.0d0)
  end module precisions
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  MODULE GlobalParam
    USE precisions
    IMPLICIT NONE
    INTEGER  ::  ImpreE, ImpreF, Nv_Prim, argunit = 6, isave
    INTEGER  ::  Nx, Ny, iterfinal, cond_lim
    INTEGER  ::  H_pv, U_pv, V_pv, P11_pv,  P12_pv, P22_pv, Hp_pv
    REAL (KIND = DP)  ::  X0, Y0, H_0, amplitude,  frottcoeff, disscoeff
    REAL (KIND = DP)  ::  CFL, TIMEOUT, period_time, pi, angle, g, phi2
    REAL (KIND = DP)  ::  lambda, gamma, beta
    REAL (KIND = DP), PARAMETER ::  EPS = 1.d-8
  END MODULE GlobalParam
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      

MODULE ModeleInterface
  USE GlobalParam
  USE precisions
  IMPLICIT NONE

CONTAINS
!--------------------------------------------------------
  
 FUNCTION InternalEn(h, p11, p22) RESULT(InternalE)
 REAL (KIND = DP)      :: h, p11, p22, InternalE
  InternalE =  0.5d0*(g*h + p11 + p22)
 END FUNCTION InternalEn
!--------------------------------------------------------

 FUNCTION Press(h, p11) RESULT(pres)
 REAL (KIND = DP)      :: h, p11, pres
  pres =  0.5d0*g*h**2.d0 + h*p11
 END FUNCTION Press
!--------------------------------------------------------

FUNCTION Sound_a_x(h, p11) RESULT(Sound_a)
 REAL (KIND = DP)      :: h, p11, Sound_a
  Sound_a =  dsqrt(g*h + 3.d0*p11)
 END FUNCTION Sound_a_x
!--------------------------------------------------------

FUNCTION Sound_a_y(h, p22) RESULT(sound_ath)
 REAL (KIND = DP)      :: h, p22, sound_ath
  sound_ath =  dsqrt(g*h + 3.d0*p22)
 END FUNCTION Sound_a_y
!--------------------------------------------------------

FUNCTION Sound_b_x( p11) RESULT(sound_b)
 REAL (KIND = DP)      :: h, p11, sound_b
  sound_b =  dsqrt(p11)
 END FUNCTION Sound_b_x
!--------------------------------------------------------

FUNCTION Sound_b_y( p22) RESULT(sound_bth)
 REAL (KIND = DP)      :: h, p22, sound_bth
  sound_bth =  dsqrt(p22)
 END FUNCTION Sound_b_y
!--------------------------------------------------------

END MODULE ModeleInterface


  PROGRAM code2D
  USE precisions
  use GlobalParam
  USE ModeleInterface
    IMPLICIT NONE  
    INTEGER :: iv, I, IT, ix, iy
    REAL (KIND = DP), ALLOCATABLE :: Ein(:,:), Sound_ax(:,:),  Pression(:,:)
    !REAL (KIND = DP), ALLOCATABLE :: RhoPowder(:,:)
    REAL (KIND = DP), ALLOCATABLE :: Sound_bx(:,:), Sound_ay(:,:),Sound_by(:,:)
    REAL (KIND = DP), ALLOCATABLE :: CONS(:,:,:), FLUX(:,:,:), Prim(:,:,:)
    REAL (KIND = DP), ALLOCATABLE :: UmaxTampon(:,:), VmaxTampon(:,:)
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:) :: MaxVP, MinVp, X, Y
    REAL (KIND = DP)  :: T1_CPU, T2_CPU, TIME,  TIME2
    REAL (KIND = DP) :: Lx, Ly,DX, Dy, Dh, DT, dt2
    REAL (KIND = DP) :: Hstar, Pstar, Ustar, Estar, Cmax
    REAL(KIND=DP):: xmax, XMIN, XMAX1, YMIN, YMAX, UMAX, VMAX, cxmax, cymax 
    
    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: NamesOfPrimVar
  !----------------------------------------------------------------------------------------------------
    pi = 4.0d0*ATAN(1.0d0)
    Nv_Prim  = 7
    H_pv     = 1
    U_pv     = 2
    V_pv     = 3
    P11_pv   = 4
    P12_pv   = 5
    P22_pv   = 6
    Hp_pv=7
  !---------------------------------------------------------------------------------------------------- 
    CALL CPU_TIME(T1_CPU)
  !----------------------------------------------------------------------------------------------------
    CALL LECTURE_DONNEES(Lx, Ly)
  !----------------------------------------------------------------------------------------------------
    DX = Lx/DFLOAT(Nx)
    Dy = Ly/DFLOAT(Ny)
    Dh =  dMIN1(Dx,Dy)
    X0 = 0.D0
    Y0 = 0.d0
    isave = -1
    print*, 'Nx =', Nx,'Ny =', Ny, 'DX =', DX, 'DY =', DY
   !----------------------------------------------------------------------------------------------------
   ALLOCATE( Ein(0:Nx+1,0:Ny+1), Sound_ax(0:Nx+1,0:Ny+1),Sound_bx(0:Nx+1,0:Ny+1), &
    & Sound_ay(0:Nx+1,0:Ny+1),Sound_by(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1) )  
   ALLOCATE( Prim(Nv_Prim,0:Nx+1,0:Ny+1), CONS(8,1:Nx,1:Ny), FLUX(8,0:Nx,0:Ny))  
   ALLOCATE( NamesOfPrimVar(Nv_Prim), UmaxTampon(1:Nx,1:Ny), VmaxTampon(1:Nx,1:Ny) )
   ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim), X(1:Nx), Y(1:Ny))
   flux(:,:,:)=0.d0; cons = 0.d0; prim = 0.d0; Sound_ay = 0.d0; Sound_by = 0.d0
   X = 0.d0; Y = 0.D0
  !----------------------------------------------------------------------------------------------------
    NamesOfPrimVar(H_pv)         = "Depht Lenght"
    NamesOfPrimVar(U_pv)         = "Velocity (x)"
    NamesOfPrimVar(V_pv)         = "Velocity (y)"
    NamesOfPrimVar(P11_pv)       = " tensor P11"
    NamesOfPrimVar(P12_pv)       = " tensor P12"
    NamesOfPrimVar(P22_pv)       = " tensor P22"
    NamesOfPrimVar(Hp_pv) = " True Height"
  !----------------------------------------------------------------------------------------------------
    DO iv = 1, Nv_Prim
       WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx, Ny
    END DO
    WRITE(6,*) " >>> End of LECTURE_DONNEES"
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
  CALL INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax,SOUND_bx,SOUND_ay, SOUND_by,Ein,CONS, Pression,  Lx)
  !----------------------------------------------------------------------------------------------------
   TIME = 0.D0
   period_time = 1.d0;
   TIME2  = period_time ;
   IT = 1
  CALL PutonScreen()
  !----------------------------------------------------------------------------------------------------
   XMIN  = MINVAL(X(:))
   XMAX1 = MaxVal(X(:))

   YMIN = MINVAL(Y(:))
   YMAX = MaxVal(Y(:))
    
  WRITE(6,'( 2(A10,E15.6))')  " Xmin = ", XMIN, &
       &                     " Xmax = ", XMAX1
  WRITE(6,'( 2(A10,E15.6))') " Ymin = ", YMIN, &
       &                     " Ymax = ", YMAX
   !----------------------------------------------------------------------------------------------------

    CALL Ecriture_donnees(X,Y, Lx, Prim, Ein, Pression, time, DX, DY )
 !stop
   !-------------------------do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT)!-------------------
   ! BOUCLE SUR LE TEMPS

   Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))
   !----------------------------------------------------------------------------------------------------
!print*, " i am here 2"
  !   do ix = 1, Nx
  !    do iy = 1, Ny
  !     UmaxTampon(ix, iy) = dmax1( Dabs(Prim(U_pv,ix,iy)+ &
  !       & Sound_ax(ix,iy)), Dabs(Prim(U_pv,ix,iy)- Sound_ax(ix,iy)) )
  !     VmaxTampon(ix, iy) =  dmax1( Dabs(Prim(v_pv,ix,iy)+ &
  !       & Sound_ay(ix,iy)), Dabs(Prim(v_pv,ix,iy)- Sound_ay(ix,iy)) )
  !    enddo
  !   enddo
  !   Umax = MAXVAL( UmaxTampon(1:Nx,1:Ny))
  !   vmax = MAXVAL(VmaxTampon(1:Nx,1:Ny))
  !  !----------------------------------------------------------------------------------------------------    
   
  !  if (ny > 1) then 
  !   DT   = dmin1(Dx/DABS(Umax), Dy/DABS(Vmax) )
  ! else
  !   DT   = Dx/DABS(Umax)
  ! endif
  !   DT   = CFL*DT
  !   dt2 = 0.5d0*dt



     Umax = MAXVAL(dabs(Prim(U_pv,1:Nx,1:Ny)))
    vmax = MAXVAL(dabs(Prim(v_pv,1:Nx,1:Ny)))
    Cxmax = MAXVAL(dabs(Sound_ax(1:Nx,1:Ny)) ) 
    Cymax = MAXVAL( dabs(Sound_ay(1:Nx,1:Ny) ) )
   !---------------------------------------------------------------------------------------------------- 
   if (ny .gt. 1) then   
    DT   = dmin1(dx/(umax + cxmax), dy/(vmax + cymax)) !Dh/( Umax + Cmax)
   else
     DT   = dx/(umax + cxmax)
    endif
    DT   = CFL*DT; dt2 = 0.5d0*dt

   !----------------------------------------------------------------------------------------------------
    IF( It == 1 ) WRITE(6,*) " Dt = ", Dt
    TIME = TIME + DT
!    !----------------------------------------------------------------------------------------------------
  ! call euler_method(DT2, time, Lx, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
   
  !----------------------------------------------------------------------------------------------------

    CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax,&
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)

   
   if (cond_lim == 2) then 
    CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  else if (cond_lim == 1) then
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by) 

  endif 
   
  !----------------------------------------------------------------------------------------------------
   call HLLC_x_sub1(prim,flux) 
  !----------------------------------------------------------------------------------------------------
   call godunov_x_sub1(cons,flux,dt,dx)
  !----------------------------------------------------------------------------------------------------
   CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, &
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------

   
  if (cond_lim == 2) then 
    CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  else if (cond_lim == 1) then
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by) 

  endif 
  
  call HLLC_x_sub2(prim,flux) 
  call godunov_x_sub2(cons,flux,dt,dx)
  CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, &
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
  if (cond_lim == 2) then 
    CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  else if (cond_lim == 1) then
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by) 

  endif 

  if (Ny .gt. 1 ) then 
   !----------------------------------------------------------------------------------------------------
   call HLLC_y_sub1(prim,flux) 
  !----------------------------------------------------------------------------------------------------
   call godunov_y_sub1(cons,flux,dt,dx)
  !----------------------------------------------------------------------------------------------------
   CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, &
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
  if (cond_lim == 2) then 
    CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  else if (cond_lim == 1) then
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by) 

  endif 
   !----------------------------------------------------------------------------------------------------
   call HLLC_y_sub2(prim,flux) 
  !----------------------------------------------------------------------------------------------------
   call godunov_y_sub2(cons,flux,dt,dx)
  !----------------------------------------------------------------------------------------------------
   CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, &
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------

   
  if (cond_lim == 2) then 
    CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  else if (cond_lim == 1) then
  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by) 

  endif 
  endif 
  !----------------------------------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------

  call euler_method(DT, time, Lx,X,Y, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
   
  !----------------------------------------------------------------------------------------------------

    CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)

  !----------------------------------------------------------------------------------------------------
    IT = IT + 1
  !----------------------------------------------------------------------------------------------------
     IF (TIME2.LE.TIME) THEN
     	CALL Ecriture_donnees(X,Y, Lx, Prim, Ein, Pression, time, DX, DY )
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
  !----------------------------------------------------------------------------------------------------
   ! FIN BOUCLE SUR LE TEMPS

   CALL Ecriture_donnees(X,Y, Lx,  Prim, Ein, Pression, time, DX, DY )
  !----------------------------------------------------------------------------------------------------
   DEALLOCATE(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, X, &
    & Y, Pression, CONS, FLUX, MinVP, MaxVp, NamesOfPrimVar)  
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

  Subroutine HLLC_x_sub1(prim,flux)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix,iy
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), FLUX(8,0:Nx,0:Ny)
  real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
  real(kind=dp) :: upstar, mpl, mpr, hpstar, ppstar
  real(kind=dp) :: EL, ER,p22l,p22r,vr,vl, HPl, HPr
  real(kind=dp) :: pstar,ustar,Estar,hstar,p12star

  do ix=0,nx
   do iy=1,ny
    !! etat gauche et droite
    ul=prim(U_pv,ix,iy);  ur=prim(U_pv,ix+1,iy)
    hl=prim(h_pv,ix,iy);  hr=prim(h_pv,ix+1,iy)
    vl=prim(v_pv,ix,iy);  vr=prim(v_pv,ix+1,iy)

    p11l=prim(p11_pv,ix,iy);  p11r=prim(p11_pv,ix+1,iy)
    p12l=prim(p12_pv,ix,iy);  p12r=prim(p12_pv,ix+1,iy)
    p22l=prim(p22_pv,ix,iy);  p22r=prim(p22_pv,ix+1,iy)

    HPl=prim(Hp_pv,ix,iy)
    HPr=prim(Hp_pv,ix+1,iy)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l);cr=dsqrt(g*hr+3.d0*p11r)

    
    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl

    ml=hl*(ul-sl)
    mpl=hpl*(ul-sl)

    mr=hr*(ur-sr)
    mpr=hpr*(ur-sr)

    ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
    pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)

    upstar=(ul*mpl-ur*mpr+pl-pr)/(mpl-mpr)
    ppstar=((ul-ur)*mpr*mpl+mpr*pl-mpl*pr)/(mpr-mpl)

  if (ustar.ge.0.d0) then
   if (sl.ge.0.d0) THEN
  !!etat gauche
   flux(1,ix,iy)=hl*ul
   flux(8,ix,iy)=HPl*ul
   flux(2,ix,iy)=hl*ul*ul+pl
   flux(3,ix,iy)=hl*ul*vl
   flux(4,ix,iy)= p11l*ul  !0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12l*ul
   flux(6,ix,iy)=hl*p22l*ul
   flux(7,ix,iy)=hl*El*ul+pl*ul
   ELSE
   !! etat star gauche
   hstar=ml/(ustar-sl)
   hpstar = mpl/(upstar-sl)
   p12star=p12l*(ul-sl)/(ustar-sl)
   !p12star = p12l*hstar/hl
   Estar=El+(pl*ul-pstar*ustar)/ml
   !! remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(8,ix,iy)=hpstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vl!hl*ul*vl + sl*vl*(hstar - hl)   
   flux(4,ix,iy)=p11l*ustar  !0.d0    !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar! + sl*(p12star - p12l)
   flux(6,ix,iy)=hstar*p22l*ustar!hl*p22l*ul + sl*p22l*(hstar - hl)
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 
   
   endif
  ELSE
   if (sr.ge.0.d0) then
  !!etat droit etoile

   hstar=mr/(ustar-sr)

   hpstar=mpr/(upstar -sr)


  ! p12star = p12r*hstar/hr
    p12star=p12r*(ur-sr)/(ustar-sr)
   Estar=Er+(pr*ur-pstar*ustar)/mr
   !remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(8,ix,iy)=hpstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vr!hl*ul*vl + sl*vl*(hstar - hl)   
   flux(4,ix,iy)=p11r*ustar  !0.d0    !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar! + sl*(p12star - p12l)
   flux(6,ix,iy)=hstar*p22r*ustar!hl*p22l*ul + sl*p22l*(hstar - hl)
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 
   
   ELSE
  !!etat droit
   flux(1,ix,iy)=hr*ur
   flux(8,ix,iy)=hpr*ur
   flux(2,ix,iy)=hr*ur*ur+pr
   flux(3,ix,iy)=hr*ur*vr
   flux(4,ix,iy)= p11r*ur  !0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12r*ur
   flux(6,ix,iy)=hr*p22r*ur
   flux(7,ix,iy)=hr*Er*ur+pr*ur
   end if
  end if

  end do
  end do
  return
  end subroutine HLLC_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Subroutine HLLC_x_sub2(prim,flux)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), FLUX(8,0:Nx,0:Ny)
  real(kind=dp) :: ul, ur, hl, hr, p11l, p11r, p12r, p12l, ml, mr, sr, sl, p22l, p22r
  real(kind=dp) :: EL, ER,vr,vl,HPl, HPr
  real(kind=dp) :: vstar,Estar,hp12star
  INTEGER :: ix,iy
  do ix=0,nx
  do iy=1,ny
    ul=prim(U_pv,ix,iy);     ur=prim(U_pv,ix+1,iy)
    hl=prim(h_pv,ix,iy);     hr=prim(h_pv,ix+1,iy)
    vl=prim(v_pv,ix,iy);     vr=prim(v_pv,ix+1,iy)
    p11l=prim(p11_pv,ix,iy); p11r=prim(p11_pv,ix+1,iy)
    p12l=prim(p12_pv,ix,iy); p12r=prim(p12_pv,ix+1,iy)
    p22l=prim(p22_pv,ix,iy); p22r=prim(p22_pv,ix+1,iy)

    HPl=prim(HP_pv,ix,iy)
    HPr=prim(Hp_pv,ix+1,iy)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0

    sl= - dsqrt(p11l)
    sr=dsqrt(p11r)

    vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
    hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))


    Flux(1,ix,iy) = 0.d0
    Flux(8,ix,iy) = 0.d0
    Flux(2,ix,iy) = 0.d0
    flux(3,ix,iy) = hp12star
    Flux(4,ix,iy) = 0.d0
    Flux(5,ix,iy) = vstar !! non conservatif
    flux(6,ix,iy) = 0.d0 !2.d0*hp12star*vstar !!inutile
    flux(7,ix,iy) = hp12star*vstar

  end do
  end do
  return
  end subroutine HLLC_x_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_x_sub1(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(8,1:nx,1:ny), FLUX(8,0:Nx,0:Ny)
  real(Kind=dp) :: hu, h, hv, he, hp22, dt, dx
  INTEGER :: k, ix, iy

  do ix=1,nx
   do iy=1,ny
!------------------------------------------------------------------------------------------------------------------------
    do k=1,8
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
   end do
!------------------------------------------------------------------------------------------------------------------------
  !! corection de la variable p11
  hu=cons(2,ix,iy);hv=cons(3,ix,iy);
  h=cons(1,ix,iy);hp22=cons(6,ix,iy);hE=cons(7,ix,iy)
  cons(4,ix,iy)=2.d0*hE-g*h*h-hp22-(hu*hu+hv*hv)/h
  end do
  end do
  end subroutine godunov_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  subroutine godunov_x_sub2(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
  implicit none
  real(KIND=dp) :: cons(8,1:nx,1:ny), FLUX(8,0:Nx,0:Ny)
  real(Kind=dp) :: hu,h,hv,he,hp11,p11,dt,dx
  INTEGER :: k,ix,iy

  do ix=1,nx
   do iy=1,ny
    p11=cons(4,ix,iy)/cons(1,ix,iy)
!------------------------------------------------------------------------------------------------------------------------
    do k=1,4
      cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
   end do
!------------------------------------------------------------------------------------------------------------------------
      cons(5,ix,iy)=cons(5,ix,iy)+dt/dx*(Flux(5,ix-1,iy)-flux(5,ix,iy))*p11 !! ici dans le flux on a stocké vstar
!------------------------------------------------------------------------------------------------------------------------
   do k=6,8
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
   end do
!------------------------------------------------------------------------------------------------------------------------
    !! corection de la variable p22
    hu=cons(2,ix,iy); hv=cons(3,ix,iy);
     h=cons(1,ix,iy); hp11=cons(4,ix,iy); hE=cons(7,ix,iy)
    cons(6,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp11
   end do
  end do
  end subroutine godunov_x_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !!! suivant y

  Subroutine HLLC_y_sub1(prim,flux)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), FLUX(8,0:Nx,0:Ny)
  real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
  real(kind=dp) :: EL, ER,p22l,p22r,vr,vl,temp, HPr, HPl
  real(kind=dp) :: pstar,ustar,Estar,hstar,p12star
  real(kind=dp) :: mpl, mpr, hpstar, upstar, ppstar 
  INTEGER :: ix,iy


  do ix=1,nx
   do iy=0,ny
    !! etat gauche et droite
    ul=prim(v_pv,ix,iy);   ur=prim(v_pv,ix,iy+1)
    hl=prim(h_pv,ix,iy);   hr=prim(h_pv,ix,iy+1)
    vl=prim(U_pv,ix,iy);   vr=prim(U_pv,ix,iy+1)

    p11l=prim(p22_pv,ix,iy); p11r=prim(p22_pv,ix,iy+1)
    p12l=prim(p12_pv,ix,iy); p12r=prim(p12_pv,ix,iy+1)
    p22l=prim(p11_pv,ix,iy); p22r=prim(p11_pv,ix,iy+1)

    HPl=prim(Hp_pv,ix,iy)
    HPr=prim(Hp_pv,ix,iy+1)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l);cr=dsqrt(g*hr+3.d0*p11r)

    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl
    ! etat star 
    ml=hl*(ul-sl)
    mr   =hr*(ur-sr)
    ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
    pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)


    mpl=hpl*(ul-sl)
    mpr   =hpr*(ur-sr)
    upstar=(ul*mpl-ur*mpr+pl-pr)/(mpl-mpr)
    ppstar=((ul-ur)*mpr*mpl+mpr*pl-mpl*pr)/(mpr-mpl)

  if (ustar.ge.0.d0) then
   if (sl.ge.0.d0) THEN
   !!etat gauche
   flux(1,ix,iy)=hl*ul
   flux(8,ix,iy)=hpl*ul
   flux(2,ix,iy)=hl*ul*ul+pl
   flux(3,ix,iy)=hl*ul*vl
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12l*ul
   flux(6,ix,iy)=hl*p22l*ul
   flux(7,ix,iy)=hl*El*ul+pl*ul
   ELSE

   !! etat star gauche
   hstar=ml/(ustar-sl)

   hpstar = mpl/(upstar -sl)
   !p12star=p12l*(ul-sl)/(ustar-sl)
    p12star = p12l*hstar/hl
   Estar=El+(pl*ul-pstar*ustar)/ml
   !! remplissage des flux
   flux(1,ix,iy)=hstar*ustar
   flux(8,ix,iy)=hpstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar
   flux(3,ix,iy)=hstar*ustar*vl
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar
   flux(6,ix,iy)=hstar*p22l*ustar
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar
   endif
  ELSE
   if (sr.ge.0.d0) then
  !!etat droit etoile
  !!etat star droit
   hstar=mr/(ustar-sr)
   p12star = p12r*hstar/hr
   Estar=Er+(pr*ur-pstar*ustar)/mr

   hpstar=mpr/(upstar-sr)


   !remplissage des flux
   flux(1,ix,iy)=hstar*ustar
   flux(8,ix,iy)=hpstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar
   flux(3,ix,iy)=hstar*ustar*vr
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar
   flux(6,ix,iy)=hstar*p22r*ustar
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar
   ELSE
  !!etat droit
   flux(1,ix,iy)=hr*ur
   flux(8,ix,iy)=hpr*ur
   flux(2,ix,iy)=hr*ur*ur+pr
   flux(3,ix,iy)=hr*ur*vr
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12r*ur
   flux(6,ix,iy)=hr*p22r*ur
   flux(7,ix,iy)=hr*Er*ur+pr*ur
   end if
  end if
  !! inversion des flux :
  temp=flux(2,ix,iy)
  flux(2,ix,iy)=flux(3,ix,iy)  
  flux(3,ix,iy)=temp
  temp=flux(6,ix,iy)
  flux(6,ix,iy)=flux(4,ix,iy)  
  flux(4,ix,iy)=temp
  end do
  end do
  return
  end subroutine HLLC_y_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Subroutine HLLC_y_sub2(prim,flux)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), FLUX(8,0:Nx,0:Ny)
  real(kind=dp) :: ul,ur,hl,hr,p11l,p11r,p12r,p12l,ml,mr,sr,p22l,p22r
  real(kind=dp) :: EL, ER,vr,vl, sl, HPl, HPr
  real(kind=dp) :: vstar,Estar,hp12star
  INTEGER :: ix,iy
  do ix=1,nx
  do iy=0,ny
    ul=prim(v_pv,ix,iy);     ur=prim(v_pv,ix,iy+1)
    hl=prim(h_pv,ix,iy);     hr=prim(h_pv,ix,iy+1)
    vl=prim(u_pv,ix,iy);     vr=prim(u_pv,ix,iy+1)
    p22l=prim(p11_pv,ix,iy); p22r=prim(p11_pv,ix,iy+1)
    p12l=prim(p12_pv,ix,iy); p12r=prim(p12_pv,ix,iy+1)
    p11l=prim(p22_pv,ix,iy); p11r=prim(p22_pv,ix,iy+1)

    HPl=prim(Hp_pv,ix,iy)
    HPr=prim(Hp_pv,ix,iy+1)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0

    sr=dsqrt(p11r)
    sl = -dsqrt(p11l)

    vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
    hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))

    Flux(1,ix,iy)=0.d0
    Flux(8,ix,iy)=0.d0
    Flux(3,ix,iy)=0.d0
    flux(2,ix,iy)=hp12star
    Flux(6,ix,iy)=0.d0
    Flux(5,ix,iy)=vstar !! non conservatif
    flux(4,ix,iy)=2.d0*hp12star*vstar !!inutile
    flux(7,ix,iy)=hp12star*vstar
  end do
  end do
  return
  end subroutine HLLC_y_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_y_sub1(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(8,1:nx,1:ny),flux(8,0:nx,0:ny)
  real(Kind=dp) :: hu,h,hv,he,hp11,dt,dx
  INTEGER :: k,ix,iy

  do ix=1,nx
   do iy=1,ny

    do k=1,Nv_Prim+1
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
   end do

  !! corection de la variable p22
  hu=cons(2,ix,iy);hv=cons(3,ix,iy);h=cons(1,ix,iy);hp11=cons(4,ix,iy);hE=cons(7,ix,iy)
  cons(6,ix,iy)=2.d0*hE-g*h*h-hp11-(hu*hu+hv*hv)/h
  end do
  end do
  end subroutine godunov_y_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_y_sub2(cons,flux,dt,dx)
  implicit none
  real(KIND=dp) :: cons(8,1:nx,1:ny), flux(8,0:nx,0:ny)
  real(Kind=dp) :: hu, h, hv, he, hp22, p22, dt, dx
  INTEGER :: k, ix, iy

do ix=1,nx
do iy=1,ny
 p22=cons(6,ix,iy)/cons(1,ix,iy)
 do k=1,4
  cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
 end do
  cons(5,ix,iy)=cons(5,ix,iy)+dt/dy*(Flux(5,ix,iy-1)-flux(5,ix,iy))*p22
 do k=6,Nv_Prim+1
  cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
 end do
  !! corection de la variable p11
  hu=cons(2,ix,iy);hv=cons(3,ix,iy);
  h=cons(1,ix,iy);hp22=cons(6,ix,iy);hE=cons(7,ix,iy)
  cons(4,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp22
end do
end do

  end subroutine godunov_y_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE PutonScreen()
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  REAL (KIND = DP)     :: MinVp(Nv_Prim), MaxVp(Nv_Prim)

  IF ( MOD(IT,ImpreE) == 0 ) THEN
    do iv = 1, Nv_Prim
    MinVp(iv)     = MINVAL(Prim(iv,1:Nx,1:Ny))
    MaxVp(iv)     = MAXVAL(Prim(iv,1:Nx,1:Ny))
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

  SUBROUTINE Ecriture_donnees(X,Y, Lx, Prim, Ein, Pression, time, DX, DY )
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy , MyUnit = 30, il
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny),  U1(1:Nx, 1:Ny)
  REAL (KIND = DP) :: h,u,v,p11,p12,p22, Lx, rhoAir, rhoP, rhoPowder
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  !REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:):: err_h, err_u,err_v,err_p11,err_p12,err_p22
  
  REAL (KIND = DP) :: error_h,error_u,error_v,error_p11,error_p12,error_p22
  REAL (KIND = DP) :: time, DX, DY
  CHARACTER(LEN=3) :: NB, Zero="000"
  

  rhoAir = 1.3d0; RhoPowder = 10.D0
  ! ECRITURE DES RESULTATS
if (Ny == 1 )  then 
        OPEN(MyUnit+1,FILE = './resu/Powder1D.out')
      endif
   ! Do ix = 1, Nx
   !   DO iy = 1, Ny
   !      Haverage(ix) = ( Haverage(ix) +1.d0/Ny*Prim(H_pv,ix, iy) )
   ! enddo; enddo 

  Do ix = 1, Nx
  DO iy = 1, Ny

    
      if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    U1(ix, iy) = 0.4d0*(dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e


    ! if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    !U1(ix) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    !&dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    !& (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e
   
   else

     U1(ix, iy) = 0.D0
   endif 

   rhoP = rhoAir + Prim(h_pv,ix, iy)*RhoPowder/Prim(Hp_pv, ix, iy)

    if (ny == 1) then   
      WRITE(MyUnit+1,'(11(E20.13,1X))') X(ix), Y(iy), Prim(h_pv,ix, iy), Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),&
      & prim(p11_pv,ix,iy), Prim(p12_pv,ix, iy), Prim(p22_pv,ix, iy), &
      &U1(ix, iy), Prim(Hp_pv, ix, iy), rhoP
    endif 
  END DO
  END DO
       WRITE(MyUnit+1,*)
       print*, 'Fr =', MAXVAL(dabs(Prim(U_pv,1:Nx, 1:Ny))/Sound_ax(1:Nx, 1:Ny))
    

 if (ny > 1) then

        isave = isave + 1
       WRITE(unit=NB, fmt="(I3)") isave
       NB    = ADJUSTL(NB)
       il    = LEN_TRIM(NB) 


      WRITE(6,*) " FILE = ", "./resu/Powder2D_"//Zero(1:3-il)//TRIM(NB)//".vtk"
      OPEN(UNIT=MyUnit, FILE="./resu/Powder2D_"//Zero(1:3-il)//TRIM(NB)//".vtk" )
      WRITE(MyUnit,'(''# vtk DataFile Version 2.0'')')
      WRITE(MyUnit,'(''Rectilinear 3D Dataset'')')
      WRITE(MyUnit,'(''ASCII'')')
      WRITE(MyUnit,'(''           '')')
      WRITE(MyUnit,'(''DATASET STRUCTURED_POINTS'')')
      WRITE(MyUnit,FMT='(''DIMENSIONS'',I8,I8,I8)') Nx+1, Ny+1, 2
      WRITE(MyUnit,FMT='(''ORIGIN'',3(E11.4,1x))') dx, 0.d0, 0.d0
      WRITE(MyUnit,FMT='(''SPACING'',3(E11.4,1x))') dx,dy, 0.0001d0
      WRITE(MyUnit,*) ' '
      WRITE(MyUnit,FMT='(''CELL_DATA '',I9)') Nx*Ny*1
      WRITE(MyUnit,*) ' '

      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'depth'
      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
     DO iy=1,Ny ; DO ix=1,Nx
        WRITE(MyUnit,'(G11.4)') Prim(H_pv,ix, iy)      
      ENDDO ; ENDDO 

      WRITE(MyUnit,FMT='(''VECTORS '',A12, '' float'')') 'Vitesse(m/s)'
      DO iy=1,Ny ; DO ix=1,Nx
        WRITE(MyUnit,'(3(E11.4,1x))') Prim(U_pv,ix, iy), Prim(V_pv,ix, iy), 0.d0
      ENDDO ; ENDDO 

      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P11'
      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
     DO iy=1,Ny ; DO ix=1,Nx
        WRITE(MyUnit,'(G11.4)') Prim(P11_pv,ix, iy)       
      ENDDO ; ENDDO 

      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P22'
      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
     DO iy=1,Ny ; DO ix=1,Nx
        WRITE(MyUnit,'(G11.4)') Prim(P22_pv,ix, iy)         
      ENDDO ; ENDDO 

      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P12'
      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
     DO iy=1,Ny ; DO ix=1,Nx
        WRITE(MyUnit,'(G11.4)') Prim(P12_pv,ix, iy)     
      ENDDO ; ENDDO 

     WRITE(MyUnit,*) ' '
   close(MyUnit)

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
 
   !  WRITE(unit=NB, fmt="(I3)") isave
   !  NB    = ADJUSTL(NB)
   !  il    = LEN_TRIM(NB) 
   !  WRITE(6,*) " FILE = ", "./resu/Powder2D_"//Zero(1:-il)//TRIM(NB)//".tp"
   !  OPEN(UNIT=MyUnit+2, FILE="./resu/Powder2D_"//Zero(1:3-il)//TRIM(NB)//".tp" )
   !  WRITE(MyUnit+2,'(A)') 'TITLE="This is a title"'  
   !  WRITE(MyUnit+2,'(A)') 'VARIABLES= "X", "Y" , "H","U", "V", "P11","P12", "P22" '
   !  WRITE(MyUnit+2,*) 'ZONE I=', NX,', J=', Ny,'DATAPACKING=POINT' 
   
   ! DO  ix=1, Nx
   !  DO  iy=1,Ny
       
   !        WRITE (MyUnit+2,'(8(E16.8,1x))') X(ix), Y(iy), Prim(H_pv,ix, iy),&
   !        & Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),prim(p11_pv,ix,iy), &
   !        &   Prim(p12_pv,ix, iy), prim(p22_pv,ix,iy)
   !     END DO
   !  END DO
 !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

      WRITE(6,*) " FILE = ", "./resu/Powder2D_"//Zero(1:3-il)//TRIM(NB)//".csv"
      OPEN(UNIT=MyUnit+6, FILE="./resu/Powder2D_"//Zero(1:3-il)//TRIM(NB)//".csv"  )
     
       WRITE(MyUnit+6,'(A)') '"X", "Y",  "H"'  
  
   DO  ix=1, Nx 
     DO  iy=1,Ny
             WRITE (MyUnit+6 ,*)  X(ix), ',' , Y(iy), ',' , &   
                 &         Prim(H_pv,ix, iy)*100
        END DO
   END DO

endif
 !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE LECTURE_DONNEES(Lx, Ly)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   REAL (KIND = DP) :: Lx, Ly
   OPEN(UNIT=21, FILE = 'dataOr1.inp', STATUS = 'OLD')
    READ(21,*) cond_lim               ! 1 box/ 2 absorbtion/ 3 batteur/ 4 jump/ 5 analyt solution
    READ(21,*) angle                  ! inclination angle
    READ(21,*) Nx, Ny                 ! NUMBER OF CELLS
    READ(21,*) Lx, Ly                 ! DOMAIN LENGTH 
    READ(21,*) TIMEOUT                ! OUTPUT TIME
    READ(21,*) iterfinal              ! Iteration final
    READ(21,*) g, CFL                 ! acceleration due to gravity
    READ(21,*) H_0                    ! stationary unstable solution (H_0, U_0, Phi_0 = 0)
    READ(21,*) phi2                   ! patit enstrophy
    READ(21,*) frottcoeff, disscoeff  ! cf, cr
    READ(21,*) amplitude              ! amplitude des perturbation
    READ(21,*) ImpreE, ImpreF
    READ(21,*) lambda, gamma, beta
    close(21)
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax, SOUND_bx,SOUND_ay, SOUND_by,Ein,CONS, Pression,  Lx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy
  REAL (KIND = DP) :: DY, DX, Lx, U_0
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim, 0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Pression(0:Nx+1,0:Ny+1), CONS(Nv_prim+1,1:Nx,1:Ny), Ein(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1), SOUND_bx(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) :: SOUND_ay(0:Nx+1,0:Ny+1), SOUND_by(0:Nx+1,0:Ny+1)

  pi = 4.0d0*ATAN(1.0d0)
  U_0 = DSQRT(g*Dtan(angle)*H_0/frottcoeff) 
 
DO ix = 1, Nx
X(ix) = 0.5D0*DX + (ix-1)*DX  
ENDDO

DO iy = 1, Ny
Y(iy) =  0.5D0*DY + (iy-1)*DY
ENDDO

DO ix = 1, Nx
DO iy = 1, Ny


if (cond_lim == 2) then 

Prim(H_pv, ix, iy) = H_0
Prim(U_pv,ix, iy) = 0.D0


Prim(Hp_pv, ix, iy) = Prim(H_pv, ix, iy) 

Prim(V_pv,ix, iy) =  0.d0
Prim(P11_pv,ix,iy) = eps
Prim(P12_pv,ix,iy) =  0.d0
Prim(P22_pv,ix,iy) =eps !H_0*H_0*phi2*(EPS)
else if (cond_lim == 1) then 
Prim(H_pv, ix, iy) = H_0*(1.d0 + amplitude*dsin(2.d0*Pi*x(ix)/Lx) )
Prim(Hp_pv, ix, iy) = Prim(H_pv, ix, iy) 
Prim(U_pv,ix, iy) = U_0
Prim(V_pv,ix, iy) =  0.d0
Prim(P11_pv,ix,iy) = Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)*phi2  
Prim(P12_pv,ix,iy) =  0.d0
Prim(P22_pv,ix,iy) = EPS

endif 

SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 &
&+ Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0


ENDDO
ENDDO


  ! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny 
      CONS(1,ix,iy) = Prim(H_pv, ix, iy)
      CONS(8,ix,iy) = Prim(Hp_pv, ix, iy)
      CONS(2,ix,iy) = Prim(H_pv, ix, iy)*Prim(U_pv, ix, iy)
      CONS(3,ix,iy) = Prim(H_pv, ix, iy)*Prim(V_pv, ix, iy)
      CONS(4,ix,iy) = Prim(H_pv, ix, iy)*Prim(P11_pv, ix, iy)
      CONS(5,ix,iy) = Prim(P12_pv, ix, iy)
      CONS(6,ix,iy) = Prim(H_pv, ix, iy)*Prim(P22_pv, ix, iy)
      CONS(7,ix,iy) = Prim(H_pv, ix, iy)*(Ein(ix,iy)+ &
        & (Prim(U_pv, ix, iy)**2.d0 + Prim(V_pv, ix, iy)**2.d0 )/2.d0)
  END DO
  END DO
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, &
    & SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   INTEGER :: ix,iy, it
  REAL (KIND = DP) :: p11, p12, p22, h
   REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), CONS(8,1:Nx,1:Ny)
  REAL (KIND = DP) ::  Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1), SOUND_bx(0:Nx+1,0:Ny+1)
      REAL (KIND = DP) ::  SOUND_ay(0:Nx+1,0:Ny+1), SOUND_by(0:Nx+1,0:Ny+1)

! CALCUL NOUVELS VARIABLES PRIMITIVES
DO ix = 1, Nx
DO iy = 1, Ny

  Prim(H_pv,ix, iy) = CONS(1,ix, iy)
  Prim(Hp_pv,ix, iy) = CONS(8,ix, iy)
  Prim(U_pv,ix, iy) = CONS(2,ix, iy)/CONS(1,ix, iy)
  Prim(V_pv,ix, iy) = CONS(3,ix, iy)/CONS(1,ix, iy)
!------------------------------------------------------------------------------------------------------------------------
  if (dabs(Prim(V_pv,ix, iy)).le.1.d-8) then
  Prim(V_pv,ix, iy)=0.d0; cons(3, ix, iy) = 0.d0
  endif
  if (dabs(Prim(U_pv,ix, iy)).le.1.d-8) then 
   Prim(U_pv,ix, iy)=0.d0; cons(2, ix, iy) = 0.d0
  endif
!------------------------------------------------------------------------------------------------------------------------
  Prim(P11_pv,ix, iy) = CONS(4,ix, iy)/CONS(1,ix, iy)
  Prim(P11_pv,ix, iy) = dmax1(Prim(P11_pv,ix, iy), 1.d-8)

  Prim(P12_pv,ix, iy) = CONS(5,ix, iy)
  Prim(P12_pv,ix, iy) = dmin1(Prim(P12_pv,ix, iy),&
  & dsqrt(Prim(P11_pv,ix, iy)*Prim(P22_pv,ix, iy)))
  Prim(P12_pv,ix, iy) = dmax1(Prim(P12_pv,ix, iy),&
  & -dsqrt(Prim(P11_pv,ix, iy)*Prim(P22_pv,ix, iy) ))
  if (dabs(Prim(p12_pv,ix, iy)).le.1.d-9) then 
     Prim(p12_pv,ix, iy)=0.d0; cons(5, ix, iy) = 0.d0
  endif


  Prim(P22_pv,ix, iy) = CONS(6,ix, iy)/CONS(1,ix, iy)
  Prim(P22_pv,ix, iy) = dmax1(Prim(P22_pv,ix, iy), 1.d-8)
!------------------------------------------------------------------------------------------------------------------------
  p11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
  h = Prim(h_pv,ix, iy)

   
!------------------------------------------------------------------------------------------------------------------------
  Ein(ix, iy) = CONS(7,ix, iy)/CONS(1,ix, iy) -&
  & ((Prim(U_pv,ix, iy))**2.d0+(Prim(V_pv,ix, iy))**2.d0 )/2.d0 !InternalEn(h, p11, p22) !
  SOUND_ax(ix, iy) = DSQRT( g*Prim(H_pv,ix, iy) + 3.d0*Prim(P11_pv,ix, iy))
  SOUND_bx(ix, iy) = DSQRT(Prim(P11_pv,ix, iy))  
  SOUND_ay(ix, iy) = DSQRT( g*Prim(H_pv,ix, iy) + 3.d0*Prim(P22_pv,ix, iy))
  SOUND_by(ix, iy) = DSQRT(Prim(P22_pv,ix, iy))  
  Pression(ix, iy) = g*(Prim(H_pv,ix, iy))**2.d0/2.d0 + Prim(P11_pv,ix, iy)*Prim(H_pv,ix, iy)
END DO  
END DO


  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1),Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) ::  SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)
!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
     Prim(:,0, iy)    = Prim(:,Nx, iy)

     Ein(0, iy)       = Ein(Nx, iy)
     SOUND_ax(0, iy)  = SOUND_ax(Nx, iy)
     SOUND_bx(0, iy)  = SOUND_bx(Nx, iy)
     SOUND_ay(0, iy)  = SOUND_ay(Nx, iy)
     SOUND_by(0, iy)  = SOUND_by(Nx, iy)
     Pression(0, iy)  = Pression(Nx, iy)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
     Prim(:,Nx+1, iy)   = Prim(:,1, iy)

     Ein(Nx+1, iy)      = Ein(1, iy)
     SOUND_ax(Nx+1, iy)  = SOUND_ax(1, iy)
     SOUND_bx(Nx+1, iy)  = SOUND_bx(1, iy)
     SOUND_ay(Nx+1, iy)  = SOUND_ay(1, iy)
     SOUND_by(Nx+1, iy)  = SOUND_by(1, iy)
     Pression(Nx+1, iy)  = Pression(1, iy)
    enddo   
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix, 0)    = Prim(:,ix, 1)
     Prim(V_pv,ix, 0)    = - Prim(V_pv,ix, 1)
     Prim(p12_pv,ix, 0)    = - Prim(p12_pv,ix, 1)

     Ein(ix, 0)       = Ein(ix, 1)
     SOUND_ax(ix, 0)  = SOUND_ax(ix, 1)
     SOUND_bx(ix, 0)  = SOUND_bx(ix, 1)
     SOUND_ay(ix, 0)  = SOUND_ay(ix, 1)
     SOUND_by(ix, 0)  = SOUND_by(ix, 1)
     Pression(ix, 0)  = Pression(ix, 1)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix, Ny+1)   = Prim(:,ix, Ny)
     Prim(V_pv,ix, Ny+1)    = - Prim(V_pv,ix, Ny)
     Prim(p12_pv,ix, Ny+1)    = - Prim(p12_pv,ix, Ny)

     Ein(ix, Ny+1)       = Ein(ix, Ny)
     SOUND_ax(ix, Ny+1)  = SOUND_ax(ix, Ny)
     SOUND_bx(ix, Ny+1)  = SOUND_bx(ix, Ny)
     SOUND_ay(ix, Ny+1)  = SOUND_ay(ix, Ny)
     SOUND_by(ix, Ny+1)  = SOUND_by(ix, Ny)
     Pression(ix, Ny+1)  = Pression(ix, Ny)
    enddo
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  
  SUBROUTINE CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) ::  SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) ::  SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  


!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
      Prim(:,0,  iy)  = Prim(:,1, iy)

      Ein(0,     iy)  = Ein(1, iy)
      SOUND_ax(0, iy) = SOUND_ax(1,  iy)
      SOUND_bx(0, iy) = SOUND_bx(1,  iy)
      SOUND_ay(0, iy) = SOUND_ay(1,  iy)
      SOUND_by(0, iy) = SOUND_by(1,  iy)
      Pression(0,iy)  = Pression(1, iy)

   ! Prim(H_pv,0, iy) = H_0

    !Prim(Hp_pv,0, iy)   = H_0
    !Prim(U_pv,0, iy)   = Prim(U_pv,1, iy) !0.D0
    !Prim(V_pv,0, iy)   = 0.d0
    !Prim(P11_pv,0, iy) = H_0*H_0*phi2 +eps 
    !Prim(P12_pv,0, iy) = Prim(P12_pv,1, iy) ! 0.d0
   ! Prim(P22_pv,0, iy) =  eps
   Ein(0, iy)         =  ( Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
   SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
   SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
   SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
   SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
   Pression(0, iy)    = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
     Prim(:,Nx+1,  iy) = Prim(:,Nx, iy)

     Ein(Nx+1,     iy)  = Ein(Nx, iy)
     SOUND_ax(Nx+1, iy) = SOUND_ax(Nx, iy)
     SOUND_bx(Nx+1, iy) = SOUND_bx(Nx, iy)
     SOUND_ay(Nx+1, iy) = SOUND_ay(Nx, iy)
     SOUND_by(Nx+1, iy) = SOUND_by(Nx, iy)
     Pression(Nx+1,iy)  = Pression(Nx,iy)
    enddo   
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  0)  = Prim(:,ix, 1)

     Ein(ix,     0)  = Ein(ix, 1)
     SOUND_ax(ix, 0) = SOUND_ax(ix,  1)
     SOUND_bx(ix, 0) = SOUND_bx(ix,  1)
     SOUND_ay(ix, 0) = SOUND_ay(ix,  1)
     SOUND_by(ix, 0) = SOUND_by(ix,  1)
     Pression(ix,0)  = Pression(ix, 1)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  Ny+1) = Prim(:,ix, Ny)

     Ein(ix,     Ny+1) = Ein(ix,     Ny)
     SOUND_ax(ix, Ny+1) = SOUND_ax(ix, Ny)
     SOUND_bx(ix, Ny+1) = SOUND_bx(ix, Ny)
     SOUND_ay(ix, Ny+1) = SOUND_ay(ix, Ny)
     SOUND_by(ix, Ny+1) = SOUND_by(ix, Ny)
     Pression(ix,Ny+1)  = Pression(ix,Ny)
    enddo
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:


  SUBROUTINE Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,TIME)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: q_0,u_0, yi, omega, TIME
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) :: SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  

    pi = 4.0d0*ATAN(1.0d0)
    omega = 6.19012d0
!------------------------------------------------------------------------------------------------------------------------
   if (angle .le. 1.d-8 .and. frottcoeff .le. 1.d-8 .and. disscoeff .le. 1.d-8) then 
      u_0 = 3.d0
   else 
      u_0 = DSQRT(H_0*g*Dtan(angle)/frottcoeff)
   end if 

   q_0 = H_0*u_0  
!------------------------------------------------------------------------------------------------------------------------
  do iy = 1, Ny
     Prim(H_pv,0, iy) = H_0*(1.d0 + amplitude*DSIN(omega*TIME))

    IF(Ny> 1 ) THEN
      yi =  REAL(iy-0.5)/REAL(Ny) 
      Prim(H_pv,0, iy) = H_0*(1.d0 + amplitude*DSIN(omega*TIME)+&
      & amplitude*SIN(2.0*pi*yi)*SIN(omega*Time/2.0))
    ENDIF
   Prim(Hp_pv,0, iy)   = Prim(H_pv,0, iy)
   Prim(U_pv,0, iy)   = q_0/Prim(H_pv,0, iy) 
   Prim(V_pv,0, iy)   = 0.d0 
   Prim(P11_pv,0, iy) = H_0*H_0*phi2*(1.d0 + EPS) 
   Prim(P12_pv,0, iy) = 0.d0
   Prim(P22_pv,0, iy) =  EPS  
   Ein(0, iy)         =  ( Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
   SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
   SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
   SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
   SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
   Pression(0, iy)    = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
  enddo
!------------------------------------------------------------------------------------------------------------------------
  do iy = 1, Ny
   Prim(H_pv,Nx+1, iy)   = Prim(H_pv,Nx, iy)
   Prim(Hp_pv,Nx+1, iy)= Prim(Hp_pv,Nx, iy)
   Prim(U_pv,Nx+1, iy)   = Prim(U_pv,Nx, iy)
   Prim(V_pv,Nx+1, iy)   = Prim(V_pv,Nx, iy) 
   Prim(P11_pv,Nx+1, iy) = Prim(P11_pv,Nx, iy)
   Prim(P12_pv,Nx+1, iy) = Prim(P12_pv,Nx, iy)
   Prim(P22_pv,Nx+1, iy) = Prim(P22_pv,Nx, iy)
   Ein(Nx+1, iy)         = Ein(Nx, iy)
   SOUND_ax(Nx+1, iy)    = SOUND_ax(Nx, iy)
   SOUND_bx(Nx+1, iy)    = SOUND_bx(Nx, iy)
    SOUND_ay(Nx+1, iy)   = SOUND_ay(Nx, iy)
   SOUND_by(Nx+1, iy)    = SOUND_by(Nx, iy)
   Pression(Nx+1, iy)    = Pression(Nx, iy)
  enddo
!------------------------------------------------------------------------------------------------------------------------
  do ix = 1, Nx
   Prim(H_pv,ix, 0)   = Prim(H_pv,ix, Ny)
   Prim(Hp_pv,ix, 0)= Prim(Hp_pv,ix, Ny)
   Prim(U_pv,ix, 0)   = Prim(U_pv,ix, Ny)
   Prim(V_pv,ix, 0)   = Prim(V_pv,ix, Ny) 
   Prim(P11_pv,ix, 0) = Prim(P11_pv,ix,Ny)
   Prim(P12_pv,ix, 0) = Prim(P12_pv,ix,Ny)
   Prim(P22_pv,ix, 0) = Prim(P22_pv,ix,Ny)
   Ein(ix, 0)         = Ein(ix, Ny)
   SOUND_ax(ix, 0)     = SOUND_ax(ix, Ny)
   SOUND_bx(ix, 0)     = SOUND_bx(ix, Ny)
   SOUND_ay(ix, 0)     = SOUND_ay(ix, Ny)
   SOUND_by(ix, 0)     = SOUND_by(ix, Ny)
   Pression(ix, 0)     = Pression(ix, Ny)

   Prim(H_pv,ix, Ny +1)   = Prim(H_pv,ix, 1)
   Prim(Hp_pv,ix, Ny +1)= Prim(Hp_pv,ix, 1)
   Prim(U_pv,ix, Ny +1)   = Prim(U_pv,ix, 1)
   Prim(V_pv,ix, Ny +1)   = Prim(V_pv,ix, 1) 
   Prim(P11_pv,ix, Ny +1) = Prim(P11_pv,ix, 1)
   Prim(P12_pv,ix, Ny +1) = Prim(P12_pv,ix, 1)
   Prim(P22_pv,ix, Ny +1) = Prim(P22_pv,ix, 1)
   Ein(ix, Ny +1)         = Ein(ix, 1)
   SOUND_ax(ix, Ny +1)    = SOUND_ax(ix, 1)
   SOUND_bx(ix, Ny +1)    = SOUND_bx(ix, 1)
    SOUND_ay(ix, Ny +1)   = SOUND_ay(ix, 1)
   SOUND_by(ix, Ny +1)    = SOUND_by(ix, 1)
   Pression(ix, Ny +1)    = Pression(ix, 1)

   end do
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE euler_method(DT, time, Lx, X,Y, CONS, Prim, Ein, &
    & SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, iy, IT,k
  REAL (KIND = DP) :: DT, traceP, H, p11, p22, p12, u, v,Q, alpha, detP
  REAL (KIND = DP) :: Vair, Vcore, rhoAir, rhoSnow, rhoPowder, rhop, Hp
  REAL (KIND = DP) :: injection_density, injection_speed, Lx
  REAL (KIND = DP) :: CONS(8,1:Nx, 1:Ny), U1(1:Nx, 1:Ny)
  REAL (KIND = DP) :: Prim(7,0:Nx+1,0:Ny+1), X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1),SOUND_bx(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ay(0:Nx+1,0:Ny+1),SOUND_by(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: fracp11,fracp22, erreur, time
  REAL (KIND = DP), ALLOCATABLE :: TS(:)

  Allocate(TS(Nv_Prim+1))
  Q = 0.d0; 
  rhoAir = 1.3d0; rhoSnow =10.d0;   TS(:) = 0.d0 
  injection_speed = 35.d0;  injection_density = 7.d0;    
  rhoPowder =10.D0 
 

 
!------------------------------------------------------------------------------------------------------------------------
  DO  ix = 1, Nx
    DO iy  = 1, Ny  
      
!------------------------------------------------------------------------------------------------------------------------ 
        if (cons(1,ix, iy).le.1.d-8) then
            print*, 'pas de l eau', cons(1,ix,iy), ix, iy, it
            stop
        end if
!------------------------------------------------------------------------------------------------------------------------

       !TERME SOURCE
  
    H = Prim(H_pv,ix, iy); hP = Prim(hp_pv, ix, iy)
    P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
    traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
    detP = p11*p22-p12**2.d0
    fracp11=p11/(P11+ P22);fracp22=p22/(P11+ P22);
   ! alpha = disscoeff/traceP
    !alpha = disscoeff*((traceP/(H**2.d0) - phi2 ))/( traceP**2.d0/(H**2.d0)  )
    alpha = disscoeff*(traceP - (H**2.d0)*phi2)/( traceP**2.d0)
    alpha = dmax1(alpha, 0.d0)
    
    if (cond_lim == 2) then 

    Vair = 0.d0 !dabs(u)*rhoAir/rhoPowder !0.2d0*dabs(u)*rhoAir/rhoPowder
    
  
! if (x(ix) >= Lx/4.d0  .and.  x(ix) .le. Lx/2.d0) then
 !   U1(ix,iy) = (2.d0*time)/(1.d0 + time**2.d0)*&
  !  &dexp(-(x(ix) - ( Lx/4.d0  + 10.d0)  -&
   ! & (35.D0)*time)**2.D0/(2.d0*(5.d0)**2.D0))  
   !endif 

  
      if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    U1(ix, iy) = 0.4d0*(dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e


    ! if ( x(ix) >= Lx/6.d0  .and.  x(ix) .le. 3.5D0*Lx/4.d0) then
    !U1(ix) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    !&dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    !& (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0)) ! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))e
   
   else

     U1(ix, iy) = 0.D0
   endif 



   rhoP = rhoAir + H*RhoPowder/hP
   
    if (  rhoP .gt. 2.d0*rhoAir .and. u .gt. 1.d-4) then 
     Vair =  (0.1d0*U + 0.9d0*dsqrt(phi2*H**2.d0 +detP)  )*H/hP
    else 
      !Vair = 0.25d0*dabs(u)*(rhoP-rhoAir)*rhoAir/rhoPowder
     Vair =  0.d0
    endif 

    if ( dabs(U1(ix, iy)) .le. eps) then 
      Vcore = 0.D0
    else

    Vcore = injection_speed*injection_density/rhoPowder
   end if 

 else 
  Vair = 0.d0
  Vcore = 0.d0
  U1(ix, iy) = 0.d0 
 endif 
    !Vcore = abs(U)*sqrt(rhoPowder*rhoSnow)/(rhoPowder+rhoSnow)

    Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)

  TS(1) = Vcore*U1(ix,iy)
  TS(8) = Vcore*U1(ix,iy) + Vair
  TS(2) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u  &
  &+ 5.d0*U1(ix,iy)**2.d0*Vcore!  ! g*Dtan(angle)*H  
  TS(3) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
  TS(4) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
  TS(5) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
  TS(6) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
  TS(7) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0)&
  & - Q + 5.d0*(U1(ix,iy)**2.d0)*Vcore*U + &
  & Vcore*U1(ix,iy)*(g*h -U**2.d0/2.d0 +3.d0/2.d0*(phi2*h**2.d0+detP))  ! g*Dtan(angle)*H*U 



    ! TS(1) = Vcore
    ! TS(8) = Vcore +Vair  !TS(8) = rhoAir*Vair + rhoPowder*Vcore
    ! TS(2) =  g*Dtan(angle)*H  - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u + U1(ix,iy)*Vcore !  g*Dtan(angle)*H 
    ! TS(3) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
    ! TS(4) = - 2.d0*alpha*p11*(dsqrt(u**2.d0 + v**2.d0 ))**(3.d0) 
    ! TS(5) =- 2.d0*alpha/H*p12*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) 
    ! TS(6) =- 2.d0*alpha*p22*(dsqrt(u**2.d0 + v**2.d0 ))**(3.d0) 
    ! TS(7) = g*Dtan(angle)*H*U  - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) &
    ! &- Q + U1(ix,iy)*Vcore*U  + Vcore*(g*h -U**2.d0/2.d0)  ! g*Dtan(angle)*H*U  
!------------------------------------------------------------------------------------------------------------------------
    do k=1,8
      CONS(k,ix, iy) = CONS(k,ix, iy) + DT*TS(k) 
    end do 


!     erreur=2.d0*CONS(7,ix, iy)/h-(CONS(2,ix, iy)**2+&
!   &CONS(3,ix, iy)**2)/h**2-g*h-(CONS(4,ix, iy)+ CONS(6,ix, iy))/h 

! !pause
! if (erreur .gt. eps) then 
!   print*, 'erreur =', erreur
!  ts(4)= - dt*fracp11*(TS(2)**2.d0 +TS(3)**2.d0 )/h                               !fracp11*erreur
!  ts(6)= - dt*fracp22*(TS(2)**2.d0 +TS(3)**2.d0 )/h                               !fracp22*erreur
!  CONS(4,ix, iy) = CONS(4,ix, iy) + dt*TS(4) !!verifier le signe
!  CONS(6,ix, iy) = CONS(6,ix, iy) + dt*TS(6) 
! endif 



!------------------------------------------------------------------------------------------------------------------------

      
  ENDDO
  ENDDO

  DEALLOCATE(TS)
  return
  END SUBROUTINE

  END PROGRAM code2D


   
