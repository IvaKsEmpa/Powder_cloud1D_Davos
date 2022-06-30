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
    INTEGER  ::  H_pv, U_pv, V_pv, Hp_pv
    REAL (KIND = DP)  ::  X0, Y0, H_0, amplitude,  frottcoeff, phi2
    REAL (KIND = DP)  ::  CFL, TIMEOUT, period_time, pi, angle, g
    REAL (KIND = DP)  ::  lambda, gamma, beta, disscoeff
    REAL (KIND = DP), PARAMETER ::  EPS = 1.d-8
    REAL (KIND = DP), PARAMETER ::  rhoAir = 1.3d0, rhoSnow =10.d0
    REAL (KIND = DP), PARAMETER ::  injection_speed = 35.d0, injection_density = 7.d0
    REAL (KIND = DP), PARAMETER ::  rhoPowder =10.D0 

  END MODULE GlobalParam
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      

MODULE ModeleInterface
  USE GlobalParam
  USE precisions
  IMPLICIT NONE

CONTAINS
!--------------------------------------------------------
  
 FUNCTION InternalEn(hm, hp)  RESULT(InternalE)
 REAL (KIND = DP)      :: hm, hp, InternalE
  InternalE =  0.5d0*(g*hm )
 END FUNCTION InternalEn
!--------------------------------------------------------

 FUNCTION Press(hm, hp) RESULT(pres)
 REAL (KIND = DP)      :: hm,hp,  pres
  pres =  0.5d0*g*hm**2.d0 
 END FUNCTION Press
!--------------------------------------------------------

FUNCTION Sound_a_x(hm, hp) RESULT(Sound_a)
 REAL (KIND = DP)      :: hm, hp, Sound_a
  Sound_a =  dsqrt(g*hm)
 END FUNCTION Sound_a_x
! !--------------------------------------------------------

END MODULE ModeleInterface


  PROGRAM code2D
  USE precisions
  use GlobalParam
  USE ModeleInterface
    IMPLICIT NONE  
    INTEGER :: iv, I, IT, ix, iy
    REAL (KIND = DP), ALLOCATABLE :: Ein(:,:), Sound_ax(:,:),  Pression(:,:)
    REAL (KIND = DP), ALLOCATABLE :: CONS(:,:,:), FLUX(:,:,:), Prim(:,:,:)
    REAL (KIND = DP), ALLOCATABLE :: UmaxTampon(:,:), VmaxTampon(:,:)
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:) :: MaxVP, MinVp, X, Y
    REAL (KIND = DP)  :: T1_CPU, T2_CPU, TIME,  TIME2
    REAL (KIND = DP) :: Lx, Ly,DX, Dy, Dh, DT, dt2
    REAL (KIND = DP) :: Hstar, Pstar, Ustar, Estar, Cmax
    REAL(KIND=DP):: xmax, XMIN, XMAX1, YMIN, YMAX, UMAX, VMAX
    
    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: NamesOfPrimVar
  !----------------------------------------------------------------------------------------------------
    pi = 4.0d0*ATAN(1.0d0)
    Nv_Prim  = 4
    H_pv     = 1
    U_pv     = 2
    V_pv     = 3
    Hp_pv=4
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
   ALLOCATE( Ein(0:Nx+1,0:Ny+1), Sound_ax(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1) )  
   ALLOCATE( Prim(Nv_Prim,0:Nx+1,0:Ny+1), CONS(Nv_Prim,1:Nx,1:Ny), FLUX(Nv_Prim,0:Nx,0:Ny))  
   ALLOCATE( NamesOfPrimVar(Nv_Prim), UmaxTampon(1:Nx,1:Ny), VmaxTampon(1:Nx,1:Ny) )
   ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim), X(1:Nx), Y(1:Ny))
   flux(:,:,:)=0.d0; cons = 0.d0; prim = 0.d0; Sound_ax = 0.d0
   X = 0.d0; Y = 0.D0
  !----------------------------------------------------------------------------------------------------
    NamesOfPrimVar(H_pv)         = "Depht Lenght"
    NamesOfPrimVar(U_pv)         = "Velocity (x)"
    NamesOfPrimVar(V_pv)         = "Velocity (y)"
    NamesOfPrimVar(Hp_pv) = " True Height"
  !----------------------------------------------------------------------------------------------------
    DO iv = 1, Nv_Prim
       WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx, Ny
    END DO
    WRITE(6,*) " >>> End of LECTURE_DONNEES"
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
  CALL INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax, Ein,CONS, Pression,  Lx)
  !----------------------------------------------------------------------------------------------------
   TIME = 0.D0
   period_time = 1.0d0;
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

    CALL Ecriture_donnees(X,Y, Lx, Prim, Pression, time, DX, DY )
 !stop
   !-------------------------do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT)!-------------------
   ! BOUCLE SUR LE TEMPS
    OPEN(argunit+1,FILE = './resu/Umax.out')


   Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))
   !----------------------------------------------------------------------------------------------------
!print*, " i am here 2"
    do ix = 1, Nx
     do iy = 1, Ny
      UmaxTampon(ix, iy) = dmax1( Dabs(Prim(U_pv,ix,iy)+ &
        & Sound_ax(ix,iy)), Dabs(Prim(U_pv,ix,iy)- Sound_ax(ix,iy)) )
     enddo
    enddo
    Umax = MAXVAL( UmaxTampon(1:Nx,1:Ny))
   
   !----------------------------------------------------------------------------------------------------    
   
   
    DT   = Dx/DABS(Umax)

    DT   = CFL*DT
    dt2 = 0.5d0*dt
   !----------------------------------------------------------------------------------------------------
    IF( It == 1 ) WRITE(6,*) " Dt = ", Dt
    TIME = TIME + DT
!    !----------------------------------------------------------------------------------------------------
   
  !----------------------------------------------------------------------------------------------------

    CALL CondLimABSORPTION(Prim, Pression, SOUND_ax)
   
  !----------------------------------------------------------------------------------------------------
   call HLLC_x_sub1(prim,flux) 
  !----------------------------------------------------------------------------------------------------
   call godunov_x_sub1(cons,flux,dt,dx)
  !----------------------------------------------------------------------------------------------------
   CALL NOUVELL_VARIABLE_PRIM(Prim, SOUND_ax, CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------
  
  call euler_method(DT, time, Lx,X,Y, CONS, Prim, SOUND_ax,Pression, it)
   
  !----------------------------------------------------------------------------------------------------

   CALL NOUVELL_VARIABLE_PRIM(Prim, SOUND_ax, CONS, Pression, it)

  !----------------------------------------------------------------------------------------------------
    IT = IT + 1
  !----------------------------------------------------------------------------------------------------
     IF (TIME2.LE.TIME) THEN
     	CALL Ecriture_donnees(X,Y, Lx, Prim, Pression, time, DX, DY )
       WRITE(argunit+1,'(2(E20.13,1X))') time, MAXVAL( Prim(U_pv,1:Nx,1:Ny))
    
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

   CALL Ecriture_donnees(X,Y, Lx,  Prim, Pression, time, DX, DY )
  !----------------------------------------------------------------------------------------------------
   DEALLOCATE(Prim, Ein, SOUND_ax, X, Y, Pression, CONS, FLUX, MinVP, MaxVp, NamesOfPrimVar)  
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
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), FLUX(Nv_Prim,0:Nx,0:Ny)
  real(kind=dp) :: ul,ur,hl,hr,pr,pl, cl,cr,ml,mr,sl,sr
  real(kind=dp) :: upstar, mpl, mpr, hpstar, ppstar
  real(kind=dp) :: vr,vl, HPl, HPr, cpl, cpr, spr, spl 
  real(kind=dp) :: pstar,ustar,hstar,p12star

  do ix=0,nx
   do iy=1,ny
    !! etat gauche et droite
    ul=prim(U_pv,ix,iy);  ur=prim(U_pv,ix+1,iy)
    hl=prim(h_pv,ix,iy);  hr=prim(h_pv,ix+1,iy)
    vl=prim(v_pv,ix,iy);  vr=prim(v_pv,ix+1,iy)

    HPl=prim(Hp_pv,ix,iy)
    HPr=prim(Hp_pv,ix+1,iy)

    !! calcul des pression
    pl=g*hl*hl*0.5d0 ! (hl/hpl)*
    pr=g*hr*hr*0.5d0  ! (hr/hpr)*
    !! calcul des vitesses du son
    cl=dsqrt(g*hl);cr=dsqrt(g*hr)  ! (hl/hpl)*  ! (hr/hpr)*

    cpl=dsqrt(g*hpl);cpr=dsqrt(g*hpr) 

    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl


    Spl=ul-cpl; if (ur-cpr < spl) spl = ur - cpr
    spr=ur+cpr; if(ul+cpl> spr) spr = ul+cpl

    ml=hl*(ul-sl)
    mpl=hpl*(ul-spl)

    mr=hr*(ur-sr)
    mpr=hpr*(ur-spr)

    ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
    pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)

    upstar=(ul*mpl-ur*mpr+pl-pr)/(mpl-mpr)
    ppstar=((ul-ur)*mpr*mpl+mpr*pl-mpl*pr)/(mpr-mpl)

  if (ustar.ge.0.d0) then
   if (sl.ge.0.d0) THEN
  !!etat gauche
   flux(1,ix,iy)=hl*ul
   flux(4,ix,iy)=HPl*ul
   flux(2,ix,iy)=hl*ul*ul+pl
   flux(3,ix,iy)=hl*ul*vl
   
   ELSE
   !! etat star gauche
   hstar=ml/(ustar-sl)
   hpstar = mpl/(upstar-spl)
  ! p12star=p12l*(ul-sl)/(ustar-sl)
 
   !! remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(4,ix,iy)=hpstar*upstar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vl!hl*ul*vl + sl*vl*(hstar - hl)   
   
   
   endif
  ELSE
   if (sr.ge.0.d0) then
  !!etat droit etoile

   hstar=mr/(ustar-sr)

   hpstar=mpr/(upstar -spr)

   !remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(4,ix,iy)=hpstar*upstar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vr!hl*ul*vl + sl*vl*(hstar - hl)   
   
   
   ELSE
  !!etat droit
   flux(1,ix,iy)=hr*ur
   flux(4,ix,iy)=hpr*ur
   flux(2,ix,iy)=hr*ur*ur+pr
   flux(3,ix,iy)=hr*ur*vr
   
   end if
  end if

  end do
  end do
  return
  end subroutine HLLC_x_sub1
  
  subroutine godunov_x_sub1(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(Nv_Prim,1:nx,1:ny), FLUX(Nv_Prim,0:Nx,0:Ny)
  real(Kind=dp) :: hu, h, hv, he, dt, dx
  INTEGER :: k, ix, iy

  do ix=1,nx
   do iy=1,ny
!------------------------------------------------------------------------------------------------------------------------
    do k=1,Nv_Prim
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
   end do
!------------------------------------------------------------------------------------------------------------------------
  end do
  end do
  end subroutine godunov_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

  SUBROUTINE Ecriture_donnees(X,Y, Lx, Prim, Pression, time, DX, DY )
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy , MyUnit = 30, il
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny),  U1(1:Nx, 1:Ny)
  REAL (KIND = DP) :: h,u,v, Lx, rhoP
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  
  REAL (KIND = DP) :: time, DX, DY
  CHARACTER(LEN=3) :: NB, Zero="000"
  
  U1 = 0.D0

  ! ECRITURE DES RESULTATS

        OPEN(MyUnit+1,FILE = './resu/SW628.txt')

  Do ix = 1, Nx
  DO iy = 1, Ny

  if ( x(ix) >= Lx/6.d0 ) then
    U1(ix,iy) = (0.4D0*dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0))! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))
  endif 

   if (dabs(U1(ix, iy)) .le. 1.d-8) then
    U1(ix, iy) = 0.D0
  endif 

  !U1(ix, iy) = 0.D0

   
   rhoP = rhoAir + Prim(h_pv,ix, iy)*RhoPowder/Prim(Hp_pv, ix, iy)

      WRITE(MyUnit+1,'(8(E20.13,1X))') X(ix), Y(iy), Prim(h_pv,ix, iy), Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),&
      &U1(ix, iy), Prim(Hp_pv, ix, iy), rhoP
  END DO
  END DO
       WRITE(MyUnit+1,*)
       print*, 'Fr =', MAXVAL(dabs(Prim(U_pv,1:Nx, 1:Ny))/Sound_ax(1:Nx, 1:Ny))
    

 ! if (ny > 1) then

 !        isave = isave + 1
 !       WRITE(unit=NB, fmt="(I3)") isave
 !       NB    = ADJUSTL(NB)
 !       il    = LEN_TRIM(NB) 


 !      WRITE(6,*) " FILE = ", "./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".vtk"
 !      OPEN(UNIT=MyUnit, FILE="./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".vtk" )
 !      WRITE(MyUnit,'(''# vtk DataFile Version 2.0'')')
 !      WRITE(MyUnit,'(''Rectilinear 3D Dataset'')')
 !      WRITE(MyUnit,'(''ASCII'')')
 !      WRITE(MyUnit,'(''           '')')
 !      WRITE(MyUnit,'(''DATASET STRUCTURED_POINTS'')')
 !      WRITE(MyUnit,FMT='(''DIMENSIONS'',I8,I8,I8)') Nx+1, Ny+1, 2
 !      WRITE(MyUnit,FMT='(''ORIGIN'',3(E11.4,1x))') dx, 0.d0, 0.d0
 !      WRITE(MyUnit,FMT='(''SPACING'',3(E11.4,1x))') dx,dy, 0.0001d0
 !      WRITE(MyUnit,*) ' '
 !      WRITE(MyUnit,FMT='(''CELL_DATA '',I9)') Nx*Ny*1
 !      WRITE(MyUnit,*) ' '

 !      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'depth'
 !      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
 !     DO iy=1,Ny ; DO ix=1,Nx
 !        WRITE(MyUnit,'(G11.4)') Prim(H_pv,ix, iy)      
 !      ENDDO ; ENDDO 

 !      WRITE(MyUnit,FMT='(''VECTORS '',A12, '' float'')') 'Vitesse(m/s)'
 !      DO iy=1,Ny ; DO ix=1,Nx
 !        WRITE(MyUnit,'(3(E11.4,1x))') Prim(U_pv,ix, iy), Prim(V_pv,ix, iy), 0.d0
 !      ENDDO ; ENDDO 

 !      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P11'
 !      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
 !     DO iy=1,Ny ; DO ix=1,Nx
 !        WRITE(MyUnit,'(G11.4)') Prim(P11_pv,ix, iy)       
 !      ENDDO ; ENDDO 

 !      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P22'
 !      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
 !     DO iy=1,Ny ; DO ix=1,Nx
 !        WRITE(MyUnit,'(G11.4)') Prim(P22_pv,ix, iy)         
 !      ENDDO ; ENDDO 

 !      WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P12'
 !      WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
 !     DO iy=1,Ny ; DO ix=1,Nx
 !        WRITE(MyUnit,'(G11.4)') Prim(P12_pv,ix, iy)     
 !      ENDDO ; ENDDO 

  !    WRITE(MyUnit,*) ' '
   !close(MyUnit)

! !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
 
!     WRITE(unit=NB, fmt="(I3)") isave
!     NB    = ADJUSTL(NB)
!     il    = LEN_TRIM(NB) 
!     WRITE(6,*) " FILE = ", "./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".tp"
!     OPEN(UNIT=MyUnit+2, FILE="./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".tp" )
!     WRITE(MyUnit+2,'(A)') 'TITLE="This is a title"'  
!     WRITE(MyUnit+2,'(A)') 'VARIABLES= "X", "Y" , "H","U", "V", "P11","P12", "P22" '
!     WRITE(MyUnit+2,*) 'ZONE I=', NX,', J=', Ny,'DATAPACKING=POINT' 
   
!    DO  ix=1, Nx
!     DO  iy=1,Ny
       
!           WRITE (MyUnit+2,'(8(E16.8,1x))') X(ix), Y(iy), Prim(H_pv,ix, iy),&
!           & Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),prim(p11_pv,ix,iy), &
!           &   Prim(p12_pv,ix, iy), prim(p22_pv,ix,iy)
!        END DO
!     END DO
!  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

!       WRITE(6,*) " FILE = ", "./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".csv"
!       OPEN(UNIT=MyUnit+6, FILE="./resu/testCFL01_"//Zero(1:3-il)//TRIM(NB)//".csv"  )
     
!        WRITE(MyUnit+6,'(A)') '"X", "Y",  "H"'  
  
!    DO  ix=1, Nx 
!      DO  iy=1,Ny
!              WRITE (MyUnit+6 ,*)  X(ix), ',' , Y(iy), ',' , &   
!                  &         Prim(H_pv,ix, iy)*100.d0
!         END DO
!    END DO

!endif
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

SUBROUTINE INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax,Ein,CONS, Pression,  Lx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy 
  REAL (KIND = DP) :: DY, DX, Lx
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim, 0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Pression(0:Nx+1,0:Ny+1), CONS(Nv_Prim,1:Nx,1:Ny), Ein(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1)
  pi = 4.0d0*ATAN(1.0d0)
 
DO ix = 1, Nx
X(ix) = 0.5D0*DX + (ix-1)*DX  
ENDDO

DO iy = 1, Ny
Y(iy) =  0.5D0*DY + (iy-1)*DY
ENDDO

DO ix = 1, Nx
DO iy = 1, Ny

 !if (x(ix) >=  Lx/2.d0 .and.  x(ix) .le.  Lx/2.d0 + 50.D0 ) then
 
Prim(H_pv, ix, iy) = H_0
Prim(U_pv,ix, iy) = 0.D0


Prim(Hp_pv, ix, iy) = Prim(H_pv, ix, iy) 

Prim(V_pv,ix, iy) =  0.d0


SOUND_ax(ix,iy) = DSQRT((g*Prim(H_pv, ix, iy)  )) ! (Prim(H_pv, ix, iy)/Prim(Hp_pv, ix, iy))& &*

Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0  ! (Prim(H_pv, ix, iy)/Prim(Hp_pv, ix, iy))* 
Ein(ix,iy) = Prim(H_pv, ix, iy)*g/2.d0  ! (Prim(H_pv, ix, iy)/Prim(Hp_pv, ix, iy))*


ENDDO
ENDDO


  ! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny 
      CONS(1,ix,iy) = Prim(H_pv, ix, iy)
      CONS(4,ix,iy) = Prim(Hp_pv, ix, iy)
      CONS(2,ix,iy) = Prim(H_pv, ix, iy)*Prim(U_pv, ix, iy)
      CONS(3,ix,iy) = Prim(H_pv, ix, iy)*Prim(V_pv, ix, iy)
     
  END DO
  END DO
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim, SOUND_ax, CONS, Pression, it)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   INTEGER :: ix,iy, it
  REAL (KIND = DP) ::  h
   REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), CONS(Nv_Prim,1:Nx,1:Ny)
  REAL (KIND = DP) ::  Pression(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1)

! CALCUL NOUVELS VARIABLES PRIMITIVES
DO ix = 1, Nx
DO iy = 1, Ny

  Prim(H_pv,ix, iy) = CONS(1,ix, iy)
  Prim(Hp_pv,ix, iy) = CONS(4,ix, iy)
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
  
  h = Prim(h_pv,ix, iy)

   
!------------------------------------------------------------------------------------------------------------------------
  
  SOUND_ax(ix, iy) = DSQRT(g*Prim(H_pv,ix, iy) ) ! (Prim(H_pv, ix, iy)/Prim(Hp_pv, ix, iy))*
   
  Pression(ix, iy) =g*(Prim(H_pv,ix, iy))**2.d0/2.d0 !  (Prim(H_pv, ix, iy)/Prim(Hp_pv, ix, iy))* 
END DO  
END DO


  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  
  SUBROUTINE CondLimABSORPTION(Prim, Pression, SOUND_ax)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) ::  SOUND_ax(0:Nx+1,0:Ny+1) 


!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
    
   Prim(H_pv,0, iy) = Prim(H_pv,1, iy)

   Prim(Hp_pv,0, iy)   =  Prim(H_pv,1, iy)
   Prim(U_pv,0, iy)   =  Prim(U_pv,1, iy) 
   Prim(V_pv,0, iy)   =  Prim(V_pv,1, iy)
  
  
   SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy))  ! ( Prim(H_pv,0, iy)/Prim(Hp_pv,0, iy))* 
  
   Pression(0, iy)    = g*(Prim(H_pv,0, iy))**2.d0/2.d0   ! ( Prim(H_pv,0, iy)/Prim(Hp_pv,0, iy))*
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do iy = 1, Ny
     Prim(:,Nx+1,  iy) = Prim(:,Nx, iy)

     SOUND_ax(Nx+1, iy) = SOUND_ax(Nx, iy)
     
     Pression(Nx+1,iy)  = Pression(Nx,iy)
    enddo   
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  0)  = Prim(:,ix, 1)

     SOUND_ax(ix, 0) = SOUND_ax(ix,  1)
     
     Pression(ix,0)  = Pression(ix, 1)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  Ny+1) = Prim(:,ix, Ny)

     SOUND_ax(ix, Ny+1) = SOUND_ax(ix, Ny)
     
     Pression(ix,Ny+1)  = Pression(ix,Ny)
    enddo
  return
  END SUBROUTINE
  
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE euler_method(DT, time, Lx, X,Y, CONS, Prim, &
    & SOUND_ax,Pression, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, iy, IT,k
  REAL (KIND = DP) :: DT,  H,  u, v,Q, alpha
  REAL (KIND = DP) :: Vair, Vcore,  rhop, Hp
  REAL (KIND = DP) ::  Lx
  REAL (KIND = DP) :: CONS(Nv_Prim,1:Nx, 1:Ny), U1(1:Nx, 1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1), X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: erreur, time
  REAL (KIND = DP) :: TS(Nv_Prim)
  H= 0.d0; u= 0.d0; v= 0.d0
  U1 =0.D0;   TS = 0.d0 
  !injection_speed = 10.d0;  injection_density = 7.d0;    rhoPowder =10.D0 
  pi = 4.0d0*ATAN(1.0d0)
!------------------------------------------------------------------------------------------------------------------------
  DO  ix = 1, Nx
    DO iy  = 1, Ny  
      
!------------------------------------------------------------------------------------------------------------------------ 
        if (cons(1,ix, iy).le.1.d-10) then
            print*, 'pas de l eau', cons(1,ix,iy), ix, iy, it
            stop
        end if
!------------------------------------------------------------------------------------------------------------------------

       !TERME SOURCE
  
    H = Prim(H_pv,ix, iy); hP = Prim(hp_pv, ix, iy)
     u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
    

 if ( x(ix) >= Lx/6.d0 ) then
    U1(ix,iy) = (dsin(0.4d0*pi*time)**2.d0)*&
    &dexp(-(x(ix) - ( Lx/6.d0  + 10.d0)  -&
    & (35.D0)*time)**2.D0/(2.d0*(14.d0)**2.D0))! (35.D0*(1.d0 - x(ix)*2.d0/Lx ))
  endif 

   if (dabs(U1(ix, iy)) .le. 1.d-8) then
    U1(ix, iy) = 0.D0
  endif 
      
   rhoP = rhoAir + H*RhoPowder/hP
   

    if (  rhoP .gt.2.d0*rhoAir .and. u .gt. eps ) then 
     Vair = u*H/hP 
    else 
      Vair =  0.d0
    endif 


    if (dabs(U1(ix, iy)) .le. eps) then 
      Vcore = 0.D0
    else

    Vcore = injection_speed*injection_density/rhoPowder
   end if 

    !Vcore = abs(U)*sqrt(rhoPowder*rhoSnow)/(rhoPowder+rhoSnow)

    TS(1) = Vcore*U1(ix,iy)
    TS(4) = Vcore*U1(ix,iy)+Vair  !TS(8) = rhoAir*Vair + rhoPowder*Vcore
    TS(2) =  g*Dtan(angle)*H  - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u + 5.D0*U1(ix,iy)*U1(ix,iy)*Vcore  !  g*Dtan(angle)*H 
    TS(3) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v + 5.d0*(U1(ix,iy)**2.d0)*Vcore*u &
    &+ Vcore*U1(ix,iy)*(g*h -U**2.d0/2.d0)  ! 
    
!------------------------------------------------------------------------------------------------------------------------
    do k=1,Nv_Prim
      CONS(k,ix, iy) = CONS(k,ix, iy) + DT*TS(k) 
    end do 
      
  ENDDO
  ENDDO

  return
  END SUBROUTINE

  END PROGRAM code2D


   
