module riemann_parameter
  implicit none
  integer, parameter :: NN = 10004  ! number of cell
  real, parameter :: ep = 1.0E-5, PROT1 = 1.0E-4  ! zero-avoid
  real, parameter :: omg = 7.0/5.0 - 1.0  ! specific heat ratio of air
  real, parameter :: TSTOP = 0.16, EPS = 0.5, TAUMAX = 10.0
end module riemann_parameter

module functions
  use riemann_parameter
  implicit none
contains
  real function G1(x)
    real, intent(in) :: x
    real :: A,B,R1,R2,omg1,r0
    A=0.
    B=0.
    R1=1.
    R2=1.
    omg1=5./3.-1.
    r0=1.
    G1=A*(1-(omg1*x)/(R1*r0))*exp(-(R1*r0)/x) + &
    B*(1-(omg1*x)/(R2*r0))*exp(-(R2*r0)/x) 
    G1=-0.
  end function G1

  real function f1(x)
    real, intent(in) :: x
    real :: A, B, R1, R2, omg1, r0
    A=0.
    B=0.
    R1=1.
    R2=1.
    omg1=5./3.-1.
    r0=1.
    f1 =-(A*(1-(omg1*x)/(R1*r0))*exp(-(R1*r0)/x) + &
    B*(1-(omg1*x)/(R2*r0))*exp(-(R2*r0)/x))/x + &
    (A*(R1*r0-omg1*x)*exp(-(R1*r0)/x) + &
    B*(R2*r0-omg1*x)*exp(-(R2*r0)/x))/(x**2)
    f1 = 0./x
  end function f1
end module functions


PROGRAM riemann_solver
  use riemann_parameter
  use functions
  implicit none
  real :: DENS(NN), XMOM(NN), ENER(NN)  ! density, momentum & energy value of each cell
  real :: FM(NN), FP(NN), FE(NN)  ! density, momentum & energy flux of each cell
  real :: DD(NN), DX(NN), DE(NN)  ! variables used in reconstruction
  real :: uul(3), uur(3), fflux(3)  ! store the left and right states and final flow for riemann solver
  integer :: I, IMAX, IMAXM, IMAXM2, IT
  real :: TIME, DXX, TAU, UBD, TBD, XI, AA
  real :: den, vel, prs, enr  ! use for initial condition
  real :: air_vel, signal_vel
  real :: DFD1, DFX1, DFE1, DFD2, DFX2, DFE2, ADE1, AXM1, AEN1, ADE2, AXM2,AEN2
  integer :: solver_number, initail_number
  character(len=20) :: arg1, arg2
  real :: start_time, end_time, elapsed_time

  call cpu_time(start_time)


  ! choose solver to use: 1 for 'HLL'; 2 for 'HLLC'; 3 for 'Roe'
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  ! 检查是否成功获取参数
  if (trim(arg1) /= '') then
    ! 根据参数值执行相应的操作
    select case(trim(arg1))
    case ("roe")
      solver_number = 3
    case ("hll")
      solver_number = 1
    case ("hllc")
      solver_number = 2
    case default
      print *, "无效的参数。请传递 'roe', 'hll', 或 'hllc'。"
    end select
  else
    print *, "请提供第一个参数。"
  endif

  if (trim(arg2) /= '') then
    select case(trim(arg2))
    case ("1")
      initail_number = 1
    case ("2")
      initail_number = 2
    case ("3")
      initail_number = 3
    case ("4")
      initail_number = 4
    case default
      print *, "无效的参数。请传递 '1', '2', '3', '4'。"
    end select
  else
    print *, "请提供第二个参数。"
  endif


  ! define file output format
  88  FORMAT(4(F14.6, 2X))
  ! initialize variables
  IMAX = NN
  DXX = 1.0 / (IMAX - 4)
  TIME = 0.0
  AA = 1.0d0
  IT = 1

  IMAXM = IMAX - 1
  IMAXM2 = IMAX - 2

  ! give initial condition
  do I = 1, IMAX
    if ( FLOAT(I-2) < IMAX/2.0 ) then
      if ( initail_number.eq.1 ) then  ! shock tube problem
        den=0.445
        vel=0.698
        prs=3.528
      else if ( initail_number.eq.2 ) then  ! Sod problem
        den=1.0
        vel=0.0
        prs=1.0
      else if ( initail_number.eq.3 ) then  ! 123 problem
        den=1.0
        vel=2.0
        prs=0.4
      else 
        den=1.0
        vel=1.0
        prs=0.4
      end if
    else 
      if ( initail_number.eq.1 ) then
        den=0.5
        vel=0.0
        prs=0.571
      else if ( initail_number.eq.2 ) then
        den=0.125
        vel=0.0
        prs=0.1
      else if ( initail_number.eq.3 ) then
        den=1.0
        vel=-2.0
        prs=0.4
      else 
        den=1.0
        vel=-1.0
        prs=0.4
      end if
    end if
    enr=(prs-G1(den))/(den*omg)
    DENS(I)=den
    XMOM(I)=den*vel
    ENER(I)=0.5*den*(vel**2)+den*enr
  end do

  ! calculate each time step
  do while (TIME <= TSTOP)

    ! find the max signal velocity for the time step length
    UBD = 0.0
    do I = 3, IMAXM2
      den=DENS(I)
      vel=XMOM(I)/DENS(I)
	    prs=(ENER(I)-0.5*XMOM(I)**2/DENS(I))*omg+G1(den)
      ! dffd=f1(den)  ! dont know what it means
      air_vel = sqrt(f1(den)+(omg+1.)*prs/den)
      signal_vel = abs(vel)+abs(air_vel)
      UBD = MAX(UBD,signal_vel)
    end do
    TBD=1.0/(UBD + PROT1)
    TAU=DXX*EPS*MIN(TAUMAX,TBD)  ! TAU is the time step length
    TIME=TIME+TAU

    ! ! ---------------------------------------------------------
    ! ! debug
    !   if ( it == 2 ) then
    !     exit
    !   end if
    ! ! ---------------------------------------------------------



    if ( TIME > TSTOP ) then
      exit  ! jump off time loop to output data
    end if
    IT=IT+1
    ! WRITE(*,*) IT,TAU,TIME  ! print the time step in terminal

    ! reconstruction: use second-order MUSCL recon-scheme
    DO I= 2,IMAXM
      DFD1=DENS(I)-DENS(I-1)
      DFX1=XMOM(I)-XMOM(I-1)
      DFE1=ENER(I)-ENER(I-1)
                  
      DFD2=DENS(I+1)-DENS(I)
      DFX2=XMOM(I+1)-XMOM(I)
      DFE2=ENER(I+1)-ENER(I)
          
      DD(I)=(SIGN(AA,DFD1)+SIGN(AA,DFD2))* &
      ((ABS(DFD1)*ABS(DFD2))/(ABS(DFD1)+ABS(DFD2)+EP))  ! /DXX
      DX(I)=(SIGN(AA,DFX1)+SIGN(AA,DFX2))* &
      ((ABS(DFX1)*ABS(DFX2))/(ABS(DFX1)+ABS(DFX2)+EP))  ! /DXX
      DE(I)=(SIGN(AA,DFE1)+SIGN(AA,DFE2))* &
      ((ABS(DFE1)*ABS(DFE2))/(ABS(DFE1)+ABS(DFE2)+EP))  ! /DXX 
    ENDDO

    ! calculate the numerical flux
    DO I= 2,IMAXM-1
      ADE1=DENS(I) +0.25*DD(I)  ! *DXX
      AXM1=XMOM(I) +0.25*DX(I)  ! *DXX
      AEN1=ENER(I) +0.25*DE(I)  ! *DXX
         
      ADE2=DENS(I+1) -0.25*DD(I+1)  ! *DXX
      AXM2=XMOM(I+1) -0.25*DX(I+1)  ! *DXX
      AEN2=ENER(I+1) -0.25*DE(I+1)  ! *DXX

      uul(1)=ADE1
	    uul(2)=AXM1/ADE1
	    uul(3)=(AEN1-0.5*AXM1**2/ADE1)*omg+G1(ADE1)
	   
      uur(1)=ADE2
	    uur(2)=AXM2/ADE2
	    uur(3)=(AEN2-0.5*AXM2**2/ADE2)*omg+G1(ADE2)
      
      ! choose hllc/hll/roe riemann solver
      if ( solver_number == 1 ) then
        CALL hll(uul,uur,fflux)
      else if ( solver_number == 2 ) then
        CALL hllc(uul,uur,fflux)
      else
        CALL roe(uul,uur,fflux)
      end if
 	   
      FM(I)=fflux(1)
      FP(I)=fflux(2)
      FE(I)=fflux(3)
    ENDDO

    ! update all state
    DO I = 3, IMAXM-1
      DENS(I)=DENS(I)+(FM(I-1)-FM(I))/DXX*TAU
      XMOM(I)=XMOM(I)+(FP(I-1)-FP(I))/DXX*TAU
      ENER(I)=ENER(I)+(FE(I-1)-FE(I))/DXX*TAU
    ENDDO

    ! rewrite edge state
    DENS(1)=DENS(4)
    DENS(2)=DENS(3)
    DENS(IMAXM)=DENS(IMAXM2)
    DENS(IMAX)=DENS(IMAXM2-1) 

  	XMOM(1)=XMOM(4)
    XMOM(2)=XMOM(3)
    XMOM(IMAXM)=XMOM(IMAXM2)
    XMOM(IMAX)=XMOM(IMAXM2-1)

	  ENER(1)=ENER(4)
    ENER(2)=ENER(3)
    ENER(IMAXM)=ENER(IMAXM2)
    ENER(IMAX)=ENER(IMAXM2-1) 

  end do

  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  print *, "Elapsed time (seconds): ", elapsed_time

  ! write to file 'riemann.plt'
  if ( solver_number == 1 ) then
    OPEN(UNIT=1,FILE='result/hll.plt',STATUS='unknown')
  else if ( solver_number == 2 ) then
    OPEN(UNIT=1,FILE='result/hllc.plt',STATUS='unknown')
  else
    OPEN(UNIT=1,FILE='result/roe.plt',STATUS='unknown')
  end if
  ! WRITE(1,*) 'IT = ', IT
  WRITE(1,*) 'title="contour"'
  WRITE(1,*) 'variables="x","den","vel","pre"'
  DO I=3,IMAXM2
    XI=DBLE(I-2)*DXX
       
    den=DENS(I)
    vel=XMOM(I)/DENS(I)
	  prs=(ENER(I)-0.5*XMOM(I)**2/DENS(I))*omg+G1(DENS(I))
    WRITE(1,88)  XI+0.045,den,vel,prs
  ENDDO 

end PROGRAM riemann_solver



subroutine hll(ul,ur,flux)
  
  use riemann_parameter
  use functions
  implicit none

  real, intent(in) :: ul(3),ur(3)
  real, intent(out) :: flux(3)
  real :: fl(3),fr(3),uf(3),uur(3),uul(3)
  real :: aal, aar, enerl, enerr, SL, SR

  aal=sqrt(f1(ul(1))+(omg+1.)*ul(3)/ul(1))
	aar=sqrt(f1(ur(1))+(omg+1.)*ur(3)/ur(1))
  SL=min(ul(2)-aal,ur(2)-aar)
  SR=max(ul(2)+aal,ur(2)+aar)

  enerl=(ul(3)-G1(ul(1)))/omg
	fl(1)=ul(1)*ul(2)
  fl(2)=ul(1)*ul(2)*ul(2)+ul(3)
  fl(3)=ul(2)*(0.5*ul(1)*ul(2)*ul(2)+enerl+ul(3))
	uul(1)=ul(1)
  uul(2)=ul(1)*ul(2)
  uul(3)=0.5*ul(1)*ul(2)*ul(2)+enerl

  enerr=(ur(3)-G1(ur(1)))/omg
	fr(1)=ur(1)*ur(2)
  fr(2)=ur(1)*ur(2)*ur(2)+ur(3)
  fr(3)=ur(2)*(0.5*ur(1)*ur(2)*ur(2)+enerr+ur(3))
  uur(1)=ur(1)
  uur(2)=ur(1)*ur(2)
  uur(3)=0.5*ur(1)*ur(2)*ur(2)+enerr

	if(0..le.SL)then
	  flux(1)=fl(1)
	  flux(2)=fl(2)
    flux(3)=fl(3)
	else if(SL.lt.0..and.0.lt.SR)then
	  flux(1)=(SR*fl(1)-SL*fr(1)-SR*SL*(uul(1)-uur(1)))/(SR-SL)
	  flux(2)=(SR*fl(2)-SL*fr(2)-SR*SL*(uul(2)-uur(2)))/(SR-SL)
    flux(3)=(SR*fl(3)-SL*fr(3)-SR*SL*(uul(3)-uur(3)))/(SR-SL)
	else
	  flux(1)=fr(1)
	  flux(2)=fr(2)
    flux(3)=fr(3)
	endif
end subroutine hll

subroutine hllc(ul,ur,flux)

  use riemann_parameter
  use functions
  implicit none

  real, intent(in) :: ul(3),ur(3)
  real, intent(out) :: flux(3)

  real :: fl(3),fr(3),uf(3),uur(3),uul(3), uur2(3),uul2(3)
  real :: aal, aar, SL, SR, SS, enerl, enerr, ad1, ad2

  aal=sqrt(f1(ul(1))+(omg+1.)*ul(3)/ul(1))
	aar=sqrt(f1(ur(1))+(omg+1.)*ur(3)/ur(1))
  SL=min(ul(2)-aal,ur(2)-aar)
  SR=max(ul(2)+aal,ur(2)+aar)
  SS=(ur(3)-ul(3)+ul(1)*ul(2)*(SL-ul(2))-ur(1)* &
  ur(2)*(SR-ur(2)))/(ul(1)*(SL-ul(2))-ur(1)*(SR-ur(2)))
  enerl=(ul(3)-G1(ul(1)))/omg
	fl(1)=ul(1)*ul(2)
  fl(2)=ul(1)*ul(2)*ul(2)+ul(3)
  fl(3)=ul(2)*(0.5*ul(1)*ul(2)*ul(2)+enerl+ul(3))

	uul(1)=ul(1)
  uul(2)=ul(1)*ul(2)
  uul(3)=0.5*ul(1)*ul(2)*ul(2)+enerl

  enerr=(ur(3)-G1(ur(1)))/omg
	fr(1)=ur(1)*ur(2)
  fr(2)=ur(1)*ur(2)*ur(2)+ur(3)
  fr(3)=ur(2)*(0.5*ur(1)*ur(2)*ur(2)+enerr+ur(3))

  uur(1)=ur(1)
  uur(2)=ur(1)*ur(2)
  uur(3)=0.5*ur(1)*ur(2)*ur(2)+enerr

  ad1=ul(1)*(SL-ul(2))/(SL-SS)
  ad2=ur(1)*(SR-ur(2))/(SR-SS)

  uul2(1)=ad1
  uul2(2)=ad1*SS
  uul2(3)=ad1*(uul(3)/uul(1)+(SS-ul(2))* &
  (SS+ul(3)/((uul(1)*(SL-ul(2))))))
  uur2(1)=ad2
  uur2(2)=ad2*SS
  uur2(3)=ad2*(uur(3)/uur(1)+(SS-ur(2))* &
  (SS+ur(3)/((uur(1)*(SR-ur(2))))))

	if(0..le.SL)then
	  flux(1)=fl(1)
	  flux(2)=fl(2)
    flux(3)=fl(3)
	else if(SL.lt.0..and.0.lt.SS)then
	  flux(1)=fl(1)+SL*(uul2(1)-uul(1))
	  flux(2)=fl(2)+SL*(uul2(2)-uul(2))
    flux(3)=fl(3)+SL*(uul2(3)-uul(3))
	else if(SS.le.0..and.0.le.SR)then
    flux(1)=fr(1)+SR*(uur2(1)-uur(1))
	  flux(2)=fr(2)+SR*(uur2(2)-uur(2))
    flux(3)=fr(3)+SR*(uur2(3)-uur(3))
	else
	  flux(1)=fr(1)
	  flux(2)=fr(2)
    flux(3)=fr(3)
	endif
end subroutine hllc

subroutine roe(ul,ur,flux)

  use riemann_parameter
  use functions
  implicit none

  real, intent(in) :: ul(3),ur(3)
  real, intent(out) :: flux(3)
  ! 局部变量
  real :: hl, hr, temp, epsilon, enerl, enerr
  real :: rho_sqrtl, rho_sqrtr, rho_avg, h_avg, a_avg, u_avg
  real :: delta(3), lambda(3), lambda_entropy_fix(3), roe_eigenvectors(3,3)
  real :: wave_strength(3), flux_l(3), flux_r(3)
  integer :: i, j

  ! 计算左侧和右侧的总焓
  hl = ((omg + 1.0) * ul(3))/(omg * ul(1)) + 0.5 * (ul(2)**2)
  hr = ((omg + 1.0) * ur(3))/(omg * ur(1)) + 0.5 * (ur(2)**2)

  ! 计算Roe平均值
  rho_sqrtl = sqrt(ul(1))
  rho_sqrtr = sqrt(ur(1))
  rho_avg = sqrt(ul(1)*ur(1))
  u_avg = (rho_sqrtl * ul(2) + rho_sqrtr * ur(2)) / (rho_sqrtl + rho_sqrtr)
  h_avg = (rho_sqrtl * hl + rho_sqrtr * hr) / (rho_sqrtl + rho_sqrtr)
  a_avg = sqrt(omg * (h_avg - 0.5 * u_avg**2))

  ! 计算特征值
  lambda(1) = u_avg - a_avg
  lambda(2) = u_avg
  lambda(3) = u_avg + a_avg

  !熵修正
  epsilon = (abs(u_avg) + abs(a_avg)) / 10.0
  do i = 1, 3
    if(abs(lambda(i)).le.epsilon) then
      lambda_entropy_fix(i) = (abs(lambda(i))**2 + epsilon**2)/(2 * epsilon)
    else
      lambda_entropy_fix(i) = abs(lambda(i))
    end if
  end do

  ! 计算通量
  enerl=(ul(3)-G1(ul(1)))/omg
	flux_l(1)=ul(1)*ul(2)
  flux_l(2)=ul(1)*ul(2)*ul(2)+ul(3)
  flux_l(3)=ul(2)*(0.5*ul(1)*ul(2)*ul(2)+enerl+ul(3))

  enerr=(ur(3)-G1(ur(1)))/omg
	flux_r(1)=ur(1)*ur(2)
  flux_r(2)=ur(1)*ur(2)*ur(2)+ur(3)
  flux_r(3)=ur(2)*(0.5*ur(1)*ur(2)*ur(2)+enerr+ur(3))

  ! 计算差分
  do i = 1, 3
    delta(i)  = ur(i) - ul(i)
  end do 

  ! 计算右特征向量矩阵
  roe_eigenvectors(1,1) = 1
  roe_eigenvectors(1,2) = u_avg - a_avg
  roe_eigenvectors(1,3) = h_avg - u_avg * a_avg
  roe_eigenvectors(2,1) = 1
  roe_eigenvectors(2,2) = u_avg
  roe_eigenvectors(2,3) = 0.5 * (u_avg**2)
  roe_eigenvectors(3,1) = 1
  roe_eigenvectors(3,2) = u_avg + a_avg
  roe_eigenvectors(3,3) = h_avg + u_avg * a_avg

  ! 计算波幅系数
  wave_strength(1) = (delta(3) - rho_avg * a_avg * delta(2))/(2 * a_avg**2)
  wave_strength(2) = delta(1) - (delta(3) / (a_avg**2))
  wave_strength(3) = (delta(3) + rho_avg * a_avg * delta(2))/(2 * a_avg**2)

  ! 得出数值通量
  do i = 1, 3
    temp = 0.0
    do j = 1, 3
      temp = temp + wave_strength(j) * lambda_entropy_fix(j) * roe_eigenvectors(j,i)
      ! temp = temp + wave_strength(j) * abs(lambda(j)) * roe_eigenvectors(j,i)
    end do
    flux(i) = 0.5*(flux_l(i)+flux_r(i)) - 0.5*temp
  end do
  
  ! ! debug
  ! if(abs(delta(2)) > 1.0) then
  !   write (*,*) "a_avg = ", a_avg
  !   write (*,*) "u_avg = ", u_avg
  !   write (*,*) "lambda = ", lambda
  !   write (*,*) "lambda_entropy_fix = ", lambda_entropy_fix
  !   write (*,*) "roe_eigenvectors = ", roe_eigenvectors
  !   write (*,*) "wave_strength = ", wave_strength
  !   write (*,*) "flux = ", flux
  !   write (*,*) "flux_l = ", flux_l
  !   write (*,*) "flux_r = ", flux_r
  ! end if
end subroutine roe
