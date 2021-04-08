module mod_mat_tab

implicit none

integer, parameter, public :: PR=selected_real_kind(8)
real(PR), parameter, public :: Pi=2.0_PR*Asin(1.0_PR)

integer :: outtabE=0
integer :: outtabr=0
integer :: outtabT=0
integer :: outtabP=0


type mat_tab
   character(len=2) :: nom
   integer :: Z
   real(PR) :: A
   real(PR) :: rhos, sig0, lamb0
   integer :: mode=1
   integer :: N_nu, N_T, N_sat, N_rho
   real(PR) :: rho0, rhof, T0, Tf, DT, nu0, nuf, Dlogr, logr0, logrf
   real(PR) :: Emin, Emax, Pmin, Pmax
   real(PR), allocatable :: nu(:), T(:), rho(:)   
   real(PR), allocatable :: P(:,:), E(:,:), Cv(:,:), Cp(:,:), g(:,:)
   real(PR), allocatable :: c(:,:), Zb(:,:), sig(:,:), lamb(:,:), Seeb(:,:), pinf(:,:)
   real(PR), allocatable :: Tm(:)
   real(PR), allocatable :: T_sat(:), P_sat(:)
   real(PR), allocatable :: nu_sat(:,:), rho_sat(:,:), E_sat(:,:), Cv_sat(:,:), Cp_sat(:,:), g_sat(:,:)
   real(PR), allocatable :: c_sat(:,:), Zb_sat(:,:), sig_sat(:,:), lamb_sat(:,:), Seeb_sat(:,:)
   real(PR), allocatable :: pinf_sat(:,:)
end type mat_tab


type(mat_tab) :: mymat


contains

subroutine read_table_NT(file_in,mat)

!N_nu,N_T,N_sat,rho0,rhof,T,nu,P,Es,Cv,Cp,gam,sound,Zbar,sigma,lambda,Seeb,P_sat,nu_sat,E_sat,Cv_sat,Cp_sat,&
!                      g_sat, c_sat, z_sat, sig_sat, lamb_sat, Seeb_sat)

   implicit none
   character(len=*) :: file_in
   type(mat_tab), intent(out) :: mat

   !integer, intent(out) :: N_nu, N_T, N_sat 
   !real(PR), intent(out) :: rho0, rhof
   !real(PR), allocatable, intent(out) :: T(:), nu(:)
   !real(PR), allocatable, intent(out) :: P_sat(:), nu_sat(:,:),E_sat(:,:),Cv_sat(:,:), Cp_sat(:,:), g_sat(:,:), c_sat(:,:), z_sat(:,:), sig_sat(:,:), lamb_sat(:,:), Seeb_sat(:,:)
   !real(PR), allocatable, intent(out) :: P(:,:), Es(:,:), Cv(:,:), Cp(:,:),&
   !gam(:,:), sound(:,:),Zbar(:,:),sigma(:,:),lambda(:,:),Seeb(:,:) 

   real(PR) :: P1, T1, E1, Cv1, Cp1, pinf, nu
   integer :: i,j,id, ierr 
   character(len=50) :: dump


   print*, ' START READ TABLE ... '
   id=12
   open(unit=id ,file=file_in,status='old',action='read',iostat=ierr)
   if(ierr.ne.0)then
     print*, " PROBLEM OPENING", file_in, " :( "
     stop
   else
     print*, " FILE", file_in, " OPENED :) "
   endif

   read(id,'(A)') dump

   read(id,'(A)') dump
   i=index(dump,'for') ; read(dump(i+4:),*) mat%nom
   print*, 'nom=', mat%nom

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%mode
   print*, 'mode=', mat%mode

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%Z
   print*, 'Z=', mat%Z

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%A
   print*, 'A=', mat%A

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%rhos
   print*, 'rhos=', mat%rhos

   if(mat%mode.eq.1)then
      read(id,'(A)') dump
      i=index(dump,'=') ; read(dump(i+1:),*) mat%sig0
      print*, 'sig0=', mat%sig0

      read(id,'(A)') dump
      i=index(dump,'=') ; read(dump(i+1:),*) mat%lamb0
      print*, 'lamb0=', mat%lamb0
   endif


   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%N_nu
   print*, 'N_nu=', mat%N_nu 
   mat%N_rho=mat%N_nu

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%rho0
   print*, 'rho0=', mat%rho0

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%rhof
   print*, 'rhof=', mat%rhof

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%N_T
   print*, 'N_T=', mat%N_T

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%T0
   print*, 'T0=', mat%T0
 
   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%Tf
   print*, 'Tf=', mat%Tf
 
   if(mat%mode.eq.1)then
     allocate(mat%T(1:mat%N_T))
     read(id,*) dump
     do i=1,mat%N_T
       read(id,'(A)') dump
       read(dump(3:),*) mat%T(i)
       !print*, dump 
     enddo
   endif

   read(id,'(A)') dump
   !print*, '1', dump

   read(id,'(A)') dump
   !print*, '2', dump

   read(id,'(A)') dump
   !print*, '3', dump

  allocate(mat%P(1:mat%N_nu,1:mat%N_T))    ; mat%P(:,:)=0.0_PR
  allocate(mat%E(1:mat%N_nu,1:mat%N_T))    ; mat%E(:,:)=0.0_PR
  allocate(mat%Cv(1:mat%N_nu,1:mat%N_T))   ; mat%Cv(:,:)=0.0_PR
  allocate(mat%g(1:mat%N_nu,1:mat%N_T))    ; mat%g(:,:)=0.0_PR
  allocate(mat%c(1:mat%N_nu,1:mat%N_T))    ; mat%c(:,:)=0.0_PR
  allocate(mat%Zb(1:mat%N_nu,1:mat%N_T))   ; mat%Zb(:,:)=0.0_PR

  if(mat%mode.eq.1)then
     allocate(mat%nu(1:mat%N_nu))             ; mat%nu(:)=0.0_PR
     allocate(mat%rho(1:mat%N_rho))           ; mat%rho(:)=0.0_PR
     allocate(mat%sig(1:mat%N_nu,1:mat%N_T))  ; mat%sig(:,:)=0.0_PR
     allocate(mat%lamb(1:mat%N_nu,1:mat%N_T)) ; mat%lamb(:,:)=0.0_PR
     allocate(mat%Seeb(1:mat%N_nu,1:mat%N_T)) ; mat%Seeb(:,:)=0.0_PR
     allocate(mat%Cp(1:mat%N_nu,1:mat%N_T))   ; mat%Cp(:,:)=0.0_PR
  else
     allocate(mat%pinf(1:mat%N_nu,1:mat%N_T))  ; mat%pinf(:,:)=0.0_PR
  endif



   do j=1,mat%N_T
     read(id,'(A)') dump
     read(id,'(A)') dump
     read(id,'(A)') dump

     if(mat%mode.eq.1)then
        do i=1,mat%N_nu
           read(id,*) mat%nu(i), P1, mat%P(i,j), E1, mat%E(i,j), mat%Cv(i,j), mat%Cp(i,j),&
                      mat%c(i,j), mat%g(i,j), mat%Zb(i,j),&
                      mat%sig(i,j), mat%lamb(i,j), mat%Seeb(i,j), pinf
        enddo
     else
        do i=1,mat%N_nu
           read(id,*) mat%P(i,j), mat%E(i,j), mat%Cv(i,j),&
                      mat%c(i,j), mat%g(i,j), mat%Zb(i,j),&
                      mat%pinf(i,j) 
        enddo
     endif
     read(id,'(A)') dump
     !print*, '7', dump 
     read(id,'(A)') dump
     !print*, '8', dump 
   enddo

   !!!---lecture de la courbe de saturation
   read(id,'(A)') dump
   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat%N_sat
 
   allocate(mat%nu_sat(1:2,1:mat%N_sat))   ; mat%nu_sat(:,:)=0.0_PR
   allocate(mat%rho_sat(1:2,1:mat%N_sat))   ; mat%rho_sat(:,:)=0.0_PR
   allocate(mat%T_sat(1:mat%N_sat))        ; mat%T_sat(:)=0.0_PR
   allocate(mat%P_sat(1:mat%N_sat))        ; mat%P_sat(:)=0.0_PR
   allocate(mat%E_sat(1:2,1:mat%N_sat))    ; mat%E_sat(:,:)=0.0_PR
   allocate(mat%Cv_sat(1:2,1:mat%N_sat))   ; mat%Cv_sat(:,:)=0.0_PR
   allocate(mat%c_sat(1:2,1:mat%N_sat))    ; mat%c_sat(:,:)=0.0_PR
   allocate(mat%g_sat(1:2,1:mat%N_sat))    ; mat%g_sat(:,:)=0.0_PR
   allocate(mat%Zb_sat(1:2,1:mat%N_sat))   ; mat%Zb_sat(:,:)=0.0_PR

   if(mat%mode.eq.1)then
      allocate(mat%Cp_sat(1:2,1:mat%N_sat))   ; mat%Cp_sat(:,:)=0.0_PR
      allocate(mat%sig_sat(1:2,1:mat%N_sat))  ; mat%sig_sat(:,:)=0.0_PR
      allocate(mat%lamb_sat(1:2,1:mat%N_sat)) ; mat%lamb_sat(:,:)=0.0_PR
      allocate(mat%Seeb_sat(1:2,1:mat%N_sat)) ; mat%Seeb_sat(:,:)=0.0_PR
   else
      allocate(mat%pinf_sat(1:2,1:mat%N_sat)) ; mat%pinf_sat(:,:)=0.0_PR
   endif

   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de bulle

   if(mat%mode.eq.1)then
      do i=1, mat%N_sat
        read(id,*) mat%T_sat(i), mat%P_sat(i), mat%nu_sat(1,i),&   !, mat%nu_sat(2,i),&
                                               mat%E_sat(1,i),&    !, mat%E_sat(2,i),&
                                               mat%Cv_sat(1,i),&   !, mat%Cv_sat(2,i),&
                                               mat%Cp_sat(1,i),&   !, mat%Cp_sat(2,i),&
                                               mat%c_sat(1,i),&    !, mat%c_sat(2,i),&
                                               mat%g_sat(1,i),&    !, mat%g_sat(2,i),&
                                               mat%Zb_sat(1,i),&   !, mat%Zb_sat(2,i),&
                                               mat%sig_sat(1,i),&  !, mat%sig_sat(2,i),&
                                               mat%lamb_sat(1,i),& !, mat%lamb_sat(2,i),&
                                               mat%Seeb_sat(1,i),& !, mat%Seeb_sat(2,i),&
                                               pinf            
      enddo
   else
      do i=1, mat%N_sat
        read(id,*) mat%T_sat(i), mat%P_sat(i), mat%nu_sat(1,i),&   !, mat%nu_sat(2,i),&
                                               mat%E_sat(1,i),&    !, mat%E_sat(2,i),&
                                               mat%Cv_sat(1,i),&   !, mat%Cv_sat(2,i),&
                                               mat%c_sat(1,i),&    !, mat%c_sat(2,i),&
                                               mat%g_sat(1,i),&    !, mat%g_sat(2,i),&
                                               mat%Zb_sat(1,i),&   !, mat%Zb_sat(2,i),&
                                               mat%pinf_sat(1,i)            
      enddo
   endif


   read(id,'(A)') dump !!! il ya N_sat+1 lignes
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de goute

   if(mat%mode.eq.1)then
      do i=1, mat%N_sat
        read(id,*) mat%T_sat(i), mat%P_sat(i), mat%nu_sat(2,i),&   !, mat%nu_sat(2,i),&
                                               mat%E_sat(2,i),&    !, mat%E_sat(2,i),&
                                               mat%Cv_sat(2,i),&   !, mat%Cv_sat(2,i),&
                                               mat%Cp_sat(2,i),&   !, mat%Cp_sat(2,i),&
                                               mat%c_sat(2,i),&    !, mat%c_sat(2,i),&
                                               mat%g_sat(2,i),&    !, mat%g_sat(2,i),&
                                               mat%Zb_sat(2,i),&   !, mat%Zb_sat(2,i),&
                                               mat%sig_sat(2,i),&  !, mat%sig_sat(2,i),&
                                               mat%lamb_sat(2,i),& !, mat%lamb_sat(2,i),&
                                               mat%Seeb_sat(2,i),& !, mat%Seeb_sat(2,i),&
                                               pinf            
      enddo
   else
      do i=1, mat%N_sat
        read(id,*) mat%T_sat(i), mat%P_sat(i), mat%nu_sat(2,i),&   !, mat%nu_sat(2,i),&
                                               mat%E_sat(2,i),&    !, mat%E_sat(2,i),&
                                               mat%Cv_sat(2,i),&   !, mat%Cv_sat(2,i),&
                                               mat%c_sat(2,i),&    !, mat%c_sat(2,i),&
                                               mat%g_sat(2,i),&    !, mat%g_sat(2,i),&
                                               mat%Zb_sat(2,i),&   !, mat%Zb_sat(2,i),&
                                               mat%pinf_sat(2,i)            
      enddo
   endif

   read(id,'(A)') dump !!! il ya N_sat+1 lignes
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de fusion

   allocate(mat%Tm(1:mat%N_nu)) ; mat%Tm(:)=0.0_PR
   do i=1, mat%N_nu
     read(id,*) nu, mat%Tm(i)
   enddo


   close(id)

   print*, ' READ FINI ! :) '

   mat%rho_sat(:,:)=1.0_PR/mat%nu_sat(:,:)

   if(mat%mode.eq.1)then
      mat%rho(:)=1.0_PR/mat%nu(:)
      open(unit=1,file='density_tab.dat',status='replace')
      do i=1,mat%N_rho
         write(1,*) mat%rho(i)
      enddo
      close(1)
   endif



    mat%nu0=1.0_PR/mat%rhof ; mat%nuf=1.0_PR/mat%rho0
    !mat%T0=mat%T(1) ; mat%Tf=mat%T(mat%N_T)
    mat%DT=(mat%Tf-mat%T0)/real(mat%N_T-1,PR) !!! N_T value=N_T-1 intervalles
    mat%logr0=log(mat%rho0) 
    mat%logrf=log(mat%rhof) 
    !mat%logr0=log(mat%rho(1))
    !mat%logrf=log(mat%rho(mat%N_rho))
    mat%Dlogr=(mat%logrf-mat%logr0)/real(mat%N_rho-1,PR)
    mat%Emin=minval(mat%E(:,:))
    mat%Emax=maxval(mat%E(:,:))
    mat%Pmin=minval(mat%P(:,:))
    mat%Pmax=maxval(mat%P(:,:))


    if(mat%mode.eq.2)then
     allocate(mat%T(1:mat%N_T))
     do i=1,mat%N_T
       mat%T(i)=mat%T0+real(i-1,PR)*mat%DT
     enddo
     allocate(mat%rho(1:mat%N_rho))
     do i=1,mat%N_rho
       mat%rho(i)=mat%rho0*(mat%rhof/mat%rho0)**(real(i-1,PR)/real(mat%N_rho-1,PR))
     enddo
    
    endif

end subroutine read_table_NT

subroutine get_index_rT(mat,rho,T,ir,iT)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho, T
   integer, intent(out) :: ir, iT
   real(PR) :: x,y
   !!! floor: integer(4) < 2e9

   call get_ir(mat,rho,ir)

   call get_iT_T(mat,T,iT)


end subroutine get_index_rT


subroutine get_index_rE(mat,rho,Es,ir,iT)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho, Es
   integer, intent(out) :: ir, iT
   real(PR) :: x, x1, x2, cr1, cr2, E1, E2
   integer :: i, iTmin, iTmax
   logical :: found 

   call get_ir(mat,rho,ir)

   if(Es.gt.mat%Emax)then

     iT=mat%N_T
     outtabE=1

   elseif(Es.lt.mat%Emin)then

     iT=1
     outtabE=-1

   else

     x=log(rho)
     x1=mat%logr0+real(ir-1,PR)*mat%Dlogr 
     !x1=log(mat%rho(ir)) 
     x2=x1+mat%Dlogr   
     x=max(x,x1) ; x=min(x,x2)
  
     cr1=(x2-x)/(x2-x1)
     cr2=(x-x1)/(x2-x1)
  
     iTmin=1
     iTmax=mat%N_T
     found=.false.
     i=0
  
     do while(.not.found)
  
        iT=(iTmin+iTmax)/2
  
        E1=cr1*mat%E(ir,iT)+cr2*mat%E(ir+1,iT)
        E2=cr1*mat%E(ir,iT+1)+cr2*mat%E(ir+1,iT+1)
  
        if(E1.gt.Es)then
          iTmax=iT 
        elseif(E2.lt.Es)then
          iTmin=iT
        else
          found=.true.
        endif
  
        i=i+1
  
        if(i.gt.990)then
            print*, iTmin, iTmax, E1, E2, Es
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM T_from_E'
            stop
        endif
  
     enddo
  
     print*, 'get_index_rE it=', i

   endif

end subroutine get_index_rE

subroutine get_index_rP(mat,rho,P,ir,iT)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho, P
   integer, intent(out) :: ir, iT
   real(PR) :: x, x1, x2, cr1, cr2, P1, P2
   integer :: i, iTmin, iTmax
   logical :: found 

   call get_ir(mat,rho,ir)

   if(P.gt.mat%Pmax)then

     iT=mat%N_T
     outtabP=1

   elseif(P.lt.mat%Pmin)then

     iT=1
     outtabP=-1

   else

     x=log(rho)
     x1=mat%logr0+real(ir-1,PR)*mat%Dlogr 
     !x1=log(mat%rho(ir)) 
     x2=x1+mat%Dlogr   
     x=max(x,x1) ; x=min(x,x2)
  
     cr1=(x2-x)/(x2-x1)
     cr2=(x-x1)/(x2-x1)
  
     iTmin=1
     iTmax=mat%N_T
     found=.false.
     i=0

     do while(.not.found)
  
        iT=(iTmin+iTmax)/2
  
        P1=cr1*mat%P(ir,iT)+cr2*mat%P(ir+1,iT)
        P2=cr1*mat%P(ir,iT+1)+cr2*mat%P(ir+1,iT+1)
  
        if(P1.gt.P)then
          iTmax=iT 
        elseif(P2.lt.P)then
          iTmin=iT
        else
          found=.true.
        endif
  
        i=i+1
  
        if(i.gt.990)then
            print*, iTmin, iTmax, P1, P2, P
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM T_from_P'
            stop
        endif
  
     enddo
  
     print*, 'get_index_rP it=', i

   endif

end subroutine get_index_rP




subroutine get_ir(mat,rho,ir)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho
   integer, intent(out) :: ir
   real(PR) :: x
   !!! floor: integer(4) < 2e9

   if(rho.gt.mat%rhof)then
    ir=mat%N_rho
    outtabr=1
   elseif(rho.lt.mat%rho0)then
    ir=1
    outtabr=-1
   else  
    x=log(rho) 
    ir=floor(min((x-mat%logr0)/mat%Dlogr,1.0e9_PR))+1
    ir=max(ir,1) ; ir=min(ir,mat%N_rho)
   endif

end subroutine get_ir

subroutine get_iT_T(mat,T,iT)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: T
   integer, intent(out) :: iT
   real(PR) :: y
   !!! floor: integer(4) < 2e9
   if(T.gt.mat%Tf)then
    iT=mat%N_T
    outtabT=1
   elseif(T.lt.mat%T0)then
    iT=1
    outtabT=-1
   else  
    y=T 
    iT=floor(min((y-mat%T0)/mat%DT,1.0e9_PR))+1
    iT=max(iT,1) ; iT=min(iT,mat%N_T)
   endif

end subroutine get_iT_T


subroutine get_from_rT_lin(mat,rho,T,ir,iT,field,val)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho, T
   integer, intent(in) :: ir, iT
   real(PR) , intent(in) :: field(1:mat%N_rho,1:mat%N_T)
   real(PR), intent(out) :: val
   real(PR) :: x,y,x1,x2,y1,y2,zNE,zNW,zSE,zSW
   integer :: ir1,ir2,iT1,iT2

   ir1=ir
   ir2=min(ir+1,mat%N_rho)
   iT1=iT
   iT2=min(iT+1,mat%N_T)

   x=log(rho) 
   x1=mat%logr0+real(ir1-1,PR)*mat%Dlogr   !!!log(mat%rho(ir1)) 
   x2=x1+mat%Dlogr   
   x=max(x,x1) ; x=min(x,x2)

   y=T 
   y1=mat%T0+real(iT1-1,PR)*mat%DT   !!!!mat%T(iT1)
   y2=y1+mat%DT
   y=max(y,y1) ; y=min(y,y2)

   zNE=field(ir2,iT2)
   zNW=field(ir1,iT2)
   zSE=field(ir2,iT1)
   zSW=field(ir1,iT1)

   val=interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y)

end subroutine get_from_rT_lin

subroutine get_from_rT_log(mat,rho,T,ir,iT,field,val)

   implicit none
   type(mat_tab), intent(in) :: mat
   real(PR), intent(in) :: rho, T
   integer, intent(in) :: ir, iT
   real(PR) , intent(in) :: field(1:mat%N_rho,1:mat%N_T)
   real(PR), intent(out) :: val
   real(PR) :: x,y,x1,x2,y1,y2,zNE,zNW,zSE,zSW
   integer :: ir1,ir2,iT1,iT2

   ir1=ir
   ir2=min(ir+1,mat%N_rho)
   iT1=iT
   iT2=min(iT+1,mat%N_T)

   x=log(rho) ; x=max(x,mat%logr0) ; x=min(x,mat%logrf)
   x1=mat%logr0+real(ir1-1,PR)*mat%Dlogr   !!!log(mat%rho(ir1)) 
   x2=x1+mat%Dlogr   

   y=T ; y=max(y,mat%T0) ; y=min(y,mat%Tf)
   y1=mat%T0+real(iT1-1,PR)*mat%DT  !!!mat%T(iT1)
   y2=y1+mat%DT

   zNE=log(field(ir2,iT2))
   zNW=log(field(ir1,iT2))
   zSE=log(field(ir2,iT1))
   zSW=log(field(ir1,iT1))

   val=exp(interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y))

end subroutine get_from_rT_log


real(PR) function interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y)

   implicit none
   real(PR), intent(in) :: x1,x2,y1,y2
   real(PR), intent(in) :: zNE,zNW,zSE,zSW
   real(PR), intent(in) :: x,y
   real(PR) :: invs, sNE, sNW, sSE, sSW

   invs=1.0_PR/((x2-x1)*(y2-y1))

   sNE=(y-y1)*(x-x1)
   sNW=(y-y1)*(x2-x)
   sSE=(y2-y)*(x-x1)
   sSW=(y2-y)*(x2-x)
   
   interplin2D=invs*(sNE*zNE+sNW*zNW+sSE*zSE+sSW*zSW)

end function interplin2D



subroutine check_tab

   implicit none

   print*, '  check tab...'

   if(outtabr.eq.1)then

     print*, '   UPPER BOUND OF TAB IN DENSITY !!!  :('

   elseif(outtabr.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN DENSITY !!!  :('

   endif

   if(outtabT.eq.1)then

     print*, '   UPPER BOUND OF TAB IN TEMPERATURE !!!  :('

   elseif(outtabT.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN TEMPERATURE !!!  :('

   endif

   if(outtabE.eq.1)then

     print*, '   UPPER BOUND OF TAB IN ENERGY !!!  :('

   elseif(outtabE.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN ENERGY !!!  :('

   endif

   if(outtabr.eq.0.and.outtabT.eq.0.and.outtabE.eq.0)then

     print*, '  -> everything is OK :) !!!'      

   endif

   if(outtabP.eq.1)then

     print*, '   UPPER BOUND OF TAB IN PRESSURE !!!  :('

   elseif(outtabP.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN PRESSURE !!!  :('

   endif



end subroutine check_tab




end module mod_mat_tab
