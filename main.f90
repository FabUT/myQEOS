program main

use mod_data 

use mod_mat_tab, only : mat_tab, mymat, read_table_NT, get_index_rT,&
                        get_from_rT_lin, get_from_rT_log,&
                        get_index_rE, check_tab,get_index_rP


!use mod_snes, only : init_snes, restart_snes, solve_snes !!!, F_test, J_test

use mod_eos, only : init_eos, P_tot, E_tot, rhos, P_i, P_e, P_b, Pcrit, Tcrit, rcrit,&
                    element_name, ddTr_E_tot, ddTp_E_tot, ddTr_P_tot,ddrT_P_tot,&
                    ddre_P,&
                    Avogadro, eV, A, Rg, solid, transport_properties,&
                    A_alpha, A_alpha2, compute_ZFS, sig0, lamb0, Tm,&
                    n, m, mfromn, cohesive_energy,&
                    bounding, electrons, ions

use thomasfermi, only : init_thomas_fermi_ionization_data, thomas_fermi_zbar

use mod_ensight, only : output_ensight

use mod_Fermi_Dirac, only : fd1h

implicit none

!integer :: N
real(PR), allocatable :: x(:), b(:)
integer :: i,j,k
real(PR) :: P, DTemp, f(1:1), dfdx(1:1)
!integer :: Z
real(PR) :: nu,nu1,nu2,Dnu,rho0,rhof,rho,T0,Tf,Tc,Pc,rhoc
real(PR) :: pmin,pmax,rhomin,rhomax,r1,r2,r3
integer :: foundP
integer :: id, N_nu, N_T, Nr, it, N_T_sat
real(PR) :: T, Te,  Et, Cvt, Cpt,ct,gt,zt,zfs,zs,sig,lamb,Sth,sigt
real(PR) :: sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu
real(PR), allocatable :: temperatures(:), Vs(:)
real(PR), allocatable  :: Pression(:,:), Es(:,:), Cv(:,:), Cp(:,:), &
                          sound(:,:),gam(:,:),w(:,:), Zbar(:,:),&
                          sigma(:,:), lambda(:,:), Seebeck(:,:),&
                          Pression_e(:,:),Pression_i(:,:),Pression_b(:,:)
real(PR) :: Ebi, Cvbi, Cpbi, cbi, gbi, zbi, sigi, Rs, pinf, pinfbi
real(PR), allocatable :: P_sat(:), E_sat(:,:), Cv_sat(:,:), Cp_sat(:,:)
real(PR), allocatable :: Pe_sat(:,:), Pi_sat(:,:), Pb_sat(:,:)
real(PR), allocatable :: c_sat(:,:), g_sat(:,:), z_sat(:,:)
real(PR), allocatable :: sig_sat(:,:),lamb_sat(:,:),Sth_sat(:,:),pinf_sat(:,:)
real(PR), allocatable :: nu_sat(:,:)
logical, allocatable :: found(:)

integer :: mode=1 !!!2 !!! 1: creation: 2: lecture

logical :: varynm=.false.
real(PR) :: facB0, facB0min, facB0max, nmin, nmax
integer :: kn, km, knm
integer :: N_n, N_m

!!!--------------TEST MYEOS------------------------------

!!!!---initialization des matériaux

IF(MODE.EQ.1)THEN

CALL INIT_EOS

Z=13  !!!26  !!29 !!!!!29 !! 13
!mat%Z=Z

IF(Z.eq.13)THEN
   !!!================================================================
   !!!                       ALUMINIUM
   !!!================================================================
   write(iout,*) ''
   write(iout,*) '-------------ALUMINIUM------------------'
   write(iout,*) ''
   Pc=Pcrit(Z) ; Tc=Tcrit(Z) ; rhoc=rcrit(Z)
   write(iout,'(A3,ES14.7,A4,ES14.7,A6,ES14.7)') 'Pc=', Pc, ' Tc=', Tc, ' rhoc=', rhoc
  
   bounding=.true. ; electrons=.true. ; ions=.true.

   N_nu=1000 ; rho0=1.0e-5_PR*rhos(Z) ; rhof=1.0_PR*rhos(Z)
 
   !N_T=400 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->10000 

   N_T=100 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 100 K 200->10000 

   do it=1,N_T
       temperatures(it)=T0+real(it-1,PR)*(Tf-T0)/real(N_T,PR)
   enddo
 
   !N_T=36 ; allocate(temperatures(1:N_T))
   !temperatures(1:36)=(/300.0_PR, 1000.0_PR, 1500.0_PR, 2000.0_PR, 2500.0_PR, 3000.0_PR,&
   !                     3500.0_PR, 4000.0_PR, 4250.0_PR, 4500.0_PR, 4625.0_PR, 4750.0_PR,&
   !                     4800.0_PR, 4850.0_PR, 4900.0_PR, 4950.0_PR, 5000.0_PR, 5050.0_PR,&
   !                     5100.0_PR, 5150.0_PR, 5200.0_PR, 5250.0_PR, 5300.0_PR, 5325.0_PR,&
   !                     5342.0_PR, 5375.0_PR, 5400.0_PR, 5425.0_PR, 5450.0_PR, 5475.0_PR,&
   !                     5500.0_PR, 6000.0_PR, 7000.0_PR, 8000.0_PR, 9000.0_PR, 10000.0_PR/)

   id=10
   call create_PE_VT_file(id,Z,N_nu,rho0,rhof,N_T,temperatures) 
 

ELSEIF(Z.eq.26)THEN
   !!!================================================================
   !!!                       FER
   !!!================================================================
   write(iout,*) ''
   write(iout,*) '----------------FER---------------------'
   write(iout,*) ''
   Pc=Pcrit(Z) ; Tc=Tcrit(Z) ; rhoc=rcrit(Z)
   Pc=10000.0_PR ; Tc=12000.0_PR ; rhoc=rcrit(Z)

   write(iout,'(A3,ES14.7,A4,ES14.7,A6,ES14.7)') 'Pc=', Pc, ' Tc=', Tc, ' rhoc=', rhoc

   N_n=1 ; N_m=1
   varynm=.false.


   !!!---------Sorties pour tables de calcul
   !N_T=401 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->10000 
   !do it=1,N_T
   !    temperatures(it)=T0+real(it-1,PR)*(Tf-T0)/real(N_T-1,PR)
   !enddo
   !!!------------------------------------------

   !!!---------Sorties pour gnuplot
   !N_T=36 ; allocate(temperatures(1:N_T))
   ! temperatures(1:36)=(/300.0_PR, 1000.0_PR, 2000.0_PR, 3000.0_PR, 3500.0_PR, 4000.0_PR, 4250.0_PR,&
   !                     4500.0_PR, 4750.0_PR, 5000.0_PR, 5100.0_PR, 5150.0_PR, 5200.0_PR, 5250.0_PR,&
   !                     5300.0_PR, 5350.0_PR, 5400.0_PR, 5450.0_PR, 5500.0_PR, 5525.0_PR, 5550.0_PR,&
   !                     5575.0_PR, 5605.0_PR, 5625.0_PR, 5650.0_PR, 5700.0_PR, 5800.0_PR, 5900.0_PR,&
   !                     6000.0_PR, 6250.0_PR, 6500.0_PR, 6750.0_PR, 7000.0_PR, 8000.0_PR, 9000.0_PR, 10000.0_PR/)
   !!!--------------------------------------

    !cohesive_energy(26)=cohesive_energy(26)*1.4_PR
    !n(26)=1.33_PR

    !N_T=100 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->10000 
    !do it=1,N_T
    !    temperatures(it)=T0+real(it-1,PR)*(Tf-T0)/real(N_T,PR)
    !enddo
 
    !N_T=36 ; allocate(temperatures(1:N_T))
    !temperatures(1:36)=(/300.0_PR, 1000.0_PR, 2000.0_PR, 3000.0_PR, 3500.0_PR, 4000.0_PR, 4250.0_PR,&
    !                    4500.0_PR, 4750.0_PR, 5000.0_PR,5250.0_PR,&
    !                    5500.0_PR,5750.0_PR, 6000.0_PR, 6100.0_PR, 6200.0_PR, 6300.0_PR, 6400.0_PR,&
    !                    6500.0_PR, 6600.0_PR, 6700.0_PR, 6800.0_PR, 6900.0_PR,&
    !                    7000.0_PR,7100.0_PR,7200.0_PR,7300.0_PR, 7400.0_PR, 7500.0_PR, 7600.0_PR,&
    !                    7700.0_PR, 7800.0_PR, 7900.0_PR, 8000.0_PR,& 
    !                    9000.0_PR, 10000.0_PR/)
  

    !temperatures(1:36)=(/300.0_PR, 1000.0_PR, 2000.0_PR, 3000.0_PR, 3500.0_PR, 4000.0_PR, 4500.0_PR,&
    !                    5000.0_PR, 5500.0_PR, 6000.0_PR, 6500.0_PR,&
    !                    7000.0_PR, 7500.0_PR, 8000.0_PR, 8500.0_PR, 9000.0_PR, 9250.0_PR, 9500.0_PR,&
    !                    9750.0_PR,  10000.0_PR, 10250.0_PR, 10500.0_PR, 10750.0_PR,&
    !                    11000.0_PR, 11250.0_PR, 11500.0_PR, 11750.0_PR, 12000.0_PR, 12100.0_PR, 12250.0_PR, 12500.0_PR,&
    !                    13000.0_PR, 13500.0_PR, 14000.0_PR,& 
    !                    14500.0_PR, 15000.0_PR/)
 

   !!!--------find Lennard Jones potential
   !mode_output=1
   !varynm=.false.
   N_T=1 ; allocate(temperatures(1:N_T))
   temperatures(1:1)=(/9375.0_PR/)
   !N_n=20 ; nmin=0.2_PR ; nmax=0.3_PR 
   !N_m=20 ; facB0min=0.5_PR ; facB0max=2.0_PR 

   
   !!!-------------------------------------

   N_nu=1000 ; rho0=1.0e-5_PR*rhos(Z) ; rhof=1.0_PR*rhos(Z)
  
   id=10
   call create_PE_VT_file(id,Z,N_nu,rho0,rhof,N_T,temperatures) 
 
ELSEIF(Z.eq.29)THEN

   !!!================================================================
   !!!                       CUIVRE
   !!!================================================================
   write(iout,*) ''
   write(iout,*) '-------------CUIVRE------------------'
   write(iout,*) ''
   Pc=Pcrit(Z) ; Tc=Tcrit(Z) ; rhoc=rcrit(Z)
   write(iout,'(A3,ES14.7,A4,ES14.7,A6,ES14.7)') 'Pc=', Pc, ' Tc=', Tc, ' rhoc=', rhoc

   mode_output=2
   N_T=401 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->10000 
   !!!N_T=1200 ; T0=200.0_PR ; Tf=30200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->30200
   do it=1,N_T
       temperatures(it)=T0+real(it-1,PR)*(Tf-T0)/real(N_T-1,PR)
   enddo
   
   !mode_output=1
   !N_T=36 ; allocate(temperatures(1:N_T))
   !temperatures(1:36)=(/300.0_PR, 1000.0_PR, 2000.0_PR, 3000.0_PR, 3500.0_PR, 4000.0_PR, 4500.0_PR,&
   !                     4750.0_PR, 5000.0_PR, 5100.0_PR, 5200.0_PR, 5300.0_PR, 5400.0_PR, 5500.0_PR,&
   !                     5600.0_PR, 5700.0_PR, 5800.0_PR, 5900.0_PR, 6000.0_PR, 6100.0_PR, 6200.0_PR,&
   !                     6300.0_PR, 6400.0_PR, 6500.0_PR, 6525.0_PR, 6550.0_PR, 6575.0_PR, 6600.0_PR,&
   !                     6625.0_PR, 6660.0_PR, 6700.0_PR, 6800.0_PR, 6900.0_PR, 7000.0_PR, 8000.0_PR, 10000.0_PR/)


   N_nu=1000 ; rho0=1.0e-5_PR*rhos(Z) ; rhof=1.5_PR*rhos(Z)
   id=10
   call create_PE_VT_file(id,Z,N_nu,rho0,rhof,N_T,temperatures) 


ENDIF !!! Z

!!!-----allocation des tableaux de sortie-----------------------------------

  allocate(Vs(1:N_nu)) ; Vs(:)=0.0_PR
  allocate(Pression(1:N_nu,1:N_T)) ; Pression(:,:)=0.0_PR
  allocate(Es(1:N_nu,1:N_T)) ; Es(:,:)=0.0_PR
  allocate(Cv(1:N_nu,1:N_T)) ; Cv(:,:)=0.0_PR
  allocate(Cp(1:N_nu,1:N_T)) ; Cp(:,:)=0.0_PR
  allocate(sound(1:N_nu,1:N_T)) ; sound(:,:)=0.0_PR
  allocate(gam(1:N_nu,1:N_T)) ; gam(:,:)=0.0_PR
  allocate(w(1:N_nu,1:N_T)) ; w(:,:)=0.0_PR
  allocate(Zbar(1:N_nu,1:N_T)) ; Zbar(:,:)=0.0_PR
  allocate(sigma(1:N_nu,1:N_T)) ; sigma(:,:)=0.0_PR
  allocate(lambda(1:N_nu,1:N_T)) ; lambda(:,:)=0.0_PR
  allocate(Seebeck(1:N_nu,1:N_T)) ; Seebeck(:,:)=0.0_PR
  allocate(Pression_e(1:N_nu,1:N_T)) ; Pression_e(:,:)=0.0_PR
  allocate(Pression_i(1:N_nu,1:N_T)) ; Pression_i(:,:)=0.0_PR
  allocate(Pression_b(1:N_nu,1:N_T)) ; Pression_b(:,:)=0.0_PR

!!!------------------------------------------------------------------------


write(iout,*) ''
write(iout,*) '==========================================='
write(iout,*) '       REMOVE VANDERWAALS LOOP !!! '
write(iout,*) '==========================================='
write(iout,*) ''

allocate(nu_sat(1:2,1:N_T)) ; allocate(P_sat(1:N_T))  ; allocate(found(1:N_T))
allocate(E_sat(1:2,1:N_T))  ; allocate(Cv_sat(1:2,1:N_T)) ; allocate(Cp_sat(1:2,1:N_T)) 
allocate(c_sat(1:2,1:N_T))  ; allocate(g_sat(1:2,1:N_T)) ; allocate(z_sat(1:2,1:N_T))
allocate(sig_sat(1:2,1:N_T)) ; allocate(lamb_sat(1:2,1:N_T)) ; allocate(Sth_sat(1:2,1:N_T)) 
allocate(pinf_sat(1:2,1:N_T))
allocate(Pe_sat(1:2,1:N_T)) ; allocate(Pi_sat(1:2,1:N_T)) ; allocate(Pb_sat(1:2,1:N_T))
nu_sat=0.0_PR ; P_sat=0.0_PR ; found=.false. ; E_sat=0.0_PR
Cv_sat=0.0_PR ; Cp_sat=0.0_PR ; c_sat=0.0_PR ; g_sat=0.0_PR ; z_sat=0.0_PR
sig_sat=0.0_PR ; lamb_sat=0.0_PR ; Sth_sat=0.0_PR ; pinf_sat=0.0_PR
Pe_sat=0.0_PR  ; Pi_sat=0.0_PR ; Pb_sat=0.0_PR


!!!------------------------- BOUCLE EN N ------------------------
!!! A désactivé en utilisation standard !!! :O

!!! knm=0
!!! 
!!! DO kn=1,N_n
!!! 
!!!   DO km=1,N_m
!!! 
!!!        knm=knm+1
!!! 
!!!         if(varynm)then
!!! 
!!!             n(Z)=nmin+real(kn-1,PR)*(nmax-nmin)/real(N_n-1,PR)
!!! 
!!!             facB0=facB0min+real(km-1,PR)*(facB0max-facB0min)/real(N_m-1,PR)
!!! 
!!!             m(Z)=mfromn(Z,n(Z),facB0)
!!! 
!!!         endif
!!! 
!!!         write(iout,*) ''
!!!         write(iout,*) '*********************************************'
!!!         !write(iout,*) ''
!!!         write(iout,*)  'RESOLUTION AVEC n=', n(Z), 'ET m=', m(Z)
!!!         write(iout,*) 'knm=', knm, 'kn=', kn, 'km=', km
!!!         write(iout,*) '*********************************************'
!!!         write(iout,*) ''

!!!------------------------- BOUCLE EN TEMPERATURE ------------------------

DO j=1,N_T

    Temp=temperatures(j)
   
    write(id,'(A)') '#----------------------------------------------------------------------------'
    if(mode_output.eq.1)then    
       write(id,'(A3,I3.3,A7,ES14.7,A50)') '# [',j,'] Temp=', temperatures(j), ' K: nu,P1,P2,Es1,Es2,Cv,Cp,c,g,Z,sig,lamb,Seeb,pinf'
    else
       write(id,'(A3,I3.3,A7,ES14.7,A21)') '# [',j,'] Temp=', temperatures(j), ' K: P,E,Cv,c,g,Z,pinf'
    endif 
    write(id,'(A)') '#----------------------------------------------------------------------------' 

    if(Temp.lt.0.8_PR*Tc)then
      rhomax=1.5_PR*rhos(Z) ; rhomin=1.0e-5_PR*rhos(Z) ; Nr=100000
    else
      rhomax=1.0_PR*rhos(Z) ; rhomin=1.0e-2_PR*rhos(Z) ; Nr=100000
    endif
   
    write(iout,'(A4,I3,A1,I3,A3,ES14.7,A2)')  ' it ', j, '/', N_T, ' T=', Temp, ' K'
   
    call MAXWELL(Pc,rhoc,rhomin,rhomax,Nr,foundP,it,P,r1,r2,r3)
   
    if(foundP.eq.1)then
      write(iout,"(A11,I4,A3,ES14.7,A6, 3(ES14.7,' '))") '     > Nit=', it, ' P=', P/pc, ' r123=', r1,r2,r3 
    elseif(foundP.eq.2)then
      write(iout,"(A11,I4,A20,ES14.7,A6,3(ES14.7,' '))") '     > Nit=', it, ' limite atteinte: P=', P, ' r123=', r1,r2,r3 
    else
      write(iout,'(A11,I4,ES14.7,A15)') '     > Nit=', it, Temp, ' P not found :('    
    endif
   
    if(foundP.eq.1.or.foundP.eq.2)then
       found(j)=.true.
    else
       found(j)=.false.
    endif
   
    !!!---calcul des grandeurs de saturation (binodales)
    if(found(j))then
       nu_sat(1,j)=1.0_PR/r1
       nu_sat(2,j)=1.0_PR/r3
       P_sat(j)=P
       Pe_sat(1,j)=P_e(Z,r1,Temp)
       Pi_sat(1,j)=P_i(Z,r1,Temp)
       Pb_sat(1,j)=P_b(Z,r1)
       Pe_sat(2,j)=P_e(Z,r3,Temp)
       Pi_sat(2,j)=P_i(Z,r3,Temp)
       Pb_sat(2,j)=P_b(Z,r3)

       E_sat(1,j)=E_tot(Z,r1,Temp)
       E_sat(2,j)=E_tot(Z,r3,Temp)
       call compute_Cv(Z,nu_sat(1,j),Temp,Cv_sat(1,j))
       call compute_Cv(Z,nu_sat(2,j),Temp,Cv_sat(2,j))
       !call compute_Cp(Z,nu_sat(1,j),Temp,Cp_sat(1,j))
       !call compute_Cp(Z,nu_sat(2,j),Temp,Cp_sat(2,j))
     
       !!!---v1: c(Cv) -> gam -> pinf 
       !call compute_soundspeed(Z,nu_sat(1,j),Temp,c_sat(1,j))
       !call compute_soundspeed(Z,nu_sat(2,j),Temp,c_sat(2,j))
       !call compute_gamma(Z,nu_sat(1,j),Temp,c_sat(1,j),g_sat(1,j))
       !call compute_gamma(Z,nu_sat(2,j),Temp,c_sat(2,j),g_sat(2,j))
       !pinf_sat(1,j)=(1.0_PR/nu_sat(1,j)*E_sat(1,j)*(g_sat(1,j)-1.0_PR)-P_sat(j))/g_sat(1,j)
       !pinf_sat(2,j)=(1.0_PR/nu_sat(2,j)*E_sat(2,j)*(g_sat(2,j)-1.0_PR)-P_sat(j))/g_sat(2,j)

       !!!---v2: gam (1/edP/dr+1) -> pinf ((re(g-1)-P)/g) -> c (sqrt(gam(p+pinf)/rho))
       call compute_gamma2(Z,nu_sat(1,j),Temp,g_sat(1,j))
       call compute_gamma2(Z,nu_sat(2,j),Temp,g_sat(2,j)) 
       call compute_pinf(Z,nu_sat(1,j),E_sat(1,j),g_sat(1,j),P,pinf_sat(1,j))
       call compute_pinf(Z,nu_sat(2,j),E_sat(2,j),g_sat(2,j),P,pinf_sat(2,j))
       call compute_soundspeed2(Z,nu_sat(1,j),P,g_sat(1,j),pinf_sat(1,j),c_sat(1,j))
       call compute_soundspeed2(Z,nu_sat(2,j),P,g_sat(2,j),pinf_sat(2,j),c_sat(2,j))
 

       if(mode_output.eq.2)then    
       call compute_Zbar(Z,nu_sat(1,j),Temp,Z_sat(1,j))
       call compute_Zbar(Z,nu_sat(2,j),Temp,Z_sat(2,j))
       elseif(mode_output.eq.1)then    
       call compute_Zbar2(Z,nu_sat(1,j),Temp,Z_sat(1,j))
       call compute_Zbar2(Z,nu_sat(2,j),Temp,Z_sat(2,j))
       endif
       call transport_properties(Z,r1,Temp,Z_sat(1,j),sig_sat(1,j),lamb_sat(1,j),Sth_sat(1,j),&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)
       call transport_properties(Z,r3,Temp,Z_sat(2,j),sig_sat(2,j),lamb_sat(2,j),Sth_sat(2,j),&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)

       call compute_Cp2(Cv_sat(1,j),g_sat(1,j),Cp_sat(1,j))
       call compute_Cp2(Cv_sat(2,j),g_sat(2,j),Cp_sat(2,j))
       

    endif
   
   
    !!!---écriture des valeurs le long de l'isoterme
    DO i=1,N_nu
   
       rho=rho0*(rhof/rho0)**(real(i-1,PR)/real(N_nu-1,PR))
       nu=1.0_PR/rho
   
       if(found(j).and.(rho.gt.r3.and.rho.lt.r1))then
   
         !call compute_Cv(Z,nu,Temp,Cvt)
         !!call compute_Cp(Z,nu,Temp,Cpt)
         !Cpt=0.0_PR
         !call compute_soundspeed(Z,nu,Temp,ct)
         !call compute_gamma(Z,nu,Temp,ct,gt)
         !call compute_Zbar(Z,nu,Temp,zt)
   
         Et=E_tot(Z,rho,Temp)
         Ebi=masslaw(nu_sat(1,j),nu,nu_sat(2,j),E_sat(1,j), E_sat(2,j))
         Cvbi=masslaw(nu_sat(1,j),nu,nu_sat(2,j),Cv_sat(1,j), Cv_sat(2,j))
         Cpbi=masslaw(nu_sat(1,j),nu,nu_sat(2,j),Cp_sat(1,j), Cp_sat(2,j))

         cbi=frozen_speed(nu_sat(1,j),nu,nu_sat(2,j),c_sat(1,j), c_sat(2,j))
         gbi=vollaw(nu_sat(1,j),nu,nu_sat(2,j),g_sat(1,j),g_sat(2,j))
         pinfbi=vollaw(nu_sat(1,j),nu,nu_sat(2,j),g_sat(1,j)*pinf_sat(1,j),g_sat(2,j)*pinf_sat(2,j))
         pinfbi=pinfbi/gbi

         zbi=masslaw(nu_sat(1,j),nu,nu_sat(2,j),z_sat(1,j), z_sat(2,j))
         sigi=sig_mixt(nu_sat(1,j),nu,nu_sat(2,j),sig_sat(1,j), sig_sat(2,j))
         lamb=sig_mixt(nu_sat(1,j),nu,nu_sat(2,j),lamb_sat(1,j), lamb_sat(2,j))
         Sth=sig_mixt(nu_sat(1,j),nu,nu_sat(2,j),Sth_sat(1,j), Sth_sat(2,j))
         Rs=nu*P/Temp

  
         if(mode_output.eq.1)then 
              write(id,"(14(ES14.7,' '))") nu, Pnu(nu), P, &
              Et,Ebi,Cvbi,Cpbi,cbi,gbi,zbi,sigi,lamb,Sth,pinfbi
         else      
              !write(id,"(7(ES14.7,' '))") P,Ebi,Cvbi,cbi,gbi,zbi,pinfbi
              write(id,"(7(ES14.7,' '))") P*1.0e6_PR,Ebi*1.0e3_PR,Cvbi*1.0e3_PR,cbi,gbi,zbi,pinfbi*1.0e6_PR

              if(P.lt.0.0_PR)then
                print*, 'ERREUR P<0:', i, j
                stop
              endif

         endif
    
         !write(id,"(15(ES14.7,' '))") nu, Pnu(nu), P, &
         !     Et,Ebi,Cvbi,Cpbi,cbi,gbi,zbi,sigi,lamb,Sth,pinfbi
    
   
         Vs(i)=nu
         Pression(i,j)=P
         Es(i,j)=Ebi
         Cv(i,j)=Cvbi
         Cp(i,j)=Cpbi
         sound(i,j)=cbi
         gam(i,j)=gbi
         w(i,j)=solid(Z,rho,Temp)
         Zbar(i,j)=zbi
         sigma(i,j)=sigi  
         lambda(i,j)=lamb 
         Seebeck(i,j)=Sth 
         Pression_e(i,j)=Pnu_e(nu)
         Pression_i(i,j)=Pnu_i(nu)
         Pression_b(i,j)=Pnu_b(nu)

       else
   
         Et=E_tot(Z,rho,Temp)

         call compute_Cv(Z,nu,Temp,Cvt)
         !call compute_Cp(Z,nu,Temp,Cpt)
                  !call compute_Zbar2(Z,nu,Temp,zt)


         !!!---v1
         !call compute_soundspeed(Z,nu,Temp,ct)
         !call compute_gamma(Z,nu,Temp,ct,gt)
         !Rs=nu*Pnu(nu)/Temp
         !pinf=(1.0_PR/nu*Et*(gt-1.0_PR)-Pnu(nu))/gt

         !!---v2
         call compute_gamma2(Z,nu,Temp,gt)
         call compute_pinf(Z,nu,Et,gt,Pnu(nu),pinf)
         call compute_soundspeed2(Z,nu,Pnu(nu),gt,pinf,ct)

         if(mode_output.eq.2)then    
         call compute_Zbar(Z,nu,Temp,zt)
         elseif(mode_output.eq.1)then    
         call compute_Zbar2(Z,nu,Temp,zt)
         endif

         call transport_properties(Z,rho,Temp,zt,sigt,lamb,Sth,&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)

         call compute_Cp2(Cvt,gt,Cpt)
         
         if(mode_output.eq.1)then 
              write(id,"(14(ES14.7,' '))") nu, Pnu(nu), Pnu(nu),&
              Et,Et,Cvt,Cpt,ct,gt,zt,sigt,lamb,Sth,pinf
         else
              write(id,"(7(ES14.7,' '))") Pnu(nu)*1.0e6_PR,Et*1.0e3_PR,Cvt*1.0e3_PR,ct,gt,zt,pinf*1.0e6_PR

              if(Pnu(nu).lt.0.0_PR)then
                print*, 'ERREUR P<0:', i, j, Pnu(nu)
                stop
              endif

         endif   

         Vs(i)=nu
         Pression(i,j)=Pnu(nu)
         Es(i,j)=Et
         Cv(i,j)=Cvt
         Cp(i,j)=Cpt
         sound(i,j)=ct
         gam(i,j)=gt
         w(i,j)=solid(Z,rho,Temp)
         Zbar(i,j)=zt
         sigma(i,j)=sigt  
         lambda(i,j)=lamb 
         Seebeck(i,j)=Sth 
         Pression_e(i,j)=Pnu_e(nu)
         Pression_i(i,j)=Pnu_i(nu)
         Pression_b(i,j)=Pnu_b(nu)

       endif
   
    ENDDO !!! nu
   
    write(id,*) ""
    write(id,*) ""
   ! write(12,*) ""
    write(iout,*) '-------------------------------------------------------------'
   
ENDDO !!!N_T


!!! ENDDO ; ENDDO !!! kn km à désactiver en utilisation standard !!! :O 

!close(10)
!close(12)


!!!---courbe de saturation

if(.not.varynm)then

write(id,'(A)') '#-------------------  Saturation curve  ----------------------'

k=0
do j=1,N_T
 if(found(j))then
 ! write(11,"(4(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(1,j), nu_sat(2,j)
  k=k+1
 endif
enddo
N_T_sat=k


print*, 'coucou Sat:', N_T_sat

write(id,'(A10,I3)') '# N_T_sat=', N_T_sat
if(mode_output.eq.1)then

   write(id,'(A)') '#----------------- Courbe de bulle ---------------------------'
   write(id,'(A)') "# T,P,nu,E,Cv,Cp,c,g,Z,sig,lambda,Seeb,pinf"
   write(id,'(A)') '#-------------------------------------------------------------'

   do j=1, N_T_sat
     write(id,"(13(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(1,j),&     !, nu_sat(2,j),&
                                                             E_sat(1,j),&      !, E_sat(2,j),&
                                                             Cv_sat(1,j),&     !, Cv_sat(2,j),&
                                                             Cp_sat(1,j),&     !, Cp_sat(2,j),&
                                                             c_sat(1,j),&      !, c_sat(2,j),&
                                                             g_sat(1,j),&      !, g_sat(2,j),&
                                                             Z_sat(1,j),&      !, Z_sat(2,j),&
                                                             sig_sat(1,j),&    !, sig_sat(2,j),&
                                                             lamb_sat(1,j),&   !, lamb_sat(2,j),&
                                                             Sth_sat(1,j),&    !, Sth_sat(2,j),&
                                                             pinf_sat(1,j)     !,pinf_sat(2,j)
      
   
   
   enddo

   j=N_T_sat
   write(id,"(13(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&       !, nu_sat(2,j),&
                                                           E_sat(2,j),&        !, E_sat(2,j),&
                                                           Cv_sat(2,j),&       !, Cv_sat(2,j),&
                                                           Cp_sat(2,j),&       !, Cp_sat(2,j),&
                                                           c_sat(2,j),&        !, c_sat(2,j),&
                                                           g_sat(2,j),&        !, g_sat(2,j),&
                                                           Z_sat(2,j),&        !, Z_sat(2,j),&
                                                           sig_sat(2,j),&      !, sig_sat(2,j),&
                                                           lamb_sat(2,j),&     !, lamb_sat(2,j),&
                                                           Sth_sat(2,j),&      !, Sth_sat(2,j),&
                                                           pinf_sat(2,j)       !,pinf_sat(2,j)

  else !!! mode_output

   write(id,'(A)') '#----------------- Courbe de bulle ---------------------------'
   write(id,'(A)') "# T,P,nu,E,Cv,c,g,Z,pinf,pe,pi,pb"
   write(id,'(A)') '#-------------------------------------------------------------'

   do j=1, N_T_sat

     write(id,"(12(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(1,j),&     !, nu_sat(2,j),&
                                                             E_sat(1,j),&      !, E_sat(2,j),&
                                                             Cv_sat(1,j),&     !, Cv_sat(2,j),&
                                                             c_sat(1,j),&      !, c_sat(2,j),&
                                                             g_sat(1,j),&      !, g_sat(2,j),&
                                                             Z_sat(1,j),&      !, Z_sat(2,j),&
                                                             pinf_sat(1,j),&     !,pinf_sat(2,j)&
                                                             Pe_sat(1,j),&
                                                             Pi_sat(1,j),&
                                                             Pb_sat(1,j)
      
   enddo

   j=N_T_sat
   write(id,"(12(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&       !, nu_sat(2,j),&
                                                           E_sat(2,j),&        !, E_sat(2,j),&
                                                           Cv_sat(2,j),&       !, Cv_sat(2,j),&
                                                           c_sat(2,j),&        !, c_sat(2,j),&
                                                           g_sat(2,j),&        !, g_sat(2,j),&
                                                           Z_sat(2,j),&        !, Z_sat(2,j),&
                                                           pinf_sat(2,j),&       !,pinf_sat(2,j)
                                                           Pe_sat(1,j),&
                                                           Pi_sat(1,j),&
                                                           Pb_sat(1,j)
  endif

write(id,*) ""
write(id,*) ""

  if(mode_output.eq.1)then

    write(id,'(A)') '#------------------- Courbe de goute -------------------------'
    write(id,'(A)') "# T,P,nu,E,Cv,Cp,c,g,Z,sig,lambda,Seeb,pinf,pe,pi,pb"
    write(id,'(A)') '#-------------------------------------------------------------'
    do j=1, N_T_sat
    write(id,"(13(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&     !, nu_sat(2,j),&
                                                          E_sat(2,j),&      !, E_sat(2,j),&
                                                          Cv_sat(2,j),&     !, Cv_sat(2,j),&
                                                          Cp_sat(2,j),&     !, Cp_sat(2,j),&
                                                          c_sat(2,j),&      !, c_sat(2,j),&
                                                          g_sat(2,j),&      !, g_sat(2,j),&
                                                          Z_sat(2,j),&      !, Z_sat(2,j),&
                                                          sig_sat(2,j),&    !, sig_sat(2,j),&
                                                          lamb_sat(2,j),&   !, lamb_sat(2,j),&
                                                          Sth_sat(2,j),&    !, Sth_sat(2,j),&
                                                          pinf_sat(2,j)     !,pinf_sat(2,j)
   

 
    enddo

    j=N_T_sat
    write(id,"(13(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&       !, nu_sat(2,j),&
                                                            E_sat(2,j),&        !, E_sat(2,j),&
                                                            Cv_sat(2,j),&       !, Cv_sat(2,j),&
                                                            Cp_sat(2,j),&       !, Cp_sat(2,j),&
                                                            c_sat(2,j),&        !, c_sat(2,j),&
                                                            g_sat(2,j),&        !, g_sat(2,j),&
                                                            Z_sat(2,j),&        !, Z_sat(2,j),&
                                                            sig_sat(2,j),&      !, sig_sat(2,j),&
                                                            lamb_sat(2,j),&     !, lamb_sat(2,j),&
                                                            Sth_sat(2,j),&      !, Sth_sat(2,j),&
                                                            pinf_sat(2,j)       !,pinf_sat(2,j)


   else !!! mode_output

    write(id,'(A)') '#------------------- Courbe de goute -------------------------'
    write(id,'(A)') "# T,P,nu,E,Cv,c,g,Z,pinf"
    write(id,'(A)') '#-------------------------------------------------------------'
    do j=1, N_T_sat
    write(id,"(12(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&     !, nu_sat(2,j),&
                                                          E_sat(2,j),&      !, E_sat(2,j),&
                                                          Cv_sat(2,j),&     !, Cv_sat(2,j),&
                                                          c_sat(2,j),&      !, c_sat(2,j),&
                                                          g_sat(2,j),&      !, g_sat(2,j),&
                                                          Z_sat(2,j),&      !, Z_sat(2,j),&
                                                          pinf_sat(2,j),&     !,pinf_sat(2,j)
                                                          Pe_sat(1,j),&
                                                          Pi_sat(1,j),&
                                                          Pb_sat(1,j)
   

 
    enddo

    j=N_T_sat
    write(id,"(12(ES14.7,' '))") temperatures(j), P_sat(j), nu_sat(2,j),&       !, nu_sat(2,j),&
                                                            E_sat(2,j),&        !, E_sat(2,j),&
                                                            Cv_sat(2,j),&       !, Cv_sat(2,j),&
                                                            c_sat(2,j),&        !, c_sat(2,j),&
                                                            g_sat(2,j),&        !, g_sat(2,j),&
                                                            Z_sat(2,j),&        !, Z_sat(2,j),&
                                                            pinf_sat(2,j),&       !,pinf_sat(2,j)
                                                            Pe_sat(2,j),&
                                                            Pi_sat(2,j),&
                                                            Pb_sat(2,j)
 



    endif

write(id,*) ""
write(id,*) ""
 
write(id,'(A)') '#------------------ Melting curve curve ----------------------'
write(id,'(A)') "# nu,Tm"
write(id,'(A)') '#-------------------------------------------------------------'

    DO i=1,N_nu
   
       rho=rho0*(rhof/rho0)**(real(i-1,PR)/real(N_nu,PR))
       nu=1.0_PR/rho

       write(id,"(2(ES14.7,' '))") nu, Tm(Z,rho) 

 
    ENDDO

endif !!!(varynm)

close(10)

!close(11)

!call init_snes(N,x,b,F_MAXWELL,J_MAXWELL) !!! N=2
!
!call solve_SNES(N,x,b)
!
!print*, 'x=', x
!print*, 'nu1=', 1.0_PR/x(1), 'nu2=', 1.0_PR/x(2)

if(mode_output.eq.1)then
call output_ensight(N_T,N_nu,nu_sat(1:2,1:N_T),temperatures(1:N_T),Vs(1:N_nu),Pression(1:N_nu,1:N_T),Es(1:N_nu,1:N_T),&
                             Cv(1:N_nu,1:N_T),Cp(1:N_nu,1:N_T),sound(1:N_nu,1:N_T),gam(1:N_nu,1:N_T),&
                             w(1:N_nu,1:N_T),Zbar(1:N_nu,1:N_T),sigma(1:N_nu,1:N_T),lambda(1:N_nu,1:N_T),&
                             Seebeck(1:N_nu,1:N_T),1)

call output_ensight(N_T,N_nu,nu_sat(1:2,1:N_T),temperatures(1:N_T),Vs(1:N_nu),Pression(1:N_nu,1:N_T),Es(1:N_nu,1:N_T),&
                             Cv(1:N_nu,1:N_T),Cp(1:N_nu,1:N_T),sound(1:N_nu,1:N_T),gam(1:N_nu,1:N_T),&
                             w(1:N_nu,1:N_T),Zbar(1:N_nu,1:N_T),sigma(1:N_nu,1:N_T),lambda(1:N_nu,1:N_T),&
                             Seebeck(1:N_nu,1:N_T),2)

endif

write(iout,*) '==========================================='
write(iout,*) '       REMOVE VANDERWAALS LOOP FINI :)'
write(iout,*) '==========================================='
write(iout,*) ""



ELSEIF(MODE.eq.2)THEN

write(iout,*) ""
write(iout,*) '==========================================='
write(iout,*) '            READ TABLE !!!                 '
write(iout,*) '==========================================='

!call read_table_NT('tab_Cu_400.dat',mymat)

call read_table_NT('tab_Fe_36.dat',mymat)

print*, 'READ TABLE FINI !!!'
print*, '--------------- test 1 -------------------'
call get_index_rT(mymat,7.8_PR,300.0_PR,i,iT)
print*, 'ir=', i, 'vs 701  iT=', iT, 'vs 125'

!print*, '--------------- test 1 -------------------'
!call get_index_rT(mymat,0.287_PR,3310.0_PR,i,iT)
!print*, 'ir=', i, 'vs 701  iT=', iT, 'vs 125'
!print*, '--------------- test 2 -------------------'
!call get_index_rT(mymat,8.93e-5_PR,201.0_PR,i,iT)
!print*, 'ir=', i, 'vs 1  iT=', iT, 'vs 1'
!print*, '--------------- test 3 -------------------'
!rho=8.80_PR ; Temp=10170.0_PR
!call get_index_rT(mymat,rho,Temp,i,iT)
!print*, 'ir=', i, 'vs 998  iT=', iT, 'vs 399'
!print*, 'P1=', mymat%P(i,iT+1), mymat%P(i+1,iT+1)
!print*, 'P2=', mymat%P(i,iT)  , mymat%P(i+1,iT)
!call get_from_rT_lin(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plin=', P
!call get_from_rT_log(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plog=', P
!print*, '--------------- test 4 -------------------'
!rho=1.0e12_PR ; Temp=1e20_PR
!call get_index_rT(mymat,rho,Temp,i,iT)
!print*, 'ir=', i, 'vs 1000  iT=', iT, 'vs 400'
!call get_from_rT_lin(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plin=', P, 'vs ', mymat%P(i,iT)
!call get_from_rT_log(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plog=', P, 'vs ', mymat%P(i,iT)
!print*, '--------------- test 5 -------------------'
!rho=1.0e-12_PR ; Temp=1e-20_PR
!call get_index_rT(mymat,rho,Temp,i,iT)
!print*, 'ir=', i, 'vs 1  iT=', iT, 'vs 1'
!call get_from_rT_lin(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plin=', P, 'vs ', mymat%P(i,iT)
!call get_from_rT_log(mymat,rho,Temp,i,iT,mymat%P,P)
!print*, 'Plog=', P, 'vs', mymat%P(i,iT)
!print*, '--------------- test 6 -------------------'
!rho=1.57e-4_PR ; Ener=6.11e3_PR  !Temp=3880.0_PR
!call get_index_rE(mymat,rho,Ener,i,iT)
!print*, 'ir=', i, 'vs 50  iT=', iT, 'vs 148'
!print*, 'T1=', mymat%T(iT)   , 'vs 3875'
!print*, 'r1=', mymat%rho(i) , 'vs 1.56895359'
!print*, '--------------- test 7 -------------------'
!rho=7.1_PR ; Ener=3.20e3_PR  !Temp=3880.0_PR
!call get_index_rE(mymat,rho,Ener,i,iT)
!print*, 'ir=', i, 'vs 980  iT=', iT, 'vs 352'
!print*, 'T1=', mymat%T(iT)   , 'vs 8975'
!print*, 'r1=', mymat%rho(i) , 'vs 7.0837749'
!print*, '--------------- test 8 -------------------'
!rho=7.1_PR ; Ener=3.20e6_PR  !Temp=3880.0_PR
!call get_index_rE(mymat,rho,Ener,i,iT)
!print*, 'ir=', i, 'vs 980  iT=', iT, 'vs 400'
!print*, 'E=', mymat%E(i,iT)   , 'vs 3635'
!print*, 'r1=', mymat%rho(i) , 'vs 7.0837749'
!print*, '--------------- test 9 -------------------'
!rho=2.0_PR ; P=15.5020635_PR  !Temp=3527.0_PR ->T1=3525,T2=3550
!call get_index_rP(mymat,rho,P,i,iT)
!print*, 'ir=', i, 'vs 870  iT=', iT, 'vs 134'
!print*, 'P(T1)=', mymat%P(i,iT)   , 'vs 15.088664'
!print*, 'P(T2)=', mymat%P(i,iT+1)   , 'vs 15.915463'
!print*, 'r1=', mymat%rho(i) , 'vs 1.9939497'
!call check_tab




call check_tab





ELSEIF(MODE.eq.3)THEN  

CALL INIT_EOS

open(unit=2, file='test_transport.dat', status='replace')

N_T=2000
T0=0.01_PR*eV
Tf=1000.0_PR*eV

!Z=13
Z=29
!rho=0.215_PR
!rho=2.5_PR
!rho=rhos(Z)/100.0_PR
N_nu=7
rpar(1:7)=(/1.0e-3_PR, 1.0e-2_PR, 1.0e-1_PR, 1.0_PR, 1.0e1_PR, 1.0e2_PR, 1.0e3_PR/)
!N_nu=2
!rpar(1:2)=(/0.215_PR, 2.5_PR/)
!rpar(1:2)=(/rhos(Z)/100.0_PR, rhos(Z)/)
N_nu=1000 ; rho0=1.0e-2_PR ; rhof=(10.0_PR)**(1.5_PR)

N_T=1
rpar(1:3)=(/1.0e4_PR, 2.0e4_PR, 3.0e4_PR/)
rpar(1:1)=(/6.0e3_PR/)

do it=1,N_T

T=rpar(it)   

do i=1,N_nu

      !rho=rpar(i)
      rho=rho0*(rhof/rho0)**(real(i-1,PR)/real(N_nu,PR))
      nu=1.0_PR/rho

   !do it=1,3 !!!N_T+1
   
      !T=T0*(Tf/T0)**(real(it-1,PR)/real(N_T,PR))
   

      call compute_Zbar(Z,nu,T,zt)
      call compute_ZFS(Z,rho,T,zt,zs,zfs)
 
      call transport_properties(Z,rho,T,ZFS,sig,lamb,Sth,&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)
   
      write(2,"(21(ES14.7,' '))") rho, T, sig, lamb, Sth, sig1, sig2, zt,zs,zfs, na, ne, nn, ni,&
                                  lnl, tau,tau1,tau2,tau_en,cross_en, mu
  
   !enddo   

    !write(2,*) ""
    !write(2,*) ""

enddo
    write(2,*) ""
    write(2,*) ""
enddo

   !!!---test_Aalpha
   !open(unit=3,file='testAalpha.dat',status='replace')
   !N_T=1000
   !T0=-10.0_PR
   !Tf=100.0_PR
   !do it=1,N_T
   !   !T=T0*(Tf/T0)**(real(it-1,PR)/real(N_T,PR))
   !   T=T0+real(it-1,PR)*(Tf-T0)/real(N_T,PR)
   !   write(3,"(4(ES14.7,' '))") T, A_alpha(T), A_alpha2(T), fd1h(T)*(1.0_PR+exp(-T)) 
   !enddo
   !close(3)

   !!!---test_Z
   Z=29
   rho=1.0e3_PR
   nu=1.0_PR/rho

   N_nu=7
   rpar(1:7)=(/1.0e-3_PR, 1.0e-2_PR, 1.0e-1_PR, 1.0_PR, 1.0e1_PR, 1.0e2_PR, 1.0e3_PR/)
   open(unit=3,file='testZ.dat',status='replace')
   N_T=1000
   T0=1.0e3_PR
   Tf=1.0e8_PR

   do i=1,N_nu

      rho=rpar(i)
      nu=1.0_PR/rho
     
      do it=1,N_T
         T=T0*(Tf/T0)**(real(it-1,PR)/real(N_T,PR))
         !T=T0+real(it-1,PR)*(Tf-T0)/real(N_T,PR)
         call compute_Zbar(Z,nu,T,zt)
         Te=T !!!10.0_PR*T  !!10000.0_PR
         call compute_ZFS(Z,rho,Te,zt,zs,zfs) 
         call transport_properties(Z,rho,Te,ZFS,sig,lamb,Sth,&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)

         write(3,"(5(ES14.7,' '))") T, zt, zs, zfs,sig
      enddo
      write(3,*) "" 
      write(3,*) "" 

   enddo

   close(3)


ELSEIF(MODE.EQ.4)THEN

CALL INIT_EOS

Z=13 !!29 !! 13
!mat%Z=Z

!!!---scaling of conductivity

rho=rhos(Z) ; nu=1.0_PR/rho ; T=300.0_PR

call compute_Zbar(Z,nu,T,zt)

call compute_ZFS(Z,rho,T,zt,zs,zfs)
 
call transport_properties(Z,rho,T,ZFS,sig,lamb,Sth,&
                sig1,sig2,K1,K2,na,ne,nn,ni,lnl,tau,tau1,tau2,tau_en,cross_en,mu)
 

 print*, 'scaling of conductivity for ', element_name(Z) 

 print*, 'sigmamax=', sig2
 print*, 'lambdamax=', K2

ELSEIF(MODE.EQ.5)THEN

print*, 'MODE=5: test Ptot >0'

CALL INIT_EOS

Z=13 !!29 !! 13

N_T=401 ; T0=200.0_PR ; Tf=10200.0_PR ; allocate(temperatures(1:N_T)) !!! 25 K 200->10000 
   do it=1,N_T
       temperatures(it)=T0+real(it-1,PR)*(Tf-T0)/real(N_T-1,PR)
   enddo
   
   N_nu=1000 ; rho0=1.0e-5_PR*rhos(Z) ; rhof=2.0_PR*rhos(Z)

   Temp=temperatures(77) 


    DO i=1,N_nu
   
       rho=rho0*(rhof/rho0)**(real(i-1,PR)/real(N_nu-1,PR))
       nu=1.0_PR/rho

       print*, P_tot(Z,1.0_PR/nu,Temp)

    ENDDO


 
ENDIF


!rhof=1.0_PR*rhos(Z) ; rho0=1.0e-2_PR*rhos(Z) ; Nr=1000
!
!allocate(Cv(N_P,N_r))
!
!
!do i=1,N_P
!
!
!do j=1,N_P
!  do i=1,N_r
!
!    rho=rho0*(rhof/rho0)**(real(i-1,PR)/real(N_nu,PR))
!    nu=1.0_PR/rho
!
!    
!
!    do while 
!










!!!------------TEST NASG----------------------------


!!!!---initialization des matériaux
!
!call init_mat
!
!mat=Cu
!
!!!!---calcul
!
!!!!---AEB
!!!!---initialization de SNES
!N=2
!allocate(x(1:N)) ; x=0.0_PR
!allocate(b(1:N)) ; b=0.0_PR
!
!x(:)=10.0_PR
!
!!call init_snes(N,x,b,F_ABE,J_ABE) !!! N=3
!!mat%A=x(1) ; mat%B=x(2) ; mat%E=x(3)
!
!call init_snes(N,x,b,F_AB,J_AB) !!! N=2
!
!call solve_SNES(N,x,b)
!
!mat%A=x(1) ; mat%B=x(2)
!
!print*, 'x=', x
!
!
!print*, log(mat%Pc)-(mat%A+(mat%B+mat%E*mat%Pc)/mat%Tc+mat%C*log(mat%Tc)+mat%D*log(mat%Pc+mat%pi))
!
!print*, log(mat%Pvap0)-(mat%A+(mat%B+mat%E*mat%Pvap0)/mat%Tvap0+mat%C*log(mat%Tvap0)+mat%D*log(mat%Pvap0+mat%pi))
!
!print*, log(mat%Pvap1)-(mat%A+(mat%B+mat%E*mat%Pvap1)/mat%Tvap1+mat%C*log(mat%Tvap1)+mat%D*log(mat%Pvap1+mat%pi))
!
!
!mat%qpg=mat%A*(mat%cpg-mat%cvg)-(mat%cpl-mat%cpg)+mat%qpl
!
!mat%ql=mat%B*(mat%cpg-mat%cvg)+mat%qg
!
!mat%bl=mat%E*(mat%cpg-mat%cvg)+mat%bg
!
!
!print*, 'qpg=', mat%qpg, 'J / kg / K'
!
!print*, 'ql=',  mat%ql, 'J / kg'
!
!print*, 'bl=', mat%bl, 'm^3 / kg'
!
!print*, 'vérif:'
!
!
!!restart_snes=.true.    
!!call init_snes(N,x,b,F_Pvap,J_Pvap)
!
!print*, 'Pc:   ', Pvap(mat%Tc), 'vs', mat%Pc
!print*, 'Pvap0:', Pvap(mat%Tvap0), 'vs', mat%Pvap0 
!print*, 'Pvap1:', Pvap(mat%Tvap1), 'vs', mat%Pvap1 
!
!
!!!!---
!
!open(unit=12, file='output.dat', status='replace')
!N_T=100
!Temp=300.0_PR
!DTemp=(mat%Tc-Temp)/real(N_T,PR)
!
!do i=0,N_T
!
!   P=Pvap(Temp)
!
!   write(12,"(4(ES14.7,' '))") Temp, P, rhol(P,Temp), rhog(P,Temp) 
!
!   Temp=Temp+DTemp
!enddo
!
!
!
!close(12)







!!!!---Pvap
!!!!---initialization de SNES
!N=1
!allocate(x(1:N)) ; x=0.0_PR
!allocate(b(1:N)) ; b=0.0_PR
!
!x(1)=1.0e10_PR
!
!
!call init_snes(N,x,b,F_Pvap,J_Pvap)
!
!N_T=1
!Temp=1000.0_PR
!DTemp=(3000.0_PR-Temp)/real(N_T,PR)
!
!do i=1,N_T
!
!   x(1)=1.0e10_PR
!   b=0.0_PR
!
!   call solve_SNES(N,x,b)
!  
!
!   write(12,*) Temp, x(1),  f, dfdx
!
!   call F_Pvap(N,x(1),f(1))
!   call J_Pvap(N,x(1),dfdx(1))
!
!   print*, Temp, x(1)*1e-5_PR,  f, dfdx
!
!   Temp=Temp+DTemp
!enddo
!
!close(12)




!!!---TEST 2 :
!
!N=2
!
!allocate(x(1:N)) ; x=0.0_PR
!allocate(b(1:N)) ; b=0.0_PR
!
!x(1)=1.0_PR
!x(2)=5.0_PR
!
!write(iout,*) 'TEST 2 solution: x = 1.07716 0.0897632'
!
!call init_snes(N,x,b,Fun,Jac)
!
!call solve_SNES(N,x,b)
!



contains

!real(PR) function rhol(P,T)
!    implicit none
!    double precision, intent(in) :: P, T
!    
!    rhol=(P+mat%pi)/((mat%gl-1.0_PR)*mat%cvl*T)
!
!end function rhol
!
!real(PR) function rhog(P,T)
!    implicit none
!    double precision, intent(in) :: P, T
!    
!    rhog=P/((mat%gg-1.0_PR)*mat%cvg*T)
!
!end function rhog




!subroutine F_AB(N,x,f)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: f(1:N)
!     real(PR) :: A,B,C,D,E,T,pinf,Tc,Pc,T1,P1,T0,P0
!
!     pinf=mat%pi
!     Tc=mat%Tc
!     Pc=mat%Pc
!     T0=mat%Tvap0
!     P0=mat%Pvap0
!     T1=mat%Tvap1
!     P1=mat%Pvap1
!     C=mat%C
!     D=mat%D
!
!     A=x(1)
!     B=x(2)
!
!     f(1)=log(Pc)-(A+B/Tc+C*log(Tc)+D*log(Pc+pinf))
!
!     f(2)=log(P0)-(A+B/T0+C*log(T0)+D*log(P0+pinf))
!
!end subroutine F_AB
!
!subroutine J_AB(N,x,J)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: J(1:N,1:N)
!     real(PR) :: A,B,C,D,E,T,pinf,Tc,Pc,T1,P1,T0,P0
!
!     pinf=mat%pi
!     Tc=mat%Tc
!     Pc=mat%Pc
!     T0=mat%Tvap0
!     P0=mat%Pvap0
!     T1=mat%Tvap1
!     P1=mat%Pvap1
!     C=mat%C
!     D=mat%D
!
!     A=x(1)
!     B=x(2)
!
!     J(1,1)=-1.0_PR
!     J(1,2)=-1.0_PR/Tc
!
!     J(2,1)=-1.0_PR
!     J(2,2)=-1.0_PR/T0
!
!end subroutine J_AB
!
!
!subroutine F_ABE(N,x,f)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: f(1:N)
!     real(PR) :: A,B,C,D,E,T,pinf,Tc,Pc,T1,P1,T0,P0
!
!     pinf=mat%pi
!     Tc=mat%Tc
!     Pc=mat%Pc
!     T0=mat%Tvap0
!     P0=mat%Pvap0
!     T1=mat%Tvap1
!     P1=mat%Pvap1
!     C=mat%C
!     D=mat%D
!
!     A=x(1)
!     B=x(2)
!     E=x(3)
!
!     f(1)=log(Pc)-(A+(B+E*Pc)/Tc+C*log(Tc)+D*log(Pc+pinf))
!
!     f(2)=log(P0)-(A+(B+E*P0)/T0+C*log(T0)+D*log(P0+pinf))
!
!     f(3)=log(P1)-(A+(B+E*P1)/T1+C*log(T1)+D*log(P1+pinf))
!
!end subroutine F_ABE
!
!subroutine J_ABE(N,x,J)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: J(1:N,1:N)
!     real(PR) :: A,B,C,D,E,T,pinf,Tc,Pc,T1,P1,T0,P0
!
!     pinf=mat%pi
!     Tc=mat%Tc
!     Pc=mat%Pc
!     T0=mat%Tvap0
!     P0=mat%Pvap0
!     T1=mat%Tvap1
!     P1=mat%Pvap1
!     C=mat%C
!     D=mat%D
!
!     A=x(1)
!     B=x(2)
!     E=x(3)
!
!     J(1,1)=-1.0_PR
!     J(1,2)=-1.0_PR/Tc
!     J(1,3)=-Pc/Tc
!
!     J(2,1)=-1.0_PR
!     J(2,2)=-1.0_PR/T0
!     J(2,3)=-P0/T0
!
!     J(3,1)=-1.0_PR
!     J(3,2)=-1.0_PR/T1
!     J(3,3)=-P1/T1
!
!end subroutine J_ABE

!subroutine F_ABEP(N,x,f)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: f(1:N)
!     real(PR) :: A,B,C,D,E,T,pinf
!
!     A=x(1)
!     B=x(2)
!     E=x(3)
!     pinf=x(4)
!
!     f(1)=log(Pc)-(A+(B+E*Pc)/Tc+C*log(Tc)+D*log(Pc+pinf))
!
!     f(2)=log(P0)-(A+(B+E*P0)/T0+C*log(T0)+D*log(P0+pinf))
!
!     f(3)=log(P1)-(A+(B+E*P1)/T1+C*log(T1)+D*log(P1+pinf))
!
!     f(4)=Cv0-(1/rho0-E*(cpg-cvg))*(p0+pinf)/((g-1)*T0)     
!   
!
!
!end subroutine F_ABEP
!

!real(PR) function Pvap(T) 
!
!    use mod_data, only : Temp, mat
!
!    implicit none
!    real(PR), intent(in) :: T
!    real(PR) :: x(1:1), b(1:1)
!    integer :: N
!
!    Temp=T    
!
!    print*, '------------------------------'
!    print*, 'T=', T
!    print*, '------------------------------'
!
!    N=1
!    x=0.0_PR ; b=0.0_PR   
! 
!    x(1)=2.0_PR*mat%Pc  !!1.0e8_PR
!
!    !x(1)=mat%Pvap0 !!1.0e8_PR
!
!    !call solve_SNES(N,x(1:1),b(1:1))
!    call My_Newton(1,x(1:1),1.0e-150_PR,3.0_PR*mat%Pc,F_Pvap,J_Pvap)
!
!    Pvap=x(1)
!
!end function Pvap
!
!
!subroutine F_Pvap(N,x,f)
!
!     use mod_data, only : mat, Temp
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: f(1:N)
!     real(PR) :: A,B,C,D,E,T,pinf
!
!     !!--- Psat
!     pinf=mat%pi
!     T=Temp
!     A=mat%A
!     B=mat%B
!     C=mat%C
!     D=mat%D
!     E=mat%E
!
!     f(1)=log(x(1))-(A+(B+E*x(1))/T+C*log(T)+D*log(x(1)+pinf))
! 
!     !!---exemple 1 : x(1)=1 ; x(2)=2
!     !f(1)=x(1)**2+x(1)*x(2)-3.0_PR
!     !f(2)=x(1)*x(2)+x(2)**2-6.0_PR
!     !!---exemple 2: result: x=1.07716 ; y=0.0897632
!     !f(1)=sin(3.0_PR*x(1))+x(2)
!     !f(2)=x(1)*x(2)-12.0_PR*x(2)**2
!    
!end subroutine F_Pvap
!
!
!subroutine J_Pvap(N,x,J)
!
!     use mod_data, only : mat, Temp
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: J(1:N,1:N)
!     real(PR) :: T, D, E, pinf
!
!     pinf=mat%pi
!     T=Temp
!     D=mat%D
!     E=mat%E
!
!     J(1,1)=1.0_PR/x(1)-(E/T+D/(x(1)+pinf))
!
!     !print*, 'J(1,1)=', J(1,1)
!
!     !!!---exemple 1
!     !J(1,1)=2.0_PR*x(1)+x(2) ; J(1,2)=x(1)
!     !J(2,1)=x(2)             ; J(2,2)=x(1)+2.0_PR*x(2)
!
!     !!!---exemple 2:
!     !J(1,1)=3.0_PR*cos(3.0_PR*x(1)) ; J(1,2)=1.0_PR
!     !J(2,1)=x(2)                    ; J(2,2)=x(1)-24.0_PR*x(2)
!
!end subroutine J_Pvap
!
!!!!===================== TRAITEMENT TABLE ===========================


real(PR) function Pnu(nu)

   implicit none
   real(PR), intent(in) :: nu
  
   Pnu=P_tot(Z,1.0_PR/nu,Temp)

end function Pnu

real(PR) function Pnu_e(nu)

   implicit none
   real(PR), intent(in) :: nu
  
   Pnu_e=P_e(Z,1.0_PR/nu,Temp)

end function Pnu_e

real(PR) function Pnu_i(nu)

   implicit none
   real(PR), intent(in) :: nu
  
   Pnu_i=P_i(Z,1.0_PR/nu,Temp)

end function Pnu_i

real(PR) function Pnu_b(nu)

   implicit none
   real(PR), intent(in) :: nu
  
   Pnu_b=P_b(Z,1.0_PR/nu)

end function Pnu_b

!subroutine F_MAXWELL(N,x,f)
!
!     use mod_eos, only : P_tot
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: f(1:N)
!
!     f(1)=Pnu(x(1))-Pnu(x(2))
!     f(2)=int_nu(Pnu,x(1),x(2))-Pnu(x(1))*(x(2)-x(1))
! 
!end subroutine F_MAXWELL
!
!subroutine J_MAXWELL(N,x,J)
!
!     implicit none
!     integer, intent(in) :: N
!     real(PR), intent(in) :: x(1:N)
!     real(PR), intent(out) :: J(1:N,1:N)
!
!     J(1,1)=ddnu(Pnu,x(1)) ; J(1,2)=-ddnu(Pnu,x(2))
!
!     J(2,1)=-ddnu(Pnu,x(1))*(x(2)-x(1)) ; J(2,2)=Pnu(x(2))-Pnu(x(1))
!
!
!end subroutine J_MAXWELL






!!!===================== ROUTINES MATHEMATIQUES =====================



!!!======================  DERIVEES  ================

real(PR) FUNCTION ddr(Z,rho,T,f)
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   interface
     double precision function fun(Z,rho,T)
     integer, intent(in) :: Z
     double precision, intent(in) :: rho,T
     end function fun
   end interface
   procedure(fun) :: f
   real(PR) :: dr
   
   dr=rho*1.0e-3_PR
   
   ddr=(f(Z,rho+dr,T)-f(Z,rho,T))/dr

end function ddr

real(PR) FUNCTION ddnu(f,nu)
   implicit none
   interface
     double precision function fun(nu)
     double precision, intent(in) :: nu
     end function fun
   end interface
   procedure(fun) :: f
   real(PR), intent(in) :: nu
   real(PR) :: dnu,nu1,nu2
  

   dnu=nu*1.0e-3_PR 
   nu1=nu
   nu2=nu+dnu

   ddnu=(f(nu2)-f(nu1))/dnu

end function ddnu

real(PR) FUNCTION int_nu(f,numin,numax)
   implicit none
   interface
     double precision function fun(nu)
     double precision, intent(in) ::nu 
     end function fun
   end interface
   procedure(fun) :: f
   real(PR), intent(in) :: numin, numax
   real(PR) :: nu,dnu,f1,f2, nu1,nu2
   integer :: i, N_nu 
 
   int_nu=0.0_PR

   N_nu=1000
   dnu=(numin-numax)/real(N_nu,PR)

   nu1=numin
   f1=f(nu1)

   do i=1,N_nu

     !!!---pas constant---------------------------------
     !nu2=nu1+dnu
     !!!---pas geo--------------------------------------
     nu2=numin*(numax/numin)**(real(i,PR)/real(N_nu,PR))
     !!!------------------------------------------------

     f2=f(nu2)

     int_nu=int_nu+(0.5_PR*(f2-f1)+f1)*(nu2-nu1)
 
     nu1=nu2
     f1=f2

   enddo

end function int_nu

real(PR) function masslaw(nu1,nu,nu2,f1,f2)

   implicit none
   real(PR), intent(in) :: nu1, nu, nu2
   real(PR), intent(in) :: f1, f2

   masslaw=( (nu2-nu)*f1 + (nu-nu1)*f2 )/(nu2-nu1)

end function masslaw


real(PR) function frozen_speed(nu1,nu,nu2,c1,c2)

   implicit none
   real(PR), intent(in) :: nu1, nu, nu2
   real(PR), intent(in) :: c1, c2

   !frozen_speed=sqrt( ( (nu2-nu)*nu/nu1*c1**2 + (nu-nu1)*nu/nu2*c2**2 )/(nu2-nu1) )
   frozen_speed=sqrt( ( (nu2-nu)*c1**2 + (nu-nu1)*c2**2 )/(nu2-nu1) )

end function frozen_speed


real(PR) function vollaw(nu1,nu,nu2,f1,f2)

   implicit none
   real(PR), intent(in) :: nu1, nu, nu2
   real(PR), intent(in) :: f1, f2
   real(PR) :: x1,x2,z1,z2

   z1=(nu2-nu)/(nu2-nu1) ; z2=(nu-nu1)/(nu2-nu1)

   x1=nu1/nu*z1 ; x2=nu2/nu*z2

   vollaw=x1*f1+x2*f2

end function vollaw


real(PR) function sig_mixt(nu1,nu,nu2,sig1,sig2)

   implicit none
   real(PR), intent(in) :: nu1, nu, nu2
   real(PR), intent(in) :: sig1, sig2
   real(PR) :: x1,x2,z1,z2
   real(PR) :: n

   !!!---Cherkas 2011 equation (20)
   z1=(nu2-nu)/(nu2-nu1) ; z2=(nu-nu1)/(nu2-nu1)

   x1=nu1/nu*z1 ; x2=nu2/nu*z2


   !!!---in serie (n=infiny)
   !sig_mixt=x1*sig1+x2*sig2

   !!!---attached phases (n=1)
   sig_mixt=sig1*sig2/(x1*sig2+x2*sig1)

   !!!---heterogeneous phase with foam-like structures
   !!! (n=3 is an exemple taken from Cherkas 2011 without justification...)
   !n=3.0_PR
   !sig_mixt=sig1*(1.0_PR+x2/(x1/n-sig1/(sig1-sig2)))
   



end function sig_mixt



!!!====================== Cp, Cv ==================

subroutine compute_Cv(Z,nu,T,Cv)
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T
   real(PR), intent(out) :: Cv
   real(PR) :: rho, E1, E2, E

   rho=1.0_PR/nu

   Cv=ddTr_E_tot(Z,rho,T)

   !dT=T*0.1_PR
   !rho=1.0_PR/nu

   !E1=E_tot(Z,rho,T-dT)
   !E2=E_tot(Z,rho,T+dT)
   !E=E_tot(Z,rho,T)


   !!Cv=(E2-E1)/dT
   ! Cv=min((E2-E)/dT,(E-E1)/dT)

end subroutine compute_Cv

subroutine compute_Cp2(Cv,gamm,Cp)

   implicit none
   real(PR), intent(in) :: Cv, gamm
   real(PR), intent(out) :: Cp

   Cp=Cv*gamm

end subroutine compute_Cp2


subroutine compute_Cp(Z,nu,T,Cp)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T
   real(PR), intent(out) :: Cp
   real(PR) :: rho, rho1, rho2, rhomin, rhomax, P, P1, P2, E1, E2,x(1:1)
   
   rho=1.0_PR/nu

   !P=P_tot(Z,rho,T)

   !!!!---find rho1
   !rhomin=rho/10.0_PR ; rhomax=rho*10.0_PR ; x(1)=rho
   !ipar(1)=Z ; rpar(1)=T1 ; rpar(2)=P

   !call My_Newton(1,x,rhomin,rhomax,F_Cp,J_Cp)

   !rho1=x(1) ; P1=P_tot(Z,rho1,T1)
   !if(abs(P-P1)/P.gt.1.0e-4_PR)then
   !  print*, 'PROBLEM Cp 1 : reduce tolerance for Newton'
   !  stop
   !endif

   !!!!---find rho2
   !rhomin=rho/10.0_PR ; rhomax=rho*10.0_PR ; x(1)=rho
   !ipar(1)=Z ; rpar(1)=T2 ; rpar(2)=P

   !call My_Newton(1,x,rhomin,rhomax,F_Cp,J_Cp)

   !rho2=x(1) ; P2=P_tot(Z,rho2,T2)
   !if(abs(P-P2)/P.gt.1.0e-4_PR)then
   !  print*, 'PROBLEM Cp 2 : reduce tolerance for Newton'
   !  stop
   !endif

   !E1=E_tot(Z,rho1,T1)
   !E2=E_tot(Z,rho2,T2)

   !Cp=(E2-E1)/dT

   Cp=ddTp_E_tot(Z,rho,T)

end subroutine compute_Cp



subroutine compute_soundspeed(Z,nu,T,c)
   !!! Formule 5.2-9 dans Hirschfelder & Curtiss p 370
   !!! attention : grandeurs molaires représentées pâr "~"
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T
   real(PR), intent(out) :: c
   real(PR) :: rho, Cv, dPdT, dPdnu

   rho=1.0_PR/nu

   Cv=ddTr_E_tot(Z,rho,T)
   dPdT=ddTr_P_tot(Z,rho,T)
   dPdnu=-ddrT_P_tot(Z,rho,T)/nu**2

   !!!---en m/s 
   if(dPdnu.lt.0.0_PR)then
     c=sqrt(1000.0_PR)*nu*sqrt(T/Cv)*dPdT*sqrt(1.0_PR-Cv/T*dPdnu/dPdT**2)
   else
     c=0.0_PR
   endif
end subroutine compute_soundspeed


subroutine compute_gamma(Z,nu,T,c,g)
   !!! Formule 5.2-9 dans Hirschfelder & Curtiss p 370
   !!! attention : grandeurs molaires représentées pâr "~"
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T, c
   real(PR), intent(out) :: g

   g=1e-3_PR*A(Z)*c**2/(T*Rg)

end subroutine compute_gamma

subroutine compute_gamma2(Z,nu,T,g)
   !!! gam=1/e*(dP/drho)e+1
   !!!    =-1/e*nu**2*(dPdnu)e+1
   !! gam -> pinf -> c
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu,T
   real(PR), intent(out) :: g
   real(PR) :: rho, e

   rho=1.0_PR/nu
   e=E_tot(Z,rho,T)
   g=1/e*ddre_P(Z,rho,T)+1.0_PR

end subroutine compute_gamma2

subroutine compute_pinf(Z,nu,e,g,p,pinf)
   !!! gam=1/e*(dP/drho)e+1
   !!!    =-1/e*nu**2*(dPdnu)e+1
   !! gam -> pinf -> c
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu,e,g,p
   real(PR), intent(out) :: pinf
   real(PR) :: rho

   rho=1.0_PR/nu
   pinf=(rho*e*(g-1.0_PR)-p)/g

end subroutine compute_pinf

subroutine compute_soundspeed2(Z,nu,p,g,pinf,c)
   !!! c=sqrt(gam*(p+pinf)/rho)
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, p, g, pinf
   real(PR), intent(out) :: c
   real(PR) :: c2

   c2=1.0e3_PR*g*(p+pinf)*nu
   
   if(c2.gt.0.0_PR)then
    c=sqrt(c2)
   else
    c=0.0_PR
    print*, 'c2<0:', g, p, pinf
    stop
   endif

end subroutine compute_soundspeed2

subroutine compute_Zbar(Z,nu,T,zbar)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T
   real(PR), intent(out) :: zbar
   real(PR) :: TeV, ndens, rho
   integer, save :: init=0

   TeV=T/eV
   rho=1.0_PR/nu
   ndens=Avogadro*rho/A(Z)

   if(init.eq.0)then
      call init_thomas_fermi_ionization_data()
      init=1
   endif
    
   zbar=thomas_fermi_zbar(ndens,TeV,real(Z,PR))
 
end subroutine compute_Zbar

subroutine compute_Zbar2(Z,nu,T,zbar)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: nu, T
   real(PR), intent(out) :: zbar
   real(PR) :: TeV, ndens, rho,zs,zfs
   integer, save :: init=0

   TeV=T/eV
   rho=1.0_PR/nu
   ndens=Avogadro*rho/A(Z)

   if(init.eq.0)then
      call init_thomas_fermi_ionization_data()
      init=1
   endif
    
   zbar=thomas_fermi_zbar(ndens,TeV,real(Z,PR))

   call compute_ZFS(Z,rho,T,zbar,zs,zfs)

   zbar=zfs
 
end subroutine compute_Zbar2

subroutine F_Cp(N,x,f)

     !use mod_data, only : rpar, ipar

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     real(PR) :: T, P0
     integer :: Z

     Z=ipar(1)
     T=rpar(1)
     P0=rpar(2)

     f(1)=P_tot(Z,x(1),T)-P0
    
end subroutine F_Cp

subroutine J_Cp(N,x,J)

     !use mod_data, only : rpar, ipar

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     real(PR) :: T, P0
     integer :: Z

     Z=ipar(1)
     T=rpar(1)
     P0=rpar(2)
  
     J(1,1)=ddr(Z,x(1),T,P_tot)
    
end subroutine J_Cp

!!!======================  NEWTON  ================

subroutine My_Newton(N,x,xmin,xmax,F,DF)
   !!! Newton une seule variable
   implicit none
   integer, intent(in) :: N
   real(PR), intent(inout) :: x(1:N)
   real(PR), intent(in) :: xmin,xmax 
   interface
       subroutine fun(N,x,f)
          integer, intent(in) :: N
          double precision, intent(in) :: x(1:N)
          double precision, intent(out) :: f(1:N)
       end subroutine fun   
       subroutine Jac(N,x,J)
          integer, intent(in) :: N
          double precision, intent(in) :: x(1:N)
          double precision, intent(out) :: J(1:N,1:N)
       end subroutine Jac
   end interface    
   procedure(fun) :: F !!! fonction
   procedure(Jac) :: DF !!! jacobienne
   real(PR) :: tol, fx(1:1), dfdx(1:1), dx
   integer :: it    

   call F(1,x,fx)
   tol=abs(fx(1))
   it=0

   do while(tol.gt.1.0e-12_PR.and.it.lt.10000)
 
      call F(1,x(1),fx(1:1))
      call DF(1,x(1),dfdx(1:1))

      dx=-fx(1)/dfdx(1)

      x(1)=x(1)+dx
      x(1)=max(x(1),xmin)
      x(1)=min(x(1),xmax)

      it=it+1

      tol=abs(fx(1))

      !print*, 'it=', it, 'x=', x(1), 'fx=', fx(1)

      !if(it.ge.9900)then
      ! print*, 'it=', it, x(1), fx(1), dx, dfdx(1)
      !endif

   enddo

   if(it.ge.10000)then
     print*, 'PROBLEM NEWTON', Temp, x(1), fx(1)
     stop 
   endif

end subroutine My_newton

subroutine foundr1r2r3(Ptarg,rhomin,rhomax,Nr,foundr1,foundr2,foundr3,r1,r2,r3)

   implicit none
   real(PR), intent(in) :: Ptarg,rhomin,rhomax
   integer, intent(in) :: Nr
   logical, intent(out) :: foundr1,foundr2,foundr3
   real(PR), intent(out) :: r1,r2,r3
   real(PR) :: P,rhop,rho,nu
   integer :: ir

   foundr1=.false. ; foundr2=.false. ; foundr3=.false.

   rho=rhomax

   !open(unit=5, file='check.dat', status='replace')

   do ir=1,Nr
    
     rhop=rho
     rho=rhomax*(rhomin/rhomax)**(real(ir,PR)/real(Nr,PR))
     nu=1.0_PR/rho
   
     P=Pnu(nu)

     if(.not.foundr1)then !!! le premier point n'a pas été trouvé

        if(P.lt.Ptarg)then
          r1=0.5_PR*(rhop+rho)  !!!r1=dicho(rhop,rho)
          foundr1=.true.
        elseif(P.eq.Ptarg)then
          r1=rho  
          foundr1=.true.
        endif
     
     elseif(foundr1.and..not.foundr2)then !!! le premier a pas été trouvé mais pas r2

        if(P.gt.Ptarg)then
          r2=0.5_PR*(rhop+rho)  !!!r1=dicho(rhop,rho)
          foundr2=.true.
        elseif(P.eq.Ptarg)then
          r2=rho  
          foundr2=.true.
        endif

     elseif(foundr1.and.foundr2.and..not.foundr3)then !!! r1 et r2 ont été trouvés mais pas r3

        if(P.lt.Ptarg)then
          r3=0.5_PR*(rhop+rho)  !!!r1=dicho(rhop,rho)
          foundr3=.true.
        elseif(P.eq.Ptarg)then
          r3=rho  
          foundr3=.true.
        endif

     endif
   
   !write(5,*) nu, ptarg, P
   
   enddo

   !close(5)

end subroutine foundr1r2r3


subroutine MAXWELL(Pc,rhoc,rhomin,rhomax,Nr,foundP,it,P,r1,r2,r3)

   implicit none
   real(PR), intent(in) :: Pc,rhoc,rhomin,rhomax
   integer, intent(in) :: Nr
   integer, intent(out) :: foundP
   integer, intent(out) :: it
   real(PR), intent(out) :: P,r1,r2,r3
   logical :: foundr123,foundr1,foundr2,foundr3,foundmin
   real(PR) :: tol, Ener, Pmin2, minP, maxP

   foundP=0
   foundr123=.true.
   Pmin2=Pnu(1.0_PR/rhomin)
   it=0

   P=Pnu(1.0_PR/rhoc)

   IF(P.ge.Pc)THEN

    foundP=3
    r1=rhoc
    r2=rhoc
    r3=rhoc

   ELSE

      call minmax(rhomin,rhomax,Nr,foundmin,minP,maxP)

      IF(foundmin)THEN

         !print*, 'foudmin minmax:', minP, maxP

         DO WHILE (foundP.eq.0.and.it.lt.100) !!! P 
  
 
            P=0.5_PR*(minP+maxP)
            
            call foundr1r2r3(P,rhomin,rhomax,Nr,foundr1,foundr2,foundr3,r1,r2,r3)
 
            !print*, foundr1,foundr2, foundr3
 
            if(foundr1.and.foundr2.and.foundr3)then
  
               foundr123=.true. 
           
               nu1=1.0_PR/r1 ; nu2=1.0_PR/r3
               
               Ener=int_nu(Pnu,nu1,nu2)-P*(nu2-nu1)
               tol=1.0e-4_PR*P*(nu2-nu1)
               
               if(abs(Ener).le.tol)then
               
                 foundP=1
               
               elseif(Ener.gt.tol)then
               
                 minP=P
               
               else
               
                 maxP=P

                 if(abs(P-Pmin2)/Pmin2.lt.1.0e-2_PR)then
                   foundP=2
                   P=Pmin2
                   r3=rhomin
                 endif

               endif


            elseif(foundr1.and.foundr2)then

               if(maxP.le.Pmin2)then
                 P=Pmin2
                 foundP=2
                 r3=rhomin
               else
                 minP=Pmin2
               endif

            else
               
               foundP=3

            endif
      
            it=it+1

          ENDDO !!! P

       ELSE

          foundP=3

       ENDIF !! foundmin

   ENDIF !! P > Pc

end subroutine MAXWELL


subroutine minmax(rhomin,rhomax,Nr,foundmin,minP,maxP)

   implicit none
   real(PR), intent(in) :: rhomin,rhomax
   integer, intent(in) :: Nr
   logical, intent(out) :: foundmin
   real(PR), intent(out) :: minP, maxP
   real(PR) :: rho,rho1,rho2,P1,P2,nu1,nu2
   integer :: ir

   foundmin=.false. 
   maxP=-1.0e30_PR

   do ir=1,Nr-1
    
     rho1=rhomax*(rhomin/rhomax)**(real(ir,PR)/real(Nr,PR))
     rho2=rhomax*(rhomin/rhomax)**(real(ir+1,PR)/real(Nr,PR))

     nu1=1.0_PR/rho1
     nu2=1.0_PR/rho2
  
     P1=Pnu(nu1)
     P2=Pnu(nu2)

     if(.not.foundmin.and.(P2.gt.P1))then
       minP=P1
       foundmin=.true.
     endif

     if(foundmin)then
       maxP=max(maxP,P1)
     endif

   enddo

end subroutine minmax


!!!=========================== TABLES =====================================


subroutine create_PE_VT_file(id,Z,N_nu,rho0,rhof,N_T,temperatures)

   implicit none
   integer, intent(in) :: id,N_nu, N_T 
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho0, rhof
   real(PR), intent(in) :: temperatures(1:N_T)
   character(len=5) :: chari

   write(chari,'(I5)') N_T

   open(unit=10, file='tab_'//element_name(Z)//'_'//trim(adjustl(chari))//'.dat', status='replace')

   print*, ' > OUTPUT table :', 'tab_'//element_name(Z)//'_'//trim(adjustl(chari))//'.dat'

   write(id,'(A)') '#=================================================================='
   write(id,'(A)') '# Phase diagram for '//element_name(Z)
   if(mode_output.eq.1)then
     write(id,'(A)') '#-------- mode= 1'
   else
     write(id,'(A)') '#-------- mode= 2' 
   endif
   write(id,'(A4,I3)') '# Z=', Z
   write(id,'(A4,ES14.7)') '# A=', A(Z)
   if(mode_output.eq.1)then
      write(id,'(A7,ES14.7)') '# rhos=', rhos(Z)
      write(id,'(A7,ES14.7)') '# sig0=', sig0(Z)
      write(id,'(A8,ES14.7)') '# lamb0=', lamb0(Z)
      write(id,'(A7,I6)') '# Nnu =', N_nu
      write(id,'(A7,ES14.7)') '# rho0=', rho0
      write(id,'(A7,ES14.7)') '# rhof=', rhof
   else
      write(id,'(A7,ES14.7)') '# rhos=', rhos(Z)*1000.0_PR
      write(id,'(A7,I6)') '# Nnu =', N_nu
      write(id,'(A7,ES14.7)') '# rho0=', rho0*1000.0_PR
      write(id,'(A7,ES14.7)') '# rhof=', rhof*1000.0_PR
   endif

   write(id,'(A7,I6)') '# N_T =', N_T
   write(id,'(A6,ES14.7)') '# T0 =', Temperatures(1)
   write(id,'(A6,ES14.7)') '# Tf =', Temperatures(N_T)

   if(mode_output.eq.1)then
      write(id,'(A)') '# Temperature list:'
      do j=1, N_T
       write(id,'(A2,ES14.7)') '# ', temperatures(j) 
      enddo
   endif

   write(id,'(A)') '#=================================================================='
   write(id,*) ""
   write(id,*) ""

end subroutine create_PE_VT_file





end program main
