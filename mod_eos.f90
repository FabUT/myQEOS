module mod_EOS

implicit none
integer, parameter, public :: PR=selected_real_kind(8)
real(PR), parameter, public :: Pi=2.0_PR*Asin(1.0_PR)

!!!=========== CONSTITUANTS
logical :: ions=.true.
logical :: electrons=.true.
logical :: bounding=.true.

!!!!========== CONSTANTES PHYSIQUES
real(PR), parameter :: Avogadro=6.0221409e23_PR
real(PR), parameter :: kb=1.38064852e-23_PR !!! J K-1
real(PR), parameter :: qe=1.602176e-19_PR !!! C
real(PR), parameter :: eV=qe/kb !!!! K
real(PR), parameter :: uma=1.660538921e-24_PR !!!g
real(PR), parameter :: mp=1.672621898e-24_PR  !!! g
real(PR), parameter :: mn=1.674927471e-24_PR  !!! g
real(PR), parameter :: me=9.10938356e-28_PR !!!g
real(PR), parameter :: meSI=9.10938356e-31_PR !!!g
real(PR), parameter :: Rg=8.31447_PR !!! J mol-1 K-1
real(PR), parameter :: h=6.626e-34_PR !!!! J s
real(PR), parameter :: hb=h/(2.0_PR*Pi) !!!! J s
real(PR), parameter :: clight=299792458.0_PR
real(PR), parameter :: eps0=8.854187817e-12_PR
real(PR), parameter :: RBohr=h**2*eps0/(pi*meSI*qe**2) !!! Bohr radius


!!!!============ MATERIAUX
real(PR), dimension(92) :: A !!! g/mol ou nucléon/atome 
real(PR), dimension(92) :: ge !!! mJ/mol/K
real(PR), dimension(92) :: n,m !!! Constants de Lennard Jones
real(PR), dimension(92) :: Tms !!! Temperature de fusion a la densité de reference K
character(len=2), dimension(92) :: symbol


!!!!============ MATERIAUX
character (len=2) :: element_name(92) !< the element abbreviated name
real(PR) :: cohesive_energy(92) !< The cohesive energy for each element, eV
real(PR) :: melting_point(92)   !< The melting point of each element, eV
real(PR) :: boiling_point(92)   !< The boiling point of each element, eV
real(PR) :: rhos(92)   !< Mass density of the element in the solid phase, g cm^-3
real(PR) :: debye_temp(92)      !< Debeye Temperature at solid density, eV
real(PR) :: gruneisen(92)       !< Grüneisen parameter at solid density
real(PR) :: bulk_modulus(92)    !< Bulk modulus, GPa
real(PR) :: bm_pderiv(92)       !< Pressure derivative of the Bulk Modulus
!!!--critical parameter : just data for post-processing
real(PR) :: Tcrit(92)       !< critical temperature 
real(PR) :: Pcrit(92)       !< critical pressure
real(PR) :: rcrit(92)     !< critical density 
!!!---Fermi energies
real(PR) :: E_Fermi(92)
real(PR) :: v_Fermi(92)
!!!---for Desjarlais Z model
real(PR) :: Eio1(92) ! First ionisation energy (eV)
real(PR) :: pola(92) ! atomic polarizability RBohr**3
real(PR) :: sigSL(92) ! ratio of conductivities sig_solid/sig_liquid at Tm
real(PR) :: sig0(92) ! reference conductivities for solids at T=300 K
real(PR) :: sigmax(92) ! max conductivity before scaling at T=300K andrho=rhosolide 
real(PR) :: lamb0(92) ! reference thermal conductivities for solids at T=300 K
real(PR) :: lambmax(92) ! max thermal conductivity before scaling 

type table_eos
   integer :: Z
   integer :: Nt
   integer :: Nr
   real(PR),allocatable :: P(:,:)
   real(PR),allocatable :: T(:)
   real(PR),allocatable :: rho(:)
   real(PR),allocatable :: Peq(:)
   integer, allocatable :: ir_eb1(:)
   integer, allocatable :: ir_eb2(:)
end type table_eos



real(PR) :: Pg, Tg, tolP,tolR
integer :: Zg

abstract interface
      subroutine func(x,fval,fderiv)
         IMPLICIT NONE
         double precision, INTENT(IN) :: x
         double precision, INTENT(OUT) :: fval,fderiv
      end subroutine func
      !subroutine pression_totale(Z,rho,T)
      !   IMPLICIT NONE
      !   double precision, INTENT(IN) :: Z
      !   double precision, INTENT(OUT) :: rho,T
      !end subroutine pression totale
end interface

!!!!!!procedure (func), pointer :: f_ptr => null ()


contains



!!!======================   CONTRIBUTION DES IONS  ====================

!!!!-------------- Modèle de Johnson
!real(PR) FUNCTION a_J(Z)
!implicit none
!integer, intent(in) :: Z
!
!a_J=1.25_PR/A(Z)**(5.0_PR/3.0_PR)
!
!END FUNCTION a_J
!
!
!real(PR) FUNCTION Tm_J(Z,rho)
!implicit none
!integer, intent(in) :: Z
!real(PR), intent(in) :: rho
!real(PR) :: a_J 
!
!
!Tm=TD(Z,rho)**2/(a_J(Z)*rho**(2.0_PR/3.0_PR))
!
!
!END FUNCTION Tm_J
!
!
!real(PR) FUNCTION TD_J(Z,rho)
!implicit none
!integer, intent(in) :: Z
!real(PR), intent(in) :: rho
!real(PR) :: G,rhoref,Tmref,TDref,C
!
!G=G_S(Z)
!
!rhoref=rhos(Z)
!Tmref=Tm(Z)
!
!TDref=sqrt(a_J(Z)*rhoref**(2.0_PR/3.0_PR)*Tmref)
!C=log(TDref)-G*rhoref
!
!TD=exp(G*rho+C)
!
!
!END FUNCTION TD_J
!
!
!real(PR) FUNCTION Ei_J(Z,rho,T)
!implicit none
!integer, intent(in) :: Z
!real(PR), intent(in) :: rho, T
!real(PR) :: G,rhoref,Tmref,TDref,C,psi
!
!G=G_S(Z)
!
!rhoref=rhos(Z)
!Tmref=Tm(Z)
!
!TDref=sqrt(a_J(Z)*rhoref**(2.0_PR/3.0_PR)*Tmref)
!C=log(TDref)-G*rhoref
!
!TD=exp(G*rho+C)
!
!
!psi=T/Tm(Z,rho)
!
!IF(psi.le.1)THEN !!! standard Debye model
!       
!
!ELSEIF(mode.eq.1)THEN
!    a1=-5.7 !!!+mixing
!    y=(201.0_PR*(1600*a1**2+2398.0_PR*(4.0_PR*a1+5.0_PR))**0.5_PR - 40.0_PR*(5.0_PR-197.0_PR*a1))/(3980.0_PR*(4.0_PR*a1+5.0_PR))
!    a3=200.0_pr
!    a2=3.0_pr/2.0_pr*(1.0_pr+a3)**3/(a3*(1.0_pr-y)*(a3*y+2.0_pr-y))
!    a4=-8.0_pr/5.0_pr*(a1+a2/(1.0_pr+a3))
!    e1=3.0_pr/2.0_pr*(psi-1.0_pr/psi)
!    e2=0.66_pr/psi
!    a1=3.0_pr/2.0_pr*(2.0_pr-psi-1.0_pr/psi)
!    a2=0.66_pr/psi-0.6_pr
!    E0=-3.0_pr/2.0_pr+3.0_pr/2.0_pr*a4*(1.0_pr-1.0_pr/(2.0_pr*psi**0.5_pr))/psi**1.5_pr+a2*(a3*y+psi**(1.0_pr-y))/(psi**y*(a3+psi**(1.0_pr-y))**2)
!    A0=1.5_pr*log(psi)+a1+a2/(psi**y*(a3+psi**(1.0_pr-y)))+a4*(1.0_pr-3.0_pr/(8.0_pr*psi**0.5_pr))/psi**1.5_pr
!
!
!  if(psi.le.1.2_PR)then
!
!    A=AD+NkT*(A0+a1)
!    E=ED+NkT*(E0+e1)
!    P=
!  else
!
!
!  endif
!
!ELSEIF(mode.eq.2)THEN
!
!
!
!
!ENDIF !!! psi
!
!END FUNCTION Ei_J
!
!








!!!---------------  Modèle de Cowan

real(PR) FUNCTION Tm(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho
real(PR) :: b, xi 
integer, save :: Zp=0
real(PR), save :: TmsCowan

b=0.6_PR*real(Z,PR)**(1.0_PR/9.0_PR)
!xi=eta(Z,rho)
xi=rho/(A(Z)/(9.0_PR*real(Z,PR)**0.3_PR))

Tm=eV*0.32_PR*( xi**(2.0_PR*b+10.0_PR/3.0_PR)/(1.0_PR+xi)**4 )

!!!!---correction de la température
if(Z.ne.Zp)then
 xi=rhos(Z)/(A(Z)/(9.0_PR*real(Z,PR)**0.3_PR))
 TmsCowan=eV*0.32_PR*( xi**(2.0_PR*b+10.0_PR/3.0_PR)/(1.0_PR+xi)**4 )
 Zp=Z
 print*, '-------------------------------------'
 print*, 'correction of Cowan temperature:'
 print*, 'Z=', Z
 print*, 'TmsCowan=', Tmscowan
 print*, 'Tms=', Tms(Z)
 print*, 'ratio=', Tms(Z)/TmsCowan
 print*, '-------------------------------------'
endif
Tm=Tm*Tms(Z)/TmsCowan


END FUNCTION Tm


real(PR) FUNCTION TD(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho
real(PR) :: b, xi

b=0.6_PR*real(Z,PR)**(1.0_PR/9.0_PR)
!xi=eta(Z,rho)
xi=rho/(A(Z)/(9.0_PR*real(Z,PR)**0.3_PR))

TD=eV*1.68_PR/(real(Z,PR)+22.0_PR)*( xi**(b+2.0_PR)/(1.0_PR+xi)**2 ) 

END FUNCTION TD

real(PR) FUNCTION G_S(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho
real(PR) :: b, xi 

!!!----Gruneisen parameter of the solid phase

b=0.6_PR*real(Z,PR)**(1.0_PR/9.0_PR)

!xi=eta(Z,rho)
xi=rho/(A(Z)/(9.0_PR*real(Z,PR)**0.3_PR))

G_S=b+2.0_PR/(1.0_PR+xi)

END FUNCTION G_S



real(PR) FUNCTION G_F(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho

!!!----Gruneisen parameter of the fluid phase

G_F=3.0_PR*G_S(Z,rho)-1.0_PR


END FUNCTION G_F


real(PR) FUNCTION E_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w

u=TD(Z,rho)/T
w=Tm(Z,rho)/T

if( w .gt. 1.0_PR )then  !!!---Solid phase
  !print*, 'coucouz1'

  E_i=3.0_PR*(kb*T/(A(Z)*uma))*( 1.0_PR + u**2.0_PR/20.0_PR - u**4.0_PR/1680.0_PR )
  
elseif( w .le. 1.0_PR )then !!!---Fluide phase
  !print*, 'coucouz2', w

  E_i=3.0_PR/2.0_PR*(kb*T/(A(Z)*uma))*( 1.0_PR + w**(1.0_PR/3.0_PR) )

endif

END FUNCTION E_i

real(PR) FUNCTION ddTr_E_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w,u1,u2,w1,w2,T1,T2,E_i1,E_i2,dT

w=Tm(Z,rho)/T
dT=0.01_PR*T

if( w .gt. 1.0_PR )then  !!!---Solid phase

   T1=T-0.5_PR*dT
   T2=min(T+0.5_PR*dt,Tm(Z,rho)) 

   u1=TD(Z,rho)/T1
   u2=TD(Z,rho)/T2

   E_i1=3.0_PR*(kb*T1/(A(Z)*uma))*( 1.0_PR + u1**2/20.0_PR - u1**4/1680.0_PR )
   E_i2=3.0_PR*(kb*T2/(A(Z)*uma))*( 1.0_PR + u2**2/20.0_PR - u2**4/1680.0_PR )

elseif( w .le. 1.0_PR )then !!!---Fluide phase

   T1=max(T-0.5_PR*dT,Tm(Z,rho))
   T2=T+0.5_PR*dT

   w1=Tm(Z,rho)/T1
   w2=Tm(Z,rho)/T2

  E_i1=3.0_PR/2.0_PR*(kb*T1/(A(Z)*uma))*( 1.0_PR + w1**(1.0_PR/3.0_PR) )
  E_i2=3.0_PR/2.0_PR*(kb*T2/(A(Z)*uma))*( 1.0_PR + w2**(1.0_PR/3.0_PR) )

endif

ddTr_E_i=(E_i2-E_i1)/(T2-T1)

END FUNCTION ddTr_E_i


real(PR) FUNCTION ddTp_E_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w,u1,u2,w1,w2,T1,T2,E_i1,E_i2,P,rho1,rho2,rhomin,rhomax,dT

P=P_tot(Z,rho,T)

w=Tm(Z,rho)/T

dT=0.01_PR*T

if( w .gt. 1.0_PR )then  !!!---Solid phase

   !T1=T-dT ; rho1=rho_TP(Z,T1,P,rho)
   !T2=T ; rho2=rho

    T1=T-0.5_PR*dT ; rho1=rho_TP(Z,T1,P,rho)
    T2=T+0.5_PR*dT ; rho2=rho_TP(Z,T2,P,rho)


   !!! T1 < T : il faut augmenter rho -> Tm augmente : on a toujours w > 1 

   u1=TD(Z,rho1)/T1
   u2=TD(Z,rho2)/T2

   E_i1=3.0_PR*(kb*T1/(A(Z)*uma))*( 1.0_PR + u1**2/20.0_PR - u1**4/1680.0_PR )
   E_i2=3.0_PR*(kb*T2/(A(Z)*uma))*( 1.0_PR + u2**2/20.0_PR - u2**4/1680.0_PR )

elseif( w .le. 1.0_PR )then !!!---Fluide phase

   !T1=T ; rho1=rho
   !T2=T+dT ; rho2=rho_TP(Z,T2,P,rho)
    T1=T-0.5_PR*dT ; rho1=rho_TP(Z,T1,P,rho)
    T2=T+0.5_PR*dT ; rho2=rho_TP(Z,T2,P,rho)

   w1=Tm(Z,rho1)/T1
   w2=Tm(Z,rho2)/T2

  E_i1=3.0_PR/2.0_PR*(kb*T1/(A(Z)*uma))*( 1.0_PR + w1**(1.0_PR/3.0_PR) )
  E_i2=3.0_PR/2.0_PR*(kb*T2/(A(Z)*uma))*( 1.0_PR + w2**(1.0_PR/3.0_PR) )

endif

ddTp_E_i=(E_i2-E_i1)/(T2-T1)

END FUNCTION ddTp_E_i

!real(PR) FUNCTION ddTr_E_i(Z,rho,T)
!implicit none
!integer, intent(in) :: Z
!real(PR), intent(in) :: rho,T
!real(PR) :: u,w
!
!
!u=TD(Z,rho)/T
!w=Tm(Z,rho)/T
!
!if( w .gt. 1.0_PR )then  !!!---Solid phase
!
!  ddTr_E_i=3.0_PR*(kb/(A(Z)*uma))*( 1.0_PR + u**2/20.0_PR - u**4/1680.0_PR )+&
!           3.0_PR*(kb/(A(Z)*uma))*( -u**2/(10.0_PR) + u**4/(420.0_PR) )
!  
!elseif( w .le. 1.0_PR )then !!!---Fluide phase
!
!  ddTr_E_i=3.0_PR/2.0_PR*(kb/(A(Z)*uma))*( 1.0_PR + w**(1.0_PR/3.0_PR) )-&
!      1.0_PR/2.0_PR*(kb/(A(Z)*uma))*w**(1.0_PR/3.0_PR) 
!
!endif
!
!END FUNCTION ddTr_E_i

real(PR) FUNCTION solid(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w

u=TD(Z,rho)/T
w=Tm(Z,rho)/T

if( w .gt. 1.0_PR )then  !!!---Solid phase

   solid=1.0_PR
  
elseif( w .le. 1.0_PR )then !!!---Fluide phase

   solid=0.0_PR

endif

END FUNCTION solid

real(PR) FUNCTION P_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w


u=TD(Z,rho)/T
w=Tm(Z,rho)/T

if( w .gt. 1.0_PR )then  !!!---Solid phase

  P_i=G_S(Z,rho)*rho*E_i(Z,rho,T)

elseif( w .le. 1.0_PR )then !!!---Fluide phase

 P_i=rho*kb*T/(A(Z)*uma)*( 1.0_PR + G_F(Z,rho)*w**(1.0_PR/3.0_PR) )

endif

END FUNCTION P_i

real(PR) FUNCTION S_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: u,w


u=TD(Z,rho)/T
w=Tm(Z,rho)/T

if( w .gt. 1.0_PR )then  !!!---Solid phase

  S_i=kb/(A(Z)*uma)*( 4.0_PR - 3.0_PR*log(u) + ( 3.0_PR*u**2/40.0_PR - 3.0_PR*u**4/2240.0_PR ) )

elseif( w .le. 1.0_PR )then !!!---Fluide phase

  S_i=kb/(A(Z)*uma)*( 7.0_PR - 3.0_PR*w**(1.0_PR/3.0_PR) + 3.0_PR/2.0_PR*log(w/u**2) )

endif

END FUNCTION S_i

real(PR) FUNCTION H_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

H_i=E_i(Z,rho,T)+P_i(Z,rho,T)/rho

END FUNCTION H_i

real(PR) FUNCTION G_i(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

 G_i=H_i(Z,rho,T)-T*S_i(Z,rho,T)

END FUNCTION G_i

real(PR) FUNCTION ddTr_P_i(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: u,w,u1,w1,u2,w2,P1,P2,dT,T1,T2
   
   dT=0.01_PR*T
   T1=T-0.5_PR*dT
   T2=T+0.5_PR*dT

   u=TD(Z,rho)/T
   w=Tm(Z,rho)/T

   u1=TD(Z,rho)/T1 ; u2=TD(Z,rho)/T2
   w1=Tm(Z,rho)/T1 ; w2=Tm(Z,rho)/T2
  
   if( w .gt. 1.0_PR )then  !!!---Solid phase

   
     P1=G_S(Z,rho)*rho*E_i(Z,rho,T1)
     P2=G_S(Z,rho)*rho*E_i(Z,rho,T2)
  
   elseif( w .le. 1.0_PR )then !!!---Fluide phase
   
     P1=rho*kb*T1/(A(Z)*uma)*( 1.0_PR + G_F(Z,rho)*w1**(1.0_PR/3.0_PR) )
     P2=rho*kb*T2/(A(Z)*uma)*( 1.0_PR + G_F(Z,rho)*w2**(1.0_PR/3.0_PR) )
  
   endif
   
   ddTr_P_i=(P2-P1)/(T2-T1)   

END FUNCTION ddTr_P_i

real(PR) FUNCTION ddrT_P_i(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: u,w,u1,w1,u2,w2,P1,P2,dr,rho1,rho2
   
   dr=0.01_PR*rho
   rho1=rho-0.5_PR*dr
   rho2=rho+0.5_PR*dr

   u=TD(Z,rho)/T
   w=Tm(Z,rho)/T

   u1=TD(Z,rho1)/T ; u2=TD(Z,rho2)/T
   w1=Tm(Z,rho1)/T ; w2=Tm(Z,rho2)/T
  
   if( w .gt. 1.0_PR )then  !!!---Solid phase

   
     P1=G_S(Z,rho1)*rho1*E_i(Z,rho1,T)
     P2=G_S(Z,rho2)*rho2*E_i(Z,rho2,T)
  
   elseif( w .le. 1.0_PR )then !!!---Fluide phase
   
     P1=rho1*kb*T/(A(Z)*uma)*( 1.0_PR + G_F(Z,rho1)*w1**(1.0_PR/3.0_PR) )
     P2=rho2*kb*T/(A(Z)*uma)*( 1.0_PR + G_F(Z,rho2)*w2**(1.0_PR/3.0_PR) )
  
   endif
   
   ddrT_P_i=(P2-P1)/(rho2-rho1)   

END FUNCTION ddrT_P_i


real(PR) FUNCTION P_gp(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

!!!---Pression de gaz parfait en MPa
P_gp=rho*Rg/A(Z)*T

END FUNCTION P_gp


!!!======================   CONTRIBUTION DES ELECTRONS  ================

real(PR) FUNCTION E_e(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: Taue,X,r,beta0,beta,TeV

!!!---- McClosekey's model

!TeV=T !!!/eV

Taue=0.6_PR

X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV

r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)

beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K^2 ->  J/g/K^2

beta=beta0*eta(Z,rho)**(-Taue)

if(beta*T/r.lt.1.0e2_PR)then
  E_e=9.0_PR*r**2/(4.0_PR*beta)*log(cosh( (2.0_PR*beta*T)/(3.0_PR*r) ))
  !print*, 'coucou1', (2.0_PR*beta*T)/(3.0_PR*r)
else
  E_e=3.0_PR/2.0_PR*r*T
  !print*, 'coucou2', (2.0_PR*beta*T)/(3.0_PR*r)
endif

END FUNCTION E_e

real(PR) FUNCTION ddTr_E_e(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: Taue,X,r,beta0,beta,TeV,dT
real(PR) :: T1, T2, E_e1, E_e2,X1,X2,r1,r2

!!!---- McClosekey's model

Taue=0.6_PR

X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV

r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)

beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K^2 ->  J/g/K^2

beta=beta0*eta(Z,rho)**(-Taue)

dT=0.01_PR*T
T1=T-0.5_PR*dT
T2=T+0.5_PR*dT

X1=real(Z,PR)**(-4.0_PR/3.0_PR)*T1/eV
X2=real(Z,PR)**(-4.0_PR/3.0_PR)*T2/eV

r1=0.85_PR*X1**0.59_PR/( 1.0_PR + 0.85_PR*X1**0.59_PR )*real(Z,PR)*Rg/A(Z)
r2=0.85_PR*X2**0.59_PR/( 1.0_PR + 0.85_PR*X2**0.59_PR )*real(Z,PR)*Rg/A(Z)

if(beta*T/r.lt.1.e3_PR)then

   if(beta*T2/r2.ge.1.e3_PR)then
      T2=T ; r2=r ; X2=X
   endif

   E_e1=9.0_PR*r1**2/(4.0_PR*beta)*log(cosh( (2.0_PR*beta*T1)/(3.0_PR*r1) ))
   E_e2=9.0_PR*r2**2/(4.0_PR*beta)*log(cosh( (2.0_PR*beta*T2)/(3.0_PR*r2) ))

else

  if(beta*T1/r1.lt.1.e3_PR)then
     T1=T ; r1=r ; X1=X
  endif

  E_e1=3.0_PR/2.0_PR*r1*T1
  E_e2=3.0_PR/2.0_PR*r2*T2

endif

ddTr_E_e=(E_e2-E_e1)/(T2-T1)

END FUNCTION ddTr_E_e

real(PR) FUNCTION ddTp_E_e(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: Taue,X,r,beta0,beta,TeV
real(PR) :: T1, T2, E_e1,E_e2,X1,X2,r1,r2,P,rho1,rho2,beta1,beta2,rhomin,rhomax,dPdr,dT

!!!---- McClosekey's model

!TeV=T !!!/eV

Taue=0.6_PR

X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV

r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)

beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K^2 ->  J/g/K^2

beta=beta0*eta(Z,rho)**(-Taue)

P=P_tot(Z,rho,T)

if(beta*T/r.lt.1.e3_PR)then

  dT=0.01_PR*T
  T1=T-0.5_PR*dT
  T2=T+0.5_PR*dT

  X1=real(Z,PR)**(-4.0_PR/3.0_PR)*T1/eV
  X2=real(Z,PR)**(-4.0_PR/3.0_PR)*T2/eV
  
  r1=0.85_PR*X1**0.59_PR/( 1.0_PR + 0.85_PR*X1**0.59_PR )*real(Z,PR)*Rg/A(Z)
  r2=0.85_PR*X2**0.59_PR/( 1.0_PR + 0.85_PR*X2**0.59_PR )*real(Z,PR)*Rg/A(Z)


  rho1=rho_TP(Z,T1,P,rho)
  rho2=rho_TP(Z,T2,P,rho)
  beta1=beta0*eta(Z,rho1)**(-Taue)
  beta2=beta0*eta(Z,rho2)**(-Taue)

  E_e1=9.0_PR*r1**2/(4.0_PR*beta1)*log(cosh( (2.0_PR*beta1*T1)/(3.0_PR*r1) ))
  E_e2=9.0_PR*r2**2/(4.0_PR*beta2)*log(cosh( (2.0_PR*beta2*T2)/(3.0_PR*r2) ))

else

  dT=0.02_PR*T
  T1=T-0.5_PR*dT
  T2=T+0.5_PR*dT
  
  X1=real(Z,PR)**(-4.0_PR/3.0_PR)*T1/eV
  X2=real(Z,PR)**(-4.0_PR/3.0_PR)*T2/eV
  
  r1=0.85_PR*X1**0.59_PR/( 1.0_PR + 0.85_PR*X1**0.59_PR )*real(Z,PR)*Rg/A(Z)
  r2=0.85_PR*X2**0.59_PR/( 1.0_PR + 0.85_PR*X2**0.59_PR )*real(Z,PR)*Rg/A(Z)

  rho1=rho_TP(Z,T1,P,rho)
  rho2=rho_TP(Z,T2,P,rho)
  beta1=beta0*eta(Z,rho1)**(-Taue)
  beta2=beta0*eta(Z,rho2)**(-Taue)

  !if(beta1*T1/r1.ge.1.0e-3_PR.and.&
  !   beta2*T2/r2.ge.1.0e-3_PR)then

     E_e1=3.0_PR/2.0_PR*r1*T1
     E_e2=3.0_PR/2.0_PR*r2*T2

  !elseif !!!! 

  !endif
  
endif

ddTp_E_e=(E_e2-E_e1)/(T2-T1)

END FUNCTION ddTp_E_e



real(PR) FUNCTION P_e(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: Taue,X,r,beta0,beta,TeV

!!!---- McClosekey's model

!TeV=T !!!!/eV

Taue=0.6_PR

X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV 

r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)

beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K ->  J/g/K

beta=beta0*eta(Z,rho)**(-Taue)


if(beta*T/r.lt.1.e3_PR)then
  P_e=rho*r**2/(Taue*beta)*log(cosh( Taue*beta*T/r ))
else
  P_e=rho*r*T
endif


!P_e=rho**(1.0_PR+Taue)*rhos(Z)**(-Taue)*r**2/(Taue*beta0)*log(cosh( Taue*beta*T/r ))

END FUNCTION P_e

real(PR) FUNCTION ddTr_P_e(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: X1,X2,r1,r2,P1,P2,dT,T1,T2
   real(PR) :: Taue,X,r,beta0,beta,TeV

   Taue=0.6_PR
   
   X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV 
   
   r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)
   
   beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K ->  J/g/K
   
   beta=beta0*eta(Z,rho)**(-Taue)

   dT=0.01_PR*T
   T1=T-0.5_PR*dT
   T2=T+0.5_PR*dT
   
   X1=real(Z,PR)**(-4.0_PR/3.0_PR)*T1/eV
   X2=real(Z,PR)**(-4.0_PR/3.0_PR)*T2/eV
   
   r1=0.85_PR*X1**0.59_PR/( 1.0_PR + 0.85_PR*X1**0.59_PR )*real(Z,PR)*Rg/A(Z)
   r2=0.85_PR*X2**0.59_PR/( 1.0_PR + 0.85_PR*X2**0.59_PR )*real(Z,PR)*Rg/A(Z)
   
   if(beta*T/r.lt.1.e3_PR)then
     P1=rho*r1**2/(Taue*beta)*log(cosh( Taue*beta*T1/r1 ))
     P2=rho*r2**2/(Taue*beta)*log(cosh( Taue*beta*T2/r2 ))
   else
     P1=rho*r1*T1
     P2=rho*r2*T2
   endif
  
   ddTr_P_e=(P2-P1)/(T2-T1)
 
end function ddTr_P_e

real(PR) FUNCTION ddrT_P_e(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: P1,P2,dr,rho1,rho2,beta1,beta2
   real(PR) :: Taue,X,r,beta0,beta,TeV

   Taue=0.6_PR
   
   X=real(Z,PR)**(-4.0_PR/3.0_PR)*T/eV 
   
   r=0.85_PR*X**0.59_PR/( 1.0_PR + 0.85_PR*X**0.59_PR )*real(Z,PR)*Rg/A(Z)
   
   beta0=ge(Z)*1.0e-3_PR/A(Z) !!! ge conversion mJ/mol/K ->  J/g/K
   
   beta=beta0*eta(Z,rho)**(-Taue)

   dr=0.01_PR*rho
   rho1=rho-0.5_PR*dr
   rho2=rho+0.5_PR*dr
   
   if(beta*T/r.lt.1.e3_PR)then
     beta1=beta0*eta(Z,rho1)**(-Taue)
     beta2=beta0*eta(Z,rho2)**(-Taue)
     P1=rho1*r**2/(Taue*beta1)*log(cosh( Taue*beta1*T/r ))
     P2=rho2*r**2/(Taue*beta2)*log(cosh( Taue*beta2*T/r ))
   else
     P1=rho1*r*T
     P2=rho2*r*T
   endif
  
   ddrT_P_e=(P2-P1)/(rho2-rho1)
 
end function ddrT_P_e

real(PR) FUNCTION ddrT_P(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: e,drho,P1,P2,rho1,rho2,T1,T2


   drho=rho*0.01_PR

   rho1=rho
   rho2=rho+drho

   P1=P_tot(Z,rho1,T)
   P2=P_tot(Z,rho2,T)

   ddrT_P=(P2-P1)/(rho2-rho1)   

end function ddrT_P

real(PR) FUNCTION ddre_P(Z,rho,T)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho,T
   real(PR) :: u,drho,P1,P2,rho1,rho2,T1,T2


   drho=rho*0.01_PR
   u=E_tot(Z,rho,T)

   rho1=rho
   rho2=rho+drho

   T1=T   
   T2=T_from_ru(Z,T,rho2,u)

   P1=P_tot(Z,rho1,T1) 
   P2=P_tot(Z,rho2,T2)

   ddre_P=(P2-P1)/(rho2-rho1)   

end function ddre_P

real(PR) function T_from_ru(Z,T0,rho,u)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: T0,rho, u
   real(PR) :: Tmin, Tmax, T, ut, tol
   integer :: it
   logical :: found, problem
   
   Tmin=1.0_PR
   Tmax=1000000000.0_PR

   found=.false.
   problem=.false.
   it=0

   do while(.not.found.and..not.problem)

     T=0.5_PR*(Tmin+Tmax)
  
     ut=E_tot(Z,rho,T)

     if(ut.lt.u)then    
        Tmin=T 
     else
        Tmax=T
     endif

     tol=2.0_PR*abs((ut-u)/(ut+u))

     if(tol.lt.1.0e-4_PR) found=.true. 

     if(it.gt.97)then
      print*, rho, T, ut, u, tol, Tmin, Tmax
     endif

     if(it.gt.100) problem=.true.

     it=it+1
 
   enddo

   if(problem)then
      !print*, 'PROBLEM T_from_ru'
      !stop
      print*, ' warning: it=100 : tol=', tol  
      T_from_ru=T
   elseif(found)then
      T_from_ru=T
   endif


end function T_from_ru




real(PR) FUNCTION H_e(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

H_e=E_e(Z,rho,T)+P_e(Z,rho,T)/rho

END FUNCTION H_e




!!!======================   ENERGIE DE LIAISON / COURBE FROIDE  ================


real(PR) FUNCTION E_b(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho
real(PR) :: A, B, inv

!!!---- Al'tshuler Lennard Jones model 1980
inv=Ecoh_J(Z)/(m(Z)-n(Z))

A=n(Z)*inv
B=m(Z)*inv

E_b = A*eta(Z,rho)**m(Z) - B*eta(Z,rho)**n(Z) + Ecoh_J(Z)


END FUNCTION E_b

real(PR) FUNCTION P_b(Z,rho)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho
   real(PR) :: A, B, inv
   
   !!!---- Al'tshuler Lennard Jones model 1980
   inv=Ecoh_J(Z)/(m(Z)-n(Z))
   
   A=n(Z)*inv
   B=m(Z)*inv
   
   P_b=rhos(Z)*( A*m(Z)*eta(Z,rho)**(m(Z)+1.0_PR) - B*n(Z)*eta(Z,rho)**(n(Z)+1.0_PR) )
   
END FUNCTION P_b

real(PR) FUNCTION ddrT_P_b(Z,rho)

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho
   real(PR) :: A, B, inv, P1, P2, dr, rho1, rho2
   
   !!!---- Al'tshuler Lennard Jones model 1980
   inv=Ecoh_J(Z)/(m(Z)-n(Z))
   
   A=n(Z)*inv
   B=m(Z)*inv
   
   dr=0.001_PR*rho
   rho1=rho-0.5_PR*dr
   rho2=rho+0.5_PR*dr
   
   P1=rhos(Z)*( A*m(Z)*eta(Z,rho1)**(m(Z)+1.0_PR) - B*n(Z)*eta(Z,rho1)**(n(Z)+1.0_PR) )
   P2=rhos(Z)*( A*m(Z)*eta(Z,rho2)**(m(Z)+1.0_PR) - B*n(Z)*eta(Z,rho2)**(n(Z)+1.0_PR) )

   ddrT_P_b=(P2-P1)/(rho2-rho1)   

END FUNCTION ddrT_P_b





real(PR) FUNCTION H_b(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

H_b=E_b(Z,rho)+P_b(Z,rho)/rho

END FUNCTION H_b

!!!======================  GRANDEURS TOTALES  ================


real(PR) FUNCTION E_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv

E_tot=0.0_PR
if(bounding) E_tot=E_tot+E_b(Z,rho)
if(ions) E_tot=E_tot+E_i(Z,rho,T)
if(electrons) E_tot=E_tot+E_e(Z,rho,T)

END FUNCTION E_tot

real(PR) FUNCTION ddTr_E_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv

ddTr_E_tot=0.0_PR
if(ions) ddTr_E_tot=ddTr_E_tot+ddTr_E_i(Z,rho,T)
if(electrons) ddTr_E_tot=ddTr_E_tot+ddTr_E_e(Z,rho,T)

END FUNCTION ddTr_E_tot

real(PR) FUNCTION ddTp_E_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv, dT

ddTp_E_tot=0.0_PR
if(ions) ddTp_E_tot=ddTp_E_tot+ddTp_E_i(Z,rho,T)
if(electrons) ddTp_E_tot=ddTp_E_tot+ddTp_E_e(Z,rho,T)

END FUNCTION ddTp_E_tot

real(PR) FUNCTION P_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv

P_tot=0.0_PR
if(bounding)  P_tot=P_tot+P_b(Z,rho)
if(ions)      P_tot=P_tot+P_i(Z,rho,T)
if(electrons) P_tot=P_tot+P_e(Z,rho,T)

END FUNCTION P_tot

real(PR) FUNCTION ddTr_P_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv

ddTr_P_tot=0.0_PR
if(ions)      ddTr_P_tot=ddTr_P_tot+ddTr_P_i(Z,rho,T)
if(electrons) ddTr_P_tot=ddTr_P_tot+ddTr_P_e(Z,rho,T)

END FUNCTION ddTr_P_tot

real(PR) FUNCTION ddrT_P_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T
real(PR) :: A, B, inv

ddrT_P_tot=0.0_PR
if(bounding)  ddrT_P_tot=ddrT_P_tot+ddrT_P_b(Z,rho)
if(ions)      ddrT_P_tot=ddrT_P_tot+ddrT_P_i(Z,rho,T)
if(electrons) ddrT_P_tot=ddrT_P_tot+ddrT_P_e(Z,rho,T)

END FUNCTION ddrT_P_tot

real(PR) FUNCTION H_tot(Z,rho,T)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho,T

H_tot=E_tot(Z,rho,T)+P_tot(Z,rho,T)/rho

END FUNCTION H_tot



!!!======================   AUTRES FONCTIONS  ================


real(PR) FUNCTION Ecoh_J(Z)
implicit none
integer, intent(in) :: Z

!Ecoh_J=Ecoh(Z)*qe*rhos(Z)/(A(Z)*uma)

!Ecoh_J=cohesive_energy(Z)*qe*rhos(Z)/(A(Z)*uma) J/cm^3
Ecoh_J=cohesive_energy(Z)*qe/(A(Z)*uma) !!! J/g

END FUNCTION Ecoh_J


real(PR) FUNCTION eta(Z,rho)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: rho

eta=rho/rhos(Z)

END FUNCTION eta


SUBROUTINE GET_R(Z,P_in,T_in,Rho_out)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: P_in, T_in
real(PR), intent(out) :: Rho_out
real(PR) :: rhomin,rhomax,P,rho
logical :: find
integer :: it,itmax

rhomin=1.0e-6_PR*rhos(Z)
rhomax=rhos(Z)*5.0_PR


Tg=T_in
Zg=Z
!tolP=0.001_PR*P_in
tolR=0.0001_PR
Pg=P_in


!!!--Newton
!Rho_out=rtnewt(P_tot_r,rhomin,rhomax,tolR)

!!!---Dichotomie sur rho
tolP=0.005_PR*P_in ; itmax=1000 ; find=.false. ; it=0
do while (.not.find.and.it.lt.itmax)
 rho=0.5_PR*(rhomin+rhomax)
 P=P_tot(Z,rho,T_in)
 if(abs(P-P_in).lt.tolP)then
   find=.true.
   rho_out=rho
 else
   if(P.lt.P_in)then
     rhomin=rho
   else 
     rhomax=rho
   endif
 endif
 it=it+1
enddo

!print*, 'GET_R: rho=', rho_out, 'P=', P_tot(Z,rho_out,T_in), 'vs', P_in 

END SUBROUTINE GET_R

SUBROUTINE P_tot_r(x,fval,fderiv)
IMPLICIT NONE
double precision, INTENT(IN) :: x
double precision, INTENT(OUT) :: fval,fderiv
real(PR) :: P, dx
logical :: find

fval=P_tot(Zg,x,Tg)

dx=tolR/2.0_PR

P=P_tot(Zg,x+dx,Tg)

!dx=0.01_PR*x
!find=.false.
!do while(.not.find)
! P=P_tot(Zg,x+dx,Tg)
! if(abs(P-fval).gt.tolP)then
!   dx=dx*0.5_PR
! else
!   find=.true.
! endif
!enddo

fderiv=( P - fval )/dx
fval=fval-Pg

END SUBROUTINE P_tot_r







SUBROUTINE GET_RT(Z,P_in,Rho_out,T_out)
implicit none
integer, intent(in) :: Z
real(PR), intent(in) :: P_in
real(PR), intent(out) :: Rho_out,T_out
real(PR) :: Tmin, Tmax, rhomin,rhomax,tol,dP,P_c,rho_c,T_c,rho,P,T
integer :: it, ir, ir_c, Nr, itmax
logical :: find

Tmin=0.0_PR
Tmax=100000.0_PR
rhomin=1.0e-4_PR*rhos(Z)
rhomax=rhos(Z)*10.0_PR



!!!---Dichotomie sur T

rho=rhomin ; dP=1.0e30_PR

it=0 ; ir_c=0 ; find=.false. ; itmax=1000

tol=0.001_PR*P_in

do while (.not.find.and.it.lt.itmax)

 T=0.5_PR*(Tmin+Tmax)
 
 do ir=1,Nr+1

  rho=rhomin+real(ir-1,PR)*(rhomax-rhomin)/real(Nr,PR)

  P=P_b(Z,rho)+P_i(Z,rho,T)+P_e(Z,rho,T)

  if(abs(P-P_in).lt.dP)then
     dP=abs(P-P_in)
     ir_c=ir
     P_c=P
     rho_c=rho
     T_c=T
  endif
 

 enddo

 if(dP.lt.tol)then
    find=.true.
    T_out=T_c
    rho_out=rho_c
 else

    if(P_c.lt.P_in)then
      Tmin=T
    else 
      Tmax=T
    endif
   
 endif

 it=it+1

enddo


!print*, 'GET_RT: it=',it, 'T=', T_out, 'rho=', rho_out, 'P=', P_c, 'Pin=', P_in 

END SUBROUTINE GET_RT

!!!========================================= MATH

REAL(PR) FUNCTION rtnewt(funcd,x1,x2,xacc) 

  !Using the Newton-Raphson method, find the root of a function known to lie in the interval
  ![x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
  !is a user-supplied subroutine that returns both the function value and the first derivative of
  !the function.
  !Parameter: MAXIT is the maximum number of iterations.

  IMPLICIT NONE
  procedure(func) :: funcd
  REAL(PR), INTENT(IN) :: x1,x2,xacc
  INTEGER, PARAMETER :: MAXIT=100
  INTEGER :: j
  REAL(PR) :: df,dx,f

    rtnewt=0.5_pr*(x1+x2) !!! Initial guess.
    do j=1,MAXIT
    call funcd(rtnewt,f,df)
    dx=f/df
    rtnewt=rtnewt-dx
    if ((x1-rtnewt)*(rtnewt-x2) < 0.0)then
          print*, 'rtnewt: values jumped out of brackets'
          RETURN
    endif
    if (abs(dx) < xacc)then
            print*, 'convergence j=',j, dx, xacc 
            RETURN  !!!!Convergence.
    endif
    end do
    print*, 'rtnewt exceeded maximum iterations'

END FUNCTION rtnewt


!
!SUBROUTINE rtnewt2(funcd,x1,x2,y1,y2,xacc,v1,v2) 
!
!  !Using the Newton-Raphson method, find the root of a function known to lie in the interval
!  ![x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
!  !is a user-supplied subroutine that returns both the function value and the first derivative of
!  !the function.
!  !Parameter: MAXIT is the maximum number of iterations.
!
!  IMPLICIT NONE
!  REAL(PR), INTENT(IN) :: x1,x2,x3,x4,xacc
!  REAL(PR), INTENT(out) :: v1,v2
!  INTERFACE
!    SUBROUTINE funcd(x,y,fval,df1,df2)
!    IMPLICIT NONE
!    double precision, INTENT(IN) :: x
!    double precision, INTENT(OUT) :: fval,df1,df2
!    END SUBROUTINE funcd
!  END INTERFACE
!  INTEGER, PARAMETER :: MAXIT=200
!  INTEGER :: j
!  REAL(PR) :: invJacobi(2:2)
!  REAL(PR) :: df,dx,f
!
!    v1=0.5_pr*(x1+x2) !!! Initial guess.
!    v2=0.5_pr*(y1+y2) !!! Initial guess.
!
!    do j=1,MAXIT
!    call funcd(v1,v2,f,df1,df2)
!
!    Jacobi
!
!    dx=f*invJacobi
!    rtnewt=rtnewt-dx
!    if ((x1-rtnewt)*(rtnewt-x2) < 0.0)&
!          print*, 'rtnewt: values jumped out of brackets'
!          RETURN
!    if (abs(dx) < xacc) RETURN  !!!!Convergence.
!    end do
!    print*, 'rtnewt exceeded maximum iterations'
!
!END SUBROUTINE rtnewt2
!


!!!========================================   MAXWELL CONSTRUCTION  =============================================










subroutine Maxwell_construct(t)
implicit none
type(table_eos), intent(inout) :: t
real(PR) :: rhoPmax,rhoPmin,rhoPmax2,Pmin,Pmax
real(PR) :: P,P1,P2,tol, intPp, intPm
logical :: find, findPmax, findPmin, findPmax2, find_eb1, find_eb2
integer :: i, ir, it , itermax,iter

print*, '' 
print*, '=============================================='
print*, '' 
print*, '          MAXWELL CONSTRUCT START !!!         '
print*, '' 
print*, '=============================================='
print*, '' 


allocate(t%ir_eb1(1:t%Nt))
allocate(t%ir_eb2(1:t%Nt))
allocate(t%Peq(1:t%Nt))



!DO IT=8,8   !!!1,t%Nt


DO IT=1,t%Nt

!!!--calcul energie totale



!IF(Enertot.gt.0.0_PR)THEN

   print*, 'Maxwell it=', it, '/', t%Nt, 'T=', t%T(it)

   findPmin=.false.
   findPmax=.false.
   findPmax2=.false.


   do ir=1,t%Nr-1
       
     if(t%P(ir+1,it).lt.t%P(ir,it).and..not.findPmax)then

      Pmax=t%P(ir,it)
      rhoPmax=t%rho(ir)  
      findPmax=.true.    
      print*, 'Pmax=', Pmax, 'rho=', 0.32_PR/rhoPmax, 'ir=', ir

     elseif(findPmax.and.((t%P(ir+1,it).gt.t%P(ir,it)).and..not.findPmin))then
      
      Pmin=t%P(ir,it)
      rhoPmin=t%rho(ir)  
      findPmin=.true.    
      print*, 'Pmin=', Pmin, 'rho=', 0.32_PR/rhoPmin, 'ir=', ir
    
      if(Pmin.lt.0.0_PR)then
        print*, 'Pmin<0 !!!'
      endif

     elseif(findPmax.and.findPmin.and.&
             t%P(ir,it).le.Pmax.and.t%P(ir+1,it).ge.Pmax)then

       rhoPmax2=t%rho(ir)
       findPmax2=.true.
       print*, 'rhoPmax2=', 0.32_PR/rhoPmax2

     endif

   enddo


!!!============= Calcul de Peq


if(findPmax.and.findPmax2.and.findPmin)then

   print*, '         Calcul de Peq',findpmin,findPmax,findPmax2

   !!!---Dichotomie sur rho
   print*, '         Dichotomie start !'
   tol=0.01_PR ; itermax=100 ; find=.false. ; iter=1
   P1=max(Pmin,0.0_PR) ; P2=Pmax

   DO WHILE (.not.find.and.iter.lt.itermax)

    P=0.5_PR*(P1+P2)
    print*, '       iter=', iter, 'P=', P


    !!!---calcul de ir_eb1
    t%ir_eb1(it)=t%Nr ; find_eb1=.false. ; ir=1
    do while(.not.find_eb1.and.ir.lt.t%Nr)
      if( t%P(ir,it).le.P.and.t%P(ir+1,it).ge.P )then
          t%ir_eb1(it)=ir+1
          find_eb1=.true.
          print*, '       eb1=', t%ir_eb1(it), 0.32_PR/t%rho(ir+1)
      endif
      ir=ir+1
    enddo
    if(.not.find_eb1) t%ir_eb1(it)=0 

    !!!---calcul de ir_eb2
    t%ir_eb2(it)=1 ;  find_eb2=.false. ; ir=t%ir_eb1(it)
    do while(.not.find_eb2.and.ir.lt.t%Nr)
      if( t%P(ir,it).le.P.and.t%P(ir+1,it).ge.P )then
          t%ir_eb2(it)=ir+1
          find_eb2=.true.
          print*, '       eb2=', t%ir_eb2(it), 0.32_PR/t%rho(ir+1) 
      endif
      ir=ir+1
    enddo
    if(.not.find_eb2) t%ir_eb2(it)=0 

    if(find_eb1.and.find_eb2)then
      !!!---integration de P
      intPp=0.0_PR ; intPm=0.0_PR
      do ir=1,t%Nr
        if(ir.gt.t%ir_eb1(it).and.ir.le.t%ir_eb2(it))then
          if(t%P(ir,it).gt.P)then
              intPp=intPp+(t%P(ir,it)-P)*abs(1.0_PR/t%rho(ir-1)-1.0_PR/t%rho(ir))
          else
              intPm=intPm+(P-t%P(ir,it))*abs(1.0_PR/t%rho(ir-1)-1.0_PR/t%rho(ir))
          endif
        endif
      enddo
      
      !!!---reussite de la Dichotomie ? 
      if( abs(intPp-intPm).lt.tol*0.5_PR*(intPp+intPm) )then
        find=.true.
        t%Peq(it)=P
        print*, 'MAXWELL: Peq trouvé ! T=', t%T(it), 'iter=', iter, 'Peq=', P, intPp, intPm 
      else
        if(intPm.lt.intPp)then
          P1=P
        else 
          P2=P
        endif
      endif
      iter=iter+1
      if(iter.eq.itermax) print*, 'MAXWELL: itermax atteint !!!'

    else
      find=.true.
      print*, 'PROBLEM :(    eb1 eb2 : ', find_eb1, find_eb2  
    endif

  ENDDO

  
   print*, '         Dichotomie fished !!!'
   

else

   print*, '         pas de Peq :O !!! ',findpmin,findPmax,findPmax2

endif

!!!-----------Modification de la table

 do ir=1,t%Nr
  if(ir.ge.t%ir_eb1(it).and.ir.lt.t%ir_eb2(it))then
    t%P(ir,it)=t%Peq(it)
  endif
 enddo


ENDDO !!! it


end subroutine Maxwell_construct

!   SUBROUTINE  Pmin(Pt,Z,rmin,dr,T)
!   implicit none
!   procedure(pression_totale) :: Pt
!   integer, intent(in) :: Z
!   real(PR), intent(in) :: rmin,rmax,T
!   real(PR) :: P, tho, Pprev
!   integer :: i,j,k,Nr,ipic
!   logical :: find
!   
!
!   Pmin=0.0_PR ; Pmax=0.0_PR
!
!
!   !!!======Recherche de r1 r2 r3 r4 r5 r6
!
!   P=0.0_PR ; rho=rmin ; ipic=0
!
!   i=1 ; find=.false.
!
!   do while(.not.find)
!
!     rhop=rho
!
!     rho=rmin+real(i,PR)*(rmax-rmin)/real(Nr,PR) 
!     
!     Pp=P
!
!     P=Pt(Z,rho,T)
!
!     if(P.lt.Pp.and.ipic.eq.0)then
!       r1=rhop
!       r2=rho
!       p1=Pp
!       p2=P 
!       ipic=1
!     elseif(P.gt.P.and.ipic.eq.1)then
!       r3=rhop
!       r4=rho
!       p3=Pp
!       p4=P 
!       ipic=2      
!     elseif(P.gt.P1.and.ipic.eq.2)then
!       r5=rhop
!       r6=rho
!       p5=Pp
!       p6=P 
!       find=.true.     
!     elseif(rho.gt.rmax.and.ipic.eq.0)then
!       find=.true. 
!     endif
!
!     i=i+1
!
!   enddo
!
!   !!!!================ Calcul de Pmin P max rhoPmin, rhoPmax ===========================
!  if(find.and.ipic.eq.2)then   
!    
!   !!!---Dichotomie sur Pmin
!   tolP=1.0e-4_PR !!! 100Pa
!   itmax=1000 ; find=.false. ; it=0
!   do while (.not.find.and.it.lt.itmax)
!    rho=0.5_PR*(r1+r2)
!    P=P_tot(Z,rho,T_in)
!    if(abs(P-P_in).lt.tolP)then
!      find=.true.
!      rho_out=rho
!    else
!      if(P.lt.P_in)then
!        rhomin=rho
!      else 
!        rhomax=rho
!      endif
!    endif
!    it=it+1
!   enddo
!
!   
!
!
!
!
!  
!
!  endif
!
!
!
!
!
!  !!!!================ Calcul de Peq ===========================
!  if(find.and.ipic.eq.2)then   
!    
!    
!
!
!
!
!  
!
!  endif
!
!
!
!END FUNCTION Pmin
!

!!!========================================   INITIALIZATION DATAS  =============================================




subroutine INIT_EOS
implicit none

!!!!!============ noms

      element_name( 1:15) = (/ ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P' /)
      element_name(16:30) = (/ ' S', ' C', 'Ar', ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn' /)
      element_name(31:45) = (/ 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh' /)
      element_name(46:60) = (/ 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd' /)
      element_name(61:75) = (/ 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re' /)
      element_name(76:86) = (/ 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn' /)


!!!!!============ Gruneisen coefficients
      gruneisen( 1: 5) = (/ 1.79146_PR, 2.12563_PR, 1.67765_PR, 1.20888_PR, 1.16175_PR /)
      gruneisen( 6:10) = (/ 1.22269_PR, 1.62934_PR, 1.37132_PR, 1.56789_PR, 1.61263_PR /)
      gruneisen(11:15) = (/ 1.86850_PR, 1.59939_PR, 2.13600_PR, 1.51886_PR, 1.67780_PR /)
      gruneisen(16:20) = (/ 1.65461_PR, 1.68219_PR, 1.83847_PR, 2.14301_PR, 1.86535_PR /)
      gruneisen(21:25) = (/ 1.59661_PR, 1.43857_PR, 1.34148_PR, 1.29223_PR, 1.29483_PR /)
      gruneisen(26:30) = (/ 1.28196_PR, 1.25924_PR, 1.25737_PR, 1.28197_PR, 1.36904_PR /)
      gruneisen(31:35) = (/ 1.46800_PR, 1.52829_PR, 1.50856_PR, 1.60737_PR, 1.82110_PR /)
      gruneisen(36:40) = (/ 1.89290_PR, 2.19694_PR, 1.94694_PR, 1.69033_PR, 1.52959_PR /)
      gruneisen(41:45) = (/ 1.42424_PR, 1.36844_PR, 1.33680_PR, 1.32271_PR, 1.32631_PR /)
      gruneisen(46:50) = (/ 1.34886_PR, 1.40126_PR, 1.49116_PR, 1.57043_PR, 1.58419_PR /)
      gruneisen(51:55) = (/ 1.63400_PR, 1.68750_PR, 1.79625_PR, 1.97781_PR, 2.28434_PR /)
      gruneisen(56:60) = (/ 2.00220_PR, 1.72999_PR, 1.69308_PR, 1.69829_PR, 1.68300_PR /)
      gruneisen(61:65) = (/ 1.66813_PR, 1.67849_PR, 1.84531_PR, 1.66448_PR, 1.65060_PR /)
      gruneisen(66:70) = (/ 1.64207_PR, 1.63544_PR, 1.62752_PR, 1.61903_PR, 1.79243_PR /)
      gruneisen(71:75) = (/ 1.60946_PR, 1.49241_PR, 1.41479_PR, 1.37144_PR, 1.34784_PR /)
      gruneisen(76:80) = (/ 1.33266_PR, 1.33626_PR, 1.36141_PR, 1.39305_PR, 1.52937_PR /)
      gruneisen(81:86) = (/ 1.59162_PR, 1.61599_PR, 1.68526_PR, 1.71311_PR, 1.84537_PR, 2.10319_PR /)

!!!!!============ Debye Temperatures eV
      debye_temp( 1: 5) = (/ 0.00395978_PR, 0.000868131_PR, 0.0233086_PR, 0.0592689_PR, 0.0741264_PR /)
      debye_temp( 6:10) = (/ 0.0871836_PR, 0.00801577_PR, 0.00875263_PR, 0.00729161_PR, 0.00466196_PR /)
      debye_temp(11:15) = (/ 0.0148618_PR, 0.0273893_PR, 0.0301668_PR, 0.0373163_PR, 0.0141415_PR /)
      debye_temp(16:20) = (/ 0.0155253_PR, 0.00990004_PR, 0.00609248_PR, 0.00967749_PR, 0.0208953_PR /)
      debye_temp(21:25) = (/ 0.031323_PR, 0.0357402_PR, 0.0403533_PR, 0.0414447_PR, 0.033887_PR /)
      debye_temp(26:30) = (/ 0.0368149_PR, 0.0365957_PR, 0.0356441_PR, 0.0302851_PR, 0.0195756_PR /)
      debye_temp(31:35) = (/ 0.0117083_PR, 0.0219696_PR, 0.0208158_PR, 0.0128074_PR, 0.00797565_PR /)
      debye_temp(36:40) = (/ 0.00493837_PR, 0.00647994_PR, 0.0139088_PR, 0.0213249_PR, 0.0257068_PR /)
      debye_temp(41:45) = (/ 0.0314038_PR, 0.0334277_PR, 0.0311502_PR, 0.0322989_PR, 0.0294256_PR /)
      debye_temp(46:50) = (/ 0.0256731_PR, 0.0198242_PR, 0.0125638_PR, 0.00990966_PR, 0.010498_PR /)
      debye_temp(51:55) = (/ 0.0133653_PR, 0.0113588_PR, 0.00761392_PR, 0.00430148_PR, 0.00468924_PR /)
      debye_temp(56:60) = (/ 0.0102879_PR, 0.013345_PR, 0.0128289_PR, 0.013399_PR, 0.0138865_PR /)
      debye_temp(61:65) = (/ 0.0142989_PR, 0.0138945_PR, 0.0110479_PR, 0.014918_PR, 0.0151129_PR /)
      debye_temp(66:70) = (/ 0.0153073_PR, 0.015503_PR, 0.0155385_PR, 0.0156902_PR, 0.0106347_PR /)
      debye_temp(71:75) = (/ 0.0159904_PR, 0.0198009_PR, 0.0241115_PR, 0.0264354_PR, 0.0259907_PR /)
      debye_temp(76:80) = (/ 0.0256124_PR, 0.0230175_PR, 0.0191656_PR, 0.0148823_PR, 0.00545389_PR /)
      debye_temp(81:86) = (/ 0.00806864_PR, 0.00800644_PR, 0.00717512_PR, 0.00685952_PR, 0.00647832_PR, 0.00320344_PR /)

!!!!!============ Solid density g / cc !!! 25°C 298.15 K
      rhos( 1: 5) = (/ 0.0760_PR, 0.1248_PR, 0.535_PR, 1.848_PR, 2.460_PR /)
      rhos( 6:10) = (/ 2.260_PR, 1.026_PR, 2.000_PR, 1.516_PR, 1.444_PR /)
      rhos(11:15) = (/ 0.968_PR, 1.738_PR, 2.698_PR, 2.330_PR, 1.823_PR /)
      rhos(16:20) = (/ 1.960_PR, 2.030_PR, 1.656_PR, 0.856_PR, 1.550_PR /)
      rhos(21:25) = (/ 2.985_PR, 4.507_PR, 6.110_PR, 7.140_PR, 7.470_PR /)
      rhos(26:30) = (/ 7.874_PR, 8.9_PR, 8.908_PR, 8.920_PR, 7.140_PR /)
      rhos(31:35) = (/ 5.904_PR, 5.323_PR, 5.727_PR, 4.819_PR, 3.120_PR /)
      rhos(36:40) = (/ 2.82_PR, 1.532_PR, 2.630_PR, 4.472_PR, 6.511_PR /)
      rhos(41:45) = (/ 8.570_PR, 10.280_PR, 11.5_PR, 12.370_PR, 12.450_PR /)
      rhos(46:50) = (/ 12.023_PR, 10.490_PR, 8.650_PR, 7.310_PR, 7.310_PR /)
      rhos(51:55) = (/ 6.697_PR, 6.240_PR, 4.940_PR, 3.540_PR, 1.879_PR /)
      rhos(56:60) = (/ 3.510_PR, 6.146_PR, 6.689_PR, 6.640_PR, 7.010_PR /)
      rhos(61:65) = (/ 7.264_PR, 7.353_PR, 5.244_PR, 7.901_PR, 8.219_PR /)
      rhos(66:70) = (/ 8.551_PR, 8.795_PR, 9.066_PR, 9.321_PR, 6.570_PR /)
      rhos(71:75) = (/ 9.841_PR, 13.310_PR, 16.650_PR, 19.250_PR, 21.020_PR /)
      rhos(76:80) = (/ 22.59_PR, 22.56_PR, 21.090_PR, 19.3_PR, 13.534_PR /)
      rhos(81:86) = (/ 11.850_PR, 11.340_PR, 9.780_PR, 9.196_PR, 7.000_PR, 4.40_PR /)

      
!!!!!============ Cohesive energy eV / atom
      cohesive_energy( 1: 5) = (/ 9.3146e-3_PR, 8.60243e-4_PR, 1.52356_PR, 3.07822_PR, 5.25474_PR /)
      cohesive_energy( 6:10) = (/ 7.41053_PR, 0.0289166_PR, 0.0353425_PR, 0.0338915_PR, 0.0181377_PR /)
  !!! cohesive_energy(11:15) = (/ 1.0126_PR, 1.32664_PR, 3.03676_PR, 3.72081_PR, 0.128518_PR /) !!!
      cohesive_energy(11:15) = (/ 1.0126_PR, 1.32664_PR, 3.39000_PR, 3.72081_PR, 0.128518_PR /) !!!
      cohesive_energy(16:20) = (/ 0.101571_PR, 0.105717_PR, 0.0673684_PR, 0.79702_PR, 1.60648_PR /)
      cohesive_energy(21:25) = (/ 3.29587_PR, 4.40486_PR, 4.69506_PR, 3.51352_PR, 2.28016_PR /)
      cohesive_energy(26:30) = (/ 3.59644_PR, 3.88664_PR, 3.91773_PR, 3.10931_PR, 1.23336_PR /)
      cohesive_energy(31:35) = (/ 2.65328_PR, 3.4617_PR, 0.335806_PR, 0.269474_PR, 0.153393_PR /)
      cohesive_energy(36:40) = (/ 0.0934867_PR, 0.746235_PR, 1.41992_PR, 3.93846_PR, 6.01134_PR /)
      cohesive_energy(41:45) = (/ 7.15142_PR, 6.21862_PR, 5.70041_PR, 6.01134_PR, 5.13036_PR /)
      cohesive_energy(46:50) = (/ 3.93846_PR, 2.64292_PR, 1.03644_PR, 2.38381_PR, 3.00567_PR /)
      cohesive_energy(51:55) = (/ 0.704777_PR, 0.49749_PR, 0.216615_PR, 0.131006_PR, 0.673684_PR /)
      cohesive_energy(56:60) = (/ 1.45101_PR, 4.14575_PR, 3.62753_PR, 3.42024_PR, 2.95385_PR /)
      cohesive_energy(61:65) = (/ 3.00567_PR, 1.81377_PR, 1.81377_PR, 3.16113_PR, 3.05749_PR /)
      cohesive_energy(66:70) = (/ 2.90202_PR, 2.74656_PR, 2.95385_PR, 2.59109_PR, 1.6583_PR /)
  !!! cohesive_energy(71:75) = (/ 4.30122_PR, 6.52956_PR, 7.61781_PR, 8.2915_PR, 7.30688_PR /)
      cohesive_energy(71:75) = (/ 4.30122_PR, 6.52956_PR, 8.10000_PR, 8.2915_PR, 7.30688_PR /)
      cohesive_energy(76:80) = (/ 6.52956_PR, 5.80405_PR, 5.07854_PR, 3.42024_PR, 0.613571_PR /)
      cohesive_energy(81:86) = (/ 1.71012_PR, 1.84486_PR, 1.6583_PR, 1.03644_PR, 0.414575_PR, 0.176194_PR /)

!!!!!============ melting point eV
      melting_point( 1: 5) = (/ 0.00121763_PR, 0.0000825656_PR, 0.0394307_PR, 0.135594_PR, 0.20408_PR /)
      melting_point( 6:10) = (/ 0.332274_PR, 0.00547975_PR, 0.00476708_PR, 0.00465409_PR, 0.00213454_PR /)
      melting_point(11:15) = (/ 0.0322327_PR, 0.0802321_PR, 0.081129_PR, 0.146632_PR, 0.0275813_PR /)
      melting_point(16:20) = (/ 0.0337528_PR, 0.0149183_PR, 0.0072875_PR, 0.0292482_PR, 0.096919_PR /)
      melting_point(21:25) = (/ 0.15767_PR, 0.168708_PR, 0.18974_PR, 0.189479_PR, 0.132031_PR /)
      melting_point(26:30) = (/ 0.157409_PR, 0.153672_PR, 0.150196_PR, 0.118005_PR, 0.0602016_PR /)
      melting_point(31:35) = (/ 0.0263263_PR, 0.105289_PR, 0.0947462_PR, 0.0429472_PR, 0.0231053_PR /)
      melting_point(36:40) = (/ 0.0100634_PR, 0.0271563_PR, 0.0912698_PR, 0.156366_PR, 0.18496_PR /)
      melting_point(41:45) = (/ 0.239019_PR, 0.251708_PR, 0.211207_PR, 0.22659_PR, 0.194433_PR /)
      melting_point(46:50) = (/ 0.158878_PR, 0.107329_PR, 0.0516444_PR, 0.0373501_PR, 0.0438971_PR /)
      melting_point(51:55) = (/ 0.0785486_PR, 0.0628072_PR, 0.0336216_PR, 0.0140231_PR, 0.0262115_PR /)
      melting_point(56:60) = (/ 0.0869242_PR, 0.103698_PR, 0.0930949_PR, 0.104654_PR, 0.112476_PR /)
      melting_point(61:65) = (/ 0.119342_PR, 0.116909_PR, 0.0951808_PR, 0.137854_PR, 0.141591_PR /)
      melting_point(66:70) = (/ 0.146458_PR, 0.151847_PR, 0.153846_PR, 0.158018_PR, 0.09492_PR /)
      melting_point(71:75) = (/ 0.168273_PR, 0.217812_PR, 0.285951_PR, 0.32115_PR, 0.300639_PR /)
      melting_point(76:80) = (/ 0.287341_PR, 0.238063_PR, 0.177425_PR, 0.116229_PR, 0.020365_PR /)
      melting_point(81:86) = (/ 0.0501608_PR, 0.0521997_PR, 0.0473188_PR, 0.0458152_PR, 0.049987_PR, 0.0175691_PR /)
 


!!!!!============ boiling point eV
      boiling_point( 1: 5) = (/ 0.00176256_PR, 0.000366765_PR, 0.140375_PR, 0.23841_PR, 0.371384_PR /)
      boiling_point( 6:10) = (/ 0.373731_PR, 0.00672345_PR, 0.00784373_PR, 0.00739006_PR, 0.00235269_PR /)
      boiling_point(11:15) = (/ 0.100482_PR, 0.118473_PR, 0.242669_PR, 0.275782_PR, 0.0481184_PR /)
      boiling_point(16:20) = (/ 0.0623909_PR, 0.0207813_PR, 0.00759169_PR, 0.0897054_PR, 0.152716_PR /)
      boiling_point(21:25) = (/ 0.269698_PR, 0.309417_PR, 0.319846_PR, 0.25588_PR, 0.202864_PR /)
      boiling_point(26:30) = (/ 0.272393_PR, 0.278129_PR, 0.276912_PR, 0.278129_PR, 0.102568_PR /)
      boiling_point(31:35) = (/ 0.215292_PR, 0.268829_PR, 0.0771033_PR, 0.0832739_PR, 0.0288675_PR /)
      boiling_point(36:40) = (/ 0.0104233_PR, 0.0835347_PR, 0.143851_PR, 0.314458_PR, 0.406931_PR /)
      boiling_point(41:45) = (/ 0.436046_PR, 0.426921_PR, 0.394416_PR, 0.384421_PR, 0.344877_PR /)
      boiling_point(46:50) = (/ 0.281258_PR, 0.211642_PR, 0.0904007_PR, 0.20382_PR, 0.249883_PR /)
      boiling_point(51:55) = (/ 0.161668_PR, 0.109608_PR, 0.0397575_PR, 0.0143534_PR, 0.0820572_PR /)
      boiling_point(56:60) = (/ 0.186264_PR, 0.3248_PR, 0.315761_PR, 0.309678_PR, 0.293164_PR /)
      boiling_point(61:65) = (/ 0.284473_PR, 0.180441_PR, 0.156453_PR, 0.306201_PR, 0.304463_PR /)
      boiling_point(66:70) = (/ 0.246841_PR, 0.2584_PR, 0.273001_PR, 0.193217_PR, 0.127686_PR /)
      boiling_point(71:75) = (/ 0.319412_PR, 0.423792_PR, 0.498101_PR, 0.506531_PR, 0.510095_PR /)
      boiling_point(76:80) = (/ 0.459339_PR, 0.408582_PR, 0.356175_PR, 0.271958_PR, 0.0547436_PR /)
      boiling_point(81:86) = (/ 0.15176_PR, 0.175747_PR, 0.159669_PR, 0.107348_PR, 0.0530289_PR, 0.0183774_PR /)

!!!!!============ Bulk Modulus GPa
      !bulk_modulus( 1: 5) = (/ 0.166_PR, 0.47_PR, 13.10_PR, 116.8_PR, 210._PR /)
      bulk_modulus( 1: 5) = (/ 0.166_PR, 0.0047_PR, 13.10_PR, 116.8_PR, 210._PR /)
      !bulk_modulus( 6:10) = (/ 444.3_PR, 301._PR, 2.97_PR, 15.0_PR, 1.07_PR /)
      bulk_modulus( 6:10) = (/ 444.3_PR, 301._PR, 2.97_PR, 1.50_PR, 1.07_PR /) ! Fix for flourine
      bulk_modulus(11:15) = (/ 6.43_PR, 33._PR, 75.8_PR, 100.1_PR, 36.0_PR /)
      bulk_modulus(16:20) = (/ 14.5_PR, 15.18_PR, 2.86_PR, 4.25_PR, 17.4_PR /)
      bulk_modulus(21:25) = (/ 60.0_PR, 110.02_PR, 154._PR, 193._PR, 158._PR /)
      bulk_modulus(26:30) = (/ 163.4_PR, 203._PR, 161._PR, 137.4_PR, 65.0_PR /)
      bulk_modulus(31:35) = (/ 55.0_PR, 87.4_PR, 78.107_PR, 48.1_PR, 14.1_PR /)
      bulk_modulus(36:40) = (/ 3.34_PR, 5.83_PR, 11.7_PR, 34.0_PR, 92.0_PR /)
      bulk_modulus(41:45) = (/ 173._PR, 265._PR, 344._PR, 267._PR, 318._PR /)
      bulk_modulus(46:50) = (/ 195._PR, 109._PR, 42.0_PR, 36.8_PR, 55.8_PR /)
      bulk_modulus(51:55) = (/ 46.7_PR, 24.0_PR, 30.38_PR, 3.63_PR, 1.723_PR /)
      bulk_modulus(56:60) = (/ 7.65_PR, 24.8_PR, 23.08_PR, 27.16_PR, 25.38_PR /)
      bulk_modulus(61:65) = (/ 38.0_PR, 30.7_PR, 11.7_PR, 34.0_PR, 37.0_PR /)
      bulk_modulus(66:70) = (/ 36.3_PR, 38.9_PR, 44.0_PR, 44.5_PR, 14.6_PR /)
      bulk_modulus(71:75) = (/ 47.0_PR, 122._PR, 194._PR, 307._PR, 372._PR /)
      bulk_modulus(76:80) = (/ 412._PR, 383._PR, 288.4_PR, 180._PR, 292._PR /)
      bulk_modulus(81:86) = (/ 35.3_PR, 43.2_PR, 54.7_PR, 39.45_PR, 30.00_PR , 3.50_PR /)

!!!!!============ Bulk Modulus derivativ with p
      bm_pderiv( 1: 5) = (/ 7.102_PR, 4.01_PR, 2.80_PR, 4.6_PR, 2.23_PR /)
      bm_pderiv( 6:10) = (/ 4.04_PR, 4.02_PR, 7.78_PR, 3.00_PR, 8.40_PR /)
      bm_pderiv(11:15) = (/ 3.84_PR, 4.32_PR, 4.109_PR, 4.16_PR, 4.50_PR /)
      bm_pderiv(16:20) = (/ 7.00_PR, 1.59_PR, 7.20_PR, 3.63_PR, 3.22_PR /)
      bm_pderiv(21:25) = (/ 2.80_PR, 3.59_PR, 4.20_PR, 4.80_PR, 4.60_PR /)
      bm_pderiv(26:30) = (/ 5.38_PR, 3.60_PR, 7.55_PR, 5.52_PR, 4.60_PR /)
      bm_pderiv(31:35) = (/ 4.70_PR, 4.31_PR, 4.317_PR, 4.33_PR, 4.60_PR /)
      bm_pderiv(36:40) = (/ 7.20_PR, 2.08_PR, 2.41_PR, 5.00_PR, 4.00_PR /)
      bm_pderiv(41:45) = (/ 4.10_PR, 4.70_PR, 4.54_PR, 4.50_PR, 4.86_PR /)
      bm_pderiv(46:50) = (/ 5.42_PR, 5.87_PR, 6.50_PR, 5.25_PR, 4.08_PR /)
      bm_pderiv(51:55) = (/ 4.90_PR, 2.30_PR, 6.13_PR, 8.08_PR, 3.87_PR /)
      bm_pderiv(56:60) = (/ 3.45_PR, 2.80_PR, 4.945_PR, 0.096_PR, 3.119_PR /)
      bm_pderiv(61:65) = (/ 1.50_PR, 2.50_PR, 3.00_PR, 4.20_PR, 2.70_PR /)
      bm_pderiv(66:70) = (/ 3.71_PR, 2.83_PR, 0.8332_PR, 2.71_PR, 1.205_PR /)
      bm_pderiv(71:75) = (/ 3.16_PR, 3.38_PR, 3.55_PR, 4.30_PR, 4.05_PR /)
      bm_pderiv(76:80) = (/ 4.40_PR, 3.10_PR, 5.05_PR, 5.61_PR, 5.50_PR /)
      bm_pderiv(81:86) = (/ 4.11_PR, 4.90_PR, 4.90_PR, 4.89_PR, 6.00_PR, 7.00_PR /)      
!!!!!============ First ionization energy (eV) NIST
   Eio1(1)=13.59843449_PR ; Eio1(2)=24.58738880_PR ;  Eio1(3)=5.39171495_PR
   Eio1(4)=9.322699_PR    ; Eio1(5)=8.298019_PR    ;  Eio1(6)=11.2602880_PR
   Eio1(7)=14.53413_PR    ; Eio1(8)=13.618055_PR   ;  Eio1(9)=17.42282_PR
   Eio1(10)=21.564540_PR  ; Eio1(11)=5.1390769_PR  ;  Eio1(12)=7.646236_PR
   Eio1(13)=5.985769_PR   ; Eio1(14)=8.15168_PR    ;  Eio1(15)=10.486686_PR
   Eio1(16)=10.36001_PR   ; Eio1(17)=12.967632_PR  ;  Eio1(18)=15.7596117_PR
   Eio1(19)=4.34066369_PR ; Eio1(20)=6.1131554_PR  ;  Eio1(21)=6.56149_PR
   Eio1(22)=6.828120_PR   ; Eio1(23)=6.746187_PR   ;  Eio1(24)=6.76651_PR
   Eio1(25)=7.4340379_PR  ; Eio1(26)=7.9024681_PR  ;  Eio1(27)=7.88101_PR
   Eio1(28)=7.639878_PR   ; Eio1(29)=7.726380_PR   ;  Eio1(30)=9.394197_PR
   Eio1(31)=5.9993020_PR  ; Eio1(32)=7.899435_PR   ;  Eio1(33)=9.78855_PR
   Eio1(34)=9.752392_PR   ; Eio1(35)=11.81381_PR   ;  Eio1(36)=13.9996053_PR
   Eio1(37)=4.1771280_PR  ; Eio1(38)=5.69486740_PR ;  Eio1(39)=6.21726_PR
   Eio1(40)=6.63412_PR    ; Eio1(41)=6.75885_PR    ;  Eio1(42)=7.09243_PR
   Eio1(43)=7.11938_PR    ; Eio1(44)=7.36050_PR    ;  Eio1(45)=7.45890_PR
   Eio1(46)=8.336839_PR   ; Eio1(47)=7.576234_PR   ;  Eio1(48)=8.993820_PR
   Eio1(49)=5.7863556_PR  ; Eio1(50)=7.343918_PR   ;  Eio1(51)=8.608389_PR
   Eio1(52)=9.00966_PR    ; Eio1(53)=10.451260_PR  ;  Eio1(54)=12.1298436_PR
   Eio1(55)=3.893905695_PR; Eio1(56)=5.2116646_PR  ;  Eio1(57)=5.5769_PR
   Eio1(58)=5.5386_PR     ; Eio1(59)=5.4702_PR     ;  Eio1(60)=5.5250_PR
   Eio1(61)=5.577_PR      ; Eio1(62)=5.64371_PR    ;  Eio1(63)=5.670385_PR
   Eio1(64)=6.14980_PR    ; Eio1(65)=5.8638_PR     ;  Eio1(66)=5.93905_PR
   Eio1(67)=6.0215_PR     ; Eio1(68)=6.1077_PR     ;  Eio1(69)=6.18431_PR
   Eio1(70)=6.254160_PR   ; Eio1(71)=5.425871_PR   ;  Eio1(72)=6.825070_PR
   Eio1(73)=7.549571_PR   ; Eio1(74)=7.86403_PR    ;  Eio1(75)=7.83352_PR
   Eio1(76)=8.43823_PR    ; Eio1(77)=8.96702_PR    ;  Eio1(78)=8.95883_PR
   Eio1(79)=9.225554_PR   ; Eio1(80)=10.437504_PR  ;  Eio1(81)=6.1082873_PR
   Eio1(82)=7.4166799_PR  ; Eio1(83)=7.285516_PR   ;  Eio1(84)=8.414_PR
   Eio1(85)=9.31751_PR    ; Eio1(86)=10.74850_PR   ;  Eio1(87)=4.0727410_PR
   Eio1(88)=5.2784239_PR  ; Eio1(89)=5.380226_PR   ;  Eio1(90)=6.30670_PR
   Eio1(91)=5.89_PR       ; Eio1(92)=6.19405_PR

!!!!!============ Polarizability Schwerdtfeger 2014 (p=alpa*E) 
!!!    "Table of experimental and calculated static dipole polarizabilities
!!!     for the electronic ground states of the neutral elements"
!!!     in atomic unit: 1 a.u. = 0.14818474 A^3 = 1.6487773e-41 C m^2 / V 

   pola(1 )=4.5_PR    ; pola(2 )=1.383191_PR ; pola(3 )=164.0_PR
   pola(4 )=37.755_PR ; pola(5 )=20.5_PR     ; pola(6 )=11.0_PR
   pola(7 )=7.6_PR    ; pola(8 )=6.04_PR     ; pola(9 )=3.76_PR
   pola(10)=2.670_PR  ; pola(11)=162.6_PR    ; pola(12)=71.7_PR
   pola(13)=46_PR     ; pola(14)=36.7_PR     ; pola(15)=24.7_PR
   pola(16)=19.6_PR   ; pola(17)=14.7_PR     ; pola(18)=11.10_PR
   pola(19)=290.6_PR  ; pola(20)=169_PR      ; pola(21)=120_PR
   pola(22)=99_PR     ; pola(23)=84_PR       ; pola(24)=78_PR
   pola(25)=63_PR     ; pola(26)=57_PR       ; pola(27)=51_PR
   pola(28)=46_PR     ; pola(29)=53.44_PR    ; pola(30)=38.8_PR
   pola(31)=54.9_PR   ; pola(32)=39.43_PR    ; pola(33)=29.1_PR
   pola(34)=26.24_PR  ; pola(35)=21.9_PR     ; pola(36)=17.075_PR
   pola(37)=316_PR    ; pola(38)=186_PR      ; pola(39)=153_PR
   pola(40)=121_PR    ; pola(41)=106_PR      ; pola(42)=86_PR
   pola(43)=77_PR     ; pola(44)=65_PR       ; pola(45)=58_PR
   pola(46)=32_PR     ; pola(47)=52.2_PR     ; pola(48)=49.65_PR
   pola(49)=68.7_PR   ; pola(50)=42.4_PR     ; pola(51)=45_PR
   pola(52)=37_PR     ; pola(53)=35.1_PR     ; pola(54)=27.815_PR
   pola(55)=401.0_PR  ; pola(56)=268_PR      ; pola(57)=210_PR
   pola(58)=200_PR    ; pola(59)=190_PR      ; pola(60)=212_PR
   pola(61)=203_PR    ; pola(62)=194_PR      ; pola(63)=187_PR
   pola(64)=159_PR    ; pola(65)=172_PR      ; pola(66)=165_PR
   pola(67)=159_PR    ; pola(68)=153_PR      ; pola(69)=147_PR
   pola(70)=142_PR    ; pola(71)=148_PR      ; pola(72)=109_PR
   pola(73)=88_PR     ; pola(74)=75_PR       ; pola(75)=65_PR
   pola(76)=57_PR     ; pola(77)=51_PR       ; pola(78)=44_PR
   pola(79)=35.1_PR   ; pola(80)=33.91_PR    ; pola(81)=51_PR
   pola(82)=47.1_PR   ; pola(83)=50_PR       ; pola(84)=46_PR
   pola(85)=45.6_PR   ; pola(86)=33.18_PR    ; pola(87)=317.8_PR
   pola(88)=246.2_PR  ; pola(89)=217_PR      ; pola(90)=217_PR
   pola(91)=171_PR    ; pola(92)=137_PR

   !!!! pola=pola*1.6487773e-41_PR !!! Cm/(V/m)
   pola=pola*0.1481847e-30_PR !!! m^3
!!!---numéro atomique A
!A(13)=26.9815386_PR
!A(26)=55.845_PR
!A(73)=180.94788_PR
!A(80)=200.59_PR
!A(82)=207.2_PR

!!!======== Mass number

A(1:92)= (/ 1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,30,40,&
            45,48,51,52,55,56,58,58,64,65,70,73,75,79,80,84,85,88,89,91,&
            93,96,98,101,103,106,108,112,115,119,122,128,127,131,133,137,139,140,141,144,&
            145,150,152,157,159,163,165,167,169,173,175,178,181,184,186,190,192,195,197,201,&
            204,207,209,209,210,222,223,226,227,232,231,238/)
            !237,244,243,247,247,251,252,257,&
            !258,259,262,261,268,263,264,269,268,272,273,277,286,289,288,292,292,293 /)


symbol(1:92)=(/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
               'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
               'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
               'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
               'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U '/)
                 !'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','&
                 ! Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ','Mt ','Ds ','Rg','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo'/)

!!!---densité,de,référence,du,solide,g,/,cc
!rhos(14)=2.6989_PR !!!! 2.773_PR 
!rhos(80)=14.38_PR !!!! 13.546,14.38,16.39 

!!!---Températures de fusion du solide K
Tms(:)=melting_point(:)*eV

!Tms(14)=934.0_PR !!!!K
!Tms(80)=234.3_PR !!!!K

!!!---énergie de cohesion
!Ecoh(14)=3.39_PR !!! eV/atom
!Ecoh(80)=0.67_PR !!! eV/atom

!!!---Electronic heat capacity constant : C = ge*T + A*T^3 Kittel p157(djvu:173) ge en mJ/mol/K**2
ge(13)=1.35_PR
ge(26)=4.98_PR
ge(29)=0.695_PR
ge(73)=5.9_PR
ge(80)=1.79_PR
ge(82)=2.98_PR

!!!---critical parameters : Ray 2017 EOS for metals in the liquid-vapor region
!!!                         MPa, K, g cm^-3
!!!---Na: 
Tcrit(11)=2448.0_PR ; rcrit(11)=0.20_PR ; Pcrit(13)=0.48e2_PR 
!!!---Al: 
Tcrit(13)=5700.0_PR ; rcrit(13)=0.32_PR ; Pcrit(13)=1.87e2_PR 
!!!---Fe:
Tcrit(26)=6900.0_PR ; rcrit(26)=1.77_PR ; Pcrit(26)=8.82e2_PR
!!!---Cu:
Tcrit(29)=7800.0_PR ; rcrit(29)=2.31_PR ; Pcrit(29)=8.94e2_PR  
!!!---Mo:
Tcrit(42)=9500.0_PR ; rcrit(42)=2.82_PR ; Pcrit(42)=9.64e2_PR  
!!!---Hg:
Tcrit(80)=1755.0_PR ; rcrit(80)=5.99_PR ; Pcrit(80)=1.65e2_PR  
!!!---Pb:
Tcrit(82)=5400.0_PR ; rcrit(82)=2.80_PR ; Pcrit(82)=4.01e2_PR  
!!!---Ta:
Tcrit(73)=9900.0_PR ; rcrit(73)=4.19_PR ; Pcrit(73)=10.6e2_PR  
!!!---W:
Tcrit(74)=12000.0_PR; rcrit(74)=4.04_PR ; Pcrit(74)=10.21e2_PR  
!!!---U:
Tcrit(92)=7000.0_PR ; rcrit(92)=8.11_PR ; Pcrit(92)=7.11e2_PR  

!!!!---Likalter 2002 (bar, K, g/cm^3) !!!!---Likalter 2002: autre ref
!Pcrit(13)=3120.0_PR                     !Pcrit(13)=4470.0_PR 
!Tcrit(13)=8860.0_PR                     !Tcrit(13)=8000.0_PR 
!rcrit(13)=0.28_PR                     !rcrit(13)=0.64_PR

!!!---Cu:
!!!---Likalter 2002 (bar, K, g/cm^3) 
!Pcrit(29)  =6510.0_PR ; Tcrit(29)  =8440.0_PR ; rcrit(29)=1.940_PR                   

!!!!---Likalter 2002: autre ref
   !Pcrit(29)=7460.0_PR ; Tcrit(29)=8390.0_PR ; !rcrit(29)=2.39.0_PR

!!!---Lennard Jones parameters for metals RAY 2017
m(11)=5.4823_PR ; n(11)=0.5898_PR !!! Na
m(13)=1.9508_PR ; n(13)=0.4531_PR !!! Al 
m(26)=1.8988_PR ; n(26)=0.6950_PR !!! Fe
m(29)=2.3857_PR ; n(29)=0.6727_PR !!! Cu
m(42)=3.2131_PR ; n(42)=0.4887_PR !!! Mo
m(73)=2.6315_PR ; n(73)=0.5568_PR !!! Ta
m(74)=2.6278_PR ; n(74)=0.4635_PR !!! W
m(80)=4.6765_PR ; n(80)=0.747_PR  !!! Hg
m(82)=1.4414_PR ; n(82)=1.3008_PR !!! Pb
m(92)=6.5917_PR ; n(92)=0.3989_PR !!! U

!m(13)=2.0_PR    ; n(13)=0.5_PR !!! Al ---RAY
!m(42)=2.0_PR    ; n(42)=0.65_PR !!! Mo
!m(73)=3.0_PR ; n(73)=0.5_PR !!! Ta
!m(80)=3.4_PR    ; n(80)=1.05_PR   !!! Hg
!m(82)=2.3_PR   ; n(82)=1.0_PR !!! Pb

!!!---Default Lenard Jones
!n(:)=m(:)=6.0_PR

!!!---Fermi energy of metals (eV) : 
E_Fermi=0.0_PR
!!! Ashcroft, N. W. and Mermin, N. D., Solid State Physics, Saunders, 1976.
E_Fermi(3 ) = 4.74_PR !!! Li 
E_Fermi(4 ) = 14.3_PR !!! Be
E_Fermi(11) = 3.24_PR !!! Na
E_Fermi(12) = 7.08_PR !!! Mg
E_Fermi(13) = 11.7_PR !!! Al
E_Fermi(19) = 2.12_PR !!! K 
E_Fermi(20) = 4.69_PR !!! Ca
E_Fermi(25) = 10.9_PR !!! Mn
E_Fermi(26) = 11.1_PR !!! Fe
E_Fermi(29) = 7.00_PR !!! Cu
E_Fermi(30) = 9.47_PR !!! Zn
E_Fermi(31) = 10.4_PR !!! Ga
E_Fermi(37) = 1.85_PR !!! Rb
E_Fermi(38) = 3.93_PR !!! Sr
E_Fermi(41) = 5.32_PR !!! Nb
E_Fermi(47) = 5.49_PR !!! Ag
E_Fermi(48) = 7.47_PR !!! Cd
E_Fermi(49) = 8.63_PR !!! In
E_Fermi(50) = 10.2_PR !!! Sn
E_Fermi(51) = 10.9_PR !!! Sb
E_Fermi(55) = 1.59_PR !!! Cs
E_Fermi(56) = 3.64_PR !!! Ba
E_Fermi(79) = 5.53_PR !!! Au
E_Fermi(80) = 7.13_PR !!! Hg
E_Fermi(81) = 8.15_PR !!! Tl
E_Fermi(82) = 9.47_PR !!! Pb
E_Fermi(83) = 9.90_PR !!! Bi
!!! Panda 2012 : Electronic structure and equilibrium properties
!!! of hcp titanium and zirconium
E_Fermi(40) = 6.9389_PR !!! Zr
E_Fermi(22) = 8.8437_PR !!! Ti


v_Fermi(:)=sqrt(2.0_PR*E_Fermi(:)*qe/(me*1.0e-3_PR))


!!!======== Ratio of conductivities sig_solid/sig_liquid from
!!! "the electronic properties of liguid metals" N. E. CUSACK
sigSL(:)=1.0_PR
sigSL(3 )=1.64_PR
sigSL(11)=1.451_PR
sigSL(12)=1.78_PR
sigSL(13)=2.20_PR
sigSL(14)=0.034_PR
sigSL(19)=1.56_PR
sigSL(25)=0.61_PR
sigSL(26)=1.09_PR
sigSL(27)=1.05_PR
sigSL(28)=1.3_PR
sigSL(29)=2.04_PR
sigSL(30)=2.24_PR
sigSL(31)=3.12_PR
sigSL(32)=0.053_PR
sigSL(34)=1.0e3_PR
sigSL(37)=1.60
sigSL(42)=1.23_PR
sigSL(47)=2.09_PR
sigSL(48)=1.97_PR
sigSL(49)=2.18_PR
sigSL(50)=2.10_PR
sigSL(51)=0.61_PR
sigSL(52)=0.091_PR
sigSL(55)=1.66_PR
sigSL(56)=1.62_PR
sigSL(74)=1.08_PR
sigSL(78)=1.40_PR
sigSL(79)=2.28_PR
sigSL(80)=4.94_PR
sigSL(81)= 2.06_PR
sigSL(82)=1.94_PR
sigSL(83)=0.35_PR

!!!======== Ref conductivity (S/m) for solids at 300 K
sig0(:)=1.0e7_PR
sig0(13)=37.7e6_PR 
sig0(26)=9.93e6_PR
sig0(29)=59.6e6_PR

!!!======== Max conductivity for solids at 300 K obtained before scaling from
!!!         Lee more algorithm
sigmax(:)=sig0(:)
sigmax(13)=86191303.0_PR
sigmax(26)=268630338.0_PR
sigmax(29)=208077178.0_PR

!!!======== Ref thermal conductivity (S/m) for solids at 300 K
lamb0(:)=1.0e2_PR
lamb0(13)=237.0_PR 
lamb0(26)=80.2_PR
lamb0(29)=401.0_PR

!!!======== Max thermal conductivity for solids at 300 K obtained before scaling from
!!!         Lee more algorithm

lambmax(:)=lamb0(:)
lambmax(13)=632.289_PR
lambmax(26)=1970.64_PR
lambmax(29)=1526.43_PR



end subroutine INIT_EOS

subroutine rhominmax_TP(Z,T,P0,rho0,dPdr,rhomin,rhomax)
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: T, P0, rho0
   real(PR), intent(out) :: rhomin, rhomax, dPdr
   real(PR) :: r1(1:3), r2(1:3), P1(1:3), P2(1:3), drho(1:3)
   real(PR) :: P, rho, rhom1, Pm1
   integer :: ir, i, N, Nr


   rhomin=1.0e-1_PR*rho0
   rhomax=1.0e1_PR*rho0
   Nr=1000
   rho=rhomin
   P=P_tot(Z,rho,T)
   N=0
   r1=0.0_PR ; r2=0.0_PR ; P1=0.0_PR ; P2=0.0_PR

   do ir=1,Nr
     rhom1=rho ; Pm1=P
     rho=rhomin+real(ir,PR)*(rhomax-rhomin)/real(Nr,PR)
     P=P_tot(Z,rho,T)
     if(Pm1.le.P0.and.P.ge.P0.or.&
        Pm1.ge.P0.and.P.le.P0)then
       N=N+1
       r1(N)=rhom1
       r2(N)=rho
       P1(N)=Pm1
       P2(N)=P
     endif     

   enddo

   if(N.eq.1)then

      rhomin=r1(1)
      rhomax=r2(1)
      dPdr=(P2(1)-P1(1))/(r2(1)-r1(1))

   elseif(N.gt.1)then

   !!!---selection du plus proche selon rho
      drho=0.0_PR
      do i=1,N
       drho(i)=min(abs(rho0-r1(i)),abs(r2(i)-rho0))
      enddo
      i=minloc(drho(1:N),1)
      rhomin=r1(i)
      rhomax=r2(i)
      dPdr=(P2(i)-P1(i))/(r2(i)-r1(i))

   elseif(N.eq.0)then

     print*, 'PROBLEM rhominmax ! :P '
     stop

   endif

end subroutine rhominmax_TP


real(PR) function rho_TP(Z,T,P0,rho0)
   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: T, P0, rho0
   real(PR) :: rhomin, rhomax, dPdr
   real(PR) :: P, rho
   logical :: problem, found
   integer :: it

   !!! V1: -----------------------------

   !rhomin=rho0*0.2_PR
   !rhomax=rho0*5.0_PR

   !!! V2 : -----------------------------

   call rhominmax_TP(Z,T,P0,rho0,dPdr,rhomin,rhomax)

   !!!----------------------------


   problem=.false.
   found=.false. 
   it=0

   do while(.not.found.and..not.problem)

      rho=0.5_PR*(rhomin+rhomax)
      P=P_tot(Z,rho,T)

      if(P.gt.P0)then
        if(dPdr.ge.0.0_PR)then
          rhomax=rho
        else
          rhomin=rho
        endif
      else
        if(dPdr.ge.0.0_PR)then
          rhomin=rho
        else
          rhomax=rho
        endif
      endif

      if(abs(P-P0)/P0.lt.1.0e-4_PR)then
        found=.true.
      endif

      it=it+1

      if(it.gt.9989)then
        print*, 'it=', it, 'rho=', rho, 'p=', p , 'vs', P0
      endif

      if(it.gt.10000)then
        problem=.true.
      endif

   enddo

   if(problem)then
      print*, ' PROBLEM IN rho_TP'
      stop
   else
      rho_TP=rho
   endif

end function rho_TP



!!!=============================  CONDUCTIVITY MODEL  ====================================
!!!
!!! Lee More 1984 + Desjarlais 2001
!!!
!!!=======================================================================================


subroutine transport_properties(Z,rho,T,Zbar,sig,K,S,&
                                sig1,sig2,K1,K2,na,ne,nn,ni,&
                                lnl,tau,tau1,tau2,tau_en,cross_en,mu)

    use mod_Fermi_Dirac, only : fd1h, fdm1h

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho, T
   real(PR), intent(in) :: Zbar
   real(PR), intent(out) :: sig,K,S,sig1,sig2,K1,K2,lnl,tau
   real(PR), intent(out) :: na, nn, ni, ne !!! density of atoms, neutrals, ions, electrons
   real(PR) :: F12, x, ln_Lambda, Te, Ti, Zi
   real(PR) :: Tmelt, R0, bmin, bmax, lambda_dense
   real(PR), intent(out) :: tau1, tau2, tau_en, cross_en
   real(PR) :: VF, vth, TF, EF, mu
   real(PR) :: lb, ldbc, ldbg, lcca
   !!!---for Redmer e-n model:
   real(PR) :: alphaD,lambda_e,Fm12,kappa,kwave_e,kar,kr0,Ak,Bk,Ck,Dk,Ek

   Te=T ; Ti=T

   na=rho/A(Z)*Avogadro*1.0e6_PR !!! m^-3

   ne=na*Zbar

   if(Zbar.ge.1.0_PR)then
     ni=na
     nn=0.0_PR
     Zi=Zbar
   else
     ni=ne
     nn=na-ni
     Zi=1.0_PR
   endif

   EF=hb**2/(2.0_PR*meSI)*(3.0_PR*Pi**2*ne)**(2.0_PR/3.0_PR) !!!J

   !mu=EF

   !print*, 'ne=', ne, 'Zbar=', Zbar, 'na=', na, 'Te=', Te, 'EF=', EF

   mu=chemical_potential(EF,Te)

   !print*, 'mu=', mu, 'Zbar=', Zbar

   TF=EF/kb

   x=mu/(kb*Te)
 
   F12=fd1h(x)

   !bmax=sqrt(1.0_PR/(&
   !         4.0_PR*Pi*ne*qe**2/(kb*sqrt(Te**2+TF**2))&
   !        +4.0_PR*Pi*ni*Zbar**2*qe**2/(kb*Ti)&
   !        ))

   !bmax=sqrt(1.0_PR/(&
   !         ne*qe**2/(eps0*kb*sqrt(Te**2+TF**2))&
   !        +ni*Zbar**2*qe**2/(eps0*kb*Ti)&
   !        ))

   bmax=sqrt(1.0_PR/(&
            ne*qe**2/(eps0*kb*sqrt(Te**2+TF**2))&
           +ni*Zi**2*qe**2/(eps0*kb*Ti)&
           ))


   !bmax=sqrt((eps0*kb*Te)/(ne*qe**2))


   !!!---electron velocity 
   vF=sqrt(2.0_PR*EF/meSI)
   !vF=v_Fermi(Z)

   vth=max(sqrt(3.0_PR*kB*Te/meSI),vF)

   !vth=max(sqrt(3.0_PR*kB*Te/meSI),0.0_PR)


   !bmin=max( Zbar*qe**2/( 4.0_PR*Pi*eps0*meSI*vth**2 ), h/(2.0_PR*meSI*vth))

   bmin=max( Zi*qe**2/( 4.0_PR*Pi*eps0*meSI*vth**2 ), h/(2.0_PR*meSI*vth))

   !bmin=sqrt( (Zbar*qe**2/( 4.0_PR*Pi*eps0*3.0_PR*kb*Te))**2 + (hb/(2.0_PR*sqrt(3.0_PR*me*kb*Te)))**2 )


   !!!---papier de Lee More 1984
   ln_lambda=0.5_PR*log(1.0_PR+bmax**2/bmin**2)

   !!!---Hayes 2016 (sur l'implementation de 84)
   !ln_lambda=0.5_PR*log(exp(1.0_PR)+bmax**2/bmin**2)

   !!!---Hayes 2016 (sur l'implementation de 84)
   !lb=exp(1.0_PR)
   !ldbc=sqrt((eps0*kb*Te)/(ne*qe**2))
   !ldbg=hb/(2.0_PR*sqrt(3.0_PR*meSI*kb*Te))
   !Lcca=1.0_PR/(4.0_PR*Pi*eps0)*Zbar*qe**2/(3.0_PR*kb*Te)

   !ln_lambda=0.5_PR*log(lb**2+Ldbc**2/(Ldbg**2+Lcca**2))

   !ln_lambda=max(0.5_PR*log(lb**2+Ldbc**2/(Ldbg**2+Lcca**2)),2.0_PR)


   lnl=ln_lambda

   !!!====Coulomb cross section for e-i with Debye-Hückel cut-off
   !!!    and degeneracy correction

  ! tau1=( (3.0_PR*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
  !        (2.0_PR*sqrt(2.0_PR)*Pi*Zbar**2*ni*qe**4*ln_Lambda) )*&
  !        ( 1.0_PR+exp(-x) )*F12

   !!!--version 0 de Lee et More
   !tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
   !       (2.0_PR*sqrt(2.0_PR)*Pi*Zi**2*ni*qe**4*max(ln_Lambda,2.0_PR)) )*&
   !       ( 1.0_PR + exp(-x) )*F12

   !!!--version de Nanagan 2015 avec correction pour e-e scattering
   tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
          (2.0_PR*sqrt(2.0_PR)*Pi*Zi**2*ni*qe**4*max(ln_Lambda,2.0_PR)) )*&
          Fc_alpha(Zi)*( 1.0_PR + exp(-x) )*F12

  


   !tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
   !       (2.0_PR*sqrt(2.0_PR)*Pi*Zbar**2*na*qe**4*max(ln_Lambda,2.0_PR)) )*&
   !       ( 1.0_PR + exp(-x) )*F12


   !tau=(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR)/&
   !    (Pi*Zbar**2*ni*qe**4*ln_lambda)

   !!!===Sphere dure ?

   if(nn.gt.0.0_PR)then

      !!!---Lee More (2e-15 cm^2 from desjarlais 2001)--------------------
      !cross_en=2.0e-19_PR
      !!!---Max value from Desjarlais (30e-15 cm^2 from desjarlais 2001)--------------------
      !cross_en=30.0e-19_PR
      !!! Desjarlais 2001: model from Redmer 1999 -------------------
      alphaD=pola(Z)
      lambda_e=sqrt(2.0_PR*Pi*hb**2/(kb*Te*meSI)) !!! thermal wavelength (cf eq 4)
      Fm12=fdm1h(x)     
      !!! eq 14 in Redmer1999
      kappa=sqrt(qe**2/(kb*Te*eps0)*2.0_PR/lambda_e**3*Fm12)     
      kwave_e=meSI*vth/hb !!! electron wave number Desjarlais 2001
      r0=(alphaD*RBohr/(2.0_PR*real(Z,PR)**(1.0_PR/3.0_PR)))**0.25_PR
      kar=kappa*r0
      kr0=kwave_e*r0

      Ak=1.0_PR+2.0_PR*kar+7.0_PR/pi**2*kar**2+pi/7.0_PR*kar**3
      Bk=exp(-18.0_PR*kar)
      Ck=(1.0_PR+22.0_PR*kar-11.3_PR*kar**2+33.0_PR*kar**4)/&
         (1.0_PR+6.0_PR*kar+4.7_PR*kar**2+2.0_PR*kar**4)
      Dk=(1.0_PR+28.0_PR*kar+13.8_PR*kar**2+3.2_PR*kar**3)/&
         (1.0_PR+8.0_PR*kar+10.0_PR*kar**2+kar**3)
      Ek=1.0_PR+0.1_PR*kar+0.3665_PR*kar**2

      cross_en=pi**3*(alphaD/(2.0_PR*r0*RBohr))**2/&
      (Ak**2+3.0_PR*Bk*kr0+7.5_PR*Ck*kr0**2-3.4_PR*Dk*kr0**3+10.6668_PR*Ek*kr0**4)

      !print*,  'cross_en=', T, cross_en, 'r0**2', r0**2
      !!!--------------------------------------------------------------
         
      tau_en=1.0_PR/(nn*vth*cross_en)
      tau=1.0_PR/(1.0_PR/tau1+1.0_PR/tau_en)
      !tau =tau1    
 
   else

      tau=tau1

   endif


   !!!====Bloch-Grüneisen - Zimann formula for solid-liquid conductivities
  
   Tmelt=Tm(Z,rho)
   R0=(1.0_PR/na)**(1.0/3.0_PR)

   if(T.le.Tmelt)then
     lambda_dense=50.0_PR*R0*(Tmelt/T)
   else
     lambda_dense=50.0_PR*R0*(Tmelt/T)/sigSL(Z)
   endif

   !tau2=max(lambda_dense/vF,R0/vF)
   tau2=lambda_dense/vF

   !taumin=min(p2a+p2b/
   !tau=max(tau,taumin)


   !!!===Final tau

   !tau=tau1
   !write(6,"(4(ES14.7,' '))"), ne, tau,A_alpha(x), x
 
   sig1=ne*qe**2*tau/meSI*A_alpha(x)
   sig2=ne*qe**2*tau2/meSI
   !!!---scaling of sig2 conductivity for solids
   sig2=sig2*sig0(Z)/sigmax(Z)


   sig=sig1+sig2   !!! max(sig1,sig2)
   sig=max(sig,0.001_PR)

   K1=ne*kb*(kb*Te)*tau/meSI*A_beta(x)
   K2=ne*kb*(kb*Te)*tau2/meSI*A_beta(x)
   !!!---scaling of K2 with reference value
   K2=K2*lamb0(Z)/lambmax(Z)

   K=K1+K2
   K=max(K,0.001_PR)

   S=(kb/qe)*A_gamma(x)

   !sigS=sig1
   !sigM=sig2

   !sigS=sig1
   !sigM=sig2


   !tau=( 3.0_PR/4.0_PR*sqrt(meSI/(2.0_PR*Pi))*(kb*Te)**(1.5_PR) )/&
   !    (Zbar**2*ni*qe**4*ln_lambda)


   !!!---Spitzer
   !bmax=sqrt(eps0*kb*Te/(ne*qe**2*(1.0_PR+Zbar)))
   !vth=sqrt(3.0_PR*kB*Te/meSI)
   !bmin=1.0_PR/(4.0_PR*Pi*eps0)*Zbar*qe**2/(meSI*vth**2)
   !ln_lambda=log(bmax/bmin)

   !bmin=Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2)
!   bmin=max(Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2), h/(2.0_PR*meSI*vth))
   !bmin=max(Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2), h/(2.0_PR*meSI*vth))

  ! ln_lambda=0.5_PR*log(1.0_PR+bmax**2/bmin**2)

   !tau=(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR)/&
   !    (Pi*Zbar**2*ni*qe**4*ln_lambda)

   !tau=3.0_PR/4.0_PR*sqrt(meSI)/sqrt(2.0_PR*pi)*(4.0_PR*Pi*eps0)**2*(kb*Te)**(1.5_PR)/(Zbar**2*ni*qe**4*ln_lambda)

   !sigS=ne*qe**2*tau/meSI*32.0_PR/(3.0_PR*Pi)


   !tau2=3.0_PR*Pi*hb**3/(4.0_PR*meSI*Zbar*qe**4*ln_lambda)


   !sigM=ne*qe**2*tau2/meSI

   !!!!!print*, 'ne=', ne, 'Zb=', Zbar, 'tau=', tau, 'lnL=', ln_lambda, 'tau2=', tau2
   !print*, 'ne=', ne, 'Zb=', Zbar, 'lnL=', ln_lambda, 'bmin/max=', bmin, bmax

end subroutine transport_properties


subroutine compute_ZFS(Z,rho,T,ZTF,ZS,ZFS)

   !!!---modification de Z par la formule (3) de Desjarlais 2001

   implicit none
   integer, intent(in) :: Z
   real(PR), intent(in) :: rho, T
   real(PR), intent(in) :: ZTF
   real(PR), intent(out) :: ZFS,ZS
   real(PR) :: na, Io, g0, g1, Ra, K, fe

   na=rho/A(Z)*Avogadro*1.0e6_PR !!! m^-3
 
   Io=Eio1(Z)*qe
   g0=2.0_PR
   g1=1.0_PR
   Ra=(3.0_PR/(4.0_PR*Pi*na))**(1.0_PR/3.0_PR)

   !!!-----Version Originale

   !K=2.0_PR*g1/g0*1.0_PR/na*(2.0_PR*Pi*me*kb*T/h**2)**1.5_PR*&
   !  exp(-Io/(kb*T)*(1.0_PR-min(0.0_PR,(1.5_PR*qe**2/(Io*4.0_PR*pi*eps0*Ra))**1.5_PR)) )

   !fe=0.5_PR*(sqrt(K**2+4.0_PR*K)-K)

   !!!!!!fe=max(fe,1.0e-10_PR)


   !ZFS=fe**(2.0_PR/ZTF**2)*ZTF+(1.0_PR-fe**(2.0_PR/ZTF**2))*fe 

   !ZS=fe

   !!!-----Version Fabien
 
   K=2.0_PR*g1/g0*(2.0_PR*Pi*me*kb*T/h**2)**1.5_PR*&
     exp(-Io/(kb*T)*(1.0_PR-min(1.0_PR,(1.5_PR*qe**2/(Io*4.0_PR*pi*eps0*Ra))**1.5_PR)) )

   fe=0.5_PR*(sqrt(K**2+4.0_PR*K*na)-K)/na
  
   ZS=fe
 
   ZFS=fe**(2.0_PR/ZTF**2)*ZTF+(1.0_PR-fe**(2.0_PR/ZTF**2))*fe

   !print*, 'fe=', fe, 'K=', K, 'ZTF=', ZTF, 'ZFS=', ZFS

end subroutine compute_ZFS

real(PR) function chemical_potential(EF,T)

   implicit none
   real(PR), intent(in) :: EF, T
   real(PR) :: Rm, denom, a(1:3), b(1:4), xi, chi
   integer :: i


   !!!---fit Rm3 de Managan 2015
   !a(1)=0.19972_PR
   !a(2)=0.17258_PR
   !a(3)=0.145_PR

   !b(1)=0.25829_PR
   !b(2)=0.28756_PR
   !b(3)=0.16842_PR
   !b(4)=0.145_PR

   !xi=sqrt(EF/(kb*T))

   !Rm=4.0_PR/sqrt(3.0_PR) 
   !do i=1,3
   !   Rm=Rm+a(i)*xi**i
   !enddo

   !denom=1.0_PR
   !do i=1,4   
   !   denom=denom+b(i)*xi**i
   !enddo

   !chi=Rm/denom*xi**3

   !chemical_potential=kb*T*log( exp(chi)-1.0_PR )

   !!!---Zimmerman's form (see Managan 2015) :

   xi=sqrt(EF/(kb*T))

   chi=(0.7531_PR+0.1679_PR*xi+0.3108*xi**2)/&
       (1.0_PR+0.2676_PR*xi+0.2280*xi**2+0.3099*xi**3)*xi**3

   !write(6,*) 'xi=', xi, 'chi=', chi
   !write(6,*) 'exp(chi)', exp(chi)


   chemical_potential=kb*T*log( max(exp(min(chi,100.0_PR))-1.0_PR,1.0e-100_PR) )

   if(chemical_potential.ne.chemical_potential)then
     print*, 'PROBLEM CHEMICAL POTENTIAL'
     stop
   endif


end function chemical_potential

real(PR) function Fc_alpha(Z)

    implicit none
    real(PR), intent(in) :: Z
    real(PR) :: x

    x=1.0_PR/Z

    Fc_alpha=0.295_PR/(1.0_PR-&
    (0.0678_PR+0.4924_PR*x+0.9760_PR*x**2+0.3008_PR*x**3)/&
    (0.0961_PR+0.7778_PR*x+1.5956_PR*x**2+1.3008_PR*x**3) )

end function Fc_alpha



real(PR) function A_alpha(x)

    use mod_Fermi_Dirac, only : fd1h, fd6h, fd4h

    implicit none
    real(PR), intent(in) :: x
    real(PR) :: F12, F3, F2

    F12=fd1h(x)
    F2=fd4h(x)
    !F3=fd6h(x)

    !A_alpha=4.0_PR/3.0_PR*F3/( (1.0_PR+exp(-1.0_PR*x))*F12**2 )

    A_alpha=4.0_PR/3.0_PR*F2/( (1.0_PR+exp(-x))*F12**2 )

end function A_alpha


real(PR) function A_alpha2(x)

    implicit none
    real(PR), intent(in) :: x
    real(PR) :: F12, F3, a(1:3), b(1:3), denom, y
    integer :: i

     y=log(1.0_PR+exp(x))
     a(1)=3.39_PR
     a(2)=0.347_PR
     a(3)=0.129_PR
     b(1)=0.0_PR
     b(2)=0.511_PR
     b(3)=0.124_PR

     A_alpha2=0.0_PR 
     do i=1,3
        A_alpha2=A_alpha2+a(i)*y**(i-1)
     enddo

     denom=1.0_PR
     do i=2,3
        denom=denom+b(i)*y**(i-1)
     enddo

     A_alpha2=A_alpha2/denom

end function A_alpha2

real(PR) function A_beta(x)

    use mod_Fermi_Dirac, only : fd1h, fd4h, fd6h, fd8h

    implicit none
    real(PR), intent(in) :: x
    real(PR) :: F12, F2, F3, F4

    F12=fd1h(x) 
    F2=fd4h(x) 
    F3=fd6h(x) 
    F4=fd8h(x)

    A_beta=20.0_PR/9.0_PR*F4*(1.0_PR-16.0_PR*F3**2/(15.0_PR*F4*F2))/&
          ( (1.0_PR+exp(-x))*F12**2 )

end function A_beta

real(PR) function A_gamma(x)

    use mod_Fermi_Dirac, only : fd1h, fd3h, fd4h, fd6h

    implicit none
    real(PR), intent(in) :: x
    real(PR) :: F12, F32, F2, F3

    F12=fd1h(x) 
    F32=fd3h(x)
    F2=fd4h(x) 
    F3=fd6h(x) 

    A_gamma=5.0_PR/3.0_PR*F32/F12-4.0_PR/3.0_PR*F3/F2

end function A_gamma

real(PR) function mfromn(Z,n,facB0)

    use mod_Fermi_Dirac, only : fd1h, fd3h, fd4h, fd6h

    implicit none
    integer, intent(in) :: Z
    real(PR), intent(in) :: n, facB0
    real(PR) :: Ecoh, rs, B0

    Ecoh=Ecoh_J(Z)
    rs=rhos(Z)
    B0=bulk_modulus(Z)*1.0e3_PR

    mfromn=n+(facB0*B0-rs*n**2*Ecoh)/(rs*n*Ecoh)


end function mfromn





end module mod_EOS
