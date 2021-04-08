module mod_data

implicit none

integer, parameter, public :: PR=selected_real_kind(8)
real(PR), parameter, public :: Pi=2.0_PR*Asin(1.0_PR)
integer :: iout=6

integer :: mode_output=2

real(PR) :: Rm=8.3144621_PR

integer :: Z
real(PR) :: Temp, Ener

real(PR) :: rpar(1:10)=0.0_PR
integer :: ipar(1:10)=0

!type material0 
!   integer :: Z
!   real(PR) :: M 
!   real(PR) :: Rs
!   real(PR) :: gl
!   real(PR) :: pi
!   real(PR) :: mu
!   real(PR) :: sigy
!   real(PR) :: rhol
!   real(PR) :: Cvl
!   real(PR) :: Cpl
!   real(PR) :: gg
!   real(PR) :: Cvg
!   real(PR) :: Cpg
!   real(PR) :: bl
!   real(PR) :: bg
!   real(PR) :: ql, qg, qpl, qpg
!   real(PR) :: A, B, C, D, E
!   real(PR) :: Tc, Pc, Tvap0, Pvap0, Tvap1, Pvap1
!end type material0

contains


!subroutine init_mat
!
!implicit none
!
!   Cu%Z=29
!   Cu%M=63.546e-3_PR
!   Cu%Rs=Rm/Cu%M
!   Cu%gl=4.22_PR
!   Cu%pi=34.2e9_PR
!   Cu%mu=9.2e10_PR
!   Cu%sigy=40.0e6_PR
!   Cu%rhol=8900.0_PR
!   Cu%Cpl=380.0_PR !!! J/kg/K
!   Cu%Cvl=Cu%Cpl/Cu%gl !!! J/kg/K
!   Cu%gg=1.667_PR
!   Cu%Cvg=Cu%Rs/(Cu%gg-1.0_PR) !!! J/kg/K
!   Cu%Cpg=Cu%gg*Cu%Cvg !!! J/kg/K
!   Cu%bl=0.0_PR
!   Cu%bg=0.0_PR
!   Cu%qpl=0.0_PR
!   Cu%qpg=0.0_PR !!! ? ref entropie of gas
!   Cu%ql=0.0_PR  !!! ? bounding energy of liquid
!   Cu%qg=0.0_PR
!   Cu%A=(Cu%Cpl-Cu%Cpg+Cu%qpg-Cu%qpl)/(Cu%Cpg-Cu%Cvg)
!   Cu%B=(Cu%ql-Cu%qg)/(Cu%Cpg-Cu%Cvg)
!   Cu%C=(Cu%Cpg-Cu%Cpl)/(Cu%Cpg-Cu%Cvg)
!   Cu%D=(Cu%Cpl-Cu%Cvl)/(Cu%Cpg-Cu%Cvg)
!   Cu%E=(Cu%bl-Cu%bg)/(Cu%Cpg-Cu%Cvg)
!
!   Cu%Pc=6510.0e5_PR
!   Cu%Tc=8440.0_PR
!   Cu%Tvap0=2835.0_PR
!   Cu%Pvap0=1.0e5_PR
!   Cu%Tvap1=1357.0_PR
!   Cu%Pvap1=0.05_PR
!
!
!
!end subroutine init_mat



end module mod_data

