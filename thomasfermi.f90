!> 
!! @file   thomasfermi.f90
!! @brief  Routines for calculating the average charge of an atom using the semi-empirical Thomas-Fermi model.
!!
!! @author Thad A. Heltemes
!! @date   23 March 2010
!!

!> Contains the data and functions necessary to calculate a Thomas-Fermi ionization state
module thomasfermi


   !use parameters, only: kind_int, ik, kind_dbl, dk

   implicit none

   integer, parameter :: kind_dbl=selected_real_kind(8)
   integer, parameter :: dk=selected_real_kind(8)


   real(kind_dbl) :: tf_coeff(9) !< The fitting coefficients
   real(kind_dbl) :: tf_alpha, tf_beta
   real(dk) ::  avogadros_number=6.0221409e23_dk

   contains

   !> Initialize the fitting coefficient data for the TF ionization model
   subroutine init_thomas_fermi_ionization_data()

      implicit none

      tf_alpha = 14.3139_dk
      tf_beta  = 0.6624_dk

      ! These are a_1 - a_4 from More's paper
      tf_coeff(1) = 0.003323_dk
      tf_coeff(2) = 0.9718_dk
      tf_coeff(3) = 0.0000926148_dk
      tf_coeff(4) = 3.10165_dk
      ! These are b_0 - b_2 from More's paper
      tf_coeff(5) = -1.7630_dk
      tf_coeff(6) = 1.43175_dk
      tf_coeff(7) = 0.31546_dk
      ! These are c_1, c_2 from More's paper
      tf_coeff(8) = -0.3666666666666667_dk
      tf_coeff(9) = 0.9833333333333334_dk

   end subroutine

   !> Calculate the Zbar of an atom based on the Thomas-Fermi model
   function thomas_fermi_zbar(ndens,temp,prot_num)

      !use constants, only: avogadros_number

      implicit none 

      real(kind_dbl) :: temp              !< gas temperature, eV
      real(kind_dbl) :: ndens             !< gas number density, cm^-3
      real(kind_dbl) :: prot_num          !< proton number (Z) of the gas
      real(kind_dbl) :: thomas_fermi_zbar !< return value, Zbar

      real(kind_dbl) :: rho1, tau1, tauf, a_fit, b_fit, c_fit, q1, fcn, x 

      rho1 = ndens/(prot_num*avogadros_number)

      if (temp==0._dk) then
         fcn = rho1;
      else
         tau1 = temp*(prot_num**(-4._dk/3._dk))
         tauf = tau1/(1._dk+tau1)
         a_fit = tf_coeff(1)*(tau1**tf_coeff(2))+tf_coeff(3)*(tau1**tf_coeff(4))
         b_fit = -dexp(tf_coeff(5)+tf_coeff(6)*tauf+tf_coeff(7)*(tauf**7._dk))
         c_fit = tf_coeff(8)*tauf+tf_coeff(9)
         q1 = a_fit*(rho1**b_fit)
         fcn = (rho1**c_fit + q1**c_fit)**(1._dk/c_fit)
      endif

      x = tf_alpha*(fcn**tf_beta)
      thomas_fermi_zbar = (prot_num*x)/(1._dk+x+dsqrt(1._dk+2._dk*x))

   end function thomas_fermi_zbar

!   !> Update the Zbar using the TF model where appropriate
!   subroutine update_tf_zbar_data()
!
!      use materials, only: mix, mat, shm_k0max, shm_fill, shm_dk0
!      use switches, only: use_busquet_nlte
!      use carlson_ionization_data, only: carlson_principal, carlson_angular
!
!      implicit none
!
!      integer(kind_int) :: i, j, imax, jmax, ii, jj
!      integer(kind_int), allocatable :: dko(:)
!      real(kind_dbl) :: tf_zbar, tfztp1, tfztm1, delta, te, beta_i
!      
!      imax = size(mat,1)
!      jmax = size(mat,2)
!   
!      ! Do the Ideal Gas calculations
!      j = 1_ik
!      do while ((mat(1,j)%z > 0_ik) .and. (j <= jmax))
!         i = 1_ik
!         do while ((mat(i,j)%z > 0_ik) .and. (i <= imax))
!            ! Check for a Thomas Fermi Zbar model
!            if (mix(j)%zmodl.eq.0_ik) then
!               allocate(dko(shm_k0max(mat(i,j)%z)))
!               ! Compute the ionization using the TF function
!
!               if (use_busquet_nlte) then
!                  beta_i = 1.34e13_dk*(mix(j)%te**3.5_dk)/(mix(j)%n*mat(i,j)%f);
!                  te = mix(j)%te/((1._dk+0.25_dk*beta_i)**0.19_dk);
!               else
!                  te = mix(j)%te
!               end if
!
!               tf_zbar = thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, te, dble(mat(i,j)%z) )
!
!               ! Now we need to correctly populate the Pj (electron) and Pk
!               ! (shell) arrays with the reduced bound electrons (occupation fractions)
!               mat(i,j)%zbar%pk = 0._dk
!               mat(i,j)%zbar%pj = 0._dk
!               mat(i,j)%zbar%pj(1:mat(i,j)%z) = 1._dk - tf_zbar/dble(mat(i,j)%z)
!
!               ! For the electron EOS calculations, we also need dPj/dt
!               delta = te*1.e-4_dk
!               tfztp1 = 1._dk - thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, te+delta, dble(mat(i,j)%z) )/dble(mat(i,j)%z)
!               tfztm1 = 1._dk - thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, te-delta, dble(mat(i,j)%z) )/dble(mat(i,j)%z)
!               mat(i,j)%zbar%dpdt(1:mat(i,j)%z) = (tfztp1-tfztm1)/(2._dk*delta)
!
!               ! And now the zero temperature data
!               mat(i,j)%zbar%pjo(1:mat(i,j)%z) = 1._dk - thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, 1.e-6_dk, dble(mat(i,j)%z) )/&
!                  dble(mat(i,j)%z)
!               mat(i,j)%zbar%dj(1:mat(i,j)%z) = mat(i,j)%zbar%pjo(1:mat(i,j)%z)
!               tfztp1 = 1._dk - thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, 1.e-6_dk+1.e-10_dk, dble(mat(i,j)%z) )/dble(mat(i,j)%z)
!               tfztm1 = 1._dk - thomas_fermi_zbar( mix(j)%n*mat(i,j)%f, 1.e-6_dk-1.e-10_dk, dble(mat(i,j)%z) )/dble(mat(i,j)%z)
!               mat(i,j)%zbar%dpdto(1:mat(i,j)%z) = (tfztp1-tfztm1)/(2._dk*delta)
!
!               dko = shm_dk0(1:shm_k0max(mat(i,j)%z))
!               dko(shm_k0max(mat(i,j)%z)) = shm_fill(mat(i,j)%z)
!               ii = shm_k0max(mat(i,j)%z)
!               do jj = 1_ik, mat(i,j)%z
!                  mat(i,j)%zbar%pk(ii) = mat(i,j)%zbar%pk(ii) + mat(i,j)%zbar%pj(jj)
!                  if (dko(ii) .eq. 1_ik) then
!                     ii = ii - 1_ik
!                  else
!                     dko(ii) = dko(ii) - 1_ik
!                  end if
!               end do
!               deallocate(dko)
!            end if
!            i = i + 1_ik
!         end do
!         j = j + 1_ik
!      end do
!
!   end subroutine

end module thomasfermi

