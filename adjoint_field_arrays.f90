module adjoint_field_arrays
  
  ! public :: phi_save
  ! public :: chi_save
  ! public :: unpert_q, pert_q
  ! public :: unpert_term, denominator, pert_term
  ! public :: omega_g, omega
  ! public :: derivative

  implicit none
  
  complex, dimension (:,:,:,:), allocatable :: phi_save, chi_save

  complex, dimension (:,:,:,:), allocatable :: unpert_q
  complex, dimension (:,:,:,:), allocatable :: pert_q

  complex, dimension (:,:), allocatable :: unpert_term, denominator
  complex, dimension (:,:), allocatable :: pert_term

  complex, dimension (:,:), allocatable :: omega_g, omega
  complex, dimension (:,:,:), allocatable :: derivative
  
end module adjoint_field_arrays
