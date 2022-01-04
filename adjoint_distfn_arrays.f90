module adjoint_distfn_arrays

  public :: unpert_l
  public :: pert_l
  public :: g_omega
  public :: gsave
  public :: lam_save

  complex, dimension (:,:,:,:,:), allocatable :: g_omega
  
  complex, dimension (:,:,:,:,:), allocatable :: unpert_l
  complex, dimension (:,:,:,:,:), allocatable :: pert_l
  
  complex, dimension (:,:,:,:,:), allocatable :: gsave
  complex, dimension (:,:,:,:,:), allocatable :: lam_save
  
end module adjoint_distfn_arrays
  
