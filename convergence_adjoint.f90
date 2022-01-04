module convergence_adjoint

  implicit none

  public :: init_convergence
  public :: omega_convergence1
  public :: omega_convergence2
  public :: adjoint_convergence
  public :: deallocate_convergence
  public :: track, track_short
  
  private
  
  real :: halfnavg
  complex, dimension (:,:,:), allocatable :: track, track_short
  
contains
  
  subroutine init_convergence
    
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_diagnostics, only: navg
    use stella_layouts, only: vmu_lo
    
    use adjoint_field_arrays, only: omega_g, omega
    
    implicit none
    
    real :: halfnavg

    halfnavg = navg/2

    if(.not. allocated(omega)) then 
       allocate(omega(naky,nakx))
       omega = 0.
    end if
    if (.not. allocated(omega_g)) then
       allocate(omega_g(naky,nakx))
       omega_g = 0.
    end if

    if(.not. allocated(track)) then
       allocate(track(navg,naky,nakx))
       track = 0.
    end if
    if(.not. allocated(track_short)) then
       allocate(track_short(nint(halfnavg),naky,nakx))
       track_short = 0.
    end if
  end subroutine init_convergence
  
  subroutine omega_convergence1 (istep, converged, adjoint)

    use kt_grids, only: nakx, naky
    use stella_time, only: code_time
    use stella_diagnostics, only: omega_vs_time, navg

    implicit none
    
    complex, dimension(:,:), allocatable :: diff
    integer, intent (in) :: istep
    integer :: ikx, iky
    real :: max_diff
    logical, intent (out) :: converged 
    logical, intent (in), optional :: adjoint
   
    converged = .False.
    
    if(istep > navg) then
       allocate(diff(naky,nakx))
       do ikx = 1, nakx
          do iky = 1, naky
             if (mod(istep,navg)+1 > 1) then
                diff(iky, ikx) = omega_vs_time(mod(istep,navg)+1,iky,ikx)- omega_vs_time(mod(istep,navg),iky,ikx)
             else
                diff(iky,ikx) = omega_vs_time(1,iky,ikx)- omega_vs_time(navg,iky,ikx)
             end if
          end do
       end do
       
       max_diff = maxval(abs(aimag(diff)))
       if(max_diff < 1E-004) converged = .True.
      
       deallocate(diff)
    end if

  end subroutine omega_convergence1

  !!Second Convergence Test

  subroutine omega_convergence2 (istep, istep_initial,  converged, adjoint)

    use stella_diagnostics, only: omega_vs_time, navg, omega_vs_time_short
    use kt_grids, only: nakx, naky
    use adjoint_field_arrays, only: omega

    implicit none
    
    integer, intent (in) :: istep, istep_initial
    
    complex, dimension (:,:), allocatable :: sum_omega, avg_omega
    complex, dimension (:,:), allocatable :: sum_omega_local, avg_omega_local
    complex, dimension(:,:), allocatable :: diff_omega
    real :: max_diff

    logical, intent (out) :: converged
    logical, intent (in), optional :: adjoint
    
    if (.not. allocated(sum_omega)) allocate(sum_omega(naky,nakx))
    allocate(sum_omega_local(naky,nakx))
    allocate(avg_omega(naky,nakx))
    allocate(avg_omega_local(naky,nakx))
    allocate(diff_omega(naky,nakx))
    
    converged = .False.

    if ((istep-istep_initial) > navg) then
       halfnavg = navg/2
       sum_omega = sum(omega_vs_time, dim=1)
       avg_omega = sum_omega/navg
       sum_omega_local = sum(omega_vs_time_short, dim=1)
       avg_omega_local = sum_omega_local/(nint(halfnavg))
       diff_omega = avg_omega - avg_omega_local       
       max_diff = maxval(abs(aimag(diff_omega)))
       if (max_diff < 1E-003) then
          converged = .True.
          write(*,*) 'CONVERGED!'
       end if
    end if
       
    deallocate(avg_omega,avg_omega_local,sum_omega,sum_omega_local,diff_omega)
    
  end subroutine omega_convergence2

  subroutine adjoint_convergence (istep, converged)

    use dist_fn_arrays, only: gnew
    use adjoint_distfn_arrays, only: g_omega
    use adjoint_field_arrays, only: omega
    
    use stella_time, only: code_dt

    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old

    use stella_diagnostics, only: navg
    use kt_grids, only: nakx, naky

    use fields, only: advance_fields, fields_updated
    use volume_averages, only: fieldline_average
    implicit none

    integer, intent (in) :: istep
    real, dimension (:,:), allocatable :: sum_omega, avg_omega
    real, dimension (:,:), allocatable :: sum_omega_local, avg_omega_local
    real, dimension (:,:), allocatable :: diff_omega

    complex, dimension (:,:), allocatable :: phiavg, phioldavg

    real :: max_diff

    logical, intent (out) :: converged
    logical :: adjoint = .true.
    real :: zero

    allocate(sum_omega(naky,nakx))
    allocate(sum_omega_local(naky,nakx))
    allocate(avg_omega(naky,nakx))
    allocate(avg_omega_local(naky,nakx))
    allocate(diff_omega(naky,nakx))

    allocate (phiavg(naky,nakx))
    allocate (phioldavg(naky,nakx))

    converged = .false.
    halfnavg = navg/2
    zero = 100.*epsilon(0.)

    fields_updated = .false.
    call advance_fields (gnew, phi, apar, dist='gbar', adjoint=adjoint)
    fields_updated = .false.
    call advance_fields (g_omega, phi_old, apar, dist='gbar', adjoint=adjoint)

    call fieldline_average (phi, phiavg)
    call fieldline_average (phi_old, phioldavg)

    where (abs(phiavg) < zero .or. abs(phioldavg) < zero)
       track(mod(istep,navg)+1,:,:) = 0.0
       track_short(mod(istep,nint(halfnavg))+1,:,:) = 0.0
    elsewhere
       track(mod(istep,navg)+1,:,:) = log(phiavg/phioldavg)/code_dt
       track_short(mod(istep,nint(halfnavg))+1,:,:) = log(phiavg/phioldavg)/code_dt
    end where

    if (istep> navg) then
       sum_omega = sum(real(track), dim=1)
       avg_omega = sum_omega/navg
       sum_omega_local = sum(real(track_short), dim=1)
       avg_omega_local = sum_omega_local/(nint(halfnavg))
       diff_omega = avg_omega - avg_omega_local

       max_diff = maxval(abs(diff_omega))

       omega = track_short(mod(istep,nint(halfnavg))+1,1,1)
       if (max_diff < 1E-006) converged = .True.
    end if

    deallocate(phiavg)
    deallocate(phioldavg)
    deallocate(avg_omega,avg_omega_local,sum_omega,sum_omega_local,diff_omega)
  end subroutine adjoint_convergence
  
  subroutine deallocate_convergence
    use adjoint_field_arrays, only: omega_g, omega

    implicit none

    deallocate(omega_g)
    deallocate(omega)
    deallocate(track)
    deallocate(track_short)
    
  end subroutine deallocate_convergence
    
end module convergence_adjoint

