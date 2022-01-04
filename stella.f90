program stella

  use mp, only: proc0
  use redistribute, only: scatter
  use job_manage, only: time_message, checkstop, job_fork
  use run_parameters, only: nstep, fphi, fapar
  use stella_time, only: update_time, code_time, code_dt
  use dist_redistribute, only: kxkyz2vmu
  use time_advance, only: advance_stella
  use stella_diagnostics, only: diagnose_stella, nsave
  use stella_save, only: stella_save_for_restart
  use dist_fn_arrays, only: gnew, gvmu
  use file_utils, only: error_unit, flush_output_file

  use fields_arrays, only: phi, apar
  
  use millerlocal, only: del
  use convergence_adjoint, only: init_convergence, deallocate_convergence
  use convergence_adjoint, only: omega_convergence1, omega_convergence2
  use convergence_adjoint, only: adjoint_convergence
  
  use adjoint_p_derivatives, only: allocate_pert, allocate_unpert
  use adjoint_p_derivatives, only: deallocate_pert, deallocate_unpert
  use adjoint_p_derivatives, only: perturb_p!, lagrangian_integrals
  use adjoint_p_derivatives, only: get_denominator, calculate_chi_lam

!  use adjoint_distfn_arrays, only: unpert_l, pert_l
  use adjoint_distfn_arrays, only: gsave
  use adjoint_distfn_arrays, only: lam_save
  use adjoint_field_arrays, only: phi_save
  use adjoint_field_arrays, only: chi_save
!  use adjoint_field_arrays, only: unpert_q, pert_q
  use adjoint_field_arrays, only: pert_term, unpert_term
  use adjoint_field_arrays, only: omega_g, omega
  use adjoint_field_arrays, only: derivative, denominator
  
  use constants, only: zi
  use zgrid, only: nzgrid, ntubes
  use stella_diagnostics, only: omega_vs_time, navg
  use fields, only: advance_fields, fields_updated
  use dist_fn_arrays, only:  g_0
  use stella_layouts, only: vmu_lo
  use stella_layouts, only: iv_idx, imu_idx, is_idx
  use species, only: spec
  use vpamu_grids, only:  maxwell_vpa, maxwell_mu, maxwell_fac
  use gyro_averages, only: gyro_average
  use kt_grids, only: nakx, naky
  use	adjoint_distfn_arrays, only: g_omega
  implicit none

  ! Add the version number and date of last change when uploading to github
  character(len=4), parameter :: VERNUM = '0.3'
  character(len=10), parameter :: VERDATE = '2021.03.26'

  logical :: debug = .false.
  logical :: stop_stella = .false.
  logical :: mpi_initialized = .false.

  integer :: istep0, istep, ierr
  integer :: istatus
  real, dimension (2) :: time_init = 0.
  real, dimension (2) :: time_diagnostics = 0.
  real, dimension (2) :: time_total = 0.

  integer :: change_p, np
  integer :: istep_initial

  logical :: stop_convergence
  logical :: do_average
  logical :: converged, convergedk
  logical :: adjoint = .True.

  integer :: ivmu,imu,iv,is, ia, iky,ikx, iz
  complex, dimension(:,:,:,:,:), allocatable :: gyro
  np = 14
  
  ! Initiate stella
  call init_stella(istep0, VERNUM, VERDATE)

  allocate(gyro(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  gyro = 0.
  call init_convergence
  call allocate_adjoint_variables
  call allocate_unpert (np)
  
  ! Add a header to the output file
  if (proc0) then
    write (*,'(A)') "############################################################"
    write (*,'(A)') "                OVERVIEW OF THE SIMULATION"
    write (*,'(A)') "############################################################"
    write (*,'(A)') " "
    write (*,'(A)') "    istep       time           dt         |phi|^2"
    write (*,'(A)') "------------------------------------------------------------"
  end if

  do_average = .True.
  converged = .False.
  stop_convergence = .False.
  
  ! Diagnose stella
  if (debug) write(*,*) 'stella::diagnose_stella'
  if (istep0.eq.0) call diagnose_stella (istep0)

  !! Run Stella as usual
  
  ! Advance stella until istep=nstep
  if (debug) write(*,*) 'stella::advance_stella'
  do istep = (istep0+1), nstep
     if (debug) write(*,*) 'istep = ', istep
     call advance_stella(istep)
     call update_time
     if (nsave > 0 .and. mod(istep,nsave)==0) then
        call scatter (kxkyz2vmu, gnew, gvmu)
        call stella_save_for_restart (gvmu, istep, code_time, code_dt, istatus)
     end if
     call time_message(.false.,time_diagnostics,' diagnostics')
     call diagnose_stella (istep)
     call time_message(.false.,time_diagnostics,' diagnostics')
     !! Run convergence test on frequency to check if phi is a single normal mode
     !! This is split into two parts; an initial check- omega_convergence1
     !! and a second convergence test - omega_convergence2
     !! fist test just checks two frequencies as adjacent time steps to see how similar
     !! second test compares two window averages - one uses the last navg time steps
     !! and one uses navg/2 time steps (rounded to be an integer)
     if (do_average) then
        call omega_convergence1 (istep, converged)
        if (converged) then
           do_average = .False.
           istep_initial = istep
        end if
     else
        call omega_convergence2 (istep, istep_initial, stop_convergence)
     end if
     if (stop_convergence) exit
     if (mod(istep,10)==0) call checkstop (stop_stella)
     if (stop_stella) exit
     ierr = error_unit()
     call flush_output_file (ierr)
  end do

  !! Store frequency: omega_g = gamma + i*omega 
  omega_g = -zi*omega_vs_time(mod(istep,navg)+1,:,:)

  ! !! store final values of g and phi as gsave and phi_save
  ! gsave = gnew
  ! fields_updated = .false.
  ! call advance_fields (gsave, phi_save, apar, dist='gbar')
  
  ! call perturb_p   
  
  deallocate(g_0)
  deallocate(gyro)
  
  call deallocate_adjoint_variables 
  call finish_stella(last_call = .True.)
  
!   !!!!  Everything Adjoint !!!!                                                                                                                                                    
!   gsave = gnew
!   fields_updated = .false.
!   call advance_fields (gsave, phi, apar, dist='gbar')
!   phi_save = phi
!   omega_g = -zi*sum(omega_vs_time,dim=1)/real(navg)
! !  call perturb_p
! !  omega_g = -zi*omega_vs_time(mod(nstep,navg)+1,:,:)

!   call finish_stella

  ! !!!! Calculate Adjoint Variable !!!! 
  ! call init_stella (istep0, VERNUM, VERDATE)

  ! if (proc0) then
  !    write (*,'(A)') "*************************** "
  !    write (*,'(A)') "**** starting adjoint ***** "
  !    write (*,'(A)') "*************************** "
  ! end if

  ! if (istep0.eq.0) call diagnose_stella (istep0)

  ! do_average = .True.
  ! converged = .False.
  ! stop_convergence = .False.

  ! do istep = (istep0+1), nstep
  !    if (debug) write(*,*) 'istep = ', istep
  !    call advance_stella(istep)
  !    call update_time
  !    if (nsave > 0 .and. mod(istep,nsave)==0) then
  !       call scatter (kxkyz2vmu, gnew, gvmu)
  !       call stella_save_for_restart (gvmu, istep, code_time, code_dt, istatus)
  !    end if
  !    call time_message(.false.,time_diagnostics,' diagnostics')
  !    call diagnose_stella (istep)
  !    call time_message(.false.,time_diagnostics,' diagnostics')
  !    call adjoint_convergence (istep, stop_convergence)
  !    if (stop_convergence) exit
  !    if (mod(istep,10)==0) call checkstop (stop_stella)
  !    if (stop_stella) exit
  !    ierr = error_unit()
  !    call flush_output_file (ierr)
  ! end do

  ! lam_save = gnew
  ! apar = 0.
  ! fields_updated = .false.
  ! call advance_fields (lam_save, chi_save, apar, dist='gbar',adjoint=adjoint)
  ! call calculate_chi_lam
  ! call perturb_p (unpert_term)
  ! call get_denominator
  ! call finish_stella  
  
  ! do change_p = 1, np
  !    write(*,*) 'change_p', change_p
  !    call init_stella (istep0, VERNUM, VERDATE, change_p)     
  !    call perturb_p (pert_term)

  !    derivative(1,1,change_p) = -(pert_term(1,1) - unpert_term(1,1))/(denominator(1,1)*del(change_p))
     
  !    open(12, file="with_nothing.txt", status="unknown",action="write",position="append")
  !    write(12,*) 'change_p', change_p, real(derivative(1,1,change_p))
  !    close(12)
  !    if (change_p .EQ. np) then
  !       call deallocate_pert
  !       call deallocate_unpert
  !       call deallocate_adjoint_variables
  !       call deallocate_convergence
  !       call finish_stella (last_call = .True.)
  !    else
  !       call deallocate_pert
  !       call finish_stella
  !       write(*,*) 'change_p', change_p, 'finish stella'
  !    end if
  ! end do
  
contains

  subroutine allocate_adjoint_variables
    use adjoint_distfn_arrays, only: gsave, g_omega
    use adjoint_distfn_arrays, only: lam_save
    use adjoint_field_arrays, only: phi_save
    use adjoint_field_arrays, only: chi_save
    
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_layouts, only: vmu_lo

    implicit none

    if(.not. allocated(g_omega)) then
       allocate(g_omega(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       g_omega = 0.
    end if
    
    if(.not. allocated(gsave)) then
       allocate(gsave(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       gsave = 0.
    end if
    if(.not. allocated(lam_save)) then
       allocate(lam_save(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       lam_save = 0.
    end if

    if(.not. allocated(phi_save)) then
       allocate(phi_save(naky,nakx,-nzgrid:nzgrid,ntubes))
       phi_save = 0.
    end if

    if(.not. allocated(chi_save)) then
       allocate(chi_save(naky,nakx,-nzgrid:nzgrid,ntubes))
       chi_save = 0.
    end if
  end subroutine allocate_adjoint_variables

  subroutine deallocate_adjoint_variables

    use adjoint_distfn_arrays, only: gsave, lam_save, g_omega
    use adjoint_field_arrays, only: phi_save, chi_save

    implicit none
    
    deallocate(g_omega)
    deallocate(gsave)
    deallocate(lam_save)
    deallocate(phi_save)
    deallocate(chi_save)

  end subroutine deallocate_adjoint_variables
  !==============================================
  !============ INITIATE STELLA =================
  !==============================================
  subroutine init_stella(istep0, VERNUM, VERDATE,change_p)
    use mp, only: init_mp, broadcast, sum_allreduce
    use mp, only: proc0,job, scope, subprocs, crossdomprocs
    use file_utils, only: init_file_utils
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: run_name, init_job_name
    use file_utils, only: flush_output_file, error_unit
    use job_manage, only: checktime, time_message, njobs
    use physics_parameters, only: init_physics_parameters, g_exb, g_exbfac
    use physics_flags, only: init_physics_flags
    use physics_flags, only: nonlinear, include_parallel_nonlinearity
    use physics_flags, only: full_flux_surface, radial_variation
    use physics_flags, only: hammett_flow_shear
    use run_parameters, only: init_run_parameters
    use run_parameters, only: avail_cpu_time, nstep, rng_seed, delt
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use run_parameters, only: delt_option_switch, delt_option_auto
    use run_parameters, only: mat_gen, mat_read
    use species, only: init_species, read_species_knobs
    use species, only: nspec, communicate_species_multibox
    use zgrid, only: init_zgrid
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: init_geometry, communicate_geo_multibox
    use stella_geometry, only: finish_init_geometry
    use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
    use response_matrix, only: init_response_matrix, read_response_matrix
    use init_g, only: ginit, init_init_g, phiinit, scale_to_phiinit
    use fields, only: init_fields, advance_fields, get_radial_correction, fields_updated
    use stella_time, only: init_tstart, init_delt
    use init_g, only: tstart
    use stella_diagnostics, only: init_stella_diagnostics
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gnew, gvmu
    use dist_fn, only: init_gxyz, init_dist_fn
    use dist_redistribute, only: init_redistribute
    use time_advance, only: init_time_advance
    use extended_zgrid, only: init_extended_zgrid
    use kt_grids, only: init_kt_grids, read_kt_grids_parameters
    use kt_grids, only: naky, nakx, ny, nx, nalpha
    use vpamu_grids, only: init_vpamu_grids, read_vpamu_grids_parameters
    use vpamu_grids, only: nvgrid, nmu
    use stella_transforms, only: init_transforms
    use stella_save, only: init_dt
    use multibox, only: read_multibox_parameters, init_multibox, rhoL, rhoR
    use multibox, only: communicate_multibox_parameters, multibox_communicate
    use ran, only: get_rnd_seed_length, init_ranf
    use volume_averages, only: init_volume_averages, volume_average

    implicit none

    integer, intent (out) :: istep0
    character(len=4), intent (in) :: VERNUM 
    character(len=10), intent (in) :: VERDATE 
    logical :: exit, list, restarted, needs_transforms
    character (500), target :: cbuff
    integer, dimension (:), allocatable  :: seed
    integer :: i, n, ierr
    real :: delt_saved, phi2, rescale

    integer, optional, intent(in) :: change_p

    ! initialize mpi message passing
    if (.not.mpi_initialized) call init_mp
    mpi_initialized = .true.

    ! initialize timer
    if (debug) write (*,*) 'stella::init_stella::check_time'
    call checktime(avail_cpu_time,exit)

    if (proc0) then
       ! write message to screen with useful info regarding start of simulation
       if (debug) write (*,*) 'stella::init_stella::write_start_message'
       call write_start_message(VERNUM, VERDATE)
       ! initialize file i/o
       if (debug) write (*,*) 'stella::init_stella::init_file_utils'
       call init_file_utils (list)
       call time_message(.false.,time_total,' Total')
       call time_message(.false.,time_init,' Initialization')
    end if

    call broadcast (list)
    call broadcast (runtype_option_switch)
    if(list) call job_fork

    debug = debug .and. proc0

    if (proc0) cbuff = trim(run_name)
    call broadcast (cbuff)
    if (.not. proc0) call init_job_name(cbuff)


    if (debug) write(6,*) "stella::init_stella::init_physics_flags"
    call init_physics_flags
    if (debug) write(6,*) "stella::init_stella::init_physics_parameters"
    call init_physics_parameters
    if (debug) write(6,*) "stella::init_stella::init_zgrid"
    call init_zgrid
    if (debug) write (6,*) "stella::init_stella::read_species_knobs"
    call read_species_knobs
    if (debug) write (6,*) "stella::init_stella::read_multibox_parameters"
    call read_kt_grids_parameters
    if (debug) write (6,*) "stella::init_stella::read_vpamu_grids_parameters"
    call read_multibox_parameters
    if (debug) write (6,*) "stella::init_stella::read_kt_grids_parameters"
    call read_vpamu_grids_parameters
    if (debug) write (6,*) "stella::init_stella::init_dist_fn_layouts"
    call init_dist_fn_layouts (nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
    if (debug) write(6,*) "stella::init_stella::init_geometry"
    if (present (change_p)) then
       call init_geometry (nalpha, change_p)
    else
       call init_geometry (nalpha)
    end if
    if (debug) write (6,*) 'stella::init_stella::init_species'
    call init_species
    if (debug) write(6,*) "stella::init_stella::init_init_g"
    call init_init_g
    if (debug) write(6,*) "stella::init_stella::init_run_parameters"
    call init_run_parameters

    if (debug) write(6,*) "stella::init_stella::init_ranf"
    n=get_rnd_seed_length()
    allocate(seed(n))
    if(rng_seed .lt. 0) then
      call init_ranf(.true.,seed,job+2)
    else
      seed = rng_seed + 37 * (/ ( i - 1, i = 1, n) /)
      call init_ranf(.false.,seed,job+2)
    endif
    deallocate(seed)

    if (debug) write (6,*) 'stella::init_stella::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (6,*) 'stella::init_stella::init_kt_grids'
    call init_kt_grids
    !if (nonlinear .or. full_flux_surface .or. include_parallel_nonlinearity & 
    !    .or. radial_variation .or. (g_exb*g_exb).gt.epsilon(0.0).or. &
    !    runtype_option_switch.eq.runtype_multibox) then
    needs_transforms = .false.
    if(nonlinear.or.include_parallel_nonlinearity) needs_transforms = .true.
    if(radial_variation.or.full_flux_surface)      needs_transforms = .true.
    if(runtype_option_switch.eq.runtype_multibox)  needs_transforms = .true.
    if(abs(g_exb*g_exbfac).gt.epsilon(0.).and..not.hammett_flow_shear) & 
      needs_transforms = .true.
    if (needs_transforms) then
       if (debug) write (*,*) "stella::init_stella::init_transforms"
       call init_transforms
    end if
    if (debug) write (6,*) 'stella::init_stella::init_multibox'
    call init_multibox
    if (proc0.and.runtype_option_switch.eq.runtype_multibox &
             .and.(job.eq.1).and.radial_variation) then
      if (debug) write (6,*) 'stella::init_stella::init_multibox_geo'
      call communicate_geo_multibox(rhoL,rhoR)
      if (debug) write (6,*) 'stella::init_stella::init_multibox_spec'
      call communicate_species_multibox(rhoL,rhoR)
    endif
    if (runtype_option_switch.eq.runtype_multibox.and.(job.eq.1)) then
      call communicate_multibox_parameters
    endif
    if (debug) write (6,*) 'stella::init_stella::finish_init_geometry'
    call finish_init_geometry
    if (debug) write (6,*) 'stella::init_stella::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write (6,*) 'stella::init_stella::init_extended_zgrid'
    call init_extended_zgrid
    if (debug) write (6,*) 'stella::init_stella::init_volume_averages'
    call init_volume_averages
    if (debug) write(6,*) "stella::init_stella::init_dist_fn"
    call init_dist_fn
    if (debug) write(6,*) "stella::init_stella::init_redistribute"
    call init_redistribute
    if (debug) write (6,*) 'stella::init_stella::init_fields'
    call init_fields
    if (debug) write(6,*) "stella::init_stella::ginit"
    call ginit (restarted,istep0)
    if (debug) write(6,*) "stella::init_stella::init_gxyz"
    call init_gxyz (restarted)

    if(restarted.and.delt_option_switch == delt_option_auto) then
      delt_saved = delt
      if (debug) write(6,*) "stella::init_stella::init_dt"
      call init_dt(delt_saved, istatus)
      if(istatus == 0) delt = delt_saved
    endif
    if (debug) write(6,*) "stella::init_stella::init_delt"
    call init_delt(delt)
    if (debug) write (6,*) 'stella::init_stella::init_time_advance'
    call init_time_advance
    if (stream_implicit .or. driftkinetic_implicit) then
       if (mat_read) then
          if (debug) write (6,*) "stella::init_stella::read_response_matrix"
          call read_response_matrix
       else
          if (debug) write (6,*) "stella::init_stella::init_response_matrix"
          call init_response_matrix
       end if
    end if

    if (debug) write (6,*) 'stella::init_stella::get_fields'
    ! get initial field from initial distribution function
    call advance_fields (gnew, phi, apar, dist='gbar')
    if(radial_variation) call get_radial_correction(gnew,phi,dist='gbar')

    if(runtype_option_switch.eq.runtype_multibox) then
      call multibox_communicate (gnew)
      if(job.eq.1) then
        fields_updated=.false.
        call advance_fields (gnew, phi, apar, dist='gbar')
      endif
    endif

    ! FLAG - the following code should probably go elsewhere
    if(.not.restarted.and.scale_to_phiinit) then
      call volume_average(phi,phi2)
      if(runtype_option_switch.eq.runtype_multibox) then
        call scope(crossdomprocs)
        call sum_allreduce(phi2)
        call scope(subprocs)
        phi2=phi2/njobs
      endif
      rescale=phiinit/sqrt(phi2)
      phi  = rescale*phi
      gnew = rescale*gnew
      gvmu = rescale*gvmu
    endif

    if (debug) write (6,*) 'stella::init_stella::init_stella_diagnostics'
    call init_stella_diagnostics (restarted,tstart)
    if (debug) write (6,*) 'stella::init_stella::init_tstart'
    call init_tstart (tstart)

    ierr = error_unit()
    if (proc0) call flush_output_file (ierr)

    if (proc0) call time_message(.false.,time_init,' Initialization')

  end subroutine init_stella

  !==============================================
  !============ WRITE START MESSAGE =============
  !==============================================
  subroutine write_start_message(VERNUM, VERDATE)
    use mp, only: proc0, nproc

    implicit none

    character(len=23) :: str
    character(len=4), intent (in) :: VERNUM 
    character(len=10), intent (in) :: VERDATE 

    if (proc0) then
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ''//achar(27)//'[32m'
      write (*,*) "            I8            ,dPYb, ,dPYb,            "
      write (*,*) "            I8            IP'`Yb IP'`Yb            "
      write (*,*) "         88888888         I8  8I I8  8I            "
      write (*,*) "            I8            I8  8' I8  8'            "
      write (*,*) "   ,g,      I8    ,ggg,   I8 dP  I8 dP    ,gggg,gg "
      write (*,*) "  ,8'8,     I8   i8' '8i  I8dP   I8dP    dP'  'Y8I "
      write (*,*) " ,8'  Yb   ,I8,  I8, ,8I  I8P    I8P    i8'    ,8I "
      write (*,*) ",8'_   8) ,d88b, `YbadP' ,d8b,_ ,d8b,_ ,d8,   ,d8b,"
      write (*,*) 'P` "YY8P8P8P""Y8888P"Y8888P`"Y888P`"Y88P"Y8888P"`Y8'
      write (*,*) ''//achar(27)//'[0m'
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) '                       Version ', VERNUM
      write (*,*) '                        ', VERDATE
      write (*,*) ' '
      write (*,*) '                   Add author names.,'
      write (*,*) '                  More author names...,'
      write (*,*) ' '
      write (*,*) '                   Add institutions...'
      write (*,*) ' '
      write (*,*) ' '
      write (*,'(A)') "############################################################"
      write (*,'(A)') "                     PARALLEL COMPUTING"
      write (*,'(A)') "############################################################"
      if (nproc==1) then
         write (str,'(I10, A)') nproc, " processor."
         write (*,'(A,A,A)') " Running on ", adjustl(trim(str))
      else
         write (str,'(I10, A)') nproc, " processors."
         write (*,'(A,A,A)') " Running on ", adjustl(trim(str))
      end if
      write(*,*)
    end if

  end subroutine write_start_message

  !==============================================
  !=============== FINISH STELLA ================
  !==============================================
  subroutine finish_stella (last_call)

    use mp, only: finish_mp
    use mp, only: proc0
    use file_utils, only: finish_file_utils
    use job_manage, only: time_message
    use physics_parameters, only: finish_physics_parameters
    use physics_flags, only: finish_physics_flags
    use run_parameters, only: finish_run_parameters, nstep
    use zgrid, only: finish_zgrid
    use species, only: finish_species
    use time_advance, only: time_gke, time_parallel_nl
    use time_advance, only: finish_time_advance
    use parallel_streaming, only: time_parallel_streaming
    use mirror_terms, only: time_mirror
    use dissipation, only: time_collisions
    use init_g, only: finish_init_g
    use dist_fn, only: finish_dist_fn
    use dist_redistribute, only: finish_redistribute
    use fields, only: finish_fields
    use fields, only: time_field_solve
    use stella_diagnostics, only: finish_stella_diagnostics
    use response_matrix, only: finish_response_matrix
    use stella_geometry, only: finish_geometry
    use extended_zgrid, only: finish_extended_zgrid
    use vpamu_grids, only: finish_vpamu_grids
    use kt_grids, only: finish_kt_grids
    use volume_averages, only: finish_volume_averages

    implicit none

    logical, intent (in), optional :: last_call

    if (debug) write (*,*) 'stella::finish_stella::finish_stella_diagnostics'
    call finish_stella_diagnostics(nstep)
    if (debug) write (*,*) 'stella::finish_stella::finish_response_matrix'
    call finish_response_matrix
    if (debug) write (*,*) 'stella::finish_stella::finish_fields'
    call finish_fields
    if (debug) write (*,*) 'stella::finish_stella::finish_time_advance'
    call finish_time_advance
    if (debug) write (*,*) 'stella::finish_stella::finish_volume_averages'
    call finish_volume_averages
    if (debug) write (*,*) 'stella::finish_stella::finish_extended_zgrid'
    call finish_extended_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_dist_fn'
    call finish_dist_fn
    if (debug) write (*,*) 'stella::finish_stella::finish_redistribute'
    call finish_redistribute
    if (debug) write (*,*) 'stella::finish_stella::finish_init_g'
    call finish_init_g
    if (debug) write (*,*) 'stella::finish_stella::finish_vpamu_grids'
    call finish_vpamu_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_kt_grids'
    call finish_kt_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_run_parameters'
    call finish_run_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_species'
    call finish_species
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_flags'
    call finish_physics_flags
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_parameters'
    call finish_physics_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_geometry'
    call finish_geometry
    if (debug) write (*,*) 'stella::finish_stella::finish_zgrid'
    call finish_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_file_utils'
    if (proc0) then
       call finish_file_utils
       call time_message(.false.,time_total,' Total')
       write (*,*)
       write (*,'(A)') "############################################################"
       write (*,'(A)') "                        ELAPSED TIME"
       write (*,'(A)') "############################################################"
       write (*,fmt=101) 'initialization:', time_init(1)/60., 'min'
       write (*,fmt=101) 'diagnostics:', time_diagnostics(1)/60., 'min'
       write (*,fmt=101) 'fields:', time_field_solve(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_field_solve(1,2)/60., 'min'
       write (*,fmt=101) 'mirror:', time_mirror(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_mirror(1,2)/60., 'min'
       write (*,fmt=101) 'stream:', time_parallel_streaming(1)/60., 'min'
       write (*,fmt=101) 'dgdx:', time_gke(1,5)/60., 'min'
       write (*,fmt=101) 'dgdy:', time_gke(1,4)/60., 'min'
       write (*,fmt=101) 'wstar:', time_gke(1,6)/60., 'min'
       write (*,fmt=101) 'collisions:', time_collisions(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_collisions(1,2)/60., 'min'
       write (*,fmt=101) 'ExB nonlin:', time_gke(1,7)/60., 'min'
       write (*,fmt=101) 'parallel nonlin:', time_parallel_nl(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_parallel_nl(1,2)/60., 'min'
       write (*,fmt=101) 'total implicit: ', time_gke(1,9)/60., 'min'
       write (*,fmt=101) 'total explicit: ', time_gke(1,8)/60., 'min'
       write (*,fmt=101) 'total:', time_total(1)/60., 'min'
       write (*,*)
    end if
101 format (a17,0pf8.2,a4)

    if (debug) write (*,*) 'stella::finish_stella::finish_mp'
    ! finish (clean up) mpi message passing
    if (present(last_call)) then
       call finish_mp
       mpi_initialized = .false.
    end if

  end subroutine finish_stella

end program stella
