!> @file
!! @brief Interpolate atmospheric and surface data from MPAS diag, init, or forecast files
!! to a target grid.
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

!> Public variables are defined below: "target" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

 module interp

 use esmf
 use netcdf
 use utils_mod

 use model_grid, only   : target_grid, &
                          input_bundle, &
                          target_bundle, & 
                          nVertLevelsPerVariable, &
                          input_names, &
                          target_names, &
                          nfields, &
                          interp_methods

 implicit none

 private
 
 public :: interp_data
 
 real(esmf_kind_r8), parameter    :: spval = 9.9E10
 integer(esmf_kind_i4), pointer   :: unmappedPtr(:)
 integer ::  isrctermprocessing = 1
 
 contains
 
 subroutine interp_data(localpet,npets)
 
    implicit none
     
    integer, intent(in)              :: localpet, npets
    integer                          :: rc, ij, i, j, n
    type(ESMF_RegridMethod_Flag)     :: method
    type(ESMF_RouteHandle)           :: rh_cons1, rh_cons2, rh_nstd, rh_patch, rh_bilin, my_rh
    type(ESMF_Field), allocatable    :: fields_input_grid(:), fields_target_grid(:)
    character(len=50)                :: interp_method
    
    call init_target_fields(localpet)

    allocate(fields_input_grid(nfields))
    allocate(fields_target_grid(nfields))
    call ESMF_FieldBundleGet(input_bundle, fieldList=fields_input_grid, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc)
    call ESMF_FieldBundleGet(target_bundle, fieldList=fields_target_grid, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc)
    if ( any(interp_methods.eq."patch")) then
       if(localpet==0)write(*,*)'making routehandle for patch regridding'
       call make_rh(fields_input_grid(1), fields_target_grid(1),ESMF_REGRIDMETHOD_PATCH, rh_patch)
    endif
    if ( any(interp_methods.eq."nearest")) then
       if(localpet==0)write(*,*)'making routehandle for nearest neighbor regridding'
       call make_rh(fields_input_grid(1), fields_target_grid(1),ESMF_REGRIDMETHOD_NEAREST_STOD, rh_nstd)
    endif
    if ( any(interp_methods.eq."conserve1")) then
       if(localpet==0)write(*,*)'making routehandle for 1st-order conservative regridding'
       call make_rh(fields_input_grid(1), fields_target_grid(1),ESMF_REGRIDMETHOD_CONSERVE, rh_cons1)
    endif
    if ( any(interp_methods.eq."conserve2")) then
       if(localpet==0)write(*,*)'making routehandle for 2nd-order conservative regridding'
       call make_rh(fields_input_grid(1), fields_target_grid(1),ESMF_REGRIDMETHOD_CONSERVE_2ND, rh_cons2)
       !fileName = 'rh_cons2_'//npets//'PETs.dat'
       !call ESMF_RouteHandleWrite(rh_cons2, fileName, rc=rc)
    endif
    if ( any(interp_methods.eq."bilinear")) then
       if(localpet==0)write(*,*)'making routehandle for bilinear regridding'
       call make_rh(fields_input_grid(1), fields_target_grid(1),ESMF_REGRIDMETHOD_BILINEAR, rh_bilin)
    endif

    do i = 1,nfields
    
       interp_method = trim(adjustl(interp_methods(i)))
       if ( interp_method.eq."patch") then
          my_rh = ESMF_RouteHandleCreate(rh_patch, rc=rc)
       else if ( interp_method.eq."nearest") then
          my_rh = ESMF_RouteHandleCreate(rh_nstd, rc=rc)
       else if ( interp_method.eq."conserve1") then
          my_rh = ESMF_RouteHandleCreate(rh_cons1, rc=rc)
       else if ( interp_method.eq."conserve2") then
          my_rh = ESMF_RouteHandleCreate(rh_cons2, rc=rc)
       else if ( interp_method.eq."bilinear") then
          my_rh = ESMF_RouteHandleCreate(rh_bilin, rc=rc)
       else
          call error_handler("Invalid interpolation method = "//interp_method, -223)
       endif
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN RouteHandleCreate", rc)
       if (localpet==0) write(*,*)"- INTERPOLATING "//trim(adjustl(input_names(i)))//" with method "//interp_method

       ! do the interpolation
       call ESMF_FieldRegrid(fields_input_grid(i), fields_target_grid(i), my_rh, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldRegrid", rc)

       ! get rid of the route handle
       call ESMF_FieldRegridRelease(routehandle=my_rh, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    enddo

    ! update the fields in target_bundle
    call ESMF_FieldBundleAddReplace(target_bundle, fields_target_grid, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN ESMF_FieldBundleAddReplace", rc)

    ! clean-up
    deallocate(fields_input_grid,fields_target_grid)

    if ( any(interp_methods.eq."patch")) then
       call ESMF_FieldRegridRelease(routehandle=rh_patch, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    endif
    if ( any(interp_methods.eq."nearest")) then
       call ESMF_FieldRegridRelease(routehandle=rh_nstd, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    endif
    if ( any(interp_methods.eq."conserve1")) then
       call ESMF_FieldRegridRelease(routehandle=rh_cons1, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    endif
    if ( any(interp_methods.eq."conserve2")) then
       call ESMF_FieldRegridRelease(routehandle=rh_cons2, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    endif
    if ( any(interp_methods.eq."bilinear")) then
       call ESMF_FieldRegridRelease(routehandle=rh_bilin, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridRelease", rc)
    endif

  end subroutine interp_data
                
  subroutine init_target_fields(localpet)
 
    implicit none
    
    integer, intent(in)          :: localpet
    integer                      :: i, rc
    integer                      :: lowBound, upperBound
    type(esmf_field),allocatable :: fields(:)
    
    allocate(fields(nfields))

    do i = 1, nfields

       if (localpet==0) print*, "- INIT FIELD ", target_names(i)

       if ( nVertLevelsPerVariable(i) .eq. 0) then
          fields(i) = ESMF_FieldCreate(target_grid, & 
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             name=target_names(i), rc=rc)
       else
          lowBound = 1
          upperBound = nVertLevelsPerVariable(i)
          fields(i) = ESMF_FieldCreate(target_grid, & 
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             name=target_names(i), &
                             ungriddedLBound=(/lowBound/), &
                             ungriddedUBound=(/upperBound/), rc=rc)
       endif
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldCreate", rc)
    enddo
    
    ! create bundle
    target_bundle = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target bundle", rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleCreate", rc)

    deallocate(fields)
 
 end subroutine init_target_fields

 subroutine make_rh(input_field,target_field,method,rh)

    implicit none

    type(ESMF_Field), intent(in)              :: input_field
    type(ESMF_Field), intent(inout)           :: target_field
    type(ESMF_RegridMethod_Flag), intent(in)  :: method
    type(ESMF_RouteHandle), intent(inout)        :: rh

    integer :: rc

    call ESMF_FieldRegridStore(input_field, target_field, &
                                regridmethod=method, &
                                routehandle=rh, &
                                srcTermProcessing=isrctermprocessing, &
                                unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN FieldRegridStore", rc)
    return 
 end subroutine make_rh

end module interp
 
 
