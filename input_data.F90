!> @file
!! @brief Read atmospheric and surface data from MPAS diag, hist, or forecast files.
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

!> Read atmospheric and/or surface input MPAS grid.
!! Supported formats include INITIalization file, forecast diag file (only 2d diagnostic
!! fields), and full 3d forecast files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

 module input_data

 use esmf
 use netcdf
 use utils_mod
 use program_setup, only          : nfiles_max, data_files_input_grid, variable_list_file, &
                                    esmf_input_mesh, esmf_input_grid

 use model_grid, only             : input_mesh, input_grid,   &
                                    nCells_input,  &
                                    nz_input, nzp1_input, &
                                    nsoil_input, strlen, &
                                    input_names, target_names, &
                                    interp_methods, vert_info, &
                                    input_bundle, &
                                    target_units, &
                                    target_longname, &
                                    nfields, nVertLevelsPerVariable, &
                                    elemIDs, nCellsPerPET, nodeIDs, &
                                    i_input, j_input, num_tiles_input_grid
 implicit none

 private
 public :: read_input_data
                                         
 contains

!> Read input data
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames CIWRO/NOAA/NSSL   
 subroutine read_input_data(localpet)

 implicit none

 include 'mpif.h'

 integer, intent(in)             :: localpet

 character(len=500)              :: the_file
 character(len=50)               :: vname, att_text
 
 integer                         :: error, ncid, rc
 integer                         :: id_dim, nfiles, tile
 integer                         :: id_var, i, j, n, nodes, nz

 type(esmf_field),allocatable    :: fields(:)
 
 real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)

 real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
 
 call init_input_fields(localpet)

 nfiles = 0
 do n = 1,nfiles_max
    if(    trim(adjustl(data_files_input_grid(n))) == '' &
      .or. trim(adjustl(data_files_input_grid(n))) == 'NULL' ) exit
    nfiles = nfiles + 1
 enddo
 if ( nfiles == 0 ) call error_handler("NO INPUT FILES", -2)

!the_file = trim(data_files_input_grid(1))
!error=nf90_open(trim(the_file),nf90_nowrite,ncid)
!call netcdf_err(error, 'opening: '//trim(the_file) )

!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for bilinear/patch interpolation
!---------------------------------------------------------------------------

 if (localpet==0) print*, "Begin reading variables"
 
 allocate(fields(nfields))
 call ESMF_FieldBundleGet(input_bundle, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN FieldBundleGet", rc)
 
 allocate(target_units(nfields))
 allocate(target_longname(nfields))

 do i = 1,nfields
    
   vname = trim(adjustl(input_names(i)))
  !call ESMF_FieldGet(fields(i), name=vname, rc=rc)
  !if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
  !       call error_handler("IN FieldGet", rc)

  !if (localpet==0) print*,"- READ ", trim(vname)
  !error=nf90_inq_varid(ncid, trim(vname), id_var)
  !call netcdf_err(error, 'reading field id' )

   ! if the grid is just on 1 tile, then we can look in many files for the desired variable
   !  if there are > 1 tiles, we are going to assume the variable is in the file. Currently,
   !  num_tiles_input_grid > 1 only for FV3-based input
   if ( num_tiles_input_grid == 1 ) then
      do n = 1,nfiles
         the_file = trim(adjustl(data_files_input_grid(n)))
         error=nf90_open(trim(the_file),nf90_nowrite,ncid)
         call netcdf_err(error, 'opening: '//trim(the_file) )
         error=nf90_inq_varid(ncid, trim(vname), id_var)
         if (error==0) then ! if true, the variable is in the file
            if (localpet==0) write(*,*)"- READ "//trim(vname)//" from "//trim(adjustl(data_files_input_grid(n)))
            exit ! if the same variable is in multiple files, "exit" means we get the variable from the earliest file in the list
         endif
         error=nf90_close(ncid) ! only occurs if the variable was not in the nth file
         if (n == nfiles) call netcdf_err(error, 'reading field id' ) ! if true, the variable is not in any file
      enddo
   else
      if (localpet < num_tiles_input_grid) then
         tile = localpet + 1
         the_file = trim(adjustl(data_files_input_grid(tile)))
         error=nf90_open(trim(the_file),nf90_nowrite,ncid)
         call netcdf_err(error, 'opening: '//trim(the_file) )
         error=nf90_inq_varid(ncid, trim(vname), id_var)
         if (error/=0) call netcdf_err(error, 'reading field id' ) ! if true, the variable is not in any file
         write(*,*)"- READ "//trim(vname)//" from "//trim(adjustl(data_files_input_grid(tile)))
      else
         ! only open the file to get id_var, which is queried by all processors
         !  to get the units and long_names attributes.  also makes nf90_close statement at
         !  end of subroutine work for all processors
         the_file = trim(adjustl(data_files_input_grid(1)))
         error=nf90_open(trim(the_file),nf90_nowrite,ncid)
         call netcdf_err(error, 'opening: '//trim(the_file) )
         error=nf90_inq_varid(ncid, trim(vname), id_var)
      endif
   endif
   error=nf90_get_att(ncid,id_var,'units',target_units(i))
   call netcdf_err(error, 'reading field units' )
   error=nf90_get_att(ncid,id_var,'long_name',target_longname(i))
   call netcdf_err(error, 'reading field long_name' )

   nz = nVertLevelsPerVariable(i) ! only used if 3-d variable

   if ( esmf_input_mesh ) then
      if ( vert_info(i).eq.'2d') then
         allocate(dummy2(nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         error=nf90_get_var(ncid, id_var, dummy2)
         call netcdf_err(error, 'reading field' )
         do j = 1, nCellsPerPET
            varptr(j) = dummy2(elemIDs(j),1)
         enddo
         nullify(varptr)
         deallocate(dummy2)
      else
         allocate(dummy3(nz,nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         error=nf90_get_var(ncid, id_var, dummy3)
         call netcdf_err(error, 'reading field' )
         do j = 1, nCellsPerPET
            varptr2(j,:) = dummy3(:,elemIDs(j),1)
         enddo
         nullify(varptr2)
         deallocate(dummy3)
      endif
   else if ( esmf_input_grid ) then
      if(localpet==0) print*,"- CALL FieldScatter FOR ",trim(adjustl(vname))
      if ( vert_info(i).eq.'2d') then
         if (localpet < num_tiles_input_grid) then
            allocate(dummy2(i_input, j_input))
            error=nf90_get_var(ncid, id_var, dummy2)
            call netcdf_err(error, 'reading field' )
         else
            allocate(dummy2(0,0))
         endif
         do tile = 1, num_tiles_input_grid
            call ESMF_FieldScatter(fields(i), dummy2, rootpet=tile-1, tile=tile, rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
               call error_handler("IN FieldScatter",rc)
         enddo
         deallocate(dummy2)
      else
         if (localpet < num_tiles_input_grid) then
            allocate(dummy3(i_input, j_input, nz))
            error=nf90_get_var(ncid, id_var, dummy3)
            call netcdf_err(error, 'reading field' )
         else
            allocate(dummy3(0,0,0))
         endif
         do tile = 1, num_tiles_input_grid
            call ESMF_FieldScatter(fields(i), dummy3, rootpet=tile-1, tile=tile, rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
               call error_handler("IN FieldScatter",rc)
         enddo
         deallocate(dummy3)
      endif
   endif
        
   if (localpet==0) print*,"- SET ON MESH ", trim(vname)

   error = nf90_close(ncid) ! file is still open, so close it

 enddo ! end loop over number of variables

 deallocate(fields)
!error = nf90_close(ncid)

 end subroutine read_input_data
 
 subroutine init_input_fields(localpet)
 
    implicit none
 
    integer, intent(in)           :: localpet   
    integer                       :: i, j, k, n, rc, lowBound, upperBound
    type(esmf_field), allocatable :: fields(:)
    
    nfields = 0 ! defined in model_grid, and set in read_varlist
    
    ! input_names, target_names,interp_methods,vert_info also defined in model_grid, and set in read_varlist
    variable_list_file = trim(adjustl(variable_list_file))
    call read_varlist(localpet,variable_list_file,nfields,input_names, target_names,interp_methods,vert_info)
    allocate(fields(nfields))
    allocate(nVertLevelsPerVariable(nfields)) ! defined in model_grid

    do i = 1, nfields
       if ( vert_info(i).eq.'2d') then
          lowBound = 0
          upperBound = 0
       else if ( vert_info(i).eq.'3d') then
          lowBound = 1
          upperBound = nz_input
       else if ( vert_info(i).eq.'3d_soil') then
          lowBound = 1
          upperBound = nsoil_input
       else if ( vert_info(i).eq.'3d_edge') then
          lowBound = 1
          upperBound = nzp1_input
       else
          call error_handler("Invalid specification of vert_info in varlist file", -123)
       endif
       nVertLevelsPerVariable(i) = upperBound
    
       if (localpet==0) print*, "- INIT FIELD ", input_names(i)

       if ( esmf_input_mesh ) then
          if ( vert_info(i).eq.'2d' ) then
             fields(i) = ESMF_FieldCreate(input_mesh, & 
                             typekind=ESMF_TYPEKIND_R8, &
                             meshloc=ESMF_MESHLOC_ELEMENT, &
                             name=input_names(i), rc=rc)
          else
             fields(i) = ESMF_FieldCreate(input_mesh, & 
                             typekind=ESMF_TYPEKIND_R8, &
                             meshloc=ESMF_MESHLOC_ELEMENT, &
                             name=input_names(i), &
                             ungriddedLBound=(/lowBound/), &
                             ungriddedUBound=(/upperBound/), rc=rc)
          endif
       else if ( esmf_input_grid ) then
          if ( vert_info(i).eq.'2d' ) then
             fields(i) = ESMF_FieldCreate(input_grid, &
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             name=input_names(i), rc=rc)
          else
             fields(i) = ESMF_FieldCreate(input_grid, &
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             name=input_names(i), &
                             ungriddedLBound=(/lowBound/), &
                             ungriddedUBound=(/upperBound/), rc=rc)
          endif
       endif
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldCreate", rc)
    enddo
    
    ! create bundle
    input_bundle = ESMF_FieldBundleCreate(fieldList=fields, & 
                                     name="input bundle", rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleCreate", rc)

    deallocate(fields)
    
 end subroutine init_input_fields

subroutine read_varlist(localpet, file,nfields,field_names_input,field_names_target,interp_methods,vert_info) ! CSS added interp_methods, vert_info

    implicit none
   
   integer, INTENT(IN)                      :: localpet 
   character(len=*), INTENT(IN)            :: file
   integer, INTENT(INOUT)                   :: nfields
   character(len=50), INTENT(INOUT), ALLOCATABLE  :: field_names_input(:), field_names_target(:)
   character(len=50), INTENT(INOUT), ALLOCATABLE  :: interp_methods(:), vert_info(:) ! CSS
   
   integer :: k, istat
   character(len=200) :: line
    
   open(14, file=trim(file), form='formatted', iostat=istat)
   if (istat /= 0) then
     call error_handler("OPENING VARLIST FILE", istat)
   endif

   nfields = 0

   !Loop over lines of file to count the number of variables
   do
     read(14, '(A)', iostat=istat) line 
     if (istat/=0) exit
     if ( trim(line) .eq. '' ) cycle
     if ( trim(line(1:1)) .eq. '!' .or. trim(line(1:1)) .eq. '#' ) cycle
     nfields = nfields+1
   enddo
   
   if (localpet==0) print*, "READING ", nfields, " FIELDS ACCORDING TO ", trim(file)
   if ( nfields == 0) call error_handler("VARLIST FILE IS EMPTY.", -1)

   allocate(field_names_input(nfields))
   allocate(field_names_target(nfields))
   allocate(interp_methods(nfields)) ! CSS added
   allocate(vert_info(nfields)) ! CSS added

   rewind(14)
    do k = 1,nfields
      read(14, *, iostat=istat) field_names_input(k), field_names_target(k), interp_methods(k), vert_info(k) 
      if (istat /= 0) call error_handler("READING VARLIST FILE", istat)
      field_names_input(k)  = trim(adjustl(field_names_input(k)))
      field_names_target(k) = trim(adjustl(field_names_target(k)))
      interp_methods(k)     = trim(adjustl(interp_methods(k)))
      vert_info(k)          = trim(adjustl(vert_info(k)))
      if (localpet==0) write(*,*) field_names_input(k), field_names_target(k), interp_methods(k), vert_info(k) 
    enddo
    
   close(14)

end subroutine read_varlist

 end module input_data
