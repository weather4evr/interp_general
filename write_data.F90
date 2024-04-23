 module write_data

 use utils_mod
 use program_setup, only   : output_file

 private

 public :: write_to_file

 contains

!> Write data on the target grid
!!
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL

 subroutine write_to_file(localpet)

 use esmf
 use netcdf
 use mpi

 use program_setup, only           : start_time, valid_time

 use model_grid, only              : target_grid, &
                                     i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     nz_input, nzp1_input, &
                                     nsoil_input, &
                                     dx, &
                                     strlen, &
                                     cen_lat, cen_lon,  &
                                     truelat1, truelat2, stand_lon, &
                                     pole_lat, pole_lon, &
                                     proj_code, map_proj_char, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     nfields, &
                                     target_bundle, &
                                     target_units, &
                                     target_longname, &
                                     vert_info

 implicit none

 integer, intent(in)              :: localpet

 character(len=128)               :: outfile
 character(len=50)                :: varname
 character(len=20)                 :: tempstr(1,19)

 integer, parameter               :: Datestrlen=19
 integer                          :: error, ncid, n, rc, i, j, k, m
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_time, dim_lon, dim_lat, dim_z, dim_zp1, dim_soil, my_dim_z
 integer                          :: dim_lonp, dim_latp, dim_str, dim_lon_stag, dim_lat_stag
 integer                          :: id_lat, id_lon, id_z, id_times, id_xtime, id_itime, id_var
 integer                          :: n2d, n3d, ndims
 integer                          :: sy,sm,sd,sh,smi,ss,  vy,vm,vd,vh,vmi,vs
 integer                          :: maxinds(2), mininds(2)

 real(esmf_kind_r8), allocatable  :: dum2d(:,:), dum2dt(:,:,:), &
                                     dum3d(:,:,:), dum3dt(:,:,:,:), &
                                     dum3dp1(:,:,:), dum3dp1t(:,:,:,:), &
                                     dumsoil(:,:,:), dumsoilt(:,:,:,:)

 type(esmf_field), allocatable    :: fields(:)

 if (localpet ==0) then
   allocate(dum2d(i_target,j_target))
   allocate(dum2dt(i_target,j_target,1))
   allocate(dum3d(i_target,j_target,nz_input))
   allocate(dum3dt(i_target,j_target,nz_input,1))
   allocate(dum3dp1(i_target,j_target,nzp1_input))
   allocate(dum3dp1t(i_target,j_target,nzp1_input,1))
   allocate(dumsoil(i_target,j_target,nsoil_input))
   allocate(dumsoilt(i_target,j_target,nsoil_input,1))
 else
   allocate(dum2d(0,0))
   allocate(dum2dt(0,0,0))
   allocate(dum3d(0,0,0))
   allocate(dum3dt(0,0,0,0))
   allocate(dum3dp1(0,0,0))
   allocate(dum3dp1t(0,0,0,0))
   allocate(dumsoil(0,0,0))
   allocate(dumsoilt(0,0,0,0))
 endif

if (localpet == 0) then

!--- open the file
   error = nf90_create(output_file, NF90_NETCDF4, ncid)
   call netcdf_err(error, 'CREATING FILE '//trim(output_file) )

!--- define dimensions
   error = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED , dim_time)
   call netcdf_err(error, 'DEFINING Time DIMENSION' )
   error = nf90_def_dim(ncid, 'west_east', i_target, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'west_east_stag', i_target+1, dim_lon_stag)
   call netcdf_err(error, 'DEFINING STAGGERED LON DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north', j_target, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north_stag', j_target+1, dim_lat_stag)
   call netcdf_err(error, 'DEFINING STAGGERED LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top', nz_input, dim_z)
   call netcdf_err(error, 'DEFINING VERTICAL DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top_stag', nzp1_input, dim_zp1)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   error = nf90_def_dim(ncid, 'soil_layers_stag', nsoil_input, dim_soil)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   error = nf90_def_dim(ncid, 'StrLen', Datestrlen, dim_str)
   call netcdf_err(error, 'DEFINING STRLEN DIMENSION' )

 !--- define global attributes
   error = nf90_put_att(ncid, NF90_GLOBAL, 'WEST-EAST_GRID_DIMENSION', i_target+1)
   call netcdf_err(error, 'DEFINING WEST-EAST GRID DIMENSION GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'SOUTH-NORTH_GRID_DIMENSION', j_target+1)
   call netcdf_err(error, 'DEFINING NORTH-SOUTH GRID DIMENSION GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'BOTTOM-TOP_GRID_DIMENSION', nz_input+1)
   call netcdf_err(error, 'DEFINING BOTTOM-TOP GRID DIMENSION GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'SIMULATION_START_DATE', start_time)
   call netcdf_err(error, 'DEFINING SUMLATION START DATE GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'START_DATE', start_time)
   call netcdf_err(error, 'DEFINING START DATE GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'VALID_TIME', valid_time)
   call netcdf_err(error, 'DEFINING START DATE GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'DX', dx)
   call netcdf_err(error, 'DEFINING DX GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'CEN_LAT', cen_lat)
   call netcdf_err(error, 'DEFINING CEN_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'CEN_LON', cen_lon)
   call netcdf_err(error, 'DEFINING CEN_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'TRUELAT1', truelat1)
   call netcdf_err(error, 'DEFINING TRUELAT1 GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'TRUELAT2', truelat2)
   call netcdf_err(error, 'DEFINING TRUELAT2 GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MOAD_CEN_LAT', cen_lat)
   call netcdf_err(error, 'DEFINING MOAD_CEN_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'STAND_LON', stand_lon)
   call netcdf_err(error, 'DEFINING STAND_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
   call netcdf_err(error, 'DEFINING POLELAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POLE_LON', pole_lon)
   call netcdf_err(error, 'DEFINING POLE_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
   call netcdf_err(error, 'DEFINING POLE_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MAP_PROJ', proj_code)
   call netcdf_err(error, 'DEFINING MAP_PROJ GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MAP_PROJ_CHAR', map_proj_char)
   call netcdf_err(error, 'DEFINING MAP_PROJ_CHAR GLOBAL ATTRIBUTE')

!--- define fields
   error = nf90_def_var(ncid, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lon, "description", "LONGITUDE, WEST IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lon, "units", "degree_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )
   error = nf90_put_att(ncid, id_lon, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lon, "coordinates", "XLONG XLAT")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lon, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lon, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "description", "LATITUDE, SOUTH IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_lat, "units", "degree_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   error = nf90_put_att(ncid, id_lat, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lat, "coordinates", "XLONG XLAT")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lat, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lat, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'Times', NF90_CHAR, (/dim_str, dim_time/), id_times)
   call netcdf_err(error, 'DEFINING Times FIELD' )
   error = nf90_put_att(ncid, id_times, "description", "Times")
   call netcdf_err(error, 'DEFINING Times NAME' )
   error = nf90_put_att(ncid, id_times, "units", "m")
   call netcdf_err(error, 'DEFINING Times UNITS' )
   error = nf90_put_att(ncid, id_times, "coordinates", "Time")
   call netcdf_err(error, 'DEFINING Times COORD' )
   error = nf90_put_att(ncid, id_times, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_times, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

endif ! localpet = 0

 allocate(fields(nfields))
 do i = 1,nfields
    call ESMF_FieldBundleGet(target_bundle, fieldList=fields, &
                       itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                       rc=error)  ! can probably move outside the loop
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
           call error_handler("IN FieldBundleGet", error)
    call ESMF_FieldGet(fields(i),name=varname,rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
       call error_handler("IN FieldGet", error)
   !call ESMF_FieldGet(fields(i), dimCount=ndims, rc=rc)
   !if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   !   call error_handler("IN FieldGet", error)

   !if (ndims==2) then
    if ( vert_info(i).eq.'2d') then
       if (localpet==0) then
          print*,"- DEFINE 2d field ON FILE TARGET GRID ", varname
          error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_var)
          call netcdf_err(error, 'DEFINING VAR' )
          error = nf90_put_att(ncid, id_var, "MemoryOrder", "XY " )
          call netcdf_err(error, 'DEFINING MEMORYORDER' )
          error = nf90_put_att(ncid, id_var, "coordinates", "XLONG XLAT XTIME")
          call netcdf_err(error, 'DEFINING COORD' )
          error = nf90_put_att(ncid, id_var, "units", target_units(i))
          call netcdf_err(error, 'DEFINING UNITS' )
          error = nf90_put_att(ncid, id_var, "description", target_longname(i))
          call netcdf_err(error, 'DEFINING LONG_NAME' )
          error = nf90_put_att(ncid, id_var, "stagger", "")
          call netcdf_err(error, 'DEFINING STAGGER' )
          error = nf90_put_att(ncid, id_var, "FieldType", 104)
          call netcdf_err(error, 'DEFINING FieldType' )
       endif

       if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname)
       call ESMF_FieldGather(fields(i), dum2d, rootPet=0, rc=error)
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGather", error)

       if (localpet==0) then
         print*, "- WRITE TO FILE ", trim(varname)
         dum2dt(:,:,1) = dum2d
         error = nf90_put_var( ncid, id_var, dum2dt,count=(/i_target,j_target,1/))
         call netcdf_err(error, 'WRITING RECORD')
       endif

    else ! 3d variable

       if ( vert_info(i).eq.'3d') then
          my_dim_z = dim_z
       else if ( vert_info(i).eq.'3d_soil') then
          my_dim_z = dim_soil
       else if ( vert_info(i).eq.'3d_edge') then
          my_dim_z = dim_zp1
       endif

       if (localpet==0) then
          print*,"- DEFINE 3d fields ON FILE TARGET GRID ", varname
          error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,my_dim_z, dim_time/), id_var)
          call netcdf_err(error, 'DEFINING VAR' )
          error = nf90_put_att(ncid, id_var, "MemoryOrder", "XYZ ")
          call netcdf_err(error, 'DEFINING MEMORYORDER' )
          error = nf90_put_att(ncid, id_var, "coordinates", "XLONG XLAT XTIME")
          call netcdf_err(error, 'DEFINING COORD' )
          error = nf90_put_att(ncid, id_var, "units", target_units(i))
          call netcdf_err(error, 'DEFINING UNITS' )
          error = nf90_put_att(ncid, id_var, "description", target_longname(i))
          call netcdf_err(error, 'DEFINING LONG_NAME' )
          error = nf90_put_att(ncid, id_var, "stagger", "")
          call netcdf_err(error, 'DEFINING STAGGER' )
          error = nf90_put_att(ncid, id_var, "FieldType", 104)
          call netcdf_err(error, 'DEFINING FieldType' )
       endif

       if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname)

       if ( vert_info(i).eq.'3d') then
          call ESMF_FieldGather(fields(i), dum3d, rootPet=0, rc=error)
          if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldGather", error)
          if (localpet==0) then
             dum3dt(:,:,:,1) = dum3d
             error = nf90_put_var( ncid, id_var, dum3dt, &
                                        count=(/i_target,j_target,nz_input,1/))
             call netcdf_err(error, 'WRITING RECORD' )
           endif
       else if ( vert_info(i).eq.'3d_soil') then
          call ESMF_FieldGather(fields(i), dumsoil, rootPet=0, rc=error)
          if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldGather", error)
          if (localpet==0) then
             dumsoilt(:,:,:,1) = dumsoil
             error = nf90_put_var( ncid, id_var, dumsoilt, &
                                        count=(/i_target,j_target,nsoil_input,1/))
             call netcdf_err(error, 'WRITING RECORD' )
           endif
       else if ( vert_info(i).eq.'3d_edge') then
          call ESMF_FieldGather(fields(i), dum3dp1, rootPet=0, rc=error)
          if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldGather", error)
          if (localpet==0) then
             dum3dp1t(:,:,:,1) = dum3dp1
             error = nf90_put_var( ncid, id_var, dum3dp1t, &
                                        count=(/i_target,j_target,nzp1_input,1/))
             call netcdf_err(error, 'WRITING RECORD' )
           endif
       endif 
    endif ! end if 3d vars
 enddo ! loop over fields

!--- write latitude, longitude

! longitude
  if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE"
  call ESMF_FieldGather(longitude_target_grid, dum2d, rootPet=0, rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGather", error)

 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lon, dum2dt, count=(/i_target,j_target,1/))
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

! latitude
  if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
  call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGather", error)
 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lat, dum2dt,count=(/i_target,j_target,1/))
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif

 !times
 if (localpet==0)  print*,"- WRITE TO FILE TARGET GRID Times"
 if (localpet ==0) then
   tempstr(1,:) = valid_time
   error = nf90_put_var( ncid, id_times, tempstr, start = (/1,1/), count=(/Datestrlen,1/))
   call netcdf_err(error, 'WRITING TIMES RECORD' )
 endif

 if (localpet==0) then
  !error = nf90_enddef(ncid, header_buffer_val,4,0,4)
  !call netcdf_err(error, 'DEFINING STUFF' )
   error = nf90_close(ncid)
   call netcdf_err(error, 'CLOSING FILE' )
 endif

 deallocate(fields)
 deallocate(dum3d, dum3dp1, dum3dt, dum3dp1t)
 deallocate(dum2d, dum2dt, dumsoil, dumsoilt)

 end subroutine write_to_file

 end module write_data
