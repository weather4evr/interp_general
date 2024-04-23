!> @file
!! @brief Set up program execution.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> This module contains code to read the setup namelist file.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module program_setup

 use esmf
 use ESMF_LogPublicMod
 use utils_mod, only : to_upper, error_handler

 implicit none

 private

 ! Public variables
 type(ESMF_LogKind_Flag), public :: LogType 
 integer, parameter, public :: nfiles_max = 200
 logical, public :: esmf_input_mesh = .false.
 logical, public :: esmf_input_grid = .false.
 
 ! Shared namelist variables
 character(len=500), public      :: variable_list_file   = "NULL" !< Full path of file containing list of variables to interpolate
 character(len=500), public      :: data_files_input_grid(nfiles_max) = "NULL" !< Full path of input history MPAS data
 character(len=500), public      :: file_target_grid = "NULL"     !<Full path of file containing target 
                                                                  !<grid information for target_grid_type='file'
 character(len=500), public      :: output_file = "NULL"          !< Full path of output file
 character(len=50), public       :: start_time !< simulation start time input data
 character(len=20), public       :: valid_time !< valid time of the forecast in WRF format
 character(len=200), public      :: input_grid_type = "NULL" !< type of input grid (mpas,fv3)
 logical, public                 :: target_grid_is_regional = .true.  !< .true. if target_grid is regional

 ! MPAS namelist variables
 character(len=500), public      :: file_input_grid_mpas = "NULL" !< Full path of MPAS file containing grid information

 ! FV3 namelist variables
 character(len=500), public      :: orog_files_input_grid(6) = "NULL" !< CSS. Input grid orography files.
 character(len=500), public      :: orog_dir_input_grid = "NULL"     !< CSS. Directory containing the input grid orography files.
 character(len=500), public      :: mosaic_file_input_grid = "NULL" !< Full path of MPAS file containing grid information

 ! Namelist Variables for latlon input grid
 character(len=500), public      :: file_input_grid_latlon = "NULL" !< Full path of file containing grid information
 character(len=500), public      :: lat_var = "NULL"      !< Name of latitude variable in $file_input_grid_latlon
 character(len=500), public      :: lon_var = "NULL"      !< Name of longitude variable in $file_input_grid_latlon
 logical, public                 :: is_regional = .true.  !< .true. if input grid is regional
 real, public                    :: dx_in_degrees = -1.0  !< Horizontal grid spacing in degrees of lat-lon grid

 ! Public subroutins
 public :: read_setup_namelist

 contains

!> Reads program configuration namelist.
!!
!! @param filename the name of the configuration file (defaults to
!! ./fort.41).
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine read_setup_namelist(unum, filename)
 implicit none

 character(len=*), intent(in), optional :: filename
 integer, intent(in), optional          :: unum
 character(:), allocatable              :: filename_to_use
 character(len=500)                     :: map_proj
 integer                                :: unit_to_use
 logical                                :: esmf_log

 integer                                :: is, ie, ierr
 
 namelist /share/ data_files_input_grid, variable_list_file, &
            file_target_grid, output_file, esmf_log, start_time, valid_time, input_grid_type, &
            target_grid_is_regional
 namelist /mpas/ file_input_grid_mpas
 namelist /fv3/ mosaic_file_input_grid, orog_files_input_grid, orog_dir_input_grid
 namelist /latlon/ file_input_grid_latlon, lat_var, lon_var, is_regional, dx_in_degrees

 !print*,"- READ SETUP NAMELIST"

 if (present(filename)) then
    filename_to_use = filename
 else
    filename_to_use = "./fort.41"
 endif

 if (present(unum)) then
     unit_to_use = unum
 else
     unit_to_use = 41
 endif

 open(unit_to_use, file=filename_to_use, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(unit_to_use, nml=share, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST SHARE.", ierr)
 read(unit_to_use, nml=mpas, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST MPAS.", ierr)
 read(unit_to_use, nml=fv3, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST FV3.", ierr)
 read(unit_to_use, nml=latlon, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST LATLON.", ierr)
 close (unit_to_use)

 input_grid_type = trim(adjustl(input_grid_type))
 input_grid_type = to_upper(input_grid_type)

 ! defaults of esmf_input_mesh and esmf_input_grid are both false.
 if ( input_grid_type.eq."MPAS" ) esmf_input_mesh = .true.
 if ( input_grid_type.eq."FV3" )  esmf_input_grid = .true.
 if ( input_grid_type.eq."LATLON" ) esmf_input_grid = .true.

 ! These should never happen, but test for them anyways
 if ( esmf_input_mesh .and. esmf_input_grid ) then
    call error_handler("BOTH esmf_input_mesh and esmf_input_grid ARE TRUE", -2)
 endif
 if ( (.not.esmf_input_mesh) .and. (.not.esmf_input_grid) ) then
    call error_handler("BOTH esmf_input_mesh and esmf_input_grid ARE FALSE", -2)
 endif
 
 if (esmf_log) then
   LogType = ESMF_LOGKIND_MULTI_ON_ERROR
 else
   LogType = ESMF_LOGKIND_NONE
 endif 
 
 end subroutine read_setup_namelist

 end module program_setup
