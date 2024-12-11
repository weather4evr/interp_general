!> @file
!! @brief Specify input and target model grids.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> Sets up the ESMF grid objects for the input data grid and target grid.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module model_grid

 use esmf
 use ESMF_LogPublicMod
 use utils_mod
 use program_setup, only       : file_input_grid_mpas, file_target_grid, input_grid_type, &
                                 esmf_input_mesh, esmf_input_grid, file_input_grid_latlon, &
                                 mosaic_file_input_grid, orog_dir_input_grid, orog_files_input_grid, &
                                 is_regional, dx_in_degrees, lat_var, lon_var, &
                                 target_grid_is_regional
 implicit none

 private

 integer, public                        :: nCells_input
                                           !< cells on input grid
 integer, public                        :: nVert_input
                                           !< vertices on input grid
 integer, public                        :: nz_input
                                           !< number of input grid atm layers
 integer, public                        :: nzp1_input
                                           !< number of input grid atm layer interfaces
 integer, public                        :: nsoil_input
                                           !< number of input soil levels
 real, public                           :: dx
                                           !< grid size (m) of target grid
 integer, public                        :: strlen
                                           !< StrLen on input file
 real, public                           :: cen_lat
                                           !< target grid projection center latitude
 real, public                           :: cen_lon
                                           !< target grid projection center longitude
!real, public                           :: moad_cen_lat
!                                          !< target grid moad center latitude
 real, public                           :: truelat1  !< First true latitude (all projections)
 real, public                           :: truelat2  !< Second true latitude (LCC only)
 real, public                           :: stand_lon !< Longitude parallel to y-axis (-180->180E)

 real, public                           :: pole_lat
                                           !< target grid projection pole latitude
 real, public                           :: pole_lon
                                           !< target grid projection pole longitude
 integer, public                        :: proj_code  !< Integer code corresponding to the requested
                                           !< target grid projection type
 character(len=500), public             :: map_proj_char

 integer, public                        :: i_target
                                           !< number of longitude points on target grid
 integer, public                        :: j_target
                                           !< number of latitude points on target grid
 integer, public                        :: ip1_target
                                           !< ip1_target plus 1
 integer, public                        :: jp1_target
                                           !< jp1_target plus 1

 integer, allocatable, public           :: elemIDs(:)
                                           !< IDs of the elements on present PET
 integer, allocatable, public           :: nodeIDs(:)
                                            !< IDs of the nodes on present PET
 integer, public                        :: nCellsPerPET
                                            !< Number of cells on this PET
 integer, public                        :: nNodesPerPET
                                           !< Number of nodes on this PET

 type(esmf_mesh),  public               :: input_mesh ! used for MPAS input
                                           !< input grid esmf mesh object
 type(esmf_grid),  public               :: input_grid ! used for FV3 input
                                           !< input grid esmf grid object
 type(esmf_grid),  public               :: target_grid
                                           !< target grid esmf grid object.

  type(esmf_field),  public              :: zgrid_input_grid
                                           !< esmf field to hold level height on input grid
  type(esmf_field),  public              :: zgrid_target_grid
                                           !< esmf field to hold level height on target grid

 type(esmf_field),  public              :: latitude_target_grid
                                           !< latitude of grid center, target grid
 type(esmf_field),  public              :: longitude_target_grid
                                           !< longitude of grid center, target grid
                                          !< number of fields read from the diag file
 type(esmf_fieldbundle), public        :: input_bundle
                                          !< bundle to hold input diag fields
 type(esmf_fieldbundle), public        :: target_bundle
                                          !< bundle to hold target diag fields
 integer, public                       :: nfields
                                          !< number of fields to regrid
 integer, allocatable, public          :: nVertLevelsPerVariable(:)
                                          !< Array to hold number of vertical levels for each variable
 character(50), allocatable, public    :: input_names(:)
                                          !< Arrays to hold field names in the input file
 character(50), allocatable, public    :: target_names(:)
                                          !< Arrays to hold target field names
 character(50), allocatable, public    :: interp_methods(:)
                                          !< Array to hold interpolation method for each variable
 character(50), allocatable, public    :: vert_info(:)
                                          !< Array to hold information about whether each variable is 2D or 3D
 character(50), allocatable, public    :: target_units(:)
                                          !< Array to hold target field units
 character(200), allocatable, public   :: target_longname(:)
                                          !< Arrays to hold target field longname
 integer, public :: i_input !< i-dimension of input grid for FV3/regular grid
                                           !! (or of each global tile)
 integer, public :: j_input !< j-dimension of input grid for FV3/regular grid
                                           !! (or of each global tile)
 integer, public :: num_tiles_input_grid !< number of tiles on the input grid for FV3

 ! public subroutines
 public :: define_target_grid
 public :: define_input_grid
 public :: cleanup_input_target_grid_data

 ! variables visible to this module, but not public
 integer :: ip1_input
 integer :: jp1_input

 contains

!> Define input grid object for MPAS input data.
!!
!! @param [in] localpet ESMF local persistent execution thread
!! @param [in] npets  Number of persistent execution threads
!! @author  Larissa Reames CIWRO/NOAA/NSSL/FRDD

 subroutine define_input_grid(localpet,npets)
    integer, intent(in)  :: localpet, npets

    num_tiles_input_grid = 1 ! set as default to 1. FV3 input will change this

    if ( input_grid_type .eq. "MPAS") then
       call define_input_grid_mpas(localpet,npets)
    else if ( input_grid_type .eq. "FV3") then
       call define_input_grid_fv3(localpet,npets)
    else if ( input_grid_type .eq. "LATLON" ) then
       call define_input_grid_latlon(localpet,npets)
    endif
 end subroutine define_input_grid

 subroutine define_input_grid_mpas(localpet,npets)

 use mpi
 use netcdf
 implicit none

 character(len=500)           :: the_file, dimname

 integer, intent(in)          :: localpet, npets

 integer                      :: error, i, j, k, rc, n
 integer                      :: ncid,id_var, id_dim, dimsize, nVertThis
 integer                      :: nCells, nVertices, maxEdges, dimids(2)
 integer                      :: actual_maxEdges
 integer                      :: cell_start, cell_end, dims(3), temp1, &
                                   temp2, temp3, temp(1), clb(1), cub(1)
 integer, allocatable         :: elemTypes2(:), vertOnCell(:,:), &
                                          nodesPET(:), nodeIDs_temp(:), &
                                          elementConn_temp(:), elementConn(:)
 integer, allocatable         :: nEdgesOnCell(:)
 real(esmf_kind_r8), allocatable       :: latCell(:), lonCell(:), &
                                          latVert(:), lonVert(:), &
                                          nodeCoords(:), &
                                          nodeCoords_temp(:), &
                                          elemCoords(:)
 real(esmf_kind_r8), parameter         :: PI=4.D0*DATAN(1.D0)

 the_file = file_input_grid_mpas

 if (localpet==0) print*,'- OPEN MPAS INPUT FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid) ! CSS
 if (error /=0) call error_handler("OPENING MPAS INPUT FILE",error)

 !Get nCells size
 if (localpet==0) print*,'- READ nCells'
 error = nf90_inq_dimid(ncid,'nCells', id_dim)
 call netcdf_err(error, 'reading nCells id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nCells)
 call netcdf_err(error, 'reading nCells')

 nCells_input = nCells

 !Get nVertices size
 if (localpet==0) print*,'- READ nVertices'
 error = nf90_inq_dimid(ncid,'nVertices',id_dim)
 call netcdf_err(error, 'reading nVertices id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nVertices)
 call netcdf_err(error, 'reading nVertices')

 nVert_input = nVertices

  ! original algorithm from Larissa, which can lead to
  !  nCellsPerPET < 0 for some meshes,
  !  probably because of using "ceiling()"

!nCellsPerPET = ceiling(real(nCells)/real(npets))
!cell_start = localpet*nCellsPerPET+1
!cell_end = min(localpet*nCellsPerPET+nCellsPerPET,nCells)
!nCellsPerPET = cell_end - cell_start + 1

 ! Get number of MPAS cells per processor, and starting/ending cells.
 call para_range(1, nCells, npets, localpet, cell_start, cell_end)
 nCellsPerPET = cell_end - cell_start + 1

 !Get nVertLevels size
 if (localpet==0) print*,'- READ nVertLevels'
 error = nf90_inq_dimid(ncid,'nVertLevels',id_dim)
 call netcdf_err(error, 'reading nVertLevels id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nz_input)
 call netcdf_err(error, 'reading nVertLevels')

 !Get nVertLevelsP1 size
 if (localpet==0) print*,'- READ nVertLevelsP1'
 error = nf90_inq_dimid(ncid,'nVertLevelsP1',id_dim)
 call netcdf_err(error, 'reading nVertLevelsP1 id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nzp1_input)
 call netcdf_err(error, 'reading nVertLevelsP1')

 !Get maxEdges size
 if (localpet==0) print*,'- READ maxEdges'
 error = nf90_inq_dimid(ncid,'maxEdges',id_dim)
 call netcdf_err(error, 'reading maxEdges id')

 error=nf90_inquire_dimension(ncid,id_dim,len=maxEdges)
 call netcdf_err(error, 'reading maxEdges')

 !Get nSoilLevels size
 if (localpet==0) print*,'- READ nSoilLevels'
 error = nf90_inq_dimid(ncid,'nSoilLevels',id_dim)
 call netcdf_err(error, 'reading nSoilLevels id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nsoil_input)
 call netcdf_err(error, 'reading nSoilLevels')
 
 allocate(latCell(cell_start:cell_end))
 allocate(lonCell(cell_start:cell_end))
 allocate(latVert(nVertices))
 allocate(lonVert(nVertices))
 allocate(vertOnCell(maxEdges,cell_start:cell_end))

 ! GET CELL CENTER LAT/LON
 if (localpet==0) print*,'- READ LONCELL ID'
 error=nf90_inq_varid(ncid, 'lonCell', id_var)
 call netcdf_err(error, 'reading lonCell id')

 if (localpet==0) print*,'- READ LONCELL'
 error=nf90_get_var(ncid, id_var, start=(/cell_start/), count=(/nCellsPerPET/),values=lonCell)
 call netcdf_err(error, 'reading lonCell')

 if (localpet==0) print*,'- READ LATCELL ID'
 error=nf90_inq_varid(ncid, 'latCell', id_var)
 call netcdf_err(error, 'reading latCell id')

 if (localpet==0) print*,'- READ LATCELL'
  error=nf90_get_var(ncid, id_var, start=(/cell_start/), count=(/nCellsPerPET/),values=latCell)
 call netcdf_err(error, 'reading latCell')

 ! GET VERTEX LAT/LON
 if (localpet==0) print*,'- READ LONVERTEX ID'
 error=nf90_inq_varid(ncid, 'lonVertex', id_var)
 call netcdf_err(error, 'reading lonVertex id')

 if (localpet==0) print*,'- READ LONVERTEX'
 error=nf90_get_var(ncid, id_var, start=(/1/),count=(/nVertices/),values=lonVert)
 call netcdf_err(error, 'reading lonVertex')
 
 if (localpet==0) print*,'- READ LATVERTEX ID'
 error=nf90_inq_varid(ncid, 'latVertex', id_var)
 call netcdf_err(error, 'reading latVertex id')

 if (localpet==0) print*,'- READ LATVERTEX'
 error=nf90_get_var(ncid, id_var,  start=(/1/),count=(/nVertices/),values=latVert)
 call netcdf_err(error, 'reading latVertex')

 if (localpet==0) print*,"- NUMBER OF CELLS ON INPUT GRID ", nCells_input
 if (localpet==0) print*,"- NUMBER OF NODES ON INPUT GRID ", nVert_input

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if (localpet==0) print*,'- READ verticesOnCell ID'
 error=nf90_inq_varid(ncid, 'verticesOnCell', id_var)
 call netcdf_err(error, 'reading verticesOnCell id')

!error = nf90_inquire_variable(ncid,id_var,dimids=dims)
!call netcdf_err(error, 'reading verticesOnCell dims')

 if (localpet==0) print*,'- READ verticesOnCell'
 error=nf90_get_var(ncid, id_var, start=(/1,cell_start/),count=(/maxEdges,nCellsPerPET/),values=vertOnCell)
 call netcdf_err(error, 'reading verticesOnCell')

 ! Start CSS...it's possible that maxEdges in the dimensions could be something like 10, but 
 !  actual number of maxEdges is really 6. Moreover, in this example, columns 7-10 could have
 !  values > 0 for vertOnCell that are meaningless, but the code below for mesh generation won't like that.
 ! So, figure out the actual number of maxEdges in the mesh.
 ! Also, below, during mesh generation, check to see if we are > nEdgesOnCell for the ith point
 allocate(nEdgesOnCell(nCells_input)) ! each processor has the whole mesh...probably not needed
 if (localpet==0) print*,'- READ nEdgesOnCell ID'
 error=nf90_inq_varid(ncid, 'nEdgesOnCell', id_var)
 call netcdf_err(error, 'reading nEdgesOnCell id')
 if (localpet==0) print*,'- READ nEdgesOnCell'
 error=nf90_get_var(ncid, id_var, start=(/1/),count=(/nCells_input/),values=nEdgesOnCell)
 call netcdf_err(error, 'reading nEdgesOnCell')
 actual_maxEdges = maxval(nEdgesOnCell)  !CSS
 if ( localpet == 0 .and. ( actual_maxEdges /= maxEdges)) then
    write(*,*)'Actual number of maxEdges = ',actual_maxEdges
 endif
 ! End CSS

 error = nf90_close(ncid)

 ! Allocate and fill element corner coordinate array.
 allocate(nodeCoords_temp(2*maxEdges*nCellsPerPET))
 allocate(nodeIDs_temp(maxEdges*nCellsPerPET))
 allocate(elementConn_temp(maxEdges*nCellsPerPET))
 allocate(elemCoords(2*nCellsPerPET))
 allocate(elemTypes2(nCellsPerPET))
 allocate(elemIDs(nCellsPerPET))

 nVertThis = 0
 k = 0

 if (localpet==0) print*, "- Create PET-local element connectivity "

 if ( 1 == 1 ) then ! newest/latest method from MPASSIT

    do i = cell_start, cell_end
       j = i - cell_start + 1
       elemIDs(j) = i
    enddo

    do i = 1,nCellsPerPET
       j = elemIDs(i)
       elemTypes2(i) = count(vertOnCell(:,j)/=0)
       nVertThis = nVertThis + elemTypes2(i)
       elemCoords(2*i-1) = lonCell(j)*180.0_esmf_kind_r8/PI
       if (elemCoords(2*i-1) > 180.0_esmf_kind_r8) then
          elemCoords(2*i-1) = elemCoords(2*i-1) - 360.0_esmf_kind_r8
       endif
       elemCoords(i*2) = latCell(j)*180.0_esmf_kind_r8/PI
       elementConn_temp(maxEdges*(i-1)+1:maxEdges*i) = vertOnCell(:,j)
    enddo
    call unique_sort(elementConn_temp,maxEdges*nCellsPerPET,nodeIDs)
    nNodesPerPET=size(nodeIDs)
    allocate(nodeCoords(2*nNodesPerPET))
    do j = 1,nNodesPerPET
       i = nodeIDs(j)
       nodeCoords(2*j-1) = lonVert(i)*180.0_esmf_kind_r8/PI
       if (nodeCoords(2*j-1) > 180.0_esmf_kind_r8) then
          nodeCoords(2*j-1) = nodeCoords(2*j-1) - 360.0_esmf_kind_r8
       endif
       nodeCoords(2*j) = latVert(i)*180.0_esmf_kind_r8/PI
    enddo
    allocate(elementConn(nVertThis))
    nVertThis = 0
    do i = 1,nCellsPerPET
       k = elemIDs(i)
       do j = 1,maxEdges
         !if(vertOnCell(j,i)/=0) then
         !temp = FINDLOC(nodeIDs,vertOnCell(j,elemIDs(i)))
          if(vertOnCell(j,k)/=0) then
             temp = FINDLOC(nodeIDs,vertOnCell(j,k))
             elementConn(nVertThis+1) = temp(1)
             nVertThis = nVertThis + 1
          endif
       enddo
    enddo

 else ! original method

 do i = cell_start,cell_end
    j = i - cell_start + 1
    elemIDs(j) = i
    elemTypes2(j) = 0
    elemCoords(2*j-1) = lonCell(i)*180.0_esmf_kind_r8/PI
    if (elemCoords(2*j-1) > 180.0_esmf_kind_r8) then
        elemCoords(2*j-1) = elemCoords(2*j-1) - 360.0_esmf_kind_r8
    endif
    elemCoords(2*j) = latCell(i)*180.0_esmf_kind_r8/PI
    do n = 1,maxEdges
        if ( n > actual_maxEdges ) cycle ! CSS
        if ( n > nEdgesOnCell(i) ) cycle ! CSS
        if (vertOnCell(n,i)>0) then
            nVertThis = nVertThis + 1

            ! Make sure we don't duplicate nodeIDs or nodeCoords on any PET
            if (.not. any(nodeIDs_temp(1:k)==vertOnCell(n,i))) then
                k = k+1
                nodeCoords_temp(2*k-1) = &
                    lonVert(vertOnCell(n,i)) *180.0_esmf_kind_r8/PI
                if (nodeCoords_temp(2*k-1) > 180.0_esmf_kind_r8) then
                    nodeCoords_temp(2*k-1)  =  &
                        nodeCoords_temp(2*k-1) - 360.0_esmf_kind_r8
                endif
                nodeCoords_temp(2*k) = &
                    latVert(vertOnCell(n,i))*180.0_esmf_kind_r8/PI

                nodeIDs_temp(k) = vertOnCell(n,i)
            endif

            elemTypes2(j) = elemTypes2(j) + 1

            temp = FINDLOC(nodeIDS_temp, vertOnCell(n,i))
            elementConn_temp(nVertThis) = temp(1)

        endif

    enddo
 enddo

 allocate(nodeCoords(2*k), nodeIDs(k), elementConn(nVertThis))
 nodeCoords = nodeCoords_temp(1:k*2)
 nodeIDs = nodeIDs_temp(1:k)
 elementConn = elementConn_temp(1:nVertThis)

 endif ! methods

 write(*,*) localpet, cell_start, cell_end, nCellsPerPET

 if (localpet==0) print*, "- CREATE MESH -"
 input_mesh = ESMF_MeshCreate(parametricDim=2, &
                     spatialDim=2, &
                     nodeIDs= nodeIDs, &
                     nodeCoords = nodeCoords, &
                     elementIDs = elemIDs, &
                     elementTypes=elemTypes2, &
                     elementConn = elementConn, &
                     elementCoords=elemCoords, &
                     coordSys=ESMF_COORDSYS_SPH_DEG, &
                     rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN MeshCreate", rc)

 ! After the creation we are through with the arrays, so they may be deallocated.
 deallocate(elemCoords,elemTypes2)
 deallocate(nodeCoords_temp, nodeCoords)
 deallocate(elementConn_temp, elementConn)
 deallocate(nodeIDs_temp)
 deallocate(lonCell, latCell, lonVert, latVert,vertOnCell)
 deallocate(nEdgesOnCell) ! CSS

 end subroutine define_input_grid_mpas

!! @param [in] localpet ESMF local persistent execution thread
!! @param [in] npets Number of persistent execution threads
!! @author chgres_cube folks...this code is taken directly from there
 subroutine define_input_grid_fv3(localpet,npets)

 use netcdf

 implicit none

 integer, intent(in)          :: localpet, npets

 character(len=500)           :: the_file

 integer                      :: id_tiles, id_dim, tile
 integer                      :: extra, error, ncid
 integer, allocatable         :: decomptile(:,:)

 if(localpet==0) print*,'- OPEN INPUT GRID MOSAIC FILE: ',trim(mosaic_file_input_grid)
 error=nf90_open(trim(mosaic_file_input_grid),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening grid mosaic file')

 if(localpet==0) print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 call netcdf_err(error, 'reading ntiles id')
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles_input_grid)
 call netcdf_err(error, 'reading ntiles')

 error = nf90_close(ncid)

 if(localpet==0) print*,'- NUMBER OF TILES, INPUT MODEL GRID IS ', num_tiles_input_grid

 if (mod(npets,num_tiles_input_grid) /= 0) then
   call error_handler("MUST RUN WITH A TASK COUNT THAT IS A MULTIPLE OF 6.", 1)
 endif

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles_input_grid

 allocate(decomptile(2,num_tiles_input_grid))

 do tile = 1, num_tiles_input_grid
   decomptile(:,tile)=(/1,extra/)
 enddo

 if(localpet==0) print*,"- CALL GridCreateMosaic FOR INPUT MODEL GRID"
 ! This defines and calculates all of the dimensions of the grid
 input_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file_input_grid), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_input_grid), &
                                  rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", error)

!-----------------------------------------------------------------------
! Get the number of lat/lons.
!-----------------------------------------------------------------------
 the_file = trim(orog_dir_input_grid) // trim(orog_files_input_grid(1))

 if(localpet==0) print*,'- OPEN FIRST INPUT GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening ororgraphy file')
 if(localpet==0) print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, 'reading lon id')
 error=nf90_inquire_dimension(ncid,id_dim,len=i_input)
 call netcdf_err(error, 'reading lon')
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, 'reading lat id')
 error=nf90_inquire_dimension(ncid,id_dim,len=j_input)
 call netcdf_err(error, 'reading lat')
 error = nf90_close(ncid)

 if(localpet==0) print*,"- I/J DIMENSIONS OF THE INPUT GRID TILES ", i_input, j_input

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 deallocate(decomptile)

 end subroutine define_input_grid_fv3

 subroutine define_input_grid_latlon(localpet,npets)

 use netcdf

 implicit none

 integer, intent(in)          :: localpet, npets

 character(len=500)           :: the_file

 integer                      :: id_dim, id_var, i, j
 integer                      :: extra, error, ncid
 integer                      :: clb(2), cub(2), starts(2), counts(2)
 real(esmf_kind_r8)           :: half_dx_in_degrees
 real(esmf_kind_r8), allocatable       :: templat(:), templon(:)
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:), &
                                          clat_src_ptr(:,:), clon_src_ptr(:,:)
 type(esmf_polekind_flag)              :: polekindflag(2)

 half_dx_in_degrees = 0.5*dx_in_degrees

 the_file = file_input_grid_latlon

 if(localpet==0) print*,'- OPEN LAT-LON FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening file')

 if(localpet==0) print*,"- READ GRID DIMENSIONS"

 error=nf90_inq_dimid(ncid, trim(adjustl(lon_var)), id_dim)
 call netcdf_err(error, 'reading lon id')
 error=nf90_inquire_dimension(ncid,id_dim,len=i_input)
 call netcdf_err(error, 'reading lon')

 error=nf90_inq_dimid(ncid, trim(adjustl(lat_var)), id_dim)
 call netcdf_err(error, 'reading lat id')
 error=nf90_inquire_dimension(ncid,id_dim,len=j_input)
 call netcdf_err(error, 'reading lat')

 if(localpet==0) print*,"- I/J DIMENSIONS OF THE INPUT GRID TILES ", i_input, j_input

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 ! basically follow subroutine define_input_grid_gaussian in chgres_cube/model_grid.F90
 if (is_regional) then
    if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR INPUT GRID"
    input_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_input,j_input/), &
                                       indexflag=ESMF_INDEX_GLOBAL, &
                                       regDecomp=(/1,npets/), &
                                       rc=error)
 else
    polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
    if (localpet==0)print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."
    input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
           maxIndex=(/i_input,j_input/), &
           polekindflag=polekindflag, &
           periodicDim=1, &
           poleDim=2,  &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           regDecomp=(/1,npets/),  &
           indexflag=ESMF_INDEX_GLOBAL, rc=error)

 endif
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridCreateNoPeriDim", error)

!-----------------------------------------------------------------------
! Create needed grid coordinates
!-----------------------------------------------------------------------
    
 if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CENTER."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
    
 if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CORNER."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
    
!----------------------------------------------------------
!  Read in coordinate values and set on grid
!----------------------------------------------------------
    
!------- Grid center coordinates
    
 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=1, &
                      farrayPtr=lon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=2, &
                      computationalLBound=clb, &
                      computationalUBound=cub, &
                      farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN GridGetCoord", error)

 allocate(templon(clb(1):cub(1)))
 allocate(templat(clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, trim(adjustl(lon_var)), id_var)
 call netcdf_err(error, 'reading longitude id')

 if (localpet==0) print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, start=(/starts(1)/),count=(/counts(1)/),values=templon)
 call netcdf_err(error, 'reading longitude')
 
 if (localpet==0) print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, trim(adjustl(lat_var)), id_var)
 call netcdf_err(error, 'reading latitude id')

 if (localpet==0) print*,'- READ LATITUDE'
 error=nf90_get_var(ncid, id_var, start=(/starts(2)/),count=(/counts(2)/),values=templat)
 call netcdf_err(error, 'reading latitude')

 do j = clb(2),cub(2)
   do i = clb(1), cub(1)
     lon_src_ptr(i,j)=real(templon(i),esmf_kind_r8)
     if (lon_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_src_ptr(i,j) = lon_src_ptr(i,j) - 360.0_esmf_kind_r8
     lat_src_ptr(i,j)=real(templat(j),esmf_kind_r8)
   enddo
 enddo

 nullify(lon_src_ptr)
 nullify(lat_src_ptr)
  
!---------- Grid corners coordinate creation
 nullify(clat_src_ptr)
 nullify(clon_src_ptr)

 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID X-COORD."
   call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=clon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID Y-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=clat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridGetCoord", error)
      
 ! For a lat/lon grid, it's easy to get the corner lat/lons from the cell center lat/lons
 ! clb, cub are bounds are for the corner points
!print *,localpet, clb(1), cub(1), clb(2), cub(2)
 do j = clb(2),cub(2)
   do i = clb(1), cub(1)
     if ( i .eq. ip1_input .and. j .eq. jp1_input ) then ! top-right corner
        clon_src_ptr(i,j)=real(templon(i_input)+half_dx_in_degrees,esmf_kind_r8)
        clat_src_ptr(i,j)=real(templat(j_input)+half_dx_in_degrees,esmf_kind_r8)
     else if ( i .eq. ip1_input ) then ! points on right edge
        clon_src_ptr(i,j)=real(templon(i_input)+half_dx_in_degrees,esmf_kind_r8)
        clat_src_ptr(i,j)=real(templat(j)-half_dx_in_degrees,esmf_kind_r8)
     else if ( j .eq. jp1_input ) then ! points on top edge
        clon_src_ptr(i,j)=real(templon(i)-half_dx_in_degrees,esmf_kind_r8)
        clat_src_ptr(i,j)=real(templat(j_input)+half_dx_in_degrees,esmf_kind_r8)
     else ! all other points
        clon_src_ptr(i,j)=real(templon(i)-half_dx_in_degrees,esmf_kind_r8)
        clat_src_ptr(i,j)=real(templat(j)-half_dx_in_degrees,esmf_kind_r8)
     endif
     if (clon_src_ptr(i,j) > 360.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) - 360.0_esmf_kind_r8
     if (clon_src_ptr(i,j) <   0.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) + 360.0_esmf_kind_r8
     if (clat_src_ptr(i,j) >  90.0_esmf_kind_r8) clat_src_ptr(i,j) = 90.0_esmf_kind_r8
     if (clat_src_ptr(i,j) < -90.0_esmf_kind_r8) clat_src_ptr(i,j) = -90.0_esmf_kind_r8
   enddo
 enddo

 nullify(clon_src_ptr)
 nullify(clat_src_ptr)
 deallocate(templat,templon)
 
 error = nf90_close(ncid)

 end subroutine define_input_grid_latlon

!> Setup the esmf grid object for the target grid.
!!
!! @param [in] localpet ESMF local persistent execution thread
!! @param [in] npets Number of persistent execution threads
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine define_target_grid(localpet, npets)
 
 use netcdf
 use mpi

 implicit none

 character(len=500)           :: the_file

 integer, intent(in)          :: localpet, npets

 integer                      :: error, extra, i, j, clb(2), cub(2), &
                                 starts(2), counts(2)

 real(esmf_kind_r8), allocatable       :: latitude(:,:), longitude(:,:), &
                                          templat(:,:), templon(:,:)
 integer                               :: ncid,id_var, id_dim, id_clon, id_clat
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:), &
                                          clat_src_ptr(:,:), clon_src_ptr(:,:)
 logical :: file_has_lat_corners = .false.
 logical :: file_has_lon_corners = .false.

 type(esmf_polekind_flag)              :: polekindflag(2)

 the_file = file_target_grid

 if (localpet==0) print*,'- OPEN WRF INPUT FILE: ',trim(the_file)
!error=nf90_open_par(trim(the_file),NF90_NOWRITE,MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid) ! CSS

 if (error /=0) call error_handler("OPENING WRF INPUT FILE",error)

 if (localpet==0) print*,'- READ WEST_EAST ID'
 error=nf90_inq_dimid(ncid, 'west_east', id_dim)
 call netcdf_err(error, 'reading west_east id')

 if (localpet==0) print*,'- READ WEST_EAST'
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 call netcdf_err(error, 'reading west_east')

 if (localpet==0) print*,'- READ SOUTH_NORTH ID'
 error=nf90_inq_dimid(ncid, 'south_north', id_dim)
 call netcdf_err(error, 'reading south_north id')

 if (localpet==0) print*,'- READ SOUTH_NORTH'
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 call netcdf_err(error, 'reading south_north')

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE DX'
 error = nf90_get_att(ncid,NF90_GLOBAL,'DX',dx)
 call netcdf_err(error, 'reading dx')

 if (localpet==0) print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

 if (localpet==0) print*, '- READING GLOBAL ATTRIBUTES'
 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LAT', cen_lat)
 call netcdf_err(error, 'GETTING CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LON', cen_lon)
 call netcdf_err(error, 'GETTING CEN_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT1', truelat1)
 call netcdf_err(error, 'GETTING TRUELAT1 GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT2', truelat2)
 call netcdf_err(error, 'GETTING TRUELAT2 GLOBAL ATTRIBUTE')

!error = nf90_get_att(ncid, NF90_GLOBAL, 'MOAD_CEN_LAT', cen_lat)
!call netcdf_err(error, 'GETTING MOAD_CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'STAND_LON', stand_lon)
 call netcdf_err(error, 'GETTING STAND_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
 call netcdf_err(error, 'GETTING POLELAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LON', pole_lon)
 call netcdf_err(error, 'GETTING POLE_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ', proj_code)
 call netcdf_err(error, 'GETTING MAP_PROJ GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ_CHAR', map_proj_char)
 if (error .ne. 0) then
   if (proj_code == 1) then
     map_proj_char = "Lambert Conformal"
   else if (proj_code == 2) then
     map_proj_char = "Polar Stereographic"
   else
     map_proj_char = "Lat/Lon"
   endif
 endif
 
!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if ( target_grid_is_regional ) then
    if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR TARGET GRID"
    target_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_target,j_target/), &
                                       indexflag=ESMF_INDEX_GLOBAL, &
                                       rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN GridCreateNoPeriDim", error)
 else
    polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
    if (localpet==0)print*,"- CALL GridCreate1PeriDim FOR TARGET GRID."
    target_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
           maxIndex=(/i_target,j_target/), &
           polekindflag=polekindflag, &
           periodicDim=1, &
           poleDim=2,  &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           regDecomp=(/1,npets/),  &
           indexflag=ESMF_INDEX_GLOBAL, rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN GridCreate1PeriDim", error)
 endif

!-----------------------------------------------------------------------
! Create needed grid coordinates
!-----------------------------------------------------------------------
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR TARGET GRID CENTER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR TARGET GRID CORNER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
    
!----------------------------------------------------------
!  Read in coordinate values and set on grid
!----------------------------------------------------------
    
!------- Grid center coordinates
    
  if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID X-COORD."
  nullify(lon_src_ptr)
  call ESMF_GridGetCoord(target_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=1, &
                      farrayPtr=lon_src_ptr, rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN GridGetCoord", error)

  if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID Y-COORD."
  nullify(lat_src_ptr)
  call ESMF_GridGetCoord(target_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=2, &
                      computationalLBound=clb, &
                      computationalUBound=cub, &
                      farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN GridGetCoord", error)

 allocate(templat(clb(1):cub(1),clb(2):cub(2)))
 allocate(templon(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'XLONG', id_var)
 if (error /= NF90_NOERR) then
   error=nf90_inq_varid(ncid, 'XLONG_M', id_var)
   call netcdf_err(error, 'reading longitude id')
 endif

 if (localpet==0) print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templon)
 call netcdf_err(error, 'reading longitude')
 
 if (localpet==0) print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'XLAT', id_var)
 if (error .ne. NF90_NOERR) then
   error=nf90_inq_varid(ncid, 'XLAT_M', id_var)
   call netcdf_err(error, 'reading latitude id')
 endif

 if (localpet==0) print*,'- READ LATITUDE'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templat)
 call netcdf_err(error, 'reading latitude')
    
 do j = clb(2),cub(2)
   do i = clb(1), cub(1)
     lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
     lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
   enddo
 enddo

 nullify(lon_src_ptr)
 nullify(lat_src_ptr)
  
 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
    call error_handler("IN FieldCreate", error)
    
 call ESMF_FieldGet(latitude_target_grid, farrayptr=lat_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
      call error_handler("IN FieldGet", error)
 call ESMF_FieldGet(longitude_target_grid, farrayptr=lon_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
      call error_handler("IN FieldGet", error)
 
 do j = clb(2),cub(2)
 do i = clb(1), cub(1)
    lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
    lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
 enddo
 enddo  

 nullify(lon_src_ptr)
 nullify(lat_src_ptr)
  
!---------- Grid corners coordinate creation

  if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET CORNERS GRID X-COORD."
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=1, &
                          farrayPtr=clon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET CORNERS GRID Y-COORD."
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=clat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)
      
  ! Determine if the corner lats/lons are in the file or if we need to 
  !  calculate them
  error=nf90_inq_varid(ncid, 'XLAT_C', id_clat)
  if (error .eq. NF90_NOERR) file_has_lat_corners = .true.
  error=nf90_inq_varid(ncid, 'XLONG_C', id_clon)
  if (error .eq. NF90_NOERR) file_has_lon_corners = .true.

  if ( file_has_lat_corners .and. file_has_lon_corners ) then
     if (localpet==0) print*,'- GETTING CORNER LAT/LONS FROM FILE'
     deallocate(templat,templon)
     ! clb, cub are for the corner grid points
     allocate(templat(clb(1):cub(1),clb(2):cub(2)))
     allocate(templon(clb(1):cub(1),clb(2):cub(2)))
     starts = (/clb(1),clb(2)/)
     counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
     error=nf90_get_var(ncid, id_clat, start=starts,count=counts,values=templat)
     call netcdf_err(error, 'reading XLAT_C')
     error=nf90_get_var(ncid, id_clon, start=starts,count=counts,values=templon)
     call netcdf_err(error, 'reading XLONG_C')
    
     do j = clb(2),cub(2)
       do i = clb(1), cub(1)
          clon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
          clat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
          if (clon_src_ptr(i,j) > 360.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) - 360.0_esmf_kind_r8
          if (clon_src_ptr(i,j) <   0.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) + 360.0_esmf_kind_r8
          if (clat_src_ptr(i,j) >  90.0_esmf_kind_r8) clat_src_ptr(i,j) = 90.0_esmf_kind_r8
          if (clat_src_ptr(i,j) < -90.0_esmf_kind_r8) clat_src_ptr(i,j) = -90.0_esmf_kind_r8
       enddo
      enddo


   else
     if (localpet==0) print*,'- CALCULATING CORNER LAT/LONS'

     allocate(latitude(i_target,j_target))
     allocate(longitude(i_target,j_target))
     call ESMF_FieldGet(latitude_target_grid, farrayPtr=lat_src_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN FieldGet", error)
      
     call ESMF_FieldGet(longitude_target_grid, farrayPtr=lon_src_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call error_handler("IN FieldGet", error)

     call get_cell_corners(lat_src_ptr, lon_src_ptr, clat_src_ptr, clon_src_ptr, real(dx,esmf_kind_r8), clb, cub)
     deallocate(latitude,longitude)
  endif
    
  nullify(lon_src_ptr)
  nullify(lat_src_ptr)
  nullify(clon_src_ptr)
  nullify(clat_src_ptr)
  deallocate(templat,templon)
 
 error = nf90_close(ncid)

 end subroutine define_target_grid

 !> For grids with equal cell sizes (e.g., lambert conformal), get lat and on of the grid
!! cell corners
!!
!! @param [in]  latitude  2d array of grid latitude
!! @param [in]  longitude 2d array of grid longitude
!! @param [in]  dx grid cell size in meters
!! @param [in]  clb  computational lower bound as returned from ESMF_GridGetCoord
!! @param [in]  cub  computational upper bound as returned from ESMF_GridGetCoord
!! @param [out] latitude_sw  returned latitude at sw corner of grid cells
!! @param [out]elongitude_sw  returned longitude at sw corner of grid cells
!! @author Larissa Reames CIWRO/NSSL/FRDD

  subroutine get_cell_corners( latitude, longitude, latitude_sw, longitude_sw, dx,clb,cub)
  implicit none

  real(esmf_kind_r8), intent(in), pointer :: latitude(:,:)
  real(esmf_kind_r8), intent(inout), pointer   :: latitude_sw(:,:)
  real(esmf_kind_r8), intent(in),pointer    :: longitude(:,:)
  real(esmf_kind_r8), intent(inout), pointer   :: longitude_sw(:,:)
  real(esmf_kind_r8), intent(in)    :: dx !grid cell side size (m)

  integer, intent(in) :: clb(2), cub(2)

  real(esmf_kind_r8)                :: lat1, lon1, lat2, lon2, d, brng


  real(esmf_kind_r8), parameter    :: pi = 3.14159265359
  real(esmf_kind_r8), parameter    :: R =  6370000.0
  real(esmf_kind_r8), parameter    :: bearingInDegrees = 135.0

  integer                           :: i, j, e

  d = sqrt((dx**2.0_esmf_kind_r8)/2.0_esmf_kind_r8)

  do j = clb(2),cub(2)
   do i = clb(1), cub(1)
         if (j == jp1_target .and. i == ip1_target) then
             lat1 = latitude(i_target,j_target)  * ( pi / 180.0_esmf_kind_r8 )
             lon1 = longitude(i_target,j_target) * ( pi / 180.0_esmf_kind_r8 )
             brng = 315.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
             lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
             lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
             latitude_sw(ip1_target,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
             longitude_sw(ip1_target,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
             cycle
         endif

     if (i == ip1_target) then
       brng = 225.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i_target,j)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i_target,j) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(ip1_target,j) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(ip1_target,j) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     if (j == jp1_target) then
       brng = 45.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i,j_target)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i,j_target) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(i,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(i,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     lat1 = latitude(i,j)  * ( pi / 180.0_esmf_kind_r8 )
     lon1 = longitude(i,j) * ( pi / 180.0_esmf_kind_r8 )

     brng = bearingInDegrees * ( pi / 180.0_esmf_kind_r8 );
     lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
     lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );

     latitude_sw(i,j) = lat2 * 180.0_esmf_kind_r8 / pi
     longitude_sw(i,j) = lon2 * 180.0_esmf_kind_r8 / pi

   enddo
 enddo

 end subroutine get_cell_corners

 subroutine cleanup_input_target_grid_data(localpet)

 implicit none
 integer, intent(in)              :: localpet
 integer                          :: rc, i
 type(esmf_field), allocatable    :: fields(:)

 if (localpet==0) print*,"- DESTROY MODEL DATA."

 if ( esmf_input_mesh ) then
    deallocate(elemIDs, nodeIDs)
 endif

 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)

 allocate(fields(nfields))
 call ESMF_FieldBundleGet(input_bundle, fieldList=fields, &
                       itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                       rc=rc)
 do i = 1, nfields
     call ESMF_FieldDestroy(fields(i), rc=rc)
 enddo
 call ESMF_FieldBundleDestroy(input_bundle)

 call ESMF_FieldBundleGet(target_bundle, fieldList=fields, &
                       itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                       rc=rc)
 do i = 1, nfields
     call ESMF_FieldDestroy(fields(i), rc=rc)
 enddo
 call ESMF_FieldBundleDestroy(target_bundle)

 deallocate(fields)

 if (esmf_input_mesh) call ESMF_MeshDestroy(input_mesh, rc=rc)
 if (esmf_input_grid) call ESMF_GridDestroy(input_grid, rc=rc)

 call ESMF_GridDestroy(target_grid, rc=rc)

 deallocate(input_names, target_names)
 deallocate(interp_methods, vert_info)
 deallocate(nVertLevelsPerVariable)
 deallocate(target_units, target_longname)

 end subroutine cleanup_input_target_grid_data

 subroutine unique_sort(val,nvals, final)
   !implicit none
    integer, intent(in) :: nvals
    integer, intent(in) :: val(nvals)
    integer, intent(inout), allocatable :: final(:)
    integer :: i = 0, min_val, max_val
    integer, dimension(nvals) :: unique

    min_val = minval(val)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        min_val = minval(val, mask=val>min_val)
        if (min_val>0) then
          i = i+1
          unique(i) = min_val
        endif
    enddo
    allocate(final(i), source=unique(1:i))   !<-- Or, just use unique(1:i)
 end subroutine unique_sort


subroutine para_range(n1, n2, nprocs, irank, ista, iend)

  integer, intent(in) :: n1, n2, nprocs, irank
  integer, intent(out) :: ista, iend

  integer :: iwork1, iwork2

  iwork1 = (n2 - n1 + 1) / nprocs
  iwork2 = mod(n2 - n1 + 1, nprocs)
  ista = irank * iwork1 + n1 + min(irank, iwork2)
  iend = ista + iwork1 - 1
  if (iwork2 > irank) iend = iend + 1
  return

end subroutine para_range


 end module model_grid
