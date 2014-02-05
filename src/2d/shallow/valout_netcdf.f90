!
! -----------------------------------------------------
!     Routine to write netcdf files in netcdf-4 (HDF5) format
!        # dlgeorge-2014.02
!        # Each file written by the fortran code has
!        # Global Attributes:
!        #           time             : time of frame (also a variable in case time attributes are desired
!        #           amr_levels       : The number of AMR levels
!        #           num_grids        : Number of grids in this frame
!        #           num_vars         : total number of variables stored in a grid (eqn + aux)
!        #           q_components     : list indicating equations stored eg. [1,1,1]
!        #           aux_components   : list indicating auxilliary variables stored eg. [1,0,0]
!        #           num_dim          : number of spatial dimensions
!        #           coord_sys        : 'Cartesian' or 'Latitude-Longitude'
!        #           coord_units      : 'meters' or 'degrees'
!        #           coord_xbounds    : the domain xlower,xupper
!        #           coord_ybounds    : the domain ylower,yupper
!        # Dimensions:
!        #           num_vars         : variables stored dimension
!        #           dimx_<gridno>    : X dimension for grid number <gridno>
!        #           dimy_<gridno>    : Y dimension for grid number <gridno>
!        # Variables:
!        #           time              : Stores the time of the frame
!        #           grid_<gridno>     : A grid of (num_vars,dimx,dimy)
!        # Attributes:
!        # (grid_<no>) gridno          : The number of this grid <grid_no>
!        #             level           : The AMR level
!        #             num_vars        : The number of variables
!        #             dim<x,y>.m      : X and Y dimensions (ie. mx,my)
!        #             dim<x,y>.low    : The lowest dimension value
!        #             dim<x,y>.d      : The distance between grid points
! -----------------------------------------------------

      subroutine valout (lst, lend, time, nvar, naux)

      use amr_module
      use netcdf
      use geoclaw_module, only: coordinate_system

      implicit none

      !arguments
      integer, intent(in) ::lst,lend,nvar,naux
      real(kind=8), intent(in) :: time

      !locals
      real(kind=8) :: xlow,ylow,dx,dy
      real, allocatable  :: grid(:,:,:)
      character(10)::  matname1, matname2
      character(4) :: gridstr
      character(40) :: dim_names
      character(40) :: coord_units
      character(40) :: coord_sys
      logical :: outaux
      integer :: i,j,ivar,iaux,iq_store,iaux_store
      integer :: iadd,iaddaux
      integer :: ipos, nstp, matunit1, matunit2,idigit
      integer :: num_vars,output_aux_num,output_q_num
      integer :: rcode,ncid
      integer :: num_vars_id,time_id,dimxid,dimyid,gridid
      integer :: level,nx,ny,loc,locaux,mitot,mjtot,mptr


 !OLD INDEXING
     !iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
     !iaddaux(i,j,iaux) = locaux + i - 1 + mitot*((iaux-1)*mjtot+j-1)
 !NEW INDEXING
      iadd(ivar,i,j)  = loc + ivar-1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*((j-1)*mitot+i-1)

   ! ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
   !graphics output of soln values for contour or surface plots.
   !modified for GeoClaw to output the surface level along with q.
   !  surface = q(i,j,1) + aux(i,j,1)
   ! :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

   !###  make the file names and open output files
   matname1 = 'fort.qxxxx'
   matname2 = 'fort.txxxx'
   matunit1 = 50
   matunit2 = 60
   nstp     = matlabu

   do ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
   enddo ! 55

   ! # how many q components requested?
   output_q_num = 0
   do i=1,nvar
      output_q_num = output_q_num + output_q_components(i)
   enddo

   ! # how many aux components requested?
   output_aux_num = 0
   do i=1,naux
      output_aux_num = output_aux_num + output_aux_components(i)
   enddo

   outaux = ((output_aux_num > 0) .and. ((.not. output_aux_onlyonce) .or. (time==t0)))
   if (.not.outaux) output_aux_num = 0

   num_vars = output_aux_num + output_q_num

         !!!!Define netcdf file and global attributes and variables
   if (output_format == 2) then
         if (coordinate_system==2) then
            coord_sys = 'Latitude-Longitude'
            coord_units = 'degrees Lat-Lon'
         else
            coord_sys = 'Cartesian'
            coord_units = 'meters'
         endif
         !prepare netcdf4 dataset (for netcdf3 change NF90_NETCDF4 to NF90_NOCLOBBER)
         rcode = nf90_create(matname1//'.nc',NF90_NETCDF4,ncid)
         if(rcode.ne.NF90_NOERR) print *,'NETCDF ERROR OPENING NETCDF FILE'
         !define dimension for number of output variables
         rcode = nf90_def_dim(ncid,'num_vars',num_vars,num_vars_id)
         !define time variable
         rcode = nf90_def_var(ncid,'time',NF90_FLOAT,time_id)
         rcode = nf90_put_att(ncid,time_id,'units','seconds')
         !assign known global attributes
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'time',time)
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'num_dim',2)
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'amr_levels',mxnest)
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'num_vars',num_vars)
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'q_components',output_q_components)
         if (outaux) then
            rcode = nf90_put_att(ncid,NF90_GLOBAL,'aux_components',output_aux_components)
         else
            rcode = nf90_put_att(ncid,NF90_GLOBAL,'aux_components',0*output_aux_components)
         endif

         rcode = nf90_put_att(ncid,NF90_GLOBAL,'coord_sys',trim(coord_sys))
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'coord_units',trim(coord_units))
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'coord_xbounds',(/xlower,xupper/))
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'coord_ybounds',(/ylower,yupper/))
         rcode=nf90_enddef(ncid)
         if(rcode.ne.NF90_NOERR) print *,'ERROR  LEAVING DEFINE MODE'
         !assign time
         rcode = nf90_put_var(ncid,time_id,time)
   endif

         level = lst
         ngrids = 0
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)

              if (output_format==2) then
                  rcode=NF90_REDEF(ncid)
                  if(rcode.ne.NF90_NOERR) print *,'NETCDF ERROR REDEFINE MODE'
                  write(gridstr,'(I4.4)') mptr
                  !define dimensions for grid
                  rcode = nf90_def_dim(ncid,'dimx_'//trim(gridstr),nx,dimxid)
                  rcode = nf90_def_dim(ncid,'dimy_'//trim(gridstr),ny,dimyid)

                  !define grid variable
                  rcode = nf90_def_var(ncid,'grid_'//trim(gridstr),NF90_FLOAT, &
                               (/num_vars_id,dimxid,dimyid/),gridid)
                  if(rcode.ne.NF90_NOERR) print *,'NETCDF ERROR DEFINING GRID'
                  !assign grid attributes
                  rcode = nf90_put_att(ncid,gridid,'gridno',mptr)
                  rcode = nf90_put_att(ncid,gridid,'level',level)
                  dim_names="['num_vars','dimx','dimy']"
                  rcode = nf90_put_att(ncid,gridid,'dim_names',TRIM(dim_names))
                  rcode = nf90_put_att(ncid,gridid,'num_vars',num_vars)
                  rcode = nf90_put_att(ncid,gridid,'dimx.m',nx)
                  rcode = nf90_put_att(ncid,gridid,'dimy.m',ny)
                  rcode = nf90_put_att(ncid,gridid,'dimx.lower',xlow)
                  rcode = nf90_put_att(ncid,gridid,'dimy.lower',ylow)
                  dx = hxposs(level)
                  dy = hyposs(level)
                  rcode = nf90_put_att(ncid,gridid,'dimx.d',dx)
                  rcode = nf90_put_att(ncid,gridid,'dimy.d',dy)
                  !grid is defined, leave define mode
                  rcode = nf90_enddef(ncid)
                  !create and fill grid array
                  allocate(grid(num_vars,nx,ny))
                  do j=nghost+1,mjtot-nghost
                  do i=nghost+1,mitot-nghost
                     iq_store = 0
                     do ivar=1,nvar
                        if (dabs(alloc(iadd(ivar,i,j))) .lt. 1d-90) then
                           alloc(iadd(ivar,i,j)) = 0.d0
                        endif
                        if (output_q_components(ivar)>0) then
                           iq_store = iq_store + 1
                           grid(iq_store,i-nghost,j-nghost)=alloc(iadd(ivar,i,j))
                        endif
                     enddo
                     iaux_store = output_q_num
                     do iaux = 1, naux
                        if (outaux.and.output_aux_components(iaux)>0) then
                           iaux_store = iaux_store+1
                           grid(iaux_store,i-nghost,j-nghost) =  alloc(iaddaux(iaux,i,j))
                        endif
                     enddo
                  enddo
                  enddo
                  !save grid
                  rcode = nf90_put_var(ncid,gridid,grid,(/1,1,1/),(/num_vars,nx,ny/))
                  if(rcode.ne.NF90_NOERR) print *,'NETCDF ERROR  Writing Grid'
                  if (allocated(grid)) deallocate(grid)
               endif

            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

      if (output_format==1) then
         write(6,601) matlabu,time
      elseif (output_format==2) then
         !assign global attribute num_grids and variable time
         rcode = nf90_put_var(ncid,time_id,time)
         rcode = nf90_put_att(ncid,NF90_GLOBAL,'num_grids',ngrids)
         rcode = nf90_close(ncid)
         write(6,602) matlabu,time
      endif
      matlabu = matlabu + 1

  601 format('GeoClaw: Frame ',i4,' output files done at time t = ', d12.6,/)
  602 format('GeoClaw: Frame ',i4,' output written to netcdf file at time t = ', d12.6,/)

      return
      end subroutine
