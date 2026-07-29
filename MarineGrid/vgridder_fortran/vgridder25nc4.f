c   vgridder.f - greate the VDatum marine grid.           

c   Has inside_poly(ipr,np,x,y,xt,yt, in,n1,iselect)
c     w/o subs. intersect, capp_conxy, ins2
c   Has check of 4 corners (plus intermediate points) in each cell
c   Has new subroutines crosses, cross_segment.
c   ver  6: Allocatable memory. New maxmin routine.
c   ver  7: small changes to make compatible with LF95
c   ver  9: Has new inside polygon routines.
c           Also, a new record in input file for output of a 
c           GTX-formatted file (Oct 2009)
c   ver 10: If no BP, the point is automatically inside, not outside. 
c           Also, include check for barriers (nobarr1), and when adding
c           extra layers, also can check for barriers (nobarr2) and whether
c           point is inside BP (icheckbp). Remove the grid size reset 
c           option (ireset). Revise how alon1,2 and alat1,2 are used to
c   ver 11: uses additional segments in BP to exclude areas from
c           being water points.
c   ver 12: Increase array sizes isgx,icrx. Exclude add'nl layers
c           in added BP segments.
c   ver 13: Corrected version, to use multiple bounding polygons (v13d.f).
c   ver 14: Revised prints for clarity.
c           Only one call to Sub. remove_barriers, after add layers
c           Remove printi, use pr_int.
c   ver 15: if added layers, new points have an index value 
c           of 2, not 1, in the marine grid. (11-18-2011)
c           Also, revised ponding.
c   ver 16: Select as inside a polygon same way for coastline, BP, NT.
c           New pril() subroutine to print integer arrays.
c           Adjust user-specified pts after init land/water determination
c           (sub. adjust)
c           Use inside_pert
c           (09.04.2012)
c   ver 17: Use EEZ polygon to exclude points outside US waters.
c           Revise rd_xyi for case when longest segment isn't first.
c           b.Exclude extension layer points that fall outside BP
c           c.cosmetic changes in print out
c   ver 18: a.remove ponds by actual points (after check for BP)
c   ver 19: correct the zeroing of points outside BP (sub. ck_bp).           
c   ver 19a:accept jfield(,) values of -1, which prevents filling.
c           ab:activate pond removal.
c   ver 20: Read ADCIRC datum file, find points that have dried.
c           a:select_water
c   ver 22: cleaned up vg21.f. revise sub. rem_barrier (Jan 2020)
c---------------------------------------------------------------------
      module saved_data
      character   ,save                       :: pname*14=
     *                                           'vgridder24.f'
c        file names
      character, save                         :: filenm(10)*100
c        input parameters
      integer, save                           :: layers,
     *                                           nobarr1,nobarr2
c        coastline and polygons
      real*8,allocatable,dimension(:),save    :: xmins,xmaxs,
     *                                           ymins,ymaxs,xminb,xmaxb
     *                                           ,yminb,ymaxb
      integer,allocatable,dimension(:),save   :: isend,isx,isy,ipt,
     *                                           iptij,npoly
      integer, save                           :: ireset,ncmx,nsegs,
     *                                           islong
c        temporary coastline
      real*8,allocatable,dimension(:),save    :: xx,yy,xi,yi
      integer,allocatable,dimension(:),save   :: ipet,ifs,ipex
      integer, save                           :: imx,nsegt
c        EEZ
      real*8,allocatable,dimension(:),save    :: xe,ye
      integer,allocatable,dimension(:),save   :: ipene,ifse
      integer, save                           :: iemax,nsege
c        coastline
      real*8,allocatable,dimension(:,:),save  :: xmc
      real*8,allocatable,dimension(:),save    :: xcs,ycs
      integer,allocatable,dimension(:),save   :: ipens,isfrst,ippx
      integer, save                           :: ismax,nsegc
c        bounding polygon
      real*8,allocatable,dimension(:,:),save  :: xmb
      real*8,allocatable,dimension(:),save    :: xbc,ybc
      integer,allocatable,dimension(:),save   :: ipenb,ibfrst
      integer, save                           :: ibmax,nsegb,npoly1
c        marine grid
      real*8,allocatable,dimension(:,:),save  :: xv,yv,isv
      real*8,allocatable,dimension(:),save    :: clat,clon
      real*8, save                            :: delx,dely,alat0,alon0,
     *                                           alat1,alat2,alon1,
     *                                           alon2,alat3,alon3,
     *                                           raddeg=0.01745329251994
      integer,allocatable,dimension(:,:),save :: jfield,idry
      integer, save                           :: imax,jmax,iset
c        ADCIRC grid                                  
      integer,save                            :: nodes,nt,ngrids
      real*8,allocatable,dimension(:),save    :: xs,ys
      real*8,allocatable,dimension(:,:),save  :: ds  ! (node,idatum)
      integer,allocatable,dimension(:,:),save :: ncellnodes,icells
      integer,save                            :: iallocsor,numnulls,
     *                                           nblock
c        adjusted points
      integer, save                           :: nadjust,ielim_ponds
      integer,allocatable,dimension(:,:),save :: iadj1,ijpon
c        print out
      real*8, save                            :: aln(2),al4(4)
      integer, dimension(20), save            :: iprt,jprt
      integer, save                           :: lm(6),lm2(6),mpr
c        coastline crossing
      parameter (icrx=200) ! potential number of crossings in one
                           ! coastline segment
      parameter (isgx=100) ! potential number of segments found for one
                           ! latitude
      end module saved_data
c-----------------------------------------------------------------------

c     vgridder.f - Select points inside lat/lon window and bounding 
c                  polygon (BP) that are water for VDatum grid. 

c======================================================================

c     programmer: Kurt Hess    	
c                 Coast Survey Development Laboratory
c                 Office of Coast Survey, National Ocean Service, NOAA
c                 1315 East-West Highway, Silver Spring, Maryland 20910
c                 kurt.hess@noaa.gov
c     date :      v.6: Jan 2006
c                 last revised Jan 2020
c     language:   Fortran 95

c=======================================================================

c     compile:   ifort vgridder.f -o vgridder.x

c     run:       vgridder.x < vgridder.in > prvg

c     input: vgridder.in
c       title
c       4, y1,y2,x1,x2                  no. of pts,lat/lon window (degrees)
c       delx,dely                       point spacing (degrees)
c       ielim_ponds                     index to remove ponds(0=no,=yes)
c       nadjust                         num points adjusted before layers added
c       n, i,j,k, i,j,k, i,j,k          points: i,j,k=new jfield(,) value
c       layers                          number of layers to be added 
c       nobarr1,nobarr2                 remove barriers,print barrier calcs
c       mpr,iprt(),jprt()               print calcs at mpr points
c       numnulls,nblock                 check for drying:
c                                       numnulls=number of nulls in 
c                                         triangular element above which 
c                                         the marine grid point is null.
c                                       nblock: if=1, prevent layers
c                                         adjacent to dried maring grid points
c       adcirc_datum_file               1.or 'none'. has datum/9.999 +++++
c       us_eez.dat                      2.the US EEZ polygon
c       coastline_file_name             3.coastline file
c       bounding_poly_file_name         4.the BP file
c       marine.dat                      5.marine grid. used for vpop.f
c       marine.gtx                      6.marine grid in GTX format             
c       water_init.gtx                  7.all initial water points
c       ponds_init.gtx                  8.points in initial, isolated ponds
c       dry_pts.gtx                     9.initial dried points

c-------------------------------------------------------------------------

c     variables:
c       clat(j) = lat. of VDatum grid point
c       clon(i) = lon. of VDatum grid point
c         ibmax = number of BP points read                    (array size)
c          icww = index to make points in a segment 
c                 counter-clockwise : 1=yes, 0=no
c          imax = number of points in vgrid in the 
c                 x (longitude) direction                     (array size)
c          ione = one (1), was number of segments in BP (see nsegb)
c         ismax = number of coastline points read (was=icmx)  (array size)
c   jfield(i,j) = marine grid index (0=land,1=water,2,3,etc=water by layers)
c          jmax = number of points in vgrid in the 
c                 y (latitude) direction                      (array size)
c        layers = number of extra cell layers to be added
c          ncmx = max number of crossings  (was lcmx)         (array size)
c       nobarr1 = index to remove barriers between points after land/water
c                 initially set. Barrier occurs when, if a cell is defined
c                 by 4 VDatum points at corners, two adjacent cells have
c                 1 or 2 nulls in the common corners. 
c                 if nobarr1==2, then do not stop if unresolved barriers
c       nobarr2 = index to print in sub. remove_barriers
c         nsegb = max number of segments in BP                (array size)
c         nsegs = max number of segments in coastline file    (array size)
c                 (was lsmx)

c     subroutines in vgridder.f:
c        prelims
c        make_layers(lu)
c        reconcile(lu)
c        maxmin(x,y,ipen,ismax,xmins,xmaxs,ymins,ymaxs,nsegs)
c        printgrd(lu)
c        pri
c        select_water_new
c        rdcoast(ipr,lu,k)
c        rdpoly(ipr,lu,k)
c        rdcon
c        setcon
c        endlons(i1,i2,iskip,nperl,k,k1max)
c        reduce_segments
c        rname(fname)
c        printi(ia,len,ifrst,ilast,jfrst,jlast,iskip,jskip,
c        coastprep2(ipr,ismax,x,y,ipen, nsegs,jccw,ilong,isend)
c        inside_poly(ipr,np,x,y,xt,yt, in,n1,iselect)
c        crosses(ipr,j,icrx,isgx)
c        cross_segment(ipr,ns,np,x,y,yt,ncross,xcross,icrx)

c     functions:
c        degs(x,y)
c        ipri(i)
c        iprj(j)
c        iprij(i,j)
c----------------------------------------------------------------------
      program vgridder
      use saved_data
      real*8 c1,c2
      integer netida1(8),netida2(8)   ! computer time variables

      write(6,"(' PROGRAM 'a14)")pname
      call date_and_time(values=netida1)

c         Read control file 
      write(6,"(/1x72(1h-)/' STEP 1. Read in data'/1x72(1h-))")
      call rdcon       

c        Compute array sizes and constants
      call setcon ; call ztime(netida1,1)

c        get EEZ 
      write(6,"(/1x72(1h-)/' STEP 2: Read EEZ, Coast, BP'/1x72(1h-))")
      write(6,"(/1x72(1h-)/'  Read in EEZ'/1x72(1h-))")
      k=2 ; iemax=0 ; ipr=0 ; lu=10
      if(filenm(k)/='none')then
        call rd_xyi(ipr,lu,k) ; allocate (xe(imx),ye(imx),ipene(imx),
     *  ifse(nsegt)) ; nsege=nsegt ; ifse(1:nsege+1)=ifs(1:nsege+1)
        iemax=imx ; xe(1:iemax)=xx(1:imx) ; ye(1:iemax)=yy(1:imx) 
        ipene(1:iemax)=ipet(1:imx) ; deallocate (xx,yy,ipet,ifs)
      endif

c        get coastline 
      write(6,"(/1x72(1h-)/'  Read in coastline'/1x72(1h-))")
      k=3 ! coast file index: 1=shoreline
      ipr=0 ; lu=10 ; ismax=0
      if(filenm(k)/='none')then
        call rd_xyi(ipr,lu,k) ; allocate (xcs(imx),ycs(imx),ipens(imx))
        ismax=imx ; xcs(1:ismax)=xx(1:imx) ; ycs(1:ismax)=yy(1:imx) 
        ipens(1:ismax)=ipet(1:imx) ; deallocate (xx,yy,ipet,ifs)
        nsegc=nsegt ; call segments_coast
      endif

c        read bounding polygon, save in xbc,ybc,ipenb
      write(6,"(/1x72(1h-)/'  Read in bounding polygon'/1x72(1h-))")
      ibmax=0 ; ipr=0 ; k=4
      if(filenm(k)/='none')then
        call rdpoly(ipr,lu,k)
      else
        write(6,"(3x'no bounding polygon data')")
      endif
      call ztime(netida1,1)

c        read the ADCIRC file for nodes that have dried  (starting with v20)
      if(numnulls>0)then
        iallocsor=0
        if(filenm(1)/='none')then
          open(unit=lu,file=filenm(1),form='formatted')
          call rd_adcirc(lu,1) ; close(lu)
        endif
      endif                                            


c        loop thru VDatum array, select land or water
      write(6,"(/1x72(1h-)/' STEP 3: Create a Marine Grid'/1x72(1h-))")
c        select all points in Vdatum grid (jfield) that are:
c         (a) within water  AND (b) within bounding polygon.
      call select_water2
c        print results
      write(6,"(/' jfield after initial selection: 0=land, 1=water')")
c      call pril2('jfield initial',jfield,imax,jmax,lm,al4)
      call ztime(netida1,1)
	  flush(6)
      k=7 ; call wrtout_gtx(k)! initial land/water

c        change user-selected points if needed (before layers added)
      write(6,"(/1x72(1h-)/' STEP X: Adjust points in Grid'/1x72(1h-)/
     * 1x'nadjust='i2)")nadjust
      if(nadjust>0)then
        call adjust
c        call pril2('jfield, adjusted',jfield,imax,jmax,lm,al4)
      endif 
	  if(nadjust<0) then
		goto 99
      endif
	  
c        first check for ponds
      if(ielim_ponds>0)then
        write(6,"(/1x72(1h-)/' STEP 3a: Check for ponds'/1x72(1h-)/3x
     * 'ielim_ponds='i2)")ielim_ponds
        call ponds(0,8)  ! io unit, set pond to, write to filenm(iout)
c        call pril2('jfield w/o ponds',jfield,imax,jmax,lm,al4)
      endif

c        check for ADCIRC nodes that are inter-tidal. set jfield=-1.
      write(6,"(/1x72(1h-)/' STEP 4: Check for dry nodes'/1x72(1h-)/
     * 1x'numnulls='i2)")numnulls
      if(numnulls>0)then 
        call ck_nodes ;call ztime(netida1,1)
        k=9 ; call wrtout_gtx(k)! land/water with dried nodes
        call ponds(-1,0)  ! io unit, set pond to, write to filenm(iout)
c        call pril2('jfield w/o dry node ponds',jfield,imax,jmax,lm,al4)
      endif

c        add layers of cells to existing ifield
      write(6,"(/1x,72(1h-)/' STEP 6: Add layers if needed'/1x72(1h-))")
      if(layers==0)then
        write(6,"(' layers='i3)")layers
      else
        c1=cos(.5*(alat0+alat1)*raddeg) ; c2=111120.
        write(6,"(' layers='i3'  layers*delx(m)='f6.1'  layers*dely(m)='
     *   f6.1)")layers,c2*layers*delx*c1,c2*layers*dely
        nx=nint(500./(c2*delx*c1)) ; ny=nint(500./(c2*dely))
        write(6,"(8x'For dx, 500 m is N='i3' layers. N layers='f6.1)")
     *  nx,c2*nx*delx*c1
        write(6,"(8x'for dy, 500 m is N='i3' layers. N layers='f6.1)")
     *  ny,c2*ny*dely
        call make_layers ;if(layers>0)call ztime(netida1,1)
      endif

c        remove points outside the BP, (e.g., in layers)
      write(6,"(/,1x,72(1h-),/,' STEP 7: Remove points outside the',
     *' BP'/1x72(1h-))")
      if(ibmax>0)call ck_bp

c        remove points outside the EEZ
      write(6,"(/,1x,72(1h-),/,' STEP 7a: Remove points outside the',
     *' EEZ'/1x72(1h-))")
      if(iemax>0)call ck_eez


c        check for and remove isolated ponds
      write(6,"(/1x72(1h-)/' STEP 8: Check for ponds'/1x72(1h-)/3x
     * 'ielim_ponds='i2)")ielim_ponds
      call ponds(0,0) ; call ztime(netida1,1)


c        Eliminate barriers: nobarr: 1=eliminate, 0=do not eliminate
      write(6,"(/,1x,72(1h-)/' STEP 9: Eliminate barriers '
     * /1x72(1h-)/3x'nobarr1='i2)")nobarr1
      if(nobarr1>=1)call remove_barriers2
      call ztime(netida1,1)

c        print out final field
      k=6 ;
      call wrtout_gtx(k)


c        print final VDatum (corner) cells
      write(6,"(/' Field, final jfield.  0=land, 1=water, 2=layer')")
c      call pril2('jfield final',jfield,imax,jmax,lm,al4)

      write(6,"(/' End of calcs. Number of ponds found='i4)")iset
      write(6,"(//1x,70('-')//1x'End Program 'a14)")pname
      call ztime(netida1,2)

 99   end
c-------------------------------------------------------------------
      subroutine adjust
c        reset jfield(,) points after inital land/water calculation
c     jfield(,) = land/water field: 0=land, 1=water, 2,3,etc=water layer
      use saved_data
      write(6,"(/' <adjust>')")
c        change after land/water determination
      do n=1,nadjust
        write(6,"(4x'n='i6'  i,j='2i8' original, revised jfield='2i3)")
     *  n,iadj1(1:2,n), jfield(iadj1(1,n),iadj1(2,n)), iadj1(3,n)
        jfield(iadj1(1,n),iadj1(2,n))=iadj1(3,n)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine ck_bp
c     if a water/layer point is outside the BP, change to land point
c       jfield(,) = land/water field: land=0, water=1, layer=2,3,etc
c     variables:
c         npoly1 = number of first (outer) polygon
c       ibfrst() = record of first point in a segment
c          ibmax = num points in polygon file
c            ncr = number of crossings in a row
c          ncmax = max number of crossings in a row
c          nsegb = number of segments in BP
c     xcr(1:ncr) = crossing points 
c          xt,yt = test longitude, latitude
c    xbc(),ybc() = points in BP
      use saved_data
      parameter (ncmax=200) ! max number of crossings in any one segment
      real*8 xt,yt,xcr(ncmax) 
      write(6,"(/' <ck_bp>'/'  ibmax='i10' nsegb='i5)")ibmax,nsegb
      write(6,"(' ibfrst='9i8)")ibfrst(1:nsegb+1)
      nch=0 ! number changed
      do j=1,jmax ; yt=clat(j) ; xt=clon(1)-delx
c        compute bounding polygon crossings at latitude 'yt'
        k=ibfrst(npoly1) ; m=ibfrst(npoly1+1)-1 ; nl=m-k+1
        call inside_pert(nl,xt,yt,xbc(k:m),ybc(k:m),ins,ncmax,ncr,xcr)
c       if(ncr==0)cycle  !++++++++++++++++++++++++++++++ Ver19 +++++
        inner: do i=1,imax
          if(jfield(i,j)==0)cycle inner
          if(ncr==0)then;jfield(i,j)=0;nch=nch+1;cycle;endif !+++ Ver19 +++
c           count number of crossings to right (=ncr)
          isum=0 ; do n=1,ncr ; if(clon(i)<xcr(n))isum=isum+1 ;enddo
          ins=0  ; if(mod(isum,2)==1)ins=1 ; if(ins==1)cycle
          jfield(i,j)=0 ; nch=nch+1
        enddo inner
      enddo
      write(6,"('  number changed='i8)")nch
      return
      end
c----------------------------------------------------------------------
      subroutine ck_eez
c     if a water point is outside the EEZ, change to land point
c       jfield(,) = land/water field: land=0, water=1, layer=2,3,etc
      use saved_data
      parameter (ncmax=200) ! max number of crossings in any one segment
      real*8 xt,yt,xcr(ncmax),xtotalc(ncmax)
      integer nsetcr(ncmax)

      write(6,"(/' <ck_eez>'/'  iemax='i10' nsege='i5)")iemax,nsege
      write(6,"('  ifse='9i8)")ifse(1:nsege+1)

c       ntot=0 ; do i=1,nsv(1,n) ; if(xt<xcro(i,n))ntot=ntot+1 ; enddo

      nch=0 ! number changed
      do j=1,jmax ; yt=clat(j) ; xt=clon(1)-delx
c         get number of crossings for each segment
        ntot=0
        do n=1,nsege
          n1=ifse(n) ; n2=ifse(n+1)-1 ; nl=n2-n1+1
          call inside_pert(nl,xt,yt,xe(n1:n2),ye(n1:n2),i,ncmax,ncr,xcr)
          if(ncr==0)cycle
          xtotalc(ntot+1:ntot+1+ncr)=xcr(1:ncr)
          nsetcr (ntot+1:ntot+1+ncr)=n
          ntot=ntot+ncr
        enddo
        if(j<=3)then
          write(6,"(' j='i4' ntot='i5)")j,ntot
          write(6,"(' xcr='8f8.2)")xtotalc(1:ntot)
          write(6,"(' nst='8i8)")nsetcr(1:ntot)
        endif

c         loop thru columns
        do i=1,imax ; xt=clon(i) 
          if(jfield(i,j)==0)cycle
c         count number of crossings to the east
          ins=0
          do n=1,nsege ! loop thru segments
            isum=0
            do k=1,ntot
              if(xt<xtotalc(k).and.nsetcr(k)==n)isum=isum+1
            enddo
           if(mod(isum,2)==1)ins=1
          enddo 
          if(ins==1)cycle
          jfield(i,j)=0 ; nch=nch+1
        enddo 
       enddo
      write(6,"('   number changed='i8)")nch

      return
      end
c----------------------------------------------------------------------
      subroutine make_layers
c      Add a number of layers of water points around water grids,
c      based on jfield(i,j). 
c      If at least one of the surounding 
c         four points [(i+1,j), (i-1,j), (i,j+1), (i,j-1)] is water,
c         turn (i,j) into water
c      Unless a point has value of -1.

c     variables:
c      ifield0(,) = index field: 0=land, 1=switched to water
c       jfield(,) = land/water field: 0=land, 1=water, 2,3,etc=water layer

      use saved_data
      integer ifield0(imax,jmax)

      write(6,"(/' <make_layers>'/'  layers       ='i3)")layers

c        find which polygon has number '1'
      nx=npoly1   ! number of the polygon of inclusion
      write(6,"('  polygon with number 1 is 'i2/)")nx

c        loop thru layers. if point is switched to water, ifield0(i,j)=1
      do l=1,layers
        jnew=1+l ! set the number for the layer
c        tag land points as water if at least 1 adjacent pt is water.
        ifield0(1:imax,1:jmax)=0      ! initialize the tag
        do i=2,imax-1
          do j=2,jmax-1
            if(jfield(i,j)>0)cycle    ! already a layer
            if(jfield(i,j)<0)cycle    ! dry point
            idx1=max(jfield(i-1,j),jfield(i+1,j),jfield(i,j-1),
     *      jfield(i,j+1))            ! look for water
            idx2=min(jfield(i-1,j),jfield(i+1,j),jfield(i,j-1),
     *      jfield(i,j+1))            ! look for dry point
            if(idx1>0)ifield0(i,j)=1
c              blocking:  neg. layers added to dry points
            if(nblock==1.and.idx2<0)ifield0(i,j)=-1
          enddo
        enddo         ! end loop thru jfield array
c          turn tagged cells [ifield0(,)] into layer
        nflip=0 ! number of points that are reset
        do i=2,imax-1
          do j=2,jmax-1
            if(ifield0(i,j)==0)cycle 
            if(ifield0(i,j)== 1)jfield(i,j)= jnew ; nflip=nflip+1
            if(ifield0(i,j)==-1)jfield(i,j)=-jnew ;! nflip=nflip+1
          enddo
        enddo
c         print istts
        write(6,"('  layer='i5' jnew='i3' num added='i12)")l,jnew,nflip
        if(nflip==0)exit
      enddo ! end of layers

c        print result
      write(6,"(/' Field after adding'i3' layer(s)')")layers
c      call pril2('jfield with +/- layers',jfield,imax,jmax,lm,al4)

c        remove negative layer points
      write(6,"(/' Field after removing dry layers')")
      do i=1,imax ; do j=1,jmax
        if(jfield(i,j)<-1)jfield(i,j)=0
      enddo ; enddo
c      call pril2('jfield with + layers',jfield,imax,jmax,lm,al4)

      return
      end
c-------------------------------------------------------------------
      subroutine multisweep(ip,ia,imax,ja,jmax,ie,je)
c        create i and j lower and upper limits for alternating
c        direction sweeps through a 2-d array.

c        will enable [for k=1-4, k=1+mod(it-1,4) ]
c          do i=i1(k),i2(k),sign(1,i2(k)-i1(k))
c          do j=j1(k),j2(k),sign(1,j2(k)-j1(k))
c        NOTE: for iteration=it, m=1+mod(it-1,4)

c     variables:
c       i1(m),i2(m) = first,last i-index for sweep 'm', m=1:4
c       j1(m),j2(m) = first j-index for sweep 'm'
c           ia,imax = first,last i-index in field
c           ja,jmax = first,last j-index in field
c                ip = print index
      integer i1(4),i2(4),j1(4),j2(4)
      integer ie(4,2),je(4,2)

      i1(1:4)=(/ia,imax,ia,imax/) ; i2(1:4)=(/imax,ia,imax,ia/)
      j1(1:4)=(/ja,jmax,jmax,ja/) ; j2(1:4)=(/jmax,ja,ja,jmax/)
      ie(1:4,1)=(/ia,imax,ia,imax/) ; ie(1:4,2)=(/imax,ia,imax,ia/)
      je(1:4,1)=(/ja,jmax,jmax,ja/) ; je(1:4,2)=(/jmax,ja,ja,jmax/)

      if(ip/=1)return
      write(6,"(3x'sweep:1=+x,+y 2=-x,-y 3=+x,-y 4=-x,+y')")
      write(6,"(/' i1(1-4)='4i10/' i2(1-4)='4i10)")i1(1:4),i2(1:4)
      write(6,"(/' j1(1-4)='4i10/' j2(1-4)='4i10)")j1(1:4),j2(1:4)
      write(6,"(' <multisweep>'4x'm='i2'  i1,i2='2i5'  j1,j2='
     *  2i5)")(m,i1(m),i2(m),j1(m),j2(m),m=1,4)

      return
      end
c-------------------------------------------------------------------------
      subroutine barrprt(idx,imax,jmax,isv,jfield,i,j)
c        print fields for barriers
c     jfield(,) = land/water field: 0=land, 1=water, 2,3,etc=water layer
      integer ii(3),isv(imax,jmax),jfield(imax,jmax)
      i1=i-1 ; i2=i+1 ; ii(1:3)=(/i1,i,i2/)
      if(idx==1)write(6,"(/'  vert barrier @ i,j='2i6' j+1='i6)")i,j,j+1
      if(idx==2)write(6,"(/'  horz barrier @ i,j='2i6' i+1='i6)")i,j,i+1
      write(6,"(16x'i='3i6,10x,3i6)")ii(1:3),ii(1:3)
      write(6,"(5x'j='i5'  isv='3i6,3x'jfield='3i6)")(jj,isv(i1:i2,jj),
     *  jfield(i1:i2,jj),jj=j+1,j-1,-1)
      return
      end
c-------------------------------------------------------------------------
      subroutine remove_barriers2

c        to add a                 to add a
c      vertical barrier(:)     horizontal barrier (..)

c       O-------O-------X       X-------X    j+1
c       |       :       |       |       |
c       |       :       |       |       |
c       |       :       |       |       |
c       X-------O-------O    i,jO.......O     j
c             i,j               |       |
c                               |       |
c       X = filled              |       |
c       O = null                X-------O    j-1


c     variables:
      use saved_data
      logical p,q

      write(6,"(/'<remove_barriers2>'/' layers='i3)")layers
      p=(nobarr2>=1)  ! print index
      allocate (isv(imax,jmax),STAT=istts)  ! point index
	  IF (istts /= 0) THEN
		write(6,'("Allocation 1 failed!")')
	  ELSE
		write(6,'("Allocation 1 successful")')
	  END IF
	  flush(6)
c        initialize
      meth=1    ! revised jfield equals layers+2
      meth=2    ! revised jfield equals max(jfield)+1
 
c        new point index: 0=land, 1=water or layer
      do i=1,imax ;do j=1,jmax 
        isv(i,j)=min(1,jfield(i,j))
      enddo ; enddo

c        begin iteration
      iter=0   ! iteration number
      ntotal=0 ! total points that have changed
      do
        iter=iter+1
        numc=0       ! number changed this iteration

c    1. look for potential vertical barrier (denoted by :)
c       j+1  X-------O-------X
c            |       :       |
c            | cell1 : cell2 |
c            |       :       |
c       j    X-------O-------X
c            i      i+1     i+2
        do i=2,imax-2 ;do j=2,jmax-2  
c           compute indices
          isum_left=isv(i  ,j)+isv(i  ,j+1)
          isum_cent=isv(i+1,j)+isv(i+1,j+1)
          isum_rght=isv(i+2,j)+isv(i+2,j+1)
          isum_bott=isv(i,j  )+isv(i+2,j)
          isum_topp=isv(i,j+1)+isv(i+2,j+1)
          icell1=isum_left+isum_cent
          icell2=isum_rght+isum_cent
c           make checks here
          if(icell1==0)cycle    ! cell1 has no valid corners
          if(icell2==0)cycle    ! cell2 has no valid corners
          if(isum_cent/=0)cycle ! barrier already exists
c           print initial situation
          numc=numc+1
          if(p)then
           write(6,"(/' iter='i3' numc='i3' potential missing ',
     *     'HORIZONTAL barrier at i='i4' j='i4)")iter,numc,i,j
           write(6,"(10x'initial:')")
           write(6,"('  j='i5' isv='3i3,5x,'jfield='3i3)")(jj,
     *     isv(i:i+2,jj),jfield(i:i+2,jj),jj=j+1,j,-1)
          endif
c           get new jfield value
          if(meth==1)then
             laynew=layers+2
          else
            laymax=0
            do ii=i,i+2 ; do jj=j,j+1
              laymax=max(laymax,jfield(ii,jj))
            enddo ; enddo
            laynew=laymax+1
          endif
c           switch isv and jfield
          do k=0,1
            if(k==0)ns=isum_bott
            if(k==1)ns=isum_topp
            if(ns>0.and.jfield(i+1,j+k)==0)then
              jfield(i+1,j+k)=laynew
              isv(i+1,j+k)=1
              ntotal=ntotal+1
            endif
          enddo
c           print final situation
          if(p)then
            write(6,"(10x'revised: laymax='i3)")laymax
            write(6,"('  j='i5' isv='3i3,5x,'jfield='3i3)")(jj,
     *      isv(i:i+2,jj),jfield(i:i+2,jj),jj=j+1,j,-1)
          endif
        enddo ; enddo


c    2. look for potential horizontal barrier (denoted by ..)
c       j+2  X-------X
c            |       |
c            | cell2 |
c            |       |
c       j+1  0.......O
c            |       |
c            | cell1 |
c            |       |
c       j    X-------X 
c            i      i+1

        do i=2,imax-2 ;do j=2,jmax-2  
c           compute indices
          isum_bott=isv(i,j  )+isv(i+1,j)
          isum_cent=isv(i,j+1)+isv(i+1,j+1)
          isum_topp=isv(i,j+2)+isv(i+1,j+2)
          isum_left=isv(i,j  )+isv(i  ,j+2)
          isum_rght=isv(i+1,j)+isv(i+1,j+2)
          icell1=isum_bott+isum_cent
          icell2=isum_topp+isum_cent
c           make checks here
          if(icell1==0)cycle    ! cell1 has no valid corners
          if(icell2==0)cycle    ! cell2 has no valid corners
          if(isum_cent/=0)cycle ! barrier already exists
c           print initial situation
          numc=numc+1
          if(p)then
            write(6,"(/' iter='i3' numc='i3' potential missing VERTICA',
     *     'L barrier at i='i4' j='i4)")iter,numc,i,j
            write(6,"(10x'initial:')")
            write(6,"('  j='i5' isv='2i3,5x,'jfield='2i3)")(jj,
     *      isv(i:i+1,jj),jfield(i:i+1,jj),jj=j+2,j,-1)
          endif
c           new layer value
          if(meth==1)then
            laynew=layers+2
          else
            laymax=0
            do ii=i,i+1 ; do jj=j,j+2
              laymax=max(laymax,jfield(ii,jj))
            enddo ; enddo
           laynew=laymax+1
          endif
c           switch isv and jfield
          do k=0,1
            if(k==0)ns=isum_left
            if(k==1)ns=isum_rght
            if(ns>0.and.jfield(i+k,j+1)==0)then
              jfield(i+k,j+1)=laynew
              ntotal=ntotal+1
              isv(i+k,j+1)=1
            endif
          enddo
c           print final situation
          if(p)then
            write(6,"(10x'revised: laynew='i3)")laynew
            write(6,"('  j='i5' isv='i3,5x,'jfield='i3)")(jj,
     *      isv(i:i+1,jj),jfield(i:i+1,jj),jj=j+2,j,-1)
          endif
        enddo ; enddo

        if(numc==0)exit

      enddo ! on iter

      write(6,"(/' iter='i3' ntotal changed='i3)")iter,ntotal
      return   
      end
c-------------------------------------------------------------------------
      subroutine ponds(jset,iout) ! io unit, set pond to, iout=1=write
c        set all water points to ipond(i,j)=0
c        set 'iset' = 1.
c        find first water point. 
c        set 'ipond' for all points connected to this point to 'iset'.
c        count the number of points with this value of 'iset'.
c        find another water point that has 'ipond' < to 'iset'
c     variables:
c   ielim_ponds = index to remove ponds (isolated water areas): >=1=yes
c          iout = if iout>0, w/o to filenm(iout)
c            is = set with largest number of points
c       isetmax = number of ponds
c     ijpon(,n) = i,j of points to set to land
c     jfield(,) = land/water field:-1=dry,0=land,1=water,2,3,etc=layer
c          jset = value to reset jfield to
c          nc() = number of points in each set
c      ipond(,) = index: 0=land, 1=in set 1, 2=in set 2, etc.
      use saved_data
      USE NETCDF        ! enables the fortran 90 interface to netcdf
	  
		   
      parameter (isetmax=100000)
      real*8 b,c,d,xt(2),yt(2),z1,d1(2),d2(2)
      integer ipond(imax,jmax),ip(4),jp(4),nc(isetmax),ie(4,2),je(4,2),
     * iex(4),isave(4,isetmax)
      character aset*2,ttl*10
      logical p
      data ip,jp/-1,1,0,0, 0,0,-1,1/,p/.false./,b/360./
	  integer ncid,xDimId,yDimId,bxId,zDims(2),bid,zid,ncerr,n
	  real*8 box(4)
	  character(len=100) :: fname
	  integer,allocatable,dimension(:,:) :: jz
	  
      fm(d)=60.*abs(d-int(d)) ! function for minutes

      write(6,"(/' <ponds>'/' ielim_ponds='i6)")ielim_ponds

c        initialize pond index
      ipond(1:imax,1:jmax)=0
      p=.true. ;  p=.false. ! print the individual ponds 

c        count wet points
      ntot=0
      do i=1,imax
        do j=1,jmax ; if(jfield(i,j)>0)ntot=ntot+1 ;enddo
      enddo ; write(6,"(4x'total water='i8)")ntot
 
c        get multi-sweep limits
      call multisweep(0,2,imax-1,2,jmax-1,ie,je) ! i1,i2,j1,j2)

c        loop thru jfield grid
      write(6,"(4x'look for isolated water')")
      iset=1         ! set to be looking for
      itmax=imax*jmax   
      do i=2,imax-1
        do j=2,jmax-1
c           check only water points
          if(jfield(i,j)<=0)cycle   
c           skip if already in a group
          if(ipond(i,j)>0.and.ipond(i,j)<=iset)cycle
c           set all adjacent to iset
          ipond(i,j)=iset
          if(iset==1)write(6,"(3x'iset='i5' first at i,j='2i8)")iset,i,j
          it=1 ;  nclast=0 !; iex(1:4)=(/imax+1,0,jmax+1,0/)
          do  ! on it
            m=1+mod(it-1,4)
            do ii=ie(m,1),ie(m,2),sign(1,ie(m,2)-ie(m,1))
              do jj=je(m,1),je(m,2),sign(1,je(m,2)-je(m,1))
                if(jfield(ii,jj)<=0.or.ipond(ii,jj)>0)cycle  ! +++++
c                  check adjacent points
                do k=1,4
                  if(jfield(ii+ip(k),jj+jp(k))<=0)cycle
                  if( ipond(ii+ip(k),jj+jp(k))==iset)then
                    ipond(ii,jj)=iset ; exit
                  endif
                enddo
              enddo
            enddo
c             count points in this pond
            nc(iset)=0
            do ii=2,imax-1 ; do jj=2,jmax-1
              if(ipond(ii,jj)==iset)nc(iset)=nc(iset)+1
            enddo ; enddo
            if(nclast==nc(iset))exit ; nclast=nc(iset)
            if(it==itmax)exit ; it=it+1
          enddo ! it
c           get pond area limits
          iex(1:4)=(/imax+1,0,jmax+1,0/)
          do ii=2,imax-1 ; do jj=2,jmax-1
            if(ipond(ii,jj)/=iset)cycle
            iex(1:2)=(/min(iex(1),ii-1),max(iex(2),ii+1)/)
            iex(3:4)=(/min(iex(3),jj-1),max(iex(4),jj+1)/)
          enddo ; enddo
c          write(6,"(3x' iset='i6' nc='i9' it='i9' limits='4i6)")iset,
c     *    nc(iset),it,iex(1:4)
          isave(1:4,iset)=iex(1:4) ! save i,j limits
c            update iset
          iset=iset+1
          if(iset>isetmax)then
            write(6,"(/' *** iset='i4' *** too large ***')")iset ; stop
          endif
        enddo ! on i
      enddo   ! on j
      iset=iset-1

c        summary
      write(6,"(/' Sumary of ponds found. number='i7)")iset
      if(iset>0)then
c         find set with most points (i.e., the main water)
        isum=0 ; is=1
        do i=1,iset
          isum=isum+nc(i) ; if(nc(i)>nc(is))is=i
        enddo
        write(6,"(4x'total of nc='i8'  set with most='i3/)")isum,is
      endif

c        write lon, lat limits for plotting
      if(iset>0.and.p)then
        c=3. ;d1(1:2)=(/-c*delx,c*delx/) ;d2(1:2)=(/-c*dely,c*dely/) 
        write(6,"(/4x'lon, lat limits for ponds (w/o main). c='f5.2)")c
        ns=0
        do i=1,iset
          if(i==is)cycle ; ns=ns+1
          xt(1:2)=clon(1)+delx*(isave(1:2,i)-1)+d1(1:2)
          yt(1:2)=clat(1)+dely*(isave(3:4,i)-1)+d2(1:2)
          write(6,"(2(i5,f6.2),2(i4,f6.2)'  pond #'i7)")int(xt(1)),
     *    fm(xt(1)),int(xt(2)),fm(xt(2)),int(yt(1)),fm(yt(1))
     *    ,int(yt(2)),fm(yt(2)),ns
        enddo
c         print out array
        do i=1,iset
          if(i<10)aset=' '//char(ichar('0')+i)
          if(i>9 )aset=char(ichar('0')+i/10)//char(ichar('0')+mod(i,10))
          if(i>99)then
            write(6,"(/' *** no. ponds too large. i='i3' ***')")i ; stop
          endif
          ttl='ipond # '//aset
          if(i==is)then
            write(6,"(/13x,a10' is the main body of water')")ttl
          else
            lm2(1:6)=(/isave(1,i),isave(2,i),1,isave(3,i),isave(4,i),1/)
c            call pril2(ttl,ipond,imax,jmax,lm2,al4)
c            call pril2('pond jfield',jfield,imax,jmax,lm2,al4)
            xt(1:2)=clon(1)+delx*(isave(1:2,i)-1)
            yt(1:2)=clat(1)+dely*(isave(3:4,i)-1)
            write(6,"(12x'xmin, xmax='2(i4,f6.2))")int(xt(1)),fm(xt(1)),
     *      int(xt(2)),fm(xt(2))
            write(6,"(12x'ymin, ymax='2(i4,f6.2))")int(yt(1)),fm(yt(1)),
     *      int(yt(2)),fm(yt(2))
          endif
        enddo ! on iset
      endif ! on p (writing arrays)

c       write out ponds (w/o largest pond) as a gtx file
      if(iout>0)then
		fname=filenm(iout)
	    if(fname/='none')then
			n=len(trim(fname))
			fname((n+1):(n+3))='.nc'	  
			box(1)=clat(1)
			box(2)=clon(1)+360.
			box(3)=dely
			box(4)=delx
			allocate (jz(imax,jmax),STAT=istts)
			IF (istts /= 0) THEN
				write(6,'("Allocation 2 failed!")')
			ELSE
				write(6,'("Allocation 2 successful")')
			END IF

			jz(1:imax,1:jmax)=-88
			
			where( ipond>0.and.ipond/=is ) jz=ipond
		  
			write(6,"(/' <wrtout_gtx>         k='i2/2xa78)")k,trim(fname)
			! Create NetCDF file
			ncerr = nf90_create(trim(fname),nf90_CLOBBER+NF90_HDF5,ncid)
			call errhandle(ncerr)
			! Define dimensions
			call errhandle(nf90_def_dim(ncid, 'box', 4, bxId))
			call errhandle(nf90_def_dim(ncid, 'longitude', imax, xDimId))
			call errhandle(nf90_def_dim(ncid, 'latitude', jmax, yDimId))
			zDims(1)=xDimId
			zDims(2)=yDimId				
			! Define variables
			ncerr=nf90_def_var(ncid,'box',nf90_DOUBLE,bxId,bid)
			call errhandle(ncerr)	
			ncerr=nf90_def_var(ncid,'field',nf90_INT,zDims,zid)
			call errhandle(ncerr)		
			ncerr=nf90_def_var_deflate(ncid,zid,0,1,8)
			call errhandle(ncerr)		
			! End definitions
			call errhandle(nf90_enddef(ncid))
			! Write data to variables
			call errhandle(nf90_put_var(ncid, bid, box))
			call errhandle(nf90_put_var(ncid, zid, jz))
			! Close NetCDF file
			call errhandle(nf90_close(ncid))

			print *, "NetCDF file created successfully"
			deallocate(jz)
        endif
      endif

      if(p.and.iset>1)then
c       write out water points (w/o layers) inside ponds as a gtx file
		fname=filenm(8)
	    if(fname/='none')then
			n=len(trim(fname))
			fname((n+1):(n+3))='.nc'	
			box(1)=clat(1)
			box(2)=clon(1)+360.
			box(3)=dely
			box(4)=delx
			allocate (jz(imax,jmax),STAT=istts)
			IF (istts /= 0) THEN
				write(6,'("Allocation 3 failed!")')
			ELSE
				write(6,'("Allocation 3 successful")')
			END IF
			flush(6)

			jz(1:imax,1:jmax)=0
			
			where( ipond>0.and.ipond/=is.and.jfield==1 ) jz=1
		  
			write(6,"(/' <wrtout_gtx>         k='i2/2xa78)")k,trim(fname)
			! Create NetCDF file
			ncerr = nf90_create(trim(fname),nf90_CLOBBER+NF90_HDF5, ncid)
			call errhandle(ncerr)
			! Define dimensions
			call errhandle(nf90_def_dim(ncid, 'box', 4, bxId))
			call errhandle(nf90_def_dim(ncid, 'longitude', imax, xDimId))
			call errhandle(nf90_def_dim(ncid, 'latitude', jmax, yDimId))
			zDims(1)=xDimId
			zDims(2)=yDimId				
			! Define variables
			ncerr=nf90_def_var(ncid,'box',nf90_DOUBLE,bxId,bid)
			call errhandle(ncerr)	
			ncerr=nf90_def_var(ncid,'field',nf90_INT,zDims,zid)
			call errhandle(ncerr)		
			ncerr=nf90_def_var_deflate(ncid,zid,0,1,8)
			call errhandle(ncerr)		
			! End definitions
			call errhandle(nf90_enddef(ncid))
			! Write data to variables
			call errhandle(nf90_put_var(ncid, bid, box))
			call errhandle(nf90_put_var(ncid, zid, jz))
			! Close NetCDF file
			call errhandle(nf90_close(ncid))

			print *, "NetCDF file created successfully"
			deallocate(jz)
		endif
      endif ! on p

c        remove ponded points
      if(ielim_ponds>0)then
        write(6,"(/' remove all ponds')")
        do i=1,imax ; do j=1,jmax
          if(ipond(i,j)>0.and.ipond(i,j)/=is)jfield(i,j)=jset
        enddo ; enddo
      endif

      return
      end
c-------------------------------------------------------------------------
      subroutine rem_prnt(idx,imax,jmax,jfield,i,j,no)
c        print the jfield around point i,j
c        idx = index:1=horiz,2=vertical
      integer jfield(imax,jmax)
      if(idx==1)then ! horiz barrier
        i1=i-1 ; i2=i+2 ; j1=j+1 ; j2=j-1
        write(6,"(/3x'horz barrier #'i5)")no
      else
        i1=i-1 ; i2=i+1 ; j1=j+2 ; j2=j-1
        write(6,"(/3x'vert barrier #'i5)")no
      endif
      write(6,"(20x'i='4i6)")(ii,ii=i1,i2)
      do k=j1,j2,-1
        write(6,"(6x'j='i6' jfield='4i6)")k,jfield(i1:i2,k)
      enddo
      return
      end
c-------------------------------------------------------------------------
      subroutine maxmin(x,y,ipen,ismax,xmins,xmaxs,ymins,ymaxs,nsegs)
c        find max and min of lats and longs in each coastline segment
c        (used for rough screening of points inside land)
c     variables:
c              ismax = number of points in the file
c              nsegs = number of segments in file
c      xmin(),xmax() = min, max of x for each segment
      real*8 x(ismax),y(ismax),xmins(nsegs),xmaxs(nsegs),ymins(nsegs),
     * ymaxs(nsegs),xrgt,xlft,ytop,ybot
      dimension ipen(ismax)

      write(6,"(' <maxmin>')")
      isum=0 ; iseg=0
      do i=1,ismax
        if(isum==0)then
          xrgt=-400. ;  xlft= 400. ; ytop=-400. ; ybot= 400.
        endif
        xrgt=max(xrgt,x(i)) ; xlft=min(xlft,x(i))
        ytop=max(ytop,y(i)) ; ybot=min(ybot,y(i))
        isum=isum+ipen(i)
        if(isum==2)then
          isum=0 ; iseg=iseg+1
          xmaxs(iseg)=xrgt ; xmins(iseg)=xlft
          ymaxs(iseg)=ytop ; ymins(iseg)=ybot
        endif
      enddo
      write(6,"('  isegs,nseg=',2i8)")iseg,nsegs
      do n=1,min(nsegs,10)
        write(6,"('  n='i6' xmin,max='2f12.6' ymin,max='2f12.6)")
     *  n,xmins(n),xmaxs(n),ymins(n),ymaxs(n)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine cross(p,indx,yt,ncmax,nsmax,nctotal,nsv,xcros) 
c     purpose: to find the crossings in each coastline segment, for a
c              single latitude 'yt'
c     variables:
c      ibfrst(,) = first point in BP segment
c      isfrst(,) = first point in coastline segment
c           indx = index: 1=BP, 2=coastline, 3=EEZ
c            nsf = number of segments found
c     nsv(1:2,n) = number of crossings, number of segment for n-th
c  xcros(1:nx,n) = crossing points for nth segment. nx=nsv(1,n)
c             yt = test latitude
c          ibmax = num points in polygon file
c          ismax = num points in coastline file
c          ncmax = max number of crossings in a row
c        nctotal = number of segments w/crossings
c          nsmax = max number of segments with crossings
c    xbc(),ybc() = points in BP
c    xcs(),ycs() = points in coastline file
c      xe(),ye() = points in EEZ file
      use saved_data
      real*8 xt,yt,xcros(ncmax,nsmax),xcr(ncmax)
      integer nsv(2,nsmax)
      logical p

c        select starting longitude
      xt=clon(1)-delx ; if(p)write(6,"(/' <cross>  xt,yt='2f12.6
     * ' indx(BP, coast, EEZ)='i1)")xt,yt,indx
c        initialize
      nsv(1:2,1:nsmax)=0 ;xcros(1:ncmax,1:nsmax)=-180.

c        crossings of bounding polygon
      nctotal=0 ! number of crossings
      if(indx==1)then
        if(p)write(6,"('  BPoly: nsegb='i6)")nsegb
        do ns=1,nsegb  ! loop thru segments
          k=ibfrst(ns) ; m=ibfrst(ns+1)-1 ; nl=m-k+1
          call inside_pert(nl,xt,yt,xbc(k:m),ybc(k:m),ins,ncmax,ncr,xcr)
          if(ncr==0)cycle
          nctotal=nctotal+1
          nsv(1:2,nctotal)=(/ncr,ns/) ; xcros(1:ncr,nctotal)=xcr(1:ncr)
          if(p)write(6,"(3x'n='i4' nseg='i6' ncr='i5)")nctotal,ns,ncr
          if(p)write(6,"(24x'xcr='f12.5)")xcr(1:ncr)
          if(p)write(6,"(24x'k,m='2i7' first,last pt in seg')")k,m
        enddo
        if(p.and.nctotal==0)write(6,"(3x'no crossings found')")

      else if(indx==2)then
        if(p)write(6,"('  Coast: nsegc='i6)")nsegc
        do ns=1,nsegc
          k=isfrst(ns) ; m=isfrst(ns+1)-1 ; nl=m-k+1
          call inside_pert(nl,xt,yt,xcs(k:m),ycs(k:m),ins,ncmax,ncr,xcr)
          if(ncr==0)cycle
          nctotal=nctotal+1
          nsv(1:2,nctotal)=(/ncr,ns/) ; xcros(1:ncr,nctotal)=xcr(1:ncr)
          if(p)write(6,"(3x'n='i4' nseg='i6' ncr='i5)")nctotal,ns,ncr
          if(p)write(6,"(24x'xcr='f12.5)")xcr(1:ncr)
          if(p)write(6,"(24x'k,m='2i7' first,last pt in seg')")k,m
          if(p)write(6,"(24x'first pt in seg: x,y='2f12.6)")xcs(k),
     *    ycs(k)
        enddo
        if(p.and.nctotal==0)write(6,"(3x'no crossings found')")
      else if(indx==3)then
      endif

c        check for array size limits
      if(nctotal>nsmax)then
        write(6,"(//' *** nctotal='i5' > nsmax='i5' ***')")nctotal,nsmax
        write(6,"(' *** at yt='f12.6)")yt ; stop
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine countup(p,indx,xt,ncmax,nsmax,nss,nsv,xcro,inn)
c     Determine whether xt,yt is inside by counting the crossings
c       inn = returned value. 1=inside BP and outside any masks, 
c             or 1=inside coast (i.e., is land)
c       indx = index: 1=BP, 2=coastline
c      ncmax = max number of crossings in any one segment
c      nsmax = max number of segments crossed
c        nss = number of segments w/crossings
c     nsv(,) = number of crossings, segment number, for set of crossings
c       ntot = number of crossings in a segment
      use saved_data
      real*8 xt,xcro(ncmax,nsmax)
      integer nsv(2,nsmax)
      logical p
      if(p)write(6,"(/' <countup>   ncmax,nsmax='2i9)")ncmax,nsmax
c        initialize the counts
      inn1=0  ! (1=yes) inside seg # 1
      inn2=0  ! (1=yes) inside seg # 2, 3, 4, etc
c        loop thru segments
      do n=1,nss
        ntot=0 ; do i=1,nsv(1,n) ; if(xt<xcro(i,n))ntot=ntot+1 ; enddo
        if(nsv(2,n)==1.and.mod(ntot,2)==1)inn1=1
        if(nsv(2,n)> 1.and.mod(ntot,2)==1)inn2=1
      enddo
      if(p)write(6,"('  inn1='i2'  inn2='i2)")inn1,inn2
c        determine inside here
      inn=0
      if(indx==1.and. inn1==1.and.inn2==0)inn=1    ! 1=inside BP
      if(indx==2.and.(inn1==1.or.inn2==1))inn=1    ! 1=land
      return
      end
c-------------------------------------------------------------------
      subroutine select_water2
c        Determine if point represents land or water.
c        Selection based on center points. Check for:
c          over water (not island). 
c        If so, switch jfield to water.

c     variables:
c          ibmax = num points in bounding polygon file
c          ismax = num points in coastline file
c      jfield(,) = land/water field: 0=land, 1=water, 2,3,etc=water layer
c          ncb/c = number of polygon/coast sets with crossings
c       nsvb/c() = saved number of polygon/coast crossings in a set
c    xbc(),ybc() = points in BP
c    xcs(),ycs() = points in coastline file
c      xcrosb(,) = for a set, the actual crossing longitudes

      use saved_data
      parameter (ncmax=200) ! max number of crossings in any one segment
      parameter (nsmax=80)  ! max number of segments crossed
      real*8 xt,yt,xcros(ncmax,nsmax),xcrosc(ncmax,nsmax),
     * xcrosb(ncmax,nsmax)
      integer nsv(2,nsmax),nsvc(2,nsmax),nsvb(2,nsmax)
      logical p,q

      write(6,"(/' <select_water2>')")
	  flush(6)
c        initialize the VDatum grid points to null (land)
      jfield(1:imax,1:jmax)=0      

c        loop thru VDatum points
      write(6,"(/' LOOP THRU THE VDATUM GRID'/'  imax,jmax='2i6/
     *'  number of cells to print at  ='i3)")imax,jmax,mpr 
      if(mpr>0)write(6,"(9x'm='i3' i,j='2i5)")(m,iprt(m),jprt(m),
     * m=1,mpr)
	  flush(6)

      do j=1,jmax 
        yt=clat(j) ; p=.false. ; if(iprj(j)==1)p=.true. !; p=.true.
        if(p)then
			write(6,"(//' MARINE GRID ROW: j='i7' yt='f11.6)")j,yt
			flush(6)
		endif

c         compute bounding polygon crossings
        if(ibmax>0)then
          call cross(p,1,yt,ncmax,nsmax,nctot,nsv,xcros) 
          if(nctot==0)cycle      ! if none in BP, can skip entire row
          ncb=nctot ; nsvb(1:2,1:nsmax)=nsv(1:2,1:nsmax)
          xcrosb(1:ncmax,1:nsmax)=xcros(1:ncmax,1:nsmax)
        endif
c         compute coastline crossings
        idx=2
        call cross(p,idx,yt,ncmax,nsmax,nctot,nsv,xcros) ! coastal crossings
        ncc=nctot ; nsvc(1:2,1:nsmax)=nsv(1:2,1:nsmax)
        xcrosc(1:ncmax,1:nsmax)=xcros(1:ncmax,1:nsmax)

c   loop across rows
        do i=1,imax
          q=.false.  ;if(p.and.ipri(i)==1)q=.true.
          xt=clon(i) ;if(q)write(6,"(/1x'COLUMN: i='i6' xt='f12.6)")i,xt
c   check if inside BP
          if(ibmax>0)then
            call countup(q,1,xt,ncmax,nsmax,ncb,nsvb,xcrosb,inn)
            if(q)write(6,"('  inn='i2' (for BP)')")inn
            if(inn==0)then
              if(q)write(6,"(' inn=0, so outside BP and skip')"); cycle
            endif
          endif
c   check if onland
          call countup(q,2,xt,ncmax,nsmax,ncc,nsvc,xcrosc,inn)
          if(q)write(6,"('  inn='i2' (for coastline)')")inn
          if(inn==1)then
            if(q)write(6,"('  keep as land')"); cycle
          endif
c   set to water
          jfield(i,j)=1 ; if(q)write(6,"('  set to water')")
        enddo
        if(p)write(6,"(/' jfield(j='i5')')")j
        if(p)write(6,"(3x70i1)")jfield(1:imax,j)

      enddo
      write(6,"(//' END OF MARINE GRID LOOP  ----------')")
c        istts report: count number of water points
      nwater=0
      do i=1,imax ; do j=1,jmax
        nwater=nwater+jfield(i,j)
      enddo ; enddo
      write(6,"( '  no. of water points ='i8)")nwater
      if(nwater==0)then
        write(6,"(//'  *** STOP: no water points found ***')")
        stop
      endif
	  flush(6)
	  
      return
      end
c-------------------------------------------------------------------
      subroutine rd_xyi(ipr,lu,k)
C      purpose - To read and select x,y data within window.
c                Move longest segment to first.
c                Find start of each segment.
c     variables -
c            imx = number of points read
c              k = index: 1=EEZ, 2=shoreline, 3=BP
c          nsegt = number of segments read in
c         islong = length of longest segment
c      xx(),yy() = points in coastline file
c         ipet() = pen commands
      use saved_data
      real*8 ala,alo

      write(6,"(/' <rd_xyi>      lu='i3' k='i3)")lu,k
      write(6,"(' file='/1xa100)")filenm(k) ;if(filenm(k)=='none')return

c        read file to find length and number of segments
      open(unit=lu,file=filenm(k),form='formatted')
      iswitch=0   ! index to switch ala and alo
      n=0         ! number of points
      isum=0      ! sum of pen up/down
      imx=0       ! number of points
      do
        read(lu,*,iostat=io)ala,alo,ipo
        if(io/=0)exit ; isum=isum+ipo ; n=n+1 ; if(n>1)cycle
        if(ala<0.0.and.alo>0.0)iswitch=1
      enddo
      imx=n ; nsegt=isum/2 ; ns=isum/2
      write(6,"(1x'imx='i9' segments='i8)")imx,nsegt
      if(imx==0)then
        write(6,"(//'     ** no points found. run stopped **')")
        stop 
      endif
      if(iswitch==1)write(6,"(1x'NOTE: interchange ala, alo')")
      allocate (xx(imx),yy(imx),ipet(imx),ifs(nsegt+1))

c        read data again, to save and find longest seg
      rewind(lu)
      ns=0        ! number of points in a segment
      nsg=0       ! number of present segment
      isum=0      ! sum of ipo
      islong=0    ! length of longest segment
      nslong=0    ! number of longest segment
      do n=1,imx
        read(lu,*)ala,alo,ipo
        ns=ns+1
c          save coastline
        xx(n)=alo ; yy(n)=ala
        if(iswitch==1)then ; xx(n)=ala ; yy(n)=alo ; endif
        ipet(n)=ipo
        if(n==1)write(6,"(2x'n='i7' xx,yy,ipet=',
     *  2f11.5,i2)")n,xx(n),yy(n),ipet(n)
        isum=isum+ipo
        if(isum==2)then
          isum=0 ; nsg=nsg+1
          if(ns>islong)then ;islong=max(islong,ns) ; nslong=nsg ;endif
          ns=0
        endif
      enddo
      close (lu)
      write(6,"(1x'islong='i7/1x'nslong='i7)")islong,nslong
      write(6,"(1x'nsg   ='i7)")nsg
      write(6,"(1x'xx(1),yx(1)='2f12.6)")xx(1),yy(1)

c           n=n+1 ; xx(n)=xcs(i) ; yy(n)=ycs(i) ; ippx(n)=ipens(i)
      if(nslong/=1)then
        write(6,"(/' move ihe longest string to first')")
        allocate (xi(imx),yi(imx),ipex(imx))
        n=0 ; isum=0
        do i=1,imx   ! save only the longest
          isum=isum+ipet(i)
          if(isum==2*nslong-1.or.isum==2*nslong)then
            n=n+1 ; xi(n)=xx(i) ; yi(n)=yy(i) ; ipex(n)=ipet(i)
          endif
        enddo
        isum=0
        do i=1,imx
          isum=isum+ipet(i)
          if(.not.(isum==2*nslong-1.or.isum==2*nslong))then
            n=n+1 ; xi(n)=xx(i) ; yi(n)=yy(i) ; ipex(n)=ipet(i)
          endif
        enddo
        xx(1:imx)=xi(1:imx);yy(1:imx)=yi(1:imx) ;ipet(1:imx)=ipex(1:imx)
        deallocate (xi,yi,ipex)
      endif

c        get first record in each segment ifs()
      isum=0 ; nsg=0
      do n=1,imx
        if(isum==0)then
          nsg=nsg+1 ; ifs(nsg)=n
        endif
        isum=isum+ipet(n)
        if(isum==2)isum=0
      enddo
      ifs(nsg+1)=imx+1
      write(6,"(1x'nsg   ='i7)")nsg !write(6,"(1x'ifs='8i8)")ifs(1:ns+1)

      return
      end 
c-------------------------------------------------------------------
      subroutine rdpoly(ipr,lu,k)
C        purpose - To read and save bounding polygon (BP). Allow
c                  multiple segments to be read in. Additional segments
c                  will EXCLUDE points from being water.
c     variables:
c           ibmax = number of points read
c         ipenb() = pen up index
c               k = index for file name
c           nsegb = number of segments in BP file
c         npoly() = number of polygon
c          npoly1 = number of first polygon
c        ibfrst() = record of first in a segment
c     xbc(),ybc() = points in BP

      use saved_data
      real*8 ala,alo,d
      integer newpoly(nsegb)
      character ffile*80

      write(6,"(/' <rdpoly>      lu='i3)")lu
      ibmax=0 ; nsegb=0
      write(6,"(' file='/1xa100)")filenm(k) ;if(filenm(k)=='none')return
      ffile=filenm(k) ; open(unit=lu,file=filenm(k),form='formatted')

      n=0         ! number of points
      isum=0
      do
        read(lu,*,iostat=io)ala,alo,ipo
        if(io/=0)exit
        n=n+1 ; isum=isum+ipo
      enddo
      ibmax=n ; nsegb=isum/2
      write(6,"(1x'ibmax='i7' nsegb='i4)")ibmax,nsegb

      allocate (xbc(ibmax),ybc(ibmax),ipenb(ibmax),ibfrst(nsegb+1),
     *  npoly(nsegb))

c        check for multiple polygons
      itypec=1
      do i=1,len_trim(ffile)-3
        if(ffile(i:i+3)=='xyij')itypec=2
      enddo
      write(6,"(' itypec='i2)")itypec

c        read data again
      rewind(lu)
      nsum=0 ; j=1
      do n=1,ibmax
        if(itypec==1)read(lu,*)ala,alo,ipo
        if(itypec==2)read(lu,*)alo,ala,ipo,j
        xbc(n)=alo ; ybc(n)=ala
        if(itypec==1.and.ala<0.0)then ;xbc(n)=ala ; ybc(n)=alo ;endif
        ipenb(n)=ipo ; nsum=nsum+ipo ; ns=(nsum+1)/2
        if(itypec==1)npoly(ns)=ns
        if(itypec==2)npoly(ns)=j
        if(ipo==1.and.mod(nsum,2)/=0)ibfrst(ns)=n ! first record in seg
        ip=0 ; if(n<=2.or.n>ibmax-2)ip=1
        if(ip==1)write(6,"(1x'n='i7' y,x,ipen,npoly='2f11.5,i2,i4
     *  )")n,ybc(n),xbc(n),ipenb(n),npoly(ns) !ala,alo,ipo
      enddo
      ibfrst(nsegb+1)=ibmax+1
      close (lu)

      npoly1=1 ! number of the main polygon
      do ns=1,nsegb
        if(npoly(ns)==1)npoly1=ns
      enddo
      write(6,"(' npoly1='i3)")npoly1
      do n=1,nsegb
        if(n<=3.or.n>nsegb-3)write(6,"(4x'ns='i6' ibfrst,last='2i8
     *  ' npoly='i4)")n,ibfrst(n),ibfrst(n+1)-1,npoly(n)
      enddo

c       if main polygon is not first, move it up to first
      if(npoly1/=1)then
        write(6,"('  move main polygon to top position')")
        allocate (xx(ibmax),yy(ibmax),ippx(ibmax))
        newpoly(1)=1 ; np=1
        m=0
        do n=ibfrst(npoly1),ibfrst(npoly1+1)-1
          m=m+1 ; xx(m)=xbc(n) ; yy(m)=ybc(n) ; ippx(m)=ipenb(n) 
        enddo
        do ns=1,nsegb
          if(npoly(ns)==1)cycle
          do n=ibfrst(ns),ibfrst(ns+1)-1
            m=m+1 ; xx(m)=xbc(n) ; yy(m)=ybc(n) ; ippx(m)=ipenb(n) 
          enddo
          np=np+1 ; newpoly(np)=npoly(ns)
        enddo
        do i=1,ibmax
          xbc(i)=xx(i) ; ybc(i)=yy(i) ; ipenb(i)=ippx(i)
        enddo
        npoly1=1 ; npoly(1:nsegb)=newpoly(1:nsegb)
        deallocate (xx,yy,ippx)
c        find start of each segment, isfrst()
        isum=0 ; nsum=0 ; ns=1 ; ibfrst(ns)=1
        do n=1,ibmax
          isum=isum+ipenb(n)
          if(isum<2)cycle
          isum=0 ; ns=ns+1 ; if(n<ibmax)ibfrst(ns)=n+1
        enddo
        ibfrst(ns)=ibmax+1
        write(6,"(4x'ns='i6' ibfrst,last='2i8' npoly='i4)")(n,ibfrst(n),
     *  ibfrst(n+1)-1,npoly(n),n=1,3)
        write(6,"(4x'ns='i6' ibfrst,last='2i8' npoly='i4)")(n,ibfrst(n),
     *  ibfrst(n+1)-1,npoly(n),n=nsegb-2,nsegb)
      endif

      n=ibfrst(2)-1
      allocate (xx(n),yy(n))
	  d=max(delx,dely)
	  call shift_polygon_outward(xbc(1:n),ybc(1:n),n,d,xx,yy)
	  do i=1,n
		if (xx(i)<alon0) xx(i)=alon0 
		if (xx(i)>alon3) xx(i)=alon3 
		if (yy(i)<alat0) yy(i)=alat0 
		if (yy(i)>alat3) yy(i)=alat3 
      end do
	  xbc(1:n)=xx(1:n)
	  ybc(1:n)=yy(1:n)
	  deallocate(xx,yy)

      return
      end 
	  
! =====================================================================
! Subroutine: shift_polygon_outward
! Shifts closed polygon coordinates outward by distance d along bisectors.
! =====================================================================
	  subroutine shift_polygon_outward(x, y, n, d, xx, yy)
	  implicit none
		
		! Arguments
	  integer, intent(in)  :: n          ! Total nodes (n-th node == 1st node)
	  real(8), intent(in)  :: x(n), y(n) ! Input coordinates
	  real(8), intent(in)  :: d          ! Shift distance
	  real(8), intent(out) :: xx(n), yy(n)! Shifted output coordinates
		
		! Local variables
	  integer :: i, i_prev
	  real(8) :: area_sum, loopdir
	  real(8) :: ax, ay, alen
	  real(8) :: n1(2), n2(2)
	  real(8) :: bisector(2), b_len, scale, dot_prod
		
		! Temporary array for edge unit normals
	  real(8), allocatable :: v(:,:)
		
	  allocate(v(n-1, 2))
		
		! -------------------------------------------------------------
		! 1. Determine winding direction using Shoelace Formula
		!    area_sum > 0 => CCW (loopdir = +1.0)
		!    area_sum < 0 => CW  (loopdir = -1.0)
		! -------------------------------------------------------------
	  area_sum = 0.0d0
	  do i = 1, n - 1
	  	area_sum = area_sum + (x(i) * y(i+1) - x(i+1) * y(i))
	  end do
		
	  if (area_sum >= 0.0d0) then
		loopdir = 1.0d0
	  else
		loopdir = -1.0d0
	  end if
		
		! -------------------------------------------------------------
		! 2. Calculate outward unit normal vectors for each segment
		!    Outward normal: [ay, -ax] for CCW, [-ay, ax] for CW
		! -------------------------------------------------------------
	  do i = 1, n - 1
		ax = x(i+1) - x(i)
		ay = y(i+1) - y(i)
		alen = sqrt(ax*ax + ay*ay)
			
		v(i, 1) = loopdir * (ay / alen)
		v(i, 2) = loopdir * (-ax / alen)
	  end do
		
		! -------------------------------------------------------------
		! 3. Calculate corner shift using angle bisector math
		! -------------------------------------------------------------
	  do i = 1, n - 1
			! Preceding edge index (wraps around to n-1 at vertex 1)
		if (i == 1) then
			i_prev = n - 1
		else
			i_prev = i - 1
	  end if
			
		n1 = v(i_prev, :)  ! Normal of incoming edge
		n2 = v(i, :)       ! Normal of outgoing edge
			
			! Bisector vector (sum of adjacent normals)
		bisector = n1 + n2
		b_len = sqrt(bisector(1)**2 + bisector(2)**2)
		bisector = bisector / b_len
			
			! Dot product: n1 . bisector
		dot_prod = n1(1) * bisector(1) + n1(2) * bisector(2)
			
			! Scale factor to maintain uniform distance d at corners
		scale = d / dot_prod
			
			! Apply shift to node i
		xx(i) = x(i) + scale * bisector(1)
		yy(i) = y(i) + scale * bisector(2)
	  end do
		
		! Ensure the closed ending node matches the shifted starting node
		xx(n) = xx(1)
		yy(n) = yy(1)
		
	  deallocate(v)
	  end subroutine shift_polygon_outward
c-----------------------------------------------------------------------
      subroutine rdcon
c        read the control file
      use saved_data
      real*8 z(8),degs 
      character namefil(10)*18,title*70 ! ++++
      data namefil/'ADCIRC grd','US EEZ','coastline data',
     * 'bounding polygon','output grid','output marine GTX',
     * 'output ponds GTX','output main GTX','init land/water',
     * 'init land/wat/dry'/

      write(6,*)'<rdcon>  ' ; lu=5

c        read title
      read(lu,"(a70)")title ; write(6,"(' title=',a70)")title

c        read lat, lon limits
      read(lu,*)n1,(z(n),n=1,n1)       ! lat,lon borders
      if(n1==8)then
        write(6,"(' y1,y2='2f8.3' y3,y4='2f8.3)")z(1:8)
        alat1=min(degs(z(1),z(2)),degs(z(3),z(4)))
        alat2=max(degs(z(1),z(2)),degs(z(3),z(4)))
        alon1=min(degs(z(5),z(6)),degs(z(7),z(8)))
        alon2=max(degs(z(5),z(6)),degs(z(7),z(8)))
      else if(n1==4)then
        write(6,"(' y1,y2='2f12.6' x1,x2='2f12.6)")z(1:4)
        alat1=min(z(1),z(2)) ; alat2=max(z(1),z(2))
        alon1=min(z(3),z(4)) ; alon2=max(z(3),z(4))
      endif
      write(6,"(1x,'regional latitude  limits=',2f12.6,/,
     * ' regional longitude limits=',2f12.6)")alat1,alat2,alon1,alon2
 
c        read VDatum grid intervals
      read(lu,*)delx,dely
      write(6,"(' delx='f10.6' deg = 'f8.4' nmi')")delx,
     *  60.*delx*cos(raddeg*.5*(alat1+alat2))
      write(6,"(' dely='f10.6' deg = 'f8.4' nmi')")dely,60.*dely

      read(lu,*)ielim_ponds
      write(6,"(' ielim_ponds    ='i10' (eliminate ponds) points')")
     *  ielim_ponds

c        read adjustment points:  i,j,new_value
      read(5,*)nadjust
      write(6,"(' nadjust'8x'='i10' (no. of pts to adjust)')")nadjust
      if(nadjust>0)then
        nt=0 ; allocate (iadj1(3,nadjust))
        do
          read(5,*)na,(iadj1(1:3,n),n=nt+1,nt+na) ! na=pts in this record
          nt=nt+na ; if(nt==nadjust)exit
        enddo
        write(6,"(5x'n='i6' i,j='2i5'   new value of jfield='i2)")(n,
     *  iadj1(1:3,n),n=1,nadjust)
      endif

c        read number of added layers 
      read(lu,*)layers
      write(6,"(' layers   ='i10' (additional layers of cells',
     * ' added to land)')")layers

c        barriers
      read(lu,*)nobarr1,nobarr2
      write(6,"(' nobarr1  ='i10' (remove barriers bet. pt)')")nobarr1
      write(6,"(' nobarr2  ='i10' (print barrier calcs)')")nobarr2

c        read cell print-out parameters
      read(lu,*)mpr,(iprt(m),jprt(m),m=1,min(20,mpr))
      mpr=min(20,mpr)
      write(6,"(/' number of cells to print at (max=20) ='6x,i3)")mpr
      if(mpr>=1)then
        mpr=min0(20,mpr)
        do m=1,mpr
          write(6,"(5x,'m=',i3,' i,j=',2i5)")m,iprt(m),jprt(m)
        enddo
      endif

c       read index for drying 
      read(lu,*)numnulls
      write(6,"(/' num. nodes to be dry ='6x,i3)")numnulls

c       file names
      write(6,"(/' read the file names')") ! 1=ADCIRC,2=EEZ,3=coast,4=bp,
      do k=1,9 !  5=marine,6=marineGTX,7=init water,8=init ponds,9=dry
        read(lu,"(a100)")filenm(k) ; call rname(filenm(k))
        write(6,"(2x'k='i2' file:'a18/4x,a100)")k,namefil(k),filenm(k)
      enddo

      write(6,"(' end <rdcon>')")

      return
      end
c-----------------------------------------------------------------------
      subroutine setcon
c        Set program constants
      use saved_data
      real*8 x1,deg,fmin,a,b,c
      data c/100./
      deg(x1)=int(x1)
      fmin(x1)=60.*abs(x1-int(x1))

c        determine computational array sizes
      write(6,"(/,' <setcon> '/' read-in parameters')")
      write(6,"('  alat1,2='2f13.6)")alat1,alat2
      write(6,"('  alon1,2='2f13.6)")alon1,alon2
      write(6,"('  delx   ='f13.6)")delx
      write(6,"('  dely   ='f13.6)")dely

c        set marine grid limits
      alon0=min(alon1,alon2) ; alon3=max(alon1,alon2)
      alat0=min(alat1,alat2) ; alat3=max(alat1,alat2)
      write(6,"(/' initial marine grid limits:')")
      write(6,"('  alon0,3='2f13.6)")alon0,alon3
      write(6,"('  alat0,3='2f13.6)")alat0,alat3

c        array sizes
      write(6,"(/' computed parameters')")
      imax=1+nint((alon3-alon0)/delx)
      jmax=1+nint((alat3-alat0)/dely)
      write(6,"('  initial i,jmax ='i6,i7)")imax,jmax

c        revise if necessary to cover upper lat,lon
      alon3=alon0+delx*(imax-1)
      alat3=alat0+dely*(jmax-1)
      if(alon3<alon2)then 
        imax=imax+1  
        alon3=alon0+delx*(imax-1) 
      endif
      if(alat3<alat2)then
        jmax=jmax+1 
        alat3=alat0+dely*(jmax-1)
      endif

      write(6,"('  final i,jmax ='i6,i7'  imax*jmax='i10)")imax,jmax,
     *  imax*jmax

c        VDatum grid origin
      write(6,"(/' final marine grid limits:')")
      write(6,"('  alon0,alon3='2f13.6)")alon0,alon3
      write(6,"('  alat0,alat3='2f13.6)")alat0,alat3

c        allocate arrays
      allocate (clon(imax),clat(jmax),jfield(imax,jmax),STAT=istts)
	  IF (istts /= 0) THEN
		write(6,'("Allocation 4 failed!")')
	  ELSE
		write(6,'("Allocation 4 successful")')
	  END IF
	  flush(6)

c        compute lat, lon of VDatum points
      ipr=0
      if(ipr==1)write(6,"(/1x'clon: imax='i6)")imax
      do i=1,imax
        clon(i)=alon0+delx*(i-1)
        if(ipr==1)write(6,"(3x'i='i6' clon='f12.6,3x'='f6.0,f8.2)")i,
     *  clon(i),deg(clon(i)),fmin(clon(i))
      enddo
      if(ipr==1)write(6,"(/1x'clat: jmax='i6)")jmax
      do j=1,jmax
        clat(j)=alat0+dely*(j-1)
        if(ipr==1)write(6,"(3x'j='i6' clat='f12.6,3x'='f6.0,f8.2)")j,
     *  clat(j),deg(clat(j)),fmin(clat(j))
      enddo

c        marine grid
      allocate(xv(imax,jmax),yv(imax,jmax),STAT=istts)
	  IF (istts /= 0) THEN
		write(6,'("Allocation 3 failed!")')
	  ELSE
		write(6,'("Allocation 3 successful")')
	  END IF
	  flush(6)
      do j=1,jmax ; xv(1:imax,j)=clon(1:imax) ; enddo
      do i=1,imax ; yv(i,1:jmax)=clat(1:jmax) ; enddo

c        longitude values for printing
      al4(1:4)=(/alon0,delx,alat0,dely/) ; aln(1:2)=(/alon0,delx/)
      lm(1:6)=(/1,imax,1,1,jmax,1/)

      return
      end
c-----------------------------------------------------------------------
      function degs(x,y)
C        PURPOSE - TO COMPUTE LONGITUDE AS A WHOLE NUMBER, GIVEN DEGREES
C                   AND MINUTES.  CAN BE NEGATIVE.
      real*8 x,y,degs
      degs=x+y/60. ;  if(x.lt.0.0)degs=x-y/60.
      RETURN
      END
c----------------------------------------------------------------------
      subroutine reduce_segments
c        eliminate a segment if it all lies outside the lat-lon window
      use saved_data
      real*8 xp(ismax),yp(ismax)
      integer ipx(ismax)
 
      write(6,"(' <reduce_segments>  ismax='i8' nsegs='i8)")ismax,nsegs
      k=0 ! points in all segments
      do n=1,nsegs
        if(alon0>xmaxs(n).or.alon3<xmins(n).or.
     *  alat0>ymaxs(n).or.alat3<ymins(n))cycle
        m1=1           ! find start (m1), end (m2) of a segment
        if(n>1)m1=isend(n-1)+1 ; m2=isend(n)
        do m=m1,m2
          k=k+1 ; xp(k)=xcs(m) ; yp(k)=ycs(m)  ; ipx(k)=ipens(m)
        enddo
      enddo

c        new number of points
      ismax=k

c        find new segments
      isum=0
      iset=0 ! number of sets
      do i=1,ismax
        if(isum==0)istart=i
        isum=isum+ipx(i)
c          check for end of segment
        if(isum==2)then
c          increment set
          iset=iset+1 ; isend(iset)=i
c            reset
          isum=0
        endif
      enddo
      do i=1,iset
c       write(6,"('  iset=',i5,' isend=',i9)")i,isend(i)
      enddo

c        reset here
      nsegs=iset
      xcs(1:ismax)=xp(1:ismax) ;  ycs(1:ismax)=yp(1:ismax)
      ipens(1:ismax)=ipx(1:ismax)
      call maxmin(xcs,ycs,ipens,ismax,xmins,xmaxs,ymins,ymaxs,nsegs)

      return
      end
c----------------------------------------------------------------------
      subroutine pril2(ttl,ia,imx,jmx,lm2,aln)
c     print integer array, with skipping and subsetting.
c       automatic 'le'. 
c     input variables:
c            ttl = title of plot
c          ia(,) = integer array
c        imx,jmx = max array dimensions
c          lm2() = i,j array limits to plot: i1,i2,i3,  j1,j2,j3
c          aln() = x0,dx
      real*8 aln(4),x1,x2,fd
      character ttl*(*),fmo*25,fc1*25,fc2*25,fl1*35,ch(jmx)*1
      integer ia(imx,jmx),icol(imx),lm2(6),ib(imx,jmx),jrow(jmx)
      data iw/80/ ! max width of page (number of characters)
      fd(x1)=60.*(abs(x1)-int(abs(x1)))
c        copy ia into ib to create subset
      ii=0 
      do i=lm2(1)+lm2(3)-1,lm2(2),lm2(3)
        ii=ii+1 ; icol(ii)=i ; jj=0
        do j=lm2(4)+lm2(6)-1,lm2(5),lm2(6)
          jj=jj+1 ; ib(ii,jj)=ia(i,j) ; jrow(jj)=j
        enddo
      enddo
      im2=ii ; jm2=jj
c        get max abs of ib(,) for le
      ia1=100000 ; ia2=-100000
      do i=1,im2 ; do j=1,jm2
        ia1=min(ia1,ib(i,j)) ; ia2=max(ia2,ib(i,j))
      enddo ; enddo
      ibmax=max(abs(ia1),abs(ia2)) ; le=1+int(1.+log10(float(ibmax))) 
      if(ibmax<10)le=1
      write(6,"(' ia1,ia2,ibmax='3i10)")ia1,ia2,ibmax
c        number of digits in longest row number[jrow(jm2)], nperl, le
      nd1=int(1.+log10(float(jrow(jm2))))
      nperl=(iw-2-nd1)/le  ! max number per line= 80 =nd1+2+nperl*len
c        number of digits in longest column number[icol(im2)],nper2,le2,is
      nd2=int(1.+log10(float(icol(im2))))
      le2=le ; nper2=nperl ; is=1
      if(nd2>=le)then ! if digits too large for column, compute skip factor
        do is=1,5
          le2=le*is ; if(le2>nd2)exit
        enddo
        nper2=(iw-2-nd1)/le2 
      endif
c        format for ib(,), header columns, trailer columns
      fmo='(i'//cc(nd1)//',1xa1,'//cc(nperl)//'i'//cc(le)//')'
      fc1='(i'//cc(nd1)//',2h: ,'//cc(nper2)//'i'//cc(le2)//')'
      fc2='(i'//cc(nd1)//',2h. ,'//cc(nper2)//'i'//cc(le2)//')'
c     write(6,"(/' fmo='a25/)")fmo
      ch=chget(jmx,aln,lm2) 
c     write(6,"(' ch='60a1)")ch(jrow(jm2):jrow(1):-1)
c     write(6,"('  jm2='i8'  jrow(jm2)='i8)")jm2,jrow(jm2)
c        loop thru rows here
      kmax=1+(im2-1)/nperl 
      do k=1,kmax
        if(k==1)write(6,"(/'<pril2>'2xa25,35x'(kmax='i4')')")ttl,kmax
        if(k> 1)write(6,"(/5xa25)")ttl
        i1=1+nperl*(k-1) ; i2=min(i1+nperl-1,im2) ! lims for ib
        i3=i1+is-1 ; i4=i2 ; i5=is   ! lims for column headers
c          print longitudes at ends
        x1=aln(1)+(icol(i1)-1)*aln(2) ; x2=aln(1)+(icol(i2)-1)*aln(2) 
        ic=nd1+le+1 ; lc=nd1+1+min(nperl,i2-i1+1)*le ! col of first, last '|'
        if(lc>iw)lc=lc-le ; ns=max(10,lc-ic-23)
        fl1='('//cc(nd1+le+1)//'x2h|<,i4,f6.2,'//cc(ns)//
     *    'x,i4,f6.2,2h>|)'
        if(ns> 10)write(6,fl1)int(x1),fd(x1),int(x2),fd(x2)
        if(ns<=10)write(6,fl1)int(x1),fd(x1)
        write(6,fc1)k,icol(i3:i4:i5) ! top column line
        js=jrow(jm2)+1 ! ++++
        do j=jm2,1,-1
          js=js-1 ; write(6,fmo)jrow(j),ch(js),ib(i1:i2,j) ! +++
        enddo
        write(6,fc2)k,icol(i3:i4:i5) ! bottom column line
      enddo
      write(6,"(3xa25' iamin,max,ibmax='3i6' ;')")ttl,ia1,ia2,ibmax
      return
      contains
       function cc(i) ! return a 2-char string for a 2-digit integer
        character cc*2
        i0=ichar('0') ; cc=char(i/10+i0)//char(mod(i,10)+i0)
       end function cc
      function chget(jmx,aln,lm2)
c        lat position like: '-75+12.3'    where + is at latitude
      real*8 aln(4),alat,fmin
      character ch(jmx)*1,cc(9)*1,chget(jmx)*1
      integer lm2(6),leng(3),istart(3)
      data leng/3,2,1/,istart/1,5,8/,intvl/10/ ! interval between '+'
      j1=lm2(5) ; j2=lm2(4) ; j3=-lm2(6)  ;            ipr=0
      i0=ichar('0') ; ch(1:jmx)=' ' ; js=j1+1
      do j=j1,j2,j3 
        alat=aln(3)+(j-1)*aln(4) ; ideg=int(alat) ; fmin=(alat-ideg)*60.
        cc(1:9)=' ' ; cc(5)='+' ; cc(8)='.' ; js=js-1
        do mm=1,3 ! degrees (4 places), minutes (2), tenths of min (1)
          if(mm==1)nump=abs(ideg) ; if(mm==2)nump=int(fmin)
          if(mm==3)nump=10.*(fmin-int(fmin)) ;numdig=0 ; le=leng(mm) 
          do l=1,le
            k=abs(nump)/10**(le-l)-10*(abs(nump)/10**(le-l+1))
            if(k>0)numdig=numdig+1;if(k==0.and.numdig>0)numdig=numdig+1
            if(numdig>0.or.l==le)cc(istart(mm)+l)=char(i0+k)
          enddo
          if(mm==1.and.ideg<0)cc(4-numdig)='-'
        enddo
        iout=0 ; if(mod(js-(j1-3),intvl)==0.and.j>4)iout=1
        if(iout==1)ch(js+4:js-4:-1)=cc(1:9)
        if(ipr==1)write(6,"(' j,js='2i5,4x'lat='f8.3,4x'deg,min='i4,f6.1
     *  ,3x'cc='9a1' io='i1)")j,js,alat,ideg,fmin,cc(1:9),iout
      enddo
      chget=ch ! to print, use js: js=j1+1; do j=j1,j2,j3 ; js=js-1
               ! do j=j1,j2,j3 ; js=js-1 ; write(6,)j,ch(js) ; enddo
      if(ipr==1)write(6,"('  j1+1='i8'  js='i8)")j1+1,js
      if(ipr==1)write(6,"('  ch='60a1)")ch(j1+1:js:-1)
      end function 
      end
c----------------------------------------------------------------------
      subroutine pri(ttl,ia,imx,jmx,lm2,aln)
c     print integer array, with skipping and subsetting, automatic 'le'. 
c     input variables:
c            ttl = title of plot
c          ia(,) = integer array
c        imx,jmx = max array dimensions
c          lm2() = i,j array limits to plot: i1,i2,i3,  j1,j2,j3
c          aln() = x0,dx
      real*8 aln(2),x1,x2,fd
      character ttl*(*)
      character*25 fmo*25,fc1*25,fc2*25,fl1*35
      integer ia(imx,jmx),icol(imx),lm2(6),ib(imx,jmx),jrow(jmx)
      data iw/80/ ! max width of page (number of characters)
      fd(x1)=60.*(abs(x1)-int(abs(x1)))
c        copy ia into ib to create subset
      ii=0 
      do i=lm2(1)+lm2(3)-1,lm2(2),lm2(3)
        ii=ii+1 ; icol(ii)=i ; jj=0
        do j=lm2(4)+lm2(6)-1,lm2(5),lm2(6)
          jj=jj+1 ; ib(ii,jj)=ia(i,j) ; jrow(jj)=j
        enddo
      enddo
      im2=ii ; jm2=jj
c        get max abs of ib(,) for le
      ia1=100000 ; ia2=-100000
      do i=1,im2 ; do j=1,jm2
        ia1=min(ia1,ib(i,j)) ; ia2=max(ia2,ib(i,j))
      enddo ; enddo
      ibmax=max(abs(ia1),abs(ia2)) ; le=1+int(1.+log10(float(ibmax))) 
c        number of digits in longest row number[jrow(jm2)], nperl, le
      if(ibmax<10)le=1
      nd1=int(1.+log10(float(jrow(jm2))))
      nperl=(iw-1-nd1)/le  ! max number per line= 80 =nd1+1+nperl*len
c        number of digits in longest column number[icol(im2)],nper2,le2,is
      nd2=int(1.+log10(float(icol(im2))))
      le2=le ; nper2=nperl ; is=1
      if(nd2>=le)then ! if digits too large for column, compute skip factor
        do is=1,5
          le2=le*is ; if(le2>nd2)exit
        enddo
        nper2=(iw-1-nd1)/le2 
      endif
c        format for ib(,), header columns, trailer columns
      fmo='(i'//cc(nd1)//',1h ,'//cc(nperl)//'i'//cc(le)//')'
      fc1='(i'//cc(nd1)//',1h:,'//cc(nper2)//'i'//cc(le2)//')'
      fc2='(i'//cc(nd1)//',1h.,'//cc(nper2)//'i'//cc(le2)//')'
c        loop thru rows here
      kmax=1+(im2-1)/nperl 
      do k=1,kmax
        if(k==1)write(6,"(/'<pri>'2xa15,45x'(kmax='i3')')")ttl,kmax
        if(k> 1)write(6,"(/7xa15)")ttl
        i1=1+nperl*(k-1) ; i2=min(i1+nperl-1,im2) ! lims for ib
        i3=i1+is-1 ; i4=i2 ; i5=is   ! lims for column headers
c          print longitudes at ends
        x1=aln(1)+(icol(i1)-1)*aln(2) ; x2=aln(1)+(icol(i2)-1)*aln(2) 
        ic=nd1+le+1 ; lc=nd1+1+min(nperl,i2-i1+1)*le ! col of first, last '|'
        if(lc>iw)lc=lc-le ; ns=max(10,lc-ic-23)
        fl1='('//cc(nd1+le)//'x2h|<,i4,f6.2,'//cc(ns)//'x,i4,f6.2,2h>|)'
        if(ns> 10)write(6,fl1)int(x1),fd(x1),int(x2),fd(x2)
        if(ns<=10)write(6,fl1)int(x1),fd(x1)
        write(6,fc1)k,icol(i3:i4:i5) ! top column line
        do j=jm2,1,-1
          write(6,fmo)jrow(j),ib(i1:i2,j)
        enddo
        write(6,fc2)k,icol(i3:i4:i5) ! bottom column line
      enddo
      write(6,"(7xa15,5x';'3x'iamin,iamax,ibmax='3i8)")ttl,ia1,ia2,ibmax
      return
      contains
       function cc(i) ! return a 2-char string for a 2-digit integer
        character*2 cc
        i0=ichar('0') ; cc=char(i/10+i0)//char(mod(i,10)+i0)
       end function cc
      end
c-----------------------------------------------------------------------
      subroutine rname(fname)
c        remove trailing, non-blank characters in file name
c        if separated by at least 1 blank space
      character fname*100
      ib=101
      do i=1,100
        if(fname(i:i)==' ')ib=i ; if(i>=ib)fname(i:i)=' '
      enddo
      return
      end    
c-----------------------------------------------------------------------
      function ipri(i)
c        determine print index: 0=no, 1=yes if i=iprt
      use saved_data
      ipri=0 ; if(mpr==0)return
      do mm=1,mpr ;if(i==iprt(mm))ipri=1 ;enddo
      return
      end
c-----------------------------------------------------------------------
      function iprj(j)
c        determine print index: 0=no, 1=yes if j=jprt
      use saved_data
      iprj=0 ; if(mpr==0)return
      do mm=1,mpr ;if(j==jprt(mm))iprj=1 ;enddo
      return
      end
c----------------------------------------------------------------------
      function iprij(i,j)
c        determine print index: 0=no, 1=yes  if i=iprt and j=jprt
      use saved_data
      iprij=0 ; if(mpr==0)return
      do mm=1,mpr ; if(i==iprt(mm).and.j==jprt(mm))iprij=1 ;enddo
      return
      end
c----------------------------------------------------------------------
      subroutine inside_pert(nx,xi,yi,x,y,i,ncx,nca,xcr)! 5.9.2013  K. Hess
c        See if point (xi,yi) is inside the polygon defined by the 'nx'
c        points (x,y). Uses a perturbation method: if a node value y(n)
c        equals yi, the node value is shifted up by the amount 'del'.
c        Point on west side is in. Point on horizontal at top is in.
c        nx = num of points in segment
c        nc = number of crossings to right of xi
c       nca = number of all crossings
c     xcr() = all crossings
c         i = result: 0=outside, 1=inside
      real*8 xi,yi,x(nx),y(nx),xc,ys(nx),del,xcr(ncx)
      character ttl(2)*7
      data ttl/'OUTSIDE','INSIDE '/,del/0.000002/
      logical p ; data p/.false./
c       store y(n) in temporary array
      ys(1:nx)=y(1:nx)
c       if a node value ys(n) equals yi, add a small amount to ys(n)
      nadd=0
      do n=1,nx
        if(yi/=ys(n))cycle
        ys(n)=ys(n)+del ; nadd=nadd+1
      enddo 
c       loop thru all points defining the polygon
      nc=0  ; nca=0 ! number of crossings
      do n=1,nx-1
c          yi value is between nodes    
        if(yi>min(ys(n),ys(n+1)).and.yi<max(ys(n),ys(n+1)))then
          if(x(n)==x(n+1))then
            xc=x(n)
          else
            xc=x(n)+(x(n+1)-x(n))*(yi-ys(n))/(ys(n+1)-ys(n)) 
          endif
          nca=nca+1 ; xcr(nca)=xc ; if(xc>xi)nc=nc+1
          if(nca>ncx)then ; write(6,"(/' ** crossings too large. nca='i5
     *    ' ncx='i5/' ** stop **')")nca,ncx ; stop ; endif
        endif
      enddo
c        inside if number of crossings to east is odd
      i=mod(nc,2) 
c        print result
      if(p)write(6,"(3x'nx='i5' nc='i3' i='i2' nadd='i5,3xa7)")nx,nc,i,
     * nadd,ttl(i+1)
      return
      end
c-----------------------------------------------------------------------
      subroutine ck_nodes
c     Find ADCIRC nodes that have gone dry and are in vgrid.
c      1.Set all points to wet: idry(i,j)=0.
c      2.If number of nodes in element containing marine grid point
c        is >= 'numnulls', set marine point to dry: idry=1.

c     variables:
c               icells(1,n) = number of cells around node n
c   icells(2:icells(1,n),n) = i-th cells around node n
c                   idry(,) = index for marine grid: 0=wet, 1=dry
c                    isv(,) = idx:-1=dry,0=land,1=water,2=1st layer, etc.
c           ncellnodes(i,n) = i-th node around cell n
c                      ndry = number of nodes within marine grid that dried

      use saved_data
      real*8 dx,dy,x1,x2,y1,y2,x(3),y(3)
      integer ijmx(4),nod(3),ia(imax,jmax)
      logical p

      write(6,"(/' <ck_nodes>')")
      allocate (idry(imax,jmax))
      dx=delx ; dy=dely

c        initialize index for having gone dry: 0=no, 1=yes
      idry(1:imax,1:jmax)=0 ; ia(1:imax,1:jmax)=0

c        search
      ndry=0
      do n=1,nodes

c         find ADCIRC nodes that went dry
        if(ds(n,1)<9.0)cycle ! bad node value is 9.999
        ill=1.+(xs(n)-xv(1,1))/dx  ;  jll=1.+(ys(n)-yv(1,1))/dy
        if(.not.(ill>0.and.ill<imax.and.jll>0.and.jll<jmax))cycle
        ndry=ndry+1
c         Are any of 4 nearest marine points in a ADCIRC element around
c          this ADCIRC node?
c         If yes, count number of ADCIRC nodes that went dry. If that 
c          number is GE numnulls, set marine grid point to dry.
        nc=icells(1,n)
        p=.false. ; if(ndry==1)p=.true.
        if(p)write(6,"(' ndry='i7/' node='i8' ds='f7.4' xs,ys='2f11.6
     *   ' ill,jll='2i5)")ndry,n,ds(n,1),xs(n),ys(n),ill,jll
        if(p)write(6,"(3x'num of cells='i2' cells='10i8)")nc,
     *   icells(2:1+nc,n)
        do i=ill,ill+1 ; do j=jll,jll+1
          if(idry(i,j)==1)cycle       ! skip if already set to dry
          x1=xv(i,j) ; y1=yv(i,j)
          if(p)write(6,"(3x'i,j='2i6' jfield='i2' x1,y1='2f12.6)")i,j,
     *    jfield(i,j),x1,y1
          do k=1,nc ! loop thru cells around the node
            kc=icells(1+k,n)            ! cell #
            nod(1:3)=ncellnodes(1:3,kc) ! nodes in cell kc
            x(1:3)=xs(nod(1:3)) ; y(1:3)=ys(nod(1:3))
            call inside_tri3(0,x,y,x1,y1,inn,3)
            if(p)write(6,"(5x'k='i2' cell='i7' nodes='3i8' inn='i2)")
     *      k,kc,nod(1:3),inn
            if(inn==0)cycle
            num9=0 ! number of nodes that went dry
            do kk=1,3 ;if(ds(nod(kk),1)>9.0)num9=num9+1 ; enddo
            if(p)write(6,"(7x'z='3f8.3' # dry='i2)")ds(nod(1:3),1),num9
            if(num9>=numnulls)then
              idry(i,j)=1 ; ia(i,j)=num9
              if(p)write(6,"(7x'set to dry')")
            else
              if(p)write(6,"(7x'not set to dry')")
            endif
          enddo
        enddo ; enddo

      enddo

c       find limits of points that have dried
      ijmx(1:4)=(/imax,0,jmax,0/)
      do i=1,imax ; do j=1,jmax
        if(idry(i,j)==1)ijmx(1:4)=(/min(ijmx(1),i),max(ijmx(2),i),
     *  min(ijmx(3),j),max(ijmx(4),j)/)
      enddo ; enddo
      write(6,"(' ndry='i8' ijmx='4i8)")ndry,ijmx(1:4)

c        inspect points. if a non-dry point is surrounded by all
c        drying, switch to drying.
      do i=2,imax-1 ; do j=2,jmax-1
        if(idry(i,j)==1)cycle ; if(idry(i-1,j)+idry(i+1,j)+idry(i,j-1)
     *  +idry(i,j+1)==4)idry(i,j)=1
      enddo ; enddo

c        reset jfield     
      do i=1,imax ; do j=1,jmax
        if(idry(i,j)==1)jfield(i,j)=-1
      enddo ; enddo     

      write(6,"(//10x' DRY POINTS IN MARINE GRID')")
c     call pril2('idry',idry,imax,jmax,lm ,al4)
c     do i=1,imax ; do j=1,jmax
c       ia(i,j)=jfield(i,j);  if(idry(i,j)==1)ia(i,j)=-1
c     enddo ; enddo
c      call pril2(' number of nulls',ia,imax,jmax,lm ,al4)
c      call pril2(' jfield, dry=*',jfield,imax,jmax,lm ,al4)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_icell
c        find the cells around node n; save in icells(10,n).
c     variables:
c                icells(1,n) = number of cells around node n
c    icells(2:icells(1,n),n) = i-th cells around node n
c            ncellnodes(i,n) = i-th node around cell n
      use saved_data

      write(6,"(/' <create_icell> ')")

c        initialize cells at nodes
      icells(1:10, 1:nodes)=0

c        loop thru cells
      do n=1,nt
        do i=1,3   ! loop thru the 3 nodes at this cell
          no=ncellnodes(i,n)          ! node number at the cell
          nn=icells(1,no)             ! present number of cells at node
          nn=nn+1                     ! update number of cells at node
		  if(nn>9) then
			write(6,"(' too many cells at node 'i8'')")no
			flush(6)
			continue
		  endif
          icells(1   ,no)=nn   ! save new number of cells
          icells(1+nn,no)=n    ! save new cell number
        enddo
      enddo

c        print cells at each node
      do n=1,nodes
        ipr=0  ; if(n<=3.or.n>=nodes-3)ipr=1
        if(ipr==1)write(6,"('  n='i8' #='i3' cells='9i8
     *  )")n,icells(1,n),(icells(1+i,n),i=1,min(9,icells(1,n)))
      enddo

      return
      end
c--------------------------------------------------------------------
      subroutine inside_tri3(ipr,x,y,xp,yp,inn,idmx)
c        See if point (xp,yp) is inside the triangle x(i),y(1), i=1,2,3.
c        If inside, in=1; otherwise, in=0.
c        If point is on one side, in=0.
      real*8 x(idmx),y(idmx),xp,yp
      r(z)=z+a180*(1.-sign(1.,z))

      if(ipr==1)then
        write(6,"(/,'  <inside_tri3>'/'   xp,yp =',2f12.6)")xp,yp
        write(6,"('   x,y   =',2f12.6)")(x(i),y(i),i=1,3)
      endif

      inn=0
      a180=atan2(0.,-1.)

c        test angles around point 1
      a12=atan2(y(2)-y(1),x(2)-x(1))
      a13=atan2(y(3)-y(1),x(3)-x(1))
      d1=r(a12-a13)
      d2=r(a13-a12)
      if(d1>a180)then
        d=d2 ; a0=a12
      else
        d=d1 ; a0=a13
      endif
      a1p=atan2(yp-y(1),xp-x(1)) ; ap=r(a1p-a0)
      if(ap>d)return

c        test angles around point 2
      a23=atan2(y(3)-y(2),x(3)-x(2)) ; a21=atan2(y(1)-y(2),x(1)-x(2))
      d1=r(a23-a21) ; d2=r(a21-a23)
      if(d1>a180)then
        d=d2 ; a0=a23
      else
        d=d1 ; a0=a21
      endif
      a2p=atan2(yp-y(2),xp-x(2)) ; ap=r(a2p-a0)
      if(ap>d)return
      inn=1

      return
      end
c--------------------------------------------------------------------
      subroutine rd_adcirc(lu,nf)
c        read a triangular grid
c     variables:
c                   iallocsor = index to allocate arrays(0=no, 1=yes)
c             ncellnodes(i,n) = i-th node around cell n
c                          nf = number of file being read
c                 icells(1,n) = number of cells around node n
c     icells(2:icells(1,n),n) = i-th cells around node n
c                   xs(),ys() = long, lat of node
c                    ds(n,nf) = datum for node 'n' and datum 'nf'
      use saved_data
      character ttle*70
      real*8 xa,xb,ya,yb

      write(6,"(/' <rd_adcirc>    lu='i4'  allocsor='i2)")lu,iallocsor

c        read header
      read(lu,"(a70)",iostat=io)ttle    
      write(6,"('  ttle='a70)")ttle

c        read number of triangles and nodes
      read(lu,*)nt,nodes
      write(6,"('  nt(# of triangles)='i9' no(# of nodes)='i9)")nt,nodes

      icallcell=0               ! index to call sub. create_icell
      if(iallocsor==0)then
        ndatums=1
        write(6,"('  allocate arrays. ndatums='i3)")ndatums
        allocate (xs(nodes),ys(nodes),
     *  ds(nodes,ndatums),ncellnodes(3,nt),icells(10,nodes),STAT=istts)
        iallocsor=1 ; icallcell=1   
	  IF (istts /= 0) THEN
		write(6,'("Allocation 10 failed!")')
	  ELSE
		write(6,'("Allocation 10 successful")')
	  END IF
	  flush(6)
		
      endif

c        read the vertices and datum elevations
      ip=2  ! print first 'ip' rows
      xa= 360. ;  ya= 360. ;  xb=-360. ; yb=-360.
      do n=1,nodes
        read(lu,*)nn,xs(n),ys(n),ds(n,nf)
        xa=min(xa,xs(n)) ; ya=min(ya,ys(n))
        xb=max(xb,xs(n)) ; yb=max(yb,ys(n))
        ipr=0 ; if(n<=ip.or.n>nodes-ip)ipr=1
        if(ipr==1)write(6,"(5x'n='i7' xs='f12.6' ys='f10.6' ds='f8.4)")
     *  nn,xs(n),ys(n),ds(n,nf)
		flush(6)
      enddo
      write(6,"(/' ADCIRC grid min,max xs='2f12.6)")xa,xb
      write(6,"(' ADCIRC grid min,max ys='2f12.6)")ya,yb
	  flush(6)
c        read nodes in each element
      write(6,"(/'  read element nodes')")
      do n=1,nt
        read(lu,*)nn,np,(ncellnodes(i,nn),i=1,np)
        if(n<=3.or.n>=nt-3+1)write(6,"(4x'n='i8' ncellnodes='6i8)")
     *  n,(ncellnodes(i,nn),i=1,np)
      enddo
      flush(6)

c        determine cells around each node
      if(icallcell==1)call create_icell

      write(6,"(/' ADCIRC grid read-in complete')")
      flush(6)
      
      return
      end
c--------------------------------------------------------------------
      subroutine segments_coast
c       divide coastline string into segments, get min & max
c       for each segment, eliminate segments outside window (saves time)
c     variables -
c          ismax = number of points (initial, final)
c          nsegc = number of segments (initial, final)
c    xcs(),ycs() = points in coastline file
c      xx(),yy() = temporary arrays
c         xmc(,) = lat, lon extrema for segment
c        ipens() = pen commands
c         ippx() = temporary array
      use saved_data
      real*8 d
      write(6,"(/' <segments_coast>')")
      write(6,"('  number of segments(nsegc)='i6)")nsegc
      allocate (isfrst(nsegc+1),xmc(4,nsegc))
      allocate (xx(ismax),yy(ismax),ippx(ismax))

c        find start of each segment, isfrst()
      isum=0 ; nsum=0 ; ns=1 ; isfrst(ns)=1
      do n=1,ismax
        isum=isum+ipens(n)
        if(isum<2)cycle
        isum=0 ; ns=ns+1 ; if(n<ismax)isfrst(ns)=n+1
      enddo
      isfrst(ns)=ismax+1
      write(6,"('  original segments: starting point')")
      write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=1,3)
      write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=nsegc-2,nsegc+1)

c       find min and max lat and lon of each segments
      do ns=1,nsegc
        xmc(1:4,ns)=(/360.,-360.,360.,-360./)
        do n=isfrst(ns),isfrst(ns+1)-1
          xmc(1:4,ns)=(/min(xmc(1,ns),xcs(n)),max(xmc(2,ns),xcs(n)),
     *    min(xmc(3,ns),ycs(n)),max(xmc(4,ns),ycs(n))/)
        enddo
      enddo
      write(6,"('  min and max long and lat')")
      write(6,"(4x'ns='i6' xmc='4f12.5)")(n,xmc(1:4,n),n=1,3)
      write(6,"(4x'ns='i6' xmc='4f12.5)")(n,xmc(1:4,n),n=nsegc-2,nsegc)

c       save only if min, max are inside the grid
      icheck_in=1
      if(icheck_in==1)then ! ++++++++++++++++++++++++++++++++++++++
      nn=0   ; nstmp=0
      write(6,"('  grid: alon0,alon3,alat0,alat3='4f11.5)")alon0,alon3,
     *  alat0,alat3
      do ns=1,nsegc
        if(xmc(1,ns)>alon3.or.xmc(2,ns)<alon0.or.xmc(3,ns)>alat3.or.
     *  xmc(4,ns)<alat0)cycle
        do n=isfrst(ns),isfrst(ns+1)-1
          nn=nn+1 ; xx(nn)=xcs(n) ; yy(nn)=ycs(n) ; ippx(nn)=ipens(n)
        enddo
        nstmp=nstmp+1
      enddo
      ismax=nn ; xcs(1:ismax)=xx(1:ismax) ; ycs(1:ismax)=yy(1:ismax)
      ipens(1:ismax)=ippx(1:ismax) ; nsegc=nstmp
      write(6,"('  revised number of points   inside grid='i8)")ismax
      write(6,"('  revised number of segments inside grid='i8)")nsegc
      if(nsegc==0)stop
      endif ! icheck_in

c        find start of each segment, isfrst()
      isum=0 ; nsum=0 ; ns=1 ; isfrst(ns)=1
      do n=1,ismax
        isum=isum+ipens(n) ; if(isum<2)cycle
        isum=0 ; ns=ns+1 ; if(n<ismax)isfrst(ns)=n+1
      enddo
      isfrst(ns)=ismax+1
      write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=1,3)
      write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=nsegc-2,nsegc+1)

c        check for segments that do not meet at ends
      write(6,"('  check for ends not meeting')")
      nn=0
      do ns=1,nsegc
        d=abs(xcs(isfrst(ns))-xcs(isfrst(ns+1)-1))+
     *    abs(ycs(isfrst(ns))-ycs(isfrst(ns+1)-1))
        if(d<0.00001)cycle
        nn=nn+1 ;write(6,"(i6' ends dont meet: n='i5' d='f10.6)")nn,ns,d
        write(6,"(10x'x,y(first)='2f12.6)")xcs(isfrst(ns)),
     *  ycs(isfrst(ns))
        write(6,"(10x'x,y(last )='2f12.6)")xcs(isfrst(ns+1)-1),
     *  ycs(isfrst(ns+1)-1)
      enddo
      if(nn>0)then
        write(6,"('  add points for ends not meeting')")
        na=0
        do ns=1,nsegc
          n1=isfrst(ns) ; n2=isfrst(ns+1)-1
          do n=n1,n2
            na=na+1 ; xcs(na)=xx(n) ; ycs(na)=yy(n) ; ipens(na)=ippx(n)
          enddo
          d=abs(xx(n1)-xx(n2))+abs(yy(n1)-yy(n2))
          if(d<0.00001)cycle
          ipens(na)=0
          na=na+1 ; ipens(na)=1 ;xcs(na)=xx(n1) ; ycs(na)=yy(n1)
        enddo
        ismax=na
        write(6,"('  revised number of points   inside grid='i8)")ismax

c        find start of each segment, isfrst()
        isum=0 ; nsum=0 ; ns=1 ; isfrst(ns)=1
        do n=1,ismax
          isum=isum+ipens(n)
          if(isum<2)cycle
          isum=0 ; ns=ns+1 ; if(n<ismax)isfrst(ns)=n+1
        enddo
        isfrst(ns)=ismax+1
        write(6,"('  revised segments: starting point')")
        write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=1,3)
        write(6,"(4x'ns='i6' isfrst='i8)")(n,isfrst(n),n=nsegc-2,
     *  nsegc+1)
      endif

c        write new coastline file to lu=55
      isum=0 ; ns=1
      do n=1,ismax
c       write(55,"(2f12.6,i2,i6)")xcs(n),ycs(n),ipens(n),ns ! ++++++
        isum=isum+ipens(n) ; if(isum<2)cycle
        isum=0 ; ns=ns+1
      enddo

      deallocate (xx,yy,ippx)
      return
      end
c----------------------------------------------------------------------
      subroutine wrtout_gtx(k)
c        write jfield in GTX format
c     jfield(,) = land/water field: -1=intertidal, 0=land, 1=water, 
c                  2,3,etc=water layer
c             k = output file number
      use saved_data
      USE NETCDF        ! enables the fortran 90 interface to netcdf
	  
 	  integer ncid,xDimId,yDimId,bxId,zDims(2),bid,zid,ncerr,n
	  real*8 box(4)
	  character(len=100) :: fname
	  integer*1,allocatable,dimension(:,:) :: jz

      if(filenm(k)=='none')return
	  fname=filenm(k)
	  n=len(trim(filenm(k)))
	  fname((n+1):(n+3))='.nc'
	  n=n+3
	  
	  allocate (jz(imax,jmax))
	  jz=jfield
	  
      write(6,"(/' <wrtout_gtx>         k='i2/2xa78)")k,fname(1:n)
	  flush(6)
		! Create NetCDF file
		ncerr = nf90_create(fname(1:n),nf90_CLOBBER+NF90_HDF5, ncid)
		call errhandle(ncerr)
		! Define dimensions
		call errhandle(nf90_def_dim(ncid, 'box', 4, bxId))
		call errhandle(nf90_def_dim(ncid, 'longitude', imax, xDimId))
		call errhandle(nf90_def_dim(ncid, 'latitude', jmax, yDimId))
		zDims(1)=xDimId
		zDims(2)=yDimId				
		! Define variables
		box(1)=clat(1)
		box(2)=clon(1)+360.
		box(3)=dely
		box(4)=delx
		ncerr=nf90_def_var(ncid,'box',nf90_DOUBLE,bxId,bid)
		call errhandle(ncerr)	
		ncerr=nf90_def_var(ncid,'field',nf90_BYTE,zDims,zid)
		call errhandle(ncerr)		
		ncerr=nf90_def_var_deflate(ncid,zid,0,1,8)
		call errhandle(ncerr)		
		! End definitions
		call errhandle(nf90_enddef(ncid))
		! Write data to variables
		call errhandle(nf90_put_var(ncid, bid, box))
		call errhandle(nf90_put_var(ncid, zid, jz))
		! Close NetCDF file
		call errhandle(nf90_close(ncid))

		write(6,"('  NetCDF file created successfully')")
		flush(6)
		deallocate(jz)

      return
      end
c----------------------------------------------------------------------
      subroutine ztime(ne1,idx)
c        write time and duration
      integer ne1(8),ne2(8)   ! computer time variables
      call date_and_time(values=ne2)
      nethr=ne2(5)-ne1(5) ; if(ne2(3)/=ne1(3))nethr=nethr+24
      ndur=3600*nethr+60*(ne2(6)-ne1(6))+ne2(7)-ne1(7)
      i=ndur/3600 ; j=(ndur-3600*i)/60 ; k=ndur-3600*i-60*j
      if(idx==1)write(6,"(/3x'duration: hr='i2' min='i2' sec='i2)")i,j,k
      if(idx==2)write(6,"(1x'date='i2'/'i2'/'i4'  hr(EST)='i2':'i2,9x'du
     *ration: hr='i2' min='i2' sec='i2)")ne2(2:3),ne2(1),ne2(5:6),i,j,k
      return
      end
c-------------------------------------------------------------------
      subroutine errhandle(istts)
      USE NETCDF        ! enables the fortran 90 interface to netcdf
	  
		integer istts
		if (istts.ne.nf90_NOERR) then
			write(6,*) nf90_strerror(istts)
			flush(6)
			stop
		endif
	  end
