c---------------------------------------------------------------
      subroutine sideset_yp(yplus,ypmin,ypmax,sslist,nss,dssopt,wdopt)
      implicit none
c  sslist - sideset list, integer array of sidesets where we need traction
c  nss - number of sidesets i.e. the size of the sslist array
c  dssopt - direct stiffness sum opt, if dssopt>1, a dsavg on traction is performed

      include 'SIZE'
      include 'TOTAL'

      real yplus(lx1,ly1,lz1,lelv), ypmin, ypmax
      integer sslist(1), nss, dssopt, wdopt

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)
      integer nij,bid,jf1,js1,jskip1,js2,jf2,jskip2
      integer e,f,i,j1,j2,ia,ntot

      real tauij(lx1,ly1,lz1,lelv,6),
     $      sij(lx1,ly1,lz1,6,lelv),
     $      ur(lxyz),us(lxyz),ut(lxyz),
     $      vr(lxyz),vs(lxyz),vt(lxyz),
     $      wr(lxyz),ws(lxyz),wt(lxyz)

      real wd(lx1,ly1,lz1,lelv)

      real visc_st_x, visc_st_y, visc_st_z, normal_st,
     $     shear_st_x, shear_st_y, shear_st_z, shear

      real dmin(lx1,ly1,lz1,lelt),emin(lx1,ly1,lz1,lelt) ! distf work arrays
      real xn(lx1,ly1,lz1,lelt),yn(lx1,ly1,lz1,lelt)
      real zn(lx1,ly1,lz1,lelt)

      real rho, mu
      integer j1_mod(6), j2_mod(6)
      real testa(lx1,ly1,lz1,lelv)
      real glmin, glmax
      integer iglsum
      real max_x,max_y,max_z
      integer n_pts_w,n_pts_lt1,n_pts_lt5

      if (wdopt.eq.1) then
        call cheap_dist(wd,1,'W  ')
      else if (wdopt.eq.2) then
        call distf(wd,1,'W  ',dmin,emin,xn,yn,zn)
      else
         if(nio.eq.0) then
           write(*,*)"Invalid wall dist option in sideset_yp...ABORT"
         endif
         call exit(0)
      endif

      nij = 2*ndim
      call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

c     NOTE: comp_sij premultiplies by 2 wrt the notation in Pope i.e. S_ij = 0.5*(du/dx + du/dx)

      do e=1,nelv
      do i=1,lxyz
         tauij(i,1,1,e,1) = sij(i,1,1,1,e)
         tauij(i,1,1,e,2) = sij(i,1,1,2,e)
         tauij(i,1,1,e,3) = sij(i,1,1,3,e)
         tauij(i,1,1,e,4) = sij(i,1,1,4,e)
         tauij(i,1,1,e,5) = sij(i,1,1,5,e)
         tauij(i,1,1,e,6) = sij(i,1,1,6,e)
      enddo
      enddo

c     performing a dssum

      if (dssopt.gt.0) then
       do i=1,6
        call dsavg(tauij(1,1,1,1,i))
       enddo
      endif

c     multiply with viscosity, assuming constant viscosity
      rho = cpfld(1,2)
      mu  = cpfld(1,1)

      call cmult(tauij,mu,lxyz*6*nelv)

      call rzero(yplus,lxyz*nelv)
      ypmin=1.0d30
      ypmax=-1.0d30
      call rzero(testa,lxyz*nelv)

c     modifiers for face indices - pick interior points behind face
c     for face =  1   2  3  4     5        6
      j1_mod = (/ 0, -1, 0, 1, lx1**2,-lx1**2/)
      j2_mod = (/ 1,  0,-1, 0,      0,      0/)

      n_pts_w = 0
      n_pts_lt1 = 0
      n_pts_lt5 = 0

      do e=1,nelv
       do f=1,2*ndim
         bid = boundaryID(f,e)
         do i=1,nss
          if(bid.eq.sslist(i)) then
c            if(nio.eq.0) write(*,*) "BOUNDARY FOUND:",bid
            call facind2(js1,jf1,jskip1,js2,jf2,jskip2,f)
            ia = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
             ia = ia + 1
             visc_st_x = tauij(j1,j2,1,e,1)*unx(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,4)*uny(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,6)*unz(ia,1,f,e)

             visc_st_y = tauij(j1,j2,1,e,4)*unx(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,2)*uny(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,5)*unz(ia,1,f,e)

             visc_st_z = tauij(j1,j2,1,e,6)*unx(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,5)*uny(ia,1,f,e)
     $                 + tauij(j1,j2,1,e,3)*unz(ia,1,f,e)

             normal_st = visc_st_x * unx(ia,1,f,e)
     $                 + visc_st_y * uny(ia,1,f,e)
     $                 + visc_st_z * unz(ia,1,f,e)

             shear_st_x = visc_st_x - normal_st * unx(ia,1,f,e)
             shear_st_y = visc_st_y - normal_st * uny(ia,1,f,e)
             shear_st_z = visc_st_z - normal_st * unz(ia,1,f,e)

             shear = sqrt(shear_st_x**2 + shear_st_y**2 +
     $                           shear_st_z**2)

             yplus(j1,j2,1,e) = wd(j1+j1_mod(f),j2+j2_mod(f),1,e)
     $                      *sqrt(shear*rho)/mu ! assumes const properties

             testa(j1+j1_mod(f),j2+j2_mod(f),1,e) = 42.0


             if(yplus(j1,j2,1,e).gt.ypmax) then
                ypmax = yplus(j1,j2,1,e)
                max_x = xm1(j1,j2,1,e)
                max_y = ym1(j1,j2,1,e)
                max_z = zm1(j1,j2,1,e)
             endif
             if(yplus(j1,j2,1,e).lt.ypmin) ypmin = yplus(j1,j2,1,e)

             n_pts_w = n_pts_w + 1
             if(yplus(j1,j2,1,e).lt.1.0) n_pts_lt1 = n_pts_lt1 + 1
             if(yplus(j1,j2,1,e).lt.5.0) n_pts_lt5 = n_pts_lt5 + 1

            enddo
            enddo
          endif
         enddo
       enddo
      enddo

      ypmin = glmin(ypmin,1)
      ypmax = glmax(ypmax,1)
      n_pts_w = iglsum(n_pts_w,1)
      n_pts_lt1 = iglsum(n_pts_lt1,1)
      n_pts_lt5 = iglsum(n_pts_lt5,1)

      if (nio.eq.0) then
       write(*,*) "MAX Y+ LOC:", max_x,max_y,max_z
       write(*,*) "Wall points:", n_pts_w
       write(*,*) "Y+ < 1:", n_pts_lt1
       write(*,*) "Y+ < 5:", n_pts_lt5
      endif

      call outpost(testa,vy,vz,pr,t,'bdt')

      return
      end
c-----------------------------------------------------------------------
      subroutine print_sideset_drag(sslist,nss,x0,scale)
      implicit none
c  sslist - sideset list, integer array of sidesets where we need traction
c  nss - number of sidesets i.e. the size of the sslist array
c  dssopt - direct stiffness sum opt, if dssopt>1, a dsavg on traction is performed

      include 'SIZE'
      include 'TOTAL'

      integer sslist(1), nss

c     for drag calculation
      real scale
      real drag
      real x0(3)

      integer iobj
      save iobj

c     drag calculation on internal obstacle

      if (istep.eq.0) then
        call create_obj(iobj,sslist,nss)
      endif

      call torque_calc(scale,x0,.false.,.false.) ! compute wall shear
      if(nio.eq.0) then
        write(*,*) "DRAGX: ", dragx(iobj)
        write(*,*) "DRAGY: ", dragy(iobj)
        write(*,*) "DRAGZ: ", dragz(iobj)
      endif

c     end of drag calculation

      return
      end
c-----------------------------------------------------------------------
      subroutine sideset_shear(shear,sslist,nss,dssopt)
      implicit none
c  sslist - sideset list, integer array of sidesets where we need traction
c  nss - number of sidesets i.e. the size of the sslist array
c  dssopt - direct stiffness sum opt, if dssopt>1, a dsavg on traction is performed

      include 'SIZE'
      include 'TOTAL'

      real shear(lx1,ly1,lz1,lelv)
      integer sslist(1), nss, dssopt

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)
      integer nij,bid,jf1,js1,jskip1,js2,jf2,jskip2
      integer e,f,i,j1,j2,ia

      real sij(lx1,ly1,lz1,6,lelv),
     $      tauij(lx1,ly1,lz1,lelv,6),
     $      ur(lxyz),us(lxyz),ut(lxyz),
     $      vr(lxyz),vs(lxyz),vt(lxyz),
     $      wr(lxyz),ws(lxyz),wt(lxyz)

      real visc_st_x, visc_st_y, visc_st_z, normal_st,
     $     shear_st_x, shear_st_y, shear_st_z

      nij = 2*ndim
      call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

c     NOTE: comp_sij premultiplies by 2 wrt the notation in Pope i.e. S_ij = 0.5*(du/dx + du/dx)

      do e=1,nelv
      do i=1,lxyz
         tauij(i,1,1,e,1) = sij(i,1,1,1,e)
         tauij(i,1,1,e,2) = sij(i,1,1,2,e)
         tauij(i,1,1,e,3) = sij(i,1,1,3,e)
         tauij(i,1,1,e,4) = sij(i,1,1,4,e)
         tauij(i,1,1,e,5) = sij(i,1,1,5,e)
         tauij(i,1,1,e,6) = sij(i,1,1,6,e)
      enddo
      enddo

c     performing a dssum

      if (dssopt.gt.0) then
       do i=1,6
        call dsavg(tauij(1,1,1,1,i))
       enddo
      endif

c     multiply with viscosity, assuming constant viscosity
      call cmult(tauij,cpfld(1,1),lxyz*6*nelv)

      call rzero(shear,lxyz*nelv)

      do e=1,nelv
       do f=1,2*ndim
         bid = boundaryID(f,e)
         do i=1,nss
          if(bid.eq.sslist(i)) then
            call facind2(js1,jf1,jskip1,js2,jf2,jskip2,f)
            ia = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
            ia = ia + 1
            visc_st_x = tauij(j1,j2,1,e,1)*unx(ia,1,f,e)
     $                + tauij(j1,j2,1,e,4)*uny(ia,1,f,e)
     $                + tauij(j1,j2,1,e,6)*unz(ia,1,f,e)

            visc_st_y = tauij(j1,j2,1,e,4)*unx(ia,1,f,e)
     $                + tauij(j1,j2,1,e,2)*uny(ia,1,f,e)
     $                + tauij(j1,j2,1,e,5)*unz(ia,1,f,e)

            visc_st_z = tauij(j1,j2,1,e,6)*unx(ia,1,f,e)
     $                + tauij(j1,j2,1,e,5)*uny(ia,1,f,e)
     $                + tauij(j1,j2,1,e,3)*unz(ia,1,f,e)

            normal_st = visc_st_x * unx(ia,1,f,e)
     $                + visc_st_y * uny(ia,1,f,e)
     $                + visc_st_z * unz(ia,1,f,e)

            shear_st_x = visc_st_x - normal_st * unx(ia,1,f,e)
            shear_st_y = visc_st_y - normal_st * uny(ia,1,f,e)
            shear_st_z = visc_st_z - normal_st * unz(ia,1,f,e)

            shear(j1,j2,1,e) = sqrt(shear_st_x**2 + shear_st_y**2 +
     $                          shear_st_z**2)

            enddo
            enddo
          endif
         enddo
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reorder_sij(tauij)
      implicit none
c  tr_x, tr_y, tr_z - traction components
c  sslist - sideset list, integer array of sidesets where we need traction
c  nss - number of sidesets i.e. the size of the sslist array
c  dssopt - direct stiffness sum opt, if dssopt>1, a dsavg on traction is performed

      include 'SIZE'
      include 'TOTAL'

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real sij(lx1,ly1,lz1,6,lelv),
     $      tauij(lx1,ly1,lz1,lelv,6),
     $      ur(lxyz),us(lxyz),ut(lxyz),
     $      vr(lxyz),vs(lxyz),vt(lxyz),
     $      wr(lxyz),ws(lxyz),wt(lxyz)

      integer e, f, i, nij

      nij = 6
      call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

      do e=1,nelv
      do i=1,lxyz
         tauij(i,1,1,e,1) = sij(i,1,1,1,e)
         tauij(i,1,1,e,2) = sij(i,1,1,2,e)
         tauij(i,1,1,e,3) = sij(i,1,1,3,e)
         tauij(i,1,1,e,4) = sij(i,1,1,4,e)
         tauij(i,1,1,e,5) = sij(i,1,1,5,e)
         tauij(i,1,1,e,6) = sij(i,1,1,6,e)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sideset_traction(tr_x,tr_y,tr_z,sslist,nss,dssopt)
      implicit none
c  tr_x, tr_y, tr_z - traction components
c  sslist - sideset list, integer array of sidesets where we need traction
c  nss - number of sidesets i.e. the size of the sslist array
c  dssopt - direct stiffness sum opt, if dssopt>1, a dsavg on traction is performed

      include 'SIZE'
      include 'TOTAL'

      integer sslist(1),nss,dssopt
      real tr_x(lx1,ly1,lz1,lelv), tr_y(lx1,ly1,lz1,lelv),
     $     tr_z(lx1,ly1,lz1,lelv)

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real tauij(lx1,ly1,lz1,6,lelv),
     $     tmp(lx1,ly1,lz1,lelv,6),
     $     ur(lxyz),us(lxyz),ut(lxyz),
     $     vr(lxyz),vs(lxyz),vt(lxyz),
     $     wr(lxyz),ws(lxyz),wt(lxyz)

      integer nij,bid,jf1,js1,jskip1,js2,jf2,jskip2
      integer e,f,i,j,j1,j2,ia,ntot
c      real ramp,ttol,fnshape,mdelta,areaf
c      real mask(lx1,ly1,lz1,lelv),ones(lx1,ly1,lz1,lelv)
      real drag_x, drag_y, drag_z, glsum

      nij = 6
      ntot = lxyz*nelv

      call rzero(tauij,lxyz*6*nelv)
      call comp_sij(tauij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

c     performing a dssum

      if (dssopt.gt.0) then

       do j=1,6
       do e=1,nelv ! reorder to match dsavg arg type
       do i=1,lxyz
          tmp(i,1,1,e,j) = tauij(i,1,1,j,e)
       enddo
       enddo
       enddo

       do i=1,6
        call dsavg(tmp(1,1,1,1,i))
       enddo

       do j=1,6
       do e=1,nelv ! reorder to return to unit stride for vectorization
       do i=1,lxyz
           tauij(i,1,1,j,e) = tmp(i,1,1,e,j)
       enddo
       enddo
       enddo

      endif

c     multiply with viscosity, assuming constant viscosity
      call cmult(tauij,-cpfld(1,1),lxyz*6*nelv) 

c     subtract pressure from diagonal terms
      do i=1,3
       do e = 1,nelv
        call add2(tauij(1,1,1,i,e),pr(1,1,1,e),lxyz)
       enddo
      enddo

c     set traction to zero except on all boundaries of interest

      call rzero(tr_x,ntot)
      call rzero(tr_y,ntot)
      call rzero(tr_z,ntot)

      drag_x = 0.0 
      drag_y = 0.0 
      drag_z = 0.0 

      do e=1,nelv
       do f=1,2*ndim
         bid = boundaryID(f,e)
         do i=1,nss
          if(bid.eq.sslist(i)) then
            call facind2(js1,jf1,jskip1,js2,jf2,jskip2,f)
            ia = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
            ia = ia + 1
            tr_x(j1,j2,1,e)=-( tauij(j1,j2,1,1,e)*unx(ia,1,f,e)
     $                        +tauij(j1,j2,1,4,e)*uny(ia,1,f,e)
     $                        +tauij(j1,j2,1,6,e)*unz(ia,1,f,e))

            tr_y(j1,j2,1,e)=-( tauij(j1,j2,1,4,e)*unx(ia,1,f,e)
     $                        +tauij(j1,j2,1,2,e)*uny(ia,1,f,e)
     $                        +tauij(j1,j2,1,5,e)*unz(ia,1,f,e))

            tr_z(j1,j2,1,e)=-( tauij(j1,j2,1,6,e)*unx(ia,1,f,e)
     $                        +tauij(j1,j2,1,5,e)*uny(ia,1,f,e)
     $                        +tauij(j1,j2,1,3,e)*unz(ia,1,f,e))

c     drag calculation to compare with torque_calc
            drag_x = drag_x + (
     $        tauij(j1,j2,1,1,e)*unx(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,4,e)*uny(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,6,e)*unz(ia,1,f,e)*area(ia,1,f,e) )

            drag_y = drag_y + (
     $        tauij(j1,j2,1,4,e)*unx(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,2,e)*uny(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,5,e)*unz(ia,1,f,e)*area(ia,1,f,e) )

            drag_z = drag_z + (
     $        tauij(j1,j2,1,6,e)*unx(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,5,e)*uny(ia,1,f,e)*area(ia,1,f,e)
     $      + tauij(j1,j2,1,3,e)*unz(ia,1,f,e)*area(ia,1,f,e) )

            enddo
            enddo
          endif
         enddo
       enddo
      enddo

      drag_x = glsum(drag_x,1)
      drag_y = glsum(drag_y,1)
      drag_z = glsum(drag_z,1)

      if (nio.eq.0) then
        write(*,71) drag_x  
        write(*,72) drag_y  
        write(*,73) drag_z  
      endif
      
71    format(1p1e20.13,': ansh drag x ')
72    format(1p1e20.13,': ansh drag y ')
73    format(1p1e20.13,': ansh drag z ')

c     create smoothing mask
c      ttol = 1.0E+06
c
c      call rzero(mask,ntot)
cc      call rone(ones,ntot)
c      fnshape = 0
c      mdelta = 0.01
c
c      do e=1,nelv
c      do f=1,6
c       bid = boundaryID(f,e)
c       do i = 1,nss
c        if(bid.eq.sslist(i)) then
cc         areaf = facint_v(ones,area,f,e)
c         call facind2(js1,jf1,jskip1,js2,jf2,jskip2,f)
c         ia = 0
c
c         do j2 = js2,jf2,jskip2
c         do j1 = js1,jf1,jskip1
c
c          ia = ia + 1
c
c          if(bid.eq.4.or.bid.eq.6) then  ! x normal faces
c            yy = ym1(j1,j2,1,e)
c            zz = zm1(j1,j2,1,e)
c            fnshape = sin(pi*yy/0.2)*sin(-pi*zz/0.2)
c          else if(bid.eq.5) then ! y normal face
c            xx = xm1(j1,j2,1,e)
c            zz = zm1(j1,j2,1,e)
c            fnshape = sin(pi*(xx-0.4)/0.1)*sin(-pi*zz/0.2)
c          else if(bid.eq.7) then ! z normal face
c            xx = xm1(j1,j2,1,e)
c            yy = ym1(j1,j2,1,e)
c            fnshape = sin(pi*(xx-0.4)/0.1)*sin(pi*yy/0.2)
c          endif
c
c          fnshape = fnshape + mdelta
c          mask(j1,j2,1,e) = fnshape  ! storing in case we need to visualise
c          tr_x(j1,j2,1,e) = tr_x(j1,j2,1,e) * fnshape
c          tr_y(j1,j2,1,e) = tr_y(j1,j2,1,e) * fnshape
c          tr_z(j1,j2,1,e) = tr_z(j1,j2,1,e) * fnshape
c         enddo
c         enddo
c        endif
c       enddo
c      enddo
c      enddo
c
c      call outpost(mask,vy,vz,pr,t,'msk')
c
c      ramp = 1.0
c      if(istep.le.1000) ramp = istep/1000.0
c
c      if (nio.eq.0) write(*,*) "RAMP FACTOR",ramp
c
c      call cmult(tr_x,ramp,ntot)
c      call cmult(tr_y,ramp,ntot)
c      call cmult(tr_z,ramp,ntot)

      return
      end
c-----------------------------------------------------------------------
