c
c  File:        $URL$
c  Package:     SAMRAI algs
c  Copyright:   (c) 1997-2014 Lawrence Livermore National Security, LLC
c  Revision:    
c  Modified:    
c  Description:    F77 routines for summing outernode data with
c                  other node or outernode data
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2014 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c



c
c***********************************************************************
c Sum node data on fine patch with outernode data on coarse patch
c and store back in node data on fine patch.
c i.e. fine_node_data += coarse_outernode_data
c***********************************************************************
c
      subroutine nodeouternodesum2d(
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,
     &  fnode,
     &  couternodelower0,couternodeupper0,
     &  couternodelower1,couternodeupper1)
c***********************************************************************
      implicit none

      integer
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1
      integer ratio(0:2-1), depth
      integer ngc0,ngc1
      double precision
     &  fnode(filo0-ngc0:fihi0+1+ngc0,
     &          filo1-ngc1:fihi1+1+ngc1,depth),
     &  couternodelower0(cilo1+1:cihi1,depth),
     &  couternodeupper0(cilo1+1:cihi1,depth),
     &  couternodelower1(cilo0:cihi0+1,depth),
     &  couternodeupper1(cilo0:cihi0+1,depth)
      integer ic0,if0,ic1,if1,id
      integer cilo0_loop,cilo1_loop,cihi0_loop,cihi1_loop
c
c***********************************************************************
c

c sum along XLOWER side
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1

      if0=cilo0*ratio(0)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            fnode(if0,if1,id) = fnode(if0,if1,id) + 
     &                          couternodelower0(ic1,id)
         enddo
      enddo

c sum along XUPPER side
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1

      if0=(cihi0+1)*ratio(0)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            fnode(if0,if1,id) = fnode(if0,if1,id) + 
     &                          couternodeupper0(ic1,id)
         enddo
      enddo

c sum along YLOWER side
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1

      if1=cilo1*ratio(1)
      do id = 1, depth
         do ic0=cilo0_loop,cihi0_loop
            if0=ic0*ratio(0)
            fnode(if0,if1,id) = fnode(if0,if1,id) + 
     &                          couternodelower1(ic0,id)
         enddo
      enddo

c sum along YUPPER side
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1

      if1=(cihi1+1)*ratio(1)
      do id = 1, depth
         do ic0=cilo0_loop,cihi0_loop
            if0=ic0*ratio(0)
            fnode(if0,if1,id) = fnode(if0,if1,id) + 
     &                          couternodeupper1(ic0,id)
         enddo
      enddo

c
      return
      end

c
c***********************************************************************
c Fill hanging nodes on fine patch (i.e. those nodes that do not overlap
c a coarse level node) along patch boundary by interpolation from 
c appropriate neighboring coarse nodes
c***********************************************************************
c
      subroutine nodehangnodeinterp2d(
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  bboxilo0,bboxilo1,bboxihi0,bboxihi1,
     &  bboxloc,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,
     &  fnode)
c***********************************************************************
      implicit none

      integer XLOWER,XUPPER,YLOWER,YUPPER
      parameter (XLOWER=0)
      parameter (XUPPER=1)
      parameter (YLOWER=2)
      parameter (YUPPER=3)

      double precision one
      parameter (one=1.0d0)
c***********************************************************************
c
      integer
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1
      integer bboxilo0,bboxilo1,bboxihi0,bboxihi1,bboxloc
      integer ratio(0:2-1),depth
      integer ngc0,ngc1
      double precision
     &  fnode(filo0-ngc0:fihi0+1+ngc0,
     &          filo1-ngc1:fihi1+1+ngc1,depth)
      integer ic0,if0,ir0,ic1,if1,ir1,id
      integer cilo0_loop,cilo1_loop,cihi0_loop,cihi1_loop
      double precision dratio0,dratio1,x0,x1
c
c***********************************************************************
c
      dratio0 = dble(ratio(0))
      dratio1 = dble(ratio(1))

      cilo0_loop=max(cilo0, bboxilo0)
      cihi0_loop=min(cihi0, bboxihi0)
      cilo1_loop=max(cilo1, bboxilo1)
      cihi1_loop=min(cihi1, bboxihi1)
c
c***********************************************************************
c 
      if (bboxloc.eq.XLOWER) then

         if0=cilo0*ratio(0)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               fnode(if0,if1+ir1,id) =
     &            fnode(if0,if1,id) * (one - x1) +
     &            fnode(if0,if1+ratio(1),id) * x1
            enddo
         enddo
      enddo

      else if (bboxloc.eq.XUPPER) then
 
         if0=(cihi0+1)*ratio(0)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               fnode(if0,if1+ir1,id) =
     &            fnode(if0,if1,id) * (one - x1) +
     &            fnode(if0,if1+ratio(1),id) * x1
            enddo
         enddo
      enddo

      else if (bboxloc.eq.YLOWER) then
      
         if1=cilo1*ratio(1)
      do id = 1, depth
         do ic0=cilo0_loop,cihi0_loop
            if0=ic0*ratio(0)
            do ir0=1,ratio(0)-1
               x0 = dble(ir0)/dratio0
               fnode(if0+ir0,if1,id) =
     &            fnode(if0,if1,id) * (one - x0) +
     &            fnode(if0+ratio(0),if1,id) * x0
            enddo
         enddo
      enddo

      else if (bboxloc.eq.YUPPER) then

         if1=(cihi1+1)*ratio(1)
      do id = 1, depth
         do ic0=cilo0_loop,cihi0_loop
            if0=ic0*ratio(0)
            do ir0=1,ratio(0)-1
               x0 = dble(ir0)/dratio0
               fnode(if0+ir0,if1,id) =
     &            fnode(if0,if1,id) * (one - x0) +
     &            fnode(if0+ratio(0),if1,id) * x0
            enddo
         enddo
      enddo
  
      endif
c
      return
      end
