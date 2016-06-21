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
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
c









c
c***********************************************************************
c Sum node data on fine patch with outernode data on coarse patch
c and store back in node data on fine patch.
c i.e. fine_node_data += coarse_outernode_data
c***********************************************************************
c
      subroutine nodeouternodesum3d(
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,ngc2,
     &  fnode,
     &  couternodelower0,couternodeupper0,
     &  couternodelower1,couternodeupper1,
     &  couternodelower2,couternodeupper2)
c***********************************************************************
      implicit none

      integer
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2
      integer ratio(0:3-1),depth
      integer ngc0,ngc1,ngc2
      double precision
     &  fnode(filo0-ngc0:fihi0+1+ngc0,
     &          filo1-ngc1:fihi1+1+ngc1,
     &          filo2-ngc2:fihi2+1+ngc2,depth),
     &  couternodelower0(cilo1+1:cihi1,
     &          cilo2+1:cihi2,depth),
     &  couternodeupper0(cilo1+1:cihi1,
     &          cilo2+1:cihi2,depth),
     &  couternodelower1(cilo0:cihi0+1,
     &          cilo2+1:cihi2,depth),
     &  couternodeupper1(cilo0:cihi0+1,
     &          cilo2+1:cihi2,depth),
     &  couternodelower2(cilo0:cihi0+1,
     &          cilo1:cihi1+1,depth),
     &  couternodeupper2(cilo0:cihi0+1,
     &          cilo1:cihi1+1,depth)
      integer ic0,if0,ic1,if1,ic2,if2,id
      integer cilo0_loop,cilo1_loop,cilo2_loop,
     &        cihi0_loop,cihi1_loop,cihi2_loop
c
c***********************************************************************
c 

c sum along XLOWER face
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if0=cilo0*ratio(0)
      do id = 1, depth
         do ic2=cilo2_loop,cihi2_loop
            if2=ic2*ratio(2)
            do ic1=cilo1_loop,cihi1_loop
               if1=ic1*ratio(1)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodelower0(ic1,ic2,id)
            enddo
         enddo
      enddo

c sum along XUPPER face
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if0=(cihi0+1)*ratio(0)
      do id = 1, depth
         do ic2=cilo2_loop,cihi2_loop
            if2=ic2*ratio(2)
            do ic1=cilo1_loop,cihi1_loop
               if1=ic1*ratio(1)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodeupper0(ic1,ic2,id)
            enddo
         enddo
      enddo

c sum along YLOWER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if1=cilo1*ratio(1)
      do id = 1, depth
         do ic2=cilo2_loop,cihi2_loop
            if2=ic2*ratio(2)
            do ic0=cilo0_loop,cihi0_loop
               if0=ic0*ratio(0)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodelower1(ic0,ic2,id)
            enddo
         enddo
      enddo

c sum along YUPPER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if1=(cihi1+1)*ratio(1)
      do id = 1, depth
         do ic2=cilo2_loop,cihi2_loop
            if2=ic2*ratio(2)
            do ic0=cilo0_loop,cihi0_loop
               if0=ic0*ratio(0)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodeupper1(ic0,ic2,id)
            enddo
         enddo
      enddo

c sum along ZLOWER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo1_loop=cilo1
      cihi1_loop=cihi1+1

      if2=cilo2*ratio(2)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            do ic0=cilo0_loop,cihi0_loop
               if0=ic0*ratio(0)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodelower2(ic0,ic1,id)
            enddo
         enddo
      enddo

c sum along ZUPPER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo1_loop=cilo1
      cihi1_loop=cihi1+1

      if2=(cihi2+1)*ratio(2)
      do id = 1, depth
         do ic1=cilo1_loop,cihi1_loop
            if1=ic1*ratio(1)
            do ic0=cilo0_loop,cihi0_loop
               if0=ic0*ratio(0)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternodeupper2(ic0,ic1,id)
            enddo
         enddo
      enddo

c
      return
      end
c

c
c***********************************************************************
c Fill hanging nodes on fine patch (i.e. those nodes that do not overlap
c a coarse level node) along patch boundary by interpolation from
c appropriate neighboring coarse nodes.
c***********************************************************************
c
      subroutine nodehangnodeinterp3d(
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2,
     &  bboxilo0,bboxilo1,bboxilo2,
     &  bboxihi0,bboxihi1,bboxihi2,
     &  bboxloc,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,ngc2,
     &  fnode)
c***********************************************************************
      implicit none

      integer XLOWER,XUPPER,YLOWER,YUPPER,ZLOWER,ZUPPER
      parameter (XLOWER=0)
      parameter (XUPPER=1)
      parameter (YLOWER=2)
      parameter (YUPPER=3)
      parameter (ZLOWER=4)
      parameter (ZUPPER=5)

      double precision one
      parameter (one=1.0d0)

c***********************************************************************
      integer
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2
      integer
     &  bboxilo0,bboxilo1,bboxilo2,
     &  bboxihi0,bboxihi1,bboxihi2,
     &  bboxloc
      integer ratio(0:3-1),depth
      integer ngc0,ngc1,ngc2
      double precision
     &  fnode(filo0-ngc0:fihi0+1+ngc0,
     &          filo1-ngc1:fihi1+1+ngc1,
     &          filo2-ngc2:fihi2+1+ngc2,depth)
      integer ic0,if0,ir0,ic1,if1,ir1,ic2,if2,ir2,id
      integer iflower0,iflower1,iflower2,ifupper0,ifupper1,ifupper2
      integer cilo0_loop,cilo1_loop,cilo2_loop,
     &        cihi0_loop,cihi1_loop,cihi2_loop
      double precision dratio0,dratio1,dratio2,x0,x1,x2
c
c***********************************************************************
c

      dratio0 = dble(ratio(0))
      dratio1 = dble(ratio(1))
      dratio2 = dble(ratio(2))

      cilo0_loop=max(cilo0, bboxilo0)
      cihi0_loop=min(cihi0, bboxihi0)
      cilo1_loop=max(cilo1, bboxilo1)
      cihi1_loop=min(cihi1, bboxihi1)
      cilo2_loop=max(cilo2, bboxilo2)
      cihi2_loop=min(cihi2, bboxihi2)

      if (bboxloc.eq.XLOWER) then
      
         if0=cilo0*ratio(0)

         do id = 1, depth

c interpolate nodes along lines in y-direction
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic2=cilo2_loop,cihi2_loop+1
                  if2=ic2*ratio(2)
                  fnode(if0,iflower1+ir1,if2,id) =
     &               fnode(if0,ifupper1,if2,id) * x1
     &           +   fnode(if0,iflower1,if2,id) * (one-x1)
               enddo
            enddo
         enddo

c interpolate nodes along lines in z-direction
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic1=cilo1_loop,cihi1_loop+1
                  if1=ic1*ratio(1)
                  fnode(if0,if1,iflower2+ir2,id) =
     &               fnode(if0,if1,ifupper2,id) * x2
     &           +   fnode(if0,if1,iflower2,id) * (one-x2)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic1=cilo1_loop,cihi1_loop
                  iflower1=ic1*ratio(1)
                  ifupper1=iflower1+ratio(1)
                  do ir1=1,ratio(1)-1
                     x1 = dble(ir1)/dratio1
                     fnode(if0,iflower1+ir1,iflower2+ir2,id) =
     &               (fnode(if0,iflower1,iflower2,id) * (one-x1) +
     &                fnode(if0,ifupper1,iflower2,id) * x1) * (one-x2) +
     &               (fnode(if0,iflower1,ifupper2,id) * (one-x1) +
     &                fnode(if0,ifupper1,ifupper2,id) * x1) * x2
                  enddo
               enddo
            enddo
         enddo

         enddo

      else if (bboxloc.eq.XUPPER) then

         if0=(cihi0+1)*ratio(0)

         do id = 1, depth

c interpolate nodes along lines in y-direction
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic2=cilo2_loop,cihi2_loop+1
                  if2=ic2*ratio(2)
                  fnode(if0,iflower1+ir1,if2,id) =
     &               fnode(if0,ifupper1,if2,id) * x1
     &           +   fnode(if0,iflower1,if2,id) * (one-x1)
               enddo
            enddo
         enddo

c interpolate nodes along lines in z-direction
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic1=cilo1_loop,cihi1_loop+1
                  if1=ic1*ratio(1)
                  fnode(if0,if1,iflower2+ir2,id) =
     &               fnode(if0,if1,ifupper2,id) * x2
     &           +   fnode(if0,if1,iflower2,id) * (one-x2)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic1=cilo1_loop,cihi1_loop
                  iflower1=ic1*ratio(1)
                  ifupper1=iflower1+ratio(1)
                  do ir1=1,ratio(1)-1
                     x1 = dble(ir1)/dratio1
                     fnode(if0,iflower1+ir1,iflower2+ir2,id) =
     &               (fnode(if0,iflower1,iflower2,id) * (one-x1) +
     &                fnode(if0,ifupper1,iflower2,id) * x1) * (one-x2) +
     &               (fnode(if0,iflower1,ifupper2,id) * (one-x1) +
     &                fnode(if0,ifupper1,ifupper2,id) * x1) * x2
                  enddo
               enddo
            enddo
         enddo

         enddo

      else if (bboxloc.eq.YLOWER) then

         if1=cilo1*ratio(1)

         do id = 1, depth

c interpolate nodes along lines in x-direction
         do ic2=cilo2_loop,cihi2_loop+1
            if2=ic2*ratio(2)
            do ic0=cilo0_loop,cihi0_loop
               iflower0=ic0*ratio(0)
               ifupper0=iflower0+ratio(0)
               do ir0=1,ratio(0)-1
                  x0 = dble(ir0)/dratio0
                  fnode(iflower0+ir0,if1,if2,id) =
     &               fnode(ifupper0,if1,if2,id) * x0
     &           +   fnode(iflower0,if1,if2,id) * (one-x0)
               enddo
            enddo
         enddo

c interpolate nodes along lines in z-direction
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic0=cilo0_loop,cihi0_loop+1
                  if0=ic0*ratio(0)
                  fnode(if0,if1,iflower2+ir2,id) =
     &               fnode(if0,if1,ifupper2,id) * x2
     &           +   fnode(if0,if1,iflower2,id) * (one-x2)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic0=cilo0_loop,cihi0_loop
                  iflower0=ic0*ratio(0)
                  ifupper0=iflower0+ratio(0)
                  do ir0=1,ratio(0)-1
                     x0 = dble(ir0)/dratio0
                     fnode(iflower0+ir0,if1,iflower2+ir2,id) =
     &               (fnode(iflower0,if1,iflower2,id) * (one-x0) +
     &                fnode(ifupper0,if1,iflower2,id) * x0) * (one-x2) +
     &               (fnode(iflower0,if1,ifupper2,id) * (one-x0) +
     &                fnode(ifupper0,if1,ifupper2,id) * x0) * x2
                  enddo
               enddo
            enddo
         enddo

         enddo

      else if (bboxloc.eq.YUPPER) then

         if1=(cihi1+1)*ratio(1)
   
         do id = 1, depth

c interpolate nodes along lines in x-direction
         do ic2=cilo2_loop,cihi2_loop+1
            if2=ic2*ratio(2)
            do ic0=cilo0_loop,cihi0_loop
               iflower0=ic0*ratio(0)
               ifupper0=iflower0+ratio(0)
               do ir0=1,ratio(0)-1
                  x0 = dble(ir0)/dratio0
                  fnode(iflower0+ir0,if1,if2,id) =
     &               fnode(ifupper0,if1,if2,id) * x0
     &           +   fnode(iflower0,if1,if2,id) * (one-x0)
               enddo
            enddo
         enddo

c interpolate nodes along lines in z-direction
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic0=cilo0_loop,cihi0_loop+1
                  if0=ic0*ratio(0)
                  fnode(if0,if1,iflower2+ir2,id) =
     &               fnode(if0,if1,ifupper2,id) * x2
     &           +   fnode(if0,if1,iflower2,id) * (one-x2)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic0=cilo0_loop,cihi0_loop
                  iflower0=ic0*ratio(0)
                  ifupper0=iflower0+ratio(0)
                  do ir0=1,ratio(0)-1
                     x0 = dble(ir0)/dratio0
                     fnode(iflower0+ir0,if1,iflower2+ir2,id) =
     &               (fnode(iflower0,if1,iflower2,id) * (one-x0) +
     &                fnode(ifupper0,if1,iflower2,id) * x0) * (one-x2) +
     &               (fnode(iflower0,if1,ifupper2,id) * (one-x0) +
     &                fnode(ifupper0,if1,ifupper2,id) * x0) * x2
                  enddo
               enddo
            enddo
         enddo

         enddo

      else if (bboxloc.eq.ZLOWER) then

         if2=cilo2*ratio(2)

         do id = 1, depth

c interpolate nodes along lines in x-direction
         do ic1=cilo1_loop,cihi1_loop+1
            if1=ic1*ratio(1)
            do ic0=cilo0_loop,cihi0_loop
               iflower0=ic0*ratio(0)
               ifupper0=iflower0+ratio(0)
               do ir0=1,ratio(0)-1
                  x0 = dble(ir0)/dratio0
                  fnode(iflower0+ir0,if1,if2,id) =
     &               fnode(ifupper0,if1,if2,id) * x0
     &           +   fnode(iflower0,if1,if2,id) * (one-x0)
               enddo
            enddo
         enddo

c interpolate nodes along lines in y-direction
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic0=cilo0_loop,cihi0_loop+1
                  if0=ic0*ratio(0)
                  fnode(if0,iflower1+ir1,if2,id) =
     &               fnode(if0,ifupper1,if2,id) * x1
     &           +   fnode(if0,iflower1,if2,id) * (one-x1)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic0=cilo0_loop,cihi0_loop
                  iflower0=ic0*ratio(0)
                  ifupper0=iflower0+ratio(0)
                  do ir0=1,ratio(0)-1
                     x0 = dble(ir0)/dratio0
                     fnode(iflower0+ir0,iflower1+ir1,if2,id) =
     &               (fnode(iflower0,iflower1,if2,id) * (one-x0) +
     &                fnode(ifupper0,iflower1,if2,id) * x0) * (one-x1) +
     &               (fnode(iflower0,ifupper1,if2,id) * (one-x0) +
     &                fnode(ifupper0,ifupper1,if2,id) * x0) * x1
                  enddo
               enddo
            enddo
         enddo

         enddo

      else if (bboxloc.eq.ZUPPER) then

         if2=(cihi2+1)*ratio(2)

         do id = 1, depth

c interpolate nodes along lines in x-direction
         do ic1=cilo1_loop,cihi1_loop+1
            if1=ic1*ratio(1)
            do ic0=cilo0_loop,cihi0_loop
               iflower0=ic0*ratio(0)
               ifupper0=iflower0+ratio(0)
               do ir0=1,ratio(0)-1
                  x0 = dble(ir0)/dratio0
                  fnode(iflower0+ir0,if1,if2,id) =
     &               fnode(ifupper0,if1,if2,id) * x0
     &           +   fnode(iflower0,if1,if2,id) * (one-x0)
               enddo
            enddo
         enddo

c interpolate nodes along lines in y-direction
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic0=cilo0_loop,cihi0_loop+1
                  if0=ic0*ratio(0)
                  fnode(if0,iflower1+ir1,if2,id) =
     &               fnode(if0,ifupper1,if2,id) * x1
     &           +   fnode(if0,iflower1,if2,id) * (one-x1)
               enddo
            enddo
         enddo

c interpolate face nodes in coarse cell interiors
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic0=cilo0_loop,cihi0_loop
                  iflower0=ic0*ratio(0)
                  ifupper0=iflower0+ratio(0)
                  do ir0=1,ratio(0)-1
                     x0 = dble(ir0)/dratio0
                     fnode(iflower0+ir0,iflower1+ir1,if2,id) =
     &               (fnode(iflower0,iflower1,if2,id) * (one-x0) +
     &                fnode(ifupper0,iflower1,if2,id) * x0) * (one-x1) +
     &               (fnode(iflower0,ifupper1,if2,id) * (one-x0) +
     &                fnode(ifupper0,ifupper1,if2,id) * x0) * x1
                  enddo
               enddo
            enddo
         enddo

         enddo

      endif
c
      return
      end
c

