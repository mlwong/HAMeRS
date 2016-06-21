c
c  File:        $URL$
c  Package:     SAMRAI mesh
c  Copyright:   (c) 1997-2014 Lawrence Livermore National Security, LLC
c  Release:     
c  Revision:    
c  Modified:    
c  Description:    F77 routines for coarsening 2d integer tag values.
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
c Constant averaging for 2d cell-centered double data.   The operation
c uses the rule that if any value on the fine mesh contained in a cell of 
c the coarse mesh is not equal to zero, the value on the coarse mesh is 
c set to one.
c***********************************************************************
c
      subroutine coarsentags2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      integer zero, one
      parameter (zero=0)
      parameter (one=1)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      integer
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  if (arrayf(if0,if1) .ne. zero)
     &               arrayc(ic0,ic1)=one
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
