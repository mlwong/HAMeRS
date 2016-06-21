c
c  File:        $URL$
c  Package:     SAMRAI algorithms
c  Copyright:   (c) 1997-2014 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: F77 routines for updating 1d flux sums from fluxes.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2014 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 1d arrays in FORTRAN routines.
c
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
      subroutine upfluxsum1d(
     &  ilo0,ihi0,
     &  flxgc0,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ihi0,
     &  flxgc0,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0),
     &  fluxsum(1)
      integer ie0
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie0 = ilo0
      else
        ie0 = ihi0+1
      endif 
 
      fluxsum(1)=fluxsum(1)+flux(ie0)
c
      return
      end
