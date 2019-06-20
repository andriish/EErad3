C=============================================================================
C---  basis functions
C=============================================================================
c---  Li2

      double complex  function cli2(z)
      implicit none
      double complex ris, z, bsli2_inside,bsli2_outside, wcli2
      double complex zlocal
      double precision zabs, pi, zeta2, border, tiny, arg

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0

      border = 0.3d0 
      tiny = 1d-14
      zabs = abs(z)
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli2(1d0/z)-zeta2-0.5d0*log(-z)**2
      elseif (zabs.le.border) then 
         ris=bsli2_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli2_outside(zlocal)
      endif

      cli2=ris
      return
      end
      
c---  recursion
      
      double complex  function wcli2(z)
      implicit none
      double complex z, cli2
      wcli2 =  cli2(z)
      return
      end

c--- Li3

      double complex  function cli3(z)
      implicit none
      double complex ris, z, bsli3_inside,bsli3_outside, wcli3
      double complex zlocal
      double precision zabs,border, pi, zeta2, zeta3,tiny,arg
      
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
     
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=wcli3(1d0/z)-log(-z)**3/6d0-zeta2*log(-z)
      elseif (zabs.le.border) then 
         ris=bsli3_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli3_outside(zlocal)
      endif

      cli3=ris
      return
      end
      
c---  recursion

      double complex  function wcli3(z)
      implicit none
      double complex z, cli3
      wcli3 =  cli3(z)
      return
      end

c--- Li4

      double complex  function cli4(z)
      implicit none
      double complex ris, z, bsli4_outside, bsli4_inside, wcli4
      double complex zlocal
      double precision zabs, pi, zeta2, zeta3, zeta4, border,tiny,arg
     
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli4(1d0/z) -log(-z)**4/24d0 - 7d0*zeta4/4d0 
     &           - zeta2*log(-z)**2/2d0
      elseif (zabs.le.border) then 
         ris=bsli4_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli4_outside(zlocal)
      endif

      
      cli4=ris
      return
      end
      
c     --- recursion for li4
      
      double complex  function wcli4(z)
      implicit none
      double complex z, cli4
      wcli4 = cli4(z)
      return
      end
c --- the case Li4(1-z^2) needs some special treatment because of its branch cut structure
c --- (that's what 'sbc' stands for: special branch cut)

c      double complex function cli4_sbc(z)
c      implicit none
c      double complex ris, z, cli4, myi,basis14
c      double complex ll1,ll2,ll3
c      double precision pi,zabs,zreal
c      integer s
      
c      pi=3.1415926535897932385D0
c      zabs = abs(z)
c      zreal = dreal(z)
c      myi = dcmplx(0d0,1d0)
           
c      if (zabs.le.1d0) then !normal li4
c         if (zreal.gt.0d0) then
c            ris = cli4(1d0 - z**2)
c         else if (zreal.eq.0d0 .and. s(z).eq.1) then !also normal li4
c            ris = cli4(1d0 - z**2 - dcmplx(0d0,1d-60))
c         else                   ! special branch cut configuration
c            ris = cli4(1d0 - z**2- dcmplx(0d0,1d-60))
c     &           - myi*pi*s(z)/3d0*(log(1d0 - z)+log(1d0+z))**3 
c         endif
c      else 
c         ll1=log(1d0/z)
c         ll2=log(1d0 - 1d0/z)
c         ll3=log(1d0 + 1d0/z)
c         ris = -2d0/3d0*ll1**4 + 4d0/3d0*ll2
c     &        *ll1**3 + 4d0/3d0*ll3*ll1**3 
c     &        - ll2**2*ll1**2 - ll3**2
c     &        *ll1**2 - 2d0*ll2*ll3
c     &        *ll1**2 - pi**2/3d0*ll1**2 + 1d0/3d0
c     &        *ll2**3*ll1 + 1d0/3d0
c     &        *ll3**3*ll1 + ll2
c     &        *ll3**2*ll1 + pi**2/3d0*ll1
c     &        *ll2 + ll3*ll2**2
c     &        *ll1 + pi**2/3d0*ll1*ll3 
c     &        - 1d0/24d0*ll2**4 - 1d0/24d0
c     &        *ll3**4 - 1d0/6d0*ll2
c     &        *ll3**3 - pi**2/12d0*ll2**2 
c     &        - pi**2/12d0*ll3**2 
c     &        - 1d0/4d0*ll2**2*ll3**2 
c     &        - 1d0/6d0*ll2**3*ll3 
c     &        - pi**2/6d0*ll2*ll3 
c     &        - 7*pi**4/360d0 
c     &        -basis14(1d0/z)
c      endif
      
c      cli4_sbc = ris
c      return
c      end

cc --- the case Li4(4z/(1+z)^2) also needs some special treatment because of its branch cut structure

c      double complex function cli4_sbc_2(z)
c      implicit none
c      double complex ris, z, cli4, myi, wcli4_sbc_2
c      double complex arg,llx,ll1px,wcli4sbc2,zlocal
c      double precision pi,zabs,zreal,ll2,tiny,arg2
c      integer s
      
c      pi=3.1415926535897932385D0
c      ll2 = dlog(2d0)
c      zabs = abs(z)
c      zreal = dreal(z)
c      myi = dcmplx(0d0,1d0)
c      llx = 1
c      ll1px = 1
c      wcli4sbc2 = 1
c      tiny = 1d-14
c      zlocal=z
     
c      ris = dcmplx(0d0,0d0)

c      if (zabs.lt.1d0) then
c         ris = cli4(4d0*z/(1d0+z)**2)
c      elseif (zabs.lt.1d0+tiny) then 
c         arg2=atan2(dimag(zlocal),dreal(zlocal))
c         zlocal=dcmplx(cos(arg2),sin(arg2))
c         arg = dcmplx(dreal(4d0*zlocal/(1d0+zlocal)**2),s(zlocal)*1d-60)
c         ris = cli4(arg)
c      else    
c         wcli4sbc2 =  wcli4_sbc_2(1d0/z)
c         llx = log(1d0/z)    
c         ll1px = log(1d0+1d0/z)          
c         ris = wcli4sbc2 + myi*pi*s(z)*
c     &(4d0*ll2**2*llx - 8d0*ll2**2* ll1px 
c     &+ 2d0*ll2*llx**2 - 8d0*ll2*llx
c     &* ll1px + 8d0*ll2*ll1px**2 
c     &+ 1d0/3d0*llx**3 - 2d0* ll1px*llx**2 
c     &+ 4d0* ll1px**2*llx - 8d0/3d0* ll1px**3 
c     &+ 8d0/3d0*ll2**3)


cc             ris = wcli4_sbc_2(1d0/z) + myi*pi*s(z)*
cc     &(4d0*ll2**2*log(1d0/z) - 8d0*ll2**2*log(1d0+1d0/z) 
cc     &+ 2d0*ll2*log(1d0/z)**2 - 8d0*ll2*log(1d0/z)
cc     &*log(1d0+1d0/z) + 8d0*ll2*log(1d0+1d0/z)**2 
cc     &+ 1d0/3d0*log(1d0/z)**3 - 2d0*log(1d0+1d0/z)*log(1d0/z)**2 
cc     &+ 4d0*log(1d0+1d0/z)**2*log(1d0/z) - 8d0/3d0*log(1d0+1d0/z)**3 
cc     &+ 8d0/3d0*ll2**3)
c      endif
      
c      cli4_sbc_2 = ris
c      return
c      end

cc     --- recursion for cli4_sbc_2
      
c      double complex  function wcli4_sbc_2(z)
c      implicit none
c      double complex z, cli4_sbc_2
c      wcli4_sbc_2 =  cli4_sbc_2(z)
c      return
c      end

cC-----------------------------------------------------------------------
cC     mapping of H_2-2(z) into convergent region
      
c      double complex  function ch2m2(z)
c      implicit none
c      double complex ris,z,bsh2m2_inside,bsh2m2_outside,cli4,cli2
c      double complex HPL4,wch2m2,myi,zlocal !,cli4_sbc,cli3
c      double precision pi,zeta2,zeta3,zeta4,zabs,zreal,border,tiny,arg
c      integer s

c      pi=3.1415926535897932385D0
c      zeta2=pi**2/6d0
c      zeta3=1.20205690315959428539973816151d0
c      zeta4=pi**4/90d0
c      myi = dcmplx(0d0,1d0)
c      tiny=1d-14
      
c      border = 0.3d0
c      zabs = abs(z)
c      zreal = dreal(z)
c      zlocal=z


c      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
c         ris = bsh2m2_inside(z)            
c      elseif (zabs.lt.1d0+tiny) then
c         if (zabs.gt.1d0) then
c            arg=atan2(dimag(zlocal),dreal(zlocal))
c            zlocal=dcmplx(cos(arg),sin(arg))
c         endif
c         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
c            ris = bsh2m2_outside(zlocal)
c         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
c            ris = HPL4(0,-1,0,1,-zlocal) 
c         endif
c      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
c         ris = dcmplx(0d0,0d0) +
c     &        wch2m2(1d0/z) 
c     &        + 37d0*pi**4/720d0 
c     &        - HPL4(0,1,0,0,1d0/z) 
c     &        - log(1d0/z)**4/24d0 
c     &        - pi**2/12d0*log(1d0/z)**2 
c     &        - pi**2/6d0*cli2(1d0/z) 
c     &        - cli4(-1d0/z) 
c     &        + 3d0*zeta3*log(1d0/z)/2d0 
c     &        - pi**3*myi*s(z)*log(1d0/z)/12d0
c      endif
      
c      ch2m2=ris
c      return
c      end
      
      
cc     --- recursion for H_2-2(z)
      
c      double complex  function wch2m2(z)
c      implicit none
c      double complex z, ch2m2
           
c      wch2m2 = ch2m2(z)
c      return
c      end

cC------------------------------------------------------------------------------
cC     mapping of H21-1(z) into convergent region 
      
c      double complex  function ch21m1(z)
c      implicit none
c      double complex ris,z,bsh21m1_inside,bsh21m1_outside_1,zlocal
c      double complex bsh21m1_outside_2,cli4,cli2,HPL4,wch21m1,myi,ch2m2
c      double precision pi,zeta2,zeta3,zeta4,border,zreal,zabs,ll2,tiny
c      double precision arg
c      integer s

c      pi=3.1415926535897932385D0
c      zeta2=pi**2/6d0
c      zeta3=1.20205690315959428539973816151d0
c      zeta4=pi**4/90d0
c      ll2 = dlog(2d0)
c      border = 0.3d0
c      myi = dcmplx(0d0,1d0)
c      tiny=1d-14

c      zabs = abs(z)
c      zreal = dreal(z)
c      zlocal=z


c      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
c         ris = bsh21m1_inside(z)           
c      elseif (zabs.lt.1d0+tiny) then
c         if (zabs.gt.1d0) then
c            arg=atan2(dimag(zlocal),dreal(zlocal))
c            zlocal=dcmplx(cos(arg),sin(arg))
c         endif
c         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
c            ris = bsh21m1_outside_1(zlocal)
c         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
c            ris = bsh21m1_outside_2(zlocal)
c         endif
c      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
c         ris = -wch21m1(1d0/z)-pi**4/144d0 -ch2m2(1d0/z) 
c     &        + log(1d0/z)**4/24d0 
c     &        + pi**2*ll2**2/3d0 - ll2**4/12d0 
c     &        + 3d0*pi**2*ll2*log(1d0/z)/4d0 
c     &        + pi**2*log(1d0/z)**2/8d0 
c     &        + pi**2*cli2(1d0/z)/4d0 
c     &        - 2*cli4(dcmplx(0.5d0,0d0)) 
c     &        + cli4(-1d0/z) 
c     &        - 7d0*zeta3*log(1d0/z)/8d0 
c     &        + myi*s(z)*(pi**3*ll2/6d0 
c     &        + pi**3*log(1d0/z)/12d0 
c     &        - 0.5d0*pi*ll2**2*log(1d0/z) 
c     &        - 0.5d0*pi*ll2*log(1d0/z)**2 
c     &        - pi*ll2*cli2(1d0/z))
c     &        - HPL4(0,0,1,-1,1d0/z) 
c     &        + HPL4(0,0,1,0,1d0/z) 
c     &        + HPL4(0,1,0,0,1d0/z) 
c     &        + HPL4(0,1,1,0,1d0/z)
c      endif

c      ch21m1=ris
c      return
c      end
      

c     --- recursion for H21-1(z)
      
c      double complex  function wch21m1(z)
c      implicit none
c      double complex z, ch21m1
            
c      wch21m1 =  ch21m1(z)
c      return
c      end
