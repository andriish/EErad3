      PROGRAM EECANALYTIC
      IMPLICIT NONE
      DOUBLE PRECISION GETLO,GETNLO,al,pi
      DOUBLE PRECISION g11,g12,g21,g22,g23,g24
      DOUBLE PRECISION g25,g26,g31,g32,g33,g34,g35
      DOUBLE PRECISION A, B
      DOUBLE PRECISION bnf,bnlc,blc
      integer i,N
      parameter(pi=3.141592653589793238d0)
      al=0.1180d0/(2.0d0*pi)      

      open(11,file='doc/predictions/analytic/LO.dat')
      do i=1,180
      N=50
      write(11,'(f8.4, E15.6)')(1.0d0*i-0.5),
     .al*GETLO(1.0d0*(1.0d0*i-1.0d0),1.0d0*i,N)      
      enddo
        
      close(11)
      open(12,file='doc/predictions/analytic/NLO.dat')
      do i=1,180
      N=50
      if (i.le.1.or.i.ge.179) then
C There are numeric problems in the first and last bins
      write(12,'(f8.4, E15.6)')(1.0d0*i-0.5),1.0d0
C      N=1
C      write(12,*)(1.0d0*i-0.5),al*GETLO(1.0d0*(1.0d0*i-1.0d0),1.0d0*i,N)
C     .+al*al*GETNLO(1.0d0*(1.0d0*i-1.0d0),1.0d0*i,N)     
      else
      write(12,'(f8.4, E15.6)')(1.0d0*i-0.5),
     .al*GETLO(1.0d0*(1.0d0*i-1.0d0),1.0d0*i,N)+
     .al*al*GETNLO(1.0d0*(1.0d0*i-1.0d0),1.0d0*i,N)     
      end if
      enddo
        
      close(12)
      write(*,*)'Run test at z=0.5'
      write(*,*) 'g11',g11(0.5d0)
      write(*,*) 'g12',g12(0.5d0)
      write(*,*) 'g21',g21(0.5d0)
      write(*,*) 'g22',g22(0.5d0)
      write(*,*) 'g23',g23(0.5d0)
      write(*,*) 'g24',g24(0.5d0)

      write(*,*) 'g31',g31(0.5d0)
      write(*,*) 'g32',g32(0.5d0)
      write(*,*) 'g33',g33(0.5d0)
      write(*,*) 'g34',g34(0.5d0)
      
      write(*,*) 'A', A(0.5d0)
      write(*,*) 'B', B(0.5d0)
      write(*,*) 'Bnf', BNF(0.5d0)
      write(*,*) 'Bnlc', BNlc(0.5d0)
      write(*,*) 'Blc', Blc(0.5d0)
    
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Simple integration of analytic functions in the given interval
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION GETLO(th1,th2,N)
      IMPLICIT NONE
      DOUBLE PRECISION z,th1, th2, pi,z1,z2,dz,sum,A
      integer i,N
      parameter(pi=3.141592653589793238d0)
      z1=(1.0d0-dcos(th1/180.0d0*pi))/2.0d0
      z2=(1.0d0-dcos(th2/180.0d0*pi))/2.0d0
      dz=(Z2-z1)/(1.0d0*N+1.0d0)  
      
        sum=0.0d0
        do i=1,N
        z=z1+dz*i
        sum=sum+A(z)/2.0d0*pi*dsqrt(z*(1.0D0-z))
        enddo
        GETLO=SUM/(1.0d0*N)
      end


      DOUBLE PRECISION FUNCTION GETNLO(th1,th2,N)
      IMPLICIT NONE
      DOUBLE PRECISION  z, th1, th2, pi,z1,z2,dz,sum,B
      integer i,N
      parameter(pi=3.141592653589793238d0)
      z1=(1.0d0-dcos(th1/180.0d0*pi))/2.0d0
      z2=(1.0d0-dcos(th2/180.0d0*pi))/2.0d0
      dz=(Z2-z1)/(1.0d0*N+1.0d0)  
        sum=0.0d0
        do i=1,50
        z=z1+dz*i
        sum=sum+B(z)/2.0d0*pi*dsqrt(z*(1.0D0-z))
        enddo
        GETNLO=SUM/(1.0d0*N)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The functions bellow were copied from Dixon's paper
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
      DOUBLE PRECISION FUNCTION A(Z)
      IMPLICIT NONE
      DOUBLE PRECISION nc, cf,ca,nf,tf, z
      nc=3.0d0
      cf=(nc**2.0d0-1.0d0)/(2.0d0*nc)
      ca=nc
      nf=5.0d0
      tf=1.0d0/2.0d0
      
      A=(cf*(-3.0d0+2.0d0*z)*(3.0d0*(2.0d0-3.0d0*z)*z+
     .2.0d0*(3.0d0-6.0d0*z+2.0d0*z**2)*DLog(1.0d0-z)))/
     .(4.0d0*(-1.0d0+z)*z**5)
      end

      DOUBLE PRECISION FUNCTION BLC(Z)
      IMPLICIT NONE
      DOUBLE PRECISION nc, cf,ca,nf,tf, z
            DOUBLE PRECISION g11,g12,g21,g22,g23,g24
      DOUBLE PRECISION g25,g26,g31,g32,g33,g34,g35
      nc=3.0d0
      cf=(nc**2-1.0d0)/(2.0d0*nc)
      ca=nc
      nf=5.0d0
      tf=1.0d0/2.0d0
      
      Blc=(63298-143577*z+72305*z**2+2064*z**3-31000*z**4+157060*z**5-
     .244800*z**6+122400*z**7)/(1440*(1-z)*z**4)-((3007-9329*z
     .+11309*z**2-6201*z**3+2716*z**4-48122*z**5+283140*z**6-667280*z**7
     .+673200*z**8-244800*z**9)*g11(z))/(720*(1-z)*z**5)-((19938
     .-38295*z+17261*z**2-336*z**3+13052*z**4-126900*z**5+422480*z**6
     .-550800*z**7+244800*z**8)*g12(z))/(720*(1-z)*z**4)+((87-211*z
     .+296*z**2-96*z**3+25*z**4-17*z**5+10*z**6+4*z**7)*g21(z))/
     .(24*(1-z)*z**5)+((3323-4726*z+1126*z**2-160*z**3-320*z**4+
     .4040*z**5-28480*z**6+61200*z**7-40800*z**8)*g22(z))/(120*z**5)-
     .((1-11*z)*g23(z))/(48*z**(7.0/2.0))-((4193-10159*z+8812*z**2-
     .2246*z**3+160*z**4+60*z**5+120*z**6)*g24(z))/(120*(1-z)*z**5)-
     .2*(3-31*z+116*z**2-170*z**3+85*z**4)*g31(z)+((5-21*z+18*z**2-
     .4*z**3)*g32(z))/(6*(1-z)*z**5)+((1+z**2)*g33(z))/(12*(1-z))
      end
      
      
      DOUBLE PRECISION FUNCTION BnLc(Z)
      IMPLICIT NONE
      DOUBLE PRECISION nc, cf,ca,nf,tf,z
      DOUBLE PRECISION g11,g12,g21,g22,g23,g24
      DOUBLE PRECISION g25,g26,g31,g32,g33,g34,g35
      nc=3.0d0
      cf=(nc**2-1.0d0)/(2.0d0*nc)
      ca=nc
      nf=5.0d0
      tf=1.0d0/2.0d0
      
      Bnlc=(9320-27552*z+14966*z**2+902*z**3-17359*z**4+75748*z**5-
     .115200*z**6+57600*z**7)/(720*(1-z)*z**4)-((4880-12412*z+
     .11322*z**2-3571*z**3+3225*z**4-31035*z**5+147846*z**6-321680*z**7
     .+316800*z**8-115200*z**9)*g11(z))/(360*(1-z)*z**5)-((11424-
     .25029*z+10971*z**2-742*z**3+18696*z**4-138600*z**5+412960*z**6-
     .518400*z**7+230400*z**8)*g12(z))/(720*(1-z)*z**4)+((314-760*z+
     .721*z**2-140*z**3+15*z**4-184*z**5+235*z**6-91*z**7)*g21(z))/
     .(120*(1-z)*z**5)+((952-1431*z+315*z**2-40*z**3-340*z**4+2660*z**5
     .-14680*z**6+28800*z**7-19200*z**8)*g22(z))/(60*z**5)+((1435+547*z
     .+992*z**2-160*z**3+960*z**4)*g23(z))/(480*z**(7.0/2.0))
     .-((1266-3143*z
     .+2647*z**2-585*z**3-130*z**4+120*z**5-120*z**6)*g24(z))/
     .(60*(1-z)*z**5)+((3-42*z+318*z**2-1196*z**3+2196*z**4-1920*z**5
     .+640*z**6)*g31(z))/(4*(1-z)*z)+((1-9*z+9*z**2-z**3-z**4+3*z**5
     .-3*z**6+2*z**7)*g32(z))/(12*(1-z)*z**5)-((1-2*z)*(1-z+z**2)
     .*g34(z))/(2*(1-z)*z)-((3+z**2+2*z**3-z**4+2*z**5)*g35(z))/(4*z**4)
      end
      
      DOUBLE PRECISION FUNCTION Bnf(Z)
      IMPLICIT NONE
      DOUBLE PRECISION nc, cf,ca,nf,tf,z
      DOUBLE PRECISION g11,g12,g21,g22,g23,g24
      DOUBLE PRECISION g25,g26,g31,g32,g33,g34,g35
      nc=3.0d0
      cf=(nc**2-1.0d0)/(2.0d0*nc)
      ca=nc
      nf=5.0d0
      tf=1.0d0/2.0d0
      
      Bnf=-(2050-4115*z+1825*z**2+48*z**3-1568*z**4+8852*z**5-
     .14400*z**6+7200*z**7)/(144*(1-z)*z**4)-((1801-4801*z+3269*z**2
     .-489*z**3-100*z**4+10960*z**5-77700*z**6+193040*z**7-198000*z**8
     .+72000*z**9)*g11(z))/(360*(1-z)*z**5)+((561-939*z+428*z**2+10*z**3
     .+1190*z**4-16650*z**5+60520*z**6-81000*z**7+36000*z**8)*g12(z))/
     .(180*(1-z)*z**4)+((9-24*z+18*z**2-4*z**3-z**7)*g21(z))/
     .(6*(1-z)*z**5)-((187-222*z+72*z**2+920*z**5-7840*z**6+18000*z**7
     .-12000*z**8)*g22(z))/(60*z**5)+((1-3*z)*g23(z))/(48*z**(7.0/2.0))
     .+((7+71*z-66*z**2+8*z**3)*g24(z))/(60*(1-z)*z**5)+
     .2*(1-16*z+66*z**2-100*z**3+50*z**4)*g31(z)
      end
  
      DOUBLE PRECISION FUNCTION g11(z)
      IMPLICIT NONE
      DOUBLE PRECISION z
CC      Cg12(z)->Log(z),
      g11=DLOG(1-z)
      end
      
      DOUBLE PRECISION FUNCTION g12(z)
      DOUBLE PRECISION z
CC      Cg12(z)->Log(z),
      g12=DLOG(z)
      end
      
      DOUBLE PRECISION FUNCTION g21(z)
      DOUBLE PRECISION z,zeta,polylog
CC     Cg21(z)->Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2)),
      g21=Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))
      end
      
      DOUBLE PRECISION FUNCTION g22(z)
      DOUBLE PRECISION z,PolyLog 
Cc      Cg22(z)->PolyLog(2,1-z)-PolyLog(2,z),
      g22=PolyLog(2,1.0d0-z)-PolyLog(2,z)
      end
      
      DOUBLE PRECISION FUNCTION g23(z)
      DOUBLE PRECISION z,polylog
C      Cg23(z)->Log((1-Sqrt(z))/(1+Sqrt(z)))*Log(z)-2*PolyLog(2,-Sqrt(z))+2*PolyLog(2,Sqrt(z)),
      g23=Log((1-Sqrt(z))/(1+Sqrt(z)))*Log(z)-2*PolyLog(2,-Sqrt(z))+
     .2*PolyLog(2,Sqrt(z))
      end
      
      DOUBLE PRECISION FUNCTION g24(z)
      DOUBLE PRECISION z,zeta
C      Cg24(z)->Zeta(2),
      g24=Zeta(2)
      end

      DOUBLE PRECISION FUNCTION g26(z)
      DOUBLE PRECISION z
C      Cg26(z)->Log(1-z)*Log(z),
      g26=Log(1-z)*Log(z)
      end
      
      DOUBLE PRECISION FUNCTION g31(z)
      DOUBLE PRECISION z, zeta, polylog 
C      Cg31(z)->-((-Log(1-z)+Log(z))*(Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))))-6*(PolyLog(3,-(z/(1-z)))-Zeta(3)),
      g31=-((-Log(1-z)+Log(z))*(Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))))
     .-6*(PolyLog(3,-(z/(1-z)))-Zeta(3))
      end
      
      DOUBLE PRECISION FUNCTION g32(z)
C      Cg32(z)->Log(1-z)**3+6*Log(1-z)*PolyLog(2,z)-12*(PolyLog(3,z)+PolyLog(3,-(z/(1-z)))),
            DOUBLE PRECISION z, polylog 
      g32=Log(1-z)**3+6*Log(1-z)*PolyLog(2,z)-12*(PolyLog(3,z)+
     .PolyLog(3,-(z/(1-z))))
      end
      
      DOUBLE PRECISION FUNCTION g33(z)
      DOUBLE PRECISION z,zeta, polylog
C      Cg33(z)->Log(1-z)**3-12*PolyLog(3,z)+6*Log(1-z)*(PolyLog(2,z)-Zeta(2)),
      g33=Log(1-z)**3-12*PolyLog(3,z)+6*Log(1-z)*(PolyLog(2,z)-Zeta(2))
      end
      
      DOUBLE PRECISION FUNCTION g34(z)
      DOUBLE PRECISION z , zeta, polylog
C      Cg34(z)->PolyLog(3,-(z/(1-z)))-3*Log(z)*Zeta(2)+8*Zeta(3),
      g34=PolyLog(3,-(z/(1-z)))-3*Log(z)*Zeta(2)+8*Zeta(3)
      end
      
      DOUBLE PRECISION FUNCTION g35(z)
      DOUBLE PRECISION z, polylog, zeta
C      Cg35(z)->Log((1+Sqrt(z))/(1-Sqrt(z)))**2*Log((1-z)/z)-8*(PolyLog(3,-(Sqrt(z)/(1-Sqrt(z))))+PolyLog(3,Sqrt(z)/(1+Sqrt(z))))+2*PolyLog(3,-(z/(1-z)))+4*Log(1-z)*Zeta(2),
      g35=Log((1+Sqrt(z))/(1-Sqrt(z)))**2*Log((1-z)/z)-8*(PolyLog(3,
     .-(Sqrt(z)/(1-Sqrt(z))))+PolyLog(3,Sqrt(z)/(1+Sqrt(z))))
     .+2*PolyLog(3,-(z/(1-z)))+4*Log(1-z)*Zeta(2)
      end
      
      DOUBLE PRECISION FUNCTION g36(z)
      DOUBLE PRECISION z,  zeta
C      Cg36(z)->Log(1-z)**3-15*Log(1-z)*Zeta(2)
      g36=Log(1-z)**3-15*Log(1-z)*Zeta(2)
      end
      
      DOUBLE PRECISION FUNCTION g37(z)
      DOUBLE PRECISION z, zeta, polylog
C      Cg37(z)->Log(1-z)*(Log(1-z)*Log(z)+PolyLog(2,z)-(15*Zeta(2))/2),
      g37=Log(1-z)*(Log(1-z)*Log(z)+PolyLog(2,z)-(15*Zeta(2))/2)
      end
      
      FUNCTION g38(z)
      DOUBLE PRECISION z, zeta
C      Cg38(z)->Zeta(3),
      g38=Zeta(3)
      end

      DOUBLE PRECISION FUNCTION B(Z)
      IMPLICIT NONE
      DOUBLE PRECISION nc, cf,ca,nf,tf,z
      DOUBLE PRECISION bnf,bnlc,blc
      nc=3.0d0
      cf=(nc**2-1.0d0)/(2.0d0*nc)
      ca=nc
      nf=5.0d0
      tf=1.0d0/2.0d0
      B=cf**2*Blc(z)+cf*(ca-2*cf)*Bnlc(z)+cf*nf*tf*BNf(z)
      end      




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The functions above were copied from Dixon's paper
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     

      DOUBLE PRECISION FUNCTION PolyLog(n,z)
      IMPLICIT NONE
      integer n
      DOUBLE PRECISION z
      double  complex cli2
      double  complex cli3
      double  complex cz
      cz=complex(z,0.0d0)
      if (n.eq.1) then
      PolyLog=-log(1.0d0-z)
      return
      else if  (n.eq.2) then
      PolyLog=dble(cli2(cz))
      return
      else if  (n.eq.3) then
      PolyLog=dble(cli3(cz))
      return
      end if    
      write(*,*)'Problem in polylog'
      stop 
      end

      DOUBLE PRECISION FUNCTION ZETA(z)
      IMPLICIT NONE
      DOUBLE PRECISION pi
      integer z
      parameter(pi=3.141592653589793238d0)
      
      if (z.eq.2) then
      ZETA=pi**2/6.0d0
      return
      else if (z.eq.3) then
      ZETA=1.2020569031595942d+00
      return
      else if (z.eq.4) then
      ZETA=pi**4/90.0d0
      return
      end if
      write(*,*)'Problem in ZETA'
      STOP
      end


       

      
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Functions below are taken from CHAPLIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C=============================================================================
C---  basis expansions
C=============================================================================

C---- expansion of dilogarithm in y = - log(1-z) with Bernoulli numbers  

      double  complex function bsli2_inside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb 
      double precision bern(11)

c     bern(i+1) = BernoulliB(2i)/(2i)!
      data bern /1.D0,0.8333333333333333D-1,-0.1388888888888889D-2
     &,0.3306878306878307D-4,-0.8267195767195767D-6,0.208767569878681D-7
     &,-0.5284190138687493D-9,0.1338253653068468D-10
     &,-0.3389680296322583D-12,0.8586062056277845D-14
     &,-0.2174868698558062D-15/
      parameter (Nmax=11)       ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5d0)
      
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = -zb**2/4d0          !accounting for BernoulliB(1) = -0.5d0
      do i=1,Nmax
         ris = ris + zb**(2*i-1)*bern(i)/(2*i-1)
      enddo
      bsli2_inside = ris
      return 
      end

C---- expansion of the dilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      double  complex function bsli2_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(29),zeta0,zeta2
      
c     zeta(i) = Zeta(2-2i-1)/(2i+1)! i.e. Zeta(-1)/6, Zeta(-3)/120, Zeta(-5)/7!....
      data zeta /-0.01388888888888889d0,0.00006944444444444444d0
     &,-7.873519778281683d-7,1.148221634332745d-8,-1.897886998897100d-10
     &,3.387301370953521d-12,-6.372636443183180d-14,1.246205991295067d-
     &15,-2.510544460899955d-17,5.178258806090624d-19,-1.088735736830085
     &d-20,2.325744114302087d-22,-5.035195213147390d-24,1.10264992943812
     &2d-25,-2.438658550900734d-27,5.440142678856252d-29,-1.222834013121
     &735d-30,2.767263468967951d-32,-6.300090591832014d-34,1.44208683884
     &1848d-35,-3.317093999159543d-37,7.663913557920658d-39,-1.777871473
     &383066d-40,4.139605898234137d-42,-9.671557036081102d-44,2.26671870
     &1676613d-45,-5.327956311328254d-47,1.255724838956433d-48,-2.967000
     &542247094d-50/

      parameter (Nmax=29) ! this is half the order we want (coz even zetaval2 are zero except for 0,2)
      parameter (zeta0 = 1.644934066848226d0)
      parameter (zeta2 = -0.2500000000000000d0)

      zb = log(z)
      ris = dcmplx(zeta0, 0d0) + zb*(1d0 -log(-zb)) 
     &     + zb**2*zeta2
      do i=1,Nmax 
         ris = ris + zb**(2*i+1)*zeta(i)
      enddo
      
      bsli2_outside=ris 
      return 
      end

C---- expansion of trilogarithm in y = - log(1-z) with Bernoulli numbers  
      
      double  complex function bsli3_inside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb 
      double precision bern(21)

c     bern(n+1) = Sum(Bern(n-k)*Bern(k)/(k+1)!(n-k)!,k=0..n)
      data bern /1d0,-0.7500000000000000d0,0.2361111111111111d0,-0.0347
     &2222222222222d0,0.0006481481481481481d0,0.0004861111111111111d0,-
     &0.00002393550012597632d0,-0.00001062925170068027d0,7.794784580498
     &866d-7,2.526087595532040d-7,-2.359163915200471d-8,-6.168132746415
     &575d-9,6.824456748981078d-10,1.524285616929085d-10,-1.91690941417
     &4054d-11,-3.791718683693992d-12,5.277408409541286d-13,9.471165533
     &842511d-13,-1.432311114490360d-14,-2.372464515550457d-15,3.846565
     &792753191d-16/
      
      parameter (Nmax=21)

      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=1,Nmax
         ris = ris + zb**(i)*bern(i)/(i)
      enddo
      bsli3_inside=ris 
      return 
      end

C---- expansion of the trilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      double  complex function bsli3_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(29),zeta0,zeta1,zeta3

c     zeta(i) = Zeta(3-2i-2)/(2i+2)! i.e. Zeta(-1)/24, Zeta(-3)/6!, Zeta(-5)/8!,....
      data zeta /-0.003472222222222222d0,0.00001157407407407407d0,-9.84
     &1899722852104d-8,1.148221634332745d-9,-1.581572499080917d-11,2.41
     &9500979252515d-13,-3.982897776989488d-15,6.923366618305929d-17,-1
     &.255272230449977d-18,2.353754002768465d-20,-4.536398903458687d-22
     &,8.945169670392643d-24,-1.798284004695496d-25,3.675499764793738d-
     &27,-7.620807971564795d-29,1.600041964369486d-30,-3.39676114756037
     &6d-32,7.282272286757765d-34,-1.575022647958003d-35,3.433540092480
     &589d-37,-7.538849998089870d-39,1.666068164765360d-40,-3.703898902
     &881387d-42,8.279211796468275d-44,-1.859914814630981d-45,4.1976272
     &25327060d-47,-9.514207698800454d-49,2.165042825786954d-50,-4.9450
     &00903745158d-52/

      parameter (zeta0 = 1.202056903159594d0)
      parameter (zeta1 = 1.644934066848226d0)
      parameter (zeta3 = -0.08333333333333333d0)
      parameter (Nmax=30)       ! again, half the order, coz odd zetas are zero except 1,3
      zb = log(z)
      ris = dcmplx(zeta0, 0d0) + zb*zeta1 + zb**3*zeta3
     &     + zb**2*(1d0 + 0.5d0 -log(-zb))/2d0
      do i=2,Nmax 
         ris = ris + zb**(2*i)*zeta(i-1)
      enddo
      bsli3_outside=ris 
      return 
      end
      
      
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Functions above are taken from CHAPLIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      