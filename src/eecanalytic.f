      PROGRAM EECANALYTIC
      IMPLICIT NONE
      real*16 GETLO,GETNLO,GETLOMOM,GETNLOMOM,al,pi
      real*16 GETLOINT
      real*16 g11,g12,g21,g22,g23,g24
      real*16 g25,g26,g31,g32,g33,g34,g35
      real*16 A, B
      real*16 bnf,bnlc,blc
      real*16 as1zcut,as2zcut,as3zcut
      integer i,N
      parameter(pi=3.14159265358979323846264338328q0)
      real*16 z
      al=0.1180q0/(2.0q0*pi)      

      open(10,file='doc/predictions/analytic/LOexpanded.dat')
      do i=1,180
      N=50000
      z=(1.0q0-cos(1.0q0*i/180.0*Pi))/2.0q0
      write(10,'(f8.4, E15.6, E15.6)')(1.0q0*i-0.5),
     .al*GETLOINT(1.0q0*i,179.9q0,N),
     .al*as1zcut(z)
      enddo
      close(10)
C      return 


      open(11,file='doc/predictions/analytic/LO.dat')
      do i=1,180
      N=500
      write(11,'(f8.4, E15.6)')(1.0q0*i-0.5),
     .al*GETLO(1.0q0*(1.0q0*i-1.0q0),1.0q0*i,N)      
      enddo
        
      close(11)
      open(12,file='doc/predictions/analytic/NLO.dat')
      do i=1,180
      N=1000
      write(12,'(f8.4, E15.6)')(1.0q0*i-0.5),
     .al*GETLO(1.0q0*(1.0q0*i-1.0q0),1.0q0*i,N)+
     .al*al*GETNLO(1.0q0*(1.0q0*i-1.0q0),1.0q0*i,N)     
      enddo
        
      close(12)
      write(*,*)'Run test at z=0.5'
      write(*,*) 'g11',g11(0.5q0)
      write(*,*) 'g12',g12(0.5q0)
      write(*,*) 'g21',g21(0.5q0)
      write(*,*) 'g22',g22(0.5q0)
      write(*,*) 'g23',g23(0.5q0)
      write(*,*) 'g24',g24(0.5q0)

      write(*,*) 'g31',g31(0.5q0)
      write(*,*) 'g32',g32(0.5q0)
      write(*,*) 'g33',g33(0.5q0)
      write(*,*) 'g34',g34(0.5q0)
      
      write(*,*) 'A', A(0.5q0)
      write(*,*) 'B', B(0.5q0)
      write(*,*) 'Bnf', BNF(0.5q0)
      write(*,*) 'Bnlc', BNlc(0.5q0)
      write(*,*) 'Blc', Blc(0.5q0)
      write(*,*)pi**2/6.0q0
      write(*,*)'EEC LO MOMENTS NUMERIC vs ANALYTIC'
      write(*,*)GETLOMOM(1,2000000),
     .64.0q0*pi/35.0q0
      write(*,*)GETLOMOM(2,200000),
     .4.0q0/9.0q0*(-33.0q0+4.0q0*pi*pi)
      write(*,*)GETLOMOM(3,200000),
     .32.0q0/15.0q0*pi*(-69.0q0+100.0q0*log(2.0q0))
      write(*,*)GETLOMOM(4,200000),
     .-8.0q0/3.0q0*(-277.0q0+28.0q0*pi*pi)
      write(*,*)GETLOMOM(5,200000),
     .-16.0q0/9.0q0*pi*(-857.0q0+1236.0q0*log(2.0q0))
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Simple integration of analytic functions in the given interval
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*16 FUNCTION GETLOINT(th1,th2,N)
      IMPLICIT NONE
      real*16 z,th1, th2, pi,z1,z2,dz,sum,A
      integer i,N
      parameter(pi=3.14159265358979323846264338328q0)
      z1=(1.0q0-cos(th1/180.0q0*pi))/2.0q0
      z2=(1.0q0-cos(th2/180.0q0*pi))/2.0q0
      dz=(Z2-z1)/(1.0q0*N+1.0q0)  
      
        sum=0.0q0
        do i=1,N
        z=z1+dz*i
        sum=sum+A(z)*(1-z)
        enddo
        GETLOINT=SUM/(1.0q0*N)*(z2-z1)*180.0/2.0d0
      end




      real*16 FUNCTION GETLO(th1,th2,N)
      IMPLICIT NONE
      real*16 z,th1, th2, pi,z1,z2,dz,sum,A
      integer i,N
      parameter(pi=3.14159265358979323846264338328q0)
      z1=(1.0q0-cos(th1/180.0q0*pi))/2.0q0
      z2=(1.0q0-cos(th2/180.0q0*pi))/2.0q0
      dz=(Z2-z1)/(1.0q0*N+1.0q0)  
      
        sum=0.0q0
        do i=1,N
        z=z1+dz*i
        sum=sum+A(z)
        enddo
        GETLO=SUM/(1.0q0*N)*(z2-z1)*180.0/2.0d0
      end


      real*16 FUNCTION GETNLO(th1,th2,N)
      IMPLICIT NONE
      real*16  z, th1, th2, pi,z1,z2,dz,sum,B
      integer i,N
      parameter(pi=3.14159265358979323846264338328q0)
      z1=(1.0q0-cos(th1/180.0q0*pi))/2.0q0
      z2=(1.0q0-cos(th2/180.0q0*pi))/2.0q0
      dz=(Z2-z1)/(1.0q0*N+1.0q0)  
        sum=0.0q0
        do i=1,N
        z=z1+dz*i
        sum=sum+B(z)
        enddo
        GETNLO=SUM/(1.0q0*N)*(z2-z1)*180.0/2.0d0
      end


      real*16 FUNCTION GETLOMOM(K,N)
      IMPLICIT NONE
      real*16 z,th1, th2, pi,z1,z2,dz,sum,A
      integer i,N,K
      parameter(pi=3.14159265358979323846264338328q0)
      z1=0.0q0
      z2=1.0q0
      dz=(Z2-z1)/(1.0q0*N+1.0q0)  
      
        sum=0.0q0
        do i=1,N
        z=z1+dz*i
        sum=sum+A(z)*(4*z*(1-z))**(K/2.0)
        enddo
        GETLOMOM=SUM/(1.0q0*N)*(z2-z1)
      end


      real*16 FUNCTION GETNLOMOM(K,N)
      IMPLICIT NONE
      real*16  z, th1, th2, pi,z1,z2,dz,sum,B
      integer i,N,K
      parameter(pi=3.14159265358979323846264338328q0)
      z1=0.0q0
      z2=1.0q0
      dz=(Z2-z1)/(1.0q0*N+1.0q0)  
        sum=0.0q0
        do i=1,N
        z=z1+dz*i
        sum=sum+B(z)*(4*z*(1-z))**(K/2.0)
        enddo
        GETNLOMOM=SUM/(1.0q0*N)*(z2-z1)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The functions bellow were copied from Dixon's paper
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
      real*16 FUNCTION A(Z)
      IMPLICIT NONE
      real*16 nc, cf,ca,nf,tf, z
      nc=3.0q0
      cf=(nc**2.0q0-1.0q0)/(2.0q0*nc)
      ca=nc
      nf=5.0q0
      tf=1.0q0/2.0q0
      
      A=(cf*(-3.0q0+2.0q0*z)*(3.0q0*(2.0q0-3.0q0*z)*z+
     .2.0q0*(3.0q0-6.0q0*z+2.0q0*z**2)*log(1.0q0-z)))/
     .(4.0q0*(-1.0q0+z)*z**5)
      end

      real*16 FUNCTION BLC(Z)
      IMPLICIT NONE
      real*16 nc, cf,ca,nf,tf, z
            real*16 g11,g12,g21,g22,g23,g24
      real*16 g25,g26,g31,g32,g33,g34,g35
      nc=3.0q0
      cf=(nc**2-1.0q0)/(2.0q0*nc)
      ca=nc
      nf=5.0q0
      tf=1.0q0/2.0q0
      
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
      
      
      real*16 FUNCTION BnLc(Z)
      IMPLICIT NONE
      real*16 nc, cf,ca,nf,tf,z
      real*16 g11,g12,g21,g22,g23,g24
      real*16 g25,g26,g31,g32,g33,g34,g35
      nc=3.0q0
      cf=(nc**2-1.0q0)/(2.0q0*nc)
      ca=nc
      nf=5.0q0
      tf=1.0q0/2.0q0
      
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
      
      real*16 FUNCTION Bnf(Z)
      IMPLICIT NONE
      real*16 nc, cf,ca,nf,tf,z
      real*16 g11,g12,g21,g22,g23,g24
      real*16 g25,g26,g31,g32,g33,g34,g35
      nc=3.0q0
      cf=(nc**2-1.0q0)/(2.0q0*nc)
      ca=nc
      nf=5.0q0
      tf=1.0q0/2.0q0
      
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
  
      real*16 FUNCTION g11(z)
      IMPLICIT NONE
      real*16 z
CC      Cg12(z)->Log(z),
      g11=log(1-z)
      end
      
      real*16 FUNCTION g12(z)
      real*16 z
CC      Cg12(z)->Log(z),
      g12=log(z)
      end
      
      real*16 FUNCTION g21(z)
      real*16 z,zeta,polylog
CC     Cg21(z)->Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2)),
      g21=Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))
      end
      
      real*16 FUNCTION g22(z)
      real*16 z,PolyLog 
Cc      Cg22(z)->PolyLog(2,1-z)-PolyLog(2,z),
      g22=PolyLog(2,1.0q0-z)-PolyLog(2,z)
      end
      
      real*16 FUNCTION g23(z)
      real*16 z,polylog
C      Cg23(z)->Log((1-Sqrt(z))/(1+Sqrt(z)))*Log(z)-2*PolyLog(2,-Sqrt(z))+2*PolyLog(2,Sqrt(z)),
      g23=Log((1-Sqrt(z))/(1+Sqrt(z)))*Log(z)-2*PolyLog(2,-Sqrt(z))+
     .2*PolyLog(2,Sqrt(z))
      end
      
      real*16 FUNCTION g24(z)
      real*16 z,zeta
C      Cg24(z)->Zeta(2),
      g24=Zeta(2)
      end

      real*16 FUNCTION g26(z)
      real*16 z
C      Cg26(z)->Log(1-z)*Log(z),
      g26=Log(1-z)*Log(z)
      end
      
      real*16 FUNCTION g31(z)
      real*16 z, zeta, polylog 
C      Cg31(z)->-((-Log(1-z)+Log(z))*(Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))))-6*(PolyLog(3,-(z/(1-z)))-Zeta(3)),
      g31=-((-Log(1-z)+Log(z))*(Log(1-z)**2+2*(PolyLog(2,z)+Zeta(2))))
     .-6*(PolyLog(3,-(z/(1-z)))-Zeta(3))
      end
      
      real*16 FUNCTION g32(z)
C      Cg32(z)->Log(1-z)**3+6*Log(1-z)*PolyLog(2,z)-12*(PolyLog(3,z)+PolyLog(3,-(z/(1-z)))),
            real*16 z, polylog 
      g32=Log(1-z)**3+6*Log(1-z)*PolyLog(2,z)-12*(PolyLog(3,z)+
     .PolyLog(3,-(z/(1-z))))
      end
      
      real*16 FUNCTION g33(z)
      real*16 z,zeta, polylog
C      Cg33(z)->Log(1-z)**3-12*PolyLog(3,z)+6*Log(1-z)*(PolyLog(2,z)-Zeta(2)),
      g33=Log(1-z)**3-12*PolyLog(3,z)+6*Log(1-z)*(PolyLog(2,z)-Zeta(2))
      end
      
      real*16 FUNCTION g34(z)
      real*16 z , zeta, polylog
C      Cg34(z)->PolyLog(3,-(z/(1-z)))-3*Log(z)*Zeta(2)+8*Zeta(3),
      g34=PolyLog(3,-(z/(1-z)))-3*Log(z)*Zeta(2)+8*Zeta(3)
      end
      
      real*16 FUNCTION g35(z)
      real*16 z, polylog, zeta
C      Cg35(z)->Log((1+Sqrt(z))/(1-Sqrt(z)))**2*Log((1-z)/z)-8*(PolyLog(3,-(Sqrt(z)/(1-Sqrt(z))))+PolyLog(3,Sqrt(z)/(1+Sqrt(z))))+2*PolyLog(3,-(z/(1-z)))+4*Log(1-z)*Zeta(2),
      g35=Log((1+Sqrt(z))/(1-Sqrt(z)))**2*Log((1-z)/z)-8*(PolyLog(3,
     .-(Sqrt(z)/(1-Sqrt(z))))+PolyLog(3,Sqrt(z)/(1+Sqrt(z))))
     .+2*PolyLog(3,-(z/(1-z)))+4*Log(1-z)*Zeta(2)
      end
      
      real*16 FUNCTION g36(z)
      real*16 z,  zeta
C      Cg36(z)->Log(1-z)**3-15*Log(1-z)*Zeta(2)
      g36=Log(1-z)**3-15*Log(1-z)*Zeta(2)
      end
      
      real*16 FUNCTION g37(z)
      real*16 z, zeta, polylog
C      Cg37(z)->Log(1-z)*(Log(1-z)*Log(z)+PolyLog(2,z)-(15*Zeta(2))/2),
      g37=Log(1-z)*(Log(1-z)*Log(z)+PolyLog(2,z)-(15*Zeta(2))/2)
      end
      
      FUNCTION g38(z)
      real*16 z, zeta
C      Cg38(z)->Zeta(3),
      g38=Zeta(3)
      end

      real*16 FUNCTION B(Z)
      IMPLICIT NONE
      real*16 nc, cf,ca,nf,tf,z
      real*16 bnf,bnlc,blc
      nc=3.0q0
      cf=(nc**2-1.0q0)/(2.0q0*nc)
      ca=nc
      nf=5.0q0
      tf=1.0q0/2.0q0
      B=cf**2*Blc(z)+cf*(ca-2*cf)*Bnlc(z)+cf*nf*tf*BNf(z)
      end      




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The functions above were copied from Dixon's paper
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     

      real*16 FUNCTION PolyLog(n,z)
      IMPLICIT NONE
      integer n
      real*16 z
      complex*32  cli2
      complex*32  cli3
      complex*32  cz
      cz=complex(z,0.0q0)
      if (n.eq.1) then
      PolyLog=-log(1.0q0-z)
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

      real*16 FUNCTION ZETA(z)
      IMPLICIT NONE
      real*16 pi
      integer z
      parameter(pi=3.14159265358979323846264338328q0)
      
      if (z.eq.2) then
      ZETA=pi**2/6.0q0
      return
      else if (z.eq.3) then
      ZETA=1.20205690315959428539973816151q0
      return
      else if (z.eq.4) then
      ZETA=pi**4/90.0q0
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


      complex*32  function bsli2_inside(z)
      implicit none
      integer i, Nmax
      complex*32 ris, z, zb 
      real*16 bern(11)

c     bern(i+1) = BernoulliB(2i)/(2i)!
      data bern /1.q0,0.8333333333333333q-1,-0.1388888888888889q-2
     &,0.3306878306878307q-4,-0.8267195767195767q-6,0.208767569878681q-7
     &,-0.5284190138687493q-9,0.1338253653068468q-10
     &,-0.3389680296322583q-12,0.8586062056277845q-14
     &,-0.2174868698558062q-15/
      parameter (Nmax=11)       ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5q0)
      
      zb = dcmplx(1q0,0q0)-z
      zb = -log(zb)
      ris = -zb**2/4q0          !accounting for BernoulliB(1) = -0.5q0
      do i=1,Nmax
         ris = ris + zb**(2*i-1)*bern(i)/(2*i-1)
      enddo
      bsli2_inside = ris
      return 
      end

C---- expansion of the dilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      complex*32  function bsli2_outside(z)
      implicit none
      integer i, Nmax
      complex*32 ris, z, zb
      real*16 zeta(29),zeta0,zeta2
     
      
c     zeta(i) = Zeta(2-2i-1)/(2i+1)! i.e. Zeta(-1)/6, Zeta(-3)/120, Zeta(-5)/7!....
      data zeta /-0.01388888888888889q0,0.00006944444444444444q0
     &,-7.873519778281683q-7,1.148221634332745q-8,-1.897886998897100q-10
     &,3.387301370953521q-12,-6.372636443183180q-14,1.246205991295067q-
     &15,-2.510544460899955q-17,5.178258806090624q-19,-1.088735736830085
     &q-20,2.325744114302087q-22,-5.035195213147390q-24,1.10264992943812
     &2q-25,-2.438658550900734q-27,5.440142678856252q-29,-1.222834013121
     &735q-30,2.767263468967951q-32,-6.300090591832014q-34,1.44208683884
     &1848q-35,-3.317093999159543q-37,7.663913557920658q-39,-1.777871473
     &383066q-40,4.139605898234137q-42,-9.671557036081102q-44,2.26671870
     &1676613q-45,-5.327956311328254q-47,1.255724838956433q-48,-2.967000
     &542247094q-50/

      parameter (Nmax=29) ! this is half the order we want (coz even zetaval2 are zero except for 0,2)
      parameter (zeta0 = 1.64493406684822643647241516664654584q0)
      parameter (zeta2 = -0.2500000000000000q0)

      zb = log(z)
      ris = dcmplx(zeta0, 0q0) + zb*(1q0 -log(-zb)) 
     &     + zb**2*zeta2
      do i=1,Nmax 
         ris = ris + zb**(2*i+1)*zeta(i)
      enddo
      
      bsli2_outside=ris 
      return 
      end

C---- expansion of trilogarithm in y = - log(1-z) with Bernoulli numbers  
      
      complex*32  function bsli3_inside(z)
      implicit none
      integer i, Nmax
      complex*32 ris, z, zb 
      real*16 bern(21)

c     bern(n+1) = Sum(Bern(n-k)*Bern(k)/(k+1)!(n-k)!,k=0..n)
      data bern /1q0,-0.7500000000000000q0,0.2361111111111111q0,-0.0347
     &2222222222222q0,0.0006481481481481481q0,0.0004861111111111111q0,-
     &0.00002393550012597632q0,-0.00001062925170068027q0,7.794784580498
     &866q-7,2.526087595532040q-7,-2.359163915200471q-8,-6.168132746415
     &575q-9,6.824456748981078q-10,1.524285616929085q-10,-1.91690941417
     &4054q-11,-3.791718683693992q-12,5.277408409541286q-13,9.471165533
     &842511q-13,-1.432311114490360q-14,-2.372464515550457q-15,3.846565
     &792753191d-16/
      
      parameter (Nmax=21)

      zb = dcmplx(1q0,0q0)-z
      zb = -log(zb)
      ris = dcmplx(0q0, 0q0)
      do i=1,Nmax
         ris = ris + zb**(i)*bern(i)/(i)
      enddo
      bsli3_inside=ris 
      return 
      end

C---- expansion of the trilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      complex*32  function bsli3_outside(z)
      implicit none
      integer i, Nmax
      complex*32 ris, z, zb
      real*16 zeta(29),zeta0,zeta1,zeta3

c     zeta(i) = Zeta(3-2i-2)/(2i+2)! i.e. Zeta(-1)/24, Zeta(-3)/6!, Zeta(-5)/8!,....
      data zeta /-0.003472222222222222q0,0.00001157407407407407q0,-9.84
     &1899722852104q-8,1.148221634332745q-9,-1.581572499080917q-11,2.41
     &9500979252515q-13,-3.982897776989488q-15,6.923366618305929q-17,-1
     &.255272230449977q-18,2.353754002768465q-20,-4.536398903458687q-22
     &,8.945169670392643q-24,-1.798284004695496q-25,3.675499764793738q-
     &27,-7.620807971564795q-29,1.600041964369486q-30,-3.39676114756037
     &6q-32,7.282272286757765q-34,-1.575022647958003q-35,3.433540092480
     &589q-37,-7.538849998089870q-39,1.666068164765360q-40,-3.703898902
     &881387q-42,8.279211796468275q-44,-1.859914814630981q-45,4.1976272
     &25327060q-47,-9.514207698800454q-49,2.165042825786954q-50,-4.9450
     &00903745158q-52/

      parameter (zeta0 = 1.20205690315959428539973816151q0)
      parameter (zeta1 = 1.64493406684822643647241516664654584q0)
      parameter (zeta3 = -0.0833333333333333333333333333q0)
      parameter (Nmax=30)       ! again, half the order, coz odd zetas are zero except 1,3
      zb = log(z)
      ris = dcmplx(zeta0, 0q0) + zb*zeta1 + zb**3*zeta3
     &     + zb**2*(1q0 + 0.5q0 -log(-zb))/2q0
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

      complex*32  function cli2(z)
      implicit none
      complex*32 ris, z, bsli2_inside,bsli2_outside, wcli2
      complex*32 zlocal
      real*16 zabs, pi, zeta2, border, tiny, arg

      pi=3.14159265358979323846264338328q0
      zeta2=pi**2/6q0

      border = 0.3q0 
      tiny = 1d-14
      zabs = abs(z)
      zlocal=z

      if (zabs.gt.1q0+tiny) then
         ris=-wcli2(1q0/z)-zeta2-0.5q0*log(-z)**2
      elseif (zabs.le.border) then 
         ris=bsli2_inside(z)
      else
         if (zabs.gt.1q0) then
            arg=atan2(imag(zlocal),real(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli2_outside(zlocal)
      endif

      cli2=ris
      return
      end
      
c---  recursion
      
      complex*32  function wcli2(z)
      implicit none
      complex*32 z, cli2
      wcli2 =  cli2(z)
      return
      end

c--- Li3

      complex*32  function cli3(z)
      implicit none
      complex*32 ris, z, bsli3_inside,bsli3_outside, wcli3
      complex*32 zlocal
      real*16 zabs,border, pi, zeta2, zeta3,tiny,arg
      
      pi=3.14159265358979323846264338328q0
      zeta2=pi**2/6q0
      zeta3=1.20205690315959428539973816151q0
     
      border = 0.3q0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1q0+tiny) then
         ris=wcli3(1q0/z)-log(-z)**3/6q0-zeta2*log(-z)
      elseif (zabs.le.border) then 
         ris=bsli3_inside(z)
      else
         if (zabs.gt.1q0) then
            arg=atan2(imag(zlocal),real(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli3_outside(zlocal)
      endif

      cli3=ris
      return
      end
      
c---  recursion

      complex*32  function wcli3(z)
      implicit none
      complex*32 z, cli3
      wcli3 =  cli3(z)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Functions above are taken from CHAPLIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Functions below are taken from TG's mail
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      real*16 function as3zcut(z)
      real*16  z
      as3zcut=-15.0496296296296296296296296296q0*log(z)**3
     .+518.865695060272668517465934465q0*log(z)**2
     .-2871.36018216039823956603951601q0*log(z)
            return
      end
     
      real*16 function as2zcut(z)
      real*16  z
      as2zcut=5.76666666666666666666666666667q0*log(z)**2
     .-81.4809061031774183115800436622q0*log(z)
     .+32.7562380661705430311159791007q0
            return
      end
      
       real*16 function as1zcut(z) 
       real*16 nc, cf,ca,nf,tf,z
       nc=3.0q0
       cf=(nc**2-1.0q0)/(2.0q0*nc)
       ca=nc
       nf=5.0q0
       tf=1.0q0/2.0q0
       as1zcut = -2.0q0*log(z) + 1.27777777777q0
     .-0.533333333333333333333q0 * z
      
       if (z.lt. 0.0000001.or.1.0q0-z.lt. 0.0000001)  then
       as1zcut = 0
       return
       end if
       as1zcut = cf * (
     .-3/2*log(z)
     .+1/2*(-1+z)*(3*z**3-13*z**2+23*z-9)*log(1.0q0-z)/z**4
     .-1/4*(-1+z)*(9*z**2-37*z+18)/z**3 )
      return       
       end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C QUADRUPLE PRECISION  BERNULLIB numbers.
C https://www.fulltextarchive.com/page/The-first-498-Bernoulli-Numbers/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      subroutine printBERNULLIB
C     1.0q0/6.0q0,
C     -1.0q0/30.0q0,
C     1.0q0/42.0q0,
C     -1.0q0/30.0q0,
C     5.0q0/66.0q0,
C     -691.0q0/2730.0q0,
C     7.0q0/6.0q0,
C     -3617.0q0/510.0q0,
C     43867.0q0/798.0q0,
C     -174611.0q0/330.0q0,
C     854513.0q0/138.0q0,
C     -236364091.0q0/2730.0q0,
C     8553103.0q0/6.0q0,
C     -23749461029.0q0/870.0q0,
C     8615841276005.0q0/14322.0q0,
C     -7709321041217.0q0/510.0q0,
C     2577687858367.0q0/6.0q0,
C     -26315271553053477373.0q0/1919190.0q0,
C     2929993913841559.0q0/6.0q0,
C     -261082718496449122051.0q0/13530.0q0,
C     1520097643918070802691.0q0/1806.0q0,
C     -27833269579301024235023.0q0/690.0q0,
C     596451111593912163277961.0q0/282.0q0,
C     -5609403368997817686249127547.0q0/46410.0q0,
C     495057205241079648212477525.0q0/66.0q0,
C     -801165718135489957347924991853.0q0/1590.0q0,
C     29149963634884862421418123812691.0q0/798.0q0,
C     -2479392929313226753685415739663229.0q0/870.0q0,
C     84483613348880041862046775994036021.0q0/354.0q0,
C     -1215233140483755572040304994079820246041491.0q0/56786730.0q0,
C     12300585434086858541953039857403386151.0q0/6.0q0,
C     -106783830147866529886385444979142647942017.0q0/510.0q0,
C     1472600022126335654051619428551932342241899101.0q0/64722.0q0,
C     -78773130858718728141909149208474606244347001.0q0/30.0q0,
C     1505381347333367003803076567377857208511438160235.0q0/4686.0q0,
C     -5827954961669944110438277244641067365282488301844260429.0q0/140100870.0q0,
C     34152417289221168014330073731472635186688307783087.0q0/6.0q0,
C     -24655088825935372707687196040585199904365267828865801.0q0/30.0q0,
C     414846365575400828295179035549542073492199375372400483487.0q0/3318.0q0,
C     -4603784299479457646935574969019046849794257872751288919656867.0q0/230010.0q0,
C     1677014149185145836823154509786269900207736027570253414881613.0q0/498.0q0,
C     -2024576195935290360231131160111731009989917391198090877281083932477.0q0/3404310.0q0,
C     660714619417678653573847847426261496277830686653388931761996983.0q0/6.0q0,
C     -1311426488674017507995511424019311843345750275572028644296919890574047.0q0/61410.0q0,
C     1179057279021082799884123351249215083775254949669647116231545215727922535.0q0/272118.0q0,
C     -1295585948207537527989427828538576749659341483719435143023316326829946247.0q0/1410.0q0,
C     1220813806579744469607301679413201203958508415202696621436215105284649447.0q0/6.0q0,
C     -211600449597266513097597728109824233673043954389060234150638733420050668349987259.0q0/4501770.0q0,
C     67908260672905495624051117546403605607342195728504487509073961249992947058239.0q0/6.0q0,
C     -94598037819122125295227433069493721872702841533066936133385696204311395415197247711.0q0/33330.0q0C       
C        

      
      
      
C      end
      