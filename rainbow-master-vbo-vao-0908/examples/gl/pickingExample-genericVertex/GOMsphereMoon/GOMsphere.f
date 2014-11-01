       subroutine GOMsphere(lamda,radius,theta,
     &                     Qabs,QNscat,QFscat,gN,gF,PN,PF)
c  Purpose:
c    To calculate absorption coefficienct, scattering efficiency, asymmetry factor and phase 
c    function for both near and far fields. 
c  Input variables:
c    lamda: wavelength in meter.
c    radius: radius of the spherical scatters ( snow grains)in meters.
c    theta: scattering angle.
c    1. If only need to calculate Qabs, just specify lamda and radius as inputs.
c    2. If only need to calculate Qscat and g, just specify lamda and radius as inputs.
c    3. theta must be specified only for phase (PN, PF) calculation. For other calculation,theta
c       can be given any value, which do not affect the results.
c  Output variables:
c    Qabs: absorption efficiency. It is the same for both near and far fields.
c    QNscat: scattering efficiency for near field scattering.
c    QFscat: scattering efficiency for far field scattering.
c    gN: asymmetry factor for near-filed scattering.
c    gF: asymmetry factor for far-field scattering.
c    PN: phase function for near-field scattering.
c    PF: phase function for far-field scattering.
c  Local variable:
c    n: number of roots for Gauss-Legendre n-point quadratures.
c  NOTE:
c   This program is distributed publically along with a paper published in Applied Optics. Citation is suggested as:
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c
c   Any errors or comments please direct them to Dr. Xiaobing Zhou at xzhou@mtech.edu 
c   or xzhou_2000_2001@yahoo.com. Thanks.  May, 2008.
c
       real*8 lamda, radius
       real*8 Qabs,QNscat,QFscat,gN,gF,theta,PN,PF
       integer*4 n
       
       n=20

c  To calcualte the absorption efficiency Qabs.
       call abs_efficiency(lamda,radius,n,Qabs)

c  To calcualte the scattering efficiency and asymmetry factor for both near and far fields
       call QscatAndG(lamda,radius,n,QNscat,QFscat,gN,gF)

c  To calculate the phase function for at a scattering angle.
       call phasefunction(lamda,radius,theta,QNscat,QFscat,
     &                          PN,PF)

       return 
       end

c     **************************************************************
c     *********To calculate the absorption efficiency Qabs**********
c     **************************************************************
      subroutine abs_efficiency(lamda,radius,n,Qabs)
c PURPOSE:
c     To calculate absorption efficiency, given wavelength, radius of sphere
c REFERENCE: Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
c             Orlando, FL: Academic Press,1984. p.481-483.
c INPUT:
c     lamda: wavelength
c     radius: radius of sphere
c     n: number of roots for Gaussian quadratures.
c OUTPUT:
c     Qabs: absorption efficiency
c LOCAL:
c     temp: temperature, used to obtain the complex refractive index; for VNIR,
c           m is independent on temp, but for microwave, it is.
c     m: refractive index
c     beta: absorption coefficient
c     eps: tolerance for truncation. = truncation error.
c     NN: number of interfaces of reflection for ray inside a sphere.
      real*8 lamda,radius,temp,Qabs
      integer*4 n,i,NN
      real*8 eps
      complex*16 m
      real*8 re,im,pi,xi
      real*8 x(n),w(n)
      integer*4 kind
      real*8 endpts(2), bb(n)
      real*8 a, b, G, absFF
      integer*4 mm,k
      real*8 xm,xk,dx
      eps=1.0d-15
      pi = dacos(-1.0d0)
      a=0.0d0
      b=1.0d0
      xm=0.5d0*(b+a)
      xk=0.5d0*(b-a)
      mm=(n+1)/2
      k=2*mm-n+1
      Qabs=0.0d0
      call GaussLeg(n,x,w)
c  To make in_angle lie in 0 to pi/2, x(i) must be selected as positive.
      do i=mm+1,n
         dx=xk*x(i)
         Qabs=Qabs+w(i)*(absFF(xm+dx,lamda,radius)
     &                  +absFF(xm-dx,lamda,radius))
      end do
      Qabs=Qabs*xk
      return
      end

       real*8 function absFF(x,lamda,radius)
c  PURPOSE:
c      To calculate the integrand function f(x) in the expression of absorption efficiency
c      
c  REFERENCE: 
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c  INPUT VARIABLES:
c      x: varible, cosine of the incident angle
c      lamda: wavelength in meter
c      radius: radius of scatter in meter
c  OUTPUT VARIABLES:
c      absFF: value of the integrand function for absorption efficiency calculation. eq.(21b)
c  LOCAL VARIABLES:
c      j: loop varible
c      beta: absorption coefficient
c      xi: path length between two reflection events inside the sphere
c      temp: snow temperature in K, only as an input to refice,no effect in VNIR
c      m: complex refractive index
c      in_angle:incident angle in radian
c      TT: transmission for natural light
c      RR: reflectance for natural light
c      Trp: transmission for parallel-polarized light
c      Trv: transmission for vertical-polarized light
c      Rep: reflection for parallel-polarized light
c      Rev: reflection for vertical-polarized light
c      tp: complex, amplitude coefficient of transmission for parallel-polarization
c      tv: complex, amplitude coefficient of transmission for vertical-polarization
c      rp: complex, amplitude coefficient of reflection for parallel-polarization
c      rv: complex, amplitude coefficient of reflection for vertical-polarization
c      t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c      r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c      u: mcos(theta)=u+iv
c      v:

      real*8 x,lamda,radius
      complex*16 m
      real*8 beta,xi,temp,in_angle,TT, RR,Trp,Trv
      complex*16 Rep,Rev,tp,tv,rp,rv,t,r
      real *8 u,v
      
      temp=293.0d0
c     temp in Kelvin where 0C = 273 K, the above temperature makes the refractive index
c     indepent of temperature, because the index data has temperature dependency on when
c     temperature is below 0C.
 
c  To obtain refractive index m
      call refice( lamda, temp, m )
      
c  To obtain absorption coefficient beta
      call abs_coef(lamda, m, beta)
      
c  To obtain a path length xi
      in_angle=dacos(x)
      call path_length(radius,in_angle,m,xi)

      call TransReflec(in_angle,m,TT,RR,Trp,Trv,
     &                        Rep,Rev,tp,tv,rp,rv,t,r,u,v)
      absFF=2.0*TT/(1.-RR*dexp(-beta*xi))*(1.0d0-dexp(-beta*xi))*x

      return
      end
       
c     ************************************************************************************
c     *********To calculate the scattering efficiency Qscat and asymmetry factor**********
c     ************************************************************************************
      subroutine QscatAndG(lamda,radius,n,QNscat,QFscat,gN,gF)
c PURPOSE:
c     To calculate scattering efficiency and asymmetry factor for both near filed (QNscat) 
c     and far filed (QFscat), given wavelength, radius of sphere, number of roots for 
c     Gaussian quadratures. Eq.(26c)+eq.(41). 
c REFERENCE: 
c     Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
c     Orlando, FL: Academic Press,1984. p.481-483.
c INPUT:
c     lamda: wavelength
c     radius: radius of sphere
c     n: number of roots for Gaussian quadratures.
c OUTPUT:
c     QFscat: scattering efficiency for far field
c     QNscat: scattering efficiency for near field
c     gN: asymmetry factor for near field
c     gF: asymmetry factor for far field
c LOCAL:
c     temp: temperature, used to obtain the complex refractive index; for VNIR,
c           m is independent on temp, but for microwave, it is.
c     m: refractive index
c     beta: absorption coefficient
c     eps: tolerance for truncation. = truncation error. 10^-10
c     NN: number of interfaces of reflection for ray inside a sphere.
c     theta: scattering angle
c     in_angle: real, incident angle
c     Trans: transmission
c     Reflec: reflectance
c     Trp: transmission for parallel-polarized light
c     Trv: transmission for vertical-polarized light
c     Rep: reflection for parallel-polarized light
c     Rev: reflection for vertical-polarized light
c     tp: real, amplitude coefficient of transmission for parallel-polarization
c     tv: real, amplitude coefficient of transmission for vertical-polarization
c     rp: real, amplitude coefficient of reflection for parallel-polarization
c     rv: real, amplitude coefficient of reflection for vertical-polarization
c     t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c     r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c     i: loop variable
c     j: jth interface
c     gNQscat: g*Qscat for near field
c     gFQscat: g*Qsact for far field
c     scatF1,scatF2: value of integrand function in Qscat expression
c     asymNF1,asymNF2: value of integrand function in g*Qscat for near field (see Eq.(48b))
c     asymFF1,asymFF2: value of integrand function in g*Qscat for far field (see Eq.(48b))
c     u: mcos(theta)=u+iv
c     v:

      real*8 lamda,radius
      real*8 QFscat,QNscat,gNQscat,gFQscat,gN,gF
      real*8 Fscat,Nscat
      real*8 scatfN1,scatfN2,scatfFB1,scatfFB2
      real*8 asymNF1,asymNF2,asymFFB1,asymFFB2
      integer*4 n
      real*8 temp, theta, in_angle,Trans, Reflec,Trp,Trv
      complex*16 Rep,Rev,tp,tv,rp,rv,t,r,thetat
      real*8 beta,eps,re,im,pi,xi,u,v
      complex*16 m
      real*8 x(n),w(n),scatF,xm,xk,dx, xmp,xkp,dxp
      integer*4 mm,k,i,NN
      real*8 a,b,bp,xacc,x1,x2,epsilon
 
       parameter (maxdiv=10)
       real*8 D
       integer*4 nbb, N2
       real*8 xb1(maxdiv),xb2(maxdiv), rts(maxdiv)
       
       nsegments=maxdiv
       nb=maxdiv
       pi = dacos(-1.0d0)
       a=0.0d0
       b=pi
       bp=pi/2.0
       xacc = 1.0d-15
       x1=0.0d0
       x2=1.0d0
       epsilon = 1.0e-8

       temp=270.0d0
       eps=1.0d-15
c  To obtain refractive index m: lamda(wavelength) and temp(temperature) are inputs. 
       call refice( lamda, temp, m )
c       
c  To obtain the abscissas x(n) and weights w(n) of the Gauss-Legendre n-point quadrature
c     formula: n (number of roots for Gaussian quadratures) is input.
       call GaussLeg(n,x,w)
c       
       xm=0.5d0*(b+a)
       xk=0.5d0*(b-a)
       xmp=0.5d0*(bp+a)
       xkp=0.5d0*(bp-a)
       mm=(n+1)/2
       k=2*mm-n+1
       QNscat=0.0d0
       QFscat=0.0d0
       Nscat=0.0d0
       Fscat=0.0d0
c  Product of g and Q for both near and far fields
       gNQscat=0.0d0
       gFQscat=0.0d0
c  Calculating scattering efficiency:
          call scat_efficiency(lamda,radius,m,QNscat,QFscat)
c          write(*,*)'QNscat=',QNscat,'QFscat=',QFscat

c  To make theta lie in 0 to pi, x(i) must be selected as positive.
       do i=mm+1,n
          dx=xk*x(i)
          dxp=xkp*x(i)
c  Calculating scattering asymmetry factor g:
          call QtimesG(xm+dx,m,lamda,radius,
     &                            scatfN1,scatfFB1,asymNF1,asymFFB1)
          call QtimesG(xm-dx,m,lamda,radius,
     &                            scatfN2,scatfFB2,asymNF2,asymFFB2)
          Nscat=Nscat+w(i)*(scatfN1 +scatfN2)
          gNQscat=gNQscat+w(i)*(asymNF1 +asymNF2)
          
          call QtimesG(xmp+dxp,m,lamda,radius,
     &                            scatfN1,scatfFB1,asymNF1,asymFFB1)
          call QtimesG(xmp-dxp,m,lamda,radius,
     &                            scatfN2,scatfFB2,asymNF2,asymFFB2)
c          write(*,*)'scatfFB1=',scatfFB1
c          write(*,*)'scatfFB2=',scatfFB2
          Fscat=Fscat+w(i)*(scatfFB1 +scatfFB2)
          gFQscat=gFQscat+w(i)*(asymFFB1 +asymFFB2)
       end do

       Nscat=Nscat*xk
       gNQscat=gNQscat*xk
       gN=gNQscat/Nscat
       
       Fscat=Fscat*xkp
       gFQscat=gFQscat*xkp
       Fscat=Fscat+Nscat
       gFQscat=gFQscat+gNQscat
       gF=gFQscat/Fscat

       return
       end


      subroutine scat_efficiency(lamda,radius,m,QNscat,QFscat)
c PURPOSE:
c     To calculate scattering efficiency for both near filed (QNscat) and far filed
c     (QFscat), given wavelength, radius of sphere, number of roots for Gaussian 
c     quadratures.
c REFERENCE: 
c     Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
c     Orlando, FL: Academic Press,1984. p.481-483.
c INPUT:
c     lamda: wavelength
c     radius: radius of sphere
c     m: refractive index
c OUTPUT:
c     QFscat: scattering efficiency for far field
c     QNscat: scattering efficiency for near field
c LOCAL:
c     a: real, lower limit of integration
c     b: real, higher limit of integration
c     n: number of roots for Gaussian quadratures.
c     temp: temperature, used to obtain the complex refractive index; for VNIR,
c           m is independent on temp, but for microwave, it is.
c     beta: scatorption coefficient
c     eps: tolerance for truncation. = truncation error. 10^-10
c     NN: number of interfaces of reflection for ray inside a sphere.
c     in_angle: real, incident angle
c     Trans: transmission
c     Reflec: reflectance
c     Trp: transmission for parallel-polarized light
c     Trv: transmission for vertical-polarized light
c     Rep: reflection for parallel-polarized light
c     Rev: reflection for vertical-polarized light
c     tp: real, amplitude coefficient of transmission for parallel-polarization
c     tv: real, amplitude coefficient of transmission for vertical-polarization
c     rp: real, amplitude coefficient of reflection for parallel-polarization
c     rv: real, amplitude coefficient of reflection for vertical-polarization
c     t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c     r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c     i: loop variable

      real*8 a,b,lamda,radius
      complex*16 m
      real*8 QFscat,QNscat
      integer*4 n
      real*8 pi
      real*8 x(20),w(20),scatFNF,xm,xk,dx
      integer*4 mm,k,i
      n=20
      pi = dacos(-1.0d0)
      a=0.0d0
      b=1.0d0

c  To obtain the abscissas and weights of the Gauss-Legendre n-point quadrature
c     formula.
       call GaussLeg(n,x,w)
c       
       xm=0.5d0*(b+a)
       xk=0.5d0*(b-a)
       mm=(n+1)/2
       k=2*mm-n+1
       QNscat=0.0d0
       QFscat=0.0d0
c  To make in_angle lie in 0 to pi/2, x(i) must be selected as positive.
       do i=mm+1,n
          dx=xk*x(i)
          
c  Calculating scattering efficiency:
          QNscat=QNscat+w(i)*(scatFNF(xm+dx,lamda,radius)
     &                   +scatFNF(xm-dx,lamda,radius))
       end do
       
       QNscat=QNscat*xk
       QFscat=QNscat+1.0d0
       
       return
       end

       real*8 function scatFNF(x,lamda,radius)
c  PURPOSE:
c      To calculate the integrand function f(x) in the expression of scattering efficiency
C      for near field.
c      
c  REFERENCE: 
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c  INPUT VARIABLES:
c      x: varible, cosine of the incident angle, also the abscissas of the Gauss-Legendre
c         n-point quadrature.
c      lamda: wavelength of incident wave in unit of meter
c      radius: radius of snow grain as a sphere
c  OUTPUT VARIABLES:
c      scatFNF: value of the integrand function for near filed scattering efficiency 
c               calculation. See eq.(23)
c  LOCAL VARIABLES:
c      NN: the number of interfaces when truncation is carried out with truncation 
c         error = 10^-10
c      m: complex, refractive index
c      beta: absorption coefficients
c      xi: real, path-length bwtween two reflecting interfaces
c      temp: snow temperature in K, only as an input to refice,no effect in VNIR
c      epspjp: real, amplitude of scattered parallel-polarized electric field
c      epspjv: real, amplitude of scattered vertically-polarized electric field
c      Trans: transmission
c      Reflec: reflectance
c      Trp: transmission for parallel-polarized light
c      Trv: transmission for vertical-polarized light
c      Rep: reflection for parallel-polarized light
c      Rev: reflection for vertical-polarized light
c      tp: real, amplitude coefficient of transmission for parallel-polarization
c      tv: real, amplitude coefficient of transmission for vertical-polarization
c      rp: real, amplitude coefficient of reflection for parallel-polarization
c      rv: real, amplitude coefficient of reflection for vertical-polarization
c      t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c      r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c      j: loop variable, jth interface
c      
      real*8 x,lamda,radius
      integer*4 NN
      complex*16 m
      real*8 temp,beta,thetai,xi,eps
      complex*16 Rep,Rev,tp,tv,rp,rv,t,r
      real*8 Trans, Reflec, Trp,Trv,u,v
      integer*4 j
      scatFNF=0.0d0
      eps=1.0d-15
      temp=270.0d0

c  To obtain refractive index m
      call refice( lamda, temp, m )
      
c  To obtain absorption coefficient beta
      call abs_coef(lamda, m, beta)
      
c  To obtain a path length xi
      thetai = dacos(x)
      call path_length(radius,thetai,m,xi)

       call TransReflec(thetai,m,Trans,Reflec,Trp,Trv,
     &                                     Rep,Rev,tp,tv,rp,rv,t,r,u,v)
c  To obtain truncation number NN
      call sumNumber3(beta,xi,t,Reflec,Trans,eps,NN)

       do j=2,NN
          scatFNF=scatFNF+(Trp*Trp*(cdabs(Rep))**(j-2)
     &                   +Trv*Trv*(cdabs(Rev))**(j-2))
     &                   *dexp(-beta*xi*dfloat(j-1))
       end do
       scatFNF=scatFNF+ 2.*Reflec
       scatFNF = scatFNF * x
       
       return
       end

       subroutine QtimesG(theta,m,lamda,radius,
     &                            scatfN,scatfFB,asymNF,asymFFB)
c  PURPOSE:
c      To calculate the integrand function f(x) in the expression of scattering
c      efficiency (eq.(26b)) and the integrand function f(x)in the expression of 
c      asymmetry factor for near field scattering and far field. 
c      For Far-field, Qsca^far=Qsca^near+1
c      
c  REFERENCE: 
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c  INPUT VARIABLES:
c      theta: varible, scattering angle in radians
c      m: complex, refractive index
c      beta: absorption coefficients
c      xi: real, path-length bwtween two interfaces
c      NN: the number of interfaces when truncation is carried out with truncation 
c         error = 10^-10
c  OUTPUT VARIABLES:
c      scatfN: value of the integrand function for near-field scattering efficiency calculation
c      asymNF: value of the integrand function for asymmetry factor calculation.
c              asymNF=gQscat for near field.
c      asymFF: value of the integrand function for asymmetry factor calculation.
c              asymFF=gQscat for far field.

c  LOCAL VARIABLES:
c      epspjp: real, amplitude of scattered parallel-polarized electric field
c      epspjv: real, amplitude of scattered vertically-polarized electric field
c      Trans: transmission
c      Reflec: reflectance
c      Trp: transmission for parallel-polarized light
c      Trv: transmission for vertical-polarized light
c      Rep: reflection for parallel-polarized light
c      Rev: reflection for vertical-polarized light
c      tp: real, amplitude coefficient of transmission for parallel-polarization
c      tv: real, amplitude coefficient of transmission for vertical-polarization
c      rp: real, amplitude coefficient of reflection for parallel-polarization
c      rv: real, amplitude coefficient of reflection for vertical-polarization
c      t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c      r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c      j: loop variable, jth interface
c      thetai: incident angle
c      x: size parameter = 2*pi*radius/wavelength
c      
       real*8 BessJ1
       
       real*8 scatfN,scatfFB,asymNF,asymFFB
       real*8 theta,beta,xi,lamda,radius
       complex*16 m,thetat
       integer*4 NN
       integer*4 nb,nsegments
       real*8 thetai,Trans, Reflec,Trp,Trv
       complex*16 Rep,Rev,tp,tv,rp,rv,t,r
       complex*16 epspjp,epspjv
       parameter (maxdiv=10)
       parameter (Nmax=1000)
       real*8 re,D
       integer*4 nbb,i,j,mm, N2
       real*8 xb1(maxdiv),xb2(maxdiv), rts(maxdiv)
       real*8 pi,u,v,a,b,xacc,x1,x2,epsilon,x
       nsegments=maxdiv
       nb=maxdiv

       a=0.0d0
       b=1.0d0
       xacc = 1.0d-15
       x1=0.0d0
       x2=1.0d0
       pi=dacos(-1.0d0)
       epsilon = 1.0e-16
       x=2.0*pi*radius/lamda
       re=dreal(m)
       j=0
       mm=0

       scatfN=0.0d0

       do while (j.le.Nmax)
          j=j+1
c  To find nb and rts(incident angles)
          call roots(j,m,theta,x1,x2,xacc, nsegments, nb, rts)
          if (nb.eq.0) go to 10
          do 9 i = 1, nb
              thetai=2.0*datan(rts(i))
              
c  subroutine TransReflec(..) is in Trans_Ref.f. To calculate tp,tv,rp,rv,t,r
              call TransReflec(thetai,m,Trans,Reflec,Trp,Trv,
     &                                     Rep,Rev,tp,tv,rp,rv,t,r,u,v)
     
c  subroutine (..) ScatFieldAmp(...) is in Trans_Ref.f. To calculate epspjp,epspjv
c  To obtain path-length xi
              call path_length(radius,thetai,m,xi)
              
c  To obtain absorption coefficient beta
              call abs_coef(lamda, m, beta)
              
              call ScatFieldAmp(thetai,tp,tv,rp,rv,beta,xi,j,u,v,
     &                         epspjp,epspjv)
              D=dsin(2.*thetai)/(4.0*dsin(theta))
     &          /dabs(dfloat(j-1)*dcos(thetai)/dsqrt(u*u+v*v)-1.0)

               scatfN=scatfN
     &                  +(cdabs(epspjp)**2+cdabs(epspjv)**2)*D
     
               mm=mm+1
9         continue              
10       continue
         end do
cccc         scatfN,scatfFB,asymNF,asymFFB
c  See eq.(26b)-(26c)
         scatfN=scatfN*dsin(theta)

c  refraction part of integrand function of far-filed scattering efficiency 
         scatfFB=2.0*BessJ1(x*dsin(theta))*BessJ1(x*dsin(theta))
     &          /dsin(theta)

c integrand function for bear-field scattering
         asymNF=scatfN*dcos(theta)
c refraction part of integrand function of far-filed asymmetry facotr (eq.50)
         asymFFB=scatfFB*dcos(theta)
 
       return
       end

       subroutine ScatFieldAmp(thetai,tp,tv,rp,rv,beta,xi,j,u,v,
     &                         epspjp,epspjv)
c PURPOSE:
c     To evaluate amplitude of scattered electric field in term of incident amplitude.
c  REFERENCE: 
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c  INPUT VARIABLES:
c      thetai: real, the incident angle (in radians)
c      m: complex, refractive index
c      tp: real, amplitude coefficient of transmission for parallel-polarization
c      tv: real, amplitude coefficient of transmission for vertical-polarization
c      rp: real, amplitude coefficient of reflection for parallel-polarization
c      rv: real, amplitude coefficient of reflection for vertical-polarization
c      beta: absorption coefficient
c      xi: path length between two reflection events inside the sphere
c      j: integer, jth interface of a ray inteacting with the sphere interface
c  OUTPUT VARIABLES:
c      epspjp: real, amplitude of scattered parallel-polarized electric field
c      epspjv: real, amplitude of scattered vertically-polarized electric field
c  LOCAL VARIABLES:
c       thetat: real, refraction angle (in radians)

       real*8 thetai,beta,xi,u,v
       complex*16 tp,tv,rp,rv,thetat
       complex*16 epspjp,epspjv
       integer*4 j

       if (j.eq.1) then
         epspjp = rp
         epspjv = rv
       else
         epspjp = complex(u,v)/dcos(thetai)*tp*tp*(-rp)**(j-2)
     &         *dexp(-beta*xi*(j-1)/2.0d0)
         epspjv = complex(u,v)/dcos(thetai)*tv*tv*(-rv)**(j-2)
     &         *dexp(-beta*xi*(j-1)/2.0d0)
       end if
       return
       end

       real*8 function BessJ1(x)
c Purpose:
c    Returns the value of the first order Bessel function J1(x) for any real x
c Reference: 
c    Press W H , S A Teukolsky, W T Vetterling and B P Flannery,"Numerical 
c    recipes in C: the art of scientific computing", Second edition, Cambridge
c    University Press,1992. p.233.
       real*8 ax,z
       real*8 x,y,ans1,ans2,xx
       ax=dabs(x)
       if(ax.lt.8.0) then
          y=x*x
          ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
     &         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))))
          ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
     &         +y*(99447.43394+y*(376.9991397+y*1.0))))
          BessJ1=ans1/ans2
        else 
          z=8.0/ax
          y=z*z
          xx=ax-2.356194491
          ans1=1.0+y*(0.183105d-2+y*(-0.3516396496d-4
     &         +y*(0.2457520174d-5+y*(-0.240337019d-6))))
          ans2=0.04687499995+y*(-0.2002690873d-3
     &        +y*(0.8449199096d-5+y*(-0.88228987d-6+y*0.105787412d-6)))
          BessJ1=dsqrt(0.636619772/ax)*(dcos(xx)*ans1-z*dsin(xx)*ans2)
          if(x.lt.0.0d0) BessJ1=-BessJ1
        endif

        return
        end
        
c     *************************************************************************
c     *********To calculate the phase function for a scattering angle**********
c     *************************************************************************
       subroutine phasefunction(lamda,radius,theta,QNscat,QFscat,
     &                          Nphase,Fphase)
c  Purpose:
c     To obtain phase function for given wavelength and radius of a sphere
c  Input:
c     lamda: wavelength, in meters
c     radius: radius of the scatter, in meters
c     theta: scattering angle
c     QFscat: scattering efficiency for far field used in scat_efficiency(...)in Q_scat2.f
c     QNscat: scattering efficiency for near field used in scat_efficiency(...)in Q_scat2.f
c  Output:
c     Nphase: phase function for near field scattering
c     Fphase: phase function for far field scattering
c  Local:
c     a: real, lower limit of integration, used in scat_efficiency(...)in Q_scat.f
c     b: real, higher limit of integration, used in scat_efficiency(...)in Q_scat.f
c     TEMP: = temperature (K) ( for WAVLEN.GT.167 microns only )
c                            (range:  213.16 to 272.16) used in refice(...) in ref_ice.f
c     theta: outgoing (scattering) angle
c     xacc: accuracy that is used to find the root.Used in roots(...) in Pphase.f

c ......in subroutine refice( WAVMET, TEMP, ref_ice ) in ref_ice.f
c .    epsilon: real, tolerance for cutoff in the summation series.
c .    N: integer, maximum number of summation of for absorption efficiency
c .    r: real, amplitude coefficient of reflection
c .    t: real, amplitude coefficient of transmission for unpolarized wave
c ......in subroutine refice( WAVMET, TEMP, ref_ice ) in ref_ice.f
c
c ......in subroutine roots(j,m,theta,x1,x2,xacc, nsegments, nb, rts)
c .     j: jth interface. 
c .     m: complex refractive index.
c .     theta: outgoing angle (scattering angle).
c .     x1,x2: limits that all the roots are to be found. 
c .     xacc: accuracy that is used to find the root.
c .     nb: as input, it is the maximum number of roots sought
c .     nsegments: number that the inteval [x1,x2] is equally subdivided in order to bracket = maximum
c .        possible roots
c .     nb: as output, it is the number of different roots that are found.
c .     rts: the different roots found that lies within [x1,x2]
c ......in subroutine roots(j,m,theta,x1,x2,xacc, nsegments, nb, rts)
c
c ......in subroutine TransReflec(thetai,m,Trans,Reflec,Trp,Trv,
c .     &                                     Rep,Rev,tp,tv,rp,rv,t,r) in Trans_Ref.f
c .   thetai: incident angle
c .   m: complex refractive index
c .   Trans: transmission for natural light
c .   Reflec: reflectance for natural light
c .   Trp: transmission for parallel-polarized light
c .   Trv: transmission for vertical-polarized light
c .   Rep: reflection for parallel-polarized light
c .   Rev: reflection for vertical-polarized light
c .    tp: real, amplitude coefficient of transmission for parallel-polarization
c .    tv: real, amplitude coefficient of transmission for vertical-polarization
c .    rp: real, amplitude coefficient of reflection for parallel-polarization
c .    rv: real, amplitude coefficient of reflection for vertical-polarization
c .    t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c .    r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c ......in subroutine TransReflec(thetai,m,Trans,Reflec,Trp,Trv,
c .     &                                     Rep,Rev,tp,tv,rp,rv,t,r) in Trans_Ref.f
c
c ......in subroutine ScatFieldAmp(thetai,m,tp,tv,rp,rv,beta,xi,j,
c .    &                         epspjp,epspjv) in Q_scat.f
c .     thetai: real, the incident angle (in radians)
c .     m: complex, refractive index
c .     tp: real, amplitude coefficient of transmission for parallel-polarization
c .     tv: real, amplitude coefficient of transmission for vertical-polarization
c .     rp: real, amplitude coefficient of reflection for parallel-polarization
c .     rv: real, amplitude coefficient of reflection for vertical-polarization
c .     beta: absorption coefficient
c .     xi: path length between two reflection events inside the sphere
c .     j: integer, jth interface of a ray inteacting with the sphere interface
c .     epspjp: real, amplitude of scattered parallel-polarized electric field
c .     epspjv: real, amplitude of scattered vertically-polarized electric field
c ......in subroutine ScatFieldAmp(thetai,m,tp,tv,rp,rv,beta,xi,j,
c .    &                         epspjp,epspjv) in Q_scat.f

c     re: real part of complex refractive index
c     D: divergence factor
c      x: size parameter. = 2*pi*radius/wavelength
c      PQN: P*Q for near field scattering
c      PQF: P*Q for far field scattering
       real*8 lamda,radius
       real*8 theta,Nphase,Fphase,PQN,PQF
       real *8 a,b
       integer*4 N,j
       real*8 QFscat,QNscat,TEMP,xacc,epsilon
       real*8 x1,x2
       integer*4 nb,nsegments,maxdiv,Nmax
       real*8 thetai,Trans,Reflec,Trp,Trv
       complex*16 Rep,Rev,tp,tv,rp,rv,t,r,thetat
       real*8 beta,xi,pi,x,u,v,BessJ1
       complex*16 m, epspjp,epspjv
       
       parameter (maxdiv=10)
c       parameter (Nmax=1000)
       parameter (Nmax=100)
       real*8 re,D
       integer*4 nbb,i,mm, N2
       real*8 xb1(maxdiv),xb2(maxdiv), rts(maxdiv)
       nsegments=maxdiv
       nb=maxdiv

       xacc = 1.0d-15
       x1=0.0d0
       x2=1.0d0
       pi=dacos(-1.0d0)
       epsilon = 1.0e-8

c To obtain complex refractive index (m). subroutine refice( ... ) is in ref_ice.f
       TEMP=270.0
       call refice( lamda, TEMP, m )

       re=dreal(m)
       Nphase=0.0d0
       Fphase=0.0d0
       PQN=0.0d0
       PQF=0.0d0
       mm=0
       j=0
       do while (j.le.Nmax)
          j=j+1

          call roots(j,m,theta,x1,x2,xacc, nsegments, nb, rts)
          
          if (nb.eq.0) go to 10
          do 9 i = 1, nb
              thetai=2.0*datan(rts(i))
c To obtain absorption coefficient              
              call abs_coef(lamda, m, beta)
c To obtain path length between two consecutive incident interfaces
              call path_length(radius,thetai,m,xi)
c To obtain interface number 
              call sumNumber3(beta,xi,t,Reflec,Trans,epsilon,N2)
              
c  subroutine TransReflec(..) is in Trans_Ref.f. To calculate tp,tv,rp,rv,t,r
              call TransReflec(thetai,m,Trans,Reflec,Trp,Trv,
     &                                     Rep,Rev,tp,tv,rp,rv,t,r,u,v)
     
c  subroutine (..) ScatFieldAmp(...) is in Trans_Ref.f. To calculate epspjp,epspjv
              call ScatFieldAmp(thetai,tp,tv,rp,rv,beta,xi,j,u,v,
     &                         epspjp,epspjv)
              D=dsin(2.*thetai)/(4.0*dsin(theta))
     &          /dabs(dfloat(j-1)*dcos(thetai)/dsqrt(u*u+v*v)-1.0)
               PQN=PQN+2.*(cdabs(epspjp)**2+cdabs(epspjv)**2)*D
              mm=mm+1
 9         continue              
10       continue
         end do
         x=2.*pi*radius/lamda
         if (theta.lt.pi/2.) then
            PQF=PQN+4.*BessJ1(x*dsin(theta))*BessJ1(x*dsin(theta))
     &         /dsin(theta)/dsin(theta)
         else
            PQF=PQN
         endif
            Nphase=PQN/QNscat
            Fphase=PQF/QFscat

       return
       end


c     ******************************************************************************
c     *********To calculate the complex refractive index for a wavelength **********
c     ******************************************************************************
      subroutine refice( WAVMET, TEMP, ref_ice )
c        Calculates complex refractive index of Ice 1H for wavelengths
c        between 45 nm and 8.6 m.  For wavelengths above 167 microns,
c        temperature dependence is included for temperatures between
c        213 and 272K.  Mainly intended for applications in Earth ice
c        clouds and snow, not other planets or interstellar space;
c        the temperature dependence or crystalline form of ice may be
c        incorrect for these latter applications.
c      I N P U T :  
c                   WAVMET = wavelength (meters)
c                   WAVLEN = wavelength (microns)
c                            (range:  0.0443 to 8.600E+06)
c                   TEMP   = temperature (K) ( for WAVLEN.GT.167 only )
c                            (range:  213.16 to 272.16)
c      O U T P U T :  ref_ice = complex refractive index
c                              ( with positive imaginary part )
c      (WARNING:  input out of range will print a warning message and
c                 return ref_ice=(0.,0.) in order not to unnecessarily 
c                 halt the calling program;  the calling program should
c                 test the real part of ref_ice to catch these errors)
c      METHOD :  Tabular interpolation, assuming
c                (1) real index is linear in log(wavelength)
c                    and linear in temperature
c                (2) log(imag. index) is linear in log(wavelength)
c                    and linear in temperature
c     AUTHORS OF subroutine refice( WAVMET, TEMP, ref_ice ):
c             Stephen Warren, Univ. of Washington (1983)
c               (sgw@cloudy.atmos.washington.edu)
c             Bo-Cai Gao, JCESS, Univ. of Maryland (1995)
c               (gao@imagecube.gsfc.nasa.gov)
c             Warren Wiscombe, NASA Goddard (1995)
c               (wiscombe@climate.gsfc.nasa.gov)
c     MODIFICATIONS IN 1995 :
c       Gao, Warren, and (to a small extent) Wiscombe modified the
c       original Warren refice program from 1984 to change values of
c       imaginary refractive index in the 0.161-0.410 and 1.445-2.50
c       micron regions.  The values in 0.161-0.410 were incorrect and
c       the values in 1.445-2.50 were among the most uncertain in 1984.
c       New measurements have made it possible to improve both regions.

c       No changes were made to real refractive indices (by re-doing a
c       Kramers-Kronig analysis), because the values of imaginary
c       index MIM involved are so small (below 0.001) that the
c       resulting changes in real index MRE would be in the third
c       decimal place at most.  (MIM has negligible effect on MRE
c       when MIM << MRE.)

c       The 0.161-0.410 micron region was changed using data provided
c       by Warren, which correct his misinterpretation of Minton's
c       measurements for 0.181-0.185 micron, and incorporate new
c       measurements of Perovich and Govoni (1991) for 0.250-0.400
c       micron.  Warren (1984) correctly represented UV measurements
c       of Seki et al. and visible measurements of Grenfell/Perovich,
c       but he plotted Minton's measurements a factor of 2.3 too low
c       because he misinterpreted base-10 as base-e.  (The UV
c       measurements of Dressler/Schnepp and Shibaguchi et al are also
c       probably expressed as absorption coefficients on base-10;
c       therefore those values also were probably plotted a factor of
c       2.3 too low in Warren's (1984) Figure 2.)

c       The details of how the present imaginary index data for
c       0.161-0.410 micron is obtained are as follows.  Point A in
c       Warren's Figure 2 at 161 nm is joined with a straight line to
c       Minton's corrected point B at 181 nm.  Minton's reported
c       values for 181-185 nm have been smoothed within his stated
c       uncertainty.  Now a smooth curve is drawn to join Minton at
c       185 nm to Perovich/Govoni (PG) at 250 nm.  PG's values from
c       their Table 1 show some unrealistic wiggles that are smaller
c       than their error bars, so a smooth curve was fitted through
c       them and values were taken from the smoothed curve at 10-nm
c       intervals.  PG ends at 400 nm, where Grenfell/Perovich (GP)
c       starts.  At 400 nm we take imaginary index=2.82E-9, the
c       average of PG (2.93E-9) and GP (2.71E-9).

c       The Warren (1984) values of imaginary index in the 1.445-2.50
c       micron region were replaced by those of Kou et al.(1994).  In
c       order to remove the resulting discontinuities near 1.445 and
c       2.5 micron, the Warren values at 1.43 and 1.44 micron were
c       changed to 0.9E-04 and 1.3E-04 respectively, and his values at
c       2.52, 2.55, and 2.565 micron were changed to 8.255E-04,
c       8.578E-04, and 8.739E-04, respectively. The latter change
c       eliminated a small local maximum at 2.5 micron which was not
c       realistic and has never been seen in spectra of snow bracketing
c       that wavelength.


c     REFERENCES :

c       Warren, S., 1984: Optical Constants of Ice from the Ultraviolet
c          to the Microwave, Appl. Opt. 23, 1206-1225

c       Kou, L., D. Labrie, and P. Chylek, 1994: Refractive indices
c          of water and ice in the 0.65- to 2.5-micron spectral range,
c          Appl. Opt. 32, 3531-3540

c       Perovich, D., and J. Govoni, 1991: Absorption Coefficients
c          of Ice from 250 to 400 nm, Geophys. Res. Lett. 18, 1233-1235
c ======================================================================

      IMPLICIT NONE

c     .. Parameters ..

      INTEGER   NWL, NWLT
      PARAMETER ( NWL = 574, NWLT = 62)
c     ..
c     .. Scalar Arguments ..

      REAL*8 TEMP, WAVLEN,WAVMET
      complex*16 ref_ice
c     ..
c     .. Local Scalars ..

      CHARACTER MESSAG*40
      LOGICAL   PASS1
      INTEGER   I, L
      REAL*8      FRAC, MIM, MRE, YHI, YLO
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC LOG, CMPLX, EXP
c     ..
c                                        ** Refractive index table

      REAL*8      TABIM( NWL ), TABIMT( NWLT, 4 ), TABRE( NWL ),
     &          TABRET( NWLT, 4), TEMREF(4), WL( NWL ), WLT( NWLT )

      COMMON / ICEREF / WL, WLT, TABRE, TABRET, TABIM, TABIMT, TEMREF

      DATA  PASS1 / .True. /

       WAVLEN = WAVMET*1.0d+6
       IF( PASS1 )  THEN

         PASS1 = .False.
c                         ** Superficially test if main table messed up

         IF( NWL.LT.100 ) CALL ERRMSG('refice--NWL value bad',.True.)
         IF( WL(1).GT.0.045 ) CALL ERRMSG('refice--WL(1) bad',.True.)
         IF( WL(NWL).LT.166.) CALL ERRMSG('refice--WL(NWL) bad',.True.)

         DO 1  I = 1, NWL

            IF( WL(I).LE.0.0 .OR. TABRE(I).LE.0.5 .OR. TABRE(I).GT.2.0
     &          .OR. TABIM(I).LT.0.0  .OR. TABIM(I).GT.10.0 )  THEN
               WRITE( MESSAG, '(A,I5,A)' )  'refice--table value ', I,
     &                ' out of bounds '
               CALL ERRMSG( MESSAG, .True. )
            END IF

            IF( I.GT.1 .AND. WL(I).LE.WL(I-1) )  THEN
               WRITE( MESSAG, '(A,I5,A)' )  'refice--table WL(', I,
     &                ') not increasing  '
               CALL ERRMSG( MESSAG, .True. )
            END IF

    1    CONTINUE

      END IF


      IF( WAVLEN.LT.WL(1) .OR. WAVLEN.GT.WLT(NWLT) ) THEN

         CALL ERRMSG('refice--wavelength outside table boundaries',
     &                .False.)
         ref_ice = (0.,0.)
         RETURN

      END IF


      IF( WAVLEN.LE.167.) THEN
c                                  ** Wavelength between 0.045 and 167
c                                  ** microns. No temperature dependence
         DO 10 I = 2, NWL
            IF( WAVLEN.LE.WL(I)) GO TO 20
   10    CONTINUE

   20    CONTINUE
         FRAC   = LOG( WAVLEN / WL(I-1) ) / 
     &            LOG( WL(I)  / WL(I-1) )
         MRE    = TABRE(I-1) + FRAC * ( TABRE(I) - TABRE(I-1) )
         MIM    = TABIM(I-1) * ( TABIM(I) / TABIM(I-1) )**FRAC


      ELSE

c               ** Wavelength greater than 167 microns
c               ** (temperature-dependent case)

c        write (*,*)'temp=',TEMP

        IF( TEMP.LT.TEMREF(4) .OR. TEMP.GT.TEMREF(1) ) THEN

            CALL ERRMSG('refice--temperature outside table boundaries',
     &                  .False.)
            ref_ice = (0.,0.)
            RETURN

         END IF
c                         ** Find position in temperature array
         DO 30 L = 2, 4
            IF( TEMP.GE.TEMREF(L) ) GO TO 40
   30    CONTINUE
c                         ** Find position in wavelength array
   40    CONTINUE
         DO 50 I = 2, NWLT
            IF( WAVLEN.LE.WLT(I) ) GO TO 60
   50    CONTINUE

   60    CONTINUE
         FRAC   = LOG( WAVLEN / WLT(I-1) ) /
     &            LOG( WLT(I) / WLT(I-1) )

         YLO    = TABRET(I-1, L) +
     &            FRAC*( TABRET(I, L) - TABRET(I-1, L) )

         YHI    = TABRET( I-1, L-1) +
     &            FRAC*( TABRET(I, L-1) - TABRET(I-1, L-1) )

         MRE    = YLO + ( YHI - YLO) * ( TEMP - TEMREF(L) ) /
     &                                 ( TEMREF(L-1) - TEMREF(L) )

         YLO    = LOG( TABIMT(I-1, L)) +
     &            FRAC*LOG( TABIMT(I, L) / TABIMT(I-1, L) )

         YHI    = LOG( TABIMT(I-1, L-1) ) +
     &            FRAC*LOG( TABIMT(I, L-1) / TABIMT(I-1, L-1) )

         MIM    = EXP( YLO + (YHI - YLO) * (TEMP - TEMREF(L)) /
     &                                     (TEMREF(L-1) - TEMREF(L)) )

      END IF


      ref_ice = CMPLX( MRE, MIM)

      END

      BLOCK DATA ICECON

c        Ice-refractive-index vs. wavelength table for subroutine refice

      IMPLICIT NONE

c     .. Parameters ..

      INTEGER   NWL, NWLT
      PARAMETER ( NWL = 574, NWLT = 62 )
c     ..
c     .. Local Scalars ..

      INTEGER   I
c     ..

      REAL*8      TABIM(NWL), TABIMT( NWLT, 4 ), TABRE(NWL),
     &          TABRET( NWLT, 4 ), TEMREF(4), WL(NWL), WLT(NWLT)

      COMMON / ICEREF / WL, WLT, TABRE, TABRET, TABIM, TABIMT, TEMREF


c       WL, WLT           wavelengths (microns) for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TABRE, TABRET     real refractive indices for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TABIM, TABIMT     imaginary refractive indices for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TEMREF            reference temperatures (-1,-5,-20,-60 deg C)
c                         for TABRET,TABIMT


      DATA  TEMREF / 272.16, 268.16, 253.16, 213.16 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ), I = 1, 10 ) /
     & 4.43000E-02, 8.34410E-01, 1.64000E-01,
     & 4.51000E-02, 8.36760E-01, 1.73000E-01,
     & 4.59000E-02, 8.37290E-01, 1.83000E-01,
     & 4.68000E-02, 8.37710E-01, 1.95000E-01,
     & 4.77000E-02, 8.38270E-01, 2.08000E-01,
     & 4.86000E-02, 8.40380E-01, 2.23000E-01,
     & 4.96000E-02, 8.47190E-01, 2.40000E-01,
     & 5.06000E-02, 8.55220E-01, 2.50000E-01,
     & 5.17000E-02, 8.60470E-01, 2.59000E-01,
     & 5.28000E-02, 8.62480E-01, 2.68000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 11, 20 ) /
     & 5.39000E-02, 8.61570E-01, 2.79000E-01,
     & 5.51000E-02, 8.60930E-01, 2.97000E-01,
     & 5.64000E-02, 8.64190E-01, 3.19000E-01,
     & 5.77000E-02, 8.69160E-01, 3.40000E-01,
     & 5.90000E-02, 8.77640E-01, 3.66000E-01,
     & 6.05000E-02, 8.92960E-01, 3.92000E-01,
     & 6.20000E-02, 9.10410E-01, 4.16000E-01,
     & 6.36000E-02, 9.30890E-01, 4.40000E-01,
     & 6.53000E-02, 9.53730E-01, 4.64000E-01,
     & 6.70000E-02, 9.81880E-01, 4.92000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 21, 30 ) /
     & 6.89000E-02, 1.02334, 5.17000E-01,
     & 7.08000E-02, 1.06735, 5.28000E-01,
     & 7.29000E-02, 1.11197, 5.33000E-01,
     & 7.38000E-02, 1.13134, 5.34000E-01,
     & 7.51000E-02, 1.15747, 5.31000E-01,
     & 7.75000E-02, 1.20045, 5.24000E-01,
     & 8.00000E-02, 1.23840, 5.10000E-01,
     & 8.27000E-02, 1.27325, 5.00000E-01,
     & 8.55000E-02, 1.32157, 4.99000E-01,
     & 8.86000E-02, 1.38958, 4.68000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 31, 40 ) /
     & 9.18000E-02, 1.41644, 3.80000E-01,
     & 9.30000E-02, 1.40906, 3.60000E-01,
     & 9.54000E-02, 1.40063, 3.39000E-01,
     & 9.92000E-02, 1.40169, 3.18000E-01,
     & 1.03300E-01, 1.40934, 2.91000E-01,
     & 1.07800E-01, 1.40221, 2.51000E-01,
     & 1.10000E-01, 1.39240, 2.44000E-01,
     & 1.12700E-01, 1.38424, 2.39000E-01,
     & 1.14000E-01, 1.38075, 2.39000E-01,
     & 1.18100E-01, 1.38186, 2.44000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 41, 50 ) /
     & 1.21000E-01, 1.39634, 2.47000E-01,
     & 1.24000E-01, 1.40918, 2.24000E-01,
     & 1.27200E-01, 1.40256, 1.95000E-01,
     & 1.29500E-01, 1.38013, 1.74000E-01,
     & 1.30500E-01, 1.36303, 1.72000E-01,
     & 1.31900E-01, 1.34144, 1.80000E-01,
     & 1.33300E-01, 1.32377, 1.94000E-01,
     & 1.34800E-01, 1.30605, 2.13000E-01,
     & 1.36200E-01, 1.29054, 2.43000E-01,
     & 1.37000E-01, 1.28890, 2.71000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 51, 60 ) /
     & 1.37800E-01, 1.28931, 2.89000E-01,
     & 1.38700E-01, 1.30190, 3.34000E-01,
     & 1.39300E-01, 1.32025, 3.44000E-01,
     & 1.40900E-01, 1.36302, 3.82000E-01,
     & 1.42500E-01, 1.41872, 4.01000E-01,
     & 1.43500E-01, 1.45834, 4.06500E-01,
     & 1.44200E-01, 1.49028, 4.05000E-01,
     & 1.45000E-01, 1.52128, 3.89000E-01,
     & 1.45900E-01, 1.55376, 3.77000E-01,
     & 1.46800E-01, 1.57782, 3.45000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 61, 70 ) /
     & 1.47600E-01, 1.59636, 3.32000E-01,
     & 1.48000E-01, 1.60652, 3.15000E-01,
     & 1.48500E-01, 1.61172, 2.98000E-01,
     & 1.49400E-01, 1.61919, 2.74000E-01,
     & 1.51200E-01, 1.62522, 2.28000E-01,
     & 1.53100E-01, 1.63404, 1.98000E-01,
     & 1.54000E-01, 1.63689, 1.72000E-01,
     & 1.55000E-01, 1.63833, 1.56000E-01,
     & 1.56900E-01, 1.63720, 1.10000E-01,
     & 1.58000E-01, 1.63233, 8.30000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 71, 80 ) /
     & 1.58900E-01, 1.62222, 5.80000E-02,
     & 1.61000E-01, 1.58270, 2.20000E-02,
     & 1.62700E-01, 1.55360, 1.00000E-02,
     & 1.65200E-01, 1.52040, 3.00000E-03,
     & 1.67500E-01, 1.49840, 1.00000E-03,
     & 1.70000E-01, 1.48010, 3.00000E-04,
     & 1.72300E-01, 1.46710, 1.00000E-04,
     & 1.74800E-01, 1.45510, 3.00000E-05,
     & 1.77100E-01, 1.44580, 1.00000E-05,
     & 1.79600E-01, 1.43700, 3.00000E-06 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 81, 90 ) /
     & 1.81000E-01, 1.43250, 1.56700E-06,
     & 1.82000E-01, 1.42950, 9.32500E-07,
     & 1.83000E-01, 1.42680, 5.39700E-07,
     & 1.84000E-01, 1.42410, 3.12200E-07,
     & 1.85000E-01, 1.42140, 1.72500E-07,
     & 1.88000E-01, 1.41430, 1.00000E-07,
     & 1.90000E-01, 1.41010, 8.20000E-08,
     & 1.95000E-01, 1.40070, 5.10000E-08,
     & 2.00000E-01, 1.39360, 3.81000E-08,
     & 2.05000E-01, 1.38670, 3.05000E-08 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 91, 100 ) /
     & 2.10000E-01, 1.38000, 2.51000E-08,
     & 2.15000E-01, 1.37610, 2.18000E-08,
     & 2.20000E-01, 1.37230, 1.98000E-08,
     & 2.25000E-01, 1.36850, 1.78000E-08,
     & 2.30000E-01, 1.36480, 1.62000E-08,
     & 2.35000E-01, 1.36120, 1.50000E-08,
     & 2.40000E-01, 1.35770, 1.43000E-08,
     & 2.45000E-01, 1.35420, 1.37000E-08,
     & 2.50000E-01, 1.35080, 1.33000E-08,
     & 2.60000E-01, 1.34720, 1.32000E-08 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 101, 110 ) /
     & 2.70000E-01, 1.34370, 1.30000E-08,
     & 2.80000E-01, 1.34030, 1.26000E-08,
     & 2.90000E-01, 1.33710, 1.18000E-08,
     & 3.00000E-01, 1.33390, 1.10000E-08,
     & 3.10000E-01, 1.33200, 9.28000E-09,
     & 3.20000E-01, 1.33020, 8.25000E-09,
     & 3.30000E-01, 1.32840, 7.65000E-09,
     & 3.40000E-01, 1.32660, 7.00000E-09,
     & 3.50000E-01, 1.32490, 6.15000E-09,
     & 3.60000E-01, 1.32380, 5.10000E-09 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 111, 120 ) /
     & 3.70000E-01, 1.32260, 4.13000E-09,
     & 3.80000E-01, 1.32150, 3.43000E-09,
     & 3.90000E-01, 1.32040, 3.12000E-09,
     & 4.00000E-01, 1.31940, 2.82000E-09,
     & 4.10000E-01, 1.31850, 2.51000E-09,
     & 4.20000E-01, 1.31775, 2.26000E-09,
     & 4.30000E-01, 1.31702, 2.08000E-09,
     & 4.40000E-01, 1.31633, 1.91000E-09,
     & 4.50000E-01, 1.31569, 1.54000E-09,
     & 4.60000E-01, 1.31509, 1.53000E-09 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 121, 130 ) /
     & 4.70000E-01, 1.31452, 1.55000E-09,
     & 4.80000E-01, 1.31399, 1.64000E-09,
     & 4.90000E-01, 1.31349, 1.78000E-09,
     & 5.00000E-01, 1.31302, 1.91000E-09,
     & 5.10000E-01, 1.31257, 2.14000E-09,
     & 5.20000E-01, 1.31215, 2.26000E-09,
     & 5.30000E-01, 1.31175, 2.54000E-09,
     & 5.40000E-01, 1.31136, 2.93000E-09,
     & 5.50000E-01, 1.31099, 3.11000E-09,
     & 5.60000E-01, 1.31064, 3.29000E-09 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 131, 140 ) /
     & 5.70000E-01, 1.31031, 3.52000E-09,
     & 5.80000E-01, 1.30999, 4.04000E-09,
     & 5.90000E-01, 1.30968, 4.88000E-09,
     & 6.00000E-01, 1.30938, 5.73000E-09,
     & 6.10000E-01, 1.30909, 6.89000E-09,
     & 6.20000E-01, 1.30882, 8.58000E-09,
     & 6.30000E-01, 1.30855, 1.04000E-08,
     & 6.40000E-01, 1.30829, 1.22000E-08,
     & 6.50000E-01, 1.30804, 1.43000E-08,
     & 6.60000E-01, 1.30780, 1.66000E-08 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 141, 150 ) /
     & 6.70000E-01, 1.30756, 1.89000E-08,
     & 6.80000E-01, 1.30733, 2.09000E-08,
     & 6.90000E-01, 1.30710, 2.40000E-08,
     & 7.00000E-01, 1.30688, 2.90000E-08,
     & 7.10000E-01, 1.30667, 3.44000E-08,
     & 7.20000E-01, 1.30646, 4.03000E-08,
     & 7.30000E-01, 1.30625, 4.30000E-08,
     & 7.40000E-01, 1.30605, 4.92000E-08,
     & 7.50000E-01, 1.30585, 5.87000E-08,
     & 7.60000E-01, 1.30566, 7.08000E-08 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 151, 160 ) /
     & 7.70000E-01, 1.30547, 8.58000E-08,
     & 7.80000E-01, 1.30528, 1.02000E-07,
     & 7.90000E-01, 1.30509, 1.18000E-07,
     & 8.00000E-01, 1.30491, 1.34000E-07,
     & 8.10000E-01, 1.30473, 1.40000E-07,
     & 8.20000E-01, 1.30455, 1.43000E-07,
     & 8.30000E-01, 1.30437, 1.45000E-07,
     & 8.40000E-01, 1.30419, 1.51000E-07,
     & 8.50000E-01, 1.30402, 1.83000E-07,
     & 8.60000E-01, 1.30385, 2.15000E-07 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 161, 170 ) /
     & 8.70000E-01, 1.30367, 2.65000E-07,
     & 8.80000E-01, 1.30350, 3.35000E-07,
     & 8.90000E-01, 1.30333, 3.92000E-07,
     & 9.00000E-01, 1.30316, 4.20000E-07,
     & 9.10000E-01, 1.30299, 4.44000E-07,
     & 9.20000E-01, 1.30283, 4.74000E-07,
     & 9.30000E-01, 1.30266, 5.11000E-07,
     & 9.40000E-01, 1.30249, 5.53000E-07,
     & 9.50000E-01, 1.30232, 6.02000E-07,
     & 9.60000E-01, 1.30216, 7.55000E-07 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 171, 180 ) /
     & 9.70000E-01, 1.30199, 9.26000E-07,
     & 9.80000E-01, 1.30182, 1.12000E-06,
     & 9.90000E-01, 1.30166, 1.33000E-06,
     & 1.00000, 1.30149, 1.62000E-06,
     & 1.01000, 1.30132, 2.00000E-06,
     & 1.02000, 1.30116, 2.25000E-06,
     & 1.03000, 1.30099, 2.33000E-06,
     & 1.04000, 1.30082, 2.33000E-06,
     & 1.05000, 1.30065, 2.17000E-06,
     & 1.06000, 1.30048, 1.96000E-06 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 181, 190 ) /
     & 1.07000, 1.30031, 1.81000E-06,
     & 1.08000, 1.30014, 1.74000E-06,
     & 1.09000, 1.29997, 1.73000E-06,
     & 1.10000, 1.29979, 1.70000E-06,
     & 1.11000, 1.29962, 1.76000E-06,
     & 1.12000, 1.29945, 1.82000E-06,
     & 1.13000, 1.29927, 2.04000E-06,
     & 1.14000, 1.29909, 2.25000E-06,
     & 1.15000, 1.29891, 2.29000E-06,
     & 1.16000, 1.29873, 3.04000E-06 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 191, 200 ) /
     & 1.17000, 1.29855, 3.84000E-06,
     & 1.18000, 1.29837, 4.77000E-06,
     & 1.19000, 1.29818, 5.76000E-06,
     & 1.20000, 1.29800, 6.71000E-06,
     & 1.21000, 1.29781, 8.66000E-06,
     & 1.22000, 1.29762, 1.02000E-05,
     & 1.23000, 1.29743, 1.13000E-05,
     & 1.24000, 1.29724, 1.22000E-05,
     & 1.25000, 1.29705, 1.29000E-05,
     & 1.26000, 1.29686, 1.32000E-05 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 201, 210 ) /
     & 1.27000, 1.29666, 1.35000E-05,
     & 1.28000, 1.29646, 1.33000E-05,
     & 1.29000, 1.29626, 1.32000E-05,
     & 1.30000, 1.29605, 1.32000E-05,
     & 1.31000, 1.29584, 1.31000E-05,
     & 1.32000, 1.29563, 1.32000E-05,
     & 1.33000, 1.29542, 1.32000E-05,
     & 1.34000, 1.29521, 1.34000E-05,
     & 1.35000, 1.29499, 1.39000E-05,
     & 1.36000, 1.29476, 1.42000E-05 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 211, 220 ) /
     & 1.37000, 1.29453, 1.48000E-05,
     & 1.38000, 1.29430, 1.58000E-05,
     & 1.39000, 1.29406, 1.74000E-05,
     & 1.40000, 1.29381, 1.98000E-05,
     & 1.41000, 1.29355, 2.50000E-05,
     & 1.42000, 1.29327, 5.40000E-05,
     & 1.43000, 1.29299, 9.00000E-05,
     & 1.44000, 1.29272, 1.30000E-04,
     & 1.44510, 1.29260, 1.60000E-04,
     & 1.44930, 1.29250, 1.89000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 221, 230 ) /
     & 1.45350, 1.29240, 2.20000E-04,
     & 1.45770, 1.29230, 2.58000E-04,
     & 1.46200, 1.29230, 2.99000E-04,
     & 1.46630, 1.29220, 3.42000E-04,
     & 1.47060, 1.29210, 3.86000E-04,
     & 1.47490, 1.29200, 4.33000E-04,
     & 1.47930, 1.29190, 4.71000E-04,
     & 1.48370, 1.29180, 5.03000E-04,
     & 1.48810, 1.29180, 5.27000E-04,
     & 1.49250, 1.29170, 5.34000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 231, 240 ) /
     & 1.49700, 1.29160, 5.38000E-04,
     & 1.50150, 1.29160, 5.33000E-04,
     & 1.50600, 1.29150, 5.27000E-04,
     & 1.51060, 1.29140, 5.17000E-04,
     & 1.51520, 1.29130, 5.09000E-04,
     & 1.51980, 1.29120, 4.98000E-04,
     & 1.52440, 1.29110, 4.87000E-04,
     & 1.52910, 1.29100, 4.75000E-04,
     & 1.53370, 1.29090, 4.62000E-04,
     & 1.53850, 1.29080, 4.51000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 241, 250 ) /
     & 1.54320, 1.29070, 4.37000E-04,
     & 1.54800, 1.29060, 4.25000E-04,
     & 1.55280, 1.29050, 4.08000E-04,
     & 1.55760, 1.29040, 3.92000E-04,
     & 1.56250, 1.29030, 3.75000E-04,
     & 1.56740, 1.29020, 3.58000E-04,
     & 1.57230, 1.29000, 3.42000E-04,
     & 1.57730, 1.28990, 3.24000E-04,
     & 1.58230, 1.28970, 3.11000E-04,
     & 1.58730, 1.28960, 2.97000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 251, 260 ) /
     & 1.59240, 1.28950, 2.84000E-04,
     & 1.59740, 1.28930, 2.74000E-04,
     & 1.60260, 1.28920, 2.66000E-04,
     & 1.60770, 1.28900, 2.59000E-04,
     & 1.61290, 1.28890, 2.53000E-04,
     & 1.61810, 1.28880, 2.49000E-04,
     & 1.62340, 1.28860, 2.47000E-04,
     & 1.62870, 1.28840, 2.44000E-04,
     & 1.63400, 1.28830, 2.40000E-04,
     & 1.63930, 1.28810, 2.35000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 261, 270 ) /
     & 1.64470, 1.28800, 2.29000E-04,
     & 1.65020, 1.28780, 2.23000E-04,
     & 1.65560, 1.28760, 2.16000E-04,
     & 1.66110, 1.28750, 2.11000E-04,
     & 1.66670, 1.28730, 2.05000E-04,
     & 1.67220, 1.28710, 1.97000E-04,
     & 1.67790, 1.28700, 1.92000E-04,
     & 1.68350, 1.28680, 1.87000E-04,
     & 1.68920, 1.28660, 1.82000E-04,
     & 1.69490, 1.28640, 1.78000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 271, 280 ) /
     & 1.70070, 1.28620, 1.72000E-04,
     & 1.70650, 1.28600, 1.67000E-04,
     & 1.71230, 1.28580, 1.63000E-04,
     & 1.71820, 1.28560, 1.59000E-04,
     & 1.72410, 1.28540, 1.54000E-04,
     & 1.73010, 1.28520, 1.52000E-04,
     & 1.73610, 1.28500, 1.49000E-04,
     & 1.74220, 1.28470, 1.48000E-04,
     & 1.74830, 1.28450, 1.45000E-04,
     & 1.75440, 1.28430, 1.42000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 281, 290 ) /
     & 1.76060, 1.28410, 1.40000E-04,
     & 1.76680, 1.28380, 1.39000E-04,
     & 1.77300, 1.28360, 1.38000E-04,
     & 1.77940, 1.28340, 1.37000E-04,
     & 1.78570, 1.28310, 1.35000E-04,
     & 1.79210, 1.28290, 1.33000E-04,
     & 1.79860, 1.28260, 1.31000E-04,
     & 1.80510, 1.28240, 1.30000E-04,
     & 1.81160, 1.28210, 1.28000E-04,
     & 1.81820, 1.28190, 1.26000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 291, 300 ) /
     & 1.82480, 1.28160, 1.25000E-04,
     & 1.83150, 1.28130, 1.23000E-04,
     & 1.83820, 1.28100, 1.23000E-04,
     & 1.84500, 1.28070, 1.26000E-04,
     & 1.85190, 1.28040, 1.26000E-04,
     & 1.85870, 1.28010, 1.33000E-04,
     & 1.86570, 1.27970, 1.42000E-04,
     & 1.87270, 1.27940, 1.56000E-04,
     & 1.87970, 1.27900, 1.87000E-04,
     & 1.88680, 1.27870, 2.33000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 301, 310 ) /
     & 1.89390, 1.27830, 2.98000E-04,
     & 1.90110, 1.27790, 4.00000E-04,
     & 1.90840, 1.27750, 4.84000E-04,
     & 1.91570, 1.27720, 5.53000E-04,
     & 1.92310, 1.27690, 6.82000E-04,
     & 1.93050, 1.27660, 8.38000E-04,
     & 1.93800, 1.27630, 9.91000E-04,
     & 1.94550, 1.27600, 1.16000E-03,
     & 1.95310, 1.27570, 1.31000E-03,
     & 1.96080, 1.27540, 1.42000E-03 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 311, 320 ) /
     & 1.96850, 1.27510, 1.52000E-03,
     & 1.97630, 1.27480, 1.59000E-03,
     & 1.98410, 1.27460, 1.63000E-03,
     & 1.99200, 1.27430, 1.64000E-03,
     & 2.00000, 1.27400, 1.64000E-03,
     & 2.00800, 1.27370, 1.63000E-03,
     & 2.01610, 1.27340, 1.59000E-03,
     & 2.02430, 1.27310, 1.55000E-03,
     & 2.03250, 1.27280, 1.52000E-03,
     & 2.04080, 1.27240, 1.45000E-03 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 321, 330 ) /
     & 2.04920, 1.27200, 1.40000E-03,
     & 2.05760, 1.27160, 1.32000E-03,
     & 2.06610, 1.27120, 1.23000E-03,
     & 2.07470, 1.27070, 1.13000E-03,
     & 2.08330, 1.27030, 1.02000E-03,
     & 2.09210, 1.26980, 9.30000E-04,
     & 2.10080, 1.26930, 8.20000E-04,
     & 2.10970, 1.26870, 7.19000E-04,
     & 2.11860, 1.26810, 6.28000E-04,
     & 2.12770, 1.26750, 5.50000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 331, 340 ) /
     & 2.13680, 1.26690, 4.85000E-04,
     & 2.14590, 1.26620, 4.28000E-04,
     & 2.15520, 1.26550, 3.83000E-04,
     & 2.16450, 1.26480, 3.48000E-04,
     & 2.17390, 1.26410, 3.17000E-04,
     & 2.18340, 1.26330, 2.93000E-04,
     & 2.19300, 1.26260, 2.71000E-04,
     & 2.20260, 1.26180, 2.54000E-04,
     & 2.21240, 1.26100, 2.38000E-04,
     & 2.22220, 1.26020, 2.27000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 341, 350 ) /
     & 2.23210, 1.25930, 2.17000E-04,
     & 2.24220, 1.25840, 2.12000E-04,
     & 2.25230, 1.25750, 2.10000E-04,
     & 2.26240, 1.25650, 2.16000E-04,
     & 2.27270, 1.25550, 2.28000E-04,
     & 2.28310, 1.25450, 2.46000E-04,
     & 2.29360, 1.25340, 2.72000E-04,
     & 2.30410, 1.25240, 3.03000E-04,
     & 2.31480, 1.25130, 3.40000E-04,
     & 2.32560, 1.25020, 3.81000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 351, 360 ) /
     & 2.33640, 1.24900, 4.20000E-04,
     & 2.34740, 1.24780, 4.61000E-04,
     & 2.35850, 1.24650, 4.97000E-04,
     & 2.36970, 1.24510, 5.27000E-04,
     & 2.38100, 1.24380, 5.54000E-04,
     & 2.39230, 1.24240, 5.77000E-04,
     & 2.40380, 1.24080, 5.99000E-04,
     & 2.41550, 1.23930, 6.21000E-04,
     & 2.42720, 1.23770, 6.40000E-04,
     & 2.43900, 1.23590, 6.61000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 361, 370 ) /
     & 2.45100, 1.23410, 6.91000E-04,
     & 2.46310, 1.23220, 7.20000E-04,
     & 2.47520, 1.23010, 7.44000E-04,
     & 2.48760, 1.22790, 7.72000E-04,
     & 2.50000, 1.22580, 8.04000E-04,
     & 2.52000, 1.22198, 8.25500E-04,
     & 2.55000, 1.21548, 8.57800E-04,
     & 2.56500, 1.21184, 8.73900E-04,
     & 2.58000, 1.20790, 8.90000E-04,
     & 2.59000, 1.20507, 9.30000E-04 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 371, 380 ) /
     & 2.60000, 1.20209, 1.01000E-03,
     & 2.62000, 1.19566, 1.35000E-03,
     & 2.67500, 1.17411, 3.42000E-03,
     & 2.72500, 1.14734, 7.92000E-03,
     & 2.77800, 1.10766, 2.00000E-02,
     & 2.81700, 1.06739, 3.80000E-02,
     & 2.83300, 1.04762, 5.20000E-02,
     & 2.84900, 1.02650, 6.80000E-02,
     & 2.86500, 1.00357, 9.23000E-02,
     & 2.88200, 9.81970E-01, 1.27000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 381, 390 ) /
     & 2.89900, 9.65030E-01, 1.69000E-01,
     & 2.91500, 9.59620E-01, 2.21000E-01,
     & 2.93300, 9.72690E-01, 2.76000E-01,
     & 2.95000, 9.91720E-01, 3.12000E-01,
     & 2.96700, 1.00668, 3.47000E-01,
     & 2.98500, 1.02186, 3.88000E-01,
     & 3.00300, 1.04270, 4.38000E-01,
     & 3.02100, 1.07597, 4.93000E-01,
     & 3.04000, 1.12954, 5.54000E-01,
     & 3.05800, 1.21267, 6.12000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 391, 400 ) /
     & 3.07700, 1.32509, 6.25000E-01,
     & 3.09600, 1.42599, 5.93000E-01,
     & 3.11500, 1.49656, 5.39000E-01,
     & 3.13500, 1.55095, 4.91000E-01,
     & 3.15500, 1.59988, 4.38000E-01,
     & 3.17500, 1.63631, 3.72000E-01,
     & 3.19500, 1.65024, 3.00000E-01,
     & 3.21500, 1.64278, 2.38000E-01,
     & 3.23600, 1.62691, 1.93000E-01,
     & 3.25700, 1.61284, 1.58000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 401, 410 ) /
     & 3.27900, 1.59245, 1.21000E-01,
     & 3.30000, 1.57329, 1.03000E-01,
     & 3.32200, 1.55770, 8.36000E-02,
     & 3.34500, 1.54129, 6.68000E-02,
     & 3.36700, 1.52654, 5.40000E-02,
     & 3.39000, 1.51139, 4.22000E-02,
     & 3.41300, 1.49725, 3.42000E-02,
     & 3.43600, 1.48453, 2.74000E-02,
     & 3.46000, 1.47209, 2.20000E-02,
     & 3.48400, 1.46125, 1.86000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 411, 420 ) /
     & 3.50900, 1.45132, 1.52000E-02,
     & 3.53400, 1.44215, 1.26000E-02,
     & 3.55900, 1.43366, 1.06000E-02,
     & 3.62400, 1.41553, 8.02000E-03,
     & 3.73200, 1.39417, 6.85000E-03,
     & 3.77500, 1.38732, 6.60000E-03,
     & 3.84700, 1.37735, 6.96000E-03,
     & 3.96900, 1.36448, 9.16000E-03,
     & 4.09900, 1.35414, 1.11000E-02,
     & 4.23900, 1.34456, 1.45000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 421, 430 ) /
     & 4.34800, 1.33882, 2.00000E-02,
     & 4.38700, 1.33807, 2.30000E-02,
     & 4.44400, 1.33847, 2.60000E-02,
     & 4.50500, 1.34053, 2.90000E-02,
     & 4.54700, 1.34287, 2.93000E-02,
     & 4.56000, 1.34418, 3.00000E-02,
     & 4.58000, 1.34634, 2.85000E-02,
     & 4.71900, 1.34422, 1.73000E-02,
     & 4.90400, 1.33453, 1.29000E-02,
     & 5.00000, 1.32897, 1.20000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 431, 440 ) /
     & 5.10000, 1.32333, 1.25000E-02,
     & 5.20000, 1.31800, 1.34000E-02,
     & 5.26300, 1.31432, 1.40000E-02,
     & 5.40000, 1.30623, 1.75000E-02,
     & 5.55600, 1.29722, 2.40000E-02,
     & 5.71400, 1.28898, 3.50000E-02,
     & 5.74700, 1.28730, 3.80000E-02,
     & 5.78000, 1.28603, 4.20000E-02,
     & 5.81400, 1.28509, 4.60000E-02,
     & 5.84800, 1.28535, 5.20000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 441, 450 ) /
     & 5.88200, 1.28813, 5.70000E-02,
     & 6.06100, 1.30156, 6.90000E-02,
     & 6.13500, 1.30901, 7.00000E-02,
     & 6.25000, 1.31720, 6.70000E-02,
     & 6.28900, 1.31893, 6.50000E-02,
     & 6.32900, 1.32039, 6.40000E-02,
     & 6.36900, 1.32201, 6.20000E-02,
     & 6.41000, 1.32239, 5.90000E-02,
     & 6.45200, 1.32149, 5.70000E-02,
     & 6.49400, 1.32036, 5.60000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 451, 460 ) /
     & 6.57900, 1.31814, 5.50000E-02,
     & 6.66700, 1.31705, 5.70000E-02,
     & 6.75700, 1.31807, 5.80000E-02,
     & 6.89700, 1.31953, 5.70000E-02,
     & 7.04200, 1.31933, 5.50000E-02,
     & 7.14300, 1.31896, 5.50000E-02,
     & 7.24600, 1.31909, 5.40000E-02,
     & 7.35300, 1.31796, 5.20000E-02,
     & 7.46300, 1.31631, 5.20000E-02,
     & 7.57600, 1.31542, 5.20000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 461, 470 ) /
     & 7.69200, 1.31540, 5.20000E-02,
     & 7.81200, 1.31552, 5.00000E-02,
     & 7.93700, 1.31455, 4.70000E-02,
     & 8.06500, 1.31193, 4.30000E-02,
     & 8.19700, 1.30677, 3.90000E-02,
     & 8.33300, 1.29934, 3.70000E-02,
     & 8.47500, 1.29253, 3.90000E-02,
     & 8.69600, 1.28389, 4.00000E-02,
     & 8.92900, 1.27401, 4.20000E-02,
     & 9.09100, 1.26724, 4.40000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 471, 480 ) /
     & 9.25900, 1.25990, 4.50000E-02,
     & 9.52400, 1.24510, 4.60000E-02,
     & 9.80400, 1.22241, 4.70000E-02,
     & 1.00000E+01, 1.19913, 5.10000E-02,
     & 1.02000E+01, 1.17150, 6.50000E-02,
     & 1.03100E+01, 1.15528, 7.50000E-02,
     & 1.04200E+01, 1.13700, 8.80000E-02,
     & 1.05300E+01, 1.11808, 1.08000E-01,
     & 1.06400E+01, 1.10134, 1.34000E-01,
     & 1.07500E+01, 1.09083, 1.68000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 481, 490 ) /
     & 1.08700E+01, 1.08734, 2.04000E-01,
     & 1.10000E+01, 1.09254, 2.48000E-01,
     & 1.11100E+01, 1.10654, 2.80000E-01,
     & 1.13600E+01, 1.14779, 3.41000E-01,
     & 1.16300E+01, 1.20202, 3.79000E-01,
     & 1.19000E+01, 1.25825, 4.09000E-01,
     & 1.22000E+01, 1.32305, 4.22000E-01,
     & 1.25000E+01, 1.38574, 4.22000E-01,
     & 1.28200E+01, 1.44478, 4.03000E-01,
     & 1.29900E+01, 1.47170, 3.89000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 491, 500 ) /
     & 1.31600E+01, 1.49619, 3.74000E-01,
     & 1.33300E+01, 1.51652, 3.54000E-01,
     & 1.35100E+01, 1.53328, 3.35000E-01,
     & 1.37000E+01, 1.54900, 3.15000E-01,
     & 1.38900E+01, 1.56276, 2.94000E-01,
     & 1.40800E+01, 1.57317, 2.71000E-01,
     & 1.42900E+01, 1.58028, 2.46000E-01,
     & 1.47100E+01, 1.57918, 1.98000E-01,
     & 1.51500E+01, 1.56672, 1.64000E-01,
     & 1.53800E+01, 1.55869, 1.52000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 501, 510 ) /
     & 1.56300E+01, 1.55081, 1.42000E-01,
     & 1.61300E+01, 1.53807, 1.28000E-01,
     & 1.63900E+01, 1.53296, 1.25000E-01,
     & 1.66700E+01, 1.53220, 1.23000E-01,
     & 1.69500E+01, 1.53340, 1.16000E-01,
     & 1.72400E+01, 1.53289, 1.07000E-01,
     & 1.81800E+01, 1.51705, 7.90000E-02,
     & 1.88700E+01, 1.50097, 7.20000E-02,
     & 1.92300E+01, 1.49681, 7.60000E-02,
     & 1.96100E+01, 1.49928, 7.50000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 511, 520 ) /
     & 2.00000E+01, 1.50153, 6.70000E-02,
     & 2.04100E+01, 1.49856, 5.50000E-02,
     & 2.08300E+01, 1.49053, 4.50000E-02,
     & 2.22200E+01, 1.46070, 2.90000E-02,
     & 2.26000E+01, 1.45182, 2.75000E-02,
     & 2.30500E+01, 1.44223, 2.70000E-02,
     & 2.36000E+01, 1.43158, 2.73000E-02,
     & 2.46000E+01, 1.41385, 2.89000E-02,
     & 2.50000E+01, 1.40676, 3.00000E-02,
     & 2.60000E+01, 1.38955, 3.40000E-02 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 521, 530 ) /
     & 2.85700E+01, 1.34894, 5.30000E-02,
     & 3.10000E+01, 1.31039, 7.55000E-02,
     & 3.33300E+01, 1.26420, 1.06000E-01,
     & 3.44800E+01, 1.23656, 1.35000E-01,
     & 3.56400E+01, 1.21663, 1.76100E-01,
     & 3.70000E+01, 1.20233, 2.22900E-01,
     & 3.82400E+01, 1.19640, 2.74600E-01,
     & 3.96000E+01, 1.19969, 3.28000E-01,
     & 4.11400E+01, 1.20860, 3.90600E-01,
     & 4.27600E+01, 1.22173, 4.64200E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 531, 540 ) /
     & 4.35800E+01, 1.24166, 5.24700E-01,
     & 4.45800E+01, 1.28175, 5.73100E-01,
     & 4.55000E+01, 1.32784, 6.36200E-01,
     & 4.61500E+01, 1.38657, 6.83900E-01,
     & 4.67100E+01, 1.46486, 7.09100E-01,
     & 4.73600E+01, 1.55323, 6.79000E-01,
     & 4.80000E+01, 1.60379, 6.25000E-01,
     & 4.87800E+01, 1.61877, 5.65400E-01,
     & 5.00300E+01, 1.62963, 5.43300E-01,
     & 5.12800E+01, 1.65712, 5.29200E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 541, 550 ) /
     & 5.27500E+01, 1.69810, 5.07000E-01,
     & 5.35000E+01, 1.72065, 4.88300E-01,
     & 5.42400E+01, 1.74865, 4.70700E-01,
     & 5.50000E+01, 1.76736, 4.20300E-01,
     & 5.57400E+01, 1.76476, 3.77100E-01,
     & 5.64000E+01, 1.75011, 3.37600E-01,
     & 5.70000E+01, 1.72327, 3.05600E-01,
     & 5.74600E+01, 1.68490, 2.83500E-01,
     & 5.84000E+01, 1.62398, 3.17000E-01,
     & 5.92900E+01, 1.59596, 3.51700E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 551, 560 ) /
     & 6.00000E+01, 1.58514, 3.90200E-01,
     & 6.10000E+01, 1.59917, 4.50900E-01,
     & 6.12500E+01, 1.61405, 4.67100E-01,
     & 6.25000E+01, 1.66625, 4.77900E-01,
     & 6.37800E+01, 1.70663, 4.89000E-01,
     & 6.46700E+01, 1.73713, 4.89900E-01,
     & 6.55800E+01, 1.76860, 4.87300E-01,
     & 6.65500E+01, 1.80343, 4.76600E-01,
     & 6.76000E+01, 1.83296, 4.50800E-01,
     & 6.90000E+01, 1.85682, 4.19300E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 561, 570 ) /
     & 7.05300E+01, 1.87411, 3.88000E-01,
     & 7.30000E+01, 1.89110, 3.43300E-01,
     & 7.50000E+01, 1.89918, 3.11800E-01,
     & 7.62900E+01, 1.90432, 2.93500E-01,
     & 8.00000E+01, 1.90329, 2.35000E-01,
     & 8.29700E+01, 1.88744, 1.98100E-01,
     & 8.50000E+01, 1.87499, 1.86500E-01,
     & 8.68000E+01, 1.86702, 1.77100E-01,
     & 9.08000E+01, 1.85361, 1.62000E-01,
     & 9.51700E+01, 1.84250, 1.49000E-01 /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 571, 574 ) /
     & 1.00000E+02, 1.83225, 1.39000E-01,
     & 1.20000E+02, 1.81914, 1.20000E-01,
     & 1.50000E+02, 1.82268, 9.62000E-02,
     & 1.67000E+02, 1.82961, 8.30000E-02 /


      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 1, 10 ) /
     & 1.67000E+02, 1.82961, 8.30000E-02, 1.82961, 8.30000E-02,
     & 1.77800E+02, 1.83258, 6.90000E-02, 1.83258, 6.90000E-02,
     & 1.88400E+02, 1.83149, 5.70000E-02, 1.83149, 5.70000E-02,
     & 1.99500E+02, 1.82748, 4.56000E-02, 1.82748, 4.56000E-02,
     & 2.11300E+02, 1.82224, 3.79000E-02, 1.82224, 3.79000E-02,
     & 2.23900E+02, 1.81718, 3.14000E-02, 1.81718, 3.14000E-02,
     & 2.37100E+02, 1.81204, 2.62000E-02, 1.81204, 2.62000E-02,
     & 2.51200E+02, 1.80704, 2.24000E-02, 1.80704, 2.24000E-02,
     & 2.66100E+02, 1.80250, 1.96000E-02, 1.80250, 1.96000E-02,
     & 2.81800E+02, 1.79834, 1.76000E-02, 1.79834, 1.76000E-02 /

      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 11, 20 ) /
     & 2.98500E+02, 1.79482, 1.66500E-02, 1.79482, 1.66500E-02,
     & 3.16200E+02, 1.79214, 1.62000E-02, 1.79214, 1.60000E-02,
     & 3.54800E+02, 1.78843, 1.55000E-02, 1.78843, 1.50000E-02,
     & 3.98100E+02, 1.78601, 1.47000E-02, 1.78601, 1.40000E-02,
     & 4.46700E+02, 1.78434, 1.39000E-02, 1.78434, 1.31000E-02,
     & 5.01200E+02, 1.78322, 1.32000E-02, 1.78322, 1.23000E-02,
     & 5.62300E+02, 1.78248, 1.25000E-02, 1.78248, 1.15000E-02,
     & 6.31000E+02, 1.78201, 1.18000E-02, 1.78201, 1.08000E-02,
     & 7.94300E+02, 1.78170, 1.06000E-02, 1.78170, 9.46000E-03,
     & 1.00000E+03, 1.78160, 9.54000E-03, 1.78160, 8.29000E-03 /

      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 21, 30 ) /
     & 1.25900E+03, 1.78190, 8.56000E-03, 1.78190, 7.27000E-03,
     & 2.50000E+03, 1.78300, 6.21000E-03, 1.78300, 4.91000E-03,
     & 5.00000E+03, 1.78430, 4.49000E-03, 1.78430, 3.30000E-03,
     & 1.00000E+04, 1.78520, 3.24000E-03, 1.78520, 2.22000E-03,
     & 2.00000E+04, 1.78620, 2.34000E-03, 1.78610, 1.49000E-03,
     & 3.20000E+04, 1.78660, 1.88000E-03, 1.78630, 1.14000E-03,
     & 3.50000E+04, 1.78680, 1.74000E-03, 1.78640, 1.06000E-03,
     & 4.00000E+04, 1.78690, 1.50000E-03, 1.78650, 9.48000E-04,
     & 4.50000E+04, 1.78700, 1.32000E-03, 1.78650, 8.50000E-04,
     & 5.00000E+04, 1.78700, 1.16000E-03, 1.78650, 7.66000E-04 /

      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 31, 40 ) /
     & 6.00000E+04, 1.78710, 8.80000E-04, 1.78650, 6.30000E-04,
     & 7.00000E+04, 1.78710, 6.95000E-04, 1.78650, 5.20000E-04,
     & 9.00000E+04, 1.78720, 4.64000E-04, 1.78650, 3.84000E-04,
     & 1.11000E+05, 1.78720, 3.40000E-04, 1.78650, 2.96000E-04,
     & 1.20000E+05, 1.78720, 3.11000E-04, 1.78650, 2.70000E-04,
     & 1.30000E+05, 1.78720, 2.94000E-04, 1.78650, 2.52000E-04,
     & 1.40000E+05, 1.78720, 2.79000E-04, 1.78650, 2.44000E-04,
     & 1.50000E+05, 1.78720, 2.70000E-04, 1.78650, 2.36000E-04,
     & 1.60000E+05, 1.78720, 2.64000E-04, 1.78650, 2.30000E-04,
     & 1.70000E+05, 1.78720, 2.58000E-04, 1.78650, 2.28000E-04 /

      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 41, 50 ) /
     & 1.80000E+05, 1.78720, 2.52000E-04, 1.78650, 2.25000E-04,
     & 2.00000E+05, 1.78720, 2.49000E-04, 1.78650, 2.20000E-04,
     & 2.50000E+05, 1.78720, 2.54000E-04, 1.78650, 2.16000E-04,
     & 2.90000E+05, 1.78720, 2.64000E-04, 1.78650, 2.17000E-04,
     & 3.20000E+05, 1.78720, 2.74000E-04, 1.78650, 2.20000E-04,
     & 3.50000E+05, 1.78720, 2.89000E-04, 1.78650, 2.25000E-04,
     & 3.80000E+05, 1.78720, 3.05000E-04, 1.78650, 2.32000E-04,
     & 4.00000E+05, 1.78720, 3.15000E-04, 1.78650, 2.39000E-04,
     & 4.50000E+05, 1.78720, 3.46000E-04, 1.78650, 2.60000E-04,
     & 5.00000E+05, 1.78720, 3.82000E-04, 1.78650, 2.86000E-04 /

      DATA  ( WLT( I ), TABRET( I, 1 ), TABIMT( I, 1 ),
     &                  TABRET( I, 2 ), TABIMT( I, 2 ), I = 51, 62 ) /
     & 6.00000E+05, 1.78720, 4.62000E-04, 1.78650, 3.56000E-04,
     & 6.40000E+05, 1.78720, 5.00000E-04, 1.78650, 3.83000E-04,
     & 6.80000E+05, 1.78720, 5.50000E-04, 1.78650, 4.15000E-04,
     & 7.20000E+05, 1.78720, 5.95000E-04, 1.78650, 4.45000E-04,
     & 7.60000E+05, 1.78720, 6.47000E-04, 1.78650, 4.76000E-04,
     & 8.00000E+05, 1.78720, 6.92000E-04, 1.78650, 5.08000E-04,
     & 8.40000E+05, 1.78720, 7.42000E-04, 1.78650, 5.40000E-04,
     & 9.00000E+05, 1.78720, 8.20000E-04, 1.78650, 5.86000E-04,
     & 1.00000E+06, 1.78720, 9.70000E-04, 1.78650, 6.78000E-04,
     & 2.00000E+06, 1.78720, 1.95000E-03, 1.78650, 1.28000E-03,
     & 5.00000E+06, 1.78720, 5.78000E-03, 1.78650, 3.55000E-03,
     & 8.60000E+06, 1.78800, 9.70000E-03, 1.78720, 5.60000E-03 /


      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 1, 10 ) /
     & 1.82961, 8.30000E-02, 1.82961, 8.30000E-02,
     & 1.83258, 6.90000E-02, 1.83258, 6.90000E-02,
     & 1.83149, 5.70000E-02, 1.83149, 5.70000E-02,
     & 1.82748, 4.56000E-02, 1.82748, 4.45000E-02,
     & 1.82224, 3.79000E-02, 1.82224, 3.55000E-02,
     & 1.81718, 3.14000E-02, 1.81718, 2.91000E-02,
     & 1.81204, 2.62000E-02, 1.81204, 2.44000E-02,
     & 1.80704, 2.19000E-02, 1.80704, 1.97000E-02,
     & 1.80250, 1.88000E-02, 1.80250, 1.67000E-02,
     & 1.79834, 1.66000E-02, 1.79834, 1.40000E-02 /

      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 11, 20 ) /
     & 1.79482, 1.54000E-02, 1.79482, 1.23500E-02,
     & 1.79214, 1.47000E-02, 1.79214, 1.08000E-02,
     & 1.78843, 1.35000E-02, 1.78843, 8.90000E-03,
     & 1.78601, 1.25000E-02, 1.78601, 7.34000E-03,
     & 1.78434, 1.15000E-02, 1.78434, 6.40000E-03,
     & 1.78322, 1.06000E-02, 1.78322, 5.60000E-03,
     & 1.78248, 9.77000E-03, 1.78248, 5.00000E-03,
     & 1.78201, 9.01000E-03, 1.78201, 4.52000E-03,
     & 1.78160, 7.66000E-03, 1.78150, 3.68000E-03,
     & 1.78140, 6.52000E-03, 1.78070, 2.99000E-03 /

      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 21, 30 ) /
     & 1.78160, 5.54000E-03, 1.78010, 2.49000E-03,
     & 1.78220, 3.42000E-03, 1.77890, 1.55000E-03,
     & 1.78310, 2.10000E-03, 1.77790, 9.61000E-04,
     & 1.78380, 1.29000E-03, 1.77730, 5.95000E-04,
     & 1.78390, 7.93000E-04, 1.77720, 3.69000E-04,
     & 1.78400, 5.70000E-04, 1.77720, 2.67000E-04,
     & 1.78400, 5.35000E-04, 1.77720, 2.51000E-04,
     & 1.78400, 4.82000E-04, 1.77720, 2.29000E-04,
     & 1.78400, 4.38000E-04, 1.77720, 2.11000E-04,
     & 1.78400, 4.08000E-04, 1.77720, 1.96000E-04 /

      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 31, 40 ) /
     & 1.78390, 3.50000E-04, 1.77720, 1.73000E-04,
     & 1.78380, 3.20000E-04, 1.77720, 1.55000E-04,
     & 1.78370, 2.55000E-04, 1.77720, 1.31000E-04,
     & 1.78370, 2.12000E-04, 1.77720, 1.13000E-04,
     & 1.78370, 2.00000E-04, 1.77720, 1.06000E-04,
     & 1.78370, 1.86000E-04, 1.77720, 9.90000E-05,
     & 1.78370, 1.75000E-04, 1.77720, 9.30000E-05,
     & 1.78370, 1.66000E-04, 1.77720, 8.73000E-05,
     & 1.78370, 1.56000E-04, 1.77720, 8.30000E-05,
     & 1.78370, 1.49000E-04, 1.77720, 7.87000E-05 /

      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 41, 50 ) /
     & 1.78370, 1.44000E-04, 1.77720, 7.50000E-05,
     & 1.78370, 1.35000E-04, 1.77720, 6.83000E-05,
     & 1.78370, 1.21000E-04, 1.77720, 5.60000E-05,
     & 1.78370, 1.16000E-04, 1.77720, 4.96000E-05,
     & 1.78370, 1.16000E-04, 1.77720, 4.55000E-05,
     & 1.78370, 1.17000E-04, 1.77720, 4.21000E-05,
     & 1.78370, 1.20000E-04, 1.77720, 3.91000E-05,
     & 1.78370, 1.23000E-04, 1.77720, 3.76000E-05,
     & 1.78370, 1.32000E-04, 1.77720, 3.40000E-05,
     & 1.78370, 1.44000E-04, 1.77720, 3.10000E-05 /

      DATA ( TABRET( I, 3 ), TABIMT( I, 3 ),
     &       TABRET( I, 4 ), TABIMT( I, 4 ), I = 51, 62 ) /
     & 1.78370, 1.68000E-04, 1.77720, 2.64000E-05,
     & 1.78370, 1.80000E-04, 1.77720, 2.51000E-05,
     & 1.78370, 1.90000E-04, 1.77720, 2.43000E-05,
     & 1.78370, 2.09000E-04, 1.77720, 2.39000E-05,
     & 1.78370, 2.16000E-04, 1.77720, 2.37000E-05,
     & 1.78370, 2.29000E-04, 1.77720, 2.38000E-05,
     & 1.78370, 2.40000E-04, 1.77720, 2.40000E-05,
     & 1.78370, 2.60000E-04, 1.77720, 2.46000E-05,
     & 1.78370, 2.92000E-04, 1.77720, 2.66000E-05,
     & 1.78370, 6.10000E-04, 1.77720, 4.45000E-05,
     & 1.78400, 1.02000E-03, 1.77720, 8.70000E-05,
     & 1.78450, 1.81000E-03, 1.77800, 1.32000E-04 /

      END

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific).
c        On Cray, symbolic dump only works if have compiled with
c        z compiler option which creates debug symbol table.

      LOGICAL       FATAL, MsgLim, Cray
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /,
     &     Cray / .False. /


      IF ( FATAL )  THEN
         WRITE (*,'(//,2A,//)')  ' ******* ERROR >>>>>>  ', MESSAG
c         IF( Cray )  CALL  SYMDUMP()
c         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE (*,'(/,2A,/)')  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         MsgLim = .True.
         WRITE (*,'(//,A,//)') ' >>>>>> TOO MANY WARNING MESSAGES --  '
     &                      //'They will no longer be printed  <<<<<<<'
      ENDIF
c     convert back from micro to meter

      RETURN
      END

      subroutine abs_coef(lamda, ref_index, beta)
c calculate absorption coefficient
c INPUT:
c   lamda:wavelength (m)
c   ref_index: complex refractive index
c OUTPUT:
c   beta: absorption coefficient (1/m)
      real*8 beta, lamda,pi,im
      complex*16 ref_index
      pi=dacos(-1.0d0)
      im = dimag(ref_index)
      beta = 4.0*pi*im/lamda
      return
      end
      
      subroutine path_length(radius,in_angle,ref_index,xi)
c PURPOSE: To calculate path length between two consecutive incident interfaces
c INPUT:
c     radius: radius of the sphere
c     in_angle: incident angle in radians
c     ref_index: complex refractive index
c OUTPUT:
c     xi: path length bwtween two consecutive incident points in m
c LOCAL:
c     re: real part of ref_index
c     im: imaginary part of ref_index

      real*8 radius,in_angle
      complex*16 ref_index
      real*8 xi
      real*8 re,im
      real*8 Trans, Reflec,Trp,Trv
      complex*16 Rep,Rev,tp,tv,rp,rv,t,r
      real*8 u,v

      re = dreal(ref_index)
      im = dimag(ref_index)
      call TransReflec(in_angle,ref_index,Trans,Reflec,Trp,Trv,
     &                                     Rep,Rev,tp,tv,rp,rv,t,r,u,v)
      xi = 2.*radius/(re*re+im*im)*dsqrt((u*re+v*im)**2+(v*re-u*im)**2)
      return
      end
      
      subroutine sumNumber1(beta,xi,r,eps,sumNumber)
c PURPOSE: 
c     To calculate the maximum number of summation for the 
c     calculation of absorption efficiency in ray tracing method. The truncation
c     is based on the ratio of incident electric field amplitude at Nth interface
c     to that transmitted at the first interface.
c INPUT VARIABLES:
c     beta: the absorption coefficient, from function abs_coef(lamda, ref_index)
c     xi: path length between two consecutive incident interfaces, from function 
c         path_length(radius,in_angle,ref_index)
c     r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2)), from
c        subroutine TransReflec(in_angle,ref_index,Trans,Reflec,tp,tv,rp,rv,t,r)
c     eps: real, tolerance for cutoff in the summation series.
c OUTPUT VARIALEs:
c     sumNumber: integer, maximum number of summation for absorption efficiency, etc.
c LOCAL VARIABLES:
c     tem: real
      real*8 beta,xi, r,eps
      integer*4 sumNumber
      real*8 tem
      tem = 0.0d0
      tem = beta*xi*0.5-2.*dlog(r)-dlog(eps)
      tem = tem /( beta*xi*0.5-dlog(r) )
c conversion from double to integer*4      
      sumNumber = idint(tem)
c at least once path is made once radiation transmitted into the sphere
      if (sumNumber.lt.2) then
         sumNumber = 2
      end if
      return
      end
          
      subroutine sumNumber2(beta,xi,r,t,eps,sumNumber)
c PURPOSE: 
c     To calculate the maximum number of summation for the 
c     calculation of absorption efficiency in ray tracing method. The truncation
c     is based on the ratio of incident electric field amplitude at Nth inteface
c     to that incident at the first inteface.
c INPUT VARIABLES:
c     beta: the absorption coefficient, from function abs_coef(lamda, ref_index)
c     xi: path length between two consecutive incident interfaces, from function 
c         path_length(radius,in_angle,ref_index)
c     r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2)), from
c        subroutine TransReflec(in_angle,ref_index,Trans,Reflec,tp,tv,rp,rv,t,r)
c     t: transmission coefficient for natural light r=sqrt(0.5(rp*2+rv**2)), from
c        subroutine TransReflec(in_angle,ref_index,Trans,Reflec,tp,tv,rp,rv,t,r)
c     eps: real, tolerance for cutoff in the summation series.
c OUTPUT VARIALEs:
c     sumNumber2: integer, maximum number of summation for absorption efficiency, etc.
c LOCAL VARIABLES:
c     tem: real
      real*8 beta,xi, r,t,eps
      integer*4 sumNumber
      real*8 tem
      tem = 0.0d0
      tem = beta*xi*0.5 + dlog(t) -2.*dlog(r)-dlog(eps)
      tem = tem /( beta*xi*0.5-dlog(r) )
c conversion from double to integer*4      
      sumNumber = idint(tem)
c at least once path is made once radiation transmitted into the sphere
      if ( sumNumber.lt.2) sumNumber = 2
      return
      end


      subroutine sumNumber3(beta,xi,t,Reflec,Trans,eps,sumNumber)
c PURPOSE: 
c     To calculate the maximum number of summation for the 
c     calculation of absorption efficiency in ray tracing method. The truncation
c     is based on the ratio of incident irradiance at Nth interface
c     to that incident at the first interface.
c INPUT VARIABLES:
c     beta: the absorption coefficient, from function abs_coef(lamda, ref_index)
c     xi: path length between two consecutive incident interfaces, from function 
c         path_length(radius,in_angle,ref_index)
c     Reflec: reflection for natural light , from
c        subroutine TransReflec(in_angle,ref_index,Trans,Reflec,tp,tv,rp,rv,t,r)
c     Trans: transmission for natural light , from
c        subroutine TransReflec(in_angle,ref_index,Trans,Reflec,tp,tv,rp,rv,t,r)
c     eps: real, tolerance for cutoff in the summation series.
c OUTPUT VARIALEs:
c     sumNumber2: integer, maximum number of summation for absorption efficiency, etc.
c LOCAL VARIABLES:
c     tem: real

      real*8 beta,xi,Reflec,Trans,eps
      complex*16 t
      integer*4 sumNumber
      real*8 tem
      tem = 0.0d0
      tem = beta*xi + 2.0*dlog(cdabs(t)) -2.*dlog(Reflec)-dlog(eps)
      tem = tem /( beta*xi-dlog(Reflec) )
c conversion from double to integer*4      
      sumNumber = idint(tem)
c at least once path is made once radiation transmitted into the sphere
      if (sumNumber.lt.2) sumNumber = 2
      return
      end

      subroutine RefractionAngle(ref_index,in_angle,thetat)
c Purpose: 
c   To calculate complex refractive angle thetat 
cInput: 
c   ref_index: complex refractive index 
c   in_angle:incident angle 
c Output: 
c   thetat: complex refraction angle
      complex*16 ref_index,thetat,thetat1, thetat2,z,uniti
      real*8 in_angle,pi, re,im
      uniti=(0.0d0,1.0d0)
      pi=dacos(-1.0d0)
      re=dreal(ref_index)
      im=dimag(ref_index)
      z=complex(re,-im)*dsin(in_angle)/(re*re+im*im)
      call cdasin(z,thetat1,thetat2)
C According to physical processes, real refraction angle should be < incident angle.
      if(cdabs(thetat1).le.in_angle) then
         thetat=thetat1
      else
         thetat=thetat2
      end if
      return
      end

      subroutine cdacos(z,z1,z2)
c Purpose:
c   To calculate the inverse trigonometric function of cosine for complex*16 argument
c Input
c   z: complex*16 argument
c Output
c   z1,z2: the results of cos(z)
       complex*16 z, z1, z2, uniti
       uniti=(0.0d0,1.0d0)
       z1=-uniti*cdlog(z+uniti*cdsqrt(1.0-z*z))
       z2=-uniti*cdlog(z-uniti*cdsqrt(1.0-z*z))
       return
       end
       
      subroutine cdasin(z,z1,z2)
c Purpose:
c   To calculate the inverse trigonometric function of sine for complex*16 argument
c   arcsin(z).
c Input
c   z: complex*16 argument
c Output
c   z1,z2: the results of sin(z)
       complex*16 z, z1, z2, uniti
       uniti=(0.0d0,1.0d0)
       z1=-uniti*cdlog(uniti*z+cdsqrt(1.0-z*z))
       z2=-uniti*cdlog(uniti*z-cdsqrt(1.0-z*z))
       return
       end


       subroutine quadrature_gauss(x,w,n,a,b,G)
c  PURPOSE:
c      To calculate the integral G of the function f(x) between a and b using the 
c      abscissas and weights of the Gauss-Legendre n-point quadrature calculated
c      from subroutine GaussLeg(n,x,w). Function FF(x) is to be writen by users.
c  REFERENCE: Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
c             Orlando, FL: Academic Press,1984. p.481-483.
c  INPUT VARIABLES:
c      a: the lower limit of integration
c      b: the upper limit of integration
c      n: number of points
c      x[]: array of n elements of abscissas
c      w[]: arrayof n elements of weights
c  OUTPUT VARIABLES:
c       G: the integral value
c  LOCAL VARIABLES:
c      i: loop varible
c      xm,xk: real*8
       parameter(max_points=50)
       real*8 a, b, G, FF
       integer*4 n,i,m,k
       real*8 x(max_points),w(max_points),xm,xk,dx
       xm=0.5d0*(b+a)
       xk=0.5d0*(b-a)
       m=(n+1)/2
       k=2*m-n+1
       G=0.0d0
       do i=k,m
          dx=xk*x(i)
          G=G+w(i)*(FF(xm+dx)+FF(xm-dx))
       end do
       G=G*xk
       return
       end
       
       real*8 function FF(x)
c PURPOSE:
c     The user-defiend function
c INPUT:
c     x: variable
c OUTPUT:
c     FF: value of teh function
      real*8 x
      FF=x*x*x+2
      return
      end

       subroutine GaussLeg(n,x,w)
c  PURPOSE:
c      To calculate the abscissas and weights of the Gauss-Legendre n-point quadrature
c      formula. Integration of f(x) is from -1 to 1.
c  Reference: Press W H , S A Teukolsky, W T Vetterling and B P Flannery,"Numerical 
c             recipes in C: the art of scientific computing, Second edition, Cambridge
c             University Press,1992. p.147-152.
c  INPUT VARIABLES:
c      n: number of points
c  OUTPUT VARIABLES:
c      x[]: array of n elements of abscissas
c      w[]: arrayof n elements of weights
c  LOCAL VARIABLES:
c      eps: relative pressision, = 1.0d-10
c      z1,z,xm,xk,
c      pp: the derivative of Legendre polynomial
c      p1: Legendre polynomial
c      p2,p3: real*8
c      m: only half of the number of roots is necessary to be found due to that
c         the roots are symmetric in the interval.
c      j,i:dummy loop variables
c      
       real*8 x1, x2
       integer*4 n
       real*8 x(n),w(n)
       integer*4 m,j,i
       real*8 eps,z1,z,xm,xk,pp,p1,p2,p3
       eps=1.0d-10
       m=(n+1)/2
       x1=-1.0d0
       x2=1.0d0
       xm=0.5d0*(x2+x1)
       xk=0.5d0*(x2-x1)
       pi=dacos(-1.0d0)
       z1=0.0d0
       do i = 1,m
          z=dcos(pi*(dfloat(i)-0.25d0)/(dfloat(n)+0.5d0))
          do while (dabs(z-z1).gt.eps)
             p1=1.0d0
             p2=0.0d0
             do j=1,n
                p3=p2
                p2=p1
                p1=((2.0d0*dfloat(j)-1.0d0)*z*p2
     &               -(dfloat(j)-1.0d0)*p3)/dfloat(j)
             end do
             pp=dfloat(n)*(z*p1-p2)/(z*z-1.0d0)
             z1=z
             z=z1-p1/pp
          end do
          x(i)=xm-xk*z
          x(n+1-i)=xm+xk*z
          w(i)=2.0*xk/((1.0d0-z*z)*pp*pp)
          w(n+1-i)=w(i)
       end do
       return
       end
       
      subroutine TransReflec(in_angle,ref_index,Trans,Reflec,Trp,Trv,
     &                                     Rep,Rev,tp,tv,rp,rv,t,r,u,v)
c****************************************************************************
c       program TR 
cc to test  subroutine TransReflec(...)
c       complex*16 ref_index
c       real*8 in_angle,Trans,Reflec,Trp,Trv,Rep,Rev,tp,tv,rp,rv,t,r
c       real*8 re,im,thetat,pi
c       ref_index =(2.0d0,0.5d0) 
c       pi = dacos(-1.0d0)
c       in_angle = pi/3.
c       call TransReflec(in_angle,ref_index,Trans,Reflec,Trp,Trv
c     &                                     Rep,Rev,tp,tv,rp,rv,t,r)
c       write(*,*)'rp=', rp
c       write(*,*)'R=',Reflec
c       stop 
c       end
c****************************************************************************
c PURPOSE: To calculate the various transmission and reflection coefficients
c INPUT:
c    in_angle: incident angle
c    ref_index: complex refractive index
c OUTPUT:
c    Trans: transmission for natural light
c    Reflec: reflectance for natural light
c    Trp: transmission for parallel-polarized light
c    Trv: transmission for vertical-polarized light
c    Rep: reflection for parallel-polarized light
c    Rev: reflection for vertical-polarized light
c     tp: real, amplitude coefficient of transmission for parallel-polarization
c     tv: real, amplitude coefficient of transmission for vertical-polarization
c     rp: real, amplitude coefficient of reflection for parallel-polarization
c     rv: real, amplitude coefficient of reflection for vertical-polarization
c     t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
c     r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
c     u: mcos(theta)=u+iv
c     v:
c LOCAL:
c     thetat: real, refraction angle in redians
      real*8 in_angle,Trans, Reflec,Trp,Trv
      complex*16 Rep,Rev,tp,tv,rp,rv,t,r,thetat
      complex*16 ref_index, uniti
      real*8 re,im,pi,u,v
      real*8 delta, gama

      pi = dacos(-1.0d0)
      uniti=(0.0d0,1.0d0)
      re = dreal(ref_index)
      im = dimag(ref_index)
      delta=re*re-im*im-(dsin(in_angle))**2
      gama=re*im
      u=dsqrt(2.0d0)/2.0*dsqrt(dsqrt(delta**2+4.0*gama**2)+delta)
      v=dsqrt(2.0d0)/2.0*dsqrt(dsqrt(delta**2+4.0*gama**2)-delta)
      tp=2.*(re+uniti*im)*dcos(in_angle)
     &   /(((re*re-im*im)*dcos(in_angle)+u)
     &      +uniti*(2.*re*im*dcos(in_angle)+v))

      tv = 2.*dcos(in_angle)/(dcos(in_angle)+u+uniti*v)

      rp = (((re*re-im*im)*dcos(in_angle)-u)
     &      +uniti*(2.*re*im*dcos(in_angle)-v))
     &    /(((re*re-im*im)*dcos(in_angle)+u)
     &      +uniti*(2.*re*im*dcos(in_angle)+v))

      rv = (dcos(in_angle)-u-uniti*v)/(dcos(in_angle)+u+uniti*v)

      t = cdsqrt( 0.5 * ( tp * tp + tv * tv ) )
      r = cdsqrt( 0.5 * ( rp * rp + rv * rv ) )
      
      Trp = re*dsqrt((u*re+v*im)**2+(v*re-u*im)**2)
     &      /((re*re+im*im)*dcos(in_angle))*cdabs(tp)*cdabs(tp)
      Trv = re*dsqrt((u*re+v*im)**2+(v*re-u*im)**2)
     &      /((re*re+im*im)*dcos(in_angle))*cdabs(tv)*cdabs(tv)
     
      Rep = rp*rp
      Rev = rv*rv
      Trans = re*dsqrt((u*re+v*im)**2+(v*re-u*im)**2)
     &      /((re*re+im*im)*dcos(in_angle))*cdabs(t)*cdabs(t)

      Reflec =cdabs(r)*cdabs(r)
      return
      end

       subroutine roots(j, m, theta, x1, x2, xacc, nsegments, nb, rts)
c   Purpose:
c      To find all possible roots of a function fx(x) between x1 and x2. First, bracket 
c      the roots using subroutine brackets(j,m,theta,x1, x2, nsegments, xb1, xb2, nb), then utilize a hybrid
c      algorithm of bisection and Newton-Raphson methds within each bracketed interval, using
c      subroutine rootfind(j, m, theta, x1, x2, xacc, rt,flag).
c  Input:
c      j: jth interface. 
c      m: complex refractive index.
c      theta: outgoing angle (scattering angle).
c      x1,x2: limits that all the roots are to be found. 
c      xacc: accuracy that is used to find the root.
c      nb: as input, it is the maximum number of roots sought
c      nsegments: number that the inteval [x1,x2] is equally subdivided in order to bracket = maximum
c         possible roots
c  Output:
c      nb: as output, it is the number of different roots that are found.
c      rts: the different roots found that lies within [x1,x2]
c  Local:
c      xb1(nb),xb2(nb): bracketed pairs within which the root is to be found
c      i: integer for loop
c      rt: the root found that lies within [xb1(i), xb2(i)]
c  Subroutines:
c      brackets(j,m,theta,x1, x2, nsegments, xb1, xb2, nb): to bracket all possible roots
c      rootfind(j, m, theta, a, b, xacc, rt,flag): to find each root within the bracketed interval [a,b]
c  References:
c      None
       
       real*8 theta
       integer*4 j,nm,nb,nsegments
       complex*16 m
       real*8 x1, x2, xacc, rts(nb)
       real*8 xb1(1000), xb2(1000), rt
       integer*4 i,ll
       integer*4 nmax,flag
       
       nmax=100
c  To enlarge the number of segments in [x1,x2] to find at least one root. If the number
c  exceed a maximum number of segments, no root exists, exit.
       nb=nmax
c       n=nmax
       call brackets(j,m,theta,x1, x2, nsegments, xb1, xb2, nb)
c  nb is changed to be the real number of possible roots after calling  brackets(..)  
         n = nsegments  
       do while (nb.eq.0.and.n.le.nmax)
          n=n+1
          nb=nmax
          call brackets(j,m,theta,x1, x2, n, xb1, xb2, nb)
       end do
c If there is no root
       if (nb.eq.0.or.n.gt.nmax) then
          nb=0
          return
       end if
       
c       call rootfind(j, m, theta, xb1(1), xb2(1), xacc, rt,flag)
c       if (flag.eq.1) return    
           
c       rts(1)=rt
c       ll=1
       ll=0
       if (nb.ge.1) then
          do 10 i = 1, nb
             call rootfind(j, m, theta, xb1(i), xb2(i), xacc, rt,flag)
             if (flag.eq.1) go to 10
             if (nb.ge.2) then
                do nm=1,i-1
                   if (rt.eq.rts(nm)) go to 10
                end do
             end if
          ll = ll+1
          rts(ll)=rt
10        continue
       end if
       nb=ll

       return
       end


       subroutine brackets(j,m,theta,x1, x2, n, xb1, xb2, nb)
c   Purpose:
c      To calculate brackets for up to nb distinct intervals which each contains at least
c      one root. For a function fx defined on the interval [x1,x2], subdivide the interval
c      n equally spaced segments, and search for zero crossings of the function. 
c  Input:
c      j: jth interface. 
c      m: complex refractive index.
c      theta: outgoing angle (scattering angle).
c      x1,x2: For a function fx defined on the interval [x1,x2]
c      n: number that the inteval [x1,x2] is equally subdivided
c      nb: as input, it is the maximum number of roots sought
c  Output:
c      nb: as output, it is the number of bracketing pairs xb1[1..nb], xb2[1..nb] that are
c          found.
c      xb1(nb),xb2(nb): bracketed pairs within which the root is to be found
c  Local:
c      nbb,i:
c      x,fp,fc,dx,fx
c  Subroutines:
c      rtfunction(j,m,theta,T,f,df): user-defined subroutine to evaluate function value and 
c                          its derivative value.
c  References:
c      Press W H , S A Teukolsky, W T Vetterling and B P Flannery,"Numerical  recipes in C: 
c      the art of scientific computing", Second edition, Cambridge University Press,1992. 
c      p.352. But note that the sign change (fc*fb<0) is used as the only condition to judge
c      if there is a root in the subdivided interval is not complete. If the one of the limits
c      of the interval is exact the root, there is no way that fc*fb<0.0 is satisfied, the
c      root is lost. Therefore, we use fc*fb<=0.0 as the condition.

       real*8 theta
       integer*4 j
       complex*16 m
       integer*4 n, nb
       integer*4 nbb, i
       real*8 x1, x2
       real*8 xb1(nb), xb2(nb)
       real*8 x, fp, fc, dx, f, df
       
       nbb = 0
       dx=(x2-x1)/dfloat(n)
       x = x1
       call rtfunction(j,m,theta,x,f,df)
       fp = f
       do i =1, n
          x = x + dx
          call rtfunction(j,m,theta,x,f,df)
          fc = f
          if ((fc*fp).le.(-1.0d-15)) then
             nbb = nbb + 1
             xb1(nbb)= x-dx
             xb2(nbb)=x
c maximum number of roots is reached
             if (nbb.eq.nb) return
          end if
          fp = fc
       end do
       nb = nbb
       return
       end
       
       subroutine rootfind(j, m, theta, x1, x2, xacc, rt, flag)
c   Purpose:
c      To find the root of a function fx(x) bracketed between x1 and x2, using a 
c      combination of Newton-Raphson and bisection methods. The root is refined till
c      its accuracy is within -xacc and +xacc. 
c  Input:
c      j: jth interface. 
c      m: complex refractive index.
c      theta: outgoing angle (scattering angle).
c      x1,x2: limits that the root is bracketed. It is from subroutine 
c             brackets(j,m,theta,x1, x2, n, xb1, xb2, nb). x1=xb1(i), x2=xb2(i).
c      xacc: accuracy that is used to find the root.
c  Output:
c      rt: the root found that lies within [x1,x2]
c      flag: =0 root is found within the maximum iiterations
c            =1 root is not found within the maximum iiterations, and thus assume there is
c               no root.
c  Local:
c      j: integer for loop
c      df: the first derivative of the function fx
c      dx: the last step size
c      dxold: the step size before last
c      f: function value
c      fl: function value at lower bound (<0)
c      fh: function value at high bound (>0)
c      temp: temporary value to store rt
c      xl: lower bound corresponding to fl (<0)
c      xh: higher bound corresponding to fh (>0)
c  Subroutines:
c      rtfunction(j,m,theta,T,f,df): user-defined subroutine to evaluate function value and 
c                          its derivative value.
c  References:
c      Press W H , S A Teukolsky, W T Vetterling and B P Flannery,"Numerical  recipes in C: 
c      the art of scientific computing", Second edition, Cambridge University Press,1992. 
c      p.366. 
       
       real*8 theta
       integer*4 j
       complex*16 m
       real*8 x1, x2, xacc, rt
       integer*4 i,maxit,flag
       real*8 df, dx, dxold, f, fl, fh, temp, xl, xh
       
       maxit=10000
       flag=0
       
       call rtfunction(j,m,theta,x1,fl,df)
       call rtfunction(j,m,theta,x2,fh,df)
       if ( fl*fh.gt.0.0d0) 
     &       call ERRMSG( 'rootfind -- root must be bracketed',.TRUE.)
       if (fl.eq.0.0d0) then 
           rt = x1
           return
       end if
       
       if (fh.eq.0.0d0) then 
           rt = x2
           return
       end if
       
       if (fl.lt.0.0d0) then
           xl = x1
           xh = x2
       else
           xl = x2
           xh = x1
       end if
       
       rt = 0.5*(x1+x2)
       dxold = dabs(x2-x1)
       dx = dxold
       call rtfunction(j,m,theta,rt,f,df)
       
       do i=1,maxit
          if ((((rt-xh)*df-f)*((rt-xl)*df-f).ge.0.0d0).or.
     &         (dabs(2.0*f).gt.dabs(dxold*df))) then
             dxold = dx
             dx = 0.5*(xh-xl)
             rt = xl+dx
             if (xl.eq.rt) return
          else
             dxold = dx
             dx = f/df
             temp = rt
             rt = rt - dx
             if (temp.eq.rt) return
          end if
          if (dabs(dx).lt.xacc) return
          call rtfunction(j,m,theta,rt,f,df)
          if (f.lt.0.0d0) then
             xl = rt
          else
             xh = rt
          end if
       end do
       
       flag=1
c       call ERRMSG( 'rootfind -- Maximum number of iterations exceeded',
c     &               .TRUE.)
       return
       end
     
       subroutine rtfunction(j,m,theta,T,f,df)
c   Purpose:
c      To find roots for equation (50). For given j, m, and theta
c  Input:
c      j: jth interface. 
c      m: complex refractive index.
c      theta: outgoing angle (scattering angle).
c      T: = tan(thetai/2), thetai is the incident angle.
c  Output:
c      f: function value.
c      df: value of the derivative of the function
c  Local:
c      re: real part of complex refractive index.
c      thetat: refractive angle.
c      thetatd: derivative of refractive angle with respect to T
c  Subroutines:
c      None
c  References:
c      None
       real*8 theta, T, f, df
       integer*4 j
       complex*16 m
       real*8 re
       real*8 thetat,thetatd

       re = dreal(m)
       thetat=dasin(2.0*T/re/(1.0+T*T))
       thetatd=dsqrt(re*re*(1.0+T*T)**2-4.0*T*T)
       f=(-1.0d0)**(j-2)*((1.-6.*T*T+T**4)*dcos(2.*dfloat(j-1)*thetat)
     &                    +4.*T*(1.-T*T)*dsin(2.*dfloat(j-1)*thetat))
     &    /(1.+T*T)**2 - dcos(theta)
       
       df = 4.*(-1.d0)**(j-2)/(1.+T*T)**3
     &      *(1.+dfloat(j-1)*(T*T-1.)/thetatd)
     &      *(4.*T*(T*T-1.)*dcos(2.*dfloat(j-1)*thetat)
     &      +(1.-6.*T*T+T**4)*dsin(2.0*dfloat(j-1)*thetat))
       return
       end
