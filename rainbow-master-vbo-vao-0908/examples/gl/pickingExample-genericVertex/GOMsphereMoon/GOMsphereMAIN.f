      program testGOM
c  This main program is an example to show how to call GOMsphere subroutine to calculate  the
c  spectral variation of single scattering albedo and asymmetry factors for both near- and far-
c  field. 
c  Input variables:
c    lamda: wavelength in meter. (For this program, lamda=0.3-4.0 micrometer) [300 ~ 4000 nm] 
c    radius: radius of the spherical scatters ( snow grains)in meters.
c    theta: scattering angle.
c    1. If only need to calculate absorption efficiency Qabs, just specify lamda and radius as inputs.
c    2. If only need to calculate the scattering efficiency Qscat and asymmetry factor g, just specify 
c       lamda and radius as inputs.
c    3. theta must be specified only for phase function (PN, PF) calculations. For other calculations,theta
c       can be given any value, which do not affect the results.
c  Output variables:
c    Qabs: absorption efficiency. It is the same for both near and far fields.
c    QNscat: scattering efficiency for near-field scattering.
c    QFscat: scattering efficiency for far-field scattering.
c    gN: asymmetry factor for near-filed scattering.
c    gF: asymmetry factor for far-field scattering.
c    PN: phase function for near-field scattering.
c    PF: phase function for far-field scattering.
c    a0: single scattering albedo=scattering efficiency/(absorption+scattering) efficiency
c    a0N: single scattering albedo for near-field
c    a0F: single scattering albedo for far-field
c  Local variable:
c    n: number of roots for Gauss-Legendre n-point quadratures.
c    pi: =acos(-1.0)=3.1415926
c
c  NOTE:
c   This program is distributed publically along a paper published in Applied Optics. Citation is suggested as:
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c
c   Any errors or comments please direct them to Dr. Xiaobing Zhou at xzhou@mtech.edu 
c   or xzhou_2000_2001@yahoo.com. Thanks.  May, 2008.
c double lambdaStart = 400e-9;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um

      real*8 lambda, radius, pi, theta, stepLambda, stepRadius, 
     &       stepTheta, lambdaNano
      real*8 Qabs, QNscat,QFscat,gN,gF, PN,PF, a0, Fscat, phase
       
      real*8  lambdaStart, lambdaEnd, radiusStart, radiusEnd, 
     &         thetaStart, thetaEnd
      
      integer*4 nSpectralSampleSteps, nRadiusteps, nThetaSteps
      
      lambdaStart = 400d-9
      lambdaEnd = 700d-9
    
      nRadiusSteps = 2
      
	  nSpectralSampleSteps = 120    
      nThetaSteps = 100
      
c      radiusStart =  0.05e-3
      radiusStart =  1.0d-3
      radiusEnd = 2.0d-3

      pi = dacos( -1.0d0)

	  thetaStart = 130.0 * pi / 180.0
	  thetaEnd = 142.0  * pi / 180.0

      stepLambda = ( lambdaEnd - lambdaStart ) 
     &                              / dfloat(nSpectralSampleSteps)
            

	  stepRadius = (radiusEnd - radiusStart) /  dfloat( nRadiusSteps )
      stepTheta = (thetaEnd - thetaStart ) / dfloat(nThetaSteps)
      
      
      open (10,file='scattering2_120_fort.txt', status='unknown')
      open (12,file='phaseFunction2_120_100_fort.txt', status='unknown')
      write (10,'(a7,1x,a10,1x,a13)' )
     &      '#radius','wavelength','QFscat * area'
      write (10,'(I6,1x,I6)' )  nRadiusSteps+1, nSpectralSampleSteps+1
      
      write (12,'(a6,1x,a10,1x,a5, 1x, a5)' ) '# radius', 'wavelength', 
     &       'theta','phase'
      write (12,'(I6,1x,I6,1x,I6)' ) nRadiusSteps+1, 
     &                       nSpectralSampleSteps+1, nThetaSteps+1
       
       
c  Specify a scattering angle just as an input, any value of theta will be OK of we only
c  concern the quantities except the phase functions.
       radius = radiusStart
              
       theta= thetaStart
c      theta could be any value when computing scattering cross section 
c  Task to calculate scattering cross section
c      write(*,*): default output
c     write (*,*) 'Working on Task 1...'
c     
      do while ( radius .le. radiusEnd )
       lambda = lambdaStart   
       do while (lambda .le. lambdaEnd )
           
         call GOMsphere(lambda,radius,theta,
     &                     Qabs,QNscat,QFscat,gN,gF,PN,PF)
     
         Fscat =QFscat * ( pi * radius * radius ) 
         
         lambdaNano = lambda * 1.0e9
         
         write (10,'(e10.4, 1x, e14.5, 1x, e14.5)' ) 
     &                radius, lambdaNano, Fscat
     
c        write(*,'(e10.4, 1x, e14.5, 1x, e14.5)'  ) 
c    &              radius, lambdaNano, Fscat
         
         lambda = lambda + stepLambda
       end do
       
       radius = radius + stepRadius
      end do
            

      close(10)
       
c Task  to calculate far-field scattering efficiency QFscat versus size parameter x for lamda=2.0micrometer
      radius = radiusStart
            
      do while ( radius .le. radiusEnd )
          
       lambda = lambdaStart   
       
       do while (lambda .le. lambdaEnd )
           
          theta = thetaStart 
          

          do while (theta .le. thetaEnd)
              
           call GOMsphere(lambda,radius,theta,
     &                     Qabs,QNscat,QFscat,gN,gF,PN,PF)
              
           lambdaNano = lambda * 1.0d9
           
           write (12,'(e10.4, 1x, e14.5, 1x, e14.5, 1x, e14.5)' ) 
     &               radius, lambdaNano, theta * (180.0/pi),
     &               PF / ( 4.0 * pi)
c          we need the phase function whose integration over sphere is 1, so divide by
c          4*pi

c           write (*,'(e10.4, 1x, e14.5, 1x, e14.5, 1x, e14.5)' ) 
c     &            radius, lambdaNano, theta * (180.0/pi), PF
           
           theta = theta + stepTheta
           
          end do
         
          lambda = lambda + stepLambda
       end do
       
       radius = radius + stepRadius
      end do
      
      close(12)
      
      stop
      end

       subroutine GOMsphere(lamda,radius,theta,
     &                     Qabs,QNscat,QFscat,gN,gF,PN,PF)
c  Purpose:
c    To calculate absorption efficiency, scattering efficiency, asymmetry factor and phase 
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
c    QNscat: scattering efficiency for near-field scattering.
c    QFscat: scattering efficiency for far-field scattering.
c    gN: asymmetry factor for near-filed scattering.
c    gF: asymmetry factor for far-field scattering.
c    PN: phase function for near-field scattering.
c    PF: phase function for far-field scattering.
c  Local variable:
c    n: number of roots for Gauss-Legendre n-point quadratures.
c  NOTE:
c   This program is to be distributed publically along a paper to be submitted to Applied Optics.
c   Any errors or comments please direct them to Xiaobing Zhou at xzhou@gi.alaska.edu or xzhou@nmt.edu
c   or xzhou_2000_2001@yahoo.com. Thanks.  March, 2002; May 2003.
       real*8 lamda, radius
       real*8 Qabs,QNscat,QFscat,gN,gF,theta,PN,PF
       integer*4 n
       
       n=20

c  To calculate the absorption efficiency Qabs.
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
c      absFF: value of the integrand function for absorption efficiency calculation
c  LOCAL VARIABLES:
c      j: loop varible
c      beta: absorption coefficient
c      xi: path length between two reflection events inside the sphere
c      temp: snow temperature in K, only as an input to ,no effect in VNIR
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
      real*8 u,v
      
c      complex*16 Rep,Rev,tp,tv,rp,rv,t,r,u,v
c      u and v are defined as reaL*8 within TransReflec function, but defined to be complex*16 here

      temp=293.0d0
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
c     Gaussian quadratures. Eqs.(43)-(44) and (47-49).
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
c     asymNF1,asymNF2: value of integrand function in g*Qscat for near field (see Eq.(44b))
c     asymFF1,asymFF2: value of integrand function in g*Qscat for far field (see Eq.(44b))
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

       temp=293.0d0
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
          Fscat=Fscat+w(i)*(scatfFB1 +scatfFB2)
          gFQscat=gFQscat+w(i)*(asymFFB1 +asymFFB2)
       end do

       Nscat=Nscat*xk
       gNQscat=gNQscat*xk
       gN=gNQscat/Nscat

       Fscat=Fscat*xkp
       gFQscat=gFQscat*xkp
c   Up to now, gF is actually the asymmetry factor due only to diffraction       
       gF=gFQscat/Fscat
c   To obtain the asymmetry factor for far-field scattering (Eq.(49))
       gF = ((QNscat * gN) + gF) / (1. + QNscat)
       
       return
       end

      subroutine scat_efficiency(lamda,radius,m,QNscat,QFscat)
c PURPOSE:
c     To calculate scattering efficiency for both near filed (QNscat) and far filed
c     (QFscat), given wavelength, radius of sphere, number of roots for Gaussian 
c     quadratures. eq.(44) & eq.(48)
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
c               calculation. See eq.(43)
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
      temp=293.0d0

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
c      For Far-field, Qsca^far=Qsca^near+1                 Eq.(5.47)
c      
c  REFERENCE: 
c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large 
c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003. 
c  INPUT VARIABLES:
c      theta: varible, scattering angle in radians
c      m: complex, refractive index
c      beta: absorption coefficients
c      xi: real, path-length between two interfaces
c      NN: the number of interfaces when truncation is carried out with truncation 
c         error = 10^-10
c  OUTPUT VARIABLES:
c      scatfN: value of the integrand function for near-field scattering efficiency calculation
c      asymNF: value of the integrand function for asymmetry factor calculation.
c              asymNF=g*Qscat for near field.
c      asymFF: value of the integrand function for asymmetry factor calculation.
c              asymFF=g*Qscat for far field.

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
c      Function BessJ1 should be declared to be real*8 because it is defined to be so
c      Otherwise, BessJ1  will be declared by default to be real*4, which does not match real*8
 
       real*8  BessJ1
       
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
c  To calculate the geometrical divergence factor Eq. (5.24)
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

c integrand function for near-field scattering
         asymNF=scatfN*dcos(theta)
c refraction part of integrand function of far-field asymmetry facotr (eq.50)
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
       TEMP=293.0
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
      PARAMETER ( NWL = 57, NWLT = 62)
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

       WAVLEN = WAVMET*1.0d+9
c     WAVLEN: nanometer; refractive index table is sorted by nanometer 
 
       IF( PASS1 )  THEN

         PASS1 = .False.
c                         ** Superficially test if main table messed up

c         IF( NWL.LT.100 ) CALL ERRMSG('refice--NWL value bad',.True.)
c         IF( WL(1).GT. 0.045 ) CALL ERRMSG('refice--WL(1) bad',.True.)
c         IF( WL(NWL).LT. 166 ) CALL ERRMSG('refice--WL(NWL) bad',.True.)

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

c      write(*,*) "lambda(m)=", WAVMET, "lambda(nm)=", WAVLEN
       
      IF( WAVLEN.LT.WL(1) .OR. WAVLEN.GT. WLT(NWLT)*1.0d+3 ) THEN

         CALL ERRMSG('refice--wavelength outside table boundaries',
     &                .False.)
         ref_ice = (0.,0.)
         RETURN

      END IF


      IF( WAVLEN .LE. 167*1.0d+3  ) THEN
c                                  ** Wavelength between 0.045*10^3 nm and 167*10^3 nm
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

c        ** Wavelength greater than 167 microns => temp dependence
c        ** (temperature-dependent case)

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

c     get the refractive index
      MIM = 0.0d0;
      ref_ice = CMPLX( MRE, MIM)

      END

      BLOCK DATA ICECON

c        Ice-refractive-index vs. wavelength table for subroutine refice

      IMPLICIT NONE

c     .. Parameters ..

      INTEGER   NWL, NWLT
      PARAMETER ( NWL = 57, NWLT = 62 )
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
     &  181.78736, 1.468492, 0.0,
     & 182.61377, 1.465508, 0.0,
     & 183.65075, 1.461964, 0.0,
     & 184.57517, 1.458976, 0.0,
     & 184.96831, 1.457748, 0.0,
     & 185.34523, 1.456594, 0.0,
     & 186.71302, 1.452589, 0.0,
     & 188.30587, 1.448250, 0.0,
     & 188.95226, 1.446576, 0.0,
     & 189.90445, 1.444204, 0.0 /
     
      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 11, 20 ) /
     & 190.74940, 1.442176, 0.0,
     & 191.60818, 1.440187, 0.0,
     & 192.92449, 1.437273, 0.0,
     & 193.43690, 1.436177, 0.0,
     & 193.74245, 1.435535, 0.0,
     & 194.44617, 1.434084, 0.0,
     & 195.47439, 1.432034, 0.0,
     & 196.31429, 1.430413, 0.0,
     & 197.15374, 1.428844, 0.0,
     & 197.97647, 1.427348, 0.0  /

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 21, 30 ) /
     & 198.91056, 1.425704, 0.0,
     & 200.84066, 1.422463, 0.0,
     & 204.22306, 1.417245, 0.0,
     & 207.60622, 1.412533, 0.0,
     & 210.40121, 1.408969, 0.0,
     & 213.13806, 1.405734, 0.0,
     & 217.53669, 1.400987, 0.0,
     & 222.57015, 1.396146, 0.0,
     & 226.34187, 1.392869, 0.0,
     & 231.16733, 1.389045, 0.0/ 

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 31, 40 ) /
     & 239.02641, 1.383575, 0.0,
     & 244.08006, 1.380471, 0.0,
     & 248.79191, 1.377822, 0.0,
     & 253.99695, 1.375140, 0.0,
     & 262.88107, 1.371063, 0.0,
     & 273.47659, 1.366882, 0.0,
     & 283.11246, 1.363603, 0.0,
     & 289.44400, 1.361673, 0.0,
     & 296.814, 1.359619, 0.0,
     & 312.657, 1.355795, 0.0/
     
      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 41, 50 ) /
     & 334.244, 1.351588, 0.0,
     & 365.119, 1.346994, 0.0,
     & 404.770, 1.342710, 0.0,
     & 435.957, 1.340179, 0.0,
     & 480.126, 1.337418, 0.0,
     & 486.269, 1.337091, 0.0,
     & 546.227, 1.334435, 0.0,
     & 587.725, 1.333012, 0.0,
     & 644.025, 1.331436, 0.0,
     & 656.454, 1.331126, 0.0/

      DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 51, 57 ) /
     & 706.714, 1.329996, 0.0,
     & 780.237, 1.328583, 0.0,
     & 852.344, 1.327365, 0.0,
     & 894.596, 1.326702, 0.0,
     & 1014.26, 1.324921, 0.0,
     & 1083.33, 1.323902, 0.0,
     & 1128.95, 1.323216, 0.0/



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
