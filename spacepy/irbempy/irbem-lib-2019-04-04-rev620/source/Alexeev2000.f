!***************************************************************************************************
! Copyright 2000, 2002, Igor Alexeev
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
!
***************************************************************************
*                                					   
*									   
*                     DYNAMIC PARABOLOID MODEL OF			   
*                      THE EARTH'S MAGNETOSPHERE			   
*                           (March 2000)				   
*	              (Revision on 19 June 2002)			   
*									   
* The dynamic paraboloid model  allows to calculate the 		   
* magnetic field inside the Earth's magnetosphere. In the base of the      
* model lies CPMOD (closed paraboloid model) software. Dynamic paraboloid  
* model enable calculation the time variations of each large scale source  
* of magnetospheric magnetic field by empirical data. The possibilities    
* are provided to "switch on" and "switch off" separate sources of the     
* magnetic field, and to change the parameters defining their intensity.   
*									   
* Contacts: Igor Alexeev						   
*           Scobeltsyn Inst. of Nuclear Physics 			   
*           Moscow State University,            			   
*           Moscow, 119899, Russia					   
*           alexeev@dec1.sinp.msu.ru					   
***************************************************************************

      subroutine a2000(ro,v,bimf,dst,al,x,bm,ifail)
*--------------------------------------------------------------------------
*  Calculation of magnetic field vector Bm in point x(3) 		   
*              inside the Earth's magnetosphere 			   
*  INPUT: 								   
*           time moment UT (Universal Time) by				   
*           year=iy;							   
*           month=mo;							   
*           day in month=id.						   
*           ro, V, bz - solar wind density and velocity, IMF Bz.	   
*           dst - value of Dst index;					   
*           AL  - value 0f al index. 					   
*  OUTPUT: magnetic field components at the point x0(3) 		   
*                      in GSM coordinates, nT:				   
*         bm(i) - total magnetic field (i=1,3)(here only external field); 			   
*         bb(1,i) - geomagnetic dipole magnetic field;			   
*         bb(2,i) - ring current magnetic field;			   
*         bb(3,i) - geomagnetic tail currents magnetic field;		   
*         bb(4,i) - magnetic field of CF currents shielding dipole;	   
*         bb(5,i) - magnetic field of CF currents shielding ring current;  
*         bb(6,i) - magnetic field of Region 1 FAC;     		   
*         bb(7,i) - IMF penetrated into the magnetosphere.      	   
*									   
* WARNING: Because of the paraboloid coordinates singularity, avoid	   
*          the magnetic field calculations at the Ox axis.		   
* 									   
*          Written by I. Alexeev					   
*--------------------------------------------------------------------------
      DIMENSION x(3), bm(3), par(10), bb(7,3), bimf(3)
      REAL*8   a2000_ut
      INTEGER*4 a2000_iyear,a2000_imonth,a2000_iday,ifail
*
      COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday
*
      ut=a2000_ut
      iy=a2000_iyear
      mo=a2000_imonth
      id=a2000_iday
      ifail=0
      call submod(ut,iy,mo,id,ro,v,bimf,dst,al,par)
      r1=par(6)
       IF(x(1)+0.5/R1*(x(2)**2+x(3)**2).GT.R1) then
       do i=1,3
       bm(i)=0.
       do j=1,7
       bb(j,i)=0.
       end do
       end do
       ifail=-1
       return
       end if
      call a_field(x,par,bm,bb)
      do 4 i=1,3
4     bm(i)=bm(i)-bb(1,i)
      return
      END
C
      subroutine submod(ut,iy,mo,id,ro,v,bimf,dst,al,par)
*---------------------------------------------------------------------*
*  Calculation of the paraboloid model input parameters 	      *
*  by empirical data						      *
*  INPUT (empirical data):					      *
*           time moment UT (Universal Time)			      *
*           year=iy;						      *
*           month=mo;						      *
*           day in month=id.					      *
*           ro, V, bimf- solar wind density and velocity, IMF.	      *
*           dst - value of Dst index;				      *
*           AL  - value of al index. 				      *
*  OUTPUT (model input parameters): 				      *
*   	      par(1) - geomagnetic dipole tilt angle, degrees;        *
*   	      par(2) - dipole magnetic field at equator, nT;	      *
*   	      par(3) - magnetic flux through the tail lobes, Wb;      *
*   	      par(4) - maximum ring current intensity, nT;	      *
*   	      par(5) - the total current of Region 1 FAC, MA;	      *
*   	      par(6) - magnetopause stand-off distance, Re;	      *
*   	      par(7) - distance to the inner edge of geotail	      *
*   		       current sheet;				      *
*   	      par(8-10) -IMF penetrated components in GSM coord., nT. *
*								      *
*    Written by V. Kalegaev					      *
*---------------------------------------------------------------------*
      DIMENSION x(3), bm(3), bimf(3), par(10), bb(7,3)
      iday=IDD(iy,mo,id)                ! the day number in year
*
***>  calculation of the tilt angle, tpsi [degrees],
***>  and the dipole magnetic field at the Earth's equator, BD [nT] 
*
      call trans(ut,iday,iy,tpsi,bd)

*
***>  calculation of the magnetopause stand-off distance, R1 [Re], 
*     by Shue et al. 1997
*
      Psw=(1.67e-6)*ro*V*V              ! Solar wind pressure [nT]          
      bz=bimf(3)
      if (bz.ge.0.) then
      r1=(11.4+0.013*bz)*(psw**(-1./6.6))
      else
      r1=(11.4+0.14*bz)*(psw**(-1./6.6))
      end if
*
***>  calculation of the distance to the geotail inner edge: R2 [Re]
*
      if(dst.ge.-10.0)then
      R2=0.7*r1
      goto 3
      end if
      fie=74.9-8.6*alog10(-dst)
      fi=fie*3.1416/180.0
      R2=1.0/cos(fi)/cos(fi)
 3     al0=sqrt(1.+2*R2/R1)
*
***>  calculation of the geotail lobe magnetic flux: FLUX [Wb] 
*
      BT=-AL/7                                                                  
      FLUX=1.5708*R1*R1*AL0*BT*6.37816*6.37816*1.E3                            
      flux=flux+395.98278e6                                                      
*
***>  calculation of the ring current mag. f. in the Earth's
*     center: BR [nT] 
*
      BR=dst-10.
      if(dst.ge.0.0)br=-10.0

*
***>  calculation of the total FAC intensity: AJ0 [MA] 
*
      AJ0=0.327744
      if(bz.le.-1.61133) AJ0=-1.017*bz/5.
      AJ0=2.*AJ0*sqrt(400./v)*((5/ro)**0.125)

       par(1)=tpsi
       par(2)=BD
       par(3)=FLUX
       par(4)=BR
       par(5)=AJ0
       par(6)=R1
       par(7)=R2
      do 2 i=1,3
2     par(i+7)=bimf(i)
      return
      END
C
      SUBROUTINE TRANS (UT,IDAY,iyear,tpsi,BD)
***************************************************************
*   Calculation of the geomagnetic dipole tilt angle and      .
*   matrices of transition from the geographic coordinates    .
*   (Cartesian) to the GSM-coordinates (G2GSM(3,3) in the     .
*   COMMON BLOCK /TRAN/).                                     .
*                _            _                               .
*                Xgsm = G2GSM*Xgeogr                          .
*                                                             .
*    ALPHA1 is the angle between the Earth's axis and the     .
*           dipole moment;                                    .
*    ALPHA2 is the angle between the Earth's axis and         .
*           the normal to the ecliptic plane;                 .
*    PHIM   is the angle between the midnight meridian plane  .
*           and the meridional plane;                         .
*    PHISE  is the angle between the Earth-Sun line and       .
*           the projection of the Earth's axis on the         .
*           ecliptic plane;                                   .
*    PSI    is the tilt angle of the geomagnetic dipole;      .
*    B      is the Sun's deflection;                          .
*    UT     is universal time;                                .
*    IDAY   is the day number;                                .
*    B1     is western longitude of the noon meridian;        .
*    B2     is the angle between the noon geogr. meridian     .
*           and the Y=0 plane in the GSM-coordinates.         .
*       Written by I. Alexeev                                 .
*							      .
*   Acknowledgment:                                           .  
*       Program SUN written by Gilbert D. Mead is used	      .
*       to determine some parameters dependent on the Earth   .
*       and Sun mutual position (SLONG, sind, and cosd from   .
*       COMMON BLOCK /ddd/).                                  .
***************************************************************
      COMMON/TRAN/g2gsm(3,3)
      common /ddd/ sind,cosd
      dimension gauss(3)         
      Ihour=int(ut)
      dmin=(ut-ihour)*60.
      min=int(dmin)
      isec=int((dmin-min)*60.)
      call SUN_a2000(IYeaR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
      P= 3.1415926/180.
      pi2=3.1415926/2
          call dipgarm(iyear, gauss)
          G1=GAUSS(2)                                                           
          H1=GAUSS(3)                                                           
          PD=G1**2+H1**2                                                        
          G10=GAUSS(1)                                                           
          BD=-SQRT(G10*G10+PD)                                                    
          ALPHA1= ATAN(-1.*SQRT(PD)/G10)                                         
          PHINP=  ATAN(H1/G1)                                                   
        ALPHA2= 23.4419*P
        PHIM= P*(UT*15+phinp/p)
        PHISE= pi2-slong+9.924E-5
         
        B1=P*15.*(UT-12)
        SB= SIN(ALPHA2)*COS(PHISE)
        CB= SQRT(1-SB*SB)
        SB1=SIN(B1)
        CB1=COS(B1)
          SPSI= -SB*COS(ALPHA1) + CB*SIN(ALPHA1)*COS(PHIM)
          CPSI=  SQRT(1-SPSI*SPSI)
          psi=asin(spsi)/p
        CB2=(COS(ALPHA1)+SB*SPSI)/CPSI/CB
        SB2=SQRT(abs(1-CB2*CB2))
        IF(PHIM.LE.0..OR.PHIM.GE.3.1415926) SB2=-SB2
          tpsi=psi
      g2gsm(1,1)=cb1*cb
      g2gsm(1,2)=-sb1*cb
      g2gsm(1,3)=sb
      g2gsm(2,1)=sb1*cb2-cb1*sb*sb2
      g2gsm(2,2)=cb1*cb2+sb1*sb*sb2
      g2gsm(2,3)=cb*sb2
      g2gsm(3,1)=-sb1*sb2-cb1*sb*cb2
      g2gsm(3,2)=-cb1*sb2+sb1*sb*cb2
      g2gsm(3,3)=cb*cb2
         RETURN
         END
  
      FUNCTION IDD(iy,mo,id)
********************************************************************
*  Calculation of the day number in a year			   *
*  INPUT PARAMETERS: year (IY), month (MO), day in the month (ID)  *
*  Written by V. Kalegaev					   *
********************************************************************
      II=0
      IF (MO.EQ.1) GOTO 4
 5    DO 10 M=1,MO-1
      GOTO (1,2,1,3,1,3,1,1,3,1,3) M
 1    II=II+31
      GOTO 10
 2    II=II+28
      GOTO 10
 3    II=II+30
 10   CONTINUE
 4    II=II+ID
      if (mod(iy,100).eq.0.and.mod(iy,400).ne.0) goto 6
      IF (MOD(IY,4).EQ.0.AND.II.GT.59.AND.MO.GT.2) II=II+1
 6    IDD=II
      RETURN
      END

C     
      subroutine dipgarm(iyear,gauss)
*------------------------------------------------------------------------*
*   Calculation of the first three Gaussian coefficients for given year  *
*   year>=1900  							 *
*   Last IGRF values are for 1995.					 *
*   Written by V. Kalegaev						 *
*------------------------------------------------------------------------*
      dimension gauss(3), g(60), sg95(3)
      data g/
     *-31543., -2298.,  5922.,
     *-31464., -2298.,  5909.,
     *-31354., -2297.,  5898.,
     *-31212., -2306.,  5875.,
     *-31060., -2317.,  5845.,
     *-30926., -2318.,  5817.,
     *-30805., -2316.,  5808.,
     *-30715., -2306.,  5812.,
     *-30654., -2292.,  5821.,
     *-30594., -2285.,  5810.,
     *-30554., -2250.,  5815.,
     *-30500., -2215.,  5820.,
     *-30421., -2169.,  5791.,
     *-30334., -2119.,  5776.,
     *-30220., -2068.,  5737.,
     *-30100., -2013.,  5675.,
     *-29992., -1956.,  5604.,
     *-29873., -1905.,  5500.,
     *-29775., -1848.,  5406.,
     *-29682., -1789.,  5318./
      data sg95/17.6, 13.0, -18.3/
      i1=(iyear-1900)
      i2=i1/5
      i3=i1-5*i2
      n=(i2)*3+1
      if(i2.gt.18) goto 1
      do i=1,3
      a=g(n-1+i)
      b=g(n+2+i)
      gauss(i)=a+(b-a)/5.*i3
      end do
      return
 1    continue
      do i=1,3
      a=g(57+i)
      gauss(i)=a+sg95(i)*(iyear-1995)
      end do
      return
      end

      subroutine A_field(x0,par,bm,bdd)
******************************************************************************
*  Calculation of the magnetic field in the magnetosphere                    *
*  inthe point x0(3) by model input parameters par(10).	                     *
*      Ver. 3 on March 2000.                                                 *
*  (The model's parameters par(1-10) are moved to  COMMON/T2/ )              *
*									     *
*  INPUT PARAMETERS: x0(3) is a point where the magnetic field is being      *
*                          calculated, in GSM coordinates, Re;		     *
*                    par(1) - geomagnetic dipole tilt angle, degrees;	     *
*                    par(2) - dipole magnetic field at equator, nT;	     *
*                    par(3) - magnetic flux through the tail lobes, Wb;      *
*                    par(4) - maximum ring current intensity, nT;	     *
*                    par(5) - the total current of Region 1 FAC, MA;	     *
*                    par(6) - magnetopause stand-off distance, Re;	     *
*                    par(7) - distance to the inner edge of geotail 	     *
*                             current sheet;				     *
*                    par(8-10) -IMF penetrated components in GSM coord., nT. *
*                                                  			     *
*  OUTPUT - magnetic field components at the point x0(3) in GSM  coord., nT. *
*         bm(i) - total magnetic field (i=1,3); 			     *
*         bdd(1,i) - geomagnetic dipole magnetic field; 		     *
*         bdd(2,i) - ring current magnetic field;			     *
*         bdd(3,i) - geomagnetic tail currents magnetic field;		     *
*         bdd(4,i) - magnetic field of CF currents shielding dipole;	     *
*         bdd(5,i) - magnetic field of CF currents shielding ring current;   *
*         bdd(6,i) - magnetic field of Region 1 FAC;     		     *
*         bdd(7,i) - IMF penetrated into the magnetosphere.     	     *
*									     *
* WARNING: Because of the paraboloid coordinates singularity, avoid	     *
*          the magnetic field calculations at the Ox axis.		     *
*									     *
* Written by  V.Kalegaev	         				     *
******************************************************************************
      COMMON/TK/B1(3),B2(3),B1A(3),B1R(3),B2A(3),B2R(3),
     *BA(3),DB(3),AB,B(3)
      COMMON/BEGFC/B1CF(3),B1CFD(3),B1CFR(3),B1D(3),B1RC(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,bd,bd0
      COMMON/IMFd/Bimf(3),b3(3)
      common/fac12/ami1,ami2,tm1,tm2
      common/bfac12/bfac1(3),bfac2(3)
      DIMENSION X0(3), FF(3), bm(3), par(10), bdd(7,3)

      psi=par(1)
      bd=par(2)
      flux=par(3) ! magnetic flux through the tail lobes
      BR=par(4)   ! maximum ring current intensity
      ami1=par(5) ! maximum of Region 1 FAC intensity
      R1=par(6)   ! magnetopause stand-off distance
      R2=par(7)   ! distance to the inner edge of geotail current sheet
      do 2 i=1,3
2     bimf(i)=par(i+7) ! IMF components
        call mas1d (flux,BR,R1,R2)
        call field(X0,FF)

         do 1 i=1,3
         bm(i)=b(i)         ! total magnetic field
         bdd(1,i)=b1d(i)    ! geomagnetic dipole contribution
         bdd(2,i)=b1rc(i)   ! ring current contribution
         bdd(3,i)=b2(i)     ! geomagnetic tail contribution
         bdd(4,i)=b1cfd(i)  ! contribution of CF currents shielding dip.
         bdd(5,i)=b1cfr(i)  ! contribution of CF currents shielding RC
         bdd(6,i)=bfac1(i)  !  contribution of Region 1 FAC. 
 1       bdd(7,i)=b3(i)     !  contribution of IMF. 
      return
      END


**********************************************************************
*                             C P M O D 			     *
*								     *
*                     CLOSED PARABOLOID MODEL OF		     *
*             MAGNETIC FIELD IN THE EARTH'S MAGNETOSPHERE	     *
*                     (ver. 3, February 2000)			     *
*								     *
* The CPMOD (closed paraboloid model) software was written based     *
* on the I.I.Alexeev paraboloid model A78. The paraboloid	     *
* model allows to calculate the magnetic field inside the Earth's    *
* magnetosphere. The magnetospheric magnetic field is described by   *
* sum of magnetic field of the following sources:                    *
* (1) the geomagnetic dipole, 	                                     *
* (2) the ring current;                                              *   
* (3) the current system of magnetotail including the dawn-dusk      *
* currents across the magnetotail current sheet and the closure	     *
* currents on the magnetopause;	                                     *
* magnetopause;                                                      *   
* (4) the Chapmen-Ferraro currents on the  magnetopause (the dipole  *
* magnetopause (the dipole screening field);                         *
* (5) the currents on the magnetopause screening the ring current;   *
* (6) Region 1 FAC;                                                  *
* (7) the IMF penetrated from magnetosheath.  			     *
*								     *
*  Written by I. Alexeev.					     *
**********************************************************************

      SUBROUTINE FIELD(UF,FF)                                           
***************************************************************************
* Program calculating the magnetic field in the magnetosphere.		  *
*									  *
* INPUT PARAMETERS:  x0(3) is a point where the magnetic field is being   *
*                          calculated, in GSM coordinates, Re;		  *
* OUTPUT PARAMETERS: ff(3) is a normalized vector of the magnetic	  *
*                          field at the point x0(3).			  *
* NOTE: The total magnetic field and the magnetic fields of the 	  *
*       magnetospheric current systems are stored in COMMON BLOCKs	  *
*       /TK/ and /BEGFC/ .						  *
* WARNING: Because of the paraboloid coordinates singularity, avoid	  *
*          the magnetic field calculations at the Ox axis.		  *
*									  *
* Written by I. Alexeev. 						  *
***************************************************************************
      COMMON/TK/B1(3),B2(3),B1A(3),B1R(3),B2A(3),B2R(3),                
     *BA(3),DB(3),AB,B(3)
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/COR1/AL,BE,SQ,PQ,QA                                        
      COMMON/COR2/CFI,SFI                                               
      COMMON/COR3/R,CT,ST                                               
      COMMON/GN/V2(3)                                                   
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/IMFd/Bimf(3),b3(3)
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      COMMON/BEGF/UFCF(3)
      COMMON/BEGFC/B1CF(3),B1CFD(3),B1CFR(3),B1D(3),B1RC(3)
      common/fac12/ami1,ami2,tm1,tm2
      common/bfac12/bfac1(3),bfac2(3)
      DIMENSION UF(3),FF(3),V1(3),V3(3),                                
     *UZ(3,3),B1IJ(3,3),B2AA(3,3),ZU(3,3)                               
     *,EZ(3,3),EA(3,3),V4(3),V5(3)                                      
     *,V6(3),V7(3),V8(3),V9(3)                                          
     *,A2X(3,3),X2A(3,3)                                                
      X=UF(1)                                                           
      YX=UF(2)                                                          
      Z=UF(3)                                                           
      T=SQRT(YX*YX+Z*Z)                                                 
      R=SQRT(X*X+T*T)                                                   
      CT=X/R                                                            
      ST=T/R                                                            
      Z=Z+Z0                                                            
      K=1                                                               
2     K=K+1                                                             
      R=SQRT(X*X+YX*YX+Z*Z)                                             
      RX=R/R1                                                           
      X1=X/R1-0.5                                                       
      RY=RX**2-X1-0.25                                                  
      RY=SQRT(ABS(RY))                                                  
      BE=RY+X1                                                          
      BE=SQRT(ABS(BE))                                                  
      AL=RY-X1                                                          
      AL=SQRT(ABS(AL))                                                  
      PQ=AL*BE                                                          
c      T=R1*PQ                                                           
      T=SQRT(YX*YX+Z*Z)                                                 
      QA=AL*AL+BE*BE                                                    
      SQ=SQRT(QA)                                                       
      CFI=Z/T 
      SFI=YX/T                                                          
      
                                                                
      IF(K.EQ.3)GOTO3                                                   
      CALL DERY4D(B2A,B2AA)                                             
c      CALL COM(B2A,EA)                                                  
      CALL PRIS(A2X,X2A)                                                
      Z=Z-Z0                                                            
      IF(K.EQ.2) GOTO2                                                  
3     IF(AL1-AL)4,4,5                                                   
4     CALL FLYD(V1,B1IJ)                                                
      CALL PRIS(ZU,UZ)                                                  
      GOTO 6                                                            
5     CALL BEG(V1,B1IJ)                                                
12    EZ(1,1)=0.                                                        
      EZ(1,2)=-ST*V1(1)/R-CT*V1(2)/R                                    
      EZ(1,3)=0.                                                        
      EZ(2,1)=0.                                                        
      EZ(2,2)=SFI/R*(CT*V1(1)-ST*V1(2))                                 
      EZ(2,3)=(CFI*V1(1)+(CT*CFI*V1(2)-SFI*V1(3))/ST)/R                 
      EZ(3,1)=0.                                                        
      EZ(3,2)=CFI/R*(CT*V1(1)-ST*V1(2))                                 
      EZ(3,3)=-(SFI*V1(1)+(CT*SFI*V1(2)+CFI*V1(3))/ST)/R                
      ZU(1,1)=CT                                                        
      ZU(1,2)=-ST                                                       
      ZU(1,3)=0.                                                        
      ZU(2,1)=ST*SFI                                                    
      ZU(2,2)=CT*SFI                                                    
      ZU(2,3)=CFI                                                       
      ZU(3,1)=ST*CFI                                                    
      ZU(3,2)=CT*CFI                                                    
      ZU(3,3)=-SFI                                                      
      DO 7 I=1,3                                                        
      DO 7 J=1,3                                                        
      UZ(I,J)=ZU(J,I)                                                   
7     CONTINUE                                                          
6     CALL PERE(V1,B1,ZU) 
        call bdipc(uf,b0,b1d)
        if (al.ge.al1) then
        do i=1,3
        b1cf(i)=b1(i)-b1d(i)
        end do
        else      
         do i=1,3
        b1cf(i)=b1(i)
        end do
        end if
	call bring(ff)
        CALL PERE(ff,B1rc,ZU)
        do i=1,3
        b1d(i)=ssd*b1d(i)*bd0
        b1cfd(i)=smr*b1cf(i)*bd0
        b1cfr(i)=b1cf(i)-b1cfd(i)
        end do

         if (bk1.eq.0.) then
        do i=1,3
         b1rc(i)=0.
         b1cfr(i)=0.
        end do
         end if
	
        do i=1,3
        b1cf(i)=b1cfr(i)+b1cfd(i)
        b1(i)=b1cf(i)+b1d(i)
        end do

      CALL PERE(B2A,B2,A2X)                                             
      CONTINUE                                                          
      CALL PERE(B1,V2,X2A)                                              
      DO 8 I=1,3                                                        
      BA(I)=B2A(I)+V2(I)                                                
8     CONTINUE                                                          
      CALL PERE(BA,B,A2X)                                               

      b3(1)=0.0116*Bimf(1)*simf
      b3(2)=0.0978*Bimf(2)*simf
      b3(3)=0.0978*Bimf(3)*simf

      call FAC(uf,Bfac1,bfac2)    ! BFAC2=0 in this version

      do 11 j=1,3
c11    b(j)=b(j)+B3(j)+bfac1(j)+bfac2(j)
11    b(j)=b1(j)+b2(j)+B3(j)+bfac1(j)+bfac2(j)+b1rc(j)
       AB=SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))                            

      DO 9 I=1,3                                                        
      FF(I)=B(I)/AB                                                     
9     CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
      SUBROUTINE mas1d (flux,Br,Rs,R2)
*****************************************************************************
* Program calculating the internal coefficients for the Bessel series,	    *
* describing the magnetic field of the magnetospheric current systems in    *
* the paraboloid model. 						    *
*									    *
* INPUT PARAMETERS: 							    *
*  FLUX is the magnetic flux through the tail lobes, Wb;		    *
*  BR	is the maximum intensity of ring current, nT;			    *
*  R1	is the distance to the subsolar point of the magnetosphere, Re;     *
*  R2	is the distance to the earthward edge of geotail current sheet, Re. *
* WARNING: You should call MAS1D after each changes in the model input      *
*          parameters.							    *
* Written by I. Alexeev. 						    *
*****************************************************************************
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      COMMON/T1/A1(12)
      COMMON/A/Y(3),F(5),V(3),U(3),YF(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      common/fac12/ami1,ami2,tm1,tm2
      dimension d1(12)
      DATA D1/0.64972264,0.21646207,
     +0.043429128,-0.000846358,-0.004917225,-0.002224403,0.94028094,
     +0.4649891,0.12928167,-0.014765534,-0.016942754,-0.022559739/
      re=6.37816
      r1=rs
      ZN = 1.
      BM = -2.*bd
      AL0=SQRT(1.+2*R2/r1)
      bt= FLUX/(1.5708*R1*R1*AL0*re*re*1.E3)          
      HN=0.06875*SQRT(HN)
      tpsi=psi
      PSI=PI/180.*PSI
      CPSI=COS(PSI)
      SPSI=SIN(PSI)
      psi=tpsi
      Z0=R1*2.*SPSI*CPSI/(3.+SPSI**2)


C******************************************************************
C      definition of the TETAM, STM, and CTMf for Region 1 FAC 
C      3,912353=31.2/7.97474; 7974.74 MWb=2*31200*10**(-9)*PI*RE**2 
C******************************************************************
          STM2=3.912353*FLUX/abs(BD)*1.e-6 
	  CTM=SQRT(1.-STM2)
	  STM=SQRT(STM2)
          TM1=ASIN(STM)*180./pi
          ami1=ami1/2.*(1.+CTM)
	  ami2=0.75*ami1
	  tm2=tm1+3.                                                   

      r0=4.
      bk1=br
      If (r2.GE.6.) then
      rkm=6.
      else 
      rkm=r2
      end if
      al1=sqrt(1.+(2.*rkm+1.)/r1)

         call masring
         b0=bd+bd1*bka
         bd0=bd/b0
      CALL MAS2D
      P=B0/R1/R1
      DO 6 I=1,6
      P=P/R1
!      P1=P1/R1
      A1(I)=CPSI*(D1(I)*P)
      A1(I+6)=SPSI*(D1(I+6)*P)
    6 CONTINUE
      return
      end

      SUBROUTINE BDIPC(x,BM,B)                                           
****************************************************************************
* Program calculating the dipole magnetic field in Cartesian coordinates.  *
*									   *
* INPUT PARAMETERS:  x(3) is a coordinates of the point where the magnetic *
*                    field is being calculated, in GSM coordinates, Re;	   *
*                    BM   is dipole moment, nT*Re^3.			   *
*                    SPSI,CPSI are Sin and Cos of dipole tilt angle. 	   *
* OUTPUT PARAMETERS: B(3) is the magnetic field in GSM coordinates, nT.    *
* Written by V. Kalegaev. 						   * 
****************************************************************************
      dimension x(3),b(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      x1=x(1)
      x2=x(2)
      x3=x(3)
      r2=(x1*x1+x2*x2+x3*x3)
      r=sqrt(r2)      
      r5=r2*r2*r
      p=x3*cpsi-x1*spsi
      br=bm/r5
      b(1)=-Br*(-r*r*spsi-3.*x1*p)
      b(2)= Br*(3.*x2*p)
      b(3)=-Br*(r*r*cpsi-3.*x3*p)
      RETURN                                                            
      END                                                               

      SUBROUTINE BDIP (P,bdp)
***********************************************************************
* Program calculating the magnetic field of the geomagnetic dipole    *
* in spherical coordinates (OX is polar axes) in the current point,   *
* defined by /cor2/, /cor3/ common blocks.			      *
*								      *
* INPUT PARAMETERS: 						      *
*     BDp  is the dipole magnetic moment, nT*Re^3;		      *
*     SPSI, CPSI - are Sin and Cos of dipole tilt angle. 	      *
* OUTPUT PARAMETERS: 						      *
*     P(3) is the dipole magnetic field, nT.			      *
* Written by I. Alexeev 					      *
***********************************************************************
      COMMON/COR2/CFI,SFI
      COMMON/COR3/R,CT,ST
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,bd,bd0
      DIMENSION P(3)
      RR=R*R
      T=Bdp/R/RR
      TC=T*CPSI
      TS=T*SPSI
      CPF=TC*CFI
       P(1)= 2*(CPF*ST-TS*CT)
       P(2)= -CPF*CT-TS*ST
       P(3)= TC*SFI
      RETURN
      END

      SUBROUTINE MASRING
***************************************************************************
* Program calculating the magnetic moment and internal coefficients 	  *
* describing the magnetic field of the ring current in  		  *
* the paraboloid model. 						  *
*									  *
* INPUT PARAMETERS: 							  *
*     BK1  is the maximum intensity of ring current, nT;		  *
*     R0   is the distance to the ring current maximum, Re;		  *
*     RKM  is the distance to the ring current edge, Re.		  *
* OUTPUT PARAMETERS: 							  *
*     BD1     is the ring current magnetic moment, nT*Re^3;		  *
*     BD1*BKA is the dipole, screening ring current, magnetic moment, 	  *
*             nT*Re^3;							  *
* WARNING: You should call MASRING after each changes in the model input  *
*          parameters.							  *
* Written by I. Alexeev 					      	  *
***************************************************************************
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      RK =RKM
      RK3 = RK*RK*RK
      T = RK/R0/2
      T2 = T*T
      T5 = T2*T2*T
      V = 1+T2
      TS = SQRT(V)
      A = T5/V/V/TS
      BKA = A
      BKB = A/RK3/T2
      BKC = 4*R0*R0
      BD1 = BK1*RK3*T2/2/(T5-A)
      RETURN
      END

                              SUBROUTINE BRING (P)
C*********************************************************************
C Calculation of the ring current field                              *
C New version 29.08.2001					     *
C                                                                    *
C    Written by Igor I. Alexeev 				     *
C*********************************************************************
c      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD,bd0
       COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,BKC
      COMMON /COR2/CFI,SFI
      COMMON/COR3/R,CT,ST
      DIMENSION P(3),PD(3),UFR(3)
      T=BD1*BKA/R/R/R
      P(1)=2.*T*(CPSI*CFI*ST-SPSI*CT)
      P(2)=-T*(CPSI*CT*CFI+SPSI*ST)
      P(3)=SFI*T*CPSI
      IF (R.GT.RKM) then
      RETURN
      else
      CALL BDIP (PD,BD1)
       RR=R*R
       RKT2= RR+BKC
       RKT = SQRT(RKT2)
       T2=RR/RKT2
       T3=T2*R/RKT
       TB=BKB*RR*R
       TC=T3*BKC/RKT2

        F1 = T3-TB-BKA
        F2 = F1+3*(TB-TC)

       P(1)=P(1)+PD(1)*F1
       P(2)=P(2)+PD(2)*F2
       P(3)=P(3)+PD(3)*F2
       endif
                             RETURN
                              END


      SUBROUTINE BRING1 (P,f1,f2)
***********************************************************************
* Program calculating the magnetic field of    Bring-BDR*Bdip         *
* in spherical coordinates (BDR=Bd1*BKA=B0-BD)			      *
*								      *
* INPUT PARAMETERS: 						      *
*     R0  - distance to ring current maximum;			      *
*     BK1 - maximum ring current magnetic field; 		      *
*     BD1,RKM,BKA,BKB,bkc - are calculated in MASRING subroutine.     *
* OUTPUT PARAMETERS: 						      *
*     P(3) is the  "Bring-BDR*Bdip"  magnetic field, nT.	      *
* Written by I. Alexeev 					      *    
***********************************************************************
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      COMMON/COR3/R,CT,ST
      DIMENSION P(3),PD(3)

***Calculation of the Magnetic Fields of Geodipole and Ring Current 

       CALL BDIP (PD,BD1)
       RR=R*R
       RKT2= RR+BKC
       RKT = SQRT(RKT2)
       T2=RR/RKT2
       T3=T2*R/RKT
       TB=BKB*RR*R
       TC=T3*BKC/RKT2
        F1 = T3-TB-BKA
        F2 = F1+3*(TB-TC)
       P(1)= PD(1)*F1
       P(2)= PD(2)*F2
       P(3)= PD(3)*F2
      RETURN
      END

      SUBROUTINE BEG(UF,VV)
****************************************************************
* Calculation of the summary dipole, ring current and Chapman- *
* Ferraro currents at the magnetopause magnetic field in       *
* the inner magnetosphere (al<al1)			       *
* Written by I. Alexeev 				       *
****************************************************************
      COMMON /COR2/CFI,SFI	
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/COR3/R,CT,ST
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T1/A1(12)
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      COMMON/begf/UFCF(3)
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      DIMENSION UF(3),VV(3,3),ufr(3)
      DIMENSION EL(7),EL1(7),EL2(7)
      EL(1)=1.
      EL1(1)=1.
      EL(2)=CT
      EL1(2)=3.*CT
      EL2(1)=0.
      EL2(2)=3.
      EL2(3)=15.*CT
      DO 19 I=2,6
      P=1.-1./I
      EL(I+1)=(1.+P)*CT*EL(I)-P*EL(I-1)
19    EL1(I+1)=(3.-P)*CT*EL1(I)-(2.-P)*EL1(I-1)
      DO 20 I=3,6
      P=1.+3./(I-1)
      EL2(I+1)=(1.+P)*CT*EL2(I)-P*EL2(I-1)
   20 CONTINUE

***Calculation of the Magnetic field of Magnetopause Screening Currents 
      T=A1(6)
      U3=6.*T*EL1(7)
      V=EL1(6)
      U2=T*V
      U1=6.*U2
      U4=36.*T*V
      U5=6.*U3
      U6=T*EL2(7)
      U7=T*EL2(6)
      T=A1(12)
      V1=6.*T*EL(7)
      V2=T*V
      V3=6.*V1
      V4=6.*V2
      V5=6.*T*EL1(7)
      I=5

21    CONTINUE
      T=A1(I)
      U3=U3*R+I*T*EL1(I+1)
      U1=U1*R+I*T*EL1(I)
      U2=U2*R+T*EL1(I)
      U4=U4*R+I*T*I*EL1(I)
      U5=U5*R+I*I*T*EL1(I+1)
      U6=U6*R+T*EL2(I+1)
      U7=U7*R+T*EL2(I)
      T=A1(I+6)
      V1=V1*R+I*T*EL(I+1)
      V2=V2*R+T*EL1(I)
      V3=V3*R+I*T*EL(I+1)*I
      V4=V4*R+I*T*EL1(I)
      V5=V5*R+I*T*EL1(I+1)
      I=I-1
      IF(I.GT.0)GOTO 21
      UFCF(1)=-ST*CFI*U1+V1
      UFCF(2)=CFI*(CT*(U1+U2)-U3)-V2*ST
      UFCF(3)=SFI*U2
      do 24 i=1,3
 24   uf(i)=ufcf(i)
      RETURN
      END

      SUBROUTINE FLYD(P,BB)
******************************************************************
* Calculation of the summary dipole, ring current and Chapman-	 *
* Ferraro currents at the magnetopause magnetic field in	 *
* the geomagnetic tail (al>al1) 				 *
* Written by I. Alexeev 					 *     
******************************************************************
      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      REAL L,L0
      COMMON/T3/ L(6,5),L0(5)
      COMMON /COR1/AL,BE,SQ,PQ,QA
      COMMON /COR2/CFI,SFI
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      DIMENSION P(3),BB(3,3)
      A2=AL*AL
      B2=BE*BE
      Y1=0.
      Y2=0.
      Y3=0.
      Y4=0.
      Y5=0.
      Y6=0.
      Y7=0.
      Y8=0.
      Y9=0.
      DO 2 N=1,5
      Z=L(1,N)
      CALL BESS(1,Z*BE,U,DU)
      CALL BESK(1,Z*AL,V,DV)
      X=CF2(N)
      Y1=Y1+X*U*DV
      Y2=Y2+X*DU*V
      Y3=Y3+CF1(N)*U*V
      X=CF3(N)
      Y4=Y4+X*U*V
      Y5=Y5+X*DU*DV
      Z=L0(N)
      X=Z*BE
      Z=Z*AL
      DU=CF0(N)
      DV=CF4(N)
      U=BESK0(Z)
      Z=BESK1(Z)
      V=BESJ0(X)
      X=BESJ1(X)
      Y6=Y6+DU*Z*V
      Y7=Y7+DU*U*X
      Y8=Y8+DV*V*U
      Y9=Y9+DV*Z*X
    2 CONTINUE
      P(1)=(-Y1*CFI+Y6)/SQ
      P(2)=(-Y2*CFI+Y7)/SQ
      P(3)=Y3*SFI/PQ
      RETURN
      END
C
      SUBROUTINE BFAC (P)
C*************************************************************
C Calculation of the magnetic field of Region 1 FAC (P(3)).
C R,CTE,STE,CFIE,SFIE are the geomagnetic spherical coordinates 
C (polar axis is directed on the Morth magnetic pole)
C STM,CTM are sin(tetam), cos(tetam),
C BFACP0=0,000098*AJ0/STM , 
C AJ0 is total field aligned current in hemisphere
C BFACP1=BFACP0*(1-CTM).
C Written by I. Alexeev
C*************************************************************
      COMMON/COR3/R,CT,ST
      COMMON/COR4/CTE,STE,CFIE,SFIE
      COMMON/TFAC/STM,CTM,BFAC0,BFAC1,TETAM,AJ0
c      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD,bd0
       DIMENSION P(3)

        U= TETAM*PI/180.
        STM=SIN(U)
        CTM=COS(U)
        BFAC0=0.000098*AJ0/STM
        BFAC1=BFAC0*(1-CTM)


      IF (R.GT.1.AND.R.LT.R1) GOTO 3
       P(1)= 0.
       P(2)= 0.
       P(3)= 0.
                              RETURN
3      CONTINUE

      U=BFAC0/R
      U1=BFAC1/R
      IF (STE.GT.STM)        GOTO 2
      IF (CTE.GT.0)          GOTO 1
C***************************************************
C               South polar cap	 
C***************************************************
       U3=U/(1-CTE)
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= U3*SFIE

                              RETURN
C***************************************************
C            North polar cap 
C***************************************************
1      CONTINUE
       U3=U/(1+CTE)
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= -U3*SFIE

                              RETURN
C***************************************************
C              Inner magnetosphere
C***************************************************
2      CONTINUE
       U3=U1/STE/STE
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= U3*SFIE*CTE
                              RETURN
                              END
			    

      SUBROUTINE DERY4D(P,BB)
***************************************************************
* Calculation of the magnetic field of geomagnetic tail
* current system
* Written by I. Alexeev 					 *     
***************************************************************
      REAL L,L0
      COMMON /COR1/AL,BE,SQ,PQ,QA
      COMMON /COR2/CFI,SFI
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T3/ L(6,5),L0(5)
      COMMON/S1/CB(6,5),CB2(6,5),
     +CD(6,5),CB3(6,5),CD2(6,5),CD3(6,5)
      DIMENSION P(3),BB(3,3),U(6,5),DU(6,5)

      if (sbt.eq.1.)then 
      IF(AL-AL0)2,2,3
2     CONTINUE
      IF(AL-174.673)40,40,41
41    PRINT 42,AL
42    FORMAT(19H GRAND EXP-DERY,AL=,E12.5)
      AL=174.670
40    CONTINUE
      E4=EXP(AL)
      DO 21 K=1,6
      M=2*K-1
      DO 21 N=1,5
      X =AL *L(K,N)
      Y=X -16.118095651
      E4=EXP(Y)
      CALL BESM(M,X,Z,DZ)
      U(K,N)=Z*E4
      DU(K,N)=DZ*E4
   21 CONTINUE
      W=0.
      GO TO 4
3     DO 22 K=1,6
      M=2*K-1
      DO 22 N=1,5
      X=AL*L(K,N)
      CALL BESK(M,X,Z,DZ)
      U(K,N)=Z*E5
      DU(K,N)=DZ*E5
22    CONTINUE
      W=+SIGN(1.,CFI)*C0/AL**2
      IF (CFI.EQ.0.) W=0.
4     R=+1.
      V2=0.
      V1=0.
      V3=0.
      V4=0.
      V5=0.
      V6=0.
      V7=0.
      V8=0.
      V9=0.
      BE1=BE/BETA0
      CVI=CFI
      SVI=-SFI
      SVS=2.*SFI*CFI
      CVS=2.*CFI**2-1.
      DO 23 K=1,6
      M=2*K-1
      Y1=0.
      Y2=0.
      Y3=0.
      Y4=0.
      Y5=0.
      CMFI=CVI*CVS-SVI*SVS
      SMFI=SVI*CVS+CVI*SVS
      CVI=CMFI
      SVI=SMFI
      DO 24 N=1,5
      X=L(K,N)*BE1
      CALL BESS(M,X,Z,DZ)
      IF(AL-AL0)5,5,6
    5 X1=CB(K,N)
      X2=CB2(K,N)
      X3=CB3(K,N)
      GO TO 7
6     X1=CD(K,N)
      X2=CD2(K,N)
      X3=CD3(K,N)
7     Y1=Y1+X2*Z*DU(K,N)
      Y2=Y2+X2*DZ*U(K,N)
      Y3=Y3+X1*Z*U(K,N)
      Y4=Y4+X3*Z*U(K,N)
      Y5=Y5+X3*DZ*DU(K,N)
   24 CONTINUE
      X1=CMFI/M*R
      X2=CMFI*M*R
      X3=SMFI*R
      V1=V1+Y1*X1
      V2=V2+Y2*X1
      V3=V3+Y3*X3
      V4=V4+Y4*X1
      V5=V5+Y3*X2
      V6=V6+Y2*X1
      V7=V7+Y1*X3
      V8=V8+Y2*X3
      R=-R
   23 CONTINUE
      IF(AL-AL0)60,60,70
60    V1=V1*E5
      V2=V2*E5
      V3=V3*E5
      V4=V4*E5
      V5=V5*E5
      V6=V6*E5
      V7=V7*E5
      V8=V8*E5
70    CONTINUE
      P(1)=(-V1-W*AL)/SQ
      P(2)=-V2/SQ
      P(3)=V3/PQ
      
      else
      do 1 i=1,3
      p(i)=0.
      do 1 j=1,3
1     bb(j,i)=0.
      end if 
      RETURN
      END
C
      SUBROUTINE PERE(UF,VV,DET)
      DIMENSION UF(3),VV(3),DET(3,3)
      DO 1 I=1,3
      P=0.
      DO 2 J=1,3
      P=P+DET(I,J)*UF(J)
    2 CONTINUE
      VV(I)=P
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE PRIS(UF,VV)
      COMMON /COR1/AL,BE,SQ,PQ,QA
      COMMON /COR2/CFI,SFI
      DIMENSION UF(3,3),VV(3,3)
      UF(1,1)=-AL/SQ
      UF(1,2)=BE/SQ
      UF(1,3)=0.
      Z=SFI/SQ
      UF(2,1)=BE*Z
      UF(2,2)=AL*Z
      UF(2,3)=CFI
      Z=CFI/SQ
      UF(3,1)=BE*Z
      UF(3,2)=AL*Z
      UF(3,3)=-SFI
      DO 1 I=1,3
      DO 1 J=1,3
      VV(I,J)=UF(J,I)
    1 CONTINUE
      RETURN
      END

      BLOCK DATA
      REAL L,L0
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T3/L(6,5),L0(5)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      COMMON/S5/CI0(9),CI1(9)
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      DATA
     *L0/3.83170597,7.01558667,10.17346814,13.3236919,16.47063005/,
     *L/1.84118390,4.2011889412,6.4156163752,8.5778364889,10.711433969,
     +12.826491226,5.3314427000,8.0152365984,10.519860874,12.932386237,
     *15.286737667,17.600266557,8.5363163000,11.345924311,13.987188630,
     *16.529365884,19.004593538,21.430854238,11.706005000,14.585848286,
     *17.312842488,19.941853366,22.501398726,25.008518704,14.863588700,
     *17.788747866,20.575514521,23.268052926,25.891277276,28.460857279
     */
      DATA
     *HN/100./,R1/10./,ON/10.4372/,B0/-31200./,PI/3.1415926/,BT/40./
      DATA SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
     */1.,0.,1.,1.,1.,1.,1.,1.,1.,0./
      DATA CI0/
     *+0.003923767,-0.016476329,+0.026355372,
     *-0.020577063,+0.009162808,-0.001575649,
     *+0.002253187,+0.013285917,+0.398942280/
      DATA CI1/
     *-0.004200587,+0.017876535,-0.028953121,
     *0.022929673,-0.010315550,+0.001638014,
     *-0.003620183,-0.039880242,+0.398942280/
      END

      SUBROUTINE MAS2D
* Written by I. Alexeev 					 *     
      REAL L,L0
      COMMON/T3/L(6,5),L0(5)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/S1/ CB(6,5),CB2(6,5),
     *CD(6,5),CB3(6,5),CD2(6,5),
     +CD3(6,5)
      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      COMMON/S5/CI0(9),CI1(9)
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      DIMENSION CF00(5),CF01(5),
     *cb0(6,5),cd0(6,5)
        DATA                      CB0/
     *7.493663380000,3.119532240E-1,1.706022280E-2,1.03359770E-3,
     *6.633343340E-5,4.421213830E-6,5.777707910E-3,2.103707660E-4,
     *7.538502750E-6,3.011267810E-7,1.313975350E-8,6.13862642E-10,
     *1.867244830E-5,5.965184960E-7,1.712466940E-8,5.43331737E-10,
     *1.89762765E-11,7.17571687E-13,8.545308210E-8,2.571561030E-9,
     *6.45145558E-11,1.75563377E-12,5.24527385E-14,1.70163718E-15,
     *4.43194331E-10,1.26958726E-11,2.88578820E-13,6.99784044E-15,
     *1.85385885E-16,0.
     */,

     *                             CD0/
     *5.40613258E-5,1.17552607E-3,2.01041928E-2,0.31415060400,
     *4.67288038000,67.3370338000,2.50590102E-3,0.202158093,
     *7.58793123,216.779913,5328.09218,118727.936,
     *1.72159326E-1,21.3628656,1157.19802,45465.0611,
     *1482947.1,42673478.6,14.734532, 2354.2227,
     *161804.316,7916120.89,316102914.,10974469500.,
     *1367.45791,254342.53,20562943.2,1179904930.,
     *54901114700.,2205064660000./
      DATA CF00/-94630.534652,-10798057.06,-652953011.
     +92,-30158380850.,-1198163932400.
     */
      DATA
     *CF01/-1318.693367,-190097.98824,-96033
     +82.4047,-369794441.86,-12502247266.
     */
      P=((10./R1)**3)
      P=P*B0/(-31200.)
      DO 1 N=1,5
      CF0(N)=CF00(N)*SPSI*P
      CF1(N)=CF01(N)*CPSI*P
      CF2(N)=CF1(N)*L(1,N)
      CF3(N)=CF2(N)*L(1,N)
      CF4(N)=CF0(N)*L0(N)
    1 CONTINUE
      AL01=SQRT(2.4)
    		      
cc	print *, 'mas',
cc     *cb(1,1),cb0(1,1),ska,sk,al0,al01,bt
		      
      DO 2 K=1,6
      M=2*K-1
      DO 2 N=1,5
      ZZ=L(k,N)*AL01
      ZA=L(k,N)*AL0
      CALL BESK(M,ZZ,SK ,DSK )
      CALL BESM(M,ZZ,SI ,DSI )
      CALL BESK(M,ZA,SKA,DSKA)
      CALL BESM(M,ZA,SIA,DSIA)
      SIII=ZA-ZZ
      SIA=SIA/SI*EXP(SIII)
      CB(K,N)=BT/40.*CB0(K,N)*AL0*SKA/AL01/SK
      CD(K,N)=BT/40.*CD0(K,N)*AL0*SIA/AL01
      CB2(K,N)=L(K,N)*CB(K,N)
      CB3(K,N)=L(K,N)*CB2(K,N)
      CD2(K,N)=L(K,N)*CD(K,N)
      CD3(K,N)=L(K,N)*CD2(K,N)
    2 CONTINUE
      AL1=SQRT(2.4)
      C0=   61.96773354*BT/40.*al0/al01
      E5=10.**7
      BETA0=1.0+0.2*(PSI/30.)**2
*--------------------------------------------------------------
*********> place to insert block "test.for"<*******************
*--------------------------------------------------------------
      RETURN
      END
C
      SUBROUTINE BESK(M,X,V,DV)
* Written by I. Alexeev 					 *     
      U=BESK0(X)
      V=BESK1(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=2.*K*V/X+U
      U=V
2     V=S
4     DV=-(M*V/X+U)
      RETURN
      END
C
      FUNCTION RINT(DZ,FX,FZ)
      DX=FX*DZ/FZ
      RINT=DX
      return
      END
C
      SUBROUTINE BESM(M,X,V,DV)
* Written by I. Alexeev 					 *     
      COMMON/S5/CI0(9),CI1(9)
      IF(X.GT.3.75)GO TO 5
      XA=-X
      IF(XA-174.673)50,50,51
 51   PRINT 52,X
52    FORMAT(18H GRAND EXP-BESM,X=,E12.5)
      X=-174.670
50    CONTINUE
      E=EXP(-X)
      V=BSI(M,X)*E
      DV=0.5*E*(BSI(M-1,X)+BSI(M+1,X))
      RETURN
5     CONTINUE
      IF(X)90,91,91
90    PRINT 92,X
92    FORMAT(19H NEGATIVE X-BESM,X=,E12.5)
   91 X=ABS(X)
      U=UG(CI0(1),X)/SQRT(X)
      V=UG(CI1(1),X)/SQRT(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=U-2.*K*V/X
      U=V
2     V=S
4     DV=U-M*V/X
      RETURN
      END
C
      SUBROUTINE BESS(M,X,V,DV)
* Written by I. Alexeev 					 *     
      IF(X.GT.3.75)GO TO 5
      V=BSJ(M,X)
      DV=0.5*(BSJ(M-1,X)-BSJ(M+1,X))
      RETURN
5     U=BESJ0(X)
      V=BESJ1(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=2*K*V/X-U
      U=V
2      V=S
4     DV=U-M*V/X
      RETURN
      END
C
      FUNCTION BSJ(N,X)
* Written by I. Alexeev 					 *     
      SUM=1.
      P=1.
      Z=-X*X/4.
      DO 2 K=1,7
      P=P*Z/K/(K+N)
2     SUM=SUM+P
      IF(N.LE.0)GO TO 3
      DO 1 K=1,N
1     SUM=SUM/K
  3   CONTINUE
7     FORMAT(16H EXP NEGATIVE,N=,I3,2HX=,E12.5)
      IF(X)4,5,4
5     CONTINUE
      BSJ=0.
      PRINT 7,N,X
      GO TO 6
4     CONTINUE
      BSJ=SUM*(X/2.)**N
6     RETURN
      END
C
      FUNCTION BSI(N,X)
* Written by I. Alexeev 					 *     
      SUM=1.
      P=1.
      Z=X*X/4.
      DO 2 K=1,7
      P=P*Z/K/(K+N)
2     SUM=SUM+P
      IF(N.LE.0)GO TO 3
      DO 1 K=1,N
1     SUM=SUM/K
  3   CONTINUE
7     FORMAT(16H EXP NEGATIVE,N=,I3,2HX=,E12.5)
      IF(X)4,5,4
5     CONTINUE
      PRINT 7,N,X
      BSI=0.
      GO TO 6
4     CONTINUE
      BSI=SUM*(X/2.)**N
6     RETURN
      END
C
      REAL FUNCTION UG(V,X)
* Written by I. Alexeev 					 *     
      DIMENSION V(9)
      UG=V(1)
      DO 7 I=2,9
      UG=UG*(3.75/X)+V(I)
7     CONTINUE
      RETURN
      END
C
C
      FUNCTION BESJY(X)
C*************************************************************************
C Calculation of the Bessel functions J0(x), J1(x), Y0(x) or Y1(x)       *
C   INPUT:                                                               * 
C        x - argument of the Bessel function                             * 
C                                                                        * 
C   OUTPUT: 					         		 *
C	BESJ0(X)= J0(x)                                           	 *
C	BESJ1(X)= J1(x)                                                  *  
C	BESY0(X)= Y0(x)                                           	 *
C	BESY1(X)= Y1(x)                                                  *  
C   Special subroutine from MSU Computer Center library  		 *
C*************************************************************************

      LOGICAL L
C
      ENTRY BESJ0(X)
C
      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =           - 0.00000 00000 000008
      B = F * A     + 0.00000 00000 000413
      A = F * B - A - 0.00000 00000 019438
      B = F * A - B + 0.00000 00000 784870
      A = F * B - A - 0.00000 00026 792535
      B = F * A - B + 0.00000 00760 816359
      A = F * B - A - 0.00000 17619 469078
      B = F * A - B + 0.00003 24603 288210
      A = F * B - A - 0.00046 06261 662063
      B = F * A - B + 0.00481 91800 694676
      A = F * B - A - 0.03489 37694 114089
      B = F * A - B + 0.15806 71023 320973
      A = F * B - A - 0.37009 49938 726498
      B = F * A - B + 0.26517 86132 033368
      A = F * B - A - 0.00872 34423 528522
      A = F * A - B + 0.31545 59429 497802
      BESJY=0.5*(A-B)
      IF(L) RETURN
C
      A =           + 0.00000 00000 000016
      B = F * A     - 0.00000 00000 000875
      A = F * B - A + 0.00000 00000 040263
      B = F * A - B - 0.00000 00001 583755
      A = F * B - A + 0.00000 00052 487948
      B = F * A - B - 0.00000 01440 723327
      A = F * B - A + 0.00000 32065 325377
      B = F * A - B - 0.00005 63207 914106
      A = F * B - A + 0.00075 31135 932578
      B = F * A - B - 0.00728 79624 795521
      A = F * B - A + 0.04719 66895 957634
      B = F * A - B - 0.17730 20127 811436
      A = F * B - A + 0.26156 73462 550466
      B = F * A - B + 0.17903 43140 771827
      A = F * B - A - 0.27447 43055 297453
      A = F * A - B - 0.06629 22264 065699
      BESJY=0.636619772367581*ALOG(X)*BESJY+0.5*(A-B)
      RETURN
C
    4 F=256.0/X**2-2.0
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000051
      B = F * A - B + 0.00000 00000 000433
      A = F * B - A - 0.00000 00000 004305
      B = F * A - B + 0.00000 00000 051683
      A = F * B - A - 0.00000 00000 786409
      B = F * A - B + 0.00000 00016 306465
      A = F * B - A - 0.00000 00517 059454
      B = F * A - B + 0.00000 30751 847875
      A = F * B - A - 0.00053 65220 468132
      A = F * A - B + 1.99892 06986 950373
      P=A-B
      B =           - 0.00000 00000 000006
      A = F * B     + 0.00000 00000 000043
      B = F * A - B - 0.00000 00000 000334
      A = F * B - A + 0.00000 00000 003006
      B = F * A - B - 0.00000 00000 032067
      A = F * B - A + 0.00000 00000 422012
      B = F * A - B - 0.00000 00007 271916
      A = F * B - A + 0.00000 00179 724572
      B = F * A - B - 0.00000 07414 498411
      A = F * B - A + 0.00006 83851 994261
      A = F * A - B - 0.03111 17092 106740
      Q=8.0*(A-B)/V
      F=V-0.785398163397448
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 6
      BESJY=F*(Q*A+P*B)
      RETURN
    6 BESJY=F*(P*A-Q*B)
      RETURN
C
      ENTRY BESJ1(X)
C
      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 5
    3 F=0.0625*X**2-2.0
      B =           + 0.00000 00000 000114
      A = F * B     - 0.00000 00000 005777
      B = F * A - B + 0.00000 00000 252812
      A = F * B - A - 0.00000 00009 424213
      B = F * A - B + 0.00000 00294 970701
      A = F * B - A - 0.00000 07617 587805
      B = F * A - B + 0.00001 58870 192399
      A = F * B - A - 0.00026 04443 893486
      B = F * A - B + 0.00324 02701 826839
      A = F * B - A - 0.02917 55248 061542
      B = F * A - B + 0.17770 91172 397283
      A = F * B - A - 0.66144 39341 345433
      B = F * A - B + 1.28799 40988 576776
      A = F * B - A - 1.19180 11605 412169
      A = F * A - B + 1.29671 75412 105298
      BESJY=0.0625*(A-B)*X
      IF(L) RETURN
C
      B =           - 0.00000 00000 000244
      A = F * B     + 0.00000 00000 012114
      B = F * A - B - 0.00000 00000 517212
      A = F * B - A + 0.00000 00018 754703
      B = F * A - B - 0.00000 00568 844004
      A = F * B - A + 0.00000 14166 243645
      B = F * A - B - 0.00002 83046 401495
      A = F * B - A + 0.00044 04786 298671
      B = F * A - B - 0.00513 16411 610611
      A = F * B - A + 0.04231 91803 533369
      B = F * A - B - 0.22662 49915 567549
      A = F * B - A + 0.67561 57807 721877
      B = F * A - B - 0.76729 63628 866459
      A = F * B - A - 0.12869 73843 813500
      A = F * A - B + 0.04060 82117 718685
      BESJY=0.636619772367581*ALOG(X)*BESJY-0.636619772367581/X
     1     +0.0625*(A-B)*X
      RETURN
C
    5 F=256.0/X**2-2.0
      B =           - 0.00000 00000 000007
      A = F * B     + 0.00000 00000 000055
      B = F * A - B - 0.00000 00000 000468
      A = F * B - A + 0.00000 00000 004699
      B = F * A - B - 0.00000 00000 057049
      A = F * B - A + 0.00000 00000 881690
      B = F * A - B - 0.00000 00018 718907
      A = F * B - A + 0.00000 00617 763396
      B = F * A - B - 0.00000 39872 843005
      A = F * B - A + 0.00089 89898 330859
      A = F * A - B + 2.00180 60817 200274
      P=A-B
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000046
      B = F * A - B + 0.00000 00000 000360
      A = F * B - A - 0.00000 00000 003264
      B = F * A - B + 0.00000 00000 035152
      A = F * B - A - 0.00000 00000 468636
      B = F * A - B + 0.00000 00008 229193
      A = F * B - A - 0.00000 00209 597814
      B = F * A - B + 0.00000 09138 615258
      A = F * B - A - 0.00009 62772 354916
      A = F * A - B + 0.09355 55741 390707
      Q=8.0*(A-B)/V
      F=V-2.356194490192345
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 7
      BESJY=F*(Q*A+P*B)
      RETURN
    7 BESJY=F*(P*A-Q*B)
      IF(X .LT. 0.0) BESJY=-BESJY
      RETURN
C
      ENTRY BESY0(X)
C
      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 4
      GO TO 8
C
      ENTRY BESY1(X)
C
      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 5
      GO TO 3
C
    9 BESJY=0.
      PRINT 100,X
      RETURN
100   FORMAT(1X,32HBESJY...NON-POSITIVE ARGUMENT X=,E12.5)
C
      END


      FUNCTION BESIK(X)
C*************************************************************************
C Calculation of the Bessel functions K0(x), K1(x), I0(x) or I1(x)       *
C   INPUT:                                                               * 
C        x - argument of the Bessel function                             * 
C                                                                        * 
C   OUTPUT: 					         		 *
C	BESI0(X)= I0(x)                                           	 *
C	BESI1(X)= I1(x)                                                  *  
C	BESK0(X)= K0(x)                                           	 *
C	BESK1(X)= K1(x)                                                  *  
C	EBESI0(X)= exp(-|x|)*I0(x)                               	 *
C	EBESI1(X)= exp(-|x|)*I1(x)                                       *  
C	EBESK0(X)= exp(x)*K0(x)                                 	 *
C	EBESK1(X)= exp(x)*K1(x)                                          *  
C   Special subroutine from MSU Computer Center library         	 *
C*************************************************************************
      LOGICAL L,E
C
      ENTRY EBESI0(X)
C
      E=.TRUE.
      GO TO 1
      ENTRY BESI0(X)
      E=.FALSE.
    1 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =               0.00000 00000 00002
      B = F * A     +   0.00000 00000 00120
      A = F * B - A +   0.00000 00000 06097
      B = F * A - B +   0.00000 00002 68828
      A = F * B - A +   0.00000 00101 69727
      B = F * A - B +   0.00000 03260 91051
      A = F * B - A +   0.00000 87383 15497
      B = F * A - B +   0.00019 24693 59688
      A = F * B - A +   0.00341 63317 66012
      B = F * A - B +   0.04771 87487 98174
      A = F * B - A +   0.50949 33654 39983
      B = F * A - B +   4.01167 37601 79349
      A = F * B - A +  22.27481 92424 62231
      B = F * A - B +  82.48903 27440 24100
      A = F * B - A + 190.49432 01727 42844
      A = F * A - B + 255.46687 96243 62167
      BESIK=0.5*(A-B)
      IF(L .AND. E) BESIK=EXP(-V)*BESIK
      IF(L) RETURN
      A =           +   0.00000 00000 00003
      B = F * A     +   0.00000 00000 00159
      A = F * B - A +   0.00000 00000 07658
      B = F * A - B +   0.00000 00003 18588
      A = F * B - A +   0.00000 00112 81211
      B = F * A - B +   0.00000 03351 95256
      A = F * B - A +   0.00000 82160 25940
      B = F * A - B +   0.00016 27083 79043
      A = F * B - A +   0.00253 63081 88086
      B = F * A - B +   0.03008 07224 20512
      A = F * B - A +   0.25908 44324 34900
      B = F * A - B +   1.51153 56760 29228
      A = F * B - A +   5.28363 28668 73920
      B = F * A - B +   8.00536 88687 00334
      A = F * B - A -   4.56343 35864 48395
      A = F * A - B -  21.05766 01774 02440
      BESIK=-ALOG(0.125*X)*BESIK+0.5*(A-B)
      IF(E) BESIK=EXP(X)*BESIK
      RETURN
    4 F=32.0/V-2.0
      B =           - 0.00000 00000 00001
      A = F * B     - 0.00000 00000 00001
      B = F * A - B + 0.00000 00000 00004
      A = F * B - A + 0.00000 00000 00010
      B = F * A - B - 0.00000 00000 00024
      A = F * B - A - 0.00000 00000 00104
      B = F * A - B + 0.00000 00000 00039
      A = F * B - A + 0.00000 00000 00966
      B = F * A - B + 0.00000 00000 01800
      A = F * B - A - 0.00000 00000 04497
      B = F * A - B - 0.00000 00000 33127
      A = F * B - A - 0.00000 00000 78957
      B = F * A - B + 0.00000 00000 29802
      A = F * B - A + 0.00000 00012 38425
      B = F * A - B + 0.00000 00085 13091
      A = F * B - A + 0.00000 00568 16966
      B = F * A - B + 0.00000 05135 87727
      A = F * B - A + 0.00000 72475 91100
      B = F * A - B + 0.00017 27006 30778
      A = F * B - A + 0.00844 51226 24921
      A = F * A - B + 2.01655 84109 17480
      BESIK=0.199471140200717*(A-B)/SQRT(V)
      IF(E) RETURN
      BESIK=EXP(V)*BESIK
      RETURN
      ENTRY EBESI1(X)
      E=.TRUE.
      GO TO 2
      ENTRY BESI1(X)
      E=.FALSE.
    2 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 3
    7 F=0.0625*X**2-2.0
      A =           +   0.00000 00000 00001
      B = F * A     +   0.00000 00000 00031
      A = F * B - A +   0.00000 00000 01679
      B = F * A - B +   0.00000 00000 79291
      A = F * B - A +   0.00000 00032 27617
      B = F * A - B +   0.00000 01119 46285
      A = F * B - A +   0.00000 32641 38122
      B = F * A - B +   0.00007 87567 85754
      A = F * B - A +   0.00154 30190 15627
      B = F * A - B +   0.02399 30791 47841
      A = F * B - A +   0.28785 55118 04672
      B = F * A - B +   2.57145 99063 47755
      A = F * B - A +  16.33455 05525 22066
      B = F * A - B +  69.39591 76337 34448
      A = F * B - A + 181.31261 60405 70265
      A = F * A - B + 259.89023 78064 77292
      BESIK=0.0625*(A-B)*X
      IF(L .AND. E) BESIK=EXP(-V)*BESIK
      IF(L) RETURN
      A =           +   0.00000 00000 00001
      B = F * A     +   0.00000 00000 00042
      A = F * B - A +   0.00000 00000 02163
      B = F * A - B +   0.00000 00000 96660
      A = F * B - A +   0.00000 00036 96783
      B = F * A - B +   0.00000 01193 67971
      A = F * B - A +   0.00000 32025 10692
      B = F * A - B +   0.00007 00106 27855
      A = F * B - A +   0.00121 70569 94516
      B = F * A - B +   0.01630 00492 89816
      A = F * B - A +   0.16107 43016 56148
      B = F * A - B +   1.10146 19930 04852
      A = F * B - A +   4.66638 70268 62842
      B = F * A - B +   9.36161 78313 95389
      A = F * B - A -   1.83923 92242 86199
      A = F * A - B -  26.68809 54808 62668
      BESIK=ALOG(0.125*X)*BESIK+1.0/X-0.0625*(A-B)*X
      IF(E) BESIK=EXP(X)*BESIK
      RETURN
    3 F=32.0/V-2.0
      B =           + 0.00000 00000 00001
      A = F * B     + 0.00000 00000 00001
      B = F * A - B - 0.00000 00000 00005
      A = F * B - A - 0.00000 00000 00010
      B = F * A - B + 0.00000 00000 00026
      A = F * B - A + 0.00000 00000 00107
      B = F * A - B - 0.00000 00000 00053
      A = F * B - A - 0.00000 00000 01024
      B = F * A - B - 0.00000 00000 01804
      A = F * B - A + 0.00000 00000 05103
      B = F * A - B + 0.00000 00000 35408
      A = F * B - A + 0.00000 00000 81531
      B = F * A - B - 0.00000 00000 47563
      A = F * B - A - 0.00000 00014 01141
      B = F * A - B - 0.00000 00096 13873
      A = F * B - A - 0.00000 00659 61142
      B = F * A - B - 0.00000 06297 24239
      A = F * B - A - 0.00000 97321 46728
      B = F * A - B - 0.00027 72053 60764
      A = F * B - A - 0.02446 74429 63276
      A = F * A - B + 1.95160 12046 52572
      BESIK=0.199471140200717*(A-B)/SQRT(V)
      IF(X .LT. 0.0) BESIK=-BESIK
      IF(E) RETURN
      BESIK=EXP(V)*BESIK
      RETURN
      ENTRY EBESK0(X)
      E=.TRUE.
      GO TO 11
      ENTRY BESK0(X)
      E=.FALSE.
   11 IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 8
      F=20.0/X-2.0
      A =           - 0.00000 00000 00002
      B = F * A     + 0.00000 00000 00011
      A = F * B - A - 0.00000 00000 00079
      B = F * A - B + 0.00000 00000 00581
      A = F * B - A - 0.00000 00000 04580
      B = F * A - B + 0.00000 00000 39044
      A = F * B - A - 0.00000 00003 64547
      B = F * A - B + 0.00000 00037 92996
      A = F * B - A - 0.00000 00450 47338
      B = F * A - B + 0.00000 06325 75109
      A = F * B - A - 0.00001 11066 85197
      B = F * A - B + 0.00026 95326 12763
      A = F * B - A - 0.01131 05046 46928
      A = F * A - B + 1.97681 63484 61652
      BESIK=0.626657068657750*(A-B)/SQRT(X)
      IF(E) RETURN
      BESIK=EXP(-X)*BESIK
      RETURN
      ENTRY EBESK1(X)
      E=.TRUE.
      GO TO 12
      ENTRY BESK1(X)
      E=.FALSE.
   12 IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 7
      F=20.0/X-2.0
      A =           + 0.00000 00000 00002
      B = F * A     - 0.00000 00000 00013
      A = F * B - A + 0.00000 00000 00089
      B = F * A - B - 0.00000 00000 00663
      A = F * B - A + 0.00000 00000 05288
      B = F * A - B - 0.00000 00000 45757
      A = F * B - A + 0.00000 00004 35417
      B = F * A - B - 0.00000 00046 45555
      A = F * B - A + 0.00000 00571 32218
      B = F * A - B - 0.00000 08451 72048
      A = F * B - A + 0.00001 61850 63810
      B = F * A - B - 0.00046 84750 28167
      A = F * B - A + 0.03546 52912 43331
      A = F * A - B + 2.07190 17175 44716
      BESIK=0.626657068657750*(A-B)/SQRT(X)
      IF(E) RETURN
      BESIK=EXP(-X)*BESIK
      RETURN
    9 BESIK=0.
      PRINT 200,X
200   FORMAT(1X,32HBESIK...NON-POSITIVE ARGUMENT X=,E12.5)
      RETURN
      END
C
       subroutine FAC(x,B1,b2)
*******************************************************************
c
C version of 27.08.98
c FAC.for calculation of the magnetic field from field-aligned
c current-I and II using dipole configuration
c with equatorial current.
c     
c         ami1 - total R-I current in MA, 
c         tm1 - colatitude of  R-I current in degrees
c         ami2 - total R-II current in MA, 
c         tm2 - colatitude of  R-II current in degrees
c x(1-3) - GSM coord. of point, B1(1-3) - mag. field GSM coord. (nT),
c                               B2(1-3) - mag. field GSM coord. (nT)
c      B2=0 in this version of the paraboloid model
* Written by V. Kalegaev         					 *     
********************************************************************
      COMMON/COR3/R,CT,ST
      COMMON/COR4/CTE,STE,CFIE,SFIE
      COMMON/TFAC/STM,CTM,BFAC0,BFAC1,TETAM,AJ0
c      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD,bd0
       COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      common/fac12/ami1,ami2,tm1,tm2
      DIMENSION x(3),B2(3),xsm(3),Bsm(3),xr(3),Br(3),sm2gsm(3,3)
     *,zu(3,3),Bsm1(3),b1(3),Br1(3)
      do i=1,3
      br(i)=0.
      b1(i)=0.
      b2(i)=0.
      end do
      tetam=tm1
      tm0=tm2
      aj0=ami1*1.e6
      ami=ami2
      p=pi/180. 
      call SMtoGSM(SM2GSM)
      call PERE2(x,xsm,sm2gsm,-1)              
      r=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      cte=xsm(3)/r
      ste=sqrt(1-cte*cte)
       if (ste.eq.0.) then 
        cfie=1.
        sfie=0.
       else
        cfie=xsm(1)/r/ste
        sfie=xsm(2)/r/ste
       end if
      xr(1)=r
      xr(3)=acos(cfie)/p
      if (xsm(2).lt.0.) xr(3)=360-xr(3)
      xr(2)=acos(cte)/p
      ZU(3,1)=Cte                                                        
      ZU(3,2)=-STe                                                      
      ZU(3,3)=0.                                                        
      ZU(2,1)=STe*SFIe                                                    
      ZU(2,2)=CTe*SFIe                                                    
      ZU(2,3)=CFIe                                                       
      ZU(1,1)=STe*CFIe                                                    
      ZU(1,2)=CTe*CFIe                                                    
      ZU(1,3)=-SFIe  
      ss2=0.             ! FAC2 Currents are cancelled!                                                    
      if (ss2.eq.0.) goto1
c        call FAC2d(ami,tm0,xr,br)     
      call PERE2(br,bsm,zu,1)              
      call PERE2(bsm,b2,sm2gsm,1)
1     if (ss1.eq.0.) return
        call bFAC(br1)     
      call PERE2(br1,bsm1,zu,1)              
      call PERE2(bsm1,b1,sm2gsm,1)

      return          
      end


*********************************************************
*     AUXILLIARY SUBROUTINES
*********************************************************
C
      SUBROUTINE PERE2(A,B,T,K)
**********************************************************
*     Transition A into B vectors by T (K>0)
*                       B=T*A      
*     or T^{-1} matrices (K<=0)
*                       B=T^{-1}*A      
* Written by V. Kalegaev         					 *     
**********************************************************
      DIMENSION A(3),B(3),T(3,3)
      if (k) 1,1,2
2     b(1)=t(1,1)*a(1)+t(1,2)*a(2)+t(1,3)*a(3)
      b(2)=t(2,1)*a(1)+t(2,2)*a(2)+t(2,3)*a(3)
      b(3)=t(3,1)*a(1)+t(3,2)*a(2)+t(3,3)*a(3)
      return
1     b(1)=t(1,1)*a(1)+t(2,1)*a(2)+t(3,1)*a(3)
      b(2)=t(1,2)*a(1)+t(2,2)*a(2)+t(3,2)*a(3)
      b(3)=t(1,3)*a(1)+t(2,3)*a(2)+t(3,3)*a(3)
      return
      end
C
      SUBROUTINE SMtoGSM(SM2GSM)
****************************************************************
*     Calculation of the transition matrix from SM coordinates to 
*     GSM ones:
*                  VectGSM=(SM2GSM)*VectSM
* Written by V. Kalegaev         					 *     
****************************************************************
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      DIMENSION SM2GSM(3,3)
      SM2GSM(1,1)=CPSI
      SM2GSM(1,2)=0.
      SM2GSM(1,3)=-SPSI
      SM2GSM(2,1)=0.
      SM2GSM(2,2)=1.
      SM2GSM(2,3)=0.
      SM2GSM(3,1)=SPSI
      SM2GSM(3,2)=0.
      SM2GSM(3,3)=CPSI
      return
      END
C
      SUBROUTINE PSTATUS(X1,X2,X3,X4,X5,X6,X7)
*****************************************************************
*     Determination of the parameters providing the model tuning
* Written by V. Kalegaev         					 *     
*****************************************************************
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      SSD =x1 ! dipole field on/off (1/0)
      SSR =x2 ! RC field on/off (1/0)
      SBT =x3 ! tail current field on/off (1/0)
      SMD =x4 ! dipole shielding field on/off (1/0)
      SMR =x5 ! RC shielding field IMF on/off (1/0)
      SS1 =x6 ! Region 1 FAC field on/off (1/0)
      SIMF=x7 ! IMF on/off (1/0)
      SS2 =0. ! Region 2 FAC field off (to be zero yet!)
      return
      END
C
      SUBROUTINE SUN_a2000(IYR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  THIS SUBROUTINE HAS BEEN COMPILED FROM: RUSSELL C.T., COSM.ELECTRO-
C  DYN., 1971, V.2,PP.184-196.
C
C
C                   AUTHOR: Gilbert D. Mead
C
C
        IMPLICIT NONE

        REAL GST,SLONG,SRASN,SDEC,RAD,T,VL,G,OBLIQ,SOB,SLP,SIND,
     1       COSD,SC

        INTEGER*4 IYR,IDAY,IHOUR,MIN,ISEC
      common /ddd/ sind,cosd
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
      IF(IYR.LT.1901.OR.IYR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYR-1900)+(IYR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
    

