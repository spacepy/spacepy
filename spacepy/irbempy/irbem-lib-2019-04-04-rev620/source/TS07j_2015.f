C This is a set of replacement routines written by Jay Albert, that improve 
C execution time by reducing calls to bessj
C
C Tested by A. C. Kellerman against the original, output was identical
C for 10,000 points R = 1:10 across 24 MLT sectors
C=======================================================================
        SUBROUTINE  SHTBNORM_E_2015(K,L,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_E
        IMPLICIT  REAL * 8  (A - H, O - Z)
        DIMENSION AK(5)
        DIMENSION AJM(0:14),AJMD(0:14)
        COMMON /TSE/ TSE(80,5,4)

        AK(1)=TSE(76,K,L)
        AK(2)=TSE(77,K,L)
        AK(3)=TSE(78,K,L)
        AK(4)=TSE(79,K,L)
        AK(5)=TSE(80,K,L)

        phi=DATAN2(Y,X)
        RHO=dsqrt(X*X+Y*Y)
        if(RHO.lt.1.D-8) then
           RHOI=1.D8
        else
           RHOI=1.D0/RHO
        endif
        DPDX=-Y*RHOI*RHOI
        DPDY= X*RHOI*RHOI

        FX=0.D0
        FY=0.D0
        FZ=0.D0

        DO 2 n=1,5
           AKN=dabs(AK(n))
           AKNR=AKN*RHO
           if(AKNR.lt.1.D-8) then
              AKNRI=1.D8
           else
              AKNRI=1.D0/AKNR
           endif

           CHZ=dcosh(Z*AKN)
           SHZ=dsinh(Z*AKN)

           AJM(0)=BESSJJ_2015(14,AKNR, AJM(1)) !!! get all n in one call
           DO 3 m=1,14
              AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
  3        CONTINUE
           AJMD(0)=-AJM(1)

           DO 4 m=0,14
              CMP=DCOS(m*phi)
              SMP=DSIN(m*phi)

              HX1=-m*DPDX*CMP*SHZ*AJM(m)
              HX2=-AKN*X*RHOI*SMP*SHZ*AJMD(m)
              HX=HX1+HX2
              HY1=-m*DPDY*CMP*SHZ*AJM(m)
              HY2=-AKN*Y*RHOI*SMP*SHZ*AJMD(m)
              HY=HY1+HY2
              HZ=-AKN*SMP*CHZ*AJM(m)

              L1=n+5*m
              FX=FX+HX*TSE(L1,K,L)
              FY=FY+HY*TSE(L1,K,L)
              FZ=FZ+HZ*TSE(L1,K,L)
  4        CONTINUE
  2     CONTINUE
      RETURN
      END
C==---------------------------------------------------------------------
        SUBROUTINE  SHTBNORM_O_2015(K,L,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_O
        IMPLICIT  REAL * 8  (A - H, O - Z)
        DIMENSION AK(5)
        DIMENSION AJM(0:14),AJMD(0:14)
        COMMON /TSO/ TSO(80,5,4)

        AK(1)=TSO(76,K,L)
        AK(2)=TSO(77,K,L)
        AK(3)=TSO(78,K,L)
        AK(4)=TSO(79,K,L)
        AK(5)=TSO(80,K,L)

        phi=DATAN2(Y,X)
        RHO=dsqrt(X*X+Y*Y)
        if(RHO.lt.1.D-8) then
           RHOI=1.D8
        else
           RHOI=1.D0/RHO
        endif
        DPDX=-Y*RHOI*RHOI
        DPDY= X*RHOI*RHOI

        FX=0.D0
        FY=0.D0
        FZ=0.D0

        DO 2 n=1,5
           AKN=dabs(AK(n))
           AKNR=AKN*RHO
           if(AKNR.lt.1.D-8) then
              AKNRI=1.D8
           else
              AKNRI=1.D0/AKNR
           endif

           CHZ=dcosh(Z*AKN)
           SHZ=dsinh(Z*AKN)

           AJM(0)=BESSJJ_2015(14,AKNR, AJM(1)) !!! get all n in one call
           DO 3 m=1,14
              AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
  3        CONTINUE
           AJMD(0)=-AJM(1)

           DO 4 m=0,14
              CMP=DCOS(m*phi)
              SMP=DSIN(m*phi)

              HX1=m*DPDX*SMP*SHZ*AJM(m)
              HX2=-AKN*X*RHOI*CMP*SHZ*AJMD(m)
              HX=HX1+HX2
              HY1=m*DPDY*SMP*SHZ*AJM(m)
              HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD(m)
              HY=HY1+HY2
              HZ=-AKN*CMP*CHZ*AJM(m)

              L1=n+5*m
              FX=FX+HX*TSO(L1,K,L)
              FY=FY+HY*TSO(L1,K,L)
              FZ=FZ+HZ*TSO(L1,K,L)
  4        CONTINUE
  2     CONTINUE
        RETURN
        END
C==---------------------------------------------------------------------
        SUBROUTINE  SHTBNORM_S_2015(K,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_S
        IMPLICIT  REAL * 8  (A - H, O - Z)
        DIMENSION AK(5)
        DIMENSION AJM(0:14),AJMD(0:14)
        COMMON /TSS/ TSS(80,5)

        AK(1)=TSS(76,K)
        AK(2)=TSS(77,K)
        AK(3)=TSS(78,K)
        AK(4)=TSS(79,K)
        AK(5)=TSS(80,K)

        phi=DATAN2(Y,X)
        RHO=dsqrt(X*X+Y*Y)
        if (RHO.lt.1.D-8) then
           RHOI=1.D8
        else
           RHOI=1.D0/RHO
        endif
        DPDX=-Y*RHOI*RHOI
        DPDY= X*RHOI*RHOI

        FX=0.D0
        FY=0.D0
        FZ=0.D0

        DO 2 n=1,5
           AKN=dabs(AK(n))
           AKNR=AKN*RHO
           if (AKNR.lt.1.D-8) then
              AKNRI=1.D8
           else
              AKNRI=1.D0/AKNR
           endif

           CHZ=dcosh(Z*AKN)
           SHZ=dsinh(Z*AKN)

           AJM(0)=BESSJJ_2015(14,AKNR, AJM(1)) !!! get all n in one call
           DO 3 m=1,14
              AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
    3      CONTINUE
           AJMD(0)=-AJM(1)

           DO 4 m=0,14
              CMP=DCOS(m*phi)
              SMP=DSIN(m*phi)
   
              HX1=m*DPDX*SMP*SHZ*AJM(m)
              HX2=-AKN*X*RHOI*CMP*SHZ*AJMD(m)
              HX=HX1+HX2
              HY1=m*DPDY*SMP*SHZ*AJM(m)
              HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD(m)
              HY=HY1+HY2
              HZ=-AKN*CMP*CHZ*AJM(m)
   
              L=n+5*m
              FX=FX+HX*TSS(L,K)
              FY=FY+HY*TSS(L,K)
              FZ=FZ+HZ*TSS(L,K)
  4        CONTINUE
  2     CONTINUE
        RETURN
        END
C==---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BESSJJ_2015(n,x, bessJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension bessJ(n) ! BESSJJ_2015 returns J0, bessJ holds J1 to Jn
      PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
CU    USES bessj0,bessj1
      if(n.lt.2) then
         print *,' *** bad argument n in BESSJJ_2015'
         stop
      endif
      ax=Dabs(x)
      if(ax.eq.0.D0)then
        BESSJJ_2015=0.
      else if(ax.gt.Dfloat(n))then
        tox=2.D0/ax
        bjm=bessj0(ax)
        bj =bessj1(ax)
        BESSJJ_2015  =bjm ! J0(x)
        bessJ(1)=bj  ! J1(x)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
          bessJ(j+1)=bjp ! J_j+1(x)
11      continue
      else
         tox=2.D0/ax
         m=2*((n+int(Dsqrt(Dfloat(IACC*n))))/2)
         jsum=0
         sum=0.D0
         bjp=0.D0
         bj =1.D0
         do i=1,n
            bessJ(i)=0.
         enddo
   
         do 12 j=m,1,-1
            bjm=j*tox*bj-bjp
            bjp=bj
            bj=bjm
   
            if (Dabs(bj).gt.BIGNO) then
               bj =bj *BIGNI
               bjp=bjp*BIGNI
               sum=sum*BIGNI
   	    do i=j+1,n
   	       bessJ(i)=bessJ(i)*BIGNI
               enddo
            endif
   
            if(jsum.ne.0)sum=sum+bj
            jsum=1-jsum
            if(j.le.n) bessJ(j)=bjp ! Jj(x)
   12    continue
   
         sum=2.D0*sum-bj
         do i=1,n
            bessJ(i)=bessJ(i)/sum
         enddo
         BESSJJ_2015=bj/sum ! J0(x)
      endif

      if (x .lt. 0.D0) then
         do i=1,n,2
            bessJ(i)=-bessJ(i)
         enddo
      endif
      return
      END
C=======================================================================
