       PROGRAM MUON00
        DIMENSION I9(4),K9(4,4),M9(4,7),J9(9,10),J4(9,5,4),JX(4),ID9(5),
     1  IXX(3,20),PLX(20)                                               
        DOUBLE PRECISION F,DZA,DZA2,DREDM                               
        COMMON/LOC001/IJK,ENERG,ECONS,ECONST,D2P1SM,D2P1S               
        COMMON/LOC002/BEM(3),ZSA(3),BE(3)                               
        COMMON/LOC003/K0,K1,K2,K3                                       
        COMMON/LOC004/NN0(3),NN1(7),NN2(7),NN3(7)                       
        COMMON/LOC005/R0(3),R1(7),R2(7),R3(7)                           
        COMMON/LOC006/IP1(7),IP2(7),IP3(7),IQ1(7),IQ2(7),IQ3(7)         
        COMMON/LOC007/IC                                                
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC009/F(60),FD                                          
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC013/COEMON(30),EXPMON(30)                             
        COMMON/LOC014/COEDIP(42),EXPDIP(42)                             
        COMMON/LOC015/COEQUA(45),EXPQUA(45)                             
        COMMON/LOC016/COEOCT(45),EXPOCT(45)                             
        COMMON/LOC017/IFM(6)                                            
        COMMON/LOC018/IFD(9)                                            
        COMMON/LOC019/IFQ(10)                                           
        COMMON/LOC020/IFO(10)                                           
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC023/COEDP(9)                                          
        COMMON/LOC024/COEQ(11)                                          
        COMMON/LOC025/COEO(11)                                          
        COMMON/LOC026/COED(4)                                           
        COMMON/LOC027/LL(20)                                            
        COMMON/LOC028/M1(7),M2(7),M3(7),YC(4),IDB                       
        COMMON/LOC029/IRR,RR(18),RAU,RAD,RA(4),RD(4),RSA(4)             
        COMMON/LOC030/POP(6),JTM(6),JTD(6),JTQ(6),JTO(6)                
        COMMON/LOC031/JM(10),JD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC033/M,E(1000),AI(1000),IA(1000),ENERGY(20,40)         
        COMMON/LOC034/DZA,DZA2,DREDM                                    
        COMMON/LOC035/ICC,CD(5),EA,EB,IDR                               
        COMMON/LOC036/ZMK,ZML,ZMM,ZMKM,ZMLM,ZMMM,IVERS                  
        COMMON/LOC037/PL(20),NPOL(20),IPOL,CL1,CL2,IDE,PLN(210),IP8     
        COMMON/LOC038/A,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)               
        COMMON/LOC039/NOPT,NMAX,ALEXP                                   
        COMMON/LOC040/AMASSA,AMASSN,HBAR                                
        COMMON/LOC041/MPU,ICPU(200),IPN                                 
C  ***  NOTE -- NON-STANDARD WAY OF OVERPRINTING -- CHANGE AS NEEDED    
C  ***  STANDARD FORTRAN REQUIRES DATA FOR I9 /1H ,1H+,1H+,1H+/         
C       Hollerith constants allow putting a character into a string
C       i.e 2HRD = RD, 2H3D = 3D etc
        DATA K9/2HRD,2H1S,2H2T,2H3T,2HRD,2H1S,2H2S,2H3S,2HRD,2H**,2H2P, 
     1  2H3P,2HRD,2H**,2H**,2H3D/                                       
        DATA I9,L9/1H+,1H+,1H+,1H ,2H--/                                
        DATA J9/124,254,455,387,387,387,455,254,124,024,056,120,024,024,
     1  024,024,126,126,124,254,455,391,030,056,112,255,511,511,510,028,
     2  112,056,014,391,254,124,014,030,054,102,255,511,006,015,015,511,
     3  511,384,508,254,007,391,254,124,124,254,385,508,510,387,455,254,
     4  124,511,511,007,014,028,056,112,224,448,124,254,387,455,254,455,
     5  387,254,124,124,254,455,387,255,127,259,254,124/                
        DATA JX,JB/1HX,1HI,1HH,1HO,1H /                                 
C
     0  DO 1100 I=1,20                                                  
        PL(I) = 0.000                                                   
 1100   CONTINUE                                                        
        M = 1                                                           
C
C       This is a 20x40 matrix. This is bigger than we need, as
C       we never go above n = 20 for the start of the cascade
        DO 1150 I=1,20                                                  
        DO 1150 J=1,40                                                  
        ENERGY(I,J) = 0.000                                             
 1150   CONTINUE                                                        
 1175   M = 1                                                           
C      
C       Read in the input file
C       When XEQ is encountered, we jump back out of here
C       ready to compute the cascade
C       The program will terminate here if STO(P) is encountered
        CALL RREAD                                                      
C
C       This sets some parameters and we are now ready to do some calculating
        CALL FFIX                                                       
C
C       Loops over the maximum N and calculates the number of subshells in the
C       energy level. This corresponds to NMX in input
        DO 1500 N=1,NMAX                                                
C         K is number of subshells in N
          K = 2*N-1                                                       
C
C         Now loop over each of the subshells
          DO 1500 L=1,K                                                   
C
C           This calculates the L value for the corresponding subshell
            L1 = (L-1)/2 + 1                                                
            NJ = N                                                          
C
C           Now calculate the Dirac point energy for (N,L)
            PT = POINT(NJ,L1)                                               
            IF(ENERGY(N,L).LE.0.000.AND.MOD(IPRINT,2).EQ.0)                 
C
     1  WRITE(IW,1200)N,L,PT                                            
C
C
 1200   FORMAT(27H NO INPUT DATA FOR STATE N=,I2,7H, (LJ)=,I2,          
     1  29H,  POINT-LIKE DIRAC ENERGY = ,F10.6,5H(MEV))                 
C       This condition only occurs if we have inputted energies in
     0  IF(ENERGY(N,L).LE.0.000)  ENERGY(N,L) = PT-ENERGY(N,L)          
 1500   CONTINUE                                                        

C       ! If nmax is 20 then we are OK as PL has been allocated to the right
C       size
        IF(NMAX.EQ.20)  GO TO 1700                                      
C
C       Make sure there are zeros in the rest of the population array
        MX1 = NMAX+1                                                    
        DO 1600 I=MX1,20                                                
        PL(I) = 0.000                                                   
 1600   CONTINUE                                                        

C       NOPT defaults to 0
C       0 is mod stat, 2 is quadratic
 1700   IF(NOPT.NE.0)  GO TO 1900                                       
        DO 1800 I=1,NMAX                                                
C
C       Here is the modified stat dist for the top energy level
        PL(I) = FLOAT(2*I-1)*EXP( ALEXP*FLOAT(I-1))                     
 1800   CONTINUE                                                        
C   
C       Calculate quadratic distribution
 1900   IF(NOPT.NE.2)  GO TO 1930                                       
        DO 1925 I=1,NMAX                                                
        PL(I) = 1.000 + CL1*FLOAT(I-1) + CL2*FLOAT((I-1)**2)            
 1925   CONTINUE                                                        

C       The following distributions are implemented as in Hartmann
C       PDJ Modification here to add in a statistical distribution
C       of 2l + 1
 1930   IF(NOPT.NE.4) GO TO 1940
        DO 1935 I=1, NMAX
        PL(I) = FLOAT(2*I-1)
 1935   CONTINUE

C 
C       PDJ purely linear distribution 1 + alpha*l
 1940   IF(NOPT.NE.8) GO TO 1950
        DO 1945 I=1, NMAX
        PL(I) = 1.0 + ALEXP*FLOAT(I-1)
 1945   CONTINUE

C       PDJ log distribution shifted to be 1 at l = 0
 1950   IF(NOPT.NE.16) GO TO 1960
        DO 1955 I=1, NMAX
        PL(I) = CL1*LOG(CL2*(FLOAT(I-1) + EXP(1.0/CL1)/CL2))
 1955   CONTINUE
C
C       SS is the full population of the top energy level
C       and is used to normalise
 1960   SS = 0.000                                                      
        DO 2000 I=1,NMAX                                                
        SS = SS + PL(I)                                                 
 2000   CONTINUE                                                        

 319    FORMAT(ES17.4)

C       Catch negative population
        IF(SS.LE.0.000)  WRITE(IW,2100)(PL(I),I=1,NMAX)                 
 2100   FORMAT(47H *** ERROR *** INITIAL L DISTRIBUTION WRONG ***/      
     1  (1X,5F10.6))                                                    
C
C       If the population is negative, then carry on
     0  IF(SS.LE.0.000)  SS=1.000                                       
C
C       Now we normalise by the population to make sure that there is one muon
        DO 2200 I=1,NMAX                                                
        PL(I) = PL(I)/SS                                                
 2200   CONTINUE                                                        
C
C
        IF(MPU.GT.1.AND.IPN.EQ.0)  WRITE(IPUNCH,2300)MPU,IDE            
 2300   FORMAT(14(1H*),19H NEW FIT -- AT MOST,I4,15H LINES PUNCHED ,    
     1  20(1H*),I5)                                                     
C
C       Idk what this is doing, L9 is some big random integer
        DO 2400 I=1,4                                                   
        DO 2400 J=1,7                                                   
        M9(I,J) = L9                                                    
 2400   CONTINUE                                                        
C
C       K0 starts off at 3. So, on the first loop round, we need to do the
C       following calculation
        IF(K0.EQ.0)  GO TO 2500                                         
C
C       Loop over 1 to 3
        DO 2450 I=1,K0                                                  
C       write(4, 69)NN0(I)
        J = NN0(I)                                                      
        M9(1,I) = K9(J+1,1)                                             
 2450   CONTINUE                                                        
69      FORMAT(I9)
C
C
 2500   IF(K1.EQ.0)  GO TO 2600                                         
        DO 2550 I=1,K1                                                  
        J = NN1(I)                                                      
        K = M1(I)                                                       
        M9(2,I) = K9(J+1,K+1)                                           
 2550   CONTINUE                                                        
C
C
 2600   IF(K2.EQ.0)  GO TO 2700                                         
        DO 2650 I=1,K2                                                  
        J = NN2(I)                                                      
        K = M2(I)                                                       
        M9(3,I) = K9(J+1,K+1)                                           
 2650   CONTINUE                                                        
C
C
 2700   IF(K3.EQ.0)  GO TO 2800                                         
        DO 2750 I=1,K3                                                  
        J = NN3(I)                                                      
        K = M3(I)                                                       
        M9(4,I) = K9(J+1,K+1)                                           
 2750   CONTINUE                                                        
 2800   IZ9 = IFIX(Z+0.500)                                             
        IZ1 = MOD(IZ9,10) + 1                                           
        IZ2 = IZ9/10 + 1                                                
        ID1 = MOD(IDE,10)+1                                             
        ID2 = MOD(IDE/10,10)+1                                          
        ID5 = IDE/10000+1                                               
        ID3 = MOD(IDE/100,10) + 1                                       
        ID4 = MOD(IDE/1000,10) + 1                                      
C
        IF(IP8.EQ.0)  GO TO 3100                                        
C
        SS = 0.000                                                      
C       NU tells us the total number of states that we have
        NU = NMAX*(NMAX+1)/2                                            
        DO 2900 I=1,NU                                                  
        SS = SS + PLN(I)                                                
 2900   CONTINUE                                                        
        IF(SS.LE.0.000)  WRITE(IW,2100)                                 
        IF(SS.LE.0.000)  SS=1.000                                       
        DO 3000 I=1,NU                                                  
        PLN(I) = PLN(I)/SS                                              
 3000   CONTINUE                                                        
 3100   ID9(1) = IZ2                                                    
        ID9(2) = IZ1                                                    
        WRITE(IW,4000)I9(1)                                             
        DO 4099 I=2,4                                                   
        WRITE(IW,4010)I9(I)                                             
 4099   CONTINUE                                                        
        DO 4199 I=1,4                                                   
        WRITE(IW,4100)I9(I)                                             
 4199   CONTINUE                                                        
        DO 4249 I=1,9                                                   
        DO 4249 J=1,2                                                   
        DO 4249 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(1,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4249   CONTINUE                                                        
        WRITE(IW,4200)I9(1),((J4(I,J,1),I=1,9),J=1,2)                   
        DO 4299 K=2,4                                                   
        WRITE(IW,4210)I9(K),((J4(I,J,K),I=1,9),J=1,2)                   
 4299   CONTINUE                                                        
        DO 4349 I=1,9                                                   
        DO 4349 J=1,2                                                   
        DO 4349 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(2,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4349   CONTINUE                                                        
        WRITE(IW,4300)I9(1),((J4(I,J,1),I=1,9),J=1,2)                   
        DO 4399 K=2,4                                                   
        WRITE(IW,4310)I9(K),((J4(I,J,K),I=1,9),J=1,2)                   
 4399   CONTINUE                                                        
        DO 4449 I=1,9                                                   
        DO 4449 J=1,2                                                   
        DO 4449 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(3,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4449   CONTINUE                                                        
        WRITE(IW,4400)I9(1),((J4(I,J,1),I=1,9),J=1,2)                   
        DO 4499 K=2,4                                                   
        WRITE(IW,4410)I9(K),((J4(I,J,K),I=1,9),J=1,2)                   
 4499   CONTINUE                                                        
        DO 4549 I=1,9                                                   
        DO 4549 J=1,2                                                   
        DO 4549 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(4,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4549   CONTINUE                                                        
        WRITE(IW,4500)I9(1),((J4(I,J,1),I=1,9),J=1,2),NMAX,NOPT         
        DO 4599 K=2,4                                                   
        WRITE(IW,4510)I9(K),((J4(I,J,K),I=1,9),J=1,2),NMAX,NOPT         
 4599   CONTINUE                                                        
        DO 4649 I=1,9                                                   
        DO 4649 J=1,2                                                   
        DO 4649 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(5,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4649   CONTINUE                                                        
        WRITE(IW,4600)I9(1),((J4(I,J,1),I=1,9),J=1,2),ALEXP,CL1,CL2     
        DO 4699 K=2,4                                                   
        WRITE(IW,4610)I9(K),((J4(I,J,K),I=1,9),J=1,2),ALEXP,CL1,CL2     
 4699   CONTINUE                                                        
        DO 4749 I=1,9                                                   
        DO 4749 J=1,2                                                   
        DO 4749 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(6,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4749   CONTINUE                                                        
        IF(IP8.EQ.0)  WRITE(IW,4700)I9(1),((J4(I,J,1),I=1,9),           
     1  J=1,2),(PL(I),I=1,10)                                           
        IF(IP8.NE.0)  WRITE(IW,4720)I9(1),((J4(I,J,1),I=1,9),J=1,2)     
        DO 4799 K=2,4                                                   
        IF(IP8.EQ.0)  WRITE(IW,4710)I9(K),((J4(I,J,K),                  
     1  I=1,9),J=1,2),(PL(I),I=1,10)                                    
        IF(IP8.NE.0)  WRITE(IW,4730)I9(K),((J4(I,J,K),I=1,9),J=1,2)     
 4799   CONTINUE                                                        
        DO 4849 I=1,9                                                   
        DO 4849 J=1,2                                                   
        DO 4849 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(7,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4849   CONTINUE                                                        
        IF(IP8.EQ.0)  WRITE(IW,4800)I9(1),((J4(I,J,1),I=1,9),           
     1  J=1,2),(PL(I),I=11,20)                                          
        IF(IP8.NE.0)  WRITE(IW,4820)I9(1),((J4(I,J,1),I=1,9),J=1,2)     
        DO 4899 K=2,4                                                   
        IF(IP8.EQ.0)  WRITE(IW,4810)I9(K),((J4(I,J,K),                  
     1  I=1,9),J=1,2),(PL(I),I=11,20)                                   
        IF(IP8.NE.0)  WRITE(IW,4830)I9(K),((J4(I,J,K),I=1,9),J=1,2)     
 4899   CONTINUE                                                        
        DO 4949 I=1,9                                                   
        DO 4949 J=1,2                                                   
        DO 4949 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(8,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 4949   CONTINUE                                                        
        WRITE(IW,4900)I9(1),((J4(I,J,1),I=1,9),J=1,2)                   
        DO 4999 K=2,4                                                   
        WRITE(IW,4910)I9(K),((J4(I,J,K),I=1,9),J=1,2)                   
 4999   CONTINUE                                                        
        DO 5049 I=1,9                                                   
        DO 5049 J=1,2                                                   
        DO 5049 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(9,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 5049   CONTINUE                                                        
        WRITE(IW,5000)I9(1),((J4(I,J,1),I=1,9),J=1,2)                   
        DO 5099 K=2,4                                                   
        WRITE(IW,5010)I9(K),((J4(I,J,K),I=1,9),J=1,2)                   
 5099   CONTINUE                                                        
        WRITE(IW,5100)I9(1)                                             
        DO 5199 K=2,4                                                   
        WRITE(IW,5110)I9(K)                                             
 5199   CONTINUE                                                        
C       WRITE(IW,5200)I9(1),ZSK,ZSL,ZSM,(POP(I),I=1,3)                  
        DO 5299 K=2,4                                                   
        WRITE(IW,5210)I9(K),ZSK,ZSL,ZSM,(POP(I),I=1,3)                  
 5299   CONTINUE                                                        
        WRITE(IW,5300)I9(1),(POP(I),I=4,6)                              
        DO 5399 K=2,4                                                   
        WRITE(IW,5310)I9(K),(POP(I),I=4,6)                              
 5399   CONTINUE                                                        
        WRITE(IW,5400)I9(1)                                             
        DO 5499 K=2,4                                                   
        WRITE(IW,5410)I9(K)                                             
 5499   CONTINUE                                                        
        WRITE(IW,5500)I9(1),IPC,WIDTHK                                  
        DO 5599 K=2,4                                                   
        WRITE(IW,5510)I9(K),IPC,WIDTHK                                  
 5599   CONTINUE                                                        
        WRITE(IW,5600)I9(1),BE                                          
        DO 5699 K=2,4                                                   
        WRITE(IW,5610)I9(K),BE                                          
 5699   CONTINUE                                                        
        WRITE(IW,5700)I9(1)                                             
        DO 5799 K=2,4                                                   
        WRITE(IW,5710)I9(K)                                             
 5799   CONTINUE                                                        
        A9 = A                                                          
        IF(AMASSA.NE.0.000)  A9 = AMASSA                                
        TQM = ABS(TFM-2.3001)                                           
        IF(TQM.GE.1.0E-10)  WRITE(IW,5800)I9(1),A9,CFM,TFM              
        IF(TQM.LT.1.0E-10)  WRITE(IW,5820)I9(1),A9                      
        DO 5899 K=2,4                                                   
        IF(TQM.GE.1.0E-10)  WRITE(IW,5810)I9(K),A9,CFM,TFM              
        IF(TQM.LT.1.0E-10)  WRITE(IW,5830)I9(K),A9                      
 5899   CONTINUE                                                        
        IF(STEP+RMATCH.GT.1.0E-20)  WRITE(IW,5900)I9(1),                
     1  STEP,RMATCH                                                     
        IF(STEP+RMATCH.LE.1.0E-20)  WRITE(IW,5920)I9(1)                 
        DO 5999 K=2,4                                                   
        IF(STEP+RMATCH.GT.1.0E-20)  WRITE(IW,5910)I9(K),                
     1  STEP,RMATCH                                                     
        IF(STEP+RMATCH.LE.1.0E-20)  WRITE(IW,5930)I9(K)                 
 5999   CONTINUE                                                        
        WRITE(IW,6000)I9(1),AMASSM,AMASSE,AMASSN                        
        DO 6099 K=2,4                                                   
        WRITE(IW,6010)I9(K),AMASSM,AMASSE,AMASSN                        
 6099   CONTINUE                                                        
        WRITE(IW,6100)I9(1),D2P1S,ESP                                   
        DO 6199 K=2,4                                                   
        WRITE(IW,6110)I9(K),D2P1S,ESP                                   
 6199   CONTINUE                                                        
        WRITE(IW,6200)I9(1)                                             
        DO 6299 K=2,4                                                   
        WRITE(IW,6210)I9(K)                                             
 6299   CONTINUE                                                        
        WRITE(IW,6300)I9(1),EHIGH,CLIMIT                                
        DO 6399 K=2,4                                                   
        WRITE(IW,6310)I9(K),EHIGH,CLIMIT                                
 6399   CONTINUE                                                        
        WRITE(IW,6400)I9(1),ELOW,ICC                                    
        DO 6499 K=2,4                                                   
        WRITE(IW,6410)I9(K),ELOW,ICC                                    
 6499   CONTINUE                                                        
        EAB = (EA-99.000)**2 + (EB-99.000)**2                           
        IF(EAB.GT.1.0E-20)  WRITE(IW,6500)I9(1),ERES,EA,EB              
        IF(EAB.LE.1.0E-20)  WRITE(IW,6520)I9(1),ERES                    
        DO 6599 K=2,4                                                   
        IF(EAB.GT.1.0E-20)  WRITE(IW,6510)I9(K),ERES,EA,EB              
        IF(EAB.LE.1.0E-20)  WRITE(IW,6530)I9(K),ERES                    
 6599   CONTINUE                                                        
        WRITE(IW,6600)I9(1),CD                                          
        DO 6699 K=2,4                                                   
        WRITE(IW,6610)I9(K),CD                                          
 6699   CONTINUE                                                        
        WRITE(IW,6700)I9(1),NPOL,IPOL                                   
        DO 6799 K=2,4                                                   
        WRITE(IW,6710)I9(K),NPOL,IPOL                                   
 6799   CONTINUE                                                        
        DO 6899 K=1,4                                                   
        WRITE(IW,6800)I9(K)                                             
 6899   CONTINUE                                                        
        DO 6999 K=1,4                                                   
        WRITE(IW,6900)I9(K)                                             
 6999   CONTINUE                                                        
        DO 7099 K=1,4                                                   
        WRITE(IW,7000)I9(K)                                             
 7099   CONTINUE                                                        
        WRITE(IW,7100)I9(1)                                             
        DO 7199 K=2,4                                                   
        WRITE(IW,7110)I9(K)                                             
 7199   CONTINUE                                                        
        WRITE(IW,7200)I9(1)                                             
        DO 7299 K=2,4                                                   
        WRITE(IW,7210)I9(K)                                             
 7299   CONTINUE                                                        
        WRITE(IW,7300)I9(1),K0,K1,K2,K3,IREAD,IW,IPUNCH                 
        DO 7399 K=2,4                                                   
        WRITE(IW,7310)I9(K),K0,K1,K2,K3,IREAD,IW,IPUNCH                 
 7399   CONTINUE                                                        
        WRITE(IW,7400)I9(1)                                             
        DO 7499 K=2,4                                                   
        WRITE(IW,7410)I9(K)                                             
 7499   CONTINUE                                                        
        DO 7599 K=1,4                                                   
        WRITE(IW,7500)I9(K)                                             
 7599   CONTINUE                                                        
        WRITE(IW,7600)I9(1),IPRINT                                      
        DO 7699 K=2,4                                                   
        WRITE(IW,7610)I9(K),IPRINT                                      
 7699   CONTINUE                                                        
        WRITE(IW,7700)I9(1),(M9(1,J),J=1,3),IDB                         
        DO 7799 K=2,4                                                   
        WRITE(IW,7710)I9(K),(M9(1,J),J=1,3),IDB                         
 7799   CONTINUE                                                        
        WRITE(IW,7800)I9(1),(M9(2,J),J=1,7),IC                          
        DO 7899 K=2,4                                                   
        WRITE(IW,7810)I9(K),(M9(2,J),J=1,7),IC                          
 7899   CONTINUE                                                        
        WRITE(IW,7900)I9(1),(M9(3,J),J=1,7)                             
        DO 7999 K=2,4                                                   
        WRITE(IW,7910)I9(K),(M9(3,J),J=1,7)                             
 7999   CONTINUE                                                        
        WRITE(IW,8000)I9(1),(M9(4,J),J=1,7)                             
        DO 8099 K=2,4                                                   
        WRITE(IW,8010)I9(K),(M9(4,J),J=1,7)                             
 8099   CONTINUE                                                        
        WRITE(IW,8100)I9(1),FD                                          
        DO 8199 K=2,4                                                   
        WRITE(IW,8110)I9(K),FD                                          
 8199   CONTINUE                                                        
        MDIR = NMAX**2                                                  
        WRITE(IW,8200)I9(1),IDR,MDIR                                    
        DO 8299 K=2,4                                                   
        WRITE(IW,8210)I9(K),IDR,MDIR                                    
 8299   CONTINUE                                                        
        WRITE(IW,8300)I9(1),MPU                                         
        DO 8399 K=2,4                                                   
        WRITE(IW,8310)I9(K),MPU                                         
 8399   CONTINUE                                                        
        WRITE(IW,8400)I9(1),IPN                                         
        DO 8499 K=2,4                                                   
        WRITE(IW,8410)I9(K),IPN                                         
 8499   CONTINUE                                                        
        WRITE(IW,8500)I9(1)                                             
        DO 8599 K=2,4                                                   
        WRITE(IW,8510)I9(K)                                             
 8599   CONTINUE                                                        
        DO 8699 K=1,4                                                   
        WRITE(IW,8600)I9(K)                                             
 8699   CONTINUE                                                        
        WRITE(IW,8700)I9(1)                                             
        DO 8799 K=2,4                                                   
        WRITE(IW,8710)I9(K)                                             
 8799   CONTINUE                                                        
        WRITE(IW,8800)I9(1)                                             
        DO 8899 K=2,4                                                   
        WRITE(IW,8810)I9(K)                                             
 8899   CONTINUE                                                        
        ID9(1) = ID5                                                    
        ID9(2) = ID4                                                    
        ID9(3) = ID3                                                    
        ID9(4) = ID2                                                    
        ID9(5) = ID1                                                    
        DO 8949 I=1,9                                                   
        DO 8949 J=1,5                                                   
        DO 8949 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(1,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 8949   CONTINUE                                                        
        WRITE(IW,8900)I9(1),((J4(I,J,1),I=1,9),J=1,5)                   
        DO 8999 K=2,4                                                   
        WRITE(IW,8910)I9(K),((J4(I,J,K),I=1,9),J=1,5)                   
 8999   CONTINUE                                                        
        DO 9049 I=1,9                                                   
        DO 9049 J=1,5                                                   
        DO 9049 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(2,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9049   CONTINUE                                                        
        WRITE(IW,9000)I9(1),((J4(I,J,1),I=1,9),J=1,5)                   
        DO 9099 K=2,4                                                   
        WRITE(IW,9010)I9(K),((J4(I,J,K),I=1,9),J=1,5)                   
 9099   CONTINUE                                                        
        DO 9149 I=1,9                                                   
        DO 9149 J=1,5                                                   
        DO 9149 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(3,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9149   CONTINUE                                                        
        WRITE(IW,9100)I9(1),JTM,((J4(I,J,1),I=1,9),J=1,5)               
        DO9199 K=2,4                                                    
        WRITE(IW,9110)I9(K),JTM,((J4(I,J,K),I=1,9),J=1,5)               
 9199   CONTINUE                                                        
        DO 9249 I=1,9                                                   
        DO 9249 J=1,5                                                   
        DO 9249 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(4,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9249   CONTINUE                                                        
        WRITE(IW,9200)I9(1),JTD,IP1,((J4(I,J,1),I=1,9),J=1,5)           
        DO 9299 K=2,4                                                   
        WRITE(IW,9210)I9(K),JTD,IP1,((J4(I,J,K),I=1,9),J=1,5)           
 9299   CONTINUE                                                        
        DO 9349 I=1,9                                                   
        DO 9349 J=1,5                                                   
        DO 9349 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(5,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9349   CONTINUE                                                        
        WRITE(IW,9300)I9(1),JTQ,IP2,((J4(I,J,1),I=1,9),J=1,5)           
        DO 9399 K=2,4                                                   
        WRITE(IW,9310)I9(K),JTQ,IP2,((J4(I,J,K),I=1,9),J=1,5)           
 9399   CONTINUE                                                        
        DO 9449 I=1,9                                                   
        DO 9449 J=1,5                                                   
        DO 9449 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(6,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9449   CONTINUE                                                        
        WRITE(IW,9400)I9(1),JTO,IP3,((J4(I,J,1),I=1,9),J=1,5)           
        DO 9499 K=2,4                                                   
        WRITE(IW,9410)I9(K),JTO,IP3,((J4(I,J,K),I=1,9),J=1,5)           
 9499   CONTINUE                                                        
        DO 9549 I=1,9                                                   
        DO 9549 J=1,5                                                   
        DO 9549 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(7,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9549   CONTINUE                                                        
        WRITE(IW,9500)I9(1),((J4(I,J,1),I=1,9),J=1,5)                   
        DO 9599 K=2,4                                                   
        WRITE(IW,9510)I9(K),((J4(I,J,K),I=1,9),J=1,5)                   
 9599   CONTINUE                                                        
        DO 9649 I=1,9                                                   
        DO 9649 J=1,5                                                   
        DO 9649 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(8,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9649   CONTINUE                                                        
        DO 9699 K=1,4                                                   
        WRITE(IW,9600)I9(K),((J4(I,J,K),I=1,9),J=1,5)                   
 9699   CONTINUE                                                        
        DO 9749 I=1,9                                                   
        DO 9749 J=1,5                                                   
        DO 9749 K=1,4                                                   
        J4(I,J,K) = JB                                                  
        I99 = ID9(J)                                                    
        IF(MOD(J9(9,I99)/2**(9-I),2).EQ.1)  J4(I,J,K) = JX(K)           
 9749   CONTINUE                                                        
        WRITE(IW,9700)I9(1),YC,((J4(I,J,1),I=1,9),J=1,5)                
        DO 9799 K=2,4                                                   
        WRITE(IW,9710)I9(K),YC,((J4(I,J,K),I=1,9),J=1,5)                
 9799   CONTINUE                                                        
        DO 9899 K=1,4                                                   
        WRITE(IW,9800)I9(K)                                             
 9899   CONTINUE                                                        
        DO 9999 K=1,4                                                   
        WRITE(IW,9900)I9(K)                                             
 9999   CONTINUE                                                        
        WRITE(IW,9910)                                                  
 4000   FORMAT(1H1/A1,120(1H*))                                         
 4010   FORMAT(A1,120(1H*))                                             
 4100   FORMAT(A1,1H*,22X,1HI,95X,1H*)                                  
 4200   FORMAT(A1,2H* ,9A1,2X,9A1,34H I TABLE OF ALL INPUT PARAMETERS -,
     1  60H-- DEFAULTS (IF APPROPRIATE) FOLLOW THE VALUES, IN PARENTHES,
     2  4HES *)                                                         
 4210   FORMAT(A1,1H*,1X,9A1,2X,9A1,1X,1HI,95X,1H*)                     
 4300   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,95(1H.),1H*)                     
 4310   FORMAT(A1,1H*,1X,9A1,2X,9A1,1X,1HI,95X,1H*)                     
 4400   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,95X,1H*)                         
 4410   FORMAT(A1,1H*,1X,9A1,2X,9A1,1X,1HI,95X,1H*)                     
 4500   FORMAT(A1,2H* ,9A1,2X,9A1,15H I E12 INIT. N=,I2,12H(MAX=20,DEF=,
     1  23H15) E11 L-DIST. OPTION=,I2,30H(DEF=0 /-1=INPUTED,0=STATIST.,,
     2  14H2=QUADRATIC) *)                                              
 4510   FORMAT(A1,1H*,1X,9A1,2X,9A1,1X,1HI,13X,I2,35X,I2,43X,1H*)       
 4600   FORMAT(A1,2H* 9A1,2X,9A1,24H I E12 STAT. DIST. EXP.=,F7.5,      
     1  25H(0.0) E14 QUAD. PARAM./A=,F7.5,8H(0.0),B=,F7.5,              
     2  20H(0.0) 1+A*L+B*L**2 *)                                        
 4610   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,22X,F7.5,25X,F7.5,8X,F7.5,19X,   
     1  1H*)                                                            
 4700   FORMAT(A1,2H* ,9A1,2X,9A1,17H I E13 NORM.INIT/,F7.6,9(1X,F7.6), 
     1  2H *)                                                           
 4710   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,15X,10(F7.6,1X),1H*)             
 4720   FORMAT(A1,2H* ,9A1,2X,9A1,7H I E13 ,24H  L-DISTRIBUTION EXTENDS,
     1  60H BEYOND STARTING N.  SEE NEXT PAGE FOR COMPLETE DISTRIBUTION,
     2  1H.,5X,1H*)                                                     
 4730   FORMAT(A1,1H*,1X,9A1,2X,9A1,2H I,95X,1H*)                       
 4800   FORMAT(A1,2H* ,9A1,2X,9A1,17H I L-DIST.(0-19)/,10(F7.6,1X),1H*) 
 4810   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,15X,10(F7.6,1X),1H*)             
 4820   FORMAT(A1,2H* ,9A1,2X,9A1,3H I ,4X,8(9X,1H*),10X,1H*)           
 4830   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,95X,1H*)                         
 4900   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,95(1H.),1H*)                     
 4910   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,95X,1H*)                         
 5000   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,45X,1HI,49X,1H*)                 
 5010   FORMAT(A1,2H* ,9A1,2X,9A1,2H I,45X,1HI,49X,1H*)                 
 5100   FORMAT(A1,1H*,22X,42HI E07 EFFECTIVE CHARGE FOR ELECTRONIC SHEL,
     1  55HLS  I E08  POPUL. OF EL. SUBSHELLS (FRACTION OF FULL) *)     
 5110   FORMAT(A1,1H*,22X,1HI,45X,1HI,49X,1H*)                          
 5200   FORMAT(A1,27H*  E01  Z /  (NEEDED)  I K/,F6.3,4H, L/,F6.3,4H, M/
     1  ,F6.3,21H  (ALL NEEDED)  I 1S=,F5.3,11H(1.000) 2S=,F5.3,        
     2  11H(1.000) 2P=,F5.3,9H(1.000) *)                                
 5210   FORMAT(A1,1H*,22X,1HI,3X,F6.3,4X,F6.3,4X,F6.3,16X,1HI,4X,F5.3,  
     1  11X,F5.3,11X,F5.3,8X,1H*)                                       
 5300   FORMAT(A1,1H*,22(1H-),1H+,17(1H-),1H+,27(1H-),5H+ 3S=,F5.3,     
     1  11H(1.000) 3P=,F5.3,11H(1.000) 3D=,F5.3,9H(1.000) *)            
 5310   FORMAT(A1,1H*,22X,1H+,17X,1H+,27X,1H+,4X,F5.3,11X,F5.3,11X,F5.3,
     1  8X,1H*)                                                         
 5400   FORMAT(A1,50H*   E10 DEPLETION OF ELECTRONIC SHELLS   I E09 ELE,
     1  20HCT. 1S WIDTH IN EV +,49(1H-),1H*)                            
 5410   FORMAT(A1,1H*,40X,1HI,27X,1H+,49X,1H*)                          
 5500   FORMAT(A1,1H*,5X,2HK/,I1,3H(0),6X,2HL/,I1,3H(0),6X,2HM/,I1,3H(0)
     1  ,5X,1HI,6X,F7.3,9H(000.000),5X,29HI E02 AV. EL. BINDING ENER. F,
     2  22HOR ATOM Z-1 (NEEDED) *)                                      
 5510   FORMAT(A1,1H*,7X,I1,11X,I1,11X,I1,8X,1HI,6X,F7.3,14X,1HI,49X,1H*
     1  )                                                               
 5600   FORMAT(A1,50H* (0=YES,1=NO - IF K/0,1S WIDTH IS USED) I (EXPER.,
     1  23H OR INPUTED VALUE) I K/,F9.2,8H(EV)  L/,F8.2,8H(EV)  M/,F8.2,
     2  6H(EV) *)                                                       
 5610   FORMAT(A1,1H*,40X,1HI,27X,1HI,3X,F9.2,8X,F8.2,8X,F8.2,5X,1H*)   
 5700   FORMAT(A1,1H*,40(1H-),1H+,27(1H-),1H+,49(1H-),1H*)              
 5710   FORMAT(A1,1H*,40X,1H+,27X,1H+,49X,1H*)                          
 5800   FORMAT(A1,20H* E19 ATOMIC WEIGHT=,F6.2,21H(140.00) E23 FERMI PA,
     1  34HRAMETERS (NOT USED IN PROGRAM)  C=,F7.5,17H(FM),SKIN THICK. ,
     2  2HT=,F7.5,6H(FM) *)                                             
 5810   FORMAT(A1,1H*,19X,F6.2,55X,F7.5,19X,F7.5,5X,1H*)                
 5820   FORMAT(A1,20H* E19 ATOMIC WEIGHT=,F6.2,21H(140.00) E23 FERMI PA,
     1  60HRAMETERS (NOT USED IN PROGRAM)    * *  N O T   S P E C I F I,
     2  13H E D  * *   *)                                               
 5830   FORMAT(A1,1H*,19X,F6.2,93X,1H*)                                 
 5900   FORMAT(A1,50H* E24 /DIRAC/ PROGRAM PARAMETERS(NOT USED)/ STEP I,
     1  14HN INTEGRATION=,1PE9.3,17H MATCHING RADIUS=,E9.3,             
     2  21H(FM)(FOR REFERENCE) *)                                       
 5910   FORMAT(A1,1H*,63X,1PE9.3,17X,E9.3,20X,1H*)                      
 5920   FORMAT(A1,50H* E24 /DIRAC/ PROGRAM PARAMETERS(NOT USED)/   *  *,
     1  46H    N  O  T     S  P  E  C  I  F  I  E  D    *,8(2X,1H*))    
 5930   FORMAT(A1,1H*,118X,1H*)                                         
 6000   FORMAT(A1,23H* E30 MASSES/ PARTICLE=,F9.4,18H(206.7686)(ELEC. M,
     1  17HASSES)  ELECTRON=,F8.1,24H(511003.4)(EV)  NUCLEON=,F6.2,     
     2  15H(931.48)(MEV) *)                                             
 6010   FORMAT(A1,1H*,22X,F9.4,35X,F8.1,24X,F6.2,14X,1H*)               
 6100   FORMAT(A1,50H* E06,E07 SPECIAL EXPERIM. TRANSITION ENERGIES/  2,
     1  5HP-1S=,-6PF8.6,25H(EMPIR. FIT)(MEV)  2S-2P=,F8.6,10H(0.0/NO TR,
     2  14HANSIT.)(MEV) *)                                              
 6110   FORMAT(A1,1H*,54X,-6PF8.6,25X,F8.6,23X,1H*)                     
 6200   FORMAT(A1,1H*,12(1H-),1H+,31(1H-),1H+,55(1H-),1H+,17(1H-),1H*)  
 6210   FORMAT(A1,1H*,12X,1H+,31X,1H+,55X,1H+,17X,1H*)                  
 6300   FORMAT(A1,1H*,12X,14H/ E20 HI. CUT=,F6.3,17H(20.000)MEV I E21,  
     1  16H INTENS. CUTOFF=,1PE9.3,33H(1.000E-06)(PER PARTICLE) I E33 S,
     2  12HTAR OPTION *)                                                
 6310   FORMAT(A1,1H*,12X,1H/,13X,F6.3,12X,1HI,20X,1PE9.3,26X,1HI,17X,  
     1  1H*)                                                            
 6400   FORMAT(A1,27H* CATALOGUE  / E20 LOW CUT=,F6.3,14H( 0.040)MEV I ,
     1  60HE25 CALIBRATION PARAMETERS (CONVERSION TO CHANNEL NO) I IN C,
     2  7HATALOG/,I1,5H(0) *)                                           
 6410   FORMAT(A1,1H*,12X,1H/,13X,F6.3,12X,1HI,55X,1HI,12X,I1,4X,1H*)   
 6500   FORMAT(A1,27H* PARAMETERS / E22 RESOL. =,F6.5,13H(.00030)MEV I, 
     1  5X,2HA=,F7.2,11H (KEV) , B=,F7.3,27H  CHANNEL NO =(E-A)/B) I 0=,
     2  15HDEF, 1=READIN *)                                             
 6510   FORMAT(A1,1H*,12X,1H/,13X,F6.5,12X,1HI,7X,F7.2,11X,F7.3,23X,    
     1  1HI,17X,1H*)                                                    
 6520   FORMAT(A1,27H* PARAMETERS / E22 RESOL. =,F6.5,13H(.00030)MEV I, 
     1  5X,57H  *  *   N O T   S P E C I F I E D   *  *  *  *   I 0=DEF,
     2  12H, 1=READIN *)                                                
 6530   FORMAT(A1,1H*,12X,1H/,13X,F6.5,12X,1HI,55X,1HI,17X,1H*)         
 6600   FORMAT(A1,1H*,12X,24H/ E34 INTENSITIES /  5*=,1PE7.1,8H(.1) 4*=,
     1  E7.1,9H(.01) 3*=,E7.1,10H(.001) 2*=,E7.1,11H(.0001) 1*=,E7.1,   
     2  10H(.00001) *)                                                  
 6610   FORMAT(A1,1H*,12X,1H/,23X,1PE7.1,8X,E7.1,9X,E7.1,10X,E7.1,11X,  
     1  E7.1,9X,1H*)                                                    
 6700   FORMAT(A1,18H* E17,E18 QUAN.DEP,4(1H/,I2,1H,,I2,1H,,I2,1H,,I2,  
     1  1H,,I2),28H(ALL/-1=START RAN) DO DEPOL=,I1,13H(0/0=Y,1=N) *)    
 6710   FORMAT(A1,1H*,17X,20(1X,I2),28X,I1,12X,1H*)                     
 6800   FORMAT(A1,1H*,118X,1H*)                                         
 6900   FORMAT(A1,120(1H*))                                             
 7000   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7100   FORMAT(A1,50H* S H E L L  A N D  S U B S H E L L   C O M B I N ,
     1  60HA T I O N S  I      B O O K K E E P I N G    P A R A M E T E,
     2  10H R S     *)                                                  
 7110   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7200   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7210   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7300   FORMAT(A1,16H* E38 CASES/  M/,I1,16H(MAX=3,DEF=3),D/,I1,3H,Q/,  
     1  I1,3H,O/,I1,48H (D,Q,O/MAX=7,DEF=4) I E31 LOGICAL UNIT NO.S REA,
     2  2HD=,I1,10H(5),WRITE=,I1,10H(6),PUNCH=,I1,5H(7) *)              
 7310   FORMAT(A1,1H*,15X,I1,16X,I1,3X,I1,3X,I1,21X,1HI,28X,I1,10X,I1,  
     1  10X,I1,4X,1H*)                                                  
 7400   FORMAT(A1,1H*,62(1H.),1HI,55(1H.),1H*)                          
 7410   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7500   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 7600   FORMAT(A1,8H* E39,40,1X,2HC1,6X,2HC2,6X,2HC3,6X,2HC4,6X,2HC5,6X,
     1  2HC6,6X,2HC7,4X,20HI E35 PRINT OPTION =,I2,17H(0/PRINT ALL) SEL,
     2  18HECT CODES 0 - 63 *)                                          
 7610   FORMAT(A1,1H*,62X,1HI,19X,I2,34X,1H*)                           
 7700   FORMAT(A1,7H* MON/ ,A2,6H(1S)  ,A2,6H(2T)  ,A2,13H(3T)  ------ ,
     1  46H ------  ------  ------  I E36 DEBUG OPTION = ,I1,           
     2  35H(0)(0=N,1=Y) DEB PRINTS ALL RATES *)                         
 7710   FORMAT(A1,1H*,6X,A2,6X,A2,6X,A2,38X,1HI,20X,I1,34X,1H*)         
 7800   FORMAT(A1,7H* DIP/ ,A2,6H(RD)* ,A2,6H(1S)  ,A2,6H(2T)  ,A2,     
     1  6H(3T)  ,3(A2,6H(--)  ),21HI E32 ACCURACY CHECK=,I1,            
     2  35H(3)(0=NO, 1 OR 2=PARTIAL, 3=FULL) *)                         
 7810   FORMAT(A1,1H*,6X,7(A2,6X),1HI,20X,I1,34X,1H*)                   
 7900   FORMAT(A1,7H* QUA/ ,A2,6H(RD)* ,A2,6H(1S)  ,A2,6H(2T)  ,A2,     
     1  6H(3T)  ,3(A2,6H(--)  ),1HI,55(1H.),1H*)                        
 7910   FORMAT(A1,1H*,6X,7(A2,6X),1HI,55X,1H*)                          
 8000   FORMAT(A1,7H* OCT/ ,A2,6H(RD)* ,A2,6H(1S)  ,A2,6H(2T)  ,A2,     
     1  6H(3T)  ,3(A2,6H(--)  ),1HI,55X,1H*)                            
 8010   FORMAT(A1,1H*,6X,7(A2,6X),1HI,55X,1H*)                          
 8100   FORMAT(A1,50H* KEY/ RD=RADIAT., T=TOTAL SHELL, --=NOT APPLIC., ,
     1  30H*=MUST BE RD I E37 FACT. DIV.=,F5.2,21H(15.00) FACT. STORED ,
     2  14HFAC(N)/FD**N *)                                              
 8110   FORMAT(A1,1H*,62X,1HI,16X,F5.2,34X,1H*)                         
 8200   FORMAT(A1,1H*,62X,36HI E36 NO. OF DIRAC ENERGIES INPUTED=,I3,   
     1  8H OUT OF ,I3,7H MAX. *)                                        
 8210   FORMAT(A1,1H*,62X,1HI,35X,I3,8X,I3,6X,1H*)                      
 8300   FORMAT(A1,63(1H*),37HI E27 MAX NO. OF TRANSITIONS PUNCHED=,I3,  
     1  17H(MAX=200,DEF=0) *)                                           
 8310   FORMAT(A1,63(1H*),1HI,36X,I3,16X,1H*)                           
 8400   FORMAT(A1,1H*,62X,39HI E29 PUNCH SPECIFIED TRANSITIONS /    ,I1,
     1  17H(0 =YES, 1 =NO) *)                                           
 8410   FORMAT(A1,1H*,62X,1HI,38X,I1,16X,1H*)                           
 8500   FORMAT(A1,1H*,6X,43HS E L E C T I O N    O F    P E N E T R A T,
     1  6H I O N,7X,1HI,55(1H.),1H*)                                    
 8510   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 8600   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 8700   FORMAT(A1,50H* E42 MAX NUMBER OF TERMS IN I E41 PENETRATION SEL,
     1  60HECTION CODES I E28 PUNCHED CARD IDENTITY NO. IN COL.S 73-78 ,
     2  10H/(10000) *)                                                  
 8710   FORMAT(A1,1H*,28X,1HI,33X,1HI,55X,1H*)                          
 8800   FORMAT(A1,50H* PENETRATION (MAX=3,DEF.=1) I     (1=Y,0=N *=MUST,
     1  14H BE 0,DEF=1) I,55X,1H*)                                      
 8810   FORMAT(A1,1H*,28X,1HI,33X,1HI,55X,1H*)                          
 8900   FORMAT(A1,1H*,28X,17HI CASES AS IN E26,17X,1HI,5(1X,9A1,1X),1H*)
 8910   FORMAT(A1,1H*,28X,1HI,33X,1HI,5(1X,9A1,1X),1H*)                 
 9000   FORMAT(A1,50H*    1S  2S  2P  3S  3P  3D  I      C1  C2  C3  C4,
     1  14H  C5  C6  C7 I,5(1X,9A1,1X),1H*)                             
 9010   FORMAT(A1,1H*,28X,1HI,33X,1HI,5(1X,9A1,1X),1H*)                 
 9100   FORMAT(A1,4H* M/,1X,6(1X,I1,2X),5HI  M/,5H  -  ,6(2X,1H-,1X),1HI
     1  ,5(1X,9A1,1X),1H*)                                              
 9110   FORMAT(A1,1H*,4X,6(1X,I1,2X),1HI,33X,1HI,5(1X,9A1,1X),1H*)      
 9200   FORMAT(A1,5H* D/ ,6(1X,I1,2X),7HI  D/  ,I1,1H*,6(3X,I1),2H I,   
     1  5(1X,9A1,1X),1H*)                                               
 9210   FORMAT(A1,1H*,4X,6(1X,I1,2X),1HI,6X,I1,1X,6(3X,I1),2H I,        
     1  5(1X,9A1,1X),1H*)                                               
 9300   FORMAT(A1,5H* Q/ ,6(1X,I1,2X),7HI  Q/  ,I1,1H*,6(3X,I1),2H I,   
     1  5(1X,9A1,1X),1H*)                                               
 9310   FORMAT(A1,1H*,4X,6(1X,I1,2X),1HI,6X,I1,1X,6(3X,I1),2H I,        
     1  5(1X,9A1,1X),1H*)                                               
 9400   FORMAT(A1,5H* O/ ,6(1X,I1,2X),7HI  O/  ,I1,1H*,6(3X,I1),2H I,   
     1  5(1X,9A1,1X),1H*)                                               
 9410   FORMAT(A1,1H*,4X,6(1X,I1,2X),1HI,6X,I1,1X,6(3X,I1),2H I,        
     1  5(1X,9A1,1X),1H*)                                               
 9500   FORMAT(A1,1H*,28(1H.),1HI,33(1H.),1HI,5(1X,9A1,1X),1H*)         
 9510   FORMAT(A1,1H*,28X,1HI,33X,1HI,5(1X,9A1,1X),1H*)                 
 9600   FORMAT(A1,1H*,62X,1HI,5(1X,9A1,1X),1H*)                         
 9700   FORMAT(A1,17H* E43 Y-CUTOFF M/,F5.2,3H D/,F5.2,3H Q/,F5.2,3H O/,
     1  F5.2,18H (ALL DEF/ 1.00) I,5(1X,9A1,1X),1H*)                    
 9710   FORMAT(A1,1H*,16X,F5.2,3(3X,F5.2),17X,1HI,5(1X,9A1,1X),1H*)     
 9800   FORMAT(A1,1H*,62X,1HI,55X,1H*)                                  
 9900   FORMAT(A1,120(1H*))                                             
 9910   FORMAT(1H1)                                                     
        IF(IP8.EQ.0)  GO TO 11100                                       
        WRITE(IW,10000)                                                 
        WRITE(IW,10100)                                                 
        WRITE(IW,10200)                                                 
        DO 10199 I=1,NMAX                                               
        N = NMAX + 1 -I                                                 
        SS =0.000                                                       
        DO 10099 J=1,N                                                  
        K = N*(N-1)/2 + J                                               
        PLX(J) = PLN(K)                                                 
        SS = SS + PLX(J)                                                
        K = IFIX(PLX(J)*1.0E6+0.50)                                     
        IXX(1,J) = MOD(K/100,10)                                        
        IXX(2,J) = MOD(K/10,10)                                         
        IXX(3,J) = MOD(K,10)                                            
        IF(AMOD(1000.0*PLX(J),1.0).GT..499999)  PLX(J)=PLX(J)-.000499999
        IF(PLX(J).LT.0.000)  PLX(J) = 0.000                             
10099   CONTINUE                                                        
        IF(N.GT.2)  WRITE(IW,10300)N,SS,(PLX(J),J=1,N)                  
        IF(N.EQ.2)  WRITE(IW,10600)N,SS,(PLX(J),J=1,2)                  
        IF(N.EQ.1)  WRITE(IW,10900)N,SS,PLX(1)                          
        IF(N.GT.2)  WRITE(IW,10400)((IXX(II,J),II=1,3),J=1,N)           
        IF(N.EQ.2)  WRITE(IW,10700)((IXX(II,J),II=1,3),J=1,2)           
        IF(N.EQ.1)  WRITE(IW,11000)(IXX(II,1),II=1,3)                   
        IF(N.GT.15)  WRITE(IW,10450)                                    
        IF(N.GT.10.AND.N.LE.15)  WRITE(IW,10460)                        
        IF(N.GT.5.AND.N.LE.10)  WRITE(IW,10470)                         
        IF(N.LE.5.AND.N.GT.3)  WRITE(IW,10480)                          
        IF(N.EQ.3)  WRITE(IW,10500)                                     
        IF(N.EQ.2)  WRITE(IW,10800)                                     
10199   CONTINUE                                                        
10000   FORMAT(20X,49H*  *   N O R M A L I Z E D   I N I T I A L   L - ,
     1  30HD I S T R I B U T I O N   *  */)                             
10100   FORMAT(7X,50HTOTAL  I  L=0  L=1  L=2  L=3  L=4 I  L=5  L=6  L=7,
     1  60H  L=8  L=9 I L=10 L=11 L=12 L=13 L=14 I L=15 L=16 L=17 L=18 ,
     2  4HL=19)                                                         
10200   FORMAT(1X,13(1H-),1H+,26(1H-),1H+,26(1H-),1H+,26(1H-),1H+,26(1H-
     1  ))                                                              
10300   FORMAT(3H N=,I2,1X,F7.5,4(2H I,1X,F4.3,1X,F4.3,1X,F4.3,1X,F4.3, 
     1  1X,F4.3))                                                       
10400   FORMAT(13X,4(2H I,2X,3I1,2X,3I1,2X,3I1,2X,3I1,2X,3I1))          
10450   FORMAT(14X,1HI,26X,1HI,26X,1HI,26X,1HI)                         
10460   FORMAT(14X,1HI,26X,1HI,26X,1HI)                                 
10470   FORMAT(14X,1HI,26X,1HI)                                         
10480   FORMAT(14X,1HI)                                                 
10500   FORMAT(14X,1HI,56X,1H+,32(1H-),1H+)                             
10600   FORMAT(3H N=,I2,1X,F7.5,3H I ,F4.3,1X,F4.3,46X,1HI,32X,1HI)     
10700   FORMAT(14X,1HI,2X,3I1,2X,3I1,46X,27HI ENTRIES FOLDED  .ABC = .A,
     1  7HBCDEF I)                                                      
10800   FORMAT(14X,1HI,56X,34HI TO SAVE SPACE    DEF           I)       
10900   FORMAT(3H N=,I2,F8.5,3H I ,F4.3,51X,1HI,32X,1HI)                
11000   FORMAT(14X,1HI,2X,3I1,51X,1H+,32(1H-),1H+/1H1)                  
11100   CALL CASCAD                                                     
        CALL SORT                                                       
C
C       Now we go back to the start
        GO TO 1175                                                      
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE RREAD                                                
        INTEGER A,S,B                                                   
        DOUBLE PRECISION F                                              
C  ***  READS INPUT CARDS AND SETS PARAMETERS ACCORDING TO CODES        
        DIMENSION A(89),S(70)                                           
        COMMON/LOC001/IJK,ENERG,ECONS,ECONST,D2P1SM,D2P1S               
        COMMON/LOC002/BEM(3),ZSA(3),BE(3)                               
        COMMON/LOC003/K0,K1,K2,K3                                       
        COMMON/LOC004/NN0(3),NN1(7),NN2(7),NN3(7)                       
        COMMON/LOC006/IP1(7),IP2(7),IP3(7),IQ1(7),IQ2(7),IQ3(7)         
        COMMON/LOC007/IC                                                
        COMMON/LOC008/IR,IW,IPUNCH,IPRINT                               
        COMMON/LOC009/F(60),FD                                          
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC017/IFM(6)                                            
        COMMON/LOC018/IFD(9)                                            
        COMMON/LOC019/IFQ(10)                                           
        COMMON/LOC020/IFO(10)                                           
        COMMON/LOC028/M1(7),M2(7),M3(7),YC(4),IDB                       
        COMMON/LOC030/POP(6),JTM(6),JTD(6),JTQ(6),JTO(6)                
        COMMON/LOC031/JM(10),JD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC033/M,E(1000),AI(1000),IA(1000),ENERGY(20,40)         
        COMMON/LOC035/ICC,CD(5),EA,EB,IDR                               
        COMMON/LOC036/ZMK,ZML,ZMM,ZMKM,ZMLM,ZMMM,IVERS                  
        COMMON/LOC037/PL(20),NPOL(20),IPOL,CL1,CL2,IDE,PLN(210),IP8     
        COMMON/LOC038/AA,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)              
        COMMON/LOC039/NOPT,NMAX,ALEXP                                   
        COMMON/LOC040/AMASSA,AMASSN,HBAR                                
        COMMON/LOC041/MPU,ICPU(200),IPN                                 
        DATA T8/0.000/                                                  
        DATA A/3HIJK,3HENE,3HAMA,3HAMN,3HD21,3HBE ,3HK  ,3HNN0,3HNN1,   
     1  3HNN2,3HNN3,3HIP1,3HIP2,3HIP3,3HIC ,3HIRE,3HIWR,3HIPU,3HFD ,    
     2  3HAMM,3HAME,3HZ  ,3HZS ,3HIFM,3HIFD,3HIFQ,3HIFO,3HM1 ,3HM2 ,    
     3  3HM3 ,3HPOP,3HJM ,3HJD ,3HJQ ,3HJO ,3HEHI,3HELO,3HCLM,3HERS,    
     4  3HICC,3HCD ,3HEAB,3HSTO,3HC  ,3HXEQ,3HA  ,3HCT ,3HSTP,3HKWD,    
     5  3HNOP,3HNMX,3HPL ,3HDIR,3HIPR,3HNPL,3HIPL,3HZM ,3HZMM,3HIQ1,    
     6  3HIQ2,3HIQ3,3HIDB,3HYC ,3HIYC,3HIJ ,3HYJ ,3HJJ1,3HIPC,3HCL ,    
     7  3HESP,3HPUN,3HIDE,3HIPN,3HJTM,3HJTD,3HJTQ,3HJTO,3HPLN,3HIP ,    
     8  3H   ,3H   ,3H   ,3H   ,3H   ,3H   ,3H   ,3H   ,3H   ,3H   /    
        DATA IRC/0/                                                     
C       IRC starts as 0, so MPU = 0?
        IF(IRC.EQ.0)  MPU = 0                                           
        IF(IRC.NE.0)  WRITE(IW,50)                                      
   50   FORMAT(//51H *** NEW CASE COMING UP *** READING UNTIL XEQ ***  )
C
C       IRC is counting the number of lines in the file
  100   IRC = IRC+1                                                     
C
        READ(IR,200)B,B1,B2,S                                           
  200   FORMAT(A3,A3,A4,70A1)                                           
C      
C       Write out entire input file at top of output file
     0  IF(MOD(IPRINT/2,2).EQ.0)  WRITE(IW,300)IRC,B,B1,B2,S            
  300   FORMAT(1X,14HINPUT CARD NO.,I3,5X,3H---,5X,2A3,A4,70A1)         
     0  J = 0                                                           
        DO 400 I=1,89                                                   
C
C       Here we have to find if our keyword matches something in the data 
C       structure. If we have a match then we make a jump
C
        IF(B.EQ.A(I))  J=I                                              
        IF(J.NE.0)  GO TO 600                                           
  400   CONTINUE                                                        
        WRITE(IW,500)B                                                  
  500   FORMAT(/46H *** ERROR *** INPUT ACTION CODE ILLEGAL *** =,A3,   
     1  21H *** CARD IGNORED ***)                                       
     0  GO TO 100                                                       
C       
C       We arrive here from line 925 and do this weird modular logic
C       JT is the row of the data type which we are in, and JU is the element
  600   JT = J/10 + 1                                                   
C       write(4,*)JT, "JT"
        JU = MOD(J,10) + 1                                              
C       write(4,*)JU, "JU"
        NI = 1                                                          
C
C       Computed go to, chooses JT'th element and goes there
        GO TO(999,1999,2999,3999,4999,5999,6999,7999,8999),JT           
  900   FORMAT(53H *** WARNING *** UNIMPLEMENTED USER CODE FOR THIS VER,
     1  29HSION (NO EFFECT PRODUCED) ***)                               
  999   GO TO (100,1100,1200,1300,1400,1500,1600,1700,1800,1900),JU     
 1100   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,0,DUM,IJK)                                     
        GO TO 100                                                       
 1200   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,1,ENERG,IDUM)                                  
        GO TO 100                                                       
 1300   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,1,AMASSA,IDUM)                                 
        WRITE(IW,1350)                                                  
 1350   FORMAT(53H *** WARNING *** AMASSA HAS BEEN PROVIDED -- IF NEW A,
     1  51H NEW AMASSA MUST BE GIVEN (OR ZERO FOR DEFAULT) ***)         
     0  GO TO 100                                                       
 1400   CALL DCYPHR(S,NI,1,AMASSN,IDUM)                                 
        GO TO 100                                                       
 1500   CALL DCYPHR(S,NI,1,D2P1S,IDUM)                                  
        GO TO 100                                                       
 1600   DO 1650 K=1,3                                                   
        CALL DCYPHR(S,NI,1,BE(K),IDUM)                                  
 1650   CONTINUE                                                        
        GO TO 100                                                       
 1700   CALL DCYPHR(S,NI,0,DUM,K0)                                      
        CALL DCYPHR(S,NI,0,DUM,K1)                                      
        CALL DCYPHR(S,NI,0,DUM,K2)                                      
        CALL DCYPHR(S,NI,0,DUM,K3)                                      
        GO TO 100                                                       
 1800   DO 1850 K=1,K0                                                  
        CALL DCYPHR(S,NI,0,DUM,NN0(K))                                  
 1850   CONTINUE                                                        
        GO TO 100                                                       
 1900   DO 1950 K=1,K1                                                  
        CALL DCYPHR(S,NI,0,DUM,NN1(K))                                  
 1950   CONTINUE                                                        
        IF(NN1(1).NE.0)  WRITE(IW,1975)                                 
        GO TO 100                                                       
 1975   FORMAT(53H *** ERROR *** RADIARION IS NOT COMPUTED AS THE FIRST,
     1  60H ITEM IN THE MULTIPOLARITY *** NO ACTION TAKEN BUT RESULTS A,
     2  8HRE BAD *)                                                     
 1999   GO TO (2000,2100,2200,2300,2400,2500,2600,2700,2800,2900),JU    
 2000   DO 2050 K=1,K2                                                  
        CALL DCYPHR(S,NI,0,DUM,NN2(K))                                  
 2050   CONTINUE                                                        
        IF(NN2(1).NE.0)  WRITE(IW,1975)                                 
        GO TO 100                                                       
 2100   DO 2150 K=1,K3                                                  
        CALL DCYPHR(S,NI,0,DUM,NN3(K))                                  
 2150   CONTINUE                                                        
        IF(NN3(1).NE.0)  WRITE(IW,1975)                                 
        GO TO 100                                                       
 2200   DO 2250 K=1,K1                                                  
        CALL DCYPHR(S,NI,0,DUM,IP1(K))                                  
 2250   CONTINUE                                                        
        IF(IP1(1).NE.0)  WRITE(IW,2275)                                 
        GO TO 100                                                       
 2275   FORMAT(53H *** ERROR *** ILLEGAL CODE FOR PENETRATION IN RADIAT,
     1  60HION (1ST ENTRY MUST BE 0) *** NO ACTION TAKEN, BUT RESULTS A,
     2  8HRE BAD *)                                                     
 2300   DO 2350 K=1,K2                                                  
        CALL DCYPHR(S,NI,0,DUM,IP2(K))                                  
 2350   CONTINUE                                                        
        IF(IP2(1).NE.0)  WRITE(IW,2275)                                 
        GO TO 100                                                       
 2400   DO 2450 K=1,K3                                                  
        CALL DCYPHR(S,NI,0,DUM,IP3(K))                                  
 2450   CONTINUE                                                        
        IF(IP3(1).NE.0)  WRITE(IW,2275)                                 
        GO TO 100                                                       
 2500   CALL DCYPHR(S,NI,0,DUM,IC)                                      
        GO TO 100                                                       
 2600   CALL DCYPHR(S,NI,0,DUM,IR)                                      
        GO TO 100                                                       
 2700   CALL DCYPHR(S,NI,0,DUM,IW)                                      
        GO TO 100                                                       
 2800   CALL DCYPHR(S,NI,0,DUM,IPUNCH)                                  
        GO TO 100                                                       
 2900   CALL DCYPHR(S,NI,1,FD,IDUM)                                     
        GO TO 100                                                       
 2999   GO TO (3000,3100,3200,3300,3400,3500,3600,3700,3800,3900),JU    
 3000   CALL DCYPHR(S,NI,1,AMASSM,IDUM)                                 
        GO TO 100                                                       
 3100   CALL DCYPHR(S,NI,1,AMASSE,IDUM)                                 
        GO TO 100                                                       
C       We arrive here from JU = JT = 3, and we grab Z from S
 3200   CALL DCYPHR(S,NI,1,Z,IDUM)                                      
        GO TO 100                                                       
 3300   CALL DCYPHR(S,NI,1,ZSK,IDUM)                                    
        CALL DCYPHR(S,NI,1,ZSL,IDUM)                                    
        CALL DCYPHR(S,NI,1,ZSM,IDUM)                                    
        GO TO 100                                                       
 3400   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 3450 K=1,6                                                   
        CALL DCYPHR(S,NI,0,DUM,IFM(K))                                  
 3450   CONTINUE                                                        
        GO TO 100                                                       
 3500   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 3550 K=1,9                                                   
        CALL DCYPHR(S,NI,0,DUM,IFD(K))                                  
 3550   CONTINUE                                                        
        GO TO 100                                                       
 3600   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 3650 K=1,10                                                  
        CALL DCYPHR(S,NI,0,DUM,IFQ(K))                                  
 3650   CONTINUE                                                        
        GO TO 100                                                       
 3700   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 3750 K=1,10                                                  
        CALL DCYPHR(S,NI,0,DUM,IFO(K))                                  
 3750   CONTINUE                                                        
        GO TO 100                                                       
 3800   DO 3850 K=1,K1                                                  
        CALL DCYPHR(S,NI,0,DUM,M1(K))                                   
 3850   CONTINUE                                                        
        GO TO 100                                                       
 3900   DO 3950 K=1,K2                                                  
        CALL DCYPHR(S,NI,0,DUM,M2(K))                                   
 3950   CONTINUE                                                        
        GO TO 100                                                       
 3999   GO TO (4000,4100,4200,4300,4400,4500,4600,4700,4800,4900),JU    
 4000   DO 4050 K=1,K3                                                  
        CALL DCYPHR(S,NI,0,DUM,M3(K))                                   
 4050   CONTINUE                                                        
        GO TO 100                                                       
 4100   DO 4150 K=1,6                                                   
        CALL DCYPHR(S,NI,1,POP(K),IDUM)                                 
        IF(POP(K).LT.0.000)  WRITE(IW,4110)K,POP(K)                     
 4110   FORMAT(21H *** WARNING *** POP(,I1,3H) =,F6.3,14H *** SET TO 1.,
     1  7H000 ***)                                                      
        IF(POP(K).LT.0.000)  POP(K) = 0.000                             
        IF(POP(K).GT.1.000)  WRITE(IW,4120)K,POP(K)                     
 4120   FORMAT(21H *** WARNING *** POP(,I1,3H) =,F6.3,14H *** SET TO 0.,
     1  7H000 ***)                                                      
        IF(POP(K).GT.1.000)  POP(K) = 1.000                             
 4150   CONTINUE                                                        
        GO TO 100                                                       
 4200   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 4250 K=1,10                                                  
        CALL DCYPHR(S,NI,0,DUM,JM(K))                                   
 4250   CONTINUE                                                        
        GO TO 100                                                       
 4300   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVRES.NE.1)  GO TO 100                                       
        DO 4350 K=1,14                                                  
        CALL DCYPHR(S,NI,0,DUM,JD(K))                                   
 4350   CONTINUE                                                        
        GO TO 100                                                       
 4400   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 4450 K=1,15                                                  
        CALL DCYPHR(S,NI,0,DUM,JQ(K))                                   
 4450   CONTINUE                                                        
        GO TO 100                                                       
 4500   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 4550 K=1,15                                                  
        CALL DCYPHR(S,NI,0,DUM,JO(K))                                   
 4550   CONTINUE                                                        
        GO TO 100                                                       
 4600   CALL DCYPHR(S,NI,1,EHIGH,IDUM)                                  
        GO TO 100                                                       
 4700   CALL DCYPHR(S,NI,1,ELOW,IDUM)                                   
        GO TO 100                                                       
 4800   CALL DCYPHR(S,NI,1,CLIMIT,IDUM)                                 
        GO TO 100                                                       
 4900   CALL DCYPHR(S,NI,1,ERES,IDUM)                                   
        GO TO 100                                                       
 4999   GO TO (5000,5100,5200,5300,100,5500,5600,5700,5800,5900),JU     
 5000   CALL DCYPHR(S,NI,0,DUM,ICC)                                     
        GO TO 100                                                       
 5100   DO 5150 K=1,5                                                   
        CALL DCYPHR(S,NI,1,CD(K),IDUM)                                  
 5150   CONTINUE                                                        
        GO TO 100                                                       
 5200   CALL DCYPHR(S,NI,1,EA,IDUM)                                     
        CALL DCYPHR(S,NI,1,EB,IDUM)                                     
        GO TO 100                                                       
 5300   WRITE(IW,5350)                                                  
        STOP                                                            
 5350   FORMAT(53H *** STOP CARD ENCOUNTERED IN THE INPUT FILE *** EXEC,
     1  60HTION TERMINATED NORMALLY *** NO SUBSEQUENT INPUT CARDS (IF A,
     2  8HNY) READ)                                                     
C       We are jumping here from the top when XEQ is encountered
 5500   GO TO 10000                                                     
 5600   CALL DCYPHR(S,NI,1,AA,IDUM)                                     
        GO TO 100                                                       
 5700   CALL DCYPHR(S,NI,1,CFM,IDUM)                                    
        CALL DCYPHR(S,NI,1,TFM,IDUM)                                    
        GO TO 100                                                       
 5800   CALL DCYPHR(S,NI,1,STEP,IDUM)                                   
        CALL DCYPHR(S,NI,1,RMATCH,IDUM)                                 
        GO TO 100                                                       
 5900   CALL DCYPHR(S,NI,1,WIDTHK,IDUM)                                 
        GO TO 100                                                       
 5999   GO TO (6000,6100,6200,6300,6400,6500,6600,6700,6800,6900),JU    
 6000   CALL DCYPHR(S,NI,0,DUM,NOPT)                                    
        GO TO 100                                                       
 6100   CALL DCYPHR(S,NI,0,DUM,NMAX)                                    
        IF(NOPT.LT.1)  CALL DCYPHR(S,NI,1,ALEXP,IDUM)                   
C
C       PDJ We need to grab the parameter in the case of 1+al
        IF(NOPT.EQ.8)  CALL DCYPHR(S,NI,1,ALEXP,IDUM)                   
        GO TO 100                                                       
 6200   CALL DCYPHR(S,NI,0,DUM,L)                                       
        CALL DCYPHR(S,NI,1,PL(L+1),IDUM)                                
        IF(L.LT.0.OR.L.GT.19.OR.PL(L+1).LT.0.000)  WRITE(IW,6250)       
 6250   FORMAT(53H *** WARNING *** SPECIFICATIONS FOR INPUTED INITIAL L,
     1  60H-DISTRIBUTION WRONG *** L OUT OF RANGE OR POPULATION NEGATIV,
     2  8HE ******)                                                     
        GO TO 100                                                       
 6300   CALL DCYPHR(S,NI,0,DUM,NSTATE)                                  
        CALL DCYPHR(S,NI,0,DUM,KAPPA)                                   
        CALL DCYPHR(S,NI,1,EVACP,IDUM)                                  
        CALL DCYPHR(S,NI,1,EBIND,IDUM)                                  
        IF(NSTATE.GT.NMAX.OR.KAPPA.GT.2*NSTATE-1.OR.NSTATE.LE.0.OR.KAPPA
     1  .LE.0)  WRITE(IW,6325)                                          
 6325   FORMAT(53H *** ERROR *** DIRAC STATE INDECIES NEGATIVE, ZERO OR,
     1  60H OUT OF LIMITS *** CHECK THE Q.NUMBERS OF THE STATES INPUTED)
        IF(EVACP.LT.0.000)  EBIND = EBIND + EVACP                       
        IDR = IDR + 1                                                   
        IF(EVACP.LT.0.000)  WRITE(IW,6350)NSTATE,KAPPA                  
 6350   FORMAT(46H *** VACUUM POLARIZATION LESS THAN ZERO FOR N=,I2,    
     1  12H, AND KAPPA=,I2,4H ***)                                      
     0  ENERGY(NSTATE,KAPPA) = EBIND + EVACP                            
        IF(ABS(EBIND).LE.1.0E-20)  ENERGY(NSTATE,KAPPA) = EBIND + EVACP 
        GO TO 100                                                       
 6400   CALL DCYPHR(S,NI,0,DUM,IPRINT)                                  
        GO TO 100                                                       
 6500   DO 6550 K=1,20                                                  
        CALL DCYPHR(S,NI,0,DUM,NPOL(K))                                 
 6550   CONTINUE                                                        
        GO TO 100                                                       
 6600   CALL DCYPHR(S,NI,0,DUM,IPOL)                                    
        GO TO 100                                                       
 6700   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,1,ZMK,IDUM)                                    
        CALL DCYPHR(S,NI,1,ZML,IDUM)                                    
        CALL DCYPHR(S,NI,1,ZMM,IDUM)                                    
        GO TO 100                                                       
 6800   IF(IVRES.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,1,ZMKM,IDUM)                                   
        CALL DCYPHR(S,NI,1,ZMLM,IDUM)                                   
        CALL DCYPHR(S,NI,1,ZMMM,IDUM)                                   
        GO TO 100                                                       
 6900   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 6950 K=1,K1                                                  
        CALL DCYPHR(S,NI,0,DUM,IQ1(K))                                  
 6950   CONTINUE                                                        
        GO TO 100                                                       
 6999   GO TO(7000,7100,7200,7300,7400,7500,7600,7700,7800,7900),JU     
 7000   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 7050 K=1,K2                                                  
        CALL DCYPHR(S,NI,0,DUM,IQ2(K))                                  
 7050   CONTINUE                                                        
        GO TO 100                                                       
 7100   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 7150 K=1,K3                                                  
        CALL DCYPHR(S,NI,0,DUM,IQ3(K))                                  
 7150   CONTINUE                                                        
        GO TO 100                                                       
 7200   CALL DCYPHR(S,NI,0,DUM,IDB)                                     
        GO TO 100                                                       
 7300   DO 7350 K=1,4                                                   
        CALL DCYPHR(S,NI,1,YC(K),IDUM)                                  
 7350   CONTINUE                                                        
        GO TO 100                                                       
 7400   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        CALL DCYPHR(S,NI,0,DUM,IYC)                                     
        GO TO 100                                                       
 7500   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 7550 K=1,4                                                   
        CALL DCYPHR(S,NI,0,DUM,IJ(K))                                   
 7550   CONTINUE                                                        
        GO TO 100                                                       
 7600   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 7650 K=1,4                                                   
        CALL DCYPHR(S,NI,1,YJ(K),INT(DUM))                                   
 7650   CONTINUE                                                        
        GO TO 100                                                       
 7700   IF(IVERS.NE.1)  WRITE(IW,900)                                   
        IF(IVERS.NE.1)  GO TO 100                                       
        DO 7750 K=1,4                                                   
        CALL DCYPHR(S,NI,0,DUM,JJ1(K))                                  
 7750   CONTINUE                                                        
        GO TO 100                                                       
 7800   DO 7850 K=1,3                                                   
        CALL DCYPHR(S,NI,0,DUM,IPC(K))                                  
 7850   CONTINUE                                                        
        GO TO 100                                                       
 7900   CALL DCYPHR(S,NI,1,CL1,IDUM)                                    
        CALL DCYPHR(S,NI,1,CL2,IDUM)                                    
        GO TO 100                                                       
 7999   GO TO(8000,8100,8200,8300,8400,8500,8600,8700,8800,8900),JU     
 8000   CALL DCYPHR(S,NI,1,ESP,IDUM)                                    
        GO TO 100                                                       
 8100   MPU = MPU + 1                                                   
        CALL DCYPHR(S,NI,0,DUM,N1J)                                     
        CALL DCYPHR(S,NI,0,DUM,L1J)                                     
        CALL DCYPHR(S,NI,0,DUM,J1J)                                     
        CALL DCYPHR(S,NI,0,DUM,N2J)                                     
        CALL DCYPHR(S,NI,0,DUM,L2J)                                     
        CALL DCYPHR(S,NI,0,DUM,J2J)                                     
        CALL DCYPHR(S,NI,0,DUM,IRS)                                     
        IF(N1J.GT.20.OR.N1J.LE.1.OR.N2J.GT.N1J.OR.N2J.LE.0.OR.L1J.GE.N1J
     1  .OR.L1J.LT.0.OR.L2J.GE.N2J.OR.L2J.LT.0.OR.J1J*J1J.NE.J1J.OR.J2J 
     2  *J2J.NE.J2J.OR.IABS(L1J-L2J).GT.3.OR.IRS.GT.2.OR.IRS.LT.        
     3  0.OR.(N2J.EQ.N1J.AND.N2J.NE.2)) WRITE(IW,8150)                  
 8150   FORMAT(53H *** WARNING *** SPECIFICATIONS FOR TRANSITION TO BE ,
     1  60HPUNCHED ARE WRONG *** NO SUCH LINE OR GROUP EXISTS *** NO PU,
     2  8HNCH ****)                                                     
        ICPU(MPU) = N1J + 32*L1J + 1024*J1J + 2048*N2J + 65536*L2J +    
     1  2097152*J2J + 4194304*IRS                                       
        GO TO 100                                                       
 8200   CALL DCYPHR(S,NI,0,DUM,IDE)                                     
        GO TO 100                                                       
 8300   CALL DCYPHR(S,NI,0,DUM,IPN)                                     
        GO TO 100                                                       
 8400   DO 8450 K=1,6                                                   
        CALL DCYPHR(S,NI,0,DUM,JTM(K))                                  
 8450   CONTINUE                                                        
        GO TO 100                                                       
 8500   DO 8550 K=1,6                                                   
        CALL DCYPHR(S,NI,0,DUM,JTD(K))                                  
 8550   CONTINUE                                                        
        GO TO 100                                                       
 8600   DO 8650 K=1,6                                                   
        CALL DCYPHR(S,NI,0,DUM,JTQ(K))                                  
 8650   CONTINUE                                                        
        GO TO 100                                                       
 8700   DO 8750 K=1,6                                                   
        CALL DCYPHR(S,NI,0,DUM,JTO(K))                                  
 8750   CONTINUE                                                        
        GO TO 100                                                       
 8800   CALL DCYPHR(S,NI,0,DUM,N8)                                      
        CALL DCYPHR(S,NI,0,DUM,LD8)                                     
        CALL DCYPHR(S,NI,0,DUM,LU8)                                     
        CALL DCYPHR(S,NI,1,A8,IDUM)                                     
        L8 = LU8 - LD8 + 1                                              
        S8 = 0.000                                                      
        DO 8850 K=1,L8                                                  
        CALL DCYPHR(S,NI,1,B8,IDUM)                                     
        IF(N8.LT.1.OR.LD8.LT.0.OR.LU8.LT.LD8.OR.N8.GT.20.OR.LU8.GE.N8.OR
     1  .A8.LT.0.000.OR.B8.LT.0.000)  WRITE(IW,8880)                    
 8880   FORMAT(53H *** WARNING *** SPEC.S FOR THE INITIAL L-DISTRIBUTIO,
     1  60HN GIVEN ARE WRONG *** INTEGERS OUT OF LIMITS OR REALS NEGATI,
     2  8HVE *****)                                                     
        K8 = N8*(N8-1)/2 + LD8 + K                                      
        PLN(K8) = B8                                                    
        S8 = S8 + B8                                                    
 8850   CONTINUE                                                        
        DO 8875 K=1,L8                                                  
        K8 = N8*(N8-1)/2 + LD8 + K                                      
        IF(A8.LE.0.000)  A8 = S8                                        
        PLN(K8) = PLN(K8)*A8/AMAX1(S8,1.E-20)                           
 8875   CONTINUE                                                        
        GO TO 100                                                       
 8900   CALL DCYPHR(S,NI,0,DUM,IP8)                                     
        GO TO 100                                                       
 8999   GO TO( 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),JU     
C  ***  THIS SECTION IS DEVOTED TO THE DISCOVERY AND RECOVERY OF ERRORS.
C  ***  ALL POSSIBLE INPUT DATA ARE SCREENED FOR ERRORS (WITHIN REASON) 
C  ***  THIS PORTION MAY BE REMOVED IF YOU ARE CONFIDENT THAT YOU MAKE  
C  ***  NO MISTAKES IN THE INPUT SPECIFICATIONS......................   
C       We have arrived here when XEQ has been read in the file
C       so the cascade can begin
C       Here we do lots of error checking
10000   IF(Z.LE.0.00)  WRITE(IW,10100)                                  
        IF(Z.GT.99.00)  WRITE(IW,10200)                                 
        IF(Z.GT.137.04)  WRITE(IW,10300)                                
        IF(AMOD(Z+0.001,1.000).GT.                                      
     1  0.002)  WRITE(IW,10400)                                         
        IF(Z.LE.0.000)  STOP                                            
10100   FORMAT(53H *** ERROR *** ATOMIC NUMBER Z ZERO OR NEGATIVE *** E,
     1  60HXECUTION TERMINATED *** ....................................)
10200   FORMAT(53H *** WARNING *** ATOMIC NUMBER Z TOO BIG *** LAST TWO,
     1  60HDIGITS WILL BE PRINTED IN BLOCK LETTERS IN THE UPPER LEFT OF,
     2  8HTABLE **)                                                     
10300   FORMAT(53H *** WARNING *** ATOMIC NUMBER Z TOO BIG *** POINT LI,
     1  60HKE DIRAC FORMULAE HAVE PROBLEMS *** NO ATTEMPT TO RECTIFY PR,
     2  8HOBLEM **)                                                     
10400   FORMAT(53H *** WARNING *** ATOMIC NUMBER Z NOT CLOSE TO AN INTE,
     1  60HGER VALUE *** INTEGER PART PRINTED IN TABLE,BUT ACTUAL VALUE,
     2  8H USED **)                                                     
        DO 10500 K=1,3                                                  
        IF(BE(K).LT.0.000)  WRITE(IW,10600)                             
10500   CONTINUE                                                        
        IF(BE(1).LT.4.0*BE(2).OR.BE(2).LT.2.0*BE(3)) WRITE(             
     1  IW,10700)                                                       
        IF(BE(1).GT.15.0*Z*Z.OR.BE(1).LT.5.00*Z*Z)  WRITE(              
     2  IW,10800)                                                       
        IF(BE(3).GT.1000.0)  WRITE(IW,10900)                            
10600   FORMAT(53H *** ERROR *** BINDING ENERGY(IES) NEGATIVE OR ZERO *,
     1  60H** PROGRAM WILL HALT LATER *** NO ATTEMPT TO CORRECT PROBLEM,
     2  8H *******)                                                     
10700   FORMAT(53H *** WARNING *** BINDING ENERGIES NOT IN ANY REASONAB,
     1  60HLE PROPORTION *** POSSIBLY ENETRED OUT OF SEQUENCE (MUST BE ,
     2  8HK,L,M **)                                                     
10800   FORMAT(53H *** WARNING *** BINDING ENERGIES NOT IN ANY REASONAB,
     1  60HLE RANGE *** POSSIBLY IN THE WRONG UNITS (MUST BE IN EV) OR ,
     2  8HZ ******)                                                     
10900   FORMAT(53H *** WARNING *** BINDING ENERGY OF M SHELL UNREASONAB,
     1  60HLY HIGH *** CHECK YOUR SOURCE OR DISREGARD IF INTENTIONALLY ,
     2  8HSET ****)                                                     
        IF(K0.GT.3.OR.K1.GT.7.OR.K2.GT.7.OR.K3.GT.7)  WRITE(IW,11000)   
11000   FORMAT(53H *** ERROR *** TOO MANY CASES SPECIFIED FOR THE MULTI,
     1  60HPOLARITIES *** MAX ARE M/3 D,Q,O/7 *** RESULTS WILL BE INCOR,
     2  8HRECT ***)                                                     
        IF(K0.LT.0.OR.K1.LT.0.OR.K2.LT.0.OR.K3.LT.0)  WRITE(IW,11100)   
11100   FORMAT(53H *** ERROR *** NEGATIVE NUMBER OF CASES SPECIFIED FOR,
     1  60H THE MULTIPOLARITIES *** MUST BE POSITIVE OR ZERO TO SKIP AN,
     2  8HYONE ***)                                                     
        IF((K0.EQ.2.AND.NN0(1).EQ.NN0(2)).OR.(K0.EQ.3.AND.(NN0(1).EQ.NN0
     1  (2).OR.NN0(1).EQ.NN0(3).OR.NN0(2).EQ.NN0(3))))  WRITE(IW,11200) 
11200   FORMAT(53H *** ERROR *** DUPLICATE SHELL SPECIFICATION IN MONOP,
     1  60HOLE RATE CALCULATION *** WILL NOT ABORT BUT RESULTS WILL BE ,
     2  8HWRONG **)                                                     
        IF((K0.EQ.1.AND.NN0(1).EQ.0).OR.(K0.EQ.1.AND.NN0(1)*NN0(2).EQ.0)
     1  .OR.(K0.EQ.3.AND.NN0(1)*NN0(2)*NN0(3).EQ.0))  WRITE(IW,11300)   
11300   FORMAT(53H *** ERROR *** RADIATION SPECIFIED FOR MONOPOLE CASES,
     1  60H *** NO SUCH RATE EXISTS, SO THE PROGRAM MIGHT DO STRANGE TH,
     2  8HINGS ***)                                                     
        IF(K1.LE.1)  GO TO 12000                                        
        DO 11500 K=2,K1                                                 
        IF(NN1(K).LE.0)                                                 
     1  WRITE(IW,11600)                                                 
        IF(M1(K).LT.0)  WRITE(IW,11700)                                 
        IF(NN1(K).GT.3)  WRITE(IW,11900)                                
        IF(M1(K).GE.NN1(K).AND.NN1(K).NE.1)  WRITE(IW,11900)            
        IF(K.EQ.2)  GO TO 11500                                         
        K4 = K-1                                                        
        DO 11400 I=2,K4                                                 
        IF(NN1(I).EQ.NN1(K).AND.M1(K).EQ.M1(I))  WRITE(IW,11800)        
        IF(NN1(I).EQ.NN1(K).AND.M1(I)*M1(K).EQ.0)  WRITE(IW,11800)      
11400   CONTINUE                                                        
11500   CONTINUE                                                        
11600   FORMAT(53H *** ERROR *** NEGATIVE SHELL Q. N. SPECIFIED OR DUPL,
     1  60HICATE RADIATION CALCULATION *** NO ATTEMPT TO FIX *** CHECK ,
     2  8HNNJ ****)                                                     
11700   FORMAT(53H *** ERROR *** NEGATIVE SUBSHELL CODE SPECIFICATION *,
     1  60H** PROGRAM WILL TAKE AN UNPREDICTABLE BRANCH OR ABORT *** CH,
     2  8HECK MJ *)                                                     
11800   FORMAT(53H *** ERROR *** DUPLICATE OR OVERLAPPING SUBSHELL CODE,
     1  60H SPECIFIED *** WILL NOT ABORT, BUT RATES WILL BE DONE TWICE ,
     2  8H********)                                                     
11900   FORMAT(53H *** ERROR *** NON EXISTENT SHELL (.GT.3) OR SUBSHELL,
     1  60H COMBINATION (E.G. 2D) *** OUTCOME UNPREDICTABLE *** CHECK N,
     2  8HNJ,MJ **)                                                     
12000   IF(K2.LE.1)  GO TO 12300                                        
        DO 12200 K=2,K2                                                 
        IF(NN2(K).LE.0)                                                 
     1  WRITE(IW,11600)                                                 
        IF(M2(K).LT.0)  WRITE(IW,11700)                                 
        IF(NN2(K).GT.3)  WRITE(IW,11900)                                
        IF(M2(K).GE.NN2(K).AND.NN2(K).NE.1)  WRITE(IW,11900)            
        IF(K.EQ.2)  GO TO 12200                                         
        K4 = K-1                                                        
        DO 12100 I=2,K4                                                 
        IF(NN2(I).EQ.NN2(K).AND.M2(K).EQ.M2(I))  WRITE(IW,11800)        
        IF(NN2(I).EQ.NN2(K).AND.M2(I)*M2(K).EQ.0)  WRITE(IW,11800)      
12100   CONTINUE                                                        
12200   CONTINUE                                                        
12300   IF(K3.LE.1)  GO TO 12600                                        
        DO 12500 K=2,K3                                                 
        IF(NN3(K).LE.0)                                                 
     1  WRITE(IW,11600)                                                 
        IF(M3(K).LT.0)  WRITE(IW,11700)                                 
        IF(NN3(K).GT.3)  WRITE(IW,11900)                                
        IF(M3(K).GE.NN3(K).AND.NN3(K).NE.1)  WRITE(IW,11900)            
        IF(K.EQ.2)  GO TO 12500                                         
        K4=K-1                                                          
        DO 12400 I=2,K4                                                 
        IF(NN3(I).EQ.NN3(K).AND.M3(K).EQ.M3(I))  WRITE(IW,11800)        
        IF(NN3(I).EQ.NN3(K).AND.M3(I)*M3(K).EQ.0)  WRITE(IW,11800)      
12400   CONTINUE                                                        
12500   CONTINUE                                                        
12600   IF(K1.LE.1)  GO TO 12900                                        
        DO 12700 K=2,K1                                                 
        IF(IP1(K)*IP1(K).NE.IP1(K))  WRITE(IW,12800)                    
12700   CONTINUE                                                        
12800   FORMAT(53H *** WARNING *** PENETRATION CODES FOR SUBSHELLS NOT ,
     1  60HZERO OR ONE *** POSSIBLE ERRONEOUS RESULT *** CHECK ARRAYS I,
     2  8HPJ(K) **)                                                     
12900   IF(K2.LE.1)  GO TO 13100                                        
        DO 13000 K=2,K2                                                 
        IF(IP2(K)*IP2(K).NE.IP2(K))  WRITE(IW,12800)                    
13000   CONTINUE                                                        
13100   IF(K3.LE.1)  GO TO 13300                                        
        DO 13200 K=2,K3                                                 
        IF(IP3(K)*IP3(K).NE.IP3(K))  WRITE(IW,12800)                    
13200   CONTINUE                                                        
13300   DO 13400 I=1,6                                                  
        IF(JTM(I).GT.3.OR.JTM(I).LT.1.OR.JTD(I).LT.1                    
     1  .OR.JTD(I).GT.3.OR.JTQ(I).LT.1.OR.JTQ(I).GT.3.OR.JTO(I).GT.3.OR.
     2  JTO(I).LT.1)  WRITE(IW,13500)                                   
13400   CONTINUE                                                        
13500   FORMAT(53H *** WARNING *** NUMBER OF TERMS IN PENETRATION CALCU,
     1  60HLATION OUTSIDE RANGE (1-3) *** DISREGARD IF PENETRATION NOT ,
     2  8HUSED ***)                                                     
        IF(IC.LT.0.OR.IC.GT.3)  WRITE(IW,13600)                         
13600   FORMAT(53H *** ERROR *** ACCURACY CONTROL OPTION CODE OUTSIDE P,
     1  60HERMISSIBLE RANGE (0-3) *** UNPREDICTABLE RESULTS CAN HAPPEN ,
     2  8H********)                                                     
        IF(IR.LE.0.OR.IW.LE.0.OR.IPUNCH.LE.0)  WRITE(IW,13700)          
13700   FORMAT(53H *** WARNING *** I/O UNIT NUMBER NEGATIVE OR ZERO ***,
     1  60H UNLIKELY TO BE READ OR WRITE UNIT *** IF PUNCH...WHO KNOWS ,
     2  8HRESULT *)                                                     
        IF(IR.GT.99.OR.IW.GT.99.OR.IPUNCH.GT.99)  WRITE(IW,13800)       
13800   FORMAT(53H *** WARNING *** I/O UNIT NUMBER EXCEEDING 99 *** NON,
     1  60H STANDARD FORTRAN ASSIGNMENT *** DISREGARD IF INTENTIONALLY ,
     2  8HSET ****)                                                     
        IF(IPR.LT.0)  WRITE(IW,13900)                                   
        IF(IPR.GT.63)  WRITE(IW,14000)                                  
13900   FORMAT(53H *** ERROR *** PRINT SELECTION OPTION CODE NEGATIVE *,
     1  60H** MODULO ROUTINE WILL FIGURE OPTIONS ERRONEOUSLY *** CHECK ,
     2  8HIPR ****)                                                     
14000   FORMAT(53H *** WARNING *** PRINT SELECTION OPTION CODE .GT. 63 ,
     1  60H*** LAST 6 BITS OF NUMBER WILL BE USED IN PRINT (MOD(IPR,64),
     2  8H) ******)                                                     
        IF(FD.GT.99.99.OR.FD.LT.0.01)  WRITE(IW,14100)                  
14100   FORMAT(53H *** WARNING *** FACTORIAL DIVIDER NOT IN ANY REASONA,
     1  60HBLE RANGE (0.01-99.99) *** COULD CAUSE SEVERE ARITHMETIC PRO,
     2  8HBLEMS **)                                                     
        IF(IDB*IDB.NE.IDB)  WRITE(IW,14200)                             
14200   FORMAT(53H *** WARNING *** DEBUG OPTION SELECTION SWITCH NOT ZE,
     1  60HRO OR ONE *** COULD CAUSE ERRORS IN THE PRINTING OF DETAILED,
     2  8H RATES *)                                                     
        IF(IPN*IPN.NE.IPN)  WRITE(IW,14300)                             
14300   FORMAT(53H *** WARNING *** PUNCH SELECTION SWITCH NOT ZERO OR O,
     1  60HNE *** COULD RESULT IN UNINTENTIONAL INCLUSION OR OMISSION O,
     2  8HF PUNCH*)                                                     
        IF(IDE.LT.0.OR.IDE.GT.99999)  WRITE(IW,14400)                   
14400   FORMAT(53H *** WARNING *** PUNCH CARD IDENTIFICATION NUMBER NEG,
     1  60HATIVE OT .GT. 99999 *** WILL PUNCH 5 STARS INSTEAD, IF .GT. ,
     2  8H5 DIGITS)                                                     
        IF(WIDTHK.LT.0.0.OR.WIDTHK.GT.999.999)  WRITE(IW,14500)         
14500   FORMAT(53H *** WARNING *** REFILLING WIDTH OF K-ELECTRON SHELL ,
     1  60HNEGATIVE OR UNREASONABLY LARGE *** CHECK FOR PROPER UNITS (E,
     2  8HV) *****)                                                     
        IF(D2P1S.LT.Z*Z*AMASSM.OR.D2P1S.GT.10.2*Z*Z*AMASSM)             
     1  WRITE(IW,14600)                                                 
C       This is the current problem being encountered
14600   FORMAT(53H *** WARNING *** ENERGY OF THE 2P-1S MUONIC TRANSITIO,
     1  60HN NEGATIVE, TOO LOW OR UNREASONABLY HIGH *** CHECK FOR UNITS,
     2  8H (EV) **)                                                     
        IF(ESP.LT.0.000.OR.ESP.GT.1.0E7)  WRITE(IW,14700)               
14700   FORMAT(53H *** WARNING *** ENERGY OF THE 2S-2P MUONIC TRANSITIO,
     1  60HN NEGATIVE OR TOO LARGE *** SET TO ZERO IF TRANSITION TO BE ,
     2  8HSKIPPED*)                                                     
        IF(NOPT.LT.-1.OR.NOPT.GT.2)  WRITE(IW,14800)                    
14800   FORMAT(53H *** WARNING *** INITIAL L-DISTRIBUTION OPTION CODE U,
     1  60HNRECOGNIZABLE *** COULD CAUSE UNEXPECTED COMPLICATIONS IF OF,
     2  8H LIMITS*)                                                     
        IF(NMAX.LT.2.OR.NMAX.GT.20)  WRITE(IW,14900)                    
14900   FORMAT(53H *** ERROR *** STARTING N QUANTUM NUMBER OF THE CASCA,
     1  60HDE NOT IN THE RANGE 2-20 *** PROGRAM WILL ABORT IF .GT. 20 O,
     2  8HR .LT. 1)                                                     
        IF(NOPT.EQ.0.AND.ABS(ALEXP).GT.1.000)  WRITE(IW,15000)          
15000   FORMAT(53H *** WARNING *** MODIFIED STATISTICAL L-DISTRIBUTION ,
     1  60HEXPONENT TOO HIGH OR TOO LOW (NEG) *** COULD CAUSE ARITH. OV,
     2  8HERFLOW *)                                                     
        IF(IP8*IP8.NE.IP8)  WRITE(IW,15100)                             
15100   FORMAT(53H *** WARNING *** L-DISTRIBUTION TABLE SELECTION CODE ,
     1  60H(TOP N ONLY OR FULL N-L) NOT ZERO OR ONE *** COULD RESULT IN,
     2  8H ERRORS*)                                                     
        IF(ABS(CL1).GT.10.0.OR.ABS(CL2).GT.10.0)  WRITE(IW,15200)       
15200   FORMAT(53H *** WARNING *** INITIAL QUADRATIC L-DISTRIBUTION PAR,
     1  60HAMETERS UNREASONABLY HIGH OR LOW (ABS .GT. 10.0) *** POSSIBL,
     2  8HE ERRORS)                                                     
        IF(IPC(1)*IPC(1).NE.IPC(1).OR.IPC(2)*IPC(2).NE.IPC(2).OR.IPC(3)*
     1  IPC(3).NE.IPC(3))  WRITE(IW,15300)                              
15300   FORMAT(53H *** WARNING *** ELECTRON REFILLING CONTROL CODES NOT,
     1  60H EQUAL TO ZERO OR ONE *** REFILLING MIGHT NOT BE DONE PROPER,
     2  8HLY *****)                                                     
        DO 15400 I=1,NMAX                                               
        IF(NPOL(I).LT.-1.OR.NPOL(I).GT.NMAX-1)                          
     1  WRITE(IW,15500)                                                 
15400   CONTINUE                                                        
        IF(IPOL*IPOL.NE.IPOL)  WRITE(IW,15600)                          
15500   FORMAT(53H *** WARNING *** POLARIZATION CODE N.S FOR EACH L OUT,
     1  60H OF RANGE *** POLARIZATION MIGHT BE WRONG OR PROGRAM WILL AB,
     2  8HORT ****)                                                     
15600   FORMAT(53H *** WARNING *** POLARIZATION CALCULATION SELECTION S,
     1  60HWITCH NOT EQUAL TO ZERO OR ONE *** POSSIBLE UNDESIRED RESULT,
     2  8HS ******)                                                     
        IF(YC(1).LT.0.0.OR.YC(1).GT.20.0.OR.YC(2).LT.0.0.OR.YC(2).GT.20.
     1  0.OR.YC(3).GT.20.0.OR.YC(3).LT.0.0.OR.YC(4).LT.0.0.OR.YC(4).GT. 
     2  20.)  WRITE(IW,15700)                                           
15700   FORMAT(53H *** WARNING *** CUTOFF Y.S FOR THE MULTIPOLARITIES N,
     1  60HEGATIVE OR UNREASONABLY HIGH *** IF HIGH CHECK IF SO DESIRED,
     2  8H *******)                                                     
        IF(ABS(AMASSM-206.7686).GT.1.0E-10)  WRITE(IW,15800)            
15800   FORMAT(53H *** WARNING *** NONSTANDARD MASS OF PARTICLE (NOT MU,
     1  60HON) *** CHECK IF THIS IS INTENTIONAL AND THE MASS IS IN ELEC,
     2  8HT. M.S *)                                                     
        IF(ABS(AMASSE-511003.4).GT.1.0E-03)  WRITE(IW,15900)            
15900   FORMAT(53H *** WARNING *** NONSTANDARD MASS FOR THE ELECTRON **,
     1  60H* CHECK IF THIS IS INTENTIONAL AND THAT THE MASS IN IN ELECT,
     2  8H. VOLTS*)                                                     
        IF(ABS(AMASSN-931.48).GT.1.0E-10)  WRITE(IW,16000)              
16000   FORMAT(53H *** WARNING NONSTANDARD MASS FOR THE AVERAGE NUCLEON,
     1  60H BOUND MASS *** CHECK IF THIS IS INTENTIONAL AND THE MASS IN,
     2  8H MEV ***)                                                     
        IF(ABS(AA-140.0).GT.1.0E-10.AND.(AA.LT.1.3*Z.OR.AA.GT.2.7*Z))   
     1  WRITE(IW,16100)                                                 
16100   FORMAT(53H *** WARNING *** ATOMIC WEIGHT A NOT CHANGED FROM DEF,
     1  60HAULT OR UNREASONABLY HIGH OR LOW *** REDUCED MASS CALCULATIO,
     2  8HN OFF **)                                                     
        IF(CFM.LT.0.0.OR.CFM.GT.7.5.OR.TFM.LT.0.5.OR.TFM.GT.3.0)        
     1  WRITE(IW,16200)                                                 
16200   FORMAT(53H *** WARNING *** FERMI DISTRIBUTION PARAMETERS UNREAS,
     1  60HONABLY HIGH OR LOW *** NOT USED IN ANY CALCULATION IN THIS P,
     2  8HROGRAM *)                                                     
        IF(STEP.LT.0..OR.RMATCH.LT.0..OR.RMATCH.GT.1.E4) WRITE(IW,16300)
16300   FORMAT(53H *** WARNING *** STEP IN INTEGRATION NEGATIVE OR MATC,
     1  60HHING RADIUS NEGATIVE OR UNREASONABLY LARGE *** NOT USED IN P,
     2  8HROGRAM *)                                                     
        IF(EHIGH.GT.30.0.OR.EHIGH.LT.0.001)  WRITE(IW,16400)            
16400   FORMAT(53H *** WARNING *** HIGH CUT OF ENERGY IN X-RAY CATALOGU,
     1  60HE TOO HIGH OR TOO LOW *** CHECK UNITS (MEV) *** OK IF INTENT,
     2  8HIONAL **)                                                     
        IF(ELOW.GT.EHIGH.OR.ELOW.LT.0.001)  WRITE(IW,16500)             
16500   FORMAT(53H *** WARNING *** LOW CUT OF ENERGY IN X-RAY CATALOGUE,
     1  60H MORE THAN THE HIGH CUT OR TOO LOW *** CHECK UNITS (MEV) AND,
     2  8H EHI ***)                                                     
        IF(CLIMIT.GT.0.5.OR.CLIMIT.LT.1.0E-7)  WRITE(IW,16600)          
16600   FORMAT(53H *** WARNING *** INTENSITY LIMIT CUTOFF FOR THE X-RAY,
     1  60H CATALOGUE TOO HIGH OR TOO LOW *** TOO FOW OR TOO MANY LINES,
     2  8HWRITTEN*)                                                     
        IF(ERES.LT.1.0E-6.OR.ERES.GT.0.05)  WRITE(IW,16700)             
16700   FORMAT(53H *** WARNING *** ENERGY RESOLUTION IN THE X-RAY CATAL,
     1  60HOGUE TOO HIGH OR TOO LOW *** POSSIBLY WRONG UNITS (MUST BE I,
     2  8HN MEV **)                                                     
        IF(ICC*ICC.NE.ICC)  WRITE(IW,16800)                             
16800   FORMAT(53H *** WARNING *** INPUTED DIVIDING POINTS SWITCH IS NO,
     1  60HT ZERO OR ONE *** POSSIBLE UNWANTED DIVIDING POINTS IN CATAL,
     2  8HOGUE ***)                                                     
        IF(CD(1).GT.0.5.OR.CD(2).LE.CD(3).OR.CD(3).LE.CD(4).OR.CD(4).LT.
     1  CD(5).OR.CD(5).LT.1.2*CLIMIT)  WRITE(IW,16900)                  
16900   FORMAT(53H *** WARNING *** NEW STAR DIVIDING POINTS FOR THE X-R,
     1  60HAY CATALOGUE TOO HIGH, TOO LOW OR UNREASONABLY CLOSE SPACED ,
     2  8H********)                                                     
        IF(EB.LE.0.000)  WRITE(IW,17000)                                
17000   FORMAT(53H *** WARNING *** CALIBRATION PARAMETER B FOR THE CONV,
     1  60HERSION OF ENERGY TO CHANNEL NUMBER NEGATIVE OR 0 *** IF 0 PR,
     2  8HOG. HALT)                                                     
        IF(MPU.GT.200)  WRITE(IW,17100)                                 
17100   FORMAT(53H *** ERROR *** TOO MANY LINES SPECIFIED TO BE PUNCHED,
     1  60H *** IF YOU WANT MORE THAN 200 LINES INCREASE THE DIM OF ICP,
     2  8HU IN L41)                                                     
        IF(IPC(1)*IPC(2)*IPC(3).EQ.0.AND.IP8.EQ.1)  WRITE(IW,17200)     
        IF(IPC(1)*IPC(2)*IPC(3).EQ.0.AND.IP8.EQ.1)  STOP                
17200   FORMAT(53H *** ERROR *** FULL (N,L) L-DISTRIBUTION REQUIRES REF,
     1  60HILLING OF ALL SHELLS (IPC 1 1 1) *** EXECUTION TERMINATED **,
     2  8H********)                                                     
        IF(IPOL.EQ.0.AND.IP8.EQ.1)  WRITE(IW,17300)                     
17300   FORMAT(53H *** WARNING *** FULL N-L DISTRIBUTION SPECIFIED AND ,
     1  60HDEPOLARIZATION CALCULATION *** DEPOLARIZATION MAY BE WRONG *,
     2  8H********)                                                     
        IF(EA.GT.2.0E4.OR.EA.LT.-2.0E4)  WRITE(IW,17400)                
17400   FORMAT(53H *** WARNING *** CALIBRATION ENERGY POINT EA UNREASON,
     1  60HABLY HIGH OR LOW *** CHECK FOR UNITS (MUST BE IN KEV) ******,
     2  8H********)                                                     
        IF(Z.LT.30.0.AND.(Z-20.0)*0.1.GT.POP(6))  WRITE(IW,17500)       
        IF(Z.LT.18.0.AND.(Z-12.0)/6.0.GT.POP(5))  WRITE(IW,17500)       
        IF(Z.LT.12.0.AND.(Z-10.0)*0.5.GT.POP(4))  WRITE(IW,17500)       
        IF(Z.LT.10.0.AND.(Z-4.0)/6.0.GT.POP(3))  WRITE(IW,17500)        
        IF(Z.LT.4.0.AND.(Z-2.0)*0.5.GT.POP(2))  WRITE(IW,17500)         
        IF((Z.LE.20.0.AND.POP(6).NE.0.0).OR.(Z.LE.12.0.AND.POP(5).NE.0.0
     1  ).OR.(Z.LE.10.0.AND.POP(4).NE.0.0).OR.(Z.LE.4.0.AND.POP(3).NE.0.
     2  ).OR.(Z.LE.2.0.AND.POP(2).NE.0.0))  WRITE(IW,17500)             
17500   FORMAT(53H *** ERROR *** Z IS TOO LOW TO SUPPORT SPECIFIED POPU,
     1  60HLATION OF ELECTRONIC SHELLS *** POSSIBLE ERRONEOUS RESULTS *,
     1  8H********)                                                     
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE DCYPHR(A,NJ,ITYPE,R,I)                               
C  ***  DECYPHERS THE NUMBERS FOR ROUTINE CREAD. FINDS INPUT ERRORS     
C       A is the string from the input file, NJ is number, ITYPE ?
C       R is the return value, and I is an integer
        INTEGER A,BB,AA,BL,SE,CO,PL,DA,PO,EE                            
        DIMENSION A(70),BB(10),AA(20)                                   
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        DATA BB/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/                
        DATA BL,SE,CO,PL,DA,PO,EE/1H ,1H/,1H,,1H+,1H-,1H.,1HE/          
        IG = 0                                                          
        DO 100 J=NJ,70                                                  
        NI1 = J+1                                                       
        IF((A(J).EQ.BL.OR.A(J).EQ.SE.OR.A(J).EQ.CO).AND.IG.EQ.0)        
     1  GO TO 100                                                       
        IF((A(J).EQ.BL.OR.A(J).EQ.SE.OR.A(J).EQ.CO).AND.IG.NE.0)        
     1  GO TO 200                                                       
        IG = IG+1                                                       
        AA(IG) = A(J)                                                   
  100   CONTINUE                                                        
  200   NJ = NI1                                                        
        IF(ITYPE.EQ.1)  GO TO 1000                                      
        DO 600 J=1,IG                                                   
        IFL = 0                                                         
        IF(J.NE.1)  GO TO 300                                           
        IF(AA(J).EQ.PL.OR.AA(J).EQ.DA)  IFL=1                           
  300   DO 400 K=1,10                                                   
        IF(AA(J).EQ.BB(K))  IFL=1                                       
  400   CONTINUE                                                        
        IF(IFL.EQ.0)  WRITE(IW,500)AA(J)                                
  500   FORMAT(/38H *** ERROR *** ILLEGAL DATA IN INPUT /,A1,5H/ ***/   
     1  60H *** STANDARD FIXUP TAKEN (TAKEN AS 0), EXECUTION CONTINUING)
     0  IF(IFL.EQ.0)  AA(J)=BB(1)                                       
  600   CONTINUE                                                        
        IS = 1                                                          
        IF(AA(1).EQ.PL.OR.AA(1).EQ.DA)  IS=2                            
        N=0                                                             
        DO 700 J=IS,IG                                                  
        K = IG-J                                                        
        N = N+ID(AA(J))*10**K                                           
  700   CONTINUE                                                        
        I = N                                                           
        IF(AA(1).EQ.DA)  I=-I                                           
        RETURN                                                          
 1000   IE = 0                                                          
        DO 1100 J=1,IG                                                  
        IF(AA(J).EQ.EE)  IE=1                                           
 1100   CONTINUE                                                        
        IF(IE.EQ.1)  GO TO 2000                                         
        IH = IG+1                                                       
        IFL = 0                                                         
        DO 1200 J=1,IG                                                  
        IF(AA(J).EQ.PO)  IH = J                                         
        IF(AA(J).EQ.PO)  IFL=IFL+1                                      
 1200   CONTINUE                                                        
        IF(IFL.EQ.1)  GO TO 1500                                        
        IF(IFL.EQ.0)  WRITE(IW,1300)                                    
 1300   FORMAT(/48H *** WARNING *** NO DECIMAL POINT IN REAL NUMBER,    
     1  39H *** ASSUMED TO BE AT THE RIGHT END ***)                     
     0  IF(IFL.GT.1)  WRITE(IW,1400)                                    
 1400   FORMAT(/53H *** ERROR *** TOO MANY DECIMAL POINTS IN REAL NUMBER
     1  ,57H *** LAST ENCOUNTERED ASSUMED, OTHERS CHANGED TO ZERO ***)  
 1500  DO 1800 J=1,IG                                                  
        IF(J.EQ.IH)  GO TO 1800                                         
        IFL = 0                                                         
        IF(J.NE.1)  GO TO 1600                                          
        IF(AA(J).EQ.PL.OR.AA(J).EQ.DA)  IFL=1                           
 1600   DO 1700 K=1,10                                                  
        IF(AA(J).EQ.BB(K))  IFL=1                                       
 1700   CONTINUE                                                        
        IF(IFL.EQ.0)  WRITE(IW,500)AA(J)                                
        IF(IFL.EQ.0)  AA(J) = BB(1)                                     
 1800   CONTINUE                                                        
        IS = 1                                                          
        IF(AA(1).EQ.PL.OR.AA(1).EQ.DA)  IS=2                            
        RE = 0.000                                                      
        DO 1900 J=IS,IG                                                 
        IF(J.EQ.IH)  GO TO 1900                                         
        K = IH-J-1                                                      
        IF(K.LT.0)  K=K+1                                               
        RE = RE + FLOAT(ID(AA(J)))*10.000**K                            
 1900   CONTINUE                                                        
        R = RE                                                          
        IF(AA(1).EQ.DA)  R=-R                                           
        RETURN                                                          
 2000   IH = 0                                                          
        IFL = 0                                                         
        DO 2100 J=1,IG                                                  
        IF(AA(J).EQ.PO)  IH=J                                           
        IF(AA(J).EQ.PO)  IFL=IFL+1                                      
        IF(AA(J).EQ.EE)  IE=J                                           
 2100   CONTINUE                                                        
        IG1 = IE-1                                                      
        IF(IFL.EQ.0)  IH=IG1+1                                          
        IF(IFL.EQ.0)  WRITE(IW,1300)                                    
        IF(IFL.GT.1)  WRITE(IW,1400)                                    
        IF(IE.EQ.IG)  AA(IG+1)=BB(1)                                    
        IF(IE.EQ.IG)  IG=IG+1                                           
        DO 2400 J=1,IG1                                                 
        IF(J.EQ.IH)  GO TO 2400                                         
        IFL = 0                                                         
        IF(J.NE.1)  GO TO 2200                                          
        IF(AA(J).EQ.PL.OR.AA(J).EQ.DA)  IFL=1                           
 2200   DO 2300 K=1,10                                                  
        IF(AA(J).EQ.BB(K))  IFL=1                                       
 2300   CONTINUE                                                        
        IF(IFL.EQ.0)  WRITE(IW,500)AA(J)                                
        IF(IFL.EQ.0)  AA(J) = BB(1)                                     
 2400   CONTINUE                                                        
        IS = 1                                                          
        IF(AA(1).EQ.PL.OR.AA(1).EQ.DA)  IS=2                            
        RE = 0.000                                                      
        DO 2500 J=IS,IG1                                                
        IF(J.EQ.IH)  GO TO 2500                                         
        K = IH-J-1                                                      
        IF(K.LT.0)  K=K+1                                               
        RE = RE + FLOAT(ID(AA(J)))*10.000**K                            
 2500   CONTINUE                                                        
        IF(AA(1).EQ.DA)  RE=-RE                                         
        IS = IE+1                                                       
        IF(AA(IE+1).EQ.PL.OR.AA(IE+1).EQ.DA)  IS=IE+2                   
        N = 0                                                           
        DO 2600 J=IS,IG                                                 
        K = IG-J                                                        
        N = N + ID(AA(J))*10**K                                         
 2600   CONTINUE                                                        
        IF(AA(IE+1).EQ.DA)  N=-N                                        
        R = RE*10.000**N                                                
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE FFIX                                                 
        DOUBLE PRECISION F,DZA,DZA2,DREDM                               
C  ***  INITIALIZES VARIABLES THAT CANNOT BE SIMPLY ASSIGNED IN DATA    
        COMMON/LOC001/IJK,ENERGY,ECONS,ECONST,D2P1SM,D2P1S              
        COMMON/LOC002/BEM(3),ZSA(3),BE(3)                               
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC009/F(60),FD                                          
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC013/COEMON(30),EXPMON(30)                             
        COMMON/LOC014/COEDIP(42),EXPDIP(42)                             
        COMMON/LOC015/COEQUA(45),EXPQUA(45)                             
        COMMON/LOC016/COEOCT(45),EXPOCT(45)                             
        COMMON/LOC028/M1(7),M2(7),M3(7),YC(4),IDB                       
        COMMON/LOC030/POP(6),JTM(6),JTD(6),JTQ(6),JTO(6)                
        COMMON/LOC031/JM(10),JD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC034/DZA,DZA2,DREDM                                    
        COMMON/LOC036/ZMK,ZML,ZMM,ZMKM,ZMLM,ZMMM,IVERS                  
        COMMON/LOC038/A,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)               
        COMMON/LOC040/AMASSA,AMASSN,HBAR                                
        DATA AMASSS/206.7686/                                           
C
C       FD is the factorial divider, which is 15 by default
        F(1)=1.000D00                                                   
        DO 100 I=2,60                                                   
C       Here we are doing F(i) = ((F(i- 1)*i-1)/15
        F(I) = F(I-1)*DBLE(FLOAT(I-1))/DBLE(FD)                         
  100   CONTINUE                                                        
        IF(Z.LT.1.0E-20)  WRITE(IW,150)                                 
  150   FORMAT(52H0*** ERROR *** Z NOT GIVEN, EXECUTION TERMINATED ***) 
     0  IF(Z.LT.1.0E-20)  STOP                                          
        IF(IVERS.NE.1)  GO TO 160                                       
        IF(Z-ZSK.LT.0.000.OR.Z-ZSK.GT.ZMKM)  ZSK=Z-ZMK                  
        IF(Z-ZSL.LT.0.000.OR.Z-ZSL.GT.ZMLM)  ZSL=Z-ZML                  
        IF(Z-ZSM.LT.0.000.OR.Z-ZSM.GT.ZMMM)  ZSM=Z-ZMM                  
  160   IF(Z-ZSK.LT.0.000.OR.Z-ZSK.GT.ZMKM)  WRITE(IW,170)              
        IF(Z-ZSL.LT.0.000.OR.Z-ZSL.GT.ZMLM)  WRITE(IW,180)              
        IF(Z-ZSM.LT.0.000.OR.Z-ZSM.GT.ZMMM)  WRITE(IW,190)              
  170   FORMAT(53H *** WARNING *** EFFECTIVE CHARGE FOR K SHELL TOO HIG,
     1  36HH OR TOO LOW *** NO ACTION TAKEN ***)                        
  180   FORMAT(53H *** WARNING *** EFFECTIVE CHARGE FOR L SHELL TOO HIG,
     1  36HH OR TOO LOW *** NO ACTION TAKEN ***)                        
  190   FORMAT(53H *** WARNING *** EFFECTIVE CHARGE FOR M SHELL TOO HIG,
     1  36HH OR TOO LOW *** NO ACTION TAKEN ***)                        
        IF(ZSK.LT.1.000)  ZSK=1.000                                     
        IF(ZSL.LT.1.000)  ZSL=1.000                                     
        IF(ZSM.LT.1.000)  ZSM=1.000                                     
        ZSKZ = ZSK/Z                                                    
        ZSLZ = ZSL/Z                                                    
        ZSMZ = ZSM/Z                                                    
        AMZZ(1) = ZSKZ/AMASSM                                           
        AMZZ(2) = ZSLZ/AMASSM                                           
        AMZZ(3) = ZSMZ/AMASSM                                           
        ECONST = ECONS*Z*Z                                              
        DO 200 I=1,3                                                    
        BEM(I) = BE(I)/AMASSE                                           
  200   CONTINUE                                                        
        BM = BEM(1)*BEM(2)*BEM(3)                                       
        IF(BM.LT.1.0E-20)  WRITE(IW,250)BE                              
  250   FORMAT(53H0*** ERROR *** UNDEFINED OR ZERO BINDING ENERGIES ***,
     1  3F12.3,30H  *** EXECUTION TERMINATED ***)                       
     0  IF(BM.LT.1.0E-20)  STOP                                         
        ZSA(1) = ZSK*ALFA                                               
        ZSA(2) = ZSL*ALFA                                               
        ZSA(3) = ZSM*ALFA                                               
        D2P1SM = D2P1S/AMASSE                                           
        ESPM = ESP/AMASSE                                               
        IF(ABS(AMASSS-AMASSM).LT.1.0E-20)  GO TO 600                    
        DO 300 I=1,30                                                   
        EXPMON(I) = EXPMON(I)/AMASSM*AMASSS                             
  300   CONTINUE                                                        
        DO 400 I=1,42                                                   
        EXPDIP(I) = EXPDIP(I)/AMASSM*AMASSS                             
  400   CONTINUE                                                        
        DO 500 I=1,45                                                   
        EXPQUA(I) = EXPQUA(I)/AMASSM*AMASSS                             
        EXPOCT(I) = EXPOCT(I)/AMASSM*AMASSS                             
  500   CONTINUE                                                        
        AMASSS = AMASSM                                                 
C       This calculates Z * fine structure constant
  600   DZA = DBLE(Z*ALFA)                                              
        DZA2 = DZA**2                                                   
        AME = AMASSE*1.000E-06                                          
C
C       Nuclear mass * muon mass / nuclear mass + electron mass
        DREDM = DBLE(A*AMASSN*AMASSM*AME/(A*AMASSN + AMASSM*AME))       
        AMASST = AMASSA*AMASSN                                          
        IF(AMASSA.GT.1.0E-20)  DREDM=DBLE(AMASST*AMASSM*AME/(AMASST +   
     1  AMASSM*AME))                                                    
        IF(IYC.EQ.0)  GO TO 700                                         
        YA = 0.0297*Z**0.666667                                         
        YB = 0.0667*SQRT(Z)                                             
        YK = 0.0758*SQRT(Z)                                             
        YD = 0.0850*SQRT(Z)                                             
        YC(1) = AMIN1(YC(1),YA)                                         
        YC(2) = AMIN1(YC(2),YB)                                         
        YC(3) = AMIN1(YC(3),YK)                                         
        YC(4) = AMIN1(YC(4),YD)                                         
  700   CONTINUE                                                        
        IF(CFM.LT.1.0E-20)  CFM = 1.100*A**0.333333                     
        IF(AMASSA.GT.1.0E-20.AND.ABS(CFM-1.100*A**0.3333).LT.1.0E-20)   
     1  CFM = 1.100*AMASSA**0.333333                                    
        IF(D2P1S.GT.1.0E-20)  GO TO 800                                 
        R1 = 1.200*A**0.3333                                            
        X = 2.000E-05*Z*R1*AMASSM/0.529                                 
        D2P1SM = ECONST*(0.750 + 3.0/X**3*(X*X-4.0-X*X*X/3.0            
     1  + EXP(-X)*(X*X+4.0+4.0*X)))                                     
  800   CONTINUE                                                        
        IF(IVERS.NE.1)  GO TO 900                                       
        JM(1)=JTM(1)                                                    
        JD(1)=JTD(1)                                                    
        JQ(1)=JTQ(1)                                                    
        JO(1)=JTO(1)                                                    
        JM(2)=JTM(2)                                                    
        JD(2)=JTD(2)                                                    
        JQ(2)=JTQ(2)                                                    
        JO(2)=JTO(2)                                                    
        JM(3)=JTM(2)                                                    
        JD(3)=JTD(2)                                                    
        JQ(3)=JTQ(2)                                                    
        JO(3)=JTO(2)                                                    
        JM(4)=JTM(3)                                                    
        JD(4)=JTD(3)                                                    
        JQ(4)=JTQ(3)                                                    
        JO(4)=JTO(3)                                                    
        JD(5)=JTD(3)                                                    
        JQ(5)=JTQ(3)                                                    
        JO(5)=JTO(3)                                                    
        JM(5)=JTM(4)                                                    
        JD(6)=JTD(4)                                                    
        JQ(6)=JTQ(4)                                                    
        JO(6)=JTO(4)                                                    
        JM(6)=JTM(4)                                                    
        JD(7)=JTD(4)                                                    
        JQ(7)=JTQ(4)                                                    
        JO(7)=JTO(4)                                                    
        JM(7)=JTM(4)                                                    
        JD(8)=JTD(4)                                                    
        JQ(8)=JTQ(4)                                                    
        JO(8)=JTO(4)                                                    
        JM(8)=JTM(5)                                                    
        JD(9)=JTD(5)                                                    
        JQ(9)=JTQ(5)                                                    
        JO(9)=JTQ(5)                                                    
        JM(9)=JTM(5)                                                    
        JD(10)=JTD(5)                                                   
        JQ(10)=JTQ(5)                                                   
        JO(10)=JTO(5)                                                   
        JD(11)=JTD(5)                                                   
        JQ(11)=JTQ(5)                                                   
        JO(11)=JTO(5)                                                   
        JD(12)=JTD(5)                                                   
        JQ(12)=JTQ(5)                                                   
        JO(12)=JTO(5)                                                   
        JM(10)=JTM(6)                                                   
        JD(13)=JTD(6)                                                   
        JQ(13)=JTQ(6)                                                   
        JO(13)=JTO(6)                                                   
        JD(14)=JTD(6)                                                   
        JQ(14)=JTQ(6)                                                   
        JO(14)=JTO(6)                                                   
        JQ(15)=JTQ(6)                                                   
        JO(15)=JTO(6)                                                   
  900   CONTINUE                                                        
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE CASCAD                                               
C  ***  MAIN CASCADE ROUTINE -- DOES ALL BOOKKEEPING...                 
        DIMENSION POPT(3),P(4),PC(3,210),PNL(210),POLPOS(210),          
     1  POLNEG(210),WIDTH(210),CONVC(210),SPORB(210),RADNT(20),         
     2  ZT(130),ZK(3,130),ZR(130),ENERGY(19),U(4),ZA0(130),ZA1(130),    
     3  ZA2(130),ZA(130),PC0(210),PC1(210),PC2(210),POP1(6),POP2(6)     
        COMMON/LOC001/IJK,ENERG,ECONS,ECONST,D2P1SM,D2P1S               
        COMMON/LOC008/IREAD,IW,IPUNCH,IPR                               
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC029/IRR,RR(18),RAU,RAD,RA(4),RD(4),RSA(4)             
        COMMON/LOC030/POP(6),JTM(6),JTD(6),JTQ(6),JTO(6)                
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC037/PL(20),NPOL(20),IPOL,CL1,CL2,IDE,PLN(210),IP8     
        COMMON/LOC038/A,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)               
        COMMON/LOC039/NOPT,NMAX,ALEXP                                   
        COMMON/LOC040/AMASSA,AMASSN,HBAR                                
        COMMON/LOC041/MPU,ICPU(200),IPN                                 
        DATA U/3HM+Q,3HD+O,3H Q ,3H O /                                 

C       For the first three electron shells,
C       we have 6 subshells, i.e
C       1s, 2s, 2p, 3s, 3p, 3d
        DO 50 I=1,6                                                     
C
C       We start off with POP being filled with 1s
        POP1(I) = POP(I)                                                
   50   CONTINUE                                                        
C
C       POP1 and POP seem to be identical, and PL is our normalised population
C       for the start of the cascade. This is doing the 3 sets of shells, K, L,
C       and M? 
C       Two electrons in 1s
        POPT(1) = 2.000*POP(1)                                          
        POP2(1) = POP1(1)                                               
C 
C       Two electrons in 2s and 6 in 2p
        POPT(2) = (2.000*POP(2)+6.000*POP(3))                           
        IERROR = 0                                                      
        MAXERR = 99                                                     
C
C       Two electrons in 3s, 6 in 3p, and 10 in 3d
        POPT(3) = (2.000*POP(4)+6.000*POP(5)+10.000*POP(6))             
        MS = 3                                                          
        POP2(2) = 1.000                                                 
C
C       POPT contains the number of electrons in the first 3 energy levels
C       So at this point we are simply normalising to 1 electron
        IF(POPT(2).GT.1.0E-20) POP2(2)=POP2(2)/POPT(2)                  
C
        POP2(3) = 1.000                                                 
C
        IF(POPT(2).GT.1.0E-20) POP2(3)=POP2(3)/POPT(2)                  
        POP2(4) = 1.000                                                 
C
        IF(POPT(3).GT.1.0E-20) POP2(4)=POP2(4)/POPT(3)                  
        POP2(5) = 1.000                                                 
C
        IF(POPT(3).GT.1.0E-20) POP2(5)=POP2(5)/POPT(3)                  
        POP2(6) = 1.000                                                 
C
        IF(POPT(3).GT.1.0E-20) POP2(6)=POP2(6)/POPT(3)                  

C       POP2 contains the normalisation for the first 6 subshells
C
        RLYMAN = 0.000                                                  
C     
C       This is the K shell refilling weight changed to atomic units
        RK = WIDTHK/HBAR                                                
        NU = NMAX*(NMAX+1)/2                                            
C
C       NU is the total number of states
C       Loop over all states and set PNL to 0
        DO 100 J=1,NU                                                   
        PNL(J) = 0.000                                                  
        IF(IP8.NE.0)  PNL(J) = PLN(J)                                   
        POLPOS(J) = 0.000                                               
        POLNEG(J) = 0.000                                               
        WIDTH(J) = 0.000                                                
        PC0(J) = 0.000                                                  
        PC1(J) = 0.000                                                  
        PC2(J) = 0.000                                                  
        CONVC(J) = 0.000                                                
        SPORB(J) = 0.000                                                
C
C       
C       MS = 3
        DO 100 IS=1,MS                                                
        PC(IS,J) = 0.000                                                
  100   CONTINUE                                                        
C
C       This tells us the start of the top of the cascade
C       i.e number of IUPAC states - NMAX
        MU = NMAX*(NMAX-1)/2                                            
C
C       NMAX is the total number of energy levels i.e 20
C       This entire loop is for the top energy level
        DO 200 J=1,NMAX                                                 
C
C       This tells us where we need to start in the array
C       to fill in populations
        JJ = MU+J        
C
C       PL is our statistical distribution of choice
C       and here we fill in the top energy level of PNL
        IF(IP8.EQ.0) PNL(JJ) = PL(J)                                    
C
C       Now we start filling in the topmost energy level
C       TODO figure out what PC is 
        ! PC0 represents the probability of there being 0 electrons at the start
        ! PC1  ""                                       1 electron
        ! PC2  ""                                       2 electrons
        ! So at the start, we are saying that we have a full K-shell
        ! Makes sense as nothing has happened yet!
        PC1(JJ) = 2.000-2.000*POP(1)                                    
        PC2(JJ) = 2.000*POP(1) - 1.000                                  
C
C      
        IF(POP(1).GT.0.500)  GO TO 150                                  
C    
C       If 1s is less than half occupied, then we know
C       there are definitely not 2 electrons
        PC0(JJ)=1.000-2.000*POP(1)                                      
        PC1(JJ)=2.000*POP(1)                                            
        PC2(JJ)=0.00                                                    
C
C       If we have some population, do these calculations
  150   IF(PNL(JJ).GT.1.0E-20)  POLPOS(JJ)=(1.0 + 2.0/FLOAT(2*J-1))/3.0 
        IF(PNL(JJ).GT.1.0E-20.AND.J.NE.1)  POLNEG(JJ)=(1.0 - 2.0/FLOAT( 
     1  2*J-1))/3.0                                                     
C  
C       Now we loop over K, L, M subshells?
        DO 200 IS=1,MS                                                  

C       POPT contains the population of each subshell
C       i.e 2, 8, 18
        PC(IS,JJ) = POPT(IS)                                            
C
C
  200   CONTINUE                                                        
C 
C       At this point we have found the population in the top energy level
C
C       Printing
        IPR1 = MOD(IPR/4,2)                                             
        IPR2 = MOD(IPR/8,2)                                             
        IPR3 = MOD(IPR/16,2)                                            
        IPR4 = MOD(IPR/32,2)                                            
C
C       BEGIN MAIN LOOP
C       This loops over each energy level 
C       and we then start from the top down
        DO 4000 I1=1,NMAX                                               
        N1 = NMAX+1-I1                                                  

C       For the very last level, skip everything
        IF(N1.EQ.1)  GO TO 4000                                         
        IF(IPR1.EQ.0.AND.N1.GT.3)  WRITE(IW,300)                        
  300   FORMAT(1H1)                                                     
     0  IF(IPR1.EQ.1)  WRITE(IW,400)                                    
  400   FORMAT(/1X,120(1H*)/)                                           
     0  IF(IPR1.EQ.0.AND.N1.LE.3)  WRITE(IW,450)                        
  450   FORMAT(/////)                                                   
C    0  CALL CHECK(N1)                                                  
C
        DO 500 I=1,20                                                   
        RADNT(I) = 0.000                                                
  500   CONTINUE                                                        
C      
C       For our value of N1, we loop over it again and find the
C       corresponding l-value
        DO 3000 I2=1,N1                                                 
C
C       L1 is our subshell L value
        L1 = I2-1                                                       
C
C       This sums up all of the rates coming out of a certain place
        RATEGT = 0.000                                                  
C
C       2x rate to account for 2 electrons in the K shell (?)
        RATE0 = 2.000*RK                                                
        RATE1 = RK                                                      
C       Sum of radiative rates
        RATERD = 0.000                                                  

C       K1 is the number in the IUPAC notation i.e K3, K1 = 3
C       i.e position in the PNL array
        K1 = N1*(N1-1)/2 + I2                                           
C
C       This triggers if we are at the top energy level
C       so we skip all this
        IF(I1.LE.1)  GO TO 700                                          
C
C       If we have no population, and we aren't in the 
C       top energy level, then we skip the next things
        IF(PNL(K1).LT.1.0E-20.AND.I1.GT.1)                              
     1  GO TO 620                                                       
C
C       So, if we have some population, then we find the normalisation
C       XNORM is the reciprocal of the population in a given state
        XNORM = 1.000/PNL(K1)                                           
C
C       So we are dividing the number of electrons that fill up the shell
C       by the actual population
        PC0(K1) = PC0(K1)*XNORM                                         
        PC1(K1) = PC1(K1)*XNORM                                         
        PC2(K1) = PC2(K1)*XNORM                                         

C       Loop from 2 to 3 
        DO 600 IS=2,MS                                                  
        PC(IS,K1) = PC(IS,K1)*XNORM                                     
        IF(PC(IS,K1).LE.0.000)  PC(IS,K1) = 0.000                       
  600   CONTINUE                                                        
        GO TO 640                                                       
C
C       We have jumped here if we have no population in our state
C       Say that the electrons are full
  620   PC(2,K1) = POPT(2)                                              
        PC(3,K1) = POPT(3)                                              
C
C 
C  ***  IF STATE IS NOT POPULADED, THE ELECTRONS ARE AS IN THE BEGINNING
C       We are encountering a state with no muonic population,
C       So we say that the K-shell is full with 100% probability
        PC0(K1) = 0.000                                                 
        PC1(K1) = 2.0-2.0*POP2(1)                                       
        PC2(K1) = 2.0*                                                  
     1  POP2(1) - 1.0                                                   
C
        IF(POP2(1).GE.0.5)  GO TO 630                                   
        PC0(K1) = 1.0 -                                                 
     2  2.0*POP2(1)                                                     
        PC1(K1) = 2.0*POP2(1)                                           
        PC2(K1) = 0.000                                                 
  630   PC(1,K1) = PC1(K1) + 2.000*PC2(K1)                              
        GO TO 700                                                       
  640   PC(1,K1) = PC1(K1) + 2.000*PC2(K1)                              
        IF(PC(1,K1).LT.0.000)  PC(1,K1) = 0.000                         
        IF(PC(1,K1).GT.2.000)  PC(1,K1) = 2.000                         
        DO 650 IS=1,MS                                                  
C
C       If IPC is 1, then we have no depletion
C       In this case, the population is set back to the original
C       value of 2, 8, 18 
        IF(IPC(IS).NE.0)  PC(IS,K1)=POPT(IS)                            
  650   CONTINUE                                                        
C
C       IF IPC is 0, then we have depletion with refilling
C       for the K shell
        IF(IPC(1).EQ.0)  GO TO 675                                      
        PC2(K1) = 2.000*POP1(1) - 1.000                                 
        PC1(K1)=2.-2.*POP1(1)                                           
        PC0(K1)=0.                                                      
C
C
        IF(POP1(1).GE..5) GO TO 675                                     
        PC2(K1)=0.0                                                     
        PC1(K1)=2.0*POP1(1)                                             
        PC0(K1)=1.0-2.0*POP1(1)                                         
C
C
  675   POLPOS(K1) = POLPOS(K1)*XNORM*FLOAT(2*L1+1)/FLOAT(L1+1)         
        IF(NPOL(N1)-L1.LE.1.AND.NPOL(N1)-L1.GE.0)  POLPOS(K1) =         
     1  (1.000 + 2.000/FLOAT(2*L1+1))/3.000                             
        IF(L1.LE.0)  GO TO 700                                          
        POLNEG(K1) = POLNEG(K1)*XNORM*FLOAT(2*L1+1)/FLOAT(L1)           
        IF(NPOL(N1)-L1.LE.1.AND.NPOL(N1)-L1.GE.0)  POLNEG(K1) =         
     1  (1.000 - 2.000/FLOAT(2*L1+1))/3.000                             
C   
C       This is the point at which all of the branches converge
  700   CONTINUE                                                        
C
C       Skip if we are in the bottom energy level
        IF(N1.EQ.1)  GO TO 3000                                         
        K = 0                                                           
        IF(IPR2.EQ.0)  WRITE(IW,800)                                    
C
C       Add up these things. Related to depletion somehow
C       It is the probability of 0, 1, 2 electrons in K shell added up
C       This should add to 1
        PC012 = PC0(K1)+PC1(K1)+PC2(K1)                                 
        IF(ABS(PC012-1.000).GT.1.000E-3*FLOAT(I1))  WRITE(IW,850)N1,L1, 
     1  PC0(K1),PC1(K1),PC2(K1)                                         
        IF(ABS(PC012).LT.1.0E-20)  PC012=1.000                          

C       Divide each of these by the sum
        PC0(K1) = PC0(K1)/PC012                                         
        PC1(K1) = PC1(K1)/PC012                                         
        PC2(K1) = PC2(K1)/PC012                                         
C
  800   FORMAT(//)                                                      
  850   FORMAT(53H0*** WARNING *** PROBABILITIES OF THE POPULATION OF T,
     1  60HHE K-SHELL ARE SIGNIFICANTLY OFF FROM BEING NORMALIZED PROPE,
     2  8HRLY. ***/20H *** HAPPENED AT N1=,I2,4H L1=,I2,5H PC0=,1PE12.5,
     3  5H PC1=,E12.5,5H PC2=,E12.5,24H NORMALIZATION FORCED.../)       
     0  POP(1) = 0.500*PC(1,K1)                                         
        POP(2) = POP2(2)*PC(2,K1)                                       
        POP(3) = POP2(3)*PC(2,K1)                                       
        POP(4) = POP2(4)*PC(3,K1)                                       
        POP(5) = POP2(5)*PC(3,K1)                                       
        POP(6) = POP2(6)*PC(3,K1)                                       
C
        DO 1100 I3=1,7                                                  
C 
C       Here we only wanting to find l values that are, at most,
C       3 units away as this is the furthest transition that can occur
C       for Auger octupole transitions
        L2 = L1-4+I3                                                    
C
C       Negative l clearly isn't allowed
        IF(L2.LT.0)  GO TO 1100                                         
C
C       Loop over each IUPAC state?
        DO 1000 I4=1,N1                                                 
        N2 = N1 - I4 + 1                                                
C 
C       Catch all invalid states
        IF(N2.LE.L2)  GO TO 1000                                        
        IF(N2.EQ.N1.AND.(N2.NE.2.OR.L1.NE.0.OR.L2.NE.1))  GO TO 1000    
        IF(N1.EQ.N2.AND.ESPM.LE.1.000E-20)  GO TO 1000                  
        K = K + 1                                                       

C
C       This is checking for the 2p -> 2s transition
C       which is ignored by default
        IF(N1.NE.N2)  GO TO 900                                         
        IJK = 1                                                         
        ENERG = ESPM                                                    
        GO TO 950                                                       
C
C
  900   IJK = 0                                                         
  950   POPQ = POP(1)                                                   
        POP(1) = 1.000                                                  
C
C       Calculation of rates
C       We give pairs of states to get the rate back
C       N1 must always be greater than N2, and L1 must be 
C       greater than L2
C       This function returns us 2 common variables as well:
C       RAD and RSA
C       RSA is the individual Auger rates for each shell
C       RAD is the total radiative rate from all contributions
C       ZT is the total rate coming into the state
C       write(4,*) N1, L1, N2, L2, "PAIRS", K
        ZT(K) = RATE(N1,L1,N2,L2) + 1.000E-10                           
C      write(4,*) "R=", ZT(K), "FOR", n1, l1, n2, l2
C
C       We start summing up all of the rates coming out of the 
C       initial state
        RATEGT = RATEGT + ZT(K)                                         
        IF(POP(1).LT.1.0E-20)  POP(1)=1.000                             
C      
C       Scale K auger by the K electron population
        RSA(1)=RSA(1)/POP(1)                                            
        POP(1) = POPQ                                                   
C
C       For the top energy level,  RAD is 0
        ZR(K) = RAD                                                     
C  
C       Sum of all of the radiative transitions
        RATERD = RATERD + RAD                                           
C
C       This is K shell rates?
        ZA(K) = RSA(1)                                                  
C
C       Total rate minus K shell Auger rate
        XR = ZT(K) - RSA(1)                                             
C
C       I believe this is equation 4 from Vogel PRA 1973
C       XR is radiative width, RSA is Auger width, and RK is total width
C       of a 1s electron vacancy
        T0 = XR + 0.500*RSA(1) + RK                                     
C
C       Again this is equation 4
        G0 = XR + 2.000*RK + 1.000E-10                                  
C
C       Here we sum the T0 stuff
        RATE1 = RATE1+T0-RK +1.000E-10                                  
C
C       Here we are summing up the radiative + L and M auger contributions
        RATE0 = RATE0+XR +1.000E-10                                     
C
C       Here, the PC sum is the number of electrons in K
C       So this is the difference of total rate and K auger rate scaled by elec
C       population
        XM = ZT(K)-ZA(K)*(PC0(K1)+0.500*PC1(K1))                        
        IERR = 0                                                        
C
C
        IF(ABS(XM).LT.1.0E-10.OR.ABS(T0).LT.1.0E-10.OR.                 
     1  ABS(G0).LT.1.0E-10.OR.ABS(ZT(K)).LT.1.0E-10)  IERR = 1          
C
C
        IF(IERR.EQ.0)  GO TO 995                                        
        IERROR = IERROR + 1                                             
        WRITE(IW,960)                                                   
        WRITE(IW,990)N1,L1,N2,L2,XM,T0,G0,K,ZT(K),PC0(K1),PC1(K1),      
     1  PC2(K1),ZR(K),XR,ZA(K),ZA0(K),ZA1(K),ZA2(K),YM,XL,ZK(1,K),      
     2  ZK(2,K),ZK(3,K)                                                 
        WRITE(IW,970)                                                   
C
C
        IF(ABS(XM).LT.1.0E-10)  XM=1.0                                  
        IF(ABS(T0).LT.1.0E-10)  T0=1.0                                  
        IF(ABS(G0).LT.1.0E-10)  G0=1.0                                  
        IF(ABS(ZT(K)).LT.1.0E-10)  ZT(K)=1.0                            
        IF(IERROR.GT.MAXERR)  WRITE(IW,980)                             
        IF(IERROR.GT.MAXERR)  STOP                                      
C
C
  960   FORMAT(//51H *** INTERNAL PROGRAM ERROR IN ROUTINE CASCAD *** B,
     1  60HAD CHOICE OF PARAMETERS HAS RESULTED IN A DIVISION BY ZERO *,
     2  2H**/55H *** ONE OF THE VARIABLES /XM,T0,G0,ZT(K),RATE0,RATE1,R,
     3  58HATEGT/ IS ZERO *** NEXT LINES GIVE MORE INFORMATION... ***//)
  970   FORMAT(//51H *** STANDARD FIXUP TAKEN (ZERO VARIABLE PUT EQUAL ,
     1  35HTO ONE), EXECUTION CONTINUING ***//)                         
  980   FORMAT(//51H *** TOO MANY DIVIDE CHECKS *** EXECUTION TERMINATE,
     1  17HD --- NO DUMP ***//)                                         
  990   FORMAT(8H *** N1=,I2,4H L1=,I2,5H, N2=,I2,4H L2=,I2/5X,3HXM=,1PE
     1  12.5,4H T0=,E12.5,4H G0=,E12.5,5H ZT(,I3,2H)=,E12.5/5X,6HPC0(K1,
     2  2H)=,E12.5,9H PC1(K1)=,E12.5,9H PC2(K1)=,E12.5/5X,6HZR(K)=,E12.5
     3  ,4H XR=,E12.5,7H ZA(K)=,E12.5/5X,4HZA0=,E12.5,5H ZA1=,E12.5,    
     4  5H ZA2=,E12.5,4H YM=,E12.5,4H XL=,E12.5/5X,8HZK(1,K)=,E12.5,    
     5  9H ZK(2,K)=,E12.5,9H ZK(3,K)=,E12.5/)                           
C
C       All of the following is from Vogel PRA
C       particularly, equations 8a, 8b, 8c
C
C       
  995   ZK(3,K) = PC(3,K1) - RSA(3)/XM                                  
        IF(ZK(3,K).LT.0.0)  ZK(3,K)=0.0                                 
C
C       This is equation 8c
        ZA2(K) = XR/ZT(K)*(PC2(K1) + RK*PC1(K1)/T0 + 2.000*RK*RK*PC0(K1)
     1  /T0/G0)                                                         
C
C       This term is the probability of the K shell having 1 electron
C       when the muon has reached state the final state
C       This is equation 8b
        ZA1(K) = ZA(K)*PC2(K1)/ZT(K) + (XR + RK*ZA(K)/ZT(K))/T0*(PC1(K1)
     1  + 2.000*RK*PC0(K1)/G0)                                          
C
C       This is equation 8a
C       PC1: Electron population probability for 1 K electron at time 0
C       ZA:  Auger width
C       T0:  sum of radiative, auger, and refilling widths
C       PC0: Electron population probability for 0 K electrons at time 0
C       G0: radiative rate + 2*refilling width
C       RK: refilling width
        ZA0(K) = PC1(K1)*ZA(K)/(2.000*T0) + PC0(K1)/G0*(XR + RK*ZA(K)/  
     1  T0)                                                             
C
C       What is this?
        YM = ZA(K)*(PC2(K1)/ZT(K) + PC1(K1)*(0.500 + RK/ZT(K))/T0 + PC0 
     1  (K1)*RK/T0/G0*(1.000 + 2.000*RK/ZT(K)))                         
C
C
        XL = RK*(PC1(K1)/T0 + 2.000*PC0(K1)/T0/G0*(T0 + RK))            
C
C
        ZK(1,K) = PC(1,K1) - YM + XL                                    
        ZK(2,K) = PC(2,K1) - RSA(2)/XM - XL                             
        IF(ZK(1,K).LT.0.00)  ZK(1,K)=0.                                 
        IF(ZK(2,K).LT.0.)  ZK(2,K)=0.0                                  
C
C       n2 and l2 loop end
 1000   CONTINUE                                                        
 1100   CONTINUE                                                        
        WIDTH(K1) = (RATEGT*PC2(K1) + (RATE1-RK)*PC1(K1) +              
     1  (RATE0-2.0*RK)*PC0(K1))*HBAR                                    
        IF(WIDTH(K1).LT.0.) WIDTH(K1)=0.0                               
        IERR = 0                                                        
        IF(ABS(RATE0).LT.1.0E-10.OR.ABS(RATE1).LT.1.0E-10.OR.           
     1  ABS(RATEGT).LT.1.0E-10)  IERR = 1                               
        IF(IERR.EQ.0)  GO TO 1175                                       
        IERROR = IERROR + 1                                             
        WRITE(IW,960)                                                   
        WRITE(IW,1150)N1,L1,RATE0,RATE1,RATEGT                          
        WRITE(IW,970)                                                   
 1150   FORMAT(8H *** N1=,I2,4H L1=,I2,7H RATE0=,1PE12.5,7H RATE1=,E12.5
     1  ,8H RATEGT=,E12.5,4H ***/)                                      
        IF(ABS(RATE0).LT.1.0E-10)  RATE0 = 1.000                        
        IF(ABS(RATE1).LT.1.0E-10)  RATE1 = 1.000                        
        IF(ABS(RATEGT).LT.1.0E-10)  RATEGT = 1.000                      
        IF(IERROR.GT.MAXERR)  WRITE(IW,980)                             
        IF(IERROR.GT.MAXERR)  STOP                                      
 1175   CONVC(K1) = RATERD/(WIDTH(K1)/HBAR - RATERD + 1.000E-10)        
        IF(CONVC(K1).LT.0.000)  CONVC(K1) = HUGE(9.999E+37)
        IF(L1.NE.0)  SPORB(K1) = 0.150*Z**4/FLOAT(N1**3*L1*(L1+1))      
!
!       Reset the counter
        K = 0                                                           
C
C       Loop over second states again
        DO 2100 I3=1,7                                                  
        L2 = L1-4+I3                                                    
        IF(L2.LT.0)  GO TO 2100                                         
        DO 2000 I4=1,N1                                                 
        N2 = N1-I4+1                                                    
C
C       Same selection rules as before
        IF(N2.LE.L2)  GO TO 2000                                        
        IF(N2.EQ.N1.AND.(N2.NE.2.OR.L1.NE.0.OR.L2.NE.1))  GO TO 2000    
        IF(N2.EQ.N1.AND.ESPM.LE.1.000E-20)  GO TO 2000                  
C
C       Calculate the IUPAC number for our next state down
        K = K + 1                                                       
        K2 = N2*(N2-1)/2 + L2 + 1                                       
C
C       THIS IS WHERE THE POPULATIONS FOR EACH LEVEL ARE CALCULATED
        BNORM = PNL(K1)*(PC2(K1)*ZT(K)/RATEGT + PC1(K1)*(ZT(K) - 0.500* 
     1  ZA(K) + RK*ZT(K)/RATEGT)/RATE1 + PC0(K1)*(ZT(K) - ZA(K) + 2.000*
     2  RK/RATE1*(ZT(K) - 0.500*ZA(K) + RK*ZT(K)/RATEGT))/RATE0)        
C
C       With no depletion, BNORM is reduced to a simpler form
C       BNORM = PNL(K1)*(PC2(K1)*ZT(K)/RATEGT)
        DO 1200 IS=1,MS                                                 
        PC(IS,K2) = PC(IS,K2) + ZK(IS,K)*BNORM                          
 1200   CONTINUE                                                        
C
C       This is where our population comes into play
C       We start filling in the population for the next levels down
C       The population is summed from all possible previous states
        PNL(K2) = PNL(K2) + BNORM                                       
C
C 
C       PC are population probabilities of K, L, M electrons at muon start
C       ZA are population probabilities of K, L, M electrons at muon end
        write(4,*) "ZA2", PC2(K2), ZA2(K), BNORM
        PC2(K2) = PC2(K2) + ZA2(K)*BNORM                                
        PC1(K2) = PC1(K2) + ZA1(K)*BNORM                                
        PC0(K2) = PC0(K2) + ZA0(K)*BNORM                              
C
C       THIS IS THE RADIATIVE INTENSITIES
C       This is very similar to BNORM except it uses only the radiative
C       rates 
C       ZR is the radiative rate
C       RADINT is the intensity for the given transition
C       RATE1 is the refilling rate
        RADINT = PNL(K1)*ZR(K)*(PC2(K1)/RATEGT + PC1(K1)*(1.000 + RK/   
     1  RATEGT)/RATE1 + PC0(K1)*(1.000 + 2.000*RK/RATE1*(1.000 + RK/    
     2  RATEGT))/RATE0)                                                 
C
C       Here we sum up each contribution coming to this state.
C       which accounts for all of the different angular momentum channels
C       in order to get a total intensity
        RADNT(N2) = RADNT(N2) + RADINT                                  
C
C
C  ***  FOR PRINTOUT OF RATES ONLY, SCHROEDINGER ENERGIES ARE USED ***  
        ENERGY(N2) = ECONST*(1.000/FLOAT(N2*N2)-1.000/FLOAT(N1*N1))     
C
C       2p -> 2s energy
        IF(N1.EQ.N2)  ENERGY(N2) = ESPM                                 
C
C       2p -> 1s energy
        IF(N2.EQ.1.AND.D2P1SM.GT.1.0E-20)  ENERGY(N2)=ENERGY(N2)+D2P1SM 
     1  -0.750*ECONST                                                   
        ENERGY(N2) = ENERGY(N2)*AMASSE                                  
C
C       deltaL
        LL = IABS(L1-L2)                                                
C
C
        CALL POPJ(L1,L2,LL,P)                                           
        IF(IPOL.NE.0)  GO TO 1500                                       
        LI = LL + 1                                                     
        J1U = 2*L1 + 1                                                  
        J2U = 2*L2 + 1                                                  
        J1D = J1U - 2                                                   
        J2D = J2U - 2                                                   
        IF(LI.GT.1)  GO TO 1300                                         
        POLPOS(K2) = POLPOS(K2) + BNORM*P(1)*POLPOS(K1)                 
        POLNEG(K2) = POLNEG(K2) + BNORM*P(2)*POLNEG(K1)                 
        GO TO 1500                                                      
 1300   IF(L2.GT.L1)  GO TO 1400                                        
        POLPOS(K2) = POLPOS(K2) + BNORM*(P(1)*POLPOS(K1)*BETA(L1,J1U,L2,
     1  J2U,LL) + P(2)*POLNEG(K1)*BETA(L1,J1D,L2,J2U,LL))               
        IF(J1D.EQ.-1.OR.J2D.EQ.-1)  GO TO 1500                          
        POLNEG(K2) = POLNEG(K2) + BNORM*P(3)*POLNEG(K1)                 
     1  *BETA(L1,J1D,L2,J2D,LL)                                         
        GO TO 1500                                                      
 1400   POLPOS(K2) = POLPOS(K2) + BNORM*P(1)*POLPOS(K1)                 
     1  *BETA(L1,J1U,L2,J2U,LL)                                         
        POLNEG(K2) = POLNEG(K2) + BNORM*(P(2)*POLPOS(K1)*BETA(L1,J1U,L2,
     1  J2D,LL) + P(3)*POLNEG(K1)*BETA(L1,J1D,L2,J2D,LL))               
 1500   LI = LL                                                         
        IF(LL.EQ.0)  LI = 2                                             
        LK = LL + 1                                                     
C       
C       RADINT PASSED THROUGH
C       This routine builds up the X-ray table
        CALL CODE(N1,L1,N2,L2,LI,RADINT)                                
C
        IF(IPR2.EQ.0)  WRITE(IW,1600)N1,L1,N2,L2,U(LK),ENERGY(N2),RADINT
 1600   FORMAT(4H N1=,I2,5H, L1=,I2,5H, N2=,I2,5H, L2=,I2,6H, MUL=,A3,  
     1  4H, E=,-6PF11.8,5H(MEV),6H, RAD=,1PE12.4,10H(PER MUON))         
C
C       Radint is radiative intensities
C       This catches all transitions to K1 and sums them up
C       We want this to be unity
     0  IF(L1.EQ.1.AND.L2.EQ.0.AND.N2.EQ.1)  RLYMAN = RLYMAN + RADINT   
C
C
        IF(MPU.LE.0.OR.IPN.NE.0)  GO TO 2000                            
        DO 1800 IT=1,MPU                                                
        IF(ICPU(IT)/4194304.NE.1)  GO TO 1800                           
        N1J = MOD(ICPU(IT),32)                                          
        N2J = MOD(ICPU(IT)/2048,32)                                     
        L1J = MOD(ICPU(IT)/32,32)                                       
        L2J = MOD(ICPU(IT)/65536,32)                                    
        IF(N1J.NE.N1.OR.N2J.NE.N2.OR.L1J.NE.L1.OR.L2J.NE.L2)  GO TO 1800
        WRITE(IPUNCH,1700)N1,L1,N2,L2,ENERGY(N2),RADINT,IDE             
 1700   FORMAT(2H 1,2I3,4X,2I3,7X,-6PF9.6,1PE12.4,26X,I5)               
C 
 1800   CONTINUE                                                        
 2000   CONTINUE                                                        
 2100   CONTINUE                                                        
C    
C       end of l1 loop
 3000   CONTINUE                              
        IF(IPR3.EQ.0)  WRITE(IW,3100)                                   
 3100   FORMAT(//)                                                      
     0  NO = N1-1                                                       
        IF(IPR3.EQ.0) WRITE(IW,3200)(N1,N2,ENERGY(N2),RADNT(N2),N2=1,NO)
 3200   FORMAT(4H N1=,I2,5H, N2=,I2,4H, E=,-6PF11.8,5H(MEV),6H, RAD=,   
     1  1PE16.9,12H(NORMALIZED)/)                                       
     0  IF(MPU.EQ.0.OR.IPN.NE.0)  GO TO 4000                            
        DO 3400 IT=1,MPU                                                
        IF(ICPU(IT)/4194304.NE.2)  GO TO 3400                           
        N1J = MOD(ICPU(IT),32)                                          
        N2 = MOD(ICPU(IT)/2048,32)                                      
        IF(N1J.EQ.N1)  WRITE(IPUNCH,3300)N1,N2,ENERGY(N2),RADNT(N2),IDE 
 3300   FORMAT(2H 2,I3,7X,I3,10X,-6PF9.6,1PE12.4,26X,I5)                
 3400   CONTINUE                                                       

C       We have arrived here if we are in the bottom energy level and
C       everything has been computed
C       Now we just want to write everything out
C       END OF N1 LOOP
 4000   CONTINUE                                                        
C
C       We can write out the sum of the Lyman intensities
        WRITE(IW,4025)RLYMAN                                            
 4025   FORMAT(//1H ,60(1H*)/40H LYMAN SERIES (NP-1S) SUM OF INTENSITIES
     1  ,3H = ,F7.5,10H(PER MUON)/                                      
     2  60H DEVIATION FROM UNITY IS THE SUM OF TRANSITION INTENSITIES  /
     3  60H ENDING IN THE 1S STATE NOT THROUGH AN NP-1S RADIATIVE TRAN-/
     4  60H SITION.  EXPERIMENTALLY SET TO UNITY FOR NORMALIZATION     /
     5  1H ,60(1H*)//)                                                  
C
     0  IF(IPR4.NE.0)  RETURN                                           
        WRITE(IW,4050)                                                  
 4050   FORMAT(1H1)                                                     
     0  WRITE(IW,4100)                                                  
        PC(1,1) = PC1(1) + 2.000*PC2(1)                                 
 4100   FORMAT(//51H N1 L1  POPULATION  POLAR.UP   POLAR.DN   WID (EV) ,
     1  60H  RAD/AUG    S-O(EV)    K-ELECT    L-ELECT    M-ELECT    ***,
     2  6H****  /1H ,120(1H-))                                          
C
C       Now we want to write out the population table
C       Here we loop over energy level
     0  DO 5000 M1=1,NMAX                                               
C
C       Calculate the actual N value
        N1 = NMAX+1-M1                                                  
C
C       Now loop over each channel
        DO 4900 LL1=1,N1                                                
C
C       Convert to actual angular momentum
        L1 = LL1-1                                                      
C
C       Calculate the array index for the population
        K1 = N1*(N1-1)/2 + LL1                                          
        POLNEG(K1) = -POLNEG(K1)*(FLOAT(L1)+0.5)/(FLOAT(L1)-0.5)        
        IF(IPOL.EQ.1)  GO TO 4300                                       
        IF(N1.EQ.1)  WRITE(IW,4110)N1,L1,PNL(K1),POLPOS(1),(PC(IS,1),   
     1  IS=1,3)                                                         
        IF(N1.EQ.1)  GO TO 4900                                         
        IF(L1.EQ.0.AND.PNL(K1).GT.1.0E-20)  WRITE(IW,4120)N1,L1,PNL(K1),
     1  POLPOS(K1),WIDTH(K1),CONVC(K1),(PC(IS,K1),IS=1,3)               
        IF(L1.EQ.0.AND.PNL(K1).GT.1.0E-20)   GO TO 4900                 
        IF(PNL(K1).LE.1.0E-20)  WRITE(IW,4130)N1,L1                     
        IF(PNL(K1).LE.                                                  
     1  1.0E-20)  GO TO 4900                                            
C
C
 4110   FORMAT(1X,I2,1H,,I2,1X,1P2E11.3,4(5X,3H***,3X),3E11.3)          
 4120   FORMAT(1X,I2,1H,,I2,1X,1P2E11.3,5X,3H***,3X,2E11.3,5X,3H***,3X, 
     1  3E11.3)                                                         
 4130   FORMAT(1X,I2,1H,,I2,1X,9(5X,3H***,3X),2X,13HNOT POPULATED)      
        WRITE(IW,4200)N1,L1,PNL(K1),POLPOS(K1),POLNEG(K1),WIDTH(K1),    
     1  CONVC(K1),SPORB(K1),(PC(IS,K1),IS=1,MS)                         
        GO TO 4900                                                      
 4200   FORMAT(1X,I2,1H,,I2,1X,1P10E11.3)                               
 4300   IF(N1.EQ.1)  WRITE(IW,4400)N1,L1,PNL(K1),(PC(IS,1),IS=1,3)      
        IF(N1.EQ.1)  GO TO 4900                                         
        IF(L1.EQ.0.AND.PNL(K1).GT.1.0E-20)                              
     1  WRITE(IW,4500)N1,L1,PNL(K1),WIDTH(K1),CONVC(K1),(PC(IS,K1),IS=  
     2  1,3)                                                            
        IF(L1.EQ.0.AND.PNL(K1).GT.1.0E-20)  GO TO 4900                  
        IF(PNL(K1).LE.1.0E-20)  WRITE(IW,4130)N1,L1                     
        IF(PNL(K1).LE.1.0E-20)  GO TO 4900                              
        WRITE(IW,4600)N1,L1,PNL(K1),WIDTH(K1),CONVC(K1),SPORB(K1),      
     1  (PC(IS,K1),IS=1,MS)                                             
 4400   FORMAT(1X,I2,1H,,I2,1X,1PE11.3,5(5X,3H***,3X),3E11.3)           
 4500   FORMAT(1X,I2,1H,,I2,1X,1PE11.3,2(5X,3H***,3X),2E11.3,5X,3H***,3X
     1  ,3E11.3)                                                        
 4600   FORMAT(1X,I2,1H,,I2,1X,1PE11.3,2(5X,3H***,3X),7E11.3)           
 4900   CONTINUE                                                        
        WRITE(IW,4950)                                                  
C
C       We have written out the full subshell population for a given N
 4950   FORMAT(1H ,120(1H-))                                            
C
C       Move down to the next energy level
 5000   CONTINUE                                                        
C
C       This is full of 1s again
        DO 5100 I=1,6                                                   
        POP(I) = POP1(I)                                                
 5100   CONTINUE                                                        
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE CODE(N1,L1,N2,L2,L,CC)                               
        COMMON/LOC008/IR,IW,IP,IPR                                      
C  ***  PUTS INFORMATION FOR LINE INTENSITIES IN COMPACT FORM TO SAVE   
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC033/M,E(1000),AI(1000),IA(1000),ENERGY(20,40)         
        IF((L+L1-L2)/2*2.NE.L+L1-L2)  WRITE(IW,100)N1,L1,N2,L2,L        
  100   FORMAT(/50H *** ERROR *** ANGULAR MOMENTA DO NOT MATCH AT N1=,I2
     1  ,5H, L1=,I2,5H, N2=,I2,5H, L2=,I2,16H, MULTIPOLARITY=,I1,4H ***)
     0  IF(CC.LE.CLIMIT)  RETURN                                        
        IF(M.GT.0995)  WRITE(IW,200)                                    
  200   FORMAT(53H *** ATTENTION *** OVERFLOWING CAPABILITY OF SORTING ,
     1  50H*** PLEASE RESTRICT CRITERIA FOR LINES TO PASS ***)          
     0  E12 = 0.000E00                                                  
        E22 = 0.000E00                                                  
        E11 = ENERGY(N1,2*L1+1)                                         
        IF(L1.NE.0)  E12 = ENERGY(N1,2*L1)                              
        E21 = ENERGY(N2,2*L2+1)                                         
        IF(L2.NE.0)  E22 = ENERGY(N2,2*L2)                              
        IF(E21-E11.LE.ELOW.AND.N1.NE.N2)  RETURN                        
        IF(E21-E11.GE.EHIGH.AND.N1.NE.N2)  RETURN                       
        IF(N1.NE.N2)  GO TO 250                                         
        E11 = 0.000                                                     
        EZ = 0.500*(ENERGY(2,2)-ENERGY(2,3))                            
        E21 = ESP*1.000E-6 - EZ                                         
        E22 = E21 + 2.000*EZ                                            
  250   IA0 = 4194304*L + 65536*L2 + 2048*N2 + 32*L1 + N1               
        IF(L.LE.0.OR.L.GT.3)  WRITE(IW,300)N1,L1,N2,L2,L                
  300   FORMAT(48H *** ERROR *** UNEXPECTED QUANTUM NUMBERS AT N1=,I2,  
     1  5H, L1=,I2,5H, N2=,I2,5H, L2=,I2,4H, L=,I1,4H ***)              
     0  GO TO(1000,2000,3000),L                                         
 1000   LL = MAX0(L1,L2)                                                
        AN = CC/FLOAT(4*LL*LL-1)                                        
        E(M) = E21-E11                                                  
        IF(L1.EQ.LL)  E(M+1)=E21-E12                                    
        IF(L2.EQ.LL)  E(M+1)=E22-E11                                    
        E(M+2)=E22-E12                                                  
        IA(M) = IA0 + 2098176                                           
        IF(L1.EQ.LL)  IA(M+1) = IA0 + 2097152                           
        IF(L2.EQ.LL)  IA(M+1) = IA0 + 1024                              
C      
C       LL is the maximum angular momentum
C       AI(M) = AN * (l+1)*(2l-1)
        IA(M+2) = IA0                                                   
        AI(M) = AN*FLOAT((LL+1)*(2*LL-1))                               
        AI(M+1) = AN                                                    
        AI(M+2) = AN*FLOAT((LL-1)*(2*LL+1))                             
        M = M + 3                                                       
        IF(LL.EQ.1)  M = M-1                                            
        RETURN                                                          
 2000   IF(L1.EQ.L2)  GO TO 2500                                        
        LL = MAX0(L1,L2)                                                
        AN = CC/FLOAT((2*LL-3)*(2*LL+1))                                
        E(M) = E21-E11                                                  
        IF(L1.EQ.LL)  E(M+1)=E21-E12                                    
        IF(L2.EQ.LL)  E(M+1)=E22-E11                                    
        E(M+2) = E22-E12                                                
        IA(M) = IA0 + 2098176                                           
        IF(L1.EQ.LL)  IA(M+1) = IA0 + 2097152                           
        IF(L2.EQ.LL)  IA(M+1) = IA0 + 1024                              
        IA(M+2) = IA0                                                   
        AI(M) = AN*FLOAT((LL+1)*(2*LL-3))                               
        AI(M+1) = 2.000*AN                                              
        AI(M+2) = AN*FLOAT((LL-2)*(2*LL+1))                             
        M = M + 3                                                       
        IF(LL.LE.2)  M = M-1                                            
        RETURN                                                          
 2500   IF(L1.EQ.0)  RETURN                                             
        AN = CC/FLOAT((2*L1+1)**2)                                      
        E(M) = E21-E11                                                  
        E(M+1) = E22-E11                                                
        E(M+2) = E21-E12                                                
        E(M+3) = E22-E12                                                
        IA(M) = IA0 + 20978176                                          
        IA(M+1) = IA0 + 1024                                            
        IA(M+2) = IA0 + 2097152                                         
        IA(M+3) = IA0                                                   
        AI(M) = AN*FLOAT((L1+2)*(2*L1-1))                               
        AI(M+1) = 3.000*AN                                              
        AI(M+2) = 3.000*AN                                              
        AI(M+3) = AN*FLOAT((L1-1)*(2*L1+3))                             
        M = M + 4                                                       
        IF(L1.EQ.1)  M = M-1                                            
        RETURN                                                          
 3000   IF(IABS(L1-L2).EQ.1)  GO TO 3500                                
        LL = MAX0(L1,L2)                                                
        AN = CC/FLOAT((2*LL-5)*(2*LL+1))                                
        E(M) = E21-E11                                                  
        IF(L1.EQ.LL)  E(M+1)=E21-E12                                    
        IF(L2.EQ.LL)  E(M+1)=E22-E11                                    
        E(M+2) = E22-E12                                                
        IA(M) = IA0 + 2098176                                           
        IF(L1.EQ.LL)  IA(M+1) = IA0 + 2097152                           
        IF(L2.EQ.LL)  IA(M+1) = IA0 + 1024                              
        IA(M+2) = IA0                                                   
        AI(M) = AN*FLOAT((2*LL-5)*(LL+1))                               
        AI(M+1) = 3.000*AN                                              
        AI(M+2) = AN*FLOAT((2*LL+1)*(LL-3))                             
        M = M + 3                                                       
        IF(LL.LE.3)  M = M-1                                            
        RETURN                                                          
 3500   LL = MAX0(L1,L2)                                                
        AN = CC/FLOAT(4*LL*LL-1)                                        
        E(M) = E21-E11                                                  
        E(M+1) = E22-E11                                                
        E(M+2) = E21-E12                                                
        E(M+3) = E22-E12                                                
        IA(M) = IA0 + 20978176                                          
        IA(M+1) = IA0 + 1024                                            
        IA(M+2) = IA0 + 2097152                                         
        IA(M+3) = IA0                                                   
        AI(M) = AN*FLOAT((2*LL-3)*(LL+2))                               
        AI(M+1) = 5.000*AN                                              
        AI(M+2) = 6.000*AN                                              
        AI(M+3) = AN*FLOAT((2*LL+3)*(LL-2))                             
        M = M + 4                                                       
        IF(LL.LE.2)  M = M-1                                            
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE SORT                                                 
        DIMENSION ST(5),C(5),LL(3)                                      
C  ***  ARRANGES LINE INTENSITIES IN ENERGY FOR THE X-RAY TABLE         
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC033/M,E(1000),AI(1000),IA(1000),ENERGY(20,40)         
        COMMON/LOC035/ICC,CD(5),EA,EB,IDIR                              
        COMMON/LOC037/PL(20),NPOL(20),IPOL,CL1,CL2,IDE,PLN(210),IP8     
        COMMON/LOC038/AB,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)              
        COMMON/LOC041/MPU,ICPU(200),IPX                                 
        DATA STAR,BLANK/1H*,1H /,LL/3HDIP,3HQUA,3HOCT/                  
        ICH = 0                                                         
        IF((EA-99.)**2 + (EB-99.)**2.LT.1.000E-20)  ICH = 1             
        CC = CLIMIT**(1.000/6.000)                                      
        DO 100 I=1,5                                                    
        C(I) = CC**I                                                    
        IF(ICC.NE.0)  C(I) = CD(I)                                      
  100   CONTINUE                                                        
        M = M-1                                                         
        IF(M.LE.0)  RETURN                                              
        M1 = M                                                          
        DO 170 I=1,M                                                    
  120   IF(AI(I).GT.CLIMIT)  GO TO 170                                  
        IF(I.GT.M1)  GO TO 170                                          
        M1 = M1 - 1                                                     
        DO 140 J=I,M1                                                   
        E(J) = E(J+1)                                                   
C       AI IS ALL INTENSITIES
        AI(J) = AI(J+1)                                                 
        IA(J) = IA(J+1)                                                 
  140   CONTINUE                                                        
        GO TO 120                                                       
  170   CONTINUE                                                        
        M = M1                                                          
        WRITE(IW,200)Z,AB,M                                             
  200   FORMAT(1H1/1H ,120(1H*)/17H0ATOMIC NUMBER = ,F5.1,4X,6HATOMIC   
     1  ,10H WEIGHT = ,F7.3,4X,18HNUMBER OF LINES = ,I5//1X,120(1H*)/)  
     0  IF(M.EQ.1)  GO TO 500                                           
        MM = M-1                                                        
        DO 400 I=1,MM                                                   
        I1 = I+1                                                        
        DO 300 J=I1,M                                                   
        IF(E(J).GE.E(I))  GO TO 300                                     
        E0 = E(J)                                                       
        AI0 = AI(J)                                                     
        IA0 = IA(J)                                                     
        E(J) = E(I)                                                     
        AI(J) = AI(I)                                                   
        IA(J) = IA(I)                                                   
        E(I) = E0                                                       
        AI(I) = AI0                                                     
        IA(I) = IA0                                                     
  300   CONTINUE                                                        
  400   CONTINUE                                                        
  500   N = 0                                                           
        EE = 0.000                                                      
        AA = 0.000                                                      
        LC = 7                                                          
        DO 1500 I=1,M                                                   
  600   II=I+N                                                          
        IF(IA(II)/16777216.EQ.1)  GO TO 1500                            
        IA0 = IA(II)                                                    
        N1 = MOD(IA0,32)                                                
        L1 = MOD(IA0/32,32)                                             
        J1 = (2*L1-1) + 2*MOD(IA0/1024,2)                               
        N2 = MOD(IA0/2048,32)                                           
        L2 = MOD(IA0/65536,32)                                          
        J2 = (2*L2-1) + 2*MOD(IA0/2097152,2)                            
        L = MOD(IA0/4194304,4)                                          
        L = LL(L)                                                       
        EN = 1000.000*E(II)                                             
        A = AI(II)                                                      
        CH = (EN-EA)/AMIN1(EB,1.000E+20)                                
        DO 700 K=1,5                                                    
        ST(K)=BLANK                                                     
        IF(A.GT.C(K))  ST(K)=STAR                                       
  700   CONTINUE                                                        
        IF(ICH.EQ.0)  WRITE(IW,800)N1,L1,J1,N2,L2,J2,L,EN,A,            
     1  ST,CH                                                           
C       A HOLDS THE INTENSITIES
        IF(ICH.EQ.1)  WRITE(IW,805)N1,L1,J1,N2,L2,J2,L,EN,A,ST          
  800   FORMAT(4H N1=,I2,4H,L1=,I2,4H,J1=,I2,9H/2    N2=,I2,4H,L2=,I2,  
     1  4H,J2=,I2,8H/2    L=,A3,7H    EN=,F12.6,13H(KEV)    INT=,1PE11.4
     2  ,4X,5A1,4X,3HCH=,0PF9.3)                                        
  805   FORMAT(4H N1=,I2,4H,L1=,I2,4H,J1=,I2,9H/2    N2=,I2,4H,L2=,I2,  
     1  4H,J2=,I2,8H/2    L=,A3,7H    EN=,F12.6,13H(KEV)    INT=,1PE11.4
     2  ,4X,5A1)                                                        
     0  LC = LC + 1                                                     
        IF(MPU.LE.0)  GO TO 830                                         
        IF(IPX.NE.0)  GO TO 830                                         
        DO 820 IT=1,MPU                                                 
        IF(ICPU(IT)/4194304.NE.0)  GO TO 820                            
        N1J = MOD(ICPU(IT),32)                                          
        N2J = MOD(ICPU(IT)/2048,32)                                     
        L1J = MOD(ICPU(IT)/32,32)                                       
        L2J = MOD(ICPU(IT)/65536,32)                                    
        J1J = MOD(ICPU(IT)/1024,2)                                      
        J2J = MOD(ICPU(IT)/2097152,2)                                   
        J1K = 2*L1J-1 + 2*J1J                                           
        J2K = 2*L2J-1 + 2*J2J                                           
        IF(N1J.NE.N1.OR.N2J.NE.N2.OR.L1J.NE.L1.OR.L2J.NE.L2.OR.J1K.NE.  
     1  J1.OR.J2K.NE.J2)  GO TO 820                                     
        WRITE(IPUNCH,810)N1,L1,J1J,N2,L2,J2J,EN,A,IDE                   
  810   FORMAT(2H 0,2I3,I2,2X,2I3,I2,5X,-3PF9.6,1PE12.4,26X,I5)         
  8200  CONTINUE                                                        
  830   IF(LC.GE.60)  WRITE(IW,850)                                     
  850   FORMAT(1H1)                                                     
     0  IF(LC.GE.60)  LC=0                                              
        IA(II) = IA(II)+16777216                                        
        EE = EE+E(II)*A                                                 
        AA = AA+A                                                       
  900   N = N+1                                                         
        IPN = I+N                                                       
        IF(E(IPN)-E(I).GE.ERES)  GO TO 1000                             
        IF(IA(IPN)/16777216.EQ.1)  GO TO 900                            
        N11 = MOD(IA(IPN),32)                                           
        N22 = MOD(IA(IPN)/2048,32)                                      
        IF(N1.EQ.N11.AND.N2.EQ.N22)  GO TO 600                          
 1000   IF(ABS(A-AA).LT.1.0E-20)  GO TO 1300                            
        EN = 1000.000*EE/AA                                             
        CH = (EN-EA)/AMIN1(EB,1.000E+20)                                
        DO 1100 K=1,5                                                   
        ST(K) = BLANK                                                   
        IF(AA.GT.C(K))  ST(K)=STAR                                      
 1100   CONTINUE                                                        
        IF(ICH.EQ.0)  WRITE(IW,1200)N1,N2,EN,AA,ST,CH                   
        IF(ICH.EQ.1)  WRITE(IW,1250)N1,N2,EN,AA,ST                      
 1200   FORMAT(4H N1=,I2,1X,13(1H-),4X,3HN2=,I2,1X,13(1H-),10X,6HAV.EN=,
     1  F12.6,13H(KEV) TOT.IN=,1PE11.4,4X,5A1,7H AV.CH=,0PF9.3)         
 1250   FORMAT(4H N1=,I2,1X,13(1H-),4X,3HN2=,I2,1X,13(1H-),10X,6HAV.EN=,
     1  F12.6,13H(KEV) TOT.IN=,1PE11.4,4X,5A1)                          
     0  LC = LC + 1                                                     
        IF(LC.GE.60)  WRITE(IW,850)                                     
        IF(LC.GE.60)  LC=0                                              
 1300   WRITE(IW,1400)                                                  
 1400   FORMAT(1H )                                                     
     0  LC = LC + 1                                                     
        IF(LC.GE.60)  WRITE(IW,850)                                     
        IF(LC.GE.60)  LC=0                                              
        N = 0                                                           
        AA = 0.000                                                      
        EE = 0.000                                                      
 1500   CONTINUE                                                        
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RATE(N1,L1,N2,L2)                                      
C  ***  MASTER RATE ROUTINE -- INTERPRETS OPTIONS AND CALLS OTHER RATES 
C  ***  POINT DIRAC OR INPUTED ENERGIES USED EXCLUSIVELY IN RATES ***   
        DIMENSION Y(3),IY(3),IP1(7),IP2(7),IP3(7),IDM(3),IDD(5,5),      
     1  IDQ(4,4),IDO(4,4),IDR(18)                                       
        COMMON/LOC001/IJK,ENERGY,ECONS,ECONST,D2P1SM,D2P1S              
        COMMON/LOC002/BEM(3),ZSA(3),BE(3)                               
        COMMON/LOC003/K0,K1,K2,K3                                       
        COMMON/LOC004/NN0(3),NN1(7),NN2(7),NN3(7)                       
        COMMON/LOC005/R0(3),R1(7),R2(7),R3(7)                           
        COMMON/LOC006/JP1(7),JP2(7),JP3(7),IQ1(7),IQ2(7),IQ3(7)         
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC028/M1(7),M2(7),M3(7),YC(4),IDB                       
        COMMON/LOC029/IRR,RR(18),RAU,RAD,RA(4),RD(4),RSA(4)             
        COMMON/LOC033/M,E(1000),AI(1000),IA(1000),ENE(20,40)            
        DATA IDM/4HM-1T,4HM-2T,4HM-3T/                                  
        DATA IDD/4HD-RT,4HD-1T,4HD-2T,4HD-3T,4H****,4H****,4HD-1S,4HD-2S
     1  ,4HD-3S,4H****,2*4H****,4HD-2P,4HD-3P,4H****,3*4H****,4HD-3D,   
     2  4H****,4*4H****,4H****/                                         
        DATA IDQ/4HQ-RT,4HQ-1T,4HQ-2T,4HQ-3T,4H****,4HQ-1S,4HQ-2S,4HQ-3S
     1  ,2*4H****,4HQ-2P,4HQ-3P,3*4H****,4HQ-3D/                        
        DATA IDO/4H8-RT,4H8-1T,4H8-2T,4H8-3T,4H****,4H8-1S,4H8-2S,4H8-3S
     1  ,2*4H****,4H8-2P,4H8-3P,3*4H****,4H8-3D/                        
        IRR = 0                                                     
C
C       This seems to be the midpoint of the two energy levels rounded
C       down
        N12 = (N1+N2+1)/2                                               
        RATE = 0.000                                                    
        RAU = 0.000                                                     
        RAD = 0.000                                                     
C
C
        DO 50 I=1,4                                                     
        RA(I) = 0.000                                                   
        RD(I) = 0.000                                                   
        RSA(I) = 0.000                                                  
   50   CONTINUE                                                        
C
C       delta L + 1, add 1 for the computed GOTO I believe
        L = IABS(L1-L2)+1                                               
C
C       Transitions with delta L > 3 are not allowed
        IF(L.GT.4)  RETURN                                              
        ENEM = ENERGY                                                   
        IF(IJK.NE.0)  GO TO 100                                         
        LL1 = 1                                                         
        LL2 = 1                                                         
C
C       Indexing fix for l = 0 states
        IF(L1.EQ.0)  LL1=0                                              
        IF(L2.EQ.0)  LL2=0                                              
C
C       This block calculates the indices for a given
C       pair of fine doublets for the 
C       two energy levels
        LJ11=2*L1+1                                                     
        LJ12=LJ11-LL1                                                   
        LJ21=2*L2+1                                                     
        LJ22=LJ21-LL2                                                   
C
C       Calculate average energy of fine doublet
C       for starting and end states
C       write(4,*) Lj11, lj12, l1, l2
C       write(4,*) ENE(N2, LJ21), ENE(N2, LJ22)
        ENE1 = 0.500*(ENE(N1,LJ11) + ENE(N1,LJ12))                      
        ENE2 = 0.500*(ENE(N2,LJ21) + ENE(N2,LJ22))                      
C
C       If we have got a 2p -> 1s transition then we simply use
C       the experimental energy given in the input file
        IF(N2.EQ.1.AND.D2P1SM.GT.1.000E-20)  ENE2=D2P1SM*AMASSE*1.0E-6 +
     1  0.500*(ENE(2,2) + ENE(2,3))                                     
C
C       Convert the transition energy to dimensionless?
        ENEM = (ENE2-ENE1)*1.000E6/AMASSE                               
C       write(4,*) ENEM, "ENERGY", N1, L1, N2, L2
C   
C       Loop over K, L, M 
  100   DO 200 I=1,3                                                    
        IY(I) = 1                                                       
C
C       BEM contains the electronic binding energies of K, L and M shell
C       converted to correct units
C       Therefore T is the energy difference between the transition energy
C       and the corresponding electron binding energy
        T = ENEM-BEM(I)                                                 
C
C       This GOTO tells us that transition energy < binding energy
C       which means we don't have a valid Auger transition
        IF(T.LE.0.000)  GO TO 200                                       
C
        IY(I) = 0                                                       
C 
C       ZSA is effective charge for each shell multiplied by fine structure
C       constant
C       Y is defined in equation A.21 of Akylas PhD thesis
C       and is a parameter for the ejected continuum electron
C       We have one value for each of the K, L, and M shells,
C       making Y a 3 length array
        Y(I) = ZSA(I)/SQRT(T*T+2.000*T)                                 
        write(4,*) Y(I), BEM(I)*AMASSE, ENEM*AMASSE,I,"ALCULATED HERE"
C       if(L==2)write(4,*) Y(I), "Y VALUE", N1, L1, N2, L2
C
C
  200   CONTINUE                                                        
C
C       delta L = 1,2,3,4
C       i.e decides which multipole we have
        GO TO (1000,2000,3000,4000),L                                   
C
C       Monopole case
C       If we have no monopoles, then skip down to quadrupole
 1000   IF(K0.EQ.0)  GO TO 3000                                         
C
C
C       K0 is the number of cases we have of monopoles
C       i.e this means that we have three possible choices:
C       K, L, or M Auger electron emission
C       There are no radiative monopoles
        DO 1100 I=1,K0                                                  
        R0(I) = 0.000                                                   
C 
C       NN = 1,2, 3
C       This is K, L, M shells.
C       This lets us get the correct Auger energies and Y variable
        NN = NN0(I)                                                     
C
C       This GOTO tells us that if we haven't got enough energy to eject an
C       electron, then we skip the monopole transition
        IF(IY(NN).NE.0)  GO TO 1100                                     
C
C 
C       If y is less than YC, then we skip the penetration for the transition
        IF(Y(NN).LT.YC(NN))  GO TO 1100                                 
C
C       This calculates the monopole rate
C       which only occurs for Auger transitions
C       Essentially, we give 2 states, an electron shell,
C       and an electron continuum parameter
        R0(I) = RMON(N1,L1,N2,L2,NN,Y(NN))                              
C
C       Add up the Auger monopole rate to a running total
C       as well as a running Auger total
        RATE = RATE+R0(I)                                               
        RAU = RAU+R0(I)                                                 
C
C       Count up the number of rates actually calculated
        IRR = IRR+1                                                     
C       Store the rate
        RR(IRR) = R0(I)                                                 
C
C       NN = 0 represents radiative rates
        IF(NN.EQ.0)  GO TO 1100                                         
C
C       Add up the rate for specific shells
        RSA(NN) = RSA(NN) + R0(I)                                       
        IDR(IRR) = IDM(NN)                                              
 1100   CONTINUE                                                        
C
C       Store the monopole Auger rate 
        RA(1) = RAU                                                     
        GO TO 3000                                                      
C   
C       Dipole case. If we have no cases of it, then skip
 2000   IF(K1.EQ.0)  GO TO 4000                                         
C
C       Loop over 4 dipole cases: radiative, and Auger for K, L, M
        DO 2100 I=1,K1                                                  
C       write(4,*) I, "DIPOLE START LOOP", IY(I)
        R1(I) = 0.000                                                   
C
C       This runs from 0 to 3
        NN = NN1(I)                                                     
        MM = NN                                                         
C
C       This is the case of radiative rates
C       I believe we do this just to make a function call
C       interface correct
        IF(MM.EQ.0)  MM=1                                               
        IP1(I) = 0                                                      
C
C       Here, we just want a radiative rate!
        IF(NN.EQ.0)  GO TO 2050                                         
C
C       Not enough energy to emit an Auger electron
C       so skip the dipole transition
        IF(IY(NN).NE.0)  GO TO 2100                                     
        IF(Y(NN).GT.YC(NN))  IP1(I) = JP1(I)                            
C
C       Here we calculate the dipole transition rate
C       M1 is a list of subshells we are considering
C       IP1 = 0 means we don't consider penetration
C       ENEM is needed for radiative rates
C       We have specifically come here to calculate a radiative rate
 2050   IF(IP1(I).EQ.0)  R1(I)=RDIPU(N1,L1,N2,L2,NN,ENEM,Y(MM),M1(I))   
C       if(NN/= 0)write(4,*) "DIPAUG", R1(I), N1, L1, N2, L2
C       if(NN== 0)write(4,*) "DIPRAD", R1(I), N1, L1, N2, L2
C        write(4,*) "RATE COMPUTED FOR ", N1, L1, N2, L2, R1(I)
C
C       Skip the following things as we have just calculated our rate above
        IF(NN.EQ.0)  GO TO 2075                                         
C
C       We want penetration,
C       IQ1 seems to always be zero
C       So only one of these if statements is done
        IF(IP1(I).NE.0.AND.N12.GE.IQ1(I))  R1(I) = RDIP(N1,L1,N2,L2,NN, 
     1  M1(I),Y(MM))                                                    

        IF(IP1(I).NE.0.AND.N12.LT.IQ1(I))  R1(I) = RDIPU(N1,L1,N2,L2,NN,
     1  ENEM,Y(MM),M1(I))                                               
C       write(4,*) "AUGER RATE COMPUTED"
C
C       Add the dipole rate up
 2075   RATE = RATE+R1(I)                                               
C
C       Add up the radiative bit
        IF(NN.EQ.0)  RAD=RAD+R1(I)                                      
C     
C       Add up the Auger bit
        IF(NN.NE.0)  RAU=RAU+R1(I)                                      
        IRR = IRR+1                                                     
        RR(IRR) = R1(I)                                                 
        MN = M1(I) + 1                                                  
        IF(NN.EQ.0)  MN = 1                                             
        IDR(IRR) = IDD(NN+1,MN)                                         
        IF(NN.EQ.0)  GO TO 2100                                         
C 
C       RSA seems to sum up the Auger rates as a function of K,L,M
        RSA(NN) = RSA(NN) + R1(I)                                       
 2100   CONTINUE                                                        
C       write(4,*) "RSA", RSA, N1, L1, N2, L2
        RD(2) = RAD                                                     
        RA(2) = RAU                                                     
        GO TO 4000                                                      
C
C       If we don't have any quadrupole cases, then
C       skip
 3000   IF(K2.EQ.0.OR.L1+L2.EQ.0)  GO TO 5000                           
C
C       Loop over the quadrupole cases
        DO 3100 I=1,K2                                                  
        R2(I) = 0.000                                                   
        NN = NN2(I)                                                     
        MM = NN                                                         
        IF(MM.EQ.0)  MM=1                                               
        IP2(I) = 0                                                      
C
C       Radiative only
        IF(NN.EQ.0)  GO TO 3050                                         
C     
C       No Auger allowed
        IF(IY(NN).NE.0)  GO TO 3100                                     
        IF(Y(NN).GT.YC(NN))  IP2(I) = JP2(I)                            
C
C       Here we calculate the quadrupole transition rates
 3050   IF(IP2(I).EQ.0)  R2(I)=RQUAU(N1,L1,N2,L2,NN,ENEM,Y(MM),M2(I))   
C       if(NN== 0)write(4,*) "QUADRAD", R2(I), N1, L1, N2, L2
C       if(NN/= 0)write(4,*) "QUADAUG", R2(I), N1, L1, N2, L2

        IF(NN.EQ.0)  GO TO 3075                                         
        IF(IP2(I).NE.0.AND.N12.GE.IQ2(I))  R2(I)=RQUA(N1,L1,N2,L2,NN,   
     1  M2(I),Y(MM))                                                    
        IF(IP2(I).NE.0.AND.N12.LT.IQ2(I))  R2(I)=RQUAU(N1,L1,N2,L2,NN,  
     1  ENEM,Y(MM),M2(I))                                               
 3075   RATE = RATE+R2(I)                                               
C
C       If NN = 0, then we have done something radiative
C       and so we add it up 
        IF(NN.EQ.0)  RD(3)=RD(3)+R2(I)                                  
        IF(NN.EQ.0)  RAD=RAD+R2(I)                                      
C
C       We then add up Auger rates
        IF(NN.NE.0)  RA(3)=RA(3)+R2(I)                                  
        IF(NN.NE.0)  RAU=RAU+R2(I)                                      
C
C       Count the rate 
        IRR = IRR+1                                                     
        RR(IRR) = R2(I)                                                 
        MN = M2(I) + 1                                                  
        IF(NN.EQ.0)  MN = 1                                             
        IDR(IRR) = IDQ(NN+1,MN)                                         
        IF(NN.EQ.0)  GO TO 3100                                         
        RSA(NN) = RSA(NN) + R2(I)                                       
 3100   CONTINUE                                                        
        GO TO 5000                                                      
C
C       If we have no octupole cases allowed,
C       or if L1 =0, L2 = 1 or vice versa
 4000   IF(K3.EQ.0.OR.L1+L2.EQ.1)  GO TO 5000                           
C
C       We have our 4 cases by default
        DO 4100 I=1,K3                                                  
        R3(I) = 0.000                                                   
C
C       Figure out which rate we want
        NN = NN3(I)                                                     
        MM = NN                                                         
        IF(MM.EQ.0)  MM=1                                               
        IP3(I) = 0                                                      
C
C       Go and perform radiative rate
        IF(NN.EQ.0)  GO TO 4050                                         
        IF(IY(NN).NE.0)  GO TO 4100                                     
        IF(Y(MM).GT.YC(NN))  IP3(I) = JP3(I)                            
C
C       Here we calculate octupole transition rates
 4050   IF(IP3(I).EQ.0)  R3(I)=ROCTU(N1,L1,N2,L2,NN,ENEM,Y(MM),M3(I))   
C       if(NN==0)write(4,*) "OCTRAD", R3(I), N1, L1, N2, L2
C       if(NN/=0)write(4,*) "OCTAUG", R3(I), N1, L1, N2, L2
        IF(NN.EQ.0)  GO TO 4075                                         
        IF(IP3(I).NE.0.AND.N12.GE.IQ3(I))  R3(I)=ROCT(N1,L1,N2,L2,NN,   
     1  M3(I),Y(MM))                                                    
        IF(IP3(I).NE.0.AND.N12.LT.IQ3(I))  R3(I)=ROCTU(N1,L1,N2,L2,NN,  
     1  ENEM,Y(MM),M3(I))                                               
C
C       Radiative octupole
 4075   IF(NN.EQ.0)  RD(4)=RD(4)+R3(I)                                  
        IF(NN.EQ.0)  RAD=RAD+R3(I)                                      
C
C       Auger octupole
        IF(NN.NE.0)  RA(4)=RA(4)+R3(I)                                  
        IF(NN.NE.0)  RAU=RAU+R3(I)                                      
        RATE = RATE+R3(I)                                               
        IRR = IRR+1                                                     
        RR(IRR) = R3(I)                                                 
        MN = M3(I) + 1                                                  
        IF(NN.EQ.0)  MN = 1                                             
        IDR(IRR) = IDO(NN+1,MN)                                         
        IF(NN.EQ.0)  GO TO 4100                                         
        RSA(NN) = RSA(NN) + R3(I)                                       
 4100   CONTINUE                                                        
 5000   IF(IDB.EQ.0)  GO TO 5400                                        
        IF(IRR.EQ.0)  GO TO 5200                                        
        WRITE(IW,5100)N1,L1,N2,L2,Y,RATE,(IDR(I),RR(I),I=1,IRR)         
 5100   FORMAT(1X,I2,1H,,I2,3H - ,I2,1H,,I2,3F7.3,1PE10.3 /7(1X,A4,     
     1  1PE12.4))                                                       
 5200   IF(IRR.EQ.0)  WRITE(IW,5300)N1,L1,N2,L2,Y                       
 5300   FORMAT(1X,I2,1H,,I2,3H - ,I2,1H,,I2,3F7.3,14H ***NO RATE***)    
 5400   IF(IRR.EQ.0)  RETURN                                            
        DO 5500 I=1,IRR                                                 
        IF(RR(I).LT.0.000)  WRITE(IW,5600)N1,L1,N2,L2,I,RR(I)           
 5500   CONTINUE                                                        
 5600   FORMAT(53H *** ERROR *** IN INTERNAL CALCULATION OF TRANSITION ,
     1  12HRATES AT N1=,I2,4H L1=,I2,5H, N2=,I2,4H L2=,I2,9H RATE NO=,I2
     2  ,7H RATE =,1PE13.5,4H ***)                                      
     0  RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RMON(N1,L1,N2,L2,N,Y)                                  
        INTEGER F                                                       
        REAL MATEL                                                      
C  ***  MONOPOLE RATE ROUTINE (PENETRATION ONLY)                        
        DIMENSION JM(10)                                                
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFFF,ALFA,AMASSE   
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC013/COEMON(30),EXPMON(30)                             
        COMMON/LOC017/F(6)                                              
        COMMON/LOC030/POP(6),DUM(24)                                    
        COMMON/LOC031/KM(10),JD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        EPIY = EXP(PI*Y)                                                
        COEFF = PICOEF*2.000*EPIY/(EPIY-1.000/EPIY)                     
        DO 500 JJ=1,10                                                  
        JM(JJ) = KM(JJ)                                                 
  500   CONTINUE                                                        
        IF(IJ(1).EQ.0)  GO TO 900                                       
        DO 600 JJ=1,10                                                  
        IF(Y.GT.YJ(1))  JM(JJ)=MIN0(KM(JJ),JJ1(1))                      
  600   CONTINUE                                                        
  900   GO TO (1000,2000,3000),N                                        
 1000   A1=0.000                                                        
        IF(F(1).EQ.1)  RETURN                                           
        DO 1100 JJ=1,3                                                  
        IF(JJ.GT.JM(1))  GO TO 1100                                     
        J=JK(JJ)                                                        
        B=Y**NE(JJ)                                                     
        A1=A1+COEMON(JJ)*MATEL(N1,L1,N2,L2,J+2,EXPMON(JJ),1)/B          
 1100   CONTINUE                                                        
        RMON=A1*A1*COEFF*POP(1)                                         
        RETURN                                                          
 2000   A1=0.000                                                        
        A2=0.000                                                        
        A3=0.000                                                        
        YY=(1.000+Y*Y)/(Y*Y)                                            
        IF(F(2)+F(3).EQ.2)  RETURN                                      
        DO 2200 JJ=1,3                                                  
        J=JK(JJ)                                                        
        B=Y**NE(JJ)                                                     
        IF(F(2).EQ.1)  GO TO 2100                                       
        IF(JJ.GT.JM(2))  GO TO 2100                                     
        A1=A1+COEMON(JJ+3)*MATEL(N1,L1,N2,L2,J+2,EXPMON(JJ+3),2)/B      
 2100   IF(F(3).EQ.1)  GO TO 2200                                       
        IF(JJ.GT.JM(3))  GO TO 2150                                     
        A2=A2+COEMON(JJ+ 6)*MATEL(N1,L1,N2,L2,J+3,EXPMON(JJ+ 6),2)/B    
 2150   IF(JJ.GT.JM(4))  GO TO 2200                                     
        A3=A3+COEMON(JJ+ 9)*MATEL(N1,L1,N2,L2,J+4,EXPMON(JJ+ 9),2)/B    
 2200   CONTINUE                                                        
        RMON=((A1+A2)**2*POP(2)+A3*A3*YY*POP(3))*COEFF                  
        RETURN                                                          
 3000   A1=0.000                                                        
        A2=0.000                                                        
        A3=0.000                                                        
        A4=0.000                                                        
        A5=0.000                                                        
        A6=0.000                                                        
        YY1=(1.000+Y*Y)/(Y*Y)                                           
        YY2=(1.000+Y*Y)*(4.000+Y*Y)/Y**4                                
        IF(F(4)+F(5)+F(6).EQ.3)  RETURN                                 
        DO 3300 JJ=1,3                                                  
        J=JK(JJ)                                                        
        B=Y**NE(JJ)                                                     
        IF(F(4).EQ.1)  GO TO 3100                                       
        IF(JJ.GT.JM(5))  GO TO 3100                                     
        A1=A1+COEMON(JJ+12)*MATEL(N1,L1,N2,L2,J+2,EXPMON(JJ+12),3)/B    
 3100   IF(F(5).EQ.1)  GO TO 3200                                       
        IF(JJ.GT.JM(6))  GO TO 3150                                     
        A2=A2+COEMON(JJ+15)*MATEL(N1,L1,N2,L2,J+3,EXPMON(JJ+15),3)/B    
 3150   IF(JJ.GT.JM(7))  GO TO 3200                                     
        A3=A3+COEMON(JJ+18)*MATEL(N1,L1,N2,L2,J+4,EXPMON(JJ+18),3)/B    
 3200   IF(F(6).EQ.1)  GO TO 3300                                       
        IF(JJ.GT.JM(8))  GO TO 3230                                     
        A4=A4+COEMON(JJ+21)*MATEL(N1,L1,N2,L2,J+4,EXPMON(JJ+21),3)/B    
 3230   IF(JJ.GT.JM(9))  GO TO 3260                                     
        A5=A5+COEMON(JJ+24)*MATEL(N1,L1,N2,L2,J+5,EXPMON(JJ+24),3)/B    
 3260   IF(JJ.GT.JM(10))  GO TO 3300                                    
        A6=A6+COEMON(JJ+27)*MATEL(N1,L1,N2,L2,J+6,EXPMON(JJ+27),3)/B    
 3300   CONTINUE                                                        
        RMON=((A1+A2+A3)**2*POP(4)+YY1*(A4+A5)**2*POP(5)+               
     1       YY2*A6*A6*POP(6))*COEFF                                    
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RDIP(N1,L1,N2,L2,N,M,Y)                                
        DIMENSION A(5),JD(14)                                           
C  ***  DIPOLE PENETRATION ROUTINE                                      
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC014/COEDIP(42),EXPDIP(42)                             
        COMMON/LOC018/F(9)                                              
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC023/COEDP(9)                                          
        COMMON/LOC030/POP(6),DUM(24)                                    
        COMMON/LOC031/JM(10),KD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        DOUBLE PRECISION ETPY, EXTPY
        INTEGER F                                                       
        REAL MATEL                                                      
        write(4,*) Y, "Y", "PI"
        ! PDJ THIS IS DODGY 
        if (Y > log(huge(1.0))/(2.0*PI)) ETPY = 1e100
        if (Y < log(huge(1.0))/(2.0*PI)) ETPY = EXP(2.000*PI*Y)                                  
        flush(4)
        FLUSH(IW)
        EXTPY = ETPY/(ETPY-1.000)                                       
        FLUSH(IW)
        YY = Y*Y                                                        
        P = EXP(Y*(2.000*ATAN(Y/FLOAT(N))-PI))                          
        DO 500 I=1,5                                                    
        A(I)=0.000                                                      
  500   CONTINUE                                                        
        DO 550 JJ=1,14                                                  
        JD(JJ) = KD(JJ)                                                 
  550   CONTINUE                                                        
        IF(IJ(2).EQ.0)  GO TO 900                                       
        DO 600 JJ=1,14                                                  
        IF(Y.GT.YJ(2))  JD(JJ) = MIN0(KD(JJ),JJ1(2))                    
  600   CONTINUE                                                        
  900   GO TO (1000,2000,3000),N                                        
 1000   A1 = 0.000                                                      
        IF(F(1).EQ.1)  GO TO 1200                                       
        DO 1100 JJ=1,3                                                  
        IF(JJ.GT.JD(1))  GO TO 1100                                     
        J = JK(JJ)                                                      
        B = Y**NE(JJ)       
C       write(4,*)MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ),1), "MATEL"
        A1 = A1+COEDIP(JJ)*MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ),1)/B        
 1100   CONTINUE                                                        
 1200   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)/YY                                                
        RDIP = COEDP(1)*EXTPY*PICOEF*ANGD*YG*(YF*P*AID*AMZZ(1)-A1)**2   
        RDIP = RDIP*POP(1)                                              
        RETURN                                                          
 2000   A1 =0.000                                                       
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        IF(M.EQ.2)  GO TO 2300                                          
        IF(F(2).EQ.1)  GO TO 2200                                       
        DO 2100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JD(2))  GO TO 2050                                     
        A1 = A1+COEDIP(JJ+3)*MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ+3),2)/B    
 2050   IF(JJ.GT.JD(3))  GO TO 2100                                     
        A1 = A1+COEDIP(JJ+ 6)*MATEL(N1,L1,N2,L2,J+4,EXPDIP(JJ+ 6),2)/B  
 2100   CONTINUE                                                        
 2200   YF = YY/(4.0+YY)                                                
        YG = (1.0+YY)/YY                                                
        A(1) = COEDP(2)*POP(2)*EXTPY*YG*(YF*P*AID*AMZZ(2)-A1)**2        
 2300   IF(M.EQ.1)  GO TO 2600                                          
        IF(F(3)+F(4).EQ.2)  GO TO 2500                                  
        DO 2400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JD(4))  GO TO 2350                                     
        A2 = A2+COEDIP(JJ+ 9)*MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ+ 9),2)/B  
 2350   IF(JJ.GT.JD(5))  GO TO 2400                                     
        A3 = A3+COEDIP(JJ+12)*MATEL(N1,L1,N2,L2,J+5,EXPDIP(JJ+12),2)/B  
 2400   CONTINUE                                                        
 2500   YF = YY/(4.0+YY)                                                
        A(2)=COEDP(3)*POP(3)*EXTPY*(YF*P*AID*AMZZ(2)-A2)**2             
        YF = YY*YY/(4.0+YY)**2                                          
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        A(3) = COEDP(4)*POP(3)*EXTPY*YG*(YF*P*AID*AMZZ(2)-A3)**2        
 2600   AT = A(1)+A(2)+A(3)                                             
        RDIP = PICOEF*ANGD*AT                                           
        RETURN                                                          
 3000   A1 = 0.000                                                      
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        A4 = 0.000                                                      
        A5 = 0.000                                                      
        IF(M.GT.1)  GO TO 3300                                          
        IF(F(5).EQ.1)  GO TO 3200                                       
        DO 3100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JD(6))  GO TO 3030                                     
        A1 = A1+COEDIP(JJ+15)*MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ+15),3)/B  
 3030   IF(JJ.GT.JD(7))  GO TO 3060                                     
        A1 = A1+COEDIP(JJ+18)*MATEL(N1,L1,N2,L2,J+4,EXPDIP(JJ+18),3)/B  
 3060   IF(JJ.GT.JD(8))  GO TO 3100                                     
        A1 = A1+COEDIP(JJ+21)*MATEL(N1,L1,N2,L2,J+5,EXPDIP(JJ+21),3)/B  
 3100   CONTINUE                                                        
 3200   YF = YY*(27.0+7.0*YY)/(9.0+YY)**2                               
        YG = (1.0+YY)/YY                                                
        A(1)=COEDP(5)*POP(4)*EXTPY*YG*(YF*P*AID*AMZZ(3)-A1)**2          
 3300   IF(M.NE.0.AND.M.NE.2)  GO TO 3600                               
        IF(F(6)+F(7).EQ.2)  GO TO 3500                                  
        DO 3400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JD(9))  GO TO 3325                                     
        A2 = A2+COEDIP(JJ+24)*MATEL(N1,L1,N2,L2,J+3,EXPDIP(JJ+24),3)/B  
 3325   IF(JJ.GT.JD(10))  GO TO 3350                                    
        A2 = A2+COEDIP(JJ+27)*MATEL(N1,L1,N2,L2,J+4,EXPDIP(JJ+27),3)/B  
 3350   IF(JJ.GT.JD(11))  GO TO 3375                                    
        A3 = A3+COEDIP(JJ+30)*MATEL(N1,L1,N2,L2,J+5,EXPDIP(JJ+30),3)/B  
 3375   IF(JJ.GT.JD(12))  GO TO 3400                                    
        A3 = A3+COEDIP(JJ+33)*MATEL(N1,L1,N2,L2,J+6,EXPDIP(JJ+33),3)/B  
 3400   CONTINUE                                                        
 3500   YF = YY*(3.0+YY)/(9.0+YY)**2                                    
        A(2) = COEDP(6)*POP(5)*EXTPY*(YF*P*AID*AMZZ(3)-A2)**2           
        YF = YY*YY/(9.0+YY)**2                                          
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        A(3) = COEDP(7)*POP(5)*EXTPY*YG*(YF*P*AID*AMZZ(3)-A3)**2        
 3600   IF(M.NE.0.AND.M.NE.3)  GO TO 3900                               
        IF(F(8)+F(9).EQ.2)  GO TO 3800                                  
        DO 3700 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JD(13))  GO TO 3650                                    
        A4 = A4+COEDIP(JJ+36)*MATEL(N1,L1,N2,L2,J+5,EXPDIP(JJ+36),3)/B  
 3650   IF(JJ.GT.JD(14))  GO TO 3700                                    
        A5 = A5+COEDIP(JJ+39)*MATEL(N1,L1,N2,L2,J+7,EXPDIP(JJ+39),3)/B  
 3700   CONTINUE                                                        
 3800   YF = YY*YY/(9.0+YY)**2                                          
        YG = (1.0+YY)/YY                                                
        A(4) = COEDP(8)*POP(6)*EXTPY*YG*(YF*P*AID*AMZZ(3)-A4)**2        
        YF = (YY/(9.0+YY))**3                                           
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        A(5) = COEDP(9)*POP(5)*EXTPY*YG*(YF*P*AID*AMZZ(3)-A5)**2        
 3900   AT = A(1)+A(2)+A(3)+A(4)+A(5)                                   
        RDIP = PICOEF*ANGD*AT                                           
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RQUA(N1,L1,N2,L2,N,M,Y)                                
        DIMENSION A(6),JQ(15)                                           
C  ***  QUADRUPOLE PENETRATION ROUTINE                                  
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC015/COEQUA(45),EXPQUA(45)                             
        COMMON/LOC019/F(10)                                             
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC024/COEQ(11)                                          
        COMMON/LOC030/POP(6),DUM(24)                                    
        COMMON/LOC031/JM(10),JD(14),KQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        INTEGER F                                                       
        REAL MATEL                                                      
        if (Y > log(huge(1.0))/(2.0*PI)) ETPY = 1e100
        if (Y < log(huge(1.0))/(2.0*PI)) ETPY = EXP(2.000*PI*Y)                                  
C       ETPY = EXP(2.000*PI*Y)                                          
        EXTPY = ETPY/(ETPY-1.000)                                       
        YY = Y*Y                                                        
        P = EXP(Y*(2.000*ATAN(Y/FLOAT(N))-PI))                          
        DO 500 I=1,6                                                    
        A(I)=0.000                                                      
  500   CONTINUE                                                        
        DO 550 JJ=1,15                                                  
        JQ(JJ) = KQ(JJ)                                                 
  550   CONTINUE                                                        
        IF(IJ(3).EQ.0)  GO TO 900                                       
        DO 600 JJ=1,15                                                  
        IF(Y.GT.YJ(3))  JQ(JJ) = MIN0(KQ(JJ),JJ1(3))                    
  600   CONTINUE                                                        
  900   GO TO (1000,2000,3000),N                                        
 1000   A1 = 0.000                                                      
        IF(F(1).EQ.1)  GO TO 1200                                       
        DO 1100 JJ=1,3                                                  
        IF(JJ.GT.JQ(1))  GO TO 1100                                     
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        A1 = A1+COEQUA(JJ)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ),1)/B        
 1100   CONTINUE                                                        
 1200   YF = YY/(4.0+YY)                                                
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = 9.0*P-1.000                                                
        RQUA = COEQ(2)*EXTPY*PICOEF*ANGQ*YG*(YF*PF*AIQ*AMZZ(1)**2-A1)**2
        RQUA = RQUA*POP(1)                                              
        RETURN                                                          
 2000   A1 =0.000                                                       
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        IF(M.EQ.2)  GO TO 2300                                          
        IF(F(2).EQ.1)  GO TO 2200                                       
        DO 2100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JQ(2))  GO TO 2050                                     
        A1 = A1+COEQUA(JJ+3)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ+3),2)/B    
 2050   IF(JJ.GT.JQ(3))  GO TO 2100                                     
        A1 = A1+COEQUA(JJ+ 6)*MATEL(N1,L1,N2,L2,J+5,EXPQUA(JJ+ 6),2)/B  
 2100   CONTINUE                                                        
 2200   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = 9.0*(4.0+5.0*YY)/(4.0+YY)*P-1.000                          
        A(1) = COEQ(3)*POP(2)*EXTPY*YG*(YF*PF*AIQ*AMZZ(2)**2-A1)**2     
 2300   IF(M.EQ.1)  GO TO 2600                                          
        IF(F(3)+F(4).EQ.2)  GO TO 2500                                  
        DO 2400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JQ(4))  GO TO 2350                                     
        A2 = A2+COEQUA(JJ+ 9)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ+ 9),2)/B  
 2350   IF(JJ.GT.JQ(5))  GO TO 2400                                     
        A3 = A3+COEQUA(JJ+12)*MATEL(N1,L1,N2,L2,J+6,EXPQUA(JJ+12),2)/B  
 2400   CONTINUE                                                        
 2500   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)/YY                                                
        PF = 3.0*P+1.000                                                
        A(2) = COEQ(4)*POP(3)*EXTPY*YG*(YF*PF*AIQ*AMZZ(2)**2-A2)**2     
        YF = YY*YY/((1.0+YY)*(9.0+YY))                                  
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = (68.0+77.0*YY)/(4.0+YY)*P-1.000                            
        A(3) = COEQ(5)*POP(3)*EXTPY*YG*(YF*PF*AIQ*AMZZ(2)**2-A3)**2     
 2600   AT = A(1)+A(2)+A(3)                                             
        RQUA = PICOEF*ANGQ*AT                                           
        RETURN                                                          
 3000   A1 = 0.000                                                      
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        A4 = 0.000                                                      
        A5 = 0.000                                                      
        A6 = 0.000                                                      
        IF(M.GT.1)  GO TO 3300                                          
        IF(F(5).EQ.1)  GO TO 3200                                       
        DO 3100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JQ(6))  GO TO 3030                                     
        A1 = A1+COEQUA(JJ+15)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ+15),3)/B  
 3030   IF(JJ.GT.JQ(7))  GO TO 3060                                     
        A1 = A1+COEQUA(JJ+18)*MATEL(N1,L1,N2,L2,J+5,EXPQUA(JJ+18),3)/B  
 3060   IF(JJ.GT.JQ(8))  GO TO 3100                                     
        A1 = A1+COEQUA(JJ+21)*MATEL(N1,L1,N2,L2,J+6,EXPQUA(JJ+21),3)/B  
 3100   CONTINUE                                                        
 3200   YF = YY*(9.0+YY)/((1.0+YY)*(4.0+YY))                            
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = (729.0+1134.0*YY+277.0*YY*YY)/(9.0+YY)**2*P-1.000          
        A(1) = COEQ(6)*POP(4)*EXTPY*YG*(YF*PF*AIQ*AMZZ(3)**2-A1)**2     
 3300   IF(M.NE.0.AND.M.NE.2)  GO TO 3600                               
        IF(F(6)+F(7).EQ.2)  GO TO 3500                                  
        DO 3400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JQ(9))  GO TO 3325                                     
        A2 = A2+COEQUA(JJ+24)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ+24),3)/B  
 3325   IF(JJ.GT.JQ(10))  GO TO 3350                                    
        A2 = A2+COEQUA(JJ+27)*MATEL(N1,L1,N2,L2,J+5,EXPQUA(JJ+27),3)/B  
 3350   IF(JJ.GT.JQ(11))  GO TO 3375                                    
        A3 = A3+COEQUA(JJ+30)*MATEL(N1,L1,N2,L2,J+6,EXPQUA(JJ+30),3)/B  
 3375   IF(JJ.GT.JQ(12))  GO TO 3400                                    
        A3 = A3+COEQUA(JJ+33)*MATEL(N1,L1,N2,L2,J+7,EXPQUA(JJ+33),3)/B  
 3400   CONTINUE                                                        
 3500   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)/YY                                                
        PF = (27.0+11.0*YY)/(9.0+YY)*P+1.000                            
        A(2) = COEQ(7)*POP(5)*EXTPY*YG*(YF*PF*AIQ*AMZZ(3)**2-A2)**2     
        YF = YY*YY/((1.0+YY)*(4.0+YY))                                  
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = (1377.0+1944.0*YY+439.0*YY*YY)/(9.0+YY)**2*P-1.000         
        A(3) = COEQ(8)*POP(5)*EXTPY*YG*(YF*PF*AIQ*AMZZ(3)**2-A3)**2     
 3600   IF(M.NE.0.AND.M.NE.3)  GO TO 3900                               
        IF(F(8)+F(9)+F(10).EQ.3)  GO TO 3800                            
        DO 3700 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JQ(13))  GO TO 3630                                    
        A4 = A4+COEQUA(JJ+36)*MATEL(N1,L1,N2,L2,J+4,EXPQUA(JJ+36),3)/B  
 3630   IF(JJ.GT.JQ(14))  GO TO 3660                                    
        A5 = A5+COEQUA(JJ+39)*MATEL(N1,L1,N2,L2,J+6,EXPQUA(JJ+39),3)/B  
 3660   IF(JJ.GT.JQ(15))  GO TO 3700                                    
        A6 = A6+COEQUA(JJ+42)*MATEL(N1,L1,N2,L2,J+8,EXPQUA(JJ+42),3)/B  
 3700   CONTINUE                                                        
 3800   YF = YY/(9.0+YY)                                                
        A(4) = COEQ(9)*POP(6)*EXTPY*(YF*P*AIQ*AMZZ(3)**2-A4)**2         
        YF = YY*YY/((1.0+YY)*(4.0+YY))                                  
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = (63.0+47.0*YY)/(9.0+YY)*P+1.000                            
        A(5) = COEQ(10)*POP(6)*EXTPY*YG*(YF*PF*AIQ*AMZZ(3)**2-A5)**2    
        YF = YY**3/((1.0+YY)*(4.0+YY)*(16.0+YY))                        
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)*(16.0+YY)/YY**4                 
        PF = (10773.0+14580.0*YY+3167.0*YY*YY)/(9.0+YY)**2*P-5.000      
        A(6) = COEQ(11)*POP(6)*EXTPY*YG*(YF*PF*AIQ*AMZZ(3)**2-A6)**2    
 3900   AT = A(1)+A(2)+A(3)+A(4)+A(5)+A(6)                              
        RQUA = PICOEF*ANGQ*AT                                           
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION ROCT(N1,L1,N2,L2,N,M,Y)                                
        DIMENSION A(6),JO(15)                                           
C  ***  OCTUPOLE PENETRATION ROUTINE                                    
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC016/COEOCT(45),EXPOCT(45)                             
        COMMON/LOC020/F(10)                                             
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC025/COEO(11)                                          
        COMMON/LOC030/POP(6),DUM(24)                                    
        COMMON/LOC031/JM(10),JD(14),JQ(15),KO(15),IYC,IJ(4),YJ(4),JJ1(4)
        DOUBLE PRECISION ETPY
        INTEGER F                                                       
        REAL MATEL                                                      
        if (Y > log(huge(1.0))/(2.0*PI)) ETPY = 1e100
        if (Y < log(huge(1.0))/(2.0*PI)) ETPY = EXP(2.000*PI*Y)                                  
C       ETPY = EXP(2.000*PI*Y)                                          
        write(4,*) Y, PI, "Y VALUE"
        write(4,*) " ETPY HERE", ETPY
        flush(4)
        EXTPY = ETPY/(ETPY-1.000)                                       
        YY = Y*Y                                                        
        P = EXP(Y*(2.000*ATAN(Y/FLOAT(N))-PI))                          
        DO 500 I=1,6                                                    
        A(I)=0.000                                                      
  500   CONTINUE                                                        
        DO 550 JJ=1,15                                                  
        JO(JJ) = KO(JJ)                                                 
  550   CONTINUE                                                        
        IF(IJ(4).EQ.0)  GO TO 900                                       
        DO 600 JJ=1,15                                                  
        IF(Y.GT.YJ(4))  JO(JJ) = MIN0(KO(JJ),JJ1(4))                    
  600   CONTINUE                                                        
  900   GO TO (1000,2000,3000),N                                        
 1000   A1 = 0.000                                                      
        IF(F(1).EQ.1)  GO TO 1200                                       
        DO 1100 JJ=1,3                                                  
        IF(JJ.GT.JO(1))  GO TO 1100                                     
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        A1 = A1+COEOCT(JJ)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ),1)/B        
 1100   CONTINUE                                                        
 1200   YF = YY*(3.0+2.0*YY)/((4.0+YY)*(9.0+YY))                        
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = 15.0*(1.0+YY)/(3.0+2.0*YY)*P-1.000                         
        ROCT = COEO(2)*PICOEF*ANGO*EXTPY*YG*(YF*PF*AIO*AMZZ(1)**3-A1)**2
        ROCT = ROCT*POP(1)                                              
        RETURN                                                          
 2000   A1 =0.000                                                       
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        IF(M.EQ.2)  GO TO 2300                                          
        IF(F(2).EQ.1)  GO TO 2200                                       
        DO 2100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JO(2))  GO TO 2050                                     
        A1 = A1+COEOCT(JJ+3)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ+3),2)/B    
 2050   IF(JJ.GT.JO(3))  GO TO 2100                                     
        A1 = A1+COEOCT(JJ+ 6)*MATEL(N1,L1,N2,L2,J+6,EXPOCT(JJ+ 6),2)/B  
 2100   CONTINUE                                                        
 2200   YF = YY*(6.0+YY)/((1.0+YY)*(9.0+YY))                            
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = 15.0*(2.0+3.0*YY)/(6.0+YY)*P-1.000                         
        IF(COEO(3) * POP(2) * EXTPY * YG < 1.0e-10) GO TO 2299
        A(1) = COEO(3)*POP(2)*EXTPY*YG*(YF*PF*AIO*AMZZ(2)**3-A1)**2     
 2299   A(1) = 0.0
 2300   IF(M.EQ.1)  GO TO 2600                                          
        IF(F(3)+F(4).EQ.2)  GO TO 2500                                  
        DO 2400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JO(4))  GO TO 2350                                     
        A2 = A2+COEOCT(JJ+ 9)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ+ 9),2)/B  
 2350   IF(JJ.GT.JO(5))  GO TO 2400                                     
        A3 = A3+COEOCT(JJ+12)*MATEL(N1,L1,N2,L2,J+7,EXPOCT(JJ+12),2)/B  
 2400   CONTINUE                                                        
 2500   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = 3.0*P+1.000                                                
        A(2) = COEO(4)*POP(3)*EXTPY*YG*(YF*PF*AIO*AMZZ(2)**3-A2)**2     
        YF = YY*YY*(68.0+13.0*YY)/((1.0+YY)*(9.0+YY)*(16.0+YY))         
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)*(16.0+YY)/YY**4                 
        PF = 5.0*(116.0+149.0*YY)/(68.0+13.0*YY)*P-1.000                
        A(3) = COEO(5)*POP(3)*EXTPY*YG*(YF*PF*AIO*AMZZ(2)**3-A3)**2     
 2600   AT = A(1)+A(2)+A(3)                                             
        ROCT = PICOEF*ANGO*AT                                           
        RETURN                                                          
 3000   A1 = 0.000                                                      
        A2 = 0.000                                                      
        A3 = 0.000                                                      
        A4 = 0.000                                                      
        A5 = 0.000                                                      
        A6 = 0.000                                                      
        IF(M.GT.1)  GO TO 3300                                          
        IF(F(5).EQ.1)  GO TO 3200                                       
        DO 3100 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JO(6))  GO TO 3030                                     
        A1 = A1+COEOCT(JJ+15)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ+15),3)/B  
 3030   IF(JJ.GT.JO(7))  GO TO 3060                                     
        A1 = A1+COEOCT(JJ+18)*MATEL(N1,L1,N2,L2,J+6,EXPOCT(JJ+18),3)/B  
 3060   IF(JJ.GT.JO(8))  GO TO 3100                                     
        A1 = A1+COEOCT(JJ+21)*MATEL(N1,L1,N2,L2,J+7,EXPOCT(JJ+21),3)/B  
 3100   CONTINUE                                                        
 3200   YF = YY*(27.00+2.000*YY)/((1.0+YY)*(4.0+YY))                    
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = 2.5*(405.0+900.0*YY+254.0*YY*YY)/                          
     1       ((9.0+YY)*(27.00+2.000*YY))*P-1.000                        
        A(1) = COEO(6)*POP(4)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A1)**2     
 3300   IF(M.NE.0.AND.M.NE.2)  GO TO 3600                               
        IF(F(6)+F(7).EQ.2)  GO TO 3500                                  
        DO 3400 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JO(9))  GO TO 3325                                     
        A2 = A2+COEOCT(JJ+24)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ+24),3)/B  
 3325   IF(JJ.GT.JO(10))  GO TO 3350                                    
        A2 = A2+COEOCT(JJ+27)*MATEL(N1,L1,N2,L2,J+6,EXPOCT(JJ+27),3)/B  
 3350   IF(JJ.GT.JO(11))  GO TO 3375                                    
        A3 = A3+COEOCT(JJ+30)*MATEL(N1,L1,N2,L2,J+7,EXPOCT(JJ+30),3)/B  
 3375   IF(JJ.GT.JO(12))  GO TO 3400                                    
        A3 = A3+COEOCT(JJ+33)*MATEL(N1,L1,N2,L2,J+8,EXPOCT(JJ+33),3)/B  
 3400   CONTINUE                                                        
 3500   YF = YY*(9.0+2.0*YY)/((1.0+YY)*(4.0+YY))                        
        YG = (1.0+YY)*(4.0+YY)/(YY*YY)                                  
        PF = (27.0+13.0*YY)/(9.0+2.0*YY)*P+1.000                        
        A(2) = COEO(7)*POP(5)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A2)**2     
        YF = YY*YY*(153.0+13.0*YY)/((1.0+YY)*(4.0+YY)*(16.0+YY))        
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)*(16.0+YY)/YY**4                 
        PF = 5.0*(2349.0+3744.0*YY+947.0*YY*YY)/                        
     1       ((9.0+YY)*(153.0+13.0*YY))*P-1.000                         
        A(3) = COEO(8)*POP(5)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A3)**2     
 3600   IF(M.NE.0.AND.M.NE.3)  GO TO 3900                               
        IF(F(8)+F(9)+F(10).EQ.3)  GO TO 3800                            
        DO 3700 JJ=1,3                                                  
        J = JK(JJ)                                                      
        B = Y**NE(JJ)                                                   
        IF(JJ.GT.JO(13))  GO TO 3630                                    
        A4 = A4+COEOCT(JJ+36)*MATEL(N1,L1,N2,L2,J+5,EXPOCT(JJ+36),3)/B  
 3630   IF(JJ.GT.JO(14))  GO TO 3660                                    
        A5 = A5+COEOCT(JJ+39)*MATEL(N1,L1,N2,L2,J+7,EXPOCT(JJ+39),3)/B  
 3660   IF(JJ.GT.JO(15))  GO TO 3700                                    
        A6 = A6+COEOCT(JJ+42)*MATEL(N1,L1,N2,L2,J+9,EXPOCT(JJ+42),3)/B  
 3700   CONTINUE                                                        
 3800   YF = YY/(1.0+YY)                                                
        YG = (1.0+YY)/YY                                                
        PF = 2.0*P+1.000                                                
        A(4) = COEO(9)*POP(6)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A4)**2     
        YF = YY*YY/((1.0+YY)*(4.0+YY))                                  
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)/YY**3                           
        PF = (63.0+47.0*YY)/(9.0+YY)*P+1.000                            
        A(5) = COEO(10)*POP(6)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A5)**2    
        YF = YY**3*(11.0+YY)/((1.0+YY)*(4.0+YY)*(16.0+YY)*(25.0+YY))    
        YG = (1.0+YY)*(4.0+YY)*(9.0+YY)*(16.0+YY)*(25.0+YY)/YY**5       
        PF = (1251.0+1850.0*YY+439.0*YY*YY)/((9.0+YY)*(11.0+YY))*P-1.000
        A(6) = COEO(11)*POP(6)*EXTPY*YG*(YF*PF*AIO*AMZZ(3)**3-A6)**2    
 3900   AT = A(1)+A(2)+A(3)+A(4)+A(5)+A(6)                              
        ROCT = PICOEF*ANGO*AT                                           
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RDIPU(N1,L1,N2,L2,N,ENEM,Y,MM)                         
C  ***  DIPOLE UNPENETRATED RATES ROUTINE                               
C
C       This routine is documented in Akylas PhD thesis, Appendix A,
C       Table A.4
C       It is for fully occupied electron shells
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC026/COED(4)                                           
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC030/POP(6),D(24)                                      
        DOUBLE PRECISION MATELU, AID, AIDSQ, COEFF, ENEM, ALFA
C
C       Here we catch 1,2,3 cases
C       N = 0 means we just want the radiative rates
        IF(N.GT.0) GO TO 1000                                           
C
C       Compute the muonic matrix element
        AID = MATELU(N1,L1,N2,L2,1)                                     
C
C       Square it
        AIDSQ = AID*AID                                                 
        L = (L1+L2+1)/2                                                 
C
C       L seems to be lambda in the thesis
C       lambda is l1 + l2 + 1 / 2 for odd multipoles
C       and l1 + l2 / 2 for even multipoles
        ANGD = FLOAT((2*L2+1)*L)/FLOAT((2*L-1)*(2*L+1))                 
C
C       Radiative dipole contribution
        RDIPU = COED(1)*COEFF/(Z*Z*AMASSM*AMASSM)*(ENEM/ALFA)**3        
     1          *ANGD*AIDSQ                                             
C       if(abs(l2 - l1) == 1)write(4,*) RDIPU, n1, l1, n2, l2
        RETURN                                                          
C
C       Calculate the exponential
 1000   EPIY = EXP(PI*Y)                                                
C
C       EXPIY = 1/(2*SINH(PIY))
        EXPIY = 1.000/(EPIY-1.000/EPIY)                                 
        YY = Y*Y                                                        
C
C       Calculate the multiplier which we need to multiply by matrix element and
C       another bunch of shit
C       PDJ TODO figure out where the 4 comes from. Is it a typo from table A.4?
        P2 = EXP(Y*(4.000*ATAN(Y/FLOAT(N))-PI))                         
C      write(4,*) P2, "p", N1, L1, N2, L2, Y
        M = MM+1                                                        
C
C       This tells us if we are K, L, or M shell
        GO TO (2000,3000,4000),N                                        
C
C       Arrive at the K-shell
C       Only 1 case for this
 2000   YF = YY/(1.000+YY)                                              
C
C       COED:
C       POP: Population of subshell
C       PICOEF: pi x fine structure x speed of light / bohr radius
C       ANGD: Angular part
C       P2: Multiplier for shell
C       EXPIY: Sinh(pi*y)
C       YF: Prefactor for P 
C       AMZZ: Effective charge/total charge
        RDIPU = COED(2)*POP(1)*PICOEF*ANGD*P2*EXPIY*YF*AIDSQ*AMZZ(1)**2 
        RETURN                                                          
C
C       2s, 2p and total are 3 separate cases
 3000   GO TO (3100,3200,3300),M                                        
C
C       Total L-shell
 3100   YF = YY*(4.0+3.0*YY)*(4.0+5.0*YY)/(4.0+YY)**3*(POP(2)+3.000*    
     1  POP(3))/4.000                                                   
        GO TO 3400                                                      
C
C       Prefactor for 2s
 3200   YF = 4.000*YY*(1.0+YY)/(4.0+YY)**2*POP(2)                       
        GO TO 3400                                                      
C
C       Prefactor for 2p
 3300   YF = YY*YY*(12.0+11.0*YY)/(4.0+YY)**3*POP(3)                    
 3400   RDIPU = COED(3)*PICOEF*ANGD*P2*EXPIY*YF*AIDSQ*AMZZ(2)**2        
        RETURN                                                          
C 
C       3s, 3p, 3d and total are 4 separate cases
 4000   GO TO (4100,4200,4300,4400),M                                   
C
C       Total M-shell
 4100   YF = YY*(81.0+78.0*YY+13.0*YY*YY)*(81.0+126.0*YY+29.0*YY*YY)    
     1  /(9.0+YY)**5*(POP(4)+3.000*POP(5)+5.000*POP(6))/9.000           
        GO TO 4500                                                      
C
C       3s shell
 4200   YF = YY*(1.0+YY)*(27.0+7.0*YY)**2/(9.0+YY)**4*POP(4)            
        GO TO 4500                                                      
C
C       3p shell
 4300   YF = 8.000*YY*YY*(81.0+96.0*YY+19.0*YY*YY)/(9.0+YY)**4*POP(5)   
        GO TO 4500                                                      
C
C       3d shell
 4400   YF = 16.000*YY**3*(45.0+11.0*YY)*(1.0+YY)/(9.0+YY)**5*POP(6)    
 4500   RDIPU = COED(4)*PICOEF*ANGD*P2*EXPIY*YF*AIDSQ*AMZZ(3)**2        
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RQUAU(N1,L1,N2,L2,N,ENEM,Y,M)                          
        DIMENSION A(6)                                                  
C  ***  QUADRUPOLE UNPENETRATES RATES ROUTINE                           
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC024/COEQ(11)                                          
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC030/POP(6),D(24)                                      
        REAL MATELU                                                     
C
C       This is the radiative rate
        IF(N.GT.0) GO TO 1000                                           
C
C       Table A.1
C       Calculate matrix element
        AIQ = MATELU(N1,L1,N2,L2,2)                                     
        AIQSQ = AIQ*AIQ                                                 
C
C       Here we have an even multipole
        L = (L1+L2)/2                                                   
C
C       Chooses correct coefficient based on delta L
C       1.5 is for if deltaL = 2
        ANGQ=1.500                                                      
C       1 is for if deltaL = 0
        IF(L2.EQ.L1)  ANGQ=1.000                                        
C
C       Angular part
        ANGQ = ANGQ*FLOAT((2*L2+1)*L*(L+1))/                            
     1         FLOAT((2*L-1)*(2*L+1)*(2*L+3))                           
        RQUAU = COEQ(1)*COEFF/(Z*AMASSM)**4*(ENEM/ALFA)**5*ANGQ*AIQSQ   
C       write(4,*) RQUAU,"QUAD RADIATIVE RATE", n1, l1, n2, l2
        RETURN                                                          

 1000   ETPY = EXP(2.000*PI*Y)                                          
        EXTPY = ETPY/(ETPY-1.000)                                       
        YY = Y*Y                                                        
        P = EXP(Y*(2.000*ATAN(Y/FLOAT(N))-PI))                          

        DO 1100 I=1,6                                                   
        A(I)=0.000                                                      
 1100   CONTINUE                                                        

        ! Electron shell
        GO TO (2000,3000,4000),N                                        
C
C       K shell Auger rate
 2000   YF = (1.0+YY)/(4.0+YY)                                          
        PF = (9.0*P-1.0)**2                                             
C       write(4,*) PF, YF
        RQUAU = COEQ(2)*POP(1)*PICOEF*ANGQ*PF*EXTPY*YF*AIQSQ*AMZZ(1)**4 
C       write(4,*) "QUAD K AUGER", RQUAU, N1, L1, N2, L2
        RETURN                                                          
C       
C       M being zero means we are doing full shell
C       so we do it all
C       If M = 2, we are only considering p-shells
 3000   IF(M.EQ.2)  GO TO 3100                                          
        YF = (4.0+YY)/(1.0+YY)                                          
        PF = (9.0*(4.0+5.0*YY)/(4.0+YY)*P-1.000)**2                     
        A(1) = COEQ(3)*POP(2)*YF*PF                                     
C       
C       If M = 1, then we are only considering s shells
 3100   IF(M.EQ.1)  GO TO 3200                                          
        YF = YY/(1.0+YY)                                                
        PF = (3.0*P+1.000)**2                                           
        A(2) = COEQ(4)*POP(3)*YF*PF                                     
C
C
        YF = YY*(4.0+YY)/((1.0+YY)*(9.0+YY))                            
        PF = ((68.0+77.0*YY)/(4.0+YY)*P-1.000)**2                       
        A(3) = COEQ(5)*POP(3)*YF*PF                                     
 3200   AT = A(1)+A(2)+A(3)                                             
C       write(4,*) "AT", AT
        RQUAU = PICOEF*ANGQ*AT*EXTPY*AIQSQ*AMZZ(2)**4                   

        RETURN                                                          
 4000   IF(M.GT.1)  GO TO 4100                                          
C
C
        YF = (9.0+YY)**2/((1.0+YY)*(4.0+YY))                            
        PF = ((729.0+1134.0*YY+277.0*YY*YY)/(9.0+YY)**2*P-1.000)**2     
        A(1) = COEQ(6)*POP(4)*YF*PF                                     
C
 4100   IF(M.NE.0.AND.M.NE.2)  GO TO 4200                               
C
C
        YF = YY/(1.0+YY)                                                
        PF = ((27.0+11.0*YY)/(9.0+YY)*P+1.000)**2                       
        A(2) = COEQ(7)*POP(5)*YF*PF                                     
C
C
        YF = YY*(9.0+YY)/((1.0+YY)*(4.0+YY))                            
        PF = ((1377.0+1944.0*YY+439.0*YY*YY)/(9.0+YY)**2*P-1.000)**2    
        A(3) = COEQ(8)*POP(5)*YF*PF                                     
C
C
 4200   IF(M.NE.0.AND.M.NE.3)  GO TO 4300                               
        YF = YY*YY/(9.0+YY)**2                                          
        PF = P**2                                                       
        A(4) = COEQ(9)*POP(6)*YF*PF                                     
C
C
        YF = YY*YY/((1.0+YY)*(4.0+YY))                                  
        PF = ((63.0+47.0*YY)/(9.0+YY)*P+1.000)**2                       
        A(5) = COEQ(10)*POP(6)*YF*PF                                    
C
        YF = YY*YY*(9.0+YY)/((1.0+YY)*(4.0+YY)*(16.0+YY))               
        PF = ((10773.0+14580.0*YY+3167.0*YY*YY)/(9.0+YY)**2*P-5.000)**2 
        A(6) = COEQ(11)*POP(6)*YF*PF                                    
 4300   AT = A(1)+A(2)+A(3)+A(4)+A(5)+A(6)                              
        RQUAU = PICOEF*ANGQ*AT*EXTPY*AIQSQ*AMZZ(3)**4                   
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION ROCTU(N1,L1,N2,L2,N,ENEM,Y,M)                          
        DIMENSION A(6)                                                  
C  ***  OCTUPOLE UNPENETRATED RATES ROUTINE                             
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC025/COEO(11)                                          
        COMMON/LOC021/ANGD,ANGQ,ANGO                                    
        COMMON/LOC022/AID,AIQ,AIO,AIDSQ,AIQSQ,AIOSQ                     
        COMMON/LOC030/POP(6),D(24)                                      
        REAL MATELU                                                     
C       Skip down to Auger transitions 
        IF(N.GT.0) GO TO 1000                                           
C
C       This radiative rate can be found in Table A.1 of Akylas'
C       PhD thesis
C       Calculate matrix element and square it
        AIO = MATELU(N1,L1,N2,L2,3)                                     
        AIOSQ = AIO*AIO                                                 
C
C       This is equivalent to lambda. Octupole is odd
        L = (L1+L2+1)/2                                                 
        ANGO = 2.500                                                    
        IF(IABS(L1-L2).EQ.1)  ANGO = 1.500                              
C
C       Angular bit of radiative transition
        ANGO = ANGO*FLOAT((2*L2+1)*(L-1)*L*(L+1))/                      
     1         FLOAT((2*L-3)*(2*L-1)*(2*L+1)*(2*L+3))                   
C
C       
        ROCTU = COEO(1)*COEFF/(Z*AMASSM)**6*(ENEM/ALFA)**7*ANGO*AIOSQ   
        RETURN                                                          

 1000   ETPY = EXP(2.000*PI*Y)                                          
        EXTPY = ETPY/(ETPY-1.000)                                       
        YY = Y*Y                                                        
        P = EXP(Y*(2.000*ATAN(Y/FLOAT(N))-PI))                          
        DO 1100 I=1,6                                                   
        A(I) = 0.000                                                    
 1100   CONTINUE                                                        
C
C       Energy level of emitted electron
        GO TO (2000,3000,4000),N                                        
C
C       K shell electron emission
 2000   YF = (1.0+YY)*(3.0+2.0*YY)**2/(YY*(4.0+YY)*(9.0+YY))            
        PF = (15.0*(1.0+YY)/(3.0+2.0*YY)*P-1.000)**2                    
        ROCTU = COEO(2)*POP(1)*PICOEF*ANGO*PF*EXTPY*YF*AIOSQ*AMZZ(1)**6 
C       write(4,*) ROCTU, n1, l1, n2, l2
        RETURN                                                          

C
C       If M = 0, we do a full shell for L
 3000   IF(M.EQ.2)  GO TO 3100                                          
        YF = (4.0+YY)*(6.0+YY)**2/(YY*(1.0+YY)*(9.0+YY))                
        PF = (15.0*(2.0+3.0*YY)/(6.0+YY)*P-1.000)**2                    
        A(1) = COEO(3)*POP(2)*YF*PF                                     
C
 3100   IF(M.EQ.1)  GO TO 3200                                          
        YF = (4.0+YY)/(1.0+YY)                                          
        PF = (3.0*P+1.000)**2                                           
        A(2) = COEO(4)*POP(3)*YF*PF                                     
        YF = (4.0+YY)*(68.0+13.0*YY)**2/((1.0+YY)*(9.0+YY)*(16.0+YY))   
        PF = (5.0*(116.0+149.0*YY)/(68.0+13.0*YY)*P-1.000)**2           
        A(3) = COEO(5)*POP(3)*YF*PF                                     
 3200   AT = A(1)+A(2)+A(3)                                             
        ROCTU = PICOEF*ANGO*AT*EXTPY*AIOSQ*AMZZ(2)**6                   
C       write(4,*) ROCTU, n1, l1, n2, l2
C       write(4,*) A(1), A(2), A(3)
        RETURN                                                          

C       Now M shell
 4000   IF(M.GT.1)  GO TO 4100                                          
        YF = (9.0+YY)*(27.00+2.000*YY)**2/(YY*(1.0+YY)*(4.0+YY))        
        PF = (2.5*(405.0+900.0*YY+254.0*YY*YY)/                         
     1  ((9.0+YY)*(27.00+2.000*YY))*P-1.000)**2                         
        A(1) = COEO(6)*POP(4)*YF*PF                                     
C 
 4100   IF(M.NE.0.AND.M.NE.2)  GO TO 4200                               
C       3p with deltaL = 3
        YF = (9.0+2.0*YY)**2/((1.0+YY)*(4.0+YY))                        
        PF = ((27.0+13.0*YY)/(9.0+2.0*YY)*P+1.000)**2                   
        A(2) = COEO(7)*POP(5)*YF*PF                                     
C
C    
        YF = (9.0+YY)*(153.0+13.0*YY)**2/((1.0+YY)*(4.0+YY)*(16.0+YY))  
        PF = (5.0*(2349.0+3744.0*YY+947.0*YY*YY)/                       
     1  ((9.0+YY)*(153.0+13.0*YY))*P-1.000)**2                          
        A(3) = COEO(8)*POP(5)*YF*PF                                     
 4200   IF(M.NE.0.AND.M.NE.3)  GO TO 4300                               
C
C
        YF = YY/(1.0+YY)                                                
        PF = (2.0*P+1.000)**2                                           
        A(4) = COEO(9)*POP(6)*YF*PF                                     
C
C
        YF = YY*(9.0+YY)/((1.0+YY)*(4.0+YY))                            
        PF = ((63.0+47.0*YY)/(9.0+YY)*P+1.000)**2                       
        A(5) =COEO(10)*POP(6)*YF*PF                                     
C
C 
        YF = YY*(9.0+YY)*(11.0+YY)**2/                                  
     1  ((1.0+YY)*(4.0+YY)*(16.0+YY)*(25.0+YY))                         
        PF = ((1251.0+1850.0*YY+439.0*YY*YY)/                           
     1  ((9.0+YY)*(11.0+YY))*P-1.0000)**2                               
        A(6) = COEO(11)*POP(6)*YF*PF                                    
 4300   AT = A(1)+A(2)+A(3)+A(4)+A(5)+A(6)                              
        ROCTU = PICOEF*ANGO*AT*EXTPY*AIOSQ*AMZZ(3)**6                   
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE POPJ(L1,L2,LL,P)                                     
        DIMENSION P(4)                                                  
        LI=LL+1                                                         
C  ***  FINDS RELATIVE POPULATION OF ALL POSSIBLE J-STATES              
        L = MAX0(L1,L2)                                                 
        P(3) = 0.000                                                    
        P(4) = 0.000                                                    
        GO TO (1000,2000,3000,4000),LI                                  
 1000   P(1) = FLOAT(L+1)/FLOAT(2*L+1)                                  
        P(2) = 1.000-P(1)                                               
        RETURN                                                          
 2000   P(1) = FLOAT((L+1)*(2*L-1))/FLOAT(4*L*L-1)                      
        P(2) = 1.000/FLOAT(4*L*L-1)                                     
        P(3) = 1.000-(P(1)+P(2))                                        
        RETURN                                                          
 3000   IF(L1.EQ.L2)  GO TO 3500                                        
        P(1) = FLOAT(L+1)/FLOAT(2*L+1)                                  
        P(2) = 2.000/FLOAT((2*L-3)*(2*L+1))                             
        P(3) = 1.000-(P(1)+P(2))                                        
        RETURN                                                          
 3500   D = FLOAT((2*L+1)**2)                                           
        P(1) = FLOAT((L+2)*(2*L-1))/D                                   
        P(2) = FLOAT((L-1)*(2*L+3))/D                                   
        P(3) = 3.000/D                                                  
        P(4) = P(3)                                                     
        RETURN                                                          
 4000   IF(IABS(L1-L2).EQ.1)  GO TO 4500                                
        P(1) = FLOAT(L+1)/FLOAT(2*L+1)                                  
        P(2) = 3.000/FLOAT((2*L-5)*(2*L+1))                             
        P(3) = 1.000-(P(1)+P(2))                                        
        RETURN                                                          
 4500   D = FLOAT(4*L*L-1)                                              
        P(1) = FLOAT((2*L-3)*(L+2))/D                                   
        P(2) = 5.000/D                                                  
        P(3) = 6.000/D                                                  
        P(4) = 1.000-(P(1)+P(2)+P(3))                                   
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  SUBROUTINE CHECK(N)                                             
        DIMENSION AA(3,20)                                              
        COMMON/LOC007/IC                                                
C  ***  NUMERICAL ACCURACY CONTROL -- COMPUTES DIAGONAL M.E.            
        COMMON/LOC008/IR,IW,IP,IPRINT                                   
        COMMON/LOC027/LL(20)                                            
        REAL MATELU                                                     
        IF(IC.LE.0)  RETURN                                             
        GO TO (1000,2000,3000),IC                                       
 1000   A = (MATELU(N,0,N,0,1)/RAV(N,0,1))**2                           
        WRITE(IW,1100)N,A                                               
 1100   FORMAT(/29H *** ACCURACY CONTROL AT N = ,I2,8H, DIP = ,         
     1  F13.9,22H (SHOULD BE UNITY) ***)                                
     0  RETURN                                                          
 2000   A = (MATELU(N,0,N,0,1)/RAV(N,0,1))**2                           
        B = (MATELU(N,0,N,0,2)/RAV(N,0,2))**2                           
        C = (MATELU(N,0,N,0,3)/RAV(N,0,3))**2                           
        WRITE(IW,2100)N,A,B,C                                           
 2100   FORMAT(/29H *** ACCURACY CONTROL AT N = ,I2,8H, DIP = ,F13.9,   
     1  8H, QUA = ,F13.9,8H, OCT = ,F13.9,22H (SHOULD BE UNITY) ***)    
     0  RETURN                                                          
 3000   DO 3100 I = 1,3                                                 
        DO 3100 J = 1,N                                                 
        I1 = I                                                          
        J1 = J-1                                                        
        AA(I,J) = (MATELU(N,J1,N,J1,I1)/RAV(N,J1,I1))**2                
 3100   CONTINUE                                                        
        WRITE(IW,3200)N,(LL(J),(AA(I,J),I=1,3),J=1,N)                   
 3200   FORMAT(/29H *** ACCURACY CONTROL AT N = ,I2,16H (SHOULD BE UNIT 
     1  ,6HY) ***/41H  L    DIPOLES    QUADRUPOLES   OCTUPOLES/         
     2  (1X,I2,3F13.9))                                                 
     0  RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  REAL FUNCTION MATEL(N1,L1,N2,L2,L,A,N)                          
C  ***  GENERAL DIMENSIONLESS MUONIC MATRIX ELEMENT FOR PENETRATION     
        DOUBLE PRECISION F,P,A1,A2,S1,S2                                
        COMMON/LOC012/AMZZ(3)                                           
        COMMON/LOC009/F(60),FD                                          
        AN = FLOAT(N1+N2)                                               
C
C       A is the exponent fitted for penetration
        P = 1.000D00 + DBLE(1.000E-03*A*FLOAT(N1*N2)/AN)                
        A1 = -2.000D00*DBLE(FLOAT(N2)/AN)/P                             
        A2 = -2.000D00*DBLE(FLOAT(N1)/AN)/P                             
        M1 = N1-L1                                                      
        M2 = N2-L2                                                      
        M3 = N1+L1+1                                                    
        M4 = N2+L2+1                                                    
        M5 = 2*L1+2                                                     
        M6 = 2*L2+2                                                     
        MM = L1+L2+L+3                                                  
        S1 = 0.000D00                                                   
        DO 200 I1 = 1,M1                                                
        K1 = I1-1                                                       
        S2 = 0.000D00                                                   
        DO 100 I2 = 1,M2                                                
        K2 = I2-1                                                       
        LA = MM+K1+K2                                                   
        LB = M2-K2                                                      
        LC = M6+K2                                                      
        S2 = S2 + A2**K2*F(LA)/(F(LB)*F(LC)*F(I2))                      
  100   CONTINUE                                                        
        LD = M1-K1                                                      
        LE = M5+K1                                                      
        S1 = S1 + S2*A1**K1/(F(LD)*F(LE)*F(I1))                         
  200   CONTINUE                                                        
        AQ = 2.000/(AN*SNGL(P))                                         
        T1 = SNGL(S1*DSQRT(F(M1)*F(M2)/DBLE(AQ**(L-1))*F(M3)*F(M4)))    
        MATEL = T1*(AQ*FLOAT(N1))**(L+L2+1)*(AQ*FLOAT(N2))**(L+L1+1)    
     1          *0.500**(L+1)*FD**(L+1)/SQRT(AQ**(L-1))*(AMZZ(N))**L    
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
C       EQUATION 7 OF VOGEL PRA 1980
     0  REAL FUNCTION MATELU(N1,L1,N2,L2,L)                             
C  ***  GENERAL DIMENSIONLESS MUONIC MATRIX ELEMENT FOR NONPENETRATION  
        DOUBLE PRECISION F,A1,A2,S1,S2                                  
        COMMON/LOC009/F(60),FD                                          
C       PDJ added dble here
        AN = DBLE(N1+N2)                                               
        A1 = -2.000D00*DBLE(FLOAT(N2)/AN)                               
        A2 = -2.000D00*DBLE(FLOAT(N1)/AN)                               
        M1 = N1-L1                                                      
        M2 = N2-L2                                                      
        M3 = N1+L1+1                                                    
        M4 = N2+L2+1                                                    
        M5 = 2*L1+2                                                     
        M6 = 2*L2+2                                                     
C
C       L is the multipolarity
        MM = L1+L2+L+3                                                  
        S1 = 0.000D00                                                   
C  
C       Here we loop over l1 to n1
        DO 200 I1 = 1,M1                                                
C
C       Calculate the index
        K1 = I1-1                                                       
        S2 = 0.000D00                                                   
C
C       Loop over l2 to n2
        DO 100 I2 = 1,M2                                                
C
C       Calculate the index
        K2 = I2-1                                                       
C
C       
        LA = MM+K1+K2                                                   
        LB = M2-K2                                                      
        LC = M6+K2                                                      
C
C       F is factorial values, F(LA) = (LA-1)!
C       PDJ added dble here
        S2 = S2 + DBLE(A2**K2*F(LA)/(F(LB)*F(LC)*F(I2)))
C
C
  100   CONTINUE                                                        
        LD = M1-K1                                                      
        LE = M5+K1                                                      
        S1 = S1 + S2*A1**K1/(F(LD)*F(LE)*F(I1))                         
        write(4,*) "SUM", S1
C
  200   CONTINUE                                                        
        AQ = 2.000/AN                                                   
        T1 = (S1*DSQRT(F(M1)*F(M2)/DBLE(AQ**(L-1))*F(M3)*F(M4)))    
C
        MATELU = T1*(AQ*FLOAT(N1))**(L+L2+1)*(AQ*FLOAT(N2))**(L+L1+1)   
     1          *0.500**(L+1)*FD**(L+1)/SQRT(AQ**(L-1))                 
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION BETA(L1,J1,L2,J2,L)                                    
        A1 = 0.250*FLOAT(J1*(J1+2))                                     
C  ***  DEPOLARIZATION FACTOR                                           
        A2 = 0.250*FLOAT(J2*(J2+2))                                     
        BETA = (A2-FLOAT(L2*(L2+1))+0.750)/(A1-FLOAT(L1*(L1+1))+0.750)  
     1  *(A1+A2-FLOAT(L*(L+1)))/(2.000*A2)                              
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
C       COMPLETE
     0  FUNCTION POINT(N,J)                                             
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        DOUBLE PRECISION DN,DZA,DZA2,D,DJ,DREDM                         
C  ***  POINT-LIKE DIRAC ENERGY FUNCTION                                
        COMMON/LOC034/DZA,DZA2,DREDM                                    
        DN = DBLE(FLOAT(N))                                             
        DJ = DBLE(FLOAT(J))                                             
C
C       D = Z/(N-J + sqrt(J^2 - (Z*alpha)^2))
C       Sommerfeld fine structure energy expression! 
        D = DZA/(DN-DJ+DSQRT(DJ*DJ-DZA2))                               
        D = DREDM/DSQRT(1.000D00+D*D)                                   

C       Now convert back to single precision
        POINT = SNGL(DREDM-D)                                           
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION ID(A)                                                  
        INTEGER A,B                                                     
        DIMENSION B(10)                                                 
C  ***  USED IN DECYPHERING TO CONVERT CHARACTERS INTO NUMBERS          
        DATA B/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/                 
        ID = 0                                                          
        DO 100 J=1,10                                                   
        IF(A.EQ.B(J))  ID=J-1                                           
  100   CONTINUE                                                        
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  FUNCTION RAV(N,L,M)                                             
        GO TO (100,200,300),M                                           
C  ***  USED BY CHECK TO FIND THE EXACT EXPECTATION VALUES              
  100   RAV = 0.500*FLOAT(3*N*N-L*(L+1))                                
        RETURN                                                          
  200   RAV = 0.500*FLOAT(N*N*(5*N*N+1-3*L*(L+1)))                      
        RETURN                                                          
  300   RAV = 0.125*FLOAT(N*N*(35*N*N*(N*N-1)-30*N*N*(L+2)*(L-1)+       
     1  3*(L+2)*(L+1)*L*(L-1)))                                         
        RETURN                                                          
        END                                                             
C-----------------------------------------------------------------------
     0  BLOCK DATA                                                      
        DOUBLE PRECISION F                                              
C  ***  BLOCK DATA WITH ALL INTRINSIC PARAMETERS AND DEFAULT VALUES     
        COMMON/LOC001/IJK,ENERGY,ECONS,ECONST,D2P1SM,D2P1S              
        COMMON/LOC002/BEM(3),ZSA(3),BE(3)                               
        COMMON/LOC003/K0,K1,K2,K3                                       
        COMMON/LOC004/NN0(3),NN1(7),NN2(7),NN3(7)                       
        COMMON/LOC006/IP1(7),IP2(7),IP3(7),IQ1(7),IQ2(7),IQ3(7)         
        COMMON/LOC007/IC                                                
        COMMON/LOC008/IREAD,IW,IPUNCH,IPRINT                            
        COMMON/LOC009/F(60),FD                                          
        COMMON/LOC010/PI,PICOEF,AMASSM,NE(3),JK(3),COEFF,ALFA,AMASSE    
        COMMON/LOC011/Z,ZSK,ZSL,ZSM,ZSKZ,ZSLZ,ZSMZ                      
        COMMON/LOC013/COEMON(30),EXPMON(30)                             
        COMMON/LOC014/COEDIP(42),EXPDIP(42)                             
        COMMON/LOC015/COEQUA(45),EXPQUA(45)                             
        COMMON/LOC016/COEOCT(45),EXPOCT(45)                             
        COMMON/LOC017/IFM(6)                                            
        COMMON/LOC018/IFD(9)                                            
        COMMON/LOC019/IFQ(10)                                           
        COMMON/LOC020/IFO(10)                                           
        COMMON/LOC023/COEDP(9)                                          
        COMMON/LOC024/COEQ(11)                                          
        COMMON/LOC025/COEO(11)                                          
        COMMON/LOC026/COED(4)                                           
        COMMON/LOC027/LL(20)                                            
        COMMON/LOC028/M1(7),M2(7),M3(7),YC(4),IDB                       
        COMMON/LOC030/POP(6),JTM(6),JTD(6),JTQ(6),JTO(6)                
        COMMON/LOC031/JM(10),JD(14),JQ(15),JO(15),IYC,IJ(4),YJ(4),JJ1(4)
        COMMON/LOC032/EHIGH,ELOW,CLIMIT,ERES,ESP,ESPM                   
        COMMON/LOC035/ICC,CD(5),EA,EB,IDIR                              
        COMMON/LOC036/ZMK,ZML,ZMM,ZMKM,ZMLM,ZMMM,IVERS                  
        COMMON/LOC037/PL(20),NPOL(20),IPOL,CL1,CL2,IDE,PLN(210),IP8     
        COMMON/LOC038/A,CFM,TFM,STEP,RMATCH,WIDTHK,IPC(3)               
        COMMON/LOC039/NOPT,NMAX,ALEXP                                   
        COMMON/LOC040/AMASSA,AMASSN,HBAR                                
        COMMON/LOC041/MPU,ICPU(200),IPN                                 
        DATA IJK,ENERGY,ECONS,D2P1S/0,1.0E4,5.505355E-03,0.000E+00/     
        DATA BE/3*0.000/                                                
        DATA K0,K1,K2,K3/3, 4,4,4/                                      
        DATA NN0/1,2,3/                                                 
        DATA NN1/0,1,2,3,3,3,3/                                         
        DATA M1/0,0,0,0,1,2,3/                                          
        DATA NN2/0,1,2,3,3,3,3/                                         
        DATA M2/0,0,0,0,1,2,3/                                          
        DATA NN3/0,1,2,3,3,3,3/                                         
        DATA M3/0,0,0,0,1,2,3/                                          
        DATA IP1/0,1,1,1,1,1,1/                                         
        DATA IP2/0,1,1,1,1,1,1/                                         
        DATA IP3/0,1,1,1,1,1,1/                                         
        DATA IQ1/0,0,0,0,0,0,0/                                         
        DATA IQ2/0,0,0,0,0,0,0/                                         
        DATA IQ3/0,0,0,0,0,0,0/                                         
        DATA IDB,YC/0,1.000,1.000,1.000,1.000/                          
        DATA PLN/210*0.000/                                             
        DATA IC,IREAD,IW,IPUNCH,FD,IPRINT,IPOL/3,5,6,7,15.000,0,0/      
        DATA PI,PICOEF,AMASSM/3.1415926535,1.298778E17,206.7686/        
C       ALFA is the fine structure constant
C       COEFF IS FINE STRUCTURE X SPEED OF LIGHT / BOHR RADIUS
        DATA COEFF,ALFA,AMASSE/4.134139E16,7.297353E-03,511003.4/       
        DATA NE,JK/0,2,2,0,2,3/                                         
        DATA Z,ZSK,ZSL,ZSM/4*0.000/                                     
        DATA IFM/6*0/                                                   
        DATA IFD/9*0/                                                   
        DATA IFQ/10*0/                                                  
        DATA IFO/10*0/                                                  
        DATA CL1,CL2/0.00,0.000/                                        
        DATA JTM,JTD,JTQ,JTO/6*1,6*1,6*1,6*1/                           
        DATA COEDP/2.133333E+01,4.266667E+01,3.555556E+00,1.137778E+02, 
     1  7.111111E+00,5.688889E+01,1.02400E+03,2.275556E+01,1.228800E+03/
        DATA COEQ/6.666667E-02,8.888889E-02,6.944444E-04,1.000000E-02,  
     1  9.37500E-04,4.064421E-05,3.511660E-03,6.503074E-05,2.107000E-02,
     2  9.290105E-05,2.064468E-06/                                      
C   
C       These are coefficients that we need for calculating the functions
C       in Akylas Table A.1 for octopole
        DATA COEO/1.693122E-03,1.015873E-02,1.984127E-05,2.125850E-04,  
     1  3.985969E-07,5.734633E-08,1.47462E-05,5.461556E-09,2.654316E-05,
     2  7.646177E-07,1.365389E-06/                                      
        DATA COED/1.333333E+00,2.133333E+01,1.066667E+01,7.111111E+00/  
        DATA IPC/0,0,0/                                                 
        DATA IPN/0/                                                     
        DATA IDIR/0/                                                    
        DATA IVERS/0/                                                   
        DATA LL/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/      
        DATA POP/ 6*1.000/                                              
        DATA JM,JD,JQ,JO/10*1,14*1,15*1,15*1/                           
        DATA EHIGH,ELOW,CLIMIT,ERES/20.000,0.040,1.000E-06,0.000300/    
        DATA ICC,CD,EA,EB/0,1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-5,99.,99./ 
        DATA NPOL/20*-1/                                                
        DATA ESP/0.000/                                                 
        DATA IDE/10000/                                                 
        DATA A,CFM,TFM,STEP,RMATCH,WIDTHK/140.0,0.000,2.3001,0.000E00,  
     1  0.000E00,0.000E00/                                              
        DATA NOPT,NMAX,ALEXP/0,15,0.000/                                
        DATA AMASSA,AMASSN,HBAR/0.000,931.48,6.582173E-16/              
        DATA IP8/0/                                                     
        DATA ZMK,ZML,ZMM,ZMKM,ZMLM,ZMMM/2.0,4.0,9.0,4.0,8.0,18.0/       
        DATA IYC,IJ,YJ,JJ1/0,1,1,1,1,0.000,0.000,0.000,0.000,1,1,1,1/   
        DATA EXPMON/      4.531799,4.706453,4.736786,3.471234,3.150354, 
     1  3.052991,4.703517,3.459110,3.249048,3.150354,2.977822,2.917974, 
     2  3.312419,2.695684,2.533755,4.914818,2.970119,2.700216,6.548446, 
     3  3.199174,2.830713,2.695684,2.415030,2.323143,2.974074,2.526711, 
     4  2.402950,2.415030,2.249824,2.188693/                            
        DATA EXPDIP/      3.766185,4.397381,4.511374,2.429543,2.723240, 
     1  2.747787,2.976262,2.927947,2.886848,4.076930,3.354237,3.192531, 
     2  2.723240,2.747093,2.736171,2.064114,2.209579,2.187568,2.546045, 
     3  2.374398,2.298389,2.876576,2.498902,2.385350,3.993301,2.876576, 
     4  2.652421,5.778643,3.134388,2.797892,2.209579,2.154546,2.118914, 
     5  2.374398,2.234588,2.178575,2.876576,2.498902,2.385350,2.154546, 
     6  2.084353,2.051376/                                              
        DATA EXPQUA/      3.531162,4.240205,4.388131,2.543817,2.683365, 
     1  2.694792,2.131239,2.520683,2.579500,2.566813,2.816281,2.819550, 
     2  2.520683,2.607425,2.618563,1.726067,1.977932,1.998559,2.061514, 
     3  2.104925,2.087356,2.283507,2.199915,2.155173,2.183688,2.283507, 
     4  2.245050,2.663082,2.438142,2.347529,1.977932,1.998217,1.987446, 
     5  2.104925,2.062287,2.036313,4.340990,2.986053,2.726274,2.283507, 
     6  2.199915,2.155173,1.998217,1.972396,1.956081/                   
        DATA EXPOCT/      3.429305,4.153371,4.322206,1.982441,2.395708, 
     1  2.470619,2.341001,2.538774,2.572215,2.195239,2.571516,2.623072, 
     2  2.395708,2.512490,2.537093,1.558161,1.837921,1.876801,1.841217, 
     3  1.946305,1.953014,2.017891,2.026300,2.011369,1.778407,2.017891, 
     4  2.032372,2.117493,2.140989,2.117166,1.837921,1.892340,1.895838, 
     5  1.946305,1.948185,1.938261,2.252893,2.331797,2.283271,2.017891, 
     6  2.026300,2.011369,1.892340,1.921936,1.898511/                   
        DATA COEMON/   9.278462E-01,-4.686876E-02, 2.606690E-03,        
     1   3.258145E-01,-1.650726E-02, 9.192955E-04,-8.728775E-02,        
     2   5.557533E-03,-3.299306E-04, 1.650726E-02,-7.889860E-04,        
     3   1.644805E-05, 1.837462E-01,-9.081039E-03, 5.039261E-04,        
     4  -7.674657E-02, 4.090333E-03,-2.414850E-04, 6.976183E-03,        
     5  -3.308179E-04, 2.022719E-05, 9.886183E-03,-4.700158E-04,        
     6   9.786634E-06,-1.115213E-03, 5.895603E-05,-1.271168E-06,        
     7   4.522730E-05,-1.882362E-06, 2.089122E-08/                      
        DATA COEDIP/   9.864278E-02,-3.552196E-03, 6.908681E-05,        
     1   2.415974E-02,-8.812569E-04, 1.721286E-05,-6.819236E-03,        
     2   3.101739E-04,-6.390954E-06, 7.532858E-02,-4.445121E-03,        
     3   2.592835E-04, 2.203142E-04,-8.210685E-06, 8.801341E-08,        
     4   3.235460E-02,-1.181545E-03, 2.303800E-05,-1.225112E-02,        
     5   5.542601E-04,-1.140757E-05, 8.892663E-04,-4.580119E-05,        
     6   9.800079E-07, 1.211143E-02,-6.669497E-04, 3.872454E-05,        
     7  -1.512719E-03, 7.933719E-05,-4.808937E-06, 4.376093E-05,        
     8  -1.627056E-06, 1.743224E-08,-5.132038E-06, 2.095740E-07,        
     9  -2.313305E-09, 2.223166E-04,-1.145030E-05, 2.450020E-07,        
     A   9.039198E-08,-3.079553E-09, 2.088144E-11/                      
        DATA COEQUA/   1.411675E-01,-3.944496E-03, 3.945040E-05,        
     1   3.270497E-01,-1.133197E-02, 1.194093E-04,-2.776085E-01,        
     2   7.846003E-03,-7.866167E-05, 4.601989E-01,-1.831626E-02,        
     3   3.675330E-04, 3.487113E-03,-1.060614E-04, 6.944398E-07,        
     4   1.253120E-00,-3.541941E-02, 3.549579E-04,-4.926037E-01,        
     5   1.705249E-02,-1.796112E-04, 3.681114E-02,-1.438755E-03,        
     6   1.569441E-05, 4.626170E-01,-1.840557E-02, 3.690030E-04,        
     7  -4.608811E-02, 2.220260E-03,-4.671028E-05, 7.870980E-03,        
     8  -2.392062E-04, 1.565542E-06,-9.473607E-04, 3.134314E-05,        
     9  -2.107953E-07, 4.492433E-02,-2.608238E-03, 1.552426E-04,        
     A   9.202786E-03,-3.596887E-04, 3.923602E-06, 5.980154E-05,        
     B  -1.722622E-06, 7.889632E-09/                                    
        DATA COEOCT/   1.835095E-02,-4.182848E-04, 2.558331E-06,        
     1   1.444508E-01,-3.328672E-03, 2.037957E-05,-4.367688E-02,        
     2   1.223204E-03,-7.847840E-06, 3.021241E-01,-8.987655E-03,        
     3   9.179069E-05, 4.438300E-03,-1.140520E-04, 5.042711E-07,        
     4   1.467168E-00,-3.379181E-02, 2.068045E-04,-5.924030E-01,        
     5   1.655453E-02,-1.061699E-04, 4.508567E-02,-1.416945E-03,        
     6   9.394653E-06, 6.820167E-01,-2.028855E-02, 2.071180E-04,        
     7  -6.910950E-02, 2.487891E-03,-2.660315E-05, 2.252870E-02,        
     8  -5.785850E-04, 2.556806E-06,-2.759089E-03, 7.677244E-05,        
     9  -3.479514E-07, 3.778402E-01,-1.582449E-02, 3.229544E-04,        
     A   4.508567E-03,-1.416945E-04, 9.394653E-07, 1.285744E-06,        
     B  -3.257990E-08, 1.067759E-10/                                    
        END                                                             
