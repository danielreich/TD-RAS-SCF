c
c     below are taken from lib
c
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       
C                                                                      
      INTEGER N,NM,IERR,MATZ                                          
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)            
C                                                                   
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 RS 00060
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     RS 00070
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             RS 00080
C     OF A REAL SYMMETRIC MATRIX.                                       RS 00090
C                                                                       RS 00100
C     ON INPUT                                                          RS 00110
C                                                                       RS 00120
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    RS 00130
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            RS 00140
C        DIMENSION STATEMENT.                                           RS 00150
C                                                                       RS 00160
C        N  IS THE ORDER OF THE MATRIX  A.                              RS 00170
C                                                                       RS 00180
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.                         RS 00190
C                                                                       RS 00200
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              RS 00210
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO          RS 00220
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    RS 00230
C                                                                       RS 00240
C     ON OUTPUT                                                         RS 00250
C                                                                       RS 00260
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.                RS 00270
C                                                                       RS 00280
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.              RS 00290
C                                                                       RS 00300
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR      RS 00310
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT   RS 00320
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.              RS 00330
C                                                                       RS 00340
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   RS 00350
C                                                                       RS 00360
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    RS 00370
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY RS 00380
C                                                                       RS 00390
C     THIS VERSION DATED AUGUST 1983.                                   RS 00400
C                                                                       RS 00410
C     ------------------------------------------------------------------RS 00420
C                                                                       RS 00430
      IF (N .LE. NM) GO TO 10                                          
      IERR = 10 * N                                                     
      GO TO 50                                                          
C                                                                     
   10 IF (MATZ .NE. 0) GO TO 20                                      
C     .......... FIND EIGENVALUES ONLY ..........                       RS 00490
      CALL  TRED1(NM,N,A,W,FV1,FV2)                                    
      CALL  TQLRAT(N,W,FV2,IERR)                                       
      GO TO 50                                                        
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........      RS 00530
   20 CALL  TRED2(NM,N,A,W,FV1,Z)                                     
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                  
   50 RETURN                                                        
      END                                                          
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                 
C                                                                     
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                       
      DOUBLE PRECISION D(N),E(N),Z(NM,N)                            
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG  
C                                                                      
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     TQL00070
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     TQL00080
C     WILKINSON.                                                        TQL00090
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   TQL00100
C                                                                       TQL00110
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            TQL00120
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               TQL00130
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              TQL00140
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  TQL00150
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  TQL00160
C                                                                       TQL00170
C     ON INPUT                                                          TQL00180
C                                                                       TQL00190
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TQL00200
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TQL00210
C          DIMENSION STATEMENT.                                         TQL00220
C                                                                       TQL00230
C        N IS THE ORDER OF THE MATRIX.                                  TQL00240
C                                                                       TQL00250
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00260
C                                                                       TQL00270
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        TQL00280
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               TQL00290
C                                                                       TQL00300
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           TQL00310
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      TQL00320
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        TQL00330
C          THE IDENTITY MATRIX.                                         TQL00340
C                                                                       TQL00350
C      ON OUTPUT                                                        TQL00360
C                                                                       TQL00370
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00380
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          TQL00390
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        TQL00400
C                                                                       TQL00410
C        E HAS BEEN DESTROYED.                                          TQL00420
C                                                                       TQL00430
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           TQL00440
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     TQL00450
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       TQL00460
C          EIGENVALUES.                                                 TQL00470
C                                                                       TQL00480
C        IERR IS SET TO                                                 TQL00490
C          ZERO       FOR NORMAL RETURN,                                TQL00500
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00510
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00520
C                                                                       TQL00530
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00540
C                                                                       TQL00550
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00560
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00570
C                                                                       TQL00580
C     THIS VERSION DATED AUGUST 1983.                                   TQL00590
C                                                                       TQL00600
C     ------------------------------------------------------------------TQL00610
C                                                                       TQL00620
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                      
      DO 100 I = 2, N                                                 
  100 E(I-1) = E(I)                                                    
C                                                                       
      F = 0.0D0                                                         
      TST1 = 0.0D0                                                     
      E(N) = 0.0D0                                                      
C                                                                      
      DO 240 L = 1, N                                                 
         J = 0                                                          
         H = DABS(D(L)) + DABS(E(L))                                   
         IF (TST1 .LT. H) TST1 = H                                      
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........        
         DO 110 M = L, N                                                
            TST2 = TST1 + DABS(E(M))                                   
            IF (TST2 .EQ. TST1) GO TO 120                             
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C                THROUGH THE BOTTOM OF THE LOOP ..........             
  110    CONTINUE                                                     
C                                                                       TQL00840
  120    IF (M .EQ. L) GO TO 220                                      
  130    IF (J .EQ. 30) GO TO 1000                                     
         J = J + 1                                                      
C     .......... FORM SHIFT ..........                                  TQL00880
         L1 = L + 1                                                    
         L2 = L1 + 1                                                  
         G = D(L)                                                       
         P = (D(L1) - G) / (2.0D0 * E(L))                              
         R = PYTHAG(P,1.0D0)                                           
         D(L) = E(L) / (P + DSIGN(R,P))                              
         D(L1) = E(L) * (P + DSIGN(R,P))                               
         DL1 = D(L1)                                                  
         H = G - D(L)                                                   
         IF (L2 .GT. N) GO TO 145                                      
C                                                                     
         DO 140 I = L2, N                                          
  140    D(I) = D(I) - H                                             
C                                                                   
  145    F = F + H                                                 
C     .......... QL TRANSFORMATION ..........                     
         P = D(M)                                                
         C = 1.0D0                                                     
         C2 = C                                                       
         EL1 = E(L1)                                                 
         S = 0.0D0                                                  
         MML = M - L                                                    
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........            
         DO 200 II = 1, MML                                           
            C3 = C2                                                  
            C2 = C                                                  
            S2 = S                                                     
            I = M - II                                                
            G = C * E(I)                                             
            H = C * P                                                  
            R = PYTHAG(P,E(I))                                         
            E(I+1) = S * R                                            
            S = E(I) / R                                             
            C = P / R                                              
            P = C * D(I) - S * G                                        
            D(I+1) = H + S * (C * G + S * D(I))                         
C     .......... FORM VECTOR ..........                                
            DO 180 K = 1, N                                            
               H = Z(K,I+1)                                            
               Z(K,I+1) = S * Z(K,I) + C * H                           
               Z(K,I) = C * Z(K,I) - S * H                             
  180       CONTINUE                                                    
C                                                                       TQL01310
  200    CONTINUE                                                     
C                                                                      
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            
         E(L) = S * P                                                   
         D(L) = C * P                                                  
         TST2 = TST1 + DABS(E(L))                                       
         IF (TST2 .GT. TST1) GO TO 130                                 
  220    D(L) = D(L) + F                                               
  240 CONTINUE                                                          
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........         
      DO 300 II = 2, N                                                 
         I = II - 1                                                     
         K = I                                                         
         P = D(I)                                                       
C                                                                      
         DO 260 J = II, N                                               
            IF (D(J) .GE. P) GO TO 260                                 
            K = J                                                       
            P = D(J)                                                    
  260    CONTINUE                                                      
C                                                                    
         IF (K .EQ. I) GO TO 300                                       
         D(K) = D(I)                                                    
         D(I) = P                                                      
C                                                                       
         DO 280 J = 1, N                                              
            P = Z(J,I)                                                 
            Z(J,I) = Z(J,K)                                             
            Z(J,K) = P                                                  
  280    CONTINUE                                                       
C                                                                      
  300 CONTINUE                                                         
C                                                                       
      GO TO 1001                                                       
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ..........             
 1000 IERR = L                                                          
 1001 RETURN                                                            
      END                                                              
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                  
C                                                                       TQL00020
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  
      DOUBLE PRECISION D(N),E2(N)                                       
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG                 
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,     TQL00070
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                TQL00080
C                                                                       TQL00090
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              TQL00100
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     TQL00110
C                                                                       TQL00120
C     ON INPUT                                                          TQL00130
C                                                                       TQL00140
C        N IS THE ORDER OF THE MATRIX.                                  TQL00150
C                                                                       TQL00160
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00170
C                                                                       TQL00180
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     TQL00190
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. TQL00200
C                                                                       TQL00210
C      ON OUTPUT                                                        TQL00220
C                                                                       TQL00230
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00240
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          TQL00250
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            TQL00260
C          THE SMALLEST EIGENVALUES.                                    TQL00270
C                                                                       TQL00280
C        E2 HAS BEEN DESTROYED.                                         TQL00290
C                                                                       TQL00300
C        IERR IS SET TO                                                 TQL00310
C          ZERO       FOR NORMAL RETURN,                                TQL00320
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00330
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00340
C                                                                       TQL00350
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00360
C                                                                       TQL00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00390
C                                                                       TQL00400
C     THIS VERSION DATED AUGUST 1983.                                   TQL00410
C                                                                       TQL00420
C     ------------------------------------------------------------------TQL00430
C                                                                      
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                      
      DO 100 I = 2, N                                                  
  100 E2(I-1) = E2(I)                                                   
C                                                                      
      F = 0.0D0                                                       
      T = 0.0D0                                                        
      E2(N) = 0.0D0                                                    
C                                                                       
      DO 290 L = 1, N                                                   
         J = 0                                                          
         H = DABS(D(L)) + DSQRT(E2(L))                                  
         IF (T .GT. H) GO TO 105                                        
         T = H                                                          
         B = EPSLON(T)                                                 
         C = B * B                                                      
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT .......... 
  105    DO 110 M = L, N                                                
            IF (E2(M) .LE. C) GO TO 120                                 
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              TQL00650
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL00660
  110    CONTINUE                                                      
C                                                                       TQL00680
  120    IF (M .EQ. L) GO TO 210                                        
  130    IF (J .EQ. 30) GO TO 1000                                      
         J = J + 1                                                     
C     .......... FORM SHIFT ..........                                  TQL00720
         L1 = L + 1                                                    
         S = DSQRT(E2(L))                                               
         G = D(L)                                                      
         P = (D(L1) - G) / (2.0D0 * S)                                
         R = PYTHAG(P,1.0D0)                                            
         D(L) = S / (P + DSIGN(R,P))                                   
         H = G - D(L)                                                 
C                                                                       TQL00800
         DO 140 I = L1, N                                               
  140    D(I) = D(I) - H                                               
C                                                                     
         F = F + H                                                      
C     .......... RATIONAL QL TRANSFORMATION ..........                 TQL00850
         G = D(M)                                                      
         IF (G .EQ. 0.0D0) G = B                                        
         H = G                                                         
         S = 0.0D0                                                      
         MML = M - L                                                   
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL00910
         DO 200 II = 1, MML                                             
            I = M - II                                                 
            P = G * H                                                 
            R = P + E2(I)                                             
            E2(I+1) = S * R                                           
            S = E2(I) / R                                              
            D(I+1) = H + S * (H + D(I))                              
            G = D(I) - E2(I) / G                                        
            IF (G .EQ. 0.0D0) G = B                                    
            H = G * P / R                                               
  200    CONTINUE                                                      
C                                                                     
         E2(L) = S * G                                                 
         D(L) = H                                                       
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST .......... TQL01060
         IF (H .EQ. 0.0D0) GO TO 210                                
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                    
         E2(L) = H * E2(L)                                              
         IF (E2(L) .NE. 0.0D0) GO TO 130                              
  210    P = D(L) + F                                                   
C     .......... ORDER EIGENVALUES ..........                          
         IF (L .EQ. 1) GO TO 250                                       
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........              
         DO 230 II = 2, L                                              
            I = L + 2 - II                                             
            IF (P .GE. D(I-1)) GO TO 270                                
            D(I) = D(I-1)                                             
  230    CONTINUE                                                       
C                                                                      
  250    I = 1                                                         
  270    D(I) = P                                                       
  290 CONTINUE                                                        
C                                                                     
      GO TO 1001                                                      
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL01260
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL01270
 1000 IERR = L                                                          
 1001 RETURN                                                            
      END                                                              
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                 
C                                                                     
      INTEGER I,J,K,L,N,II,NM,JP1                                     
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)                          
      DOUBLE PRECISION F,G,H,SCALE                                      
C                                                                      
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    TRE00070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRE00080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRE00090
C                                                                       TRE00100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   TRE00110
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           TRE00120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRE00130
C                                                                       TRE00140
C     ON INPUT                                                          TRE00150
C                                                                       TRE00160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRE00170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             TRE00280
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         TRE00290
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.        TRE00300
C                                                                       TRE00310
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00320
C                                                                       TRE00330
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00340
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00350
C                                                                       TRE00360
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    TRE00370
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        TRE00380
C                                                                       TRE00390
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00400
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00410
C                                                                       TRE00420
C     THIS VERSION DATED AUGUST 1983.                                   TRE00430
C                                                                       TRE00440
C     ------------------------------------------------------------------TRE00450
C                                                                       TRE00460
      DO 100 I = 1, N                                                  
         D(I) = A(N,I)                                                  
         A(N,I) = A(I,I)                                                
  100 CONTINUE                                                         
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........              
      DO 300 II = 1, N                                                  
         I = N + 1 - II                                                 
         L = I - 1                                                      
         H = 0.0D0                                                      
         SCALE = 0.0D0                                                  
         IF (L .LT. 1) GO TO 130                                        
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRE00580
         DO 120 K = 1, L                                               
  120    SCALE = SCALE + DABS(D(K))                                     
C                                                                      
         IF (SCALE .NE. 0.0D0) GO TO 140                                
C                                                                       
         DO 125 J = 1, L                                          
            D(J) = A(L,J)                                               
            A(L,J) = A(I,J)                                             
            A(I,J) = 0.0D0                                              
  125    CONTINUE                                                       
C                                                                       
  130    E(I) = 0.0D0                                                   
         E2(I) = 0.0D0                                                  
         GO TO 300                                                      
C                                                                       
  140    DO 150 K = 1, L                                               
            D(K) = D(K) / SCALE                                         
            H = H + D(K) * D(K)                                         
  150    CONTINUE                                                       
C                                                                      
         E2(I) = SCALE * SCALE * H                                      
         F = D(L)                                                      
         G = -DSIGN(DSQRT(H),F)                                         
         E(I) = SCALE * G                                              
         H = H - F * G                                                  
         D(L) = F - G                                                  
         IF (L .EQ. 1) GO TO 285                                      
C     .......... FORM A*U ..........                                    
         DO 170 J = 1, L                                              
  170    E(J) = 0.0D0                                                 
C                                                                       
         DO 240 J = 1, L                                               
            F = D(J)                                                    
            G = E(J) + A(J,J) * F                                    
            JP1 = J + 1                                               
            IF (L .LT. JP1) GO TO 220                                   
C                                                                      
            DO 200 K = JP1, L                                           
               G = G + A(K,J) * D(K)                                  
               E(K) = E(K) + A(K,J) * F                                
  200       CONTINUE                                                  
C                                                                      
  220       E(J) = G                                               
  240    CONTINUE                                                       
C     .......... FORM P ..........                                    
         F = 0.0D0                                                      
C                                                                     
         DO 245 J = 1, L                                               
            E(J) = E(J) / H                                        
            F = F + E(J) * D(J)                                         
  245    CONTINUE                                                
C                                                                      
         H = F / (H + H)                                                
C     .......... FORM Q ..........                                      
         DO 250 J = 1, L                                              
  250    E(J) = E(J) - H * D(J)                                        
C     .......... FORM REDUCED A ..........                             
         DO 280 J = 1, L                                               
            F = D(J)                                                   
            G = E(J)                                                    
C                                                                       
            DO 260 K = J, L                                            
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)                     
C                                                                       
  280    CONTINUE                                                      
C                                                                       
  285    DO 290 J = 1, L                                                
            F = D(J)                                                   
            D(J) = A(L,J)                                              
            A(L,J) = A(I,J)                                            
            A(I,J) = F * SCALE                                          
  290    CONTINUE                                                       
C                                                                       
  300 CONTINUE                                                         
C                                                                       
      RETURN                                                           
      END                                                               
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    
C                                                                       
      INTEGER I,J,K,L,N,II,NM,JP1                                       
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)                        
      DOUBLE PRECISION F,G,H,HH,SCALE                                  
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,  
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A             
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                           
C                                                                   
C     ON INPUT                                                    
C                                                                      
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL       
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00280
C                                                                       TRE00290
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00300
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00310
C                                                                       TRE00320
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                TRE00330
C          PRODUCED IN THE REDUCTION.                                   TRE00340
C                                                                       TRE00350
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            TRE00360
C                                                                       TRE00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00390
C                                                                       TRE00400
C     THIS VERSION DATED AUGUST 1983.                                   TRE00410
C                                                                       TRE00420
C     ------------------------------------------------------------------TRE00430
C                                                                       TRE00440
      DO 100 I = 1, N                                                   
C                                                                       TRE00460
         DO 80 J = I, N                                                 
   80    Z(J,I) = A(J,I)                                               
C                                                                       
         D(I) = A(N,I)                                                  
  100 CONTINUE                                                        
C                                                                    
      IF (N .EQ. 1) GO TO 510                                         
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               
      DO 300 II = 2, N                                               
         I = N + 2 - II                                               
         L = I - 1                                                    
         H = 0.0D0                                                    
         SCALE = 0.0D0                                                 
         IF (L .LT. 2) GO TO 130                                      
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........      
         DO 120 K = 1, L                                                
  120    SCALE = SCALE + DABS(D(K))                                   
C                                                                       
         IF (SCALE .NE. 0.0D0) GO TO 140                                
  130    E(I) = D(L)                                                    
C                                                                       
         DO 135 J = 1, L                                              
            D(J) = Z(L,J)                                              
            Z(I,J) = 0.0D0                                            
            Z(J,I) = 0.0D0                                           
  135    CONTINUE                                                       
C                                                                       
         GO TO 290                                                      
C                                                                       
  140    DO 150 K = 1, L                                                
            D(K) = D(K) / SCALE                                         
            H = H + D(K) * D(K)                                         
  150    CONTINUE                                                      
C                                                                       
         F = D(L)                                                     
         G = -DSIGN(DSQRT(H),F)                                      
         E(I) = SCALE * G                                             
         H = H - F * G                                                 
         D(L) = F - G                                                
C     .......... FORM A*U ..........                                 
         DO 170 J = 1, L                                           
  170    E(J) = 0.0D0                                                 
C                                                                      
         DO 240 J = 1, L                                              
            F = D(J)                                                    
            Z(J,I) = F                                                 
            G = E(J) + Z(J,J) * F                                  
            JP1 = J + 1                                                
            IF (L .LT. JP1) GO TO 220                                  
C                                                                    
            DO 200 K = JP1, L                                         
               G = G + Z(K,J) * D(K)                                   
               E(K) = E(K) + Z(K,J) * F                                
  200       CONTINUE                                                    
C                                                                    
  220       E(J) = G                                                  
  240    CONTINUE                                                    
C     .......... FORM P ..........                                   
         F = 0.0D0                                                   
C                                                                    
         DO 245 J = 1, L                                           
            E(J) = E(J) / H                                             
            F = F + E(J) * D(J)                                        
  245    CONTINUE                                                       
C                                                                      
         HH = F / (H + H)                                             
C     .......... FORM Q ..........                                      
         DO 250 J = 1, L                                              
  250    E(J) = E(J) - HH * D(J)                                        
C     .......... FORM REDUCED A ..........                              
         DO 280 J = 1, L                                                
            F = D(J)                                                    
            G = E(J)                                                  
C                                                                      
            DO 260 K = J, L                                          
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)                      
C                                                                     
            D(J) = Z(L,J)                                             
            Z(I,J) = 0.0D0                                         
  280    CONTINUE                                                      
C                                                                     
  290    D(I) = H                                                  
  300 CONTINUE                                                          
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     TRE01300
      DO 500 I = 2, N                                                  
         L = I - 1                                                     
         Z(N,L) = Z(L,L)                                                
         Z(L,L) = 1.0D0                                               
         H = D(I)                                                      
         IF (H .EQ. 0.0D0) GO TO 380                                    
C                                                                       
         DO 330 K = 1, L                                               
  330    D(K) = Z(K,I) / H                                              
C                                                                       
         DO 360 J = 1, L                                                
            G = 0.0D0                                                   
C                                                                       
            DO 340 K = 1, L                                             
  340       G = G + Z(K,I) * Z(K,J)                                    
C                                                                       
            DO 360 K = 1, L                                            
               Z(K,J) = Z(K,J) - G * D(K)                               
  360    CONTINUE                                                      
C                                                                      
  380    DO 400 K = 1, L                                               
  400    Z(K,I) = 0.0D0                                                 
C                                                                     
  500 CONTINUE                                                         
C                                                                       
  510 DO 520 I = 1, N                                                  
         D(I) = Z(N,I)                                                  
         Z(N,I) = 0.0D0                                               
  520 CONTINUE                                                          
C                                                                      
      Z(N,N) = 1.0D0                                                
      E(1) = 0.0D0                                                      
      RETURN                                                       
      END                                                             







      DOUBLE PRECISION FUNCTION EPSLON (X)                            
      DOUBLE PRECISION X                                               
C                                                                    
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.                  
C                                                                     
      DOUBLE PRECISION A,B,C,EPS                                       
C                                                                     
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS            
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,                      
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT               
C            NUMBERS IS NOT A POWER OF THREE.                           EPS00110
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO         EPS00120
C            THE ACCURACY USED IN FLOATING POINT VARIABLES              EPS00130
C            THAT ARE STORED IN MEMORY.                                 EPS00140
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO          EPS00150
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING            EPS00160
C     ASSUMPTION 2.                                                     EPS00170
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,                  EPS00180
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,                    EPS00190
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,                   EPS00200
C            C  IS NOT EXACTLY EQUAL TO ONE,                            EPS00210
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM                   EPS00220
C                 THE NEXT LARGER FLOATING POINT NUMBER.                EPS00230
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED         EPS00240
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.            EPS00250
C                                                                       EPS00260
C     THIS VERSION DATED 4/6/83.                                        EPS00270
C                                                                       EPS00280
      A = 4.0D0/3.0D0                                                  
   10 B = A - 1.0D0                                                    
      C = B + B + B                                                     
      EPS = DABS(C-1.0D0)                                              
      IF (EPS .EQ. 0.0D0) GO TO 10                                     
      EPSLON = EPS*DABS(X)                                             
      RETURN                                                          
      END                                                            
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)                         
      DOUBLE PRECISION A,B                                              
C                                                                       
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW  
C                                                                       
      DOUBLE PRECISION P,R,S,T,U                                        
      P = DMAX1(DABS(A),DABS(B))                                        
      IF (P .EQ. 0.0D0) GO TO 20                                        
      R = (DMIN1(DABS(A),DABS(B))/P)**2                                 
   10 CONTINUE                                                          
         T = 4.0D0 + R                                                  
         IF (T .EQ. 4.0D0) GO TO 20                                     
         S = R/T                                                        
         U = 1.0D0 + 2.0D0*S                                            
         P = U*P                                                        
         R = (S/U)**2 * R                                              
      GO TO 10                                                        
   20 PYTHAG = P                                                    
      RETURN                                                       
      END     
	
	
 
