!
!
!                             G E N C L
!                             =========
!      
!     ------------------------------------------------------------------

      PROGRAM GENCLF
!
!                C O P Y R I G H T -- 1994
!
!                By C. Froese Fischer
!                   Bin Liu
!                   Vanderbilt University
!                   Nashville, TN 37235 USA
!
!                Written: May, 1983
!
!     Computer Physics Communications, Vol. 64, 406--416 (1991).
!
!     Modified   By G. Gaigalas for f-shells  (1997)
!
!     Cleaned up and checked by J. Grumer (2015) by hand and by using 
!        the PlusFORT routine SPAG. This was done to simplify translation
!        to f90
!
!     ------------------------------------------------------------------
!
!
!       This program computes all possible couplings for each member of
!   the reference set, generates all unique, possible configurations
!   and their couplings  from the active set, and for each replacement,
!   generates configurations and their couplings, then outputs the
!   configurations and couplings with given final term(s).
!     Input : (Interactive)
!        i)   Header ;
!       ii)   List of closed shells ;
!      iii)   The "reference" set ;
!       iv)   The "active" set ;
!        v)   Replacements from the reference set ;
!             Virtual Set if Replacement is 's' or 'd' or 'sd'
!       vi)   Final Term
!     Output :
!       i )   Header ;
!       ii)   List of closed shells ;
!      iii)   Configurations (FORMAT(8(1X,A3,'(',I2,')'))
!               and their couplings (FORMAT(15(1X,A3))


!   I/O allocation
!    FILE1(.)  :  Configurations by Active Set  (internal file)
!    FILE2(.)  :  Configurations by Replacement (internal file)
!    FILE3(.)  :  couplings                     (internal file)
!           6  :  Terminal output
!           7  :  File CLIST
!  FBETA(1,.)  :  Information of Beta2          (internal file )
!  FBETA(2,.)  :  Information of Beta3          (internal file )
!  FBETA(3,.)  :  Information of Beta4          (internal file )
!  FBETA(4,.)  :  Information of Beta5          (internal file )
!  FBETA(5,.)  :  Information of Beta6          (internal file )
!  FBETA(6,.)  :  Information of Beta7          (internal file )
!  FBETA(7,.)  :  Information of Beta8          (internal file )
!
!
! ---------------------------------------------------------------------
!               M A I N         P R O G R A M
! ---------------------------------------------------------------------
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS, NSHEL, NCOUPL, NCFG, NSCOUP
      PARAMETER (NELS=21,NSHEL=8,NCOUPL=2*NSHEL-1,NCFG=10000,
     &           NSCOUP=10000)
!
! COMMON variables
!
      CHARACTER*60 ACT, REF(NELS), REPL(NELS)
      INTEGER CONST, M, MA, MAXFINAL, MINFINAL, ML, MR, MS(NELS), NF,
     &        NFTM, NQ, NR, NREF, ORD0, ORD9, ORDLA, ORDLZ, ORDUA,
     &        ORDUZ, PARITYVAL, Q(NELS), Q1, Q10, Q11, Q12, Q13, Q14,
     &        Q15, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, QL(NELS), QR(NELS),
     &        QS(NELS,NELS), RL(NELS)
      CHARACTER*3 EL(NELS), ELA(NELS), ELL(NELS), ELR(NELS), ELS(NELS)
      CHARACTER*8 FBETA(7,NSCOUP)
      CHARACTER*32 FILE1(NCFG), FILE3(NSCOUP)
      CHARACTER*40 FILE2(NCFG)
      CHARACTER*2 FTM(NELS)
      CHARACTER*72 HEADER, SHELLS, VIRTUL
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK0  / ORDLA, ORDLZ, ORDUA, ORDUZ, ORD0, ORD9
      COMMON /BLK1  / HEADER, SHELLS, ACT, VIRTUL
      COMMON /BLK2  / REF, REPL, FTM
      COMMON /BLK3  / EL, ELL, ELR, ELS, ELA
      COMMON /BLK4  / Q, QL, ML, QR, MR, M, QS, MS, MA, RL, NREF, Q1,
     &                Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12,
     &                Q13, Q14, Q15
      COMMON /FILES / FBETA, FILE1, FILE2, FILE3
!
! Local variables
!
      CHARACTER CH1
      CHARACTER*2 CH2
      CHARACTER*3 CH3, COUPLE(NCOUPL), ELB(NELS), ELC(NELS), ELV(NELS)
      INTEGER DFLAG, F(NELS), H(NELS), I, IFLAG, ITYPE, J, K, L, LEFT,
     &        LL1, LL2, LLMAX, LLMIN, LMAX, LMIN, LVIR, MINQ1, MV, N,
     &        N1, N2, NC, NL, NREPL, NVIR, PL, QA(NELS), QB(NELS),
     &        RIGHT, SFLAG, VL(NELS)
      CHARACTER*24 FILENAME
      CHARACTER*60 FINALS(NELS), STRL, STRR
      INTEGER ICTOI, LVAL
      DOUBLE PRECISION S, SUMQ
      CHARACTER*72 TEMP
      EXTERNAL CONFIG, COUPLD, DECOMP, DEL, HELP, ICTOI, INPUT, LVAL,
     &         PRINT_COUPLINGS, REPLAC, STRSH, VPAIR
!
!      
!**************
!     Declaration of variable for the input data
!
!     HEADER  =  header of output
!     SHELLS  =  Closed Shells
!     REF  =  Reference set
!     ACT  =  Active Set
!     REPL  =  Replacement
!     VIRTUL  = Virtual Set
!     FTM  = Final Term (NELS)
!
      FILENAME = 'cfg.inp'

!**************

!     Obtain ASCII value for the bound of character and digit
!
      ORDLA = ICHAR('a')
      ORDLZ = ICHAR('z')
      ORDUA = ICHAR('A')
      ORDUZ = ICHAR('Z')
      ORD0 = ICHAR('0')
      ORD9 = ICHAR('9')
      NQ = 0
      IFLAG = 0
      SFLAG = 0
      DFLAG = 0

!***************
!     Disply for the beginning of the program
!

      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,*)
     &'     -----------------------------------------------------------'
      WRITE (0,*)
      WRITE (0,*) '     You are under the program GENCL'
      WRITE (0,*) '            which GENerates a Configuration List'
      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,*)
     &'     Type H (Help) to obtain information about the input format '
      WRITE (0,*)
      WRITE (0,*) '     Type <RETURN> if you already know '
      WRITE (0,*)
      WRITE (0,*)
     &   '     --------------------------------------------------------'
      WRITE (0,*)

!**************
!     If the user needs to obtain information about input format,
!     call subroutine HELP .
!
      READ (5,'(A)') STRL
      I = INDEX(STRL,'H')
      J = INDEX(STRL,'h')
      IF ( I.NE.0 .OR. J.NE.0 ) CALL HELP

!**********************************************************************
!     READ AND ANALYSIS THE INPUT DATA
!
!
!     Input Header
!
 100  WRITE (0,99001)
99001 FORMAT (T13,'Header  ?  ')
      READ (5,'(A)') HEADER

      CH1 = CH3(2:2)              ! CHECK THIS! IS CH3 DEFINED?

!
!**************
!     Input Closed Shell, check if it satisfies the FORMAT(18(1X,A3))
!
 200  WRITE (0,99002)
99002 FORMAT ('     Closed Shells  ?  ')
      READ (5,'(A)') TEMP
      I = INDEX(TEMP,'B')
      J = INDEX(TEMP,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 100
      N = 2
      SHELLS = ' '
 300  IF ( INDEX(TEMP,'     ').GT.2 ) THEN
         CALL DEL(TEMP)
         J = ICHAR(TEMP(1:1))
         IF ( J.LT.ORD0 .OR. J.GT.ORD9 ) THEN
            WRITE (0,*) '     INPUT ERROR !'
            GO TO 200
         END IF
         SHELLS(N:N+2) = TEMP(:3)
         N = N + 4
         CALL STRSH(TEMP,4)
         GO TO 300
      END IF


!**************
!     Input Reference Set, and check the input error
!
!     LEFT     =  position of the first '(' in the string
!     RIGHT     =  position of the first ')' in the string
!     SUMQ    =  total number of Q(i) for each configuration
!     NL     =  Total quantum number for each reference set
!     CONST     =  Total number of Qi for the first reference set
!     PARITYVAL  =  Parity for the first reference set
!
!
 400  WRITE (0,99003)
99003 FORMAT ('     Reference Set  ?  ')
      CALL INPUT(NELS,NREF,REF,IFLAG,SFLAG,DFLAG,*200,*400)
      DO I = 1, NREF
         SUMQ = 0
         NL = 0
         TEMP = REF(I)
 450     IF ( TEMP(1:5).NE.'     ' ) THEN
            CALL DEL(TEMP)
!
!           Error if the input has unmatched parenthesis
!
            LEFT = INDEX(TEMP,'(')
            RIGHT = INDEX(TEMP(:7),')')
            IF ( LEFT.EQ.0 .OR. RIGHT.EQ.0 ) THEN
               WRITE (0,*) '       Unmatched parenthesis!',
     &                     ' Please input again .'
               GO TO 400
            END IF
!
!           Error if number of electron is more than FULL
!
            CH3 = TEMP(LEFT+1:RIGHT-1)
            N = ICTOI(CH3)
            L = LVAL(TEMP(2:2))
            IF ( N.GT.L*4+2 ) THEN
               WRITE (0,*)
     &        '     Number of electrons in a shell is  more than FULL !'
               GO TO 400
            END IF
            SUMQ = SUMQ + N
            NL = NL + N*L
            CALL STRSH(TEMP,RIGHT+1)
            GO TO 450
         END IF

!**************
!        The first Reference Set is defined as the standard value
!
         IF ( I.EQ.1 ) THEN
            CONST = INT(SUMQ)
            PARITYVAL = MOD(NL,2)
         ELSE
!
!           Error if members of the reference set have different
!           electron numbers or parity
            IF ( INT(SUMQ).NE.CONST ) THEN
               WRITE (0,*) '     Total number of electrons is wrong!'
               GO TO 400
            END IF
            IF ( MOD(NL,2).NE.PARITYVAL ) THEN
               WRITE (0,*) '       Parity is wrong!'
               GO TO 400
            END IF
         END IF
      END DO

!******************
!     Input Active Set
!
 500  WRITE (0,99004)
99004 FORMAT ('        Active Set  ?  ')
      READ (5,'(A)') ACT
      I = INDEX(ACT,'B')
      J = INDEX(ACT,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 400
      IF ( ACT.NE.'   ' ) THEN
         WRITE (0,99005)
99005    FORMAT ('Type of set generation ?')
         READ *, ITYPE
      END IF
!*******************
!     Input Replacement
!
 600  WRITE (0,99006)
99006 FORMAT ('      Replacements  ?  ')
      IFLAG = 1
      CALL INPUT(NELS,NREPL,REPL,IFLAG,SFLAG,DFLAG,*500,*1000)

!*******************
!     If Replacement = S or D or SD, input Virtual Set
!
!     ELV  =  Electron label for virtual set
!     MV  =  Number of electron for Virtual Set
!     VL  =  Parity for each Qi in Virtual Set
!
      IF ( SFLAG.EQ.0 .AND. DFLAG.EQ.0 ) GO TO 1000

      WRITE (0,99007)
99007 FORMAT ('       Virtual Set  ?  ')
      READ (5,'(A)') VIRTUL
      TEMP = VIRTUL
      J = INDEX(TEMP,'     ')
      TEMP(J:J) = ','
      MV = 0

!
!     Decompose the input of Virtual Set
!
 700  IF ( TEMP(:5).NE.'     ' ) THEN
         MV = MV + 1
         IF ( MV.GT.(NELS) ) THEN
            WRITE (0,*) ' Virtual set too large: MAX=', NELS
            STOP
         END IF
         CALL DEL(TEMP)
         J = INDEX(TEMP,',')
         VL(MV) = MOD(LVAL(TEMP(2:2)),2)
!
!        Convert the input of uppercase to lowercase and assign value
!        ELVi
         N = ICHAR(TEMP(2:2))
         IF ( N.GE.ORDUA .AND. N.LE.ORDUZ ) TEMP(2:2)
     &        = CHAR(N-ORDUA+ORDLA)
         ELV(MV) = TEMP(:J-1)
         CALL STRSH(TEMP,J+1)
         GO TO 700
      END IF

 800  WRITE (0,99008)
99008 FORMAT ('  From which shell  ?  ')
      READ (5,'(A)') TEMP
      I = INDEX(TEMP,'B')
      J = INDEX(TEMP,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 600
      CALL DEL(TEMP)
      NVIR = ICHAR(TEMP(1:1)) - ORD0
      IF ( NVIR.LT.1 .OR. NVIR.GT.9 ) THEN
         WRITE (0,*) '    Please input digit for the shell position!'
         GO TO 800
      END IF

 900  WRITE (0,99009)
99009 FORMAT ('    To which shell  ?  ')
      READ (5,'(A)') TEMP
      I = INDEX(TEMP,'B')
      J = INDEX(TEMP,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 600
      CALL DEL(TEMP)
      LVIR = ICHAR(TEMP(:1)) - ORD0
      IF ( LVIR.LT.0 .OR. LVIR.GT.9 ) THEN
         WRITE (0,*) '     Please input digit for the shell position!'
         GO TO 900
      END IF

!******************
!     Input Final Terms
!
!     NFTM  =  Number of input for final term
!
 1000 WRITE (0,99010)
99010 FORMAT ('       Final Terms  ?  ')
      IFLAG = 0
      CALL INPUT(NELS,NFTM,FINALS,IFLAG,SFLAG,DFLAG,*600,*1100)
      DO I = 1, NFTM
         CH2 = FINALS(I)(:2)

!
!        Conver the input of lowercase to uppercase
!
         N = ICHAR(CH2(2:2))
         IF ( N.GE.ORDLA .AND. N.LE.ORDLZ ) CH2(2:2)
     &        = CHAR(N-ORDLA+ORDUA)
         FTM(I) = CH2
      END DO

!**********************************************************************
!     Open the file for CI.LST
!
 1100 OPEN (7,FILE=FILENAME,STATUS='UNKNOWN')

!**********************************************************************
!     PRINT OUT ALL INPUT DATA OF THE USER
!
      WRITE (6,99011)
99011 FORMAT (///////T5,'***************          I N P U T',
     &        '  D A T A          **********'///)
      WRITE (6,99012) HEADER
99012 FORMAT (T5,'          Header  :  ',A60/)
      WRITE (6,99013) SHELLS
99013 FORMAT (T5,'    Closed shell  :  ',A60/)
      WRITE (6,99014) REF(1)
99014 FORMAT (T5,'   Reference Set  :  ',A60/)
      DO I = 2, NREF
         WRITE (6,99033) I, REF(I)
      END DO
      WRITE (6,99015) ACT
99015 FORMAT (T5,'      Active Set  :  ',A60/)
      WRITE (6,99016) REPL(1)
99016 FORMAT (T5,'    Replacements  :  ',A60/)
      IF ( SFLAG.NE.0 .OR. DFLAG.NE.0 ) THEN
         WRITE (6,99017) VIRTUL
99017    FORMAT (T5,'     Virtual Set  :  ',A72/)
         WRITE (6,*) '             From which shell  :  ',
     &               CHAR(NVIR+ORD0)
         WRITE (6,*)
         WRITE (6,*) '               To which shell  :  ',
     &               CHAR(LVIR+ORD0)
         WRITE (6,*)
      END IF
      DO I = 2, NREPL
         WRITE (6,99033) I, REPL(I)
      END DO
      WRITE (6,99018) FINALS(1)
99018 FORMAT (T5,'     Final Terms  :  ',A60/)
      DO I = 2, NFTM
         WRITE (6,99033) I, FINALS(I)
      END DO

      WRITE (7,99019) HEADER
99019 FORMAT (' ',A60)
      WRITE (7,'(1X,A72)') SHELLS

!**********************************************************************
!     PROCESS THE REFERENCE SET
!
!     M  =  Maximin number of ELi for reference set
!     ELS  =  Storage of ELi for different Reference Set
!     QS  =  Storage of Qi for each Reference Set
!     MS  =  Storage of M for each Reference Set
!     STRL  =  Temporary string
!
!
      DO NF = 1, NREF
         M = 0
         TEMP = REF(NF)
         DO I = 1, NELS

!**************
!           Decompose the input for Reference Set
!
            IF ( TEMP(:5).EQ.'     ' ) THEN
               Q(I) = 0
               EXIT
            END IF
            CALL DEL(TEMP)
            LEFT = INDEX(TEMP,'(')
            RIGHT = INDEX(TEMP,')')
            M = M + 1

!**************
!           Convert the input of uppercase to lowercase, and assign
!           initial values for ELi and Qi
!
            CH3 = TEMP(:LEFT-1)
            N = ICHAR(CH3(2:2))
            IF ( N.GE.ORDUA .AND. N.LE.ORDUZ ) CH3(2:2)
     &           = CHAR(N-ORDUA+ORDLA)
            EL(M) = CH3
            CH3 = TEMP(LEFT+1:RIGHT-1)
            Q(M) = ICTOI(CH3)
            CALL STRSH(TEMP,RIGHT+1)
         END DO
         DO J = 1, NELS
            QS(J,NF) = 0
         END DO

!**************
!        Add new electrons to the first configuration
!
         IF ( NF.EQ.1 ) THEN
            DO J = 1, M
               ELS(J) = EL(J)
               QS(J,NF) = Q(J)
            END DO
            MS(NF) = M
         ELSE
            DO J = 1, M
               N = MS(1)
               DO I = 1, N
                  IF ( EL(J).EQ.ELS(I) ) THEN
                     QS(I,NF) = Q(J)
                     GO TO 1120
                  END IF
               END DO
               N = N + 1
               IF ( N.GT.NELS ) THEN
                  WRITE (0,*) ' Too many shells in reference set: MAX=',
     &                        NELS
                  STOP
               END IF
               MS(1) = N
               ELS(N) = EL(J)
               QS(N,NF) = Q(J)
 1120       END DO
            MS(NF) = N
         END IF

!***************
!        Generate all couplings for the reference set
!        MAX = -5  means the first time to call subroutine COUPLD
!
         K = 0
         DO I = 1, M
            IF ( Q(I).NE.0 ) THEN
               K = K + 1
               ELB(K) = EL(I)
               QB(K) = Q(I)
            END IF
         END DO

         MAXFINAL = -5
         CALL COUPLD(ELB,QB,K,ELC,NC,*1700)
         IF ( NC.NE.0 ) THEN
            WRITE (6,99020)
99020       FORMAT (//
     &'      GENERATE ALL COUPLINGS FOR EACH MEMBER  OF THE REFERENCE SE
     &T'//)
            CALL PRINT_COUPLINGS(ELC,QB,K,NC,1)
         END IF

!**************
!        Compute max and min value for the finals.
!        Rule :  MAXFINAL = 2*|S+L| ,  MINFINAL = 2*|S-L|
!
         IF ( NF.EQ.1 ) THEN
            MAXFINAL = 0
            MINFINAL = 100
            DO I = 1, NC
               N = 2*K - 1
               READ (FILE3(I),99021) (COUPLE(J),J=1,N)
99021          FORMAT (9(A3))
               CH3 = COUPLE(N)(1:1)
               CH1 = COUPLE(N)(2:2)
               S = (ICTOI(CH3)-1)/2.
               L = LVAL(CH1)
               LMIN = 2*ABS(S-L)
               LMAX = 2*ABS(S+L)
               IF ( LMIN.LT.MINFINAL ) MINFINAL = LMIN
               IF ( LMAX.GT.MAXFINAL ) MAXFINAL = LMAX
            END DO
         END IF
      END DO
      IF ( ACT(:5).EQ.'     ' ) GO TO 1500

!**********************************************************************
!     PROCESS THE ACTIVE SET
!
!     The member of the input is limited to NELS .
!
!     NQ  =  Number of configurations
!     ELA  =  ELi for Active Set
!     QA  =  Qi for Active Set
!     MA  =  Number of ELi for Active Set
!     RL  =  L-value for each shell
!     F  =  Full value for each shell
!
      MA = 0
      TEMP = ACT

!**************
!     Decompose the input of Active Set, Conver the input of uppercase
!     to lowercase, assign values to ELA,QA,RL,F .
!
 1200 IF ( MA.LT.NELS ) THEN
         IF ( TEMP(1:5).EQ.'     ' ) GO TO 1300
         CALL DEL(TEMP)
         MA = MA + 1
         N = ICHAR(TEMP(2:2))
         IF ( N.GE.ORDUA .AND. N.LE.ORDUZ ) TEMP(2:2)
     &        = CHAR(N-ORDUA+ORDLA)
         ELA(MA) = TEMP(:2)
         RL(MA) = LVAL(TEMP(2:2))
         IF ( RL(MA).GT.3 ) THEN
            F(MA) = 2
         ELSE
            F(MA) = MIN0(CONST,4*RL(MA)+2)
         END IF
         CALL STRSH(TEMP,4)
         GO TO 1200
      ELSE
         WRITE (0,*) 'Too many electrons in the active set: MAX= 15'
         STOP
      END IF

!***************
!     NELS is maximun value of active electrons
!
 1300 DO I = MA + 1, NELS
         F(I) = 0
      END DO

!**************
!     Generate other possible configurations from the given set
!     F(I) is the maximum allowed number of electrons
!     H(I) is the number of unassigned electrons or ``holes''
!
 1400 IF ( ITYPE.EQ.0 ) THEN
         MINQ1 = 0
      ELSE IF ( ITYPE.EQ.1 ) THEN
         MINQ1 = MAX0(MIN0(F(1),CONST)-1,0)
      ELSE IF ( ITYPE.EQ.2 ) THEN
         MINQ1 = MAX0(MIN0(F(1),CONST)-2,0)
      ELSE IF ( ITYPE.EQ.3 ) THEN
         MINQ1 = MAX0(MIN0(F(1),CONST)-3,0)
      ELSE
         WRITE (0,*) ' Unknown type: Re-enter'
         READ *, ITYPE
         GO TO 1400
      END IF
      DO Q1 = MIN0(F(1),CONST), MINQ1, -1
         H(1) = MAX0(0,CONST-Q1)
         DO Q2 = MIN0(F(2),H(1)), 0, -1
            H(2) = MAX0(0,H(1)-Q2)
            DO Q3 = MIN0(F(3),H(2)), 0, -1
               H(3) = MAX0(0,H(2)-Q3)
               DO Q4 = MIN0(F(4),H(3)), 0, -1
                  H(4) = MAX0(0,H(3)-Q4)
                  DO Q5 = MIN0(F(5),H(4)), 0, -1
                     H(5) = MAX0(0,H(4)-Q5)
                     DO Q6 = MIN0(F(6),H(5)), 0, -1
                        H(6) = MAX0(0,H(5)-Q6)
                        DO Q7 = MIN0(F(7),H(6)), 0, -1
                           H(7) = MAX0(0,H(6)-Q7)
                           DO Q8 = MIN0(F(8),H(7)), 0, -1
                              H(8) = MAX0(0,H(7)-Q8)
                              DO Q9 = MIN0(F(9),H(8)), 0, -1
                                 H(9) = MAX0(0,H(8)-Q9)
                                 DO Q10 = MIN0(F(10),H(9)), 0, -1
                                    H(10) = MAX0(0,H(9)-Q10)
                                    DO Q11 = MIN0(F(11),H(10)), 0, -1
                                       H(11) = MAX0(0,H(10)-Q11)
                                       DO Q12 = MIN0(F(12),H(11)), 0, -1
                                         H(12) = MAX0(0,H(11)-Q12)
                                         DO Q13 = MIN0(F(13),H(12)), 0,
     &                                      -1
                                         H(13) = MAX0(0,H(12)-Q13)
                                         DO Q14 = MIN0(F(14),H(13)), 0,
     &                                      -1
                                         H(14) = MAX0(0,H(13)-Q14)
                                         IF ( H(14).LE.F(15) ) THEN
                                         Q15 = H(14)
                                         CALL CONFIG
                                         END IF
                                         END DO
                                         END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
!**************
!     Print Header of the output file
!
      IF ( NQ.EQ.0 ) GO TO 1500
      WRITE (6,99022)
99022 FORMAT (//
     &'      GENERATE ALL POSSIBLE CONFIGURATIONS FROM THE   ACTIVE SET'
     &//)
      J = INDEX(HEADER,'         ')
      WRITE (6,99023) HEADER(1:J)
99023 FORMAT (' '/T10,'-------------       ',A,'  --------'/)
      WRITE (6,99024) SHELLS
99024 FORMAT ('          Closed Shells  :  ',A60/)
      WRITE (6,99025) ACT
99025 FORMAT ('             Active Set  :  ',A60/)

!**************
!     Print all configurations generated from Active Set
!
      WRITE (6,*) '         Configurations  :'
      DO I = 1, NQ
         READ (FILE1(I),99034) (QA(J),J=1,MA)
         IF ( ELA(1)(1:1).EQ.' ' ) THEN
            N = 27
         ELSE
            N = 28
         END IF

         K = 0
         DO J = 1, MA
            IF ( QA(J).NE.0 .OR. MA.EQ.2 ) THEN
               K = K + 1
               ELB(K) = ELA(J)
               QB(K) = QA(J)
            END IF
         END DO
         WRITE (6,99026) (ELB(J),QB(J),J=1,K)
99026    FORMAT (T28,8(1X,A3,'(',I2,')'))
      END DO

!***************
!     For each new configuration, generate all possible couplings
!
      N1 = 0
      N2 = 0
      DO I = 1, NQ
         READ (FILE1(I),99034) (QA(J),J=1,MA)

!**************
!        Omit ELA(i) if corresponding QA(i)=0
!
         K = 0
         DO J = 1, MA
            IF ( QA(J).NE.0 .OR. MA.EQ.2 ) THEN
               K = K + 1
               ELB(K) = ELA(J)
               QB(K) = QA(J)
            END IF
         END DO

!***************
!        Omit configurations which have more than 8 shells
!
         IF ( K.LE.8 ) THEN
            CALL COUPLD(ELB,QB,K,ELC,NC,*1700)
            IF ( NC.GT.0 ) THEN
               WRITE (6,99027)
99027          FORMAT (//'    LIST THE COUPLINGS FOR EACH',
     &                 ' CONFIGURATION GENERATED BY THE ACTIVE SET'//)
               CALL PRINT_COUPLINGS(ELC,QB,K,NC,2)
            ELSE
               N2 = N2 + 1
            END IF
         ELSE
            N1 = N1 + 1
         END IF
      END DO
      IF ( N1.NE.0 ) PRINT 99028, N1
99028 FORMAT (T5,'Too many occuplied shells --- ',I3,
     &        ' configuration omitted!')
      IF ( N2.NE.0 ) PRINT 99029, N2
99029 FORMAT (T5,'No final term as your selection for ',I3,
     &        ' Active set!')

      IF ( NREPL.EQ.0 ) GO TO 1700

!**********************************************************************
!     PROCESS THE REPACEMENTS
!
!     STRL  :  String of the left side of '='
!     STRR  =  String of the right side of '='
!

 1500 IF ( SFLAG.NE.0 .OR. DFLAG.NE.0 ) GO TO 1600
      N1 = 0
      DO NR = 1, NREPL
         STRR = REPL(NR)
         J = INDEX(STRR,'=')
         STRL = STRR(1:J-1)
         CALL STRSH(STRR,J+1)

!**************
!        Decompose the substring on the left of ''=''
!
!        ELL  =  Old value of ELi to be replaced
!        QL  =  Old value of Qi to be replaced
!        ML  =  Number of Old value of ELi to be replaced
!
         CALL DEL(STRL)
         CALL DECOMP(STRL,ELL,QL,ML)

!**************
!        Decompose the substring on the right of ''=''
!
!        ELR  =  New value of ELi
!        QR  =  New vaRue of Qi
!        MR  =  Number of the new value of ELi
!
         CALL DEL(STRR)
         CALL DECOMP(STRR,ELR,QR,MR)

!**************
!        For each Replacement, replace all Reference Sets
!
         DO NF = 1, NREF
            M = MS(NF)
            IF ( M.LT.ML ) CYCLE
            DO I = 1, M
               EL(I) = ELS(I)
               Q(I) = QS(I,NF)
            END DO
            CALL REPLAC(ELB,QB,K,*1550)
            CALL COUPLD(ELB,QB,K,ELC,NC,*1700)
            IF ( NC.GT.0 ) THEN
               WRITE (6,99030)
99030          FORMAT (//'      FOR EACH REPLACEMENT, GENERATE',
     &             ' CONFIGURATIONS AND COUPLINGS FOR THE REFERENCE SET'
     &             //)
               CALL PRINT_COUPLINGS(ELC,QB,K,NC,3)
            ELSE
               N1 = N1 + 1
            END IF
 1550    END DO
      END DO
      IF ( N1.NE.0 ) WRITE (0,99031) N1
99031 FORMAT (T5,'No Final Term as selected for ',I3,' Replacements!')
      GO TO 1700

!**********************************************************************
!     PROCESS THE VIRTUAL SET
!
!     RL  =  Parity of ELi of Reference Set
!     VL  =  Parity of ELi of Virtual Set
!
!
 1600 DO NF = 1, NREF

!**************
!        Set the initial values of Reference Set
!
         M = MS(NF)
         DO I = 1, NELS
            CH3 = ELS(I)
            EL(I) = CH3
            ELB(I) = CH3
            Q(I) = QS(I,NF)
            QB(I) = Q(I)
            RL(I) = MOD(LVAL(CH3(2:2)),2)
         END DO

         IF ( SFLAG.EQ.0 .AND. DFLAG.NE.0 ) GO TO 1650

!**************
!        Preparation for the Single Replacement
!
         QL(1) = 1
         ML = 1
         QR(1) = 1
         MR = 1
         DO I = NVIR, LVIR
            IF ( Q(I).EQ.0 ) CYCLE
            ELL(1) = EL(I)
            IF ( RL(I).GT.0 ) THEN
               N = 1
            ELSE
               N = 0
            END IF
            DO J = 1, MV

!**************
!              Replace Qi of Reference Set which has the same parity
!              with Qj of Virtual Set
!
               IF ( VL(J).EQ.N ) THEN
                  ELR(1) = ELV(J)
                  CALL REPLAC(ELB,QB,K,*1620)
                  CALL COUPLD(ELB,QB,K,ELC,NC,*1700)
                  IF ( NC.GT.0 ) THEN
                     WRITE (6,99032)
99032                FORMAT (//'        FOR VIRTUAL SET, ',
     &          'GENERATE CONFIGURATION AND COUPLINGS FOR S-REPLACEMENT'
     &          //)
                     NR = NELS
                     TEMP = ELL(1)//' = '//ELR(1)
                     CALL DEL(TEMP)
                     REPL(NR) = TEMP
                     CALL PRINT_COUPLINGS(ELC,QB,K,NC,4)
                  END IF
               END IF
 1620       END DO
         END DO


!*************
!        Rreplace pairs of Q(i) and Q(j) by Double Virtual Set
!
!        PL  =  Pairty for the pair Qi and Qj in Reference Set
!
 1650    IF ( DFLAG.EQ.0 ) EXIT
         ML = 2
         QL(1) = 1
         QL(2) = 1
         DO I = NVIR, LVIR - 1
            ELL(1) = EL(I)
            LL1 = LVAL(ELL(1)(2:2))
            DO J = I + 1, M
               IF ( Q(I).EQ.0 .OR. Q(J).EQ.0 ) CYCLE
               ELL(2) = EL(J)
               LL2 = LVAL(ELL(2)(2:2))
               TEMP = ELL(1)//'.'//ELL(2)//' = '
               LLMIN = IABS(LL1-LL2)
               LLMAX = LL1 + LL2
               PL = MOD(LLMAX,2)
               CALL VPAIR(ELV,MV,PL,LLMIN,LLMAX,TEMP,*1700)
            END DO
         END DO

!**************
!        Replace pairs of (Qi)=2 by Double Virtual Set
!
         ML = 1
         QL(1) = 2
         DO I = NVIR, LVIR
            IF ( Q(I).GT.1 ) THEN
               ELL(1) = EL(I)
               LL1 = LVAL(ELL(1)(2:2))
               LLMIN = 0
               LLMAX = LL1 + LL1
               TEMP = ELL(1)//'(2) = '
               PL = MOD(LLMAX,2)
               CALL VPAIR(ELV,MV,PL,LLMIN,LLMAX,TEMP,*1700)
            END IF
         END DO
      END DO

!**********************************************************************
!     THE END OF THE PROGRAM
!
 1700 WRITE (7,'(A1)') '*'
      CLOSE (7)
      WRITE (0,*)
      WRITE (0,*) '         OK !'
      WRITE (0,*) '         List of configurations and their couplings'
      WRITE (0,*) '         is in the file ', FILENAME
99033 FORMAT (T29,I2,'  :  ',A60/)
99034 FORMAT (21(I2))
      END PROGRAM GENCLF
!
!
!                  S U B R O U T I N E S
!                  =====================
!
!
! ----------------------------------------------------------------------
!               SUBROUTINE      C O N F I G
! ----------------------------------------------------------------------
!
!               Examine if the new configuration has the same
!       electron number and parity .
!
      SUBROUTINE CONFIG
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS, NCFG, NSCOUP
      PARAMETER (NELS=21,NCFG=10000,NSCOUP=10000)
!
! COMMON variables
!
      INTEGER CONST, M, MA, MAXFINAL, MINFINAL, ML, MR, MS(NELS), NF,
     &        NFTM, NQ, NR, NREF, PARITYVAL, Q(NELS), Q1, Q10, Q11, Q12,
     &        Q13, Q14, Q15, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, QA(NELS),
     &        QL(NELS), QR(NELS), QS(NELS,NELS), RL(NELS)
      CHARACTER*3 EL(NELS), ELA(NELS), ELL(NELS), ELR(NELS), ELS(NELS)
      CHARACTER*8 FBETA(7,NSCOUP)
      CHARACTER*32 FILE1(NCFG), FILE3(NSCOUP)
      CHARACTER*40 FILE2(NCFG)
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK3  / EL, ELL, ELR, ELS, ELA
      COMMON /BLK4  / Q, QL, ML, QR, MR, M, QS, MS, MA, RL, NREF, Q1,
     &                Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12,
     &                Q13, Q14, Q15
      COMMON /FILES / FBETA, FILE1, FILE2, FILE3
!
! Local variables
!
      CHARACTER*3 CH3
      INTEGER I, IQ, J, K, N, NEWP
!
!
      EQUIVALENCE (Q1,QA(1))
      EQUIVALENCE (Q2,QA(2))
      EQUIVALENCE (Q3,QA(3))
      EQUIVALENCE (Q4,QA(4))
      EQUIVALENCE (Q5,QA(5))
      EQUIVALENCE (Q6,QA(6))
      EQUIVALENCE (Q7,QA(7))
      EQUIVALENCE (Q8,QA(8))
      EQUIVALENCE (Q9,QA(9))
      EQUIVALENCE (Q10,QA(10))
      EQUIVALENCE (Q11,QA(11))
      EQUIVALENCE (Q12,QA(12))
      EQUIVALENCE (Q13,QA(13))
      EQUIVALENCE (Q14,QA(14))
      EQUIVALENCE (Q15,QA(15))

!     INPUT :
!     Q,Qi  =  Occupation number
!     L  =  L-value corresponding ELi
!     OUTPUT :
!     NQ  =  Number of configurations

!**************
!     Return if the new configuration has different pairty
!
      NEWP = 0
      DO I = 1, MA
         NEWP = NEWP + QA(I)*RL(I)
      END DO

      NEWP = MOD(NEWP,2)
      IF ( NEWP.NE.PARITYVAL ) RETURN

!**************
!     Return if the new condiguration is the same as the Reference Set
!
      DO I = 1, NREF
         M = MS(I)
         N = 0
         DO J = 1, M
            CH3 = ELS(J)
            IQ = QS(J,I)
            DO K = 1, MA
               IF ( CH3.EQ.ELA(K) .AND. IQ.EQ.QA(K) ) N = N + 1
            END DO
         END DO
         IF ( N.EQ.M ) RETURN
      END DO

!**************
!     Otherwise, write them into the configuration file
!
      NQ = NQ + 1
      WRITE (FILE1(NQ),99001) (QA(J),J=1,MA)
99001 FORMAT (15(I2))
      CONTINUE
      END SUBROUTINE CONFIG
!
! ----------------------------------------------------------------------
!               SUBROUTINE      C O U P L D
! ----------------------------------------------------------------------
!
!          This subroutine generates all possible couplings .
!       First, compute Alpha value from the given configuration,
!       then compute Beta from each value of Alpha .

      SUBROUTINE COUPLD(EL,Q,M,ELC,NC,*)

!   INPUT :
!       El  =  electron label
!              where EL(I)(1=1)  ---  blank
!                    EL(I)(2=2)  ---  n-symbol
!                    EL(I)(3=3)  ---  L-symbol
!       Q  =  Occupation number
!               0 (empty) <= Q(i)  <= 2(2Li+1) (full)
!       M  =  number of shells
!               0  <=  M  <=  5
!   OUTPUT :
!       NC  =  number of couplings
!        *  =  Return label if the maximun number of couplings > NSCOUP
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS, NSHEL, NCOUPL, NCFG, NSCOUP
      PARAMETER (NELS=21,NSHEL=8,NCOUPL=2*NSHEL-1,NCFG=10000,
     &           NSCOUP=10000)
!
! COMMON variables
!
      INTEGER CONST, MAXFINAL, MINFINAL, NF, NFTM, NQ, NR, ORD0, ORD9,
     &        ORDLA, ORDLZ, ORDUA, ORDUZ, PARITYVAL
      CHARACTER*8 FBETA(7,NSCOUP)
      CHARACTER*32 FILE1(NCFG), FILE3(NSCOUP)
      CHARACTER*40 FILE2(NCFG)
      CHARACTER*2 FTM(NELS)
      CHARACTER*60 REF(NELS), REPL(NELS)
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK0  / ORDLA, ORDLZ, ORDUA, ORDUZ, ORD0, ORD9
      COMMON /BLK2  / REF, REPL, FTM
      COMMON /FILES / FBETA, FILE1, FILE2, FILE3
!
! Dummy arguments
!
      INTEGER M, NC
      CHARACTER*3 EL(*), ELC(*)
      INTEGER Q(*)
!
! Local variables
!
      CHARACTER*2 A2, B1, CH2
      CHARACTER*3 ALFA(NSHEL,NSCOUP), CCH3, CH3, COUPLE(NCOUPL,NSCOUP)
      INTEGER BETA(NELS), CHILD, I, J, JJ, K, KS, L, L1, L2, L3, L4,
     &        LMAX, LMIN, LOCA, LOCB, LOCT, LPOSIT(10), MBETA, N, NALFA,
     &        NB, NBETA, NT, NTERM(27), PARENT, POSIT(NELS), PTR
      CHARACTER*15 CALFA
      CHARACTER CCH1, CH1
      DOUBLE PRECISION FULL, HALF, S, S1, S2, S3, S4
      INTEGER ICTOI, LVAL
      CHARACTER SYMB
      CHARACTER*65 TERM(27), TERM2(27), TERM3(27), TERM4(27), TERM5(27),
     &             TERM6(27), TERMQ3
      EXTERNAL ICTOI, LVAL, SYMB
!
!
!**************
!     Number of possible terms for configurations P(1-3),
!     D(1-5), F(1-2), G(1-2), ... M(1-2)
!
      DATA (NTERM(I),I=1,27)/1, 3, 3, 1, 5, 8, 16, 16, 1, 7, 17, 47, 73,
     &      119, 119, 1, 9, 1, 9, 1, 9, 1, 9, 1, 9, 1, 9/
!**************
!     Starting position in term table for given L
!
      DATA (LPOSIT(I),I=1,9)/1, 4, 9, 16, 18, 20, 22, 24, 26/

!**************
!     Possible terms for configurations P(1-3),D(1-5),F(1-2),G(1-2)
!
      DATA (TERM(I),I=1,27)/'2P1', '1S01D23P2', '2P12D34S3', '2D1',
     &      '1S01D21G23P23F2', '2D12P32D32F32G32H34P34F3',
     &      '1S01D21G23P23F21S41D41F41G41I43P43D43F43G43H45D4',
     &      '2D12P32D32F32G32H34P34F32S52D52F52G52I54D54G56S5', '2F1',
     &      '1S01D21G21I23P23F23H2',
     &      '2P12D12D22F12F22G12G22H12H22I12K12L14S14D14F14G14I1',
     & '1S11S21D11D21D31D41F11G11G21G31G41H11H21I11I21I31K11L11L21N13P1'
     & ,
     & '2P12P22P32P42D12D22D32D42D52F12F22F32F42F52F62F72G12G22G32G42G5'
     & ,
     & '1S11S21S31S41P01D11D21D31D41D51D61F11F21F31F41G11G21G31G41G51G6'
     & ,
     & '2S12S22P12P22P32P42P52D12D22D32D42D52D62D72F12F22F32F42F52F62F7'
     & , '2G1', '1S01D21G21I21L23P23F23H23K2', '2H1',
     & '1S01D21G21I21L23P23F23H23K2', '2I1',
     & '1S01D21G21I21L23P23F23H23K2', '2K1',
     & '1S01D21G21I21L23P23F23H23K2', '2L1',
     & '1S01D21G21I21L23P23F23H23K2', '2M1',
     & '1S01D21G21I21L23P23F23H23K2'/
      DATA (TERM2(I),I=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR',
     &      'ERROR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',
     & '3P23P33D13D23F13F23F33F43G13G23G33H13H23H33H43I13I23K13K23L13M1'
     & ,
     & '2G62H12H22H32H42H52H62H72I12I22I32I42I52K12K22K32K42K52L12L22L3'
     & ,
     & '1G71G81H11H21H31H41I11I21I31I41I51I61I71K11K21K31L11L21L31L41M1'
     & ,
     & '2F82F92FA2G12G22G32G42G52G62G72G82G92GA2H12H22H32H42H52H62H72H8'
     & , 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',
     & 'ERR', 'ERROR', 'ERR', 'ERROR'/
      DATA (TERM3(I),I=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR',
     &      'ERROR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',
     &      '5S05D15F15G15I1',
     & '2M12M22N12O04S14P14P24D14D24D34F14F24F34F44G14G24G34G44H14H24H3'
     & ,
     & '1M21N11N21Q03P13P23P33P43P53P63D13D23D33D43D53F13F23F33F43F53F6'
     & ,
     & '2H92I12I22I32I42I52I62I72I82I92K12K22K32K42K52K62K72L12L22L32L4'
     & , 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',
     & 'ERR', 'ERROR', 'ERR', 'ERROR'/
      DATA (TERM4(I),I=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR',
     &      'ERROR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR',
     &      '4I14I24I34K14K24L14M06P06F06H0',
     & '3F73F83F93G13G23G33G43G53G63G73H13H23H33H43H53H63H73H83H93I13I2'
     & ,
     & '2L52M12M22M32M42N12N22O02Q04S14S24P14P24D14D24D34D44D54D64F14F2'
     & , 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',
     & 'ERR', 'ERROR', 'ERR', 'ERROR'/
      DATA (TERM5(I),I=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR',
     &      'ERROR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR',
     &      'ERROR',
     & '3I33I43I53I63K13K23K33K43K53K63L13L23L33M13M23M33N03O07F05S05P0'
     & ,
     & '4F34F44F54G14G24G34G44G54G64G74H14H24H34H44H54I14I24I34I44I54K1'
     & , 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',
     & 'ERR', 'ERROR', 'ERR', 'ERROR'/
      DATA (TERM6(I),I=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR',
     &      'ERROR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR',
     &      'ERROR', '5D15D25D35F15F25G15G25G35H15H25I15I25K05L0',
     &      '4K24K34L14L24L34M04N06P06D06F06G06H06I08S0', 'ERR',
     &      'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',
     &      'ERR', 'ERROR', 'ERR', 'ERROR'/
      DATA TERMQ3/'1S11D11G11I13P13F13H1'/
!**********************************************************************
!     COMPUTE THE POSSIBLE VALUES FOR ALPHA
!
!     EMPTY  =  1 when Q(i)=0
!     FULL  =   2*(2L+1)
!     HALF  =  2L+1
!     ALFA  =  matrix of (NELS,NSCOUP)
!     NALFA  =  Number of ALFA
!     POSIT  =  array of NELS, store the position in table NTERM and
!     TERM corresponding Q(i). Rule  :
!     position = (L-1)*2+Q(i)        if 1 <= Qi <= HALF
!     position = (L-1)*2+(FULL-Qi)   if Qi > HALF
!

      NALFA = 1
      DO I = 1, M
         CH3 = EL(I)
         CH1 = CH3(2:2)
         IF ( CH3(3:3).EQ.' ' ) THEN
            CCH3 = CH3
            CH3 = ' '//CCH3(1:2)
         END IF
         ELC(I) = CH3
         CCH1 = CHAR(ICHAR(CH1)-ORDLA+ORDUA)
         CH1 = CCH1
         L = LVAL(CH1)
         FULL = 4*L + 2
         K = Q(I)

!**************
!        If shell is full, ALFA(i) = 1S0
!
         IF ( K.EQ.0 .OR. K.EQ.INT(FULL) ) THEN
            CH3 = '1S0'

!**************
!           If Q(i) = 1,  then ALFA(i)=2<L-symbol>1
!
         ELSE IF ( K.EQ.1 ) THEN
            CH3 = '2'//CH1//'1'

!**************
!           Otherwise, get the possible value from array NTERM and TERM
!
         ELSE
            HALF = FULL/2
            IF ( K.LE.HALF ) THEN
               POSIT(I) = LPOSIT(L) + K - 1
            ELSE
               POSIT(I) = LPOSIT(L) + (FULL-K) - 1
            END IF
            NALFA = NALFA*NTERM(POSIT(I))
            CH3 = '   '
         END IF
         IF ( NALFA.GT.NSCOUP ) THEN
            WRITE (0,*) 'Array ALFA in routine COUPLD exceeded'
            STOP
         END IF

!**************
!        CALFA is a string storing the elements in one ALFA
!
         CALFA(3*I-2:3*I) = CH3
      END DO

!**************
!     Assign values to all elements of ALFA
!
!     NT  =
!     LOCT  =  Current position in the table TERM
!     LOCA  =  Current position in the matrix ALFA
!
      NT = 1
      DO I = M, 1, -1
         CH3 = CALFA(I*3-2:I*3)
         IF ( CH3.NE.'   ' ) THEN
            DO J = 1, NALFA
               ALFA(I,J) = CH3
            END DO
         ELSE
            LOCT = POSIT(I)
            N = NTERM(LOCT)
            LOCA = 1
 20         IF ( LOCA.LE.NALFA ) THEN
               DO J = 1, N
!                 CH3 = TERM(LOCT)(J*4-3:J*4-1)
                  IF ( LOCT.EQ.10 .AND. Q(I).EQ.12 ) THEN
!                    The Case f(12)
                     CH3 = TERMQ3(J*3-2:J*3)
!
                  ELSE
!                    Other Cases
                     IF ( J.LT.22 ) THEN
                        CH3 = TERM(LOCT)(J*3-2:J*3)
                     ELSE IF ( J.LT.43 ) THEN
                        JJ = J - 21
                        CH3 = TERM2(LOCT)(JJ*3-2:JJ*3)
                     ELSE IF ( J.LT.64 ) THEN
                        JJ = J - 42
                        CH3 = TERM3(LOCT)(JJ*3-2:JJ*3)
                     ELSE IF ( J.LT.85 ) THEN
                        JJ = J - 63
                        CH3 = TERM4(LOCT)(JJ*3-2:JJ*3)
                     ELSE IF ( J.LT.106 ) THEN
                        JJ = J - 84
                        CH3 = TERM5(LOCT)(JJ*3-2:JJ*3)
                     ELSE IF ( J.LT.120 ) THEN
                        JJ = J - 105
                        CH3 = TERM6(LOCT)(JJ*3-2:JJ*3)
                     END IF
                  END IF
                  DO K = 1, NT
                     ALFA(I,LOCA) = CH3
                     LOCA = LOCA + 1
                  END DO
               END DO
               GO TO 20
            END IF
            NT = NT*N
         END IF
      END DO

!**********************************************************************
!     GENERATE POSSIBLE VALUE OF BETA FROM ALFA
!

      NC = 0
      DO NB = 1, NALFA

!**************
!        There is only one coupling if M = 1
!
         IF ( M.EQ.1 ) THEN
            BETA(1) = 1
            COUPLE(1,1) = ALFA(1,NB)
            NBETA = 1
            GO TO 100
         END IF

!**************
!        Define BETA(1)=ALFA(1), then the next basic coupling steps is :
!        S1 = (BETA(1)(1:1)-1)/2 ,   S2 = (ALFA(2)(1:1)-1)/2 ;
!        | S1-S2 | <= BETA(2)(1:1) <= | S1+S2 |
!        L1 = L-number of BETA(1)(2:2),   L2 = L-number of ALFA(2)(2:2)
!        ; Symbol(| L1-L2 |) <= BETA(2)(2:2) <= Symbol(| L1+L2 |) .
!
         B1 = ALFA(1,NB)(1:2)
         DO J = 2, M
            A2 = ALFA(J,NB)(1:2)
            PARENT = 1
            CHILD = 1
            BETA(J) = 0
 40         S1 = (ICHAR(B1(1:1))-ORD0-1)/2.
            S2 = (ICHAR(A2(1:1))-ORD0-1)/2.
            S3 = ABS(S1-S2)
            S4 = ABS(S1+S2)
            L1 = LVAL(B1(2:2))
            L2 = LVAL(A2(2:2))
            L3 = ABS(L1-L2)
            L4 = ABS(L1+L2)
            MBETA = (S4-S3+1)*(L4-L3+1)
!

!*************
!           Generate Beta from each alpha .
!           There are four scratch files for storing the information
!           about Beta(i), 1<i<6, shown as follows :
!           --------------------------------------------
!           |    Parent    |    Value    |    Child    |
!           --------------------------------------------
!           Define   Beta(i) is child of Beta(i-1) and parent of
!           beta(i+1) ; Parent is a pointer to the parent of Beta(i),
!           that is Beta(i-1) ; Value is one of the possible value for
!           BETA(i), and Child is the number of children for Beta(i) .
!
!           MBETA  =  Number of couplings generated from ALFA(i)
!           PARENT  =  Current pointer to the parent of Beta(i)
!
            DO S = S3, S4
               KS = 2*S
               CH1 = CHAR(ORD0+1+KS)
               DO L = L3, L4
                  WRITE (FBETA(J-1,CHILD),99004) PARENT, CH1//SYMB(L), 1
                  CHILD = CHILD + 1
               END DO
            END DO
            BETA(J) = BETA(J) + MBETA
            IF ( J.EQ.2 ) GO TO 60

!**************
!           Correct the number of its children for Beta(j-1)
!
            LOCB = PARENT
            DO K = J - 1, 2, -1
               READ (FBETA(K-1,LOCB),99004) PTR, CH2, N
               WRITE (FBETA(K-1,LOCB),99004) PTR, CH2, N + MBETA - 1
               LOCB = PTR
            END DO

!**************
!           If pointer to the parent is the end of file for Beta(j-1),
!           prepare to generate Beta(i+1) ; otherwise, generate next
!           Beta(j) according Alfa(j-1) and Beta(j-1) .
!
            PARENT = PARENT + 1
            IF ( PARENT.GT.BETA(J-1) ) THEN
               GO TO 60
            ELSE
               READ (FBETA(J-2,PARENT),99004) PTR, B1, N
               GO TO 40
            END IF
 60         READ (FBETA(J-1,1),99004) PTR, B1, N
         END DO

!**************
!        Assign values to the couplings forward ,
!        COUPLE(I) = Alpha(i) for coupling(j) if i <= M ;
!
         NBETA = BETA(M)
         IF ( NBETA.GT.NSCOUP ) THEN
            WRITE (0,*) 'Array COUPLE in routine COUPL exceeded'
            STOP
         END IF
         DO J = 1, NBETA
            DO I = 1, M
               COUPLE(I,J) = ALFA(I,NB)
            END DO
            DO K = 2*M - 1, 9
               COUPLE(K,J) = '   '
            END DO
         END DO

!**************
!        Assign values to the couplings backward ,
!        COUPLE(i) = Beta(i-m) for coupling(j) if i > m
!
         DO I = M, 2, -1
            N = 1
            DO J = 1, BETA(I)
               READ (FBETA(I-1,J),99004) PTR, CH2, NT
               DO K = 1, NT
                  COUPLE(M+I-1,N) = CH2//'0'
                  N = N + 1
                  IF ( N.GT.NSCOUP ) THEN
                     WRITE (0,*) 'Array COUPLE exceeded'
                     STOP
                  END IF
               END DO
            END DO
         END DO

!**************
!        Selection from generated couplings according the following
!        rules
 100     DO I = 1, NBETA
            N = 2*M - 1

!**************
!           If the first time to call COUPLD, not compute MAX and MIN
!
            IF ( MAXFINAL.EQ.-5 ) GO TO 120

!**************
!           Compute MAX and MIN value for each final term, keep it if
!           intersection is non-empty
!
            CH2 = COUPLE(N,I)(1:2)
            CH3 = CH2(1:1)
            CH1 = CH2(2:2)
            S = (ICTOI(CH3)-1)/2.
            L = LVAL(CH1)
            LMIN = 2*ABS(S-L)
            LMAX = 2*ABS(S+L)
            IF ( LMIN.GT.MAXFINAL .OR. LMAX.LT.MINFINAL ) CYCLE

!**************
!           If Final Terms are given, do selection
!
 120        IF ( NFTM.NE.0 ) THEN
               CH2 = COUPLE(N,I)(:2)
               DO K = 1, NFTM
                  IF ( CH2.EQ.FTM(K) ) GO TO 140
               END DO
               CONTINUE
               CYCLE
            END IF

!**************
!           Waining if the number of couplings > NSCOUP
!
            IF ( NC.EQ.NSCOUP ) THEN
               WRITE (0,*) '          WARNING !'
               WRITE (0,*) '          The number of couplings',
     &                     ' is greater than 2500 .',
     &                     '          Please select the Final Term .'
               RETURN 1
            END IF

!**************
!           Write configurations and couplings to CI.LST
!
 140        NC = NC + 1
            WRITE (FILE3(NC),99001) (COUPLE(J,I),J=1,N)
99001       FORMAT (9(A3))
            WRITE (7,99002) (ELC(J),Q(J),J=1,M)
99002       FORMAT (5(1X,A3,'(',I2,')'))
            WRITE (7,99003) (COUPLE(J,I),J=1,N)
99003       FORMAT (15(1X,A3))
         END DO
      END DO
      CONTINUE
99004 FORMAT (I3,A2,I3)
      END SUBROUTINE COUPLD
!
! ----------------------------------------------------------------------
!       SUBROUTINE      D E C O M P
! ----------------------------------------------------------------------
!
!       Decompose the string of Replacement
!
      SUBROUTINE DECOMP(STR,EL,Q,MR)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS
      PARAMETER (NELS=21)
!
! COMMON variables
!
      INTEGER ORD0, ORD9, ORDLA, ORDLZ, ORDUA, ORDUZ
      COMMON /BLK0  / ORDLA, ORDLZ, ORDUA, ORDUZ, ORD0, ORD9
!
! Dummy arguments
!
      INTEGER MR
      CHARACTER*60 STR
      CHARACTER*3 EL(NELS)
      INTEGER Q(NELS)
!
! Local variables
!
      CHARACTER*3 CH3
      INTEGER I, J, K, LEFT, N, RIGHT
      INTEGER ICTOI
      EXTERNAL DEL, ICTOI, STRSH
!
!
!     INPUT  :
!     STR  =  String to be decomposed
!     EL and Q
!     OUTPUT  :
!     MR  =  Number of EL to be replaced
!
      DO I = 1, 5
         IF ( STR(:5).EQ.'     ' ) THEN
            MR = I - 1
            RETURN
         END IF
         CALL DEL(STR)

         LEFT = INDEX(STR,'(')
         RIGHT = INDEX(STR,')')

!**************
!        If the Replacement is like 2s.2p = 3s.3p
!
         IF ( LEFT.EQ.0 ) THEN
            K = 1
            N = ICHAR(STR(3:3))
            IF ( N.GE.ORD0 .AND. N.LE.ORD9 ) THEN
               J = 3
            ELSE
               J = 2
            END IF
            CH3 = STR(:J)

!**************
!           If the Replacement is like 2p(2) = 3p(2)
!
         ELSE
            CH3 = STR(LEFT+1:RIGHT-1)
            K = ICTOI(CH3)
            CH3 = STR(:LEFT-1)
         END IF

!**************
!        Convert uppercase to lowercase, and assign value to ELi,Qi
!
         N = ICHAR(CH3(2:2))
         IF ( N.GE.ORDUA .AND. N.LE.ORDUZ ) CH3(2:2)
     &        = CHAR(N-ORDUA+ORDLA)
         EL(I) = CH3
         Q(I) = K
         IF ( LEFT.EQ.0 ) THEN
            CALL STRSH(STR,(J+2))
         ELSE
            CALL STRSH(STR,(RIGHT+1))
         END IF
      END DO
      CONTINUE
      END SUBROUTINE DECOMP
!
! ----------------------------------------------------------------------
!               FUNCTION        D E L
! ----------------------------------------------------------------------
!
!       Delete the leading space of the string
!
      SUBROUTINE DEL(STR)
      IMPLICIT NONE
!
!
! Dummy arguments
!
      CHARACTER*(*) STR
!
! Local variables
!
      INTEGER I, LENGTH
      CHARACTER*72 TEMP
!

      LENGTH = LEN(STR)
      I = 0
 100  IF ( STR(I+1:I+1).EQ.' ' ) THEN
         I = I + 1
         IF ( I.LT.LENGTH ) GO TO 100
      END IF
      TEMP = STR(I+1:)
      STR = TEMP
      CONTINUE
      END SUBROUTINE DEL
!
! ----------------------------------------------------------------------
!       SUBROUTINE      H E L P
! ----------------------------------------------------------------------
!
!       Explanation of the input format
!
      SUBROUTINE HELP
      IMPLICIT NONE
!
! Local variables
!
      INTEGER I, J
      CHARACTER*10 STR
!
 100  WRITE (0,99001)
99001 FORMAT (//,5X,
     &        'This program prompts for each required input.  The user'/
     &        5X,
     &        'should enter data or a RETURN after a question (?) mark')
      WRITE (0,*)
      WRITE (0,*) '     Example 1 :'
      WRITE (0,*) '    --------------------'
      WRITE (0,*) '                 Header  ?  S II ground state'
      WRITE (0,*) '          Closed shells  ?   1s 2s 2p'
      WRITE (0,*) '          Reference Set  ?  3s(2) 3p(3)'
      WRITE (0,*) '                      2  ?  RETURN'
      WRITE (0,*) '             Active Set  ?  3s,3p'
      WRITE (0,*) 'Type of set generation   ?  0'
      WRITE (0,*) '            Replacement  ?  3s(2) = 3d(2)'
      WRITE (0,*) '                      2  ?  3s = 3d'
      WRITE (0,*) '                      3  ?  3s.3p = 4s.3d'
      WRITE (0,*) '                      4  ?  <RETURN>'
      WRITE (0,*) '             Final Term  ?  4S'
      WRITE (0,*) '                      2  ?  RETURN'
      WRITE (0,99002)
99002 FORMAT (/5X,
     &        'Header and Closed Shells cannot exceed 72 characters and'
     &        /5X,'will be copied to the output file. The electrons are'
     &        /5X,'separated by blanks in the Closed Shells.')
      WRITE (0,99003)
99003 FORMAT (5X,'Press RETURN for more... ')
      READ (5,'(A)') STR

 200  WRITE (0,99004)
99004 FORMAT (///////5X,
     &       'Items are separated by a blank in the Reference Set, by a'
     &       /5X,
     &       'comma or a blank in the Active set, and by a period or a'/
     &       5X,'blank in Replacements.'//5X,
     &       'Reference Set, Replacement, and Final Term are three sets'
     &       /5X,
     &       'of input, each with 0 to 10 members.  Each member must be'
     &       /5X,'entered on a separate line.')

      WRITE (0,*) '       PRINT RETURN to terminate the input set.'
      WRITE (0,*) '       PRINT RETURN if the set is empty.'
      WRITE (0,*)
      WRITE (0,*) '     Example 2 :'
      WRITE (0,*) '    --------------------'
      WRITE (0,*)
      WRITE (0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
      WRITE (0,*) '                     2  ?  RETURN'
      WRITE (0,*) '            Active Set  ?  2s,2p,3s'
      WRITE (0,*) 'Type of set generation  ?  0'
      WRITE (0,*) '           Replacement  ?  RETURN '
      WRITE (0,*) '            Final Term  ?  RETURN'
      WRITE (0,*)
      WRITE (0,*)
     &        '     Where the Replacement and the Final Term are empty.'
      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,99010)
      READ (5,'(A)') STR
      I = INDEX(STR,'B')
      J = INDEX(STR,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 100

 300  WRITE (0,99005)
99005 FORMAT (///////5X,
     &      'By inputing "s" or "d" or "sd" you can compute the config-'
     &      /5X,'urations from the Virtual Set, where  S means Single'/5
     &      X,'Replacement, D means Double Replacement, SD means Single'
     &      /5X,'and Double Replacement.'//5X,
     &      'GENCL will give you prompts for the Virtual Set automati-'/
     &      5X,
     &      'cally, then you need to specify the range of shells that'/5
     &      X,'are to participate in the replacements. For instance, a'/
     &      5X,
     &      'response of 2 to "From which shell" and of 3 to "To which'/
     &      5X,'shell" implies that shells 2 and 3 participates in the'/
     &      5X,'replacements, and shell 1 does not enter into any',1X,
     &      'replacements.')
      WRITE (0,*)
      WRITE (0,*) '     Example 3  :'
      WRITE (0,*) '    -------------------------'
      WRITE (0,*) '                ...'
      WRITE (0,*) '         Reference Set  ?  2s(1) 2p(1) 3s(1)'
      WRITE (0,*) '            Active Set  ?  RETURN'
      WRITE (0,*) '           Replacement  ?  sd'
      WRITE (0,*) '           Virtual Set  ?  3p,3d,4s'
      WRITE (0,*) '      From which shell  ?  2'
      WRITE (0,*) '        To which shell  ?  3'
      WRITE (0,*) '            Final Term  ?  RETURN'
      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,99010)
      READ (5,'(A)') STR
      I = INDEX(STR,'B')
      J = INDEX(STR,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 200
 400  WRITE (0,99006)
99006 FORMAT (///////5X,
     &       'After terminating an input line, if you find the previous'
     &       /5X,
     &     'input to be wrong, type ''B'' or ''b'' to go back one step.'
     &     /5X,
     &     'For example, before inputing Active Set, if you find the'/5X
     &     ,'wrong spelling in the Header, type "B" and GENCL will '/5X,
     &     'prompt for the Header again.'//)
      WRITE (0,*) '     Example 4 :'
      WRITE (0,*) '    ----------------------'
      WRITE (0,*) '                Header  ?  OXYYEN'
      WRITE (0,*) '            Active Set  ?  B'
      WRITE (0,*) '                Header  ?  OXYGEN '
      WRITE (0,*) '            Active Set  ?  2s '
      WRITE (0,*)
      WRITE (0,*) '     Then the following prompts will continue.'
      WRITE (0,*)
      WRITE (0,*)
      WRITE (0,99010)
      READ (5,'(A)') STR
      I = INDEX(STR,'B')
      J = INDEX(STR,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 300

 500  WRITE (0,99007)
99007 FORMAT (///////5X,
     &        'Example 4 shows the procedure for going back four steps'/
     &        5X,'to correct the Closed Shells from  5s  to 1s  2s.'//)
      WRITE (0,*) '     Example 5 :'
      WRITE (0,*) '    -----------------------'
      WRITE (0,*)
      WRITE (0,*) '            Active Set  ?  5s'
      WRITE (0,*) 'Type of set generation  ?  0'
      WRITE (0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
      WRITE (0,*) '                     2  ?  2P(4)'
      WRITE (0,*) '                     3  ?  b '
      WRITE (0,*) '                     2  ?  b '
      WRITE (0,*) '         Reference Set  ?  b '
      WRITE (0,*) '         Closed Shells  ?    1s  2s '
      WRITE (0,*) '         Reference Set  ?  '
      WRITE (0,*)
      WRITE (0,*)
!
      WRITE (0,*) '   Then reenter the data for the Reference Set',
     &            ' and continue the input.'
      WRITE (0,*)
      WRITE (0,99010)
      READ (5,'(A)') STR
      I = INDEX(STR,'B')
      J = INDEX(STR,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 400

      WRITE (0,99008)
99008 FORMAT (///////5X,
     &        'When the following error conditions are detected,',
     &        ' a message is given.'/10X,
     &        ' 1). The parentheses are not matched'/10X,
     &       ' 2). The number of electrons in a shell is more than FULL'
     &       //15X,'For each member of the Reference Set ,'/10X,
     &       ' 3). The number of electrons is not the same '/10X,
     &       ' 4). The parity is not the same '//10X,
     &       ' 5). The number of couplings generated by a configuration'
     &       /10X,'     is more than 500'///////)
      WRITE (0,99009)
99009 FORMAT (5X,'Press ''b'' for BACK or RETURN to begin the',1X,
     &        'program.')
      READ (5,'(A)') STR
      I = INDEX(STR,'B')
      J = INDEX(STR,'b')
      IF ( I.NE.0 .OR. J.NE.0 ) GO TO 500
      DO I = 1, 30
         WRITE (0,*)
      END DO
      CONTINUE
99010 FORMAT (5X,'Press ''b'' for BACK or RETURN for more... ')
      END SUBROUTINE HELP
!
! ----------------------------------------------------------------------
!               FUNCTION        I C T O I
! ----------------------------------------------------------------------
!
!       Convert character string into its corresponding integer
!
      INTEGER FUNCTION ICTOI(STR)
      IMPLICIT NONE
!
! Dummy arguments
!
      CHARACTER*(*) STR
!
! Local variables
!
      INTEGER N
!
      N = ICHAR(STR(1:1)) - ICHAR('0')
      IF ( STR(2:2).NE.' ' ) N = N*10 + ICHAR(STR(2:2)) - ICHAR('0')
      ICTOI = N
      CONTINUE
      END FUNCTION ICTOI
!
! ----------------------------------------------------------------------
!               SUBROUTINE      I N P U T
! ----------------------------------------------------------------------
!
!       Process the input set and check the input error
!
      SUBROUTINE INPUT(MAXSET,NSET,SET,MARK,SFLAG,DFLAG,*,*)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER DFLAG, MARK, MAXSET, NSET, SFLAG
      CHARACTER*60 SET(*)
!
! Local variables
!
      CHARACTER*60 CH1, TEMP
      INTEGER I, J
      EXTERNAL DEL
!

!     MAXSET  =  Maximum number of input elements
!     NSET  =  Number of members in the set
!     SET  =  Character array with NSET elements
!     *1  =  Return label if input is 'B'
!     *2  =  Return label if the set is empty
!     MARK  =  1 if input is Replacement, 0 otherwise

      NSET = 0
      SET(1) = '    '
 100  IF ( NSET.LT.MAXSET ) THEN
         READ (5,'(A60)') TEMP

!**************
!        If input is 'B' or 'b', go back one step
!
         I = INDEX(TEMP,'B')
         J = INDEX(TEMP,'b')
         IF ( I.NE.0 .OR. J.NE.0 ) THEN
            IF ( NSET.EQ.0 ) THEN
               RETURN 1
            ELSE
               NSET = NSET - 1
               GO TO 150
            END IF
         END IF

!**************
!        Go to for next input if the input is empty
!        Return if the input is finished
!
         IF ( TEMP(1:5).EQ.'     ' ) THEN
            IF ( NSET.EQ.0 ) THEN
               RETURN 2
            ELSE
               RETURN
            END IF
         END IF
         CALL DEL(TEMP)

!**************
!        If Replacement is 's' or 'd' or 'sd', set single and
!        double flag for Virtual Set
!
         IF ( MARK.NE.0 ) THEN
            CH1 = TEMP(:1)
            IF ( CH1.EQ.'S' .OR. CH1.EQ.'s' ) SFLAG = 1
            I = INDEX(TEMP,'SD')
            J = INDEX(TEMP,'sd')
            IF ( CH1.EQ.'D' .OR. CH1.EQ.'d' .OR. I.NE.0 .OR. J.NE.0 )
     &           DFLAG = 1
            IF ( SFLAG.NE.0 .OR. DFLAG.NE.0 ) THEN
               NSET = 1
               SET(1) = TEMP
               RETURN
            END IF
         END IF

!**************
!        READ the input and delete the repeated member
!
         DO I = 1, NSET
            IF ( SET(I).EQ.TEMP ) THEN
               WRITE (0,*) '     You give a repeated input!'
               GO TO 150
            END IF
         END DO
         NSET = NSET + 1
         SET(NSET) = TEMP
 150     WRITE (0,99001) NSET + 1, '  ?  '
99001    FORMAT (T7,I10,A)
         GO TO 100
      END IF
      CONTINUE
      END SUBROUTINE INPUT
!
! ----------------------------------------------------------------------
!               FUNCTION        L V A L
! ----------------------------------------------------------------------
!
!       convert the symbol into its corresponding quantum number
!
      INTEGER FUNCTION LVAL(SYMBOL)
      IMPLICIT NONE
!
! Dummy arguments
!
      CHARACTER SYMBOL
!
! Local variables
!
      INTEGER LOCATE
      CHARACTER*40 SET
!
      DATA SET/'spdfghiklmnopqrstuvwSPDFGHIKLMNOPQRSTUVW'/

      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE.LE.20 ) THEN
         LVAL = LOCATE - 1
      ELSE
         LVAL = LOCATE - 21
      END IF
      CONTINUE
      END FUNCTION LVAL
!
! ----------------------------------------------------------------------
!               SUBROUTINE      P R I N T
! ----------------------------------------------------------------------
!
!       Print out the values of couplings
!
      SUBROUTINE PRINT_COUPLINGS(EL,Q,M,NC,MARK)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS, NSHEL, NCOUPL, NCFG, NSCOUP
      PARAMETER (NELS=21,NSHEL=8,NCOUPL=2*NSHEL-1,NCFG=10000,
     &           NSCOUP=10000)
!
! COMMON variables
!
      CHARACTER*60 ACT, REF(NELS), REPL(NELS)
      INTEGER CONST, MAXFINAL, MINFINAL, NF, NFTM, NQ, NR, PARITYVAL
      CHARACTER*8 FBETA(7,NSCOUP)
      CHARACTER*32 FILE1(NCFG), FILE3(NSCOUP)
      CHARACTER*40 FILE2(NCFG)
      CHARACTER*2 FTM(NELS)
      CHARACTER*72 HEADER, SHELLS, VIRTUL
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK1  / HEADER, SHELLS, ACT, VIRTUL
      COMMON /BLK2  / REF, REPL, FTM
      COMMON /FILES / FBETA, FILE1, FILE2, FILE3
!
! Dummy arguments
!
      INTEGER M, MARK, NC
      CHARACTER*3 EL(NELS)
      INTEGER Q(NELS)
!
! Local variables
!
      CHARACTER CH1
      CHARACTER*3 COUPLE(NCOUPL)
      INTEGER I, J, K
!
!     MARK  =  1 for Reference Set
!           =  2 for Reference Set, Active Set
!           =  3 for Reference Set, Replacement
!           =  4 for Reference Set, Virtual Set
!
!**************
!     Print the input of Header
!
      J = INDEX(HEADER,'          ')
      WRITE (6,99001) HEADER(1:J)
99001 FORMAT (' '/T10,'-------------       ',A,'  --------',/)

!**************
!     Print the input of Closed Shells
!
      J = INDEX(SHELLS,'     ')
      WRITE (6,99002) SHELLS(:J)
99002 FORMAT ('          Closed Shells  :  ',A/)
      IF ( MARK.EQ.1 .OR. MARK.EQ.2 ) GO TO 200

!**************
!     Print the input of Reference Set
!
      WRITE (6,99003) REF(NF)
99003 FORMAT ('          Reference set  :  ',A60/)
      IF ( MARK.EQ.3 ) GO TO 100

!**************
!     Print the input of Virtual Set
!
      IF ( MARK.EQ.4 ) CH1 = 'S'
      IF ( MARK.EQ.5 ) CH1 = 'D'
      IF ( MARK.EQ.4 .OR. MARK.EQ.5 ) THEN
         WRITE (6,99004) VIRTUL
99004    FORMAT ('            Virtual Set  :  ',A60/)
         WRITE (6,99005) CH1, REPL(NELS)
99005    FORMAT ('          ',A1,'-Replacement  :  ',A60/)
         GO TO 200
      END IF

!**************
!     Print the input of Active Set
!
      WRITE (6,99006) ACT
99006 FORMAT ('             Active Set  :  ',A60/)
      IF ( MARK.EQ.2 ) GO TO 200

!**************
!     Print the input of Replacement
!
 100  WRITE (6,99007) REPL(NR)
99007 FORMAT ('            Replacement  :  ',A60/)

!**************
!     Print the new configuration by Replacement
!
 200  IF ( EL(1)(1:1).EQ.' ' ) THEN
         K = 1
      ELSE
         K = 2
      END IF
      WRITE (6,99008) (EL(J),Q(J),J=1,M)
99008 FORMAT ('          Configuration  :',5(2X,A3,'(',I2,')')/)

!**************
!     Print couplings generated from the configuration
!
      K = 2*M - 1
      IF ( M.LE.3 ) THEN
         READ (FILE3(1),99013) (COUPLE(J),J=1,K)
         WRITE (6,99009) (COUPLE(J),J=1,K)
99009    FORMAT ('        Their couplings  :  ',9(5X,A3))
         DO I = 2, NC
            READ (FILE3(I),99013) (COUPLE(J),J=1,K)
            WRITE (6,99010) (COUPLE(J),J=1,K)
99010       FORMAT (T29,9(5X,A3))
         END DO
      ELSE
         WRITE (6,99011)
99011    FORMAT ('         Their couplings  :'/)
         DO I = 1, NC
            READ (FILE3(I),99013) (COUPLE(J),J=1,K)
            WRITE (6,99012) (COUPLE(J),J=1,K)
99012       FORMAT (T5,9(5X,A3))
         END DO
      END IF
      CONTINUE
99013 FORMAT (9(A3))
      END SUBROUTINE PRINT_COUPLINGS
!
! ----------------------------------------------------------------------
!       SUBROUTINE      R E P L A C E
! ----------------------------------------------------------------------
!
      SUBROUTINE REPLAC(ELB,QB,MB,*)

!   OUTPUT :
!       ELB,QB,MB
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS, NCFG, NSCOUP
      PARAMETER (NELS=21,NCFG=10000,NSCOUP=10000)
!
! COMMON variables
!
      INTEGER CONST, M, MA, MAXFINAL, MINFINAL, ML, MR, MS(NELS), NF,
     &        NFTM, NQ, NR, NREF, PARITYVAL, Q(NELS), Q1, Q10, Q11, Q12,
     &        Q13, Q14, Q15, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, QL(NELS),
     &        QR(NELS), QS(NELS,NELS), RL(NELS)
      CHARACTER*3 EL(NELS), ELA(NELS), ELL(NELS), ELR(NELS), ELS(NELS)
      CHARACTER*8 FBETA(7,NSCOUP)
      CHARACTER*32 FILE1(NCFG), FILE3(NSCOUP)
      CHARACTER*40 FILE2(NCFG)
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK3  / EL, ELL, ELR, ELS, ELA
      COMMON /BLK4  / Q, QL, ML, QR, MR, M, QS, MS, MA, RL, NREF, Q1,
     &                Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12,
     &                Q13, Q14, Q15
      COMMON /FILES / FBETA, FILE1, FILE2, FILE3
!
! Dummy arguments
!
      INTEGER MB
      CHARACTER*3 ELB(NELS)
      INTEGER QB(NELS)
!
! Local variables
!
      CHARACTER*3 CH3, ELC(NELS)
      INTEGER I, J, K, L, MARK, MC, N, NP, PP(NCFG), QA(NELS), QC(NELS),
     &        QQ
      INTEGER LVAL
      EXTERNAL LVAL
!


      DO I = 1, NELS
         ELC(I) = EL(I)
         QC(I) = Q(I)
      END DO

!
!*********     Correct Q(i) by subtraction
!
      DO I = 1, ML
         CH3 = ELL(I)
         MARK = 1
         DO J = 1, M
            IF ( CH3.EQ.ELC(J) ) THEN
               QC(J) = QC(J) - QL(I)
               IF ( QC(J).GE.0 ) MARK = 0
            END IF
         END DO
      END DO
      IF ( MARK.NE.0 ) RETURN 1

!
!*********     Correct QC(i) by adding
!
      MC = M
      DO I = 1, MR
         CH3 = ELR(I)
         MARK = 0
         DO J = 1, MC
            IF ( CH3.EQ.ELC(J) ) THEN
               QC(J) = QC(J) + QR(I)
               MARK = 1
            END IF
         END DO
         IF ( MARK.EQ.0 ) THEN
            MC = MC + 1
            ELC(MC) = CH3
            QC(MC) = QR(I)
         END IF
      END DO

!**************
!     Delete EL(i) if Q(i) = 0
!
      MB = 0
      DO I = 1, MC
         IF ( QC(I).NE.0 .OR. MC.EQ.2 ) THEN
            MB = MB + 1
            ELB(MB) = ELC(I)
            QB(MB) = QC(I)
         END IF
      END DO

!
!*********     Check the input error after replacement
!
      J = 0
      K = 0
      DO I = 1, MB
         CH3 = ELB(I)
         L = LVAL(CH3(2:2))
         IF ( QB(I).GT.L*4+2 ) RETURN 1
         J = J + QB(I)
         K = K + QB(I)*L
      END DO

      IF ( J.NE.CONST ) RETURN 1

      IF ( MOD(K,2).NE.PARITYVAL ) RETURN 1

!
!*********    If the replacement duplicates a configuration in the
!     active set, it should not be sent to CI.LST
!
      DO I = 1, NQ
         READ (FILE1(I),99001) (QA(J),J=1,MA)
99001    FORMAT (15(I2))
         N = 0
         DO J = 1, MB
            CH3 = ELB(J)
            DO K = 1, MA
               IF ( CH3.EQ.ELA(K) .AND. QB(K).EQ.QA(K) ) N = N + 1
            END DO
         END DO
         IF ( N.EQ.MB ) RETURN 1
      END DO

!**************
!     If the replacement duplicates a configuration in the
!     Reference Set, it should not be sent to CI.LST
!
      DO I = 1, NREF
         L = MS(I)
         N = 0
         DO J = 1, L
            CH3 = ELS(J)
            QQ = QS(J,I)
            DO K = 1, MB
               IF ( CH3.EQ.ELB(K) .AND. QQ.EQ.QB(K) ) N = N + 1
            END DO
         END DO
         IF ( N.EQ.MB ) RETURN 1
      END DO

!**************
!     If the replacement duplicates a configuration in the previous
!     replacement, it should not sent to CI.LST
!
      DO I = 1, NP
         L = PP(I)
         READ (FILE2(I),99002) (ELC(J),QC(J),J=1,L)
99002    FORMAT (8(A3,I2))
         N = 0
         DO J = 1, L
            CH3 = ELC(J)
            DO K = 1, MB
               IF ( CH3.EQ.ELB(K) .AND. QC(K).EQ.QB(K) ) N = N + 1
            END DO
         END DO
         IF ( N.EQ.MB ) RETURN 1
      END DO

      NP = NP + 1
      PP(NP) = MB
      WRITE (FILE2(NP),99003) (ELB(J),QB(J),J=1,MB)
99003 FORMAT (8(A3,I2))
      CONTINUE
      END SUBROUTINE REPLAC
!
! ----------------------------------------------------------------------
!               SUBROUTINE        S T R S H
! ----------------------------------------------------------------------
!
!       Shift the string left
!
      SUBROUTINE STRSH(STR,I)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER I
      CHARACTER*(*) STR
!
! Local variables
!
      CHARACTER*72 TEMP
!
      TEMP = STR(I:)
      STR = TEMP
      CONTINUE
      END SUBROUTINE STRSH
!
! ----------------------------------------------------------------------
!               FUNCTION        S Y M B
! ----------------------------------------------------------------------
!
!       Convert the quantum number into its corresponding symbol
!
      CHARACTER FUNCTION SYMB(L)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER L
!
! Local variables
!
      CHARACTER*20 SET
!
      DATA SET/'SPDFGHIKLMNOPQRSTUVW'/

      SYMB = SET(L+1:L+1)
      CONTINUE
      END FUNCTION SYMB
!
! ----------------------------------------------------------------------
!       SUBROUTINE      V P A I R
! ----------------------------------------------------------------------
!
!       Generate occupied or virtual pairs for D-Replacement
!
      SUBROUTINE VPAIR(ELV,MV,PL,LLMIN,LLMAX,STR,*)

!   INPUT :
!       ELV  =  ELi for Virtual Set
!        MV  =  Number of ELi for Virtual Set
!        PL  =  Parity of Qi for Reference Set
!     LLMIN  =  Minimum angular coupling of the pair
!     LLMAX  =  Maximum angular coupling of the pair
!       STR  =  String to be packed as output for Replacement
!

      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER NELS
      PARAMETER (NELS=21)
!
! COMMON variables
!
      INTEGER CONST, M, MA, MAXFINAL, MINFINAL, ML, MR, MS(NELS), NF,
     &        NFTM, NQ, NR, NREF, PARITYVAL, Q(NELS), Q1, Q10, Q11, Q12,
     &        Q13, Q14, Q15, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, QL(NELS),
     &        QR(NELS), QS(NELS,NELS), RL(NELS)
      CHARACTER*3 EL(NELS), ELA(NELS), ELL(NELS), ELR(NELS), ELS(NELS)
      CHARACTER*2 FTM(NELS)
      CHARACTER*60 REF(NELS), REPL(NELS)
      COMMON  NF, NR, NFTM, MAXFINAL, MINFINAL, PARITYVAL, CONST, NQ
      COMMON /BLK2  / REF, REPL, FTM
      COMMON /BLK3  / EL, ELL, ELR, ELS, ELA
      COMMON /BLK4  / Q, QL, ML, QR, MR, M, QS, MS, MA, RL, NREF, Q1,
     &                Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12,
     &                Q13, Q14, Q15
!
! Dummy arguments
!
      INTEGER LLMAX, LLMIN, MV, PL
      CHARACTER*72 STR
      CHARACTER*3 ELV(NELS)
!
! Local variables
!
      CHARACTER*3 ELB(NELS), ELC(NELS)
      INTEGER I, J, K, LR1, LR2, LRMAX, LRMIN, N, NC, PR, QB(NELS)
      INTEGER LVAL
      CHARACTER*72 SSTR
      EXTERNAL COUPLD, DEL, LVAL, PRINT_COUPLINGS, REPLAC
!

      N = INDEX(STR,'=')
      NR = NELS

!**************
!     D-Replacement for the pair of two single ELi
!
      MR = 2
      QR(1) = 1
      QR(2) = 1
      DO I = 1, MV - 1
         ELR(1) = ELV(I)
         LR1 = LVAL(ELR(1)(2:2))
         DO J = I + 1, MV
            ELR(2) = ELV(J)
            LR2 = LVAL(ELR(2)(2:2))
            LRMIN = IABS(LR1-LR2)
            LRMAX = LR1 + LR2
            PR = MOD(LRMAX,2)
!
!           If the pair has the same parity with the left side, and the
!           angular coupling of the two pairs have values in common,
!           replace them, then generate couplings for the new
!           configuration
            IF ( PL.NE.PR .OR. LRMIN.GT.LLMAX .OR. LRMAX.LT.LLMIN )
     &           CYCLE
            CALL REPLAC(ELB,QB,K,*50)
            CALL COUPLD(ELB,QB,K,ELC,NC,*200)
            IF ( NC.GT.0 ) THEN
               WRITE (6,99001)
99001          FORMAT (//'   FOR VIRTUAL SET, GENERATE',
     &                 ' CONFIGURATION AND COUPLINGS FOR D-REPLACEMENT'
     &                 //)
               SSTR = STR(:N)//ELR(1)//'.'//ELR(2)
               CALL DEL(SSTR)
               REPL(NR) = SSTR
               CALL PRINT_COUPLINGS(ELC,QB,K,NC,4)
            END IF
 50      END DO
      END DO

!**************
!     D-Replacement for the pairs which has the value Qi=2
!
      MR = 1
      QR(1) = 2
      DO I = 1, MV
         ELR(1) = ELV(I)
         LR1 = LVAL(ELR(1)(2:2))
         LRMAX = LR1 + LR1
         PR = MOD(LRMAX,2)
!
!        If it has the same parity with the left side, replace them,
!        then generate couplings for the new configuration
!
         IF ( PL.EQ.PR .AND. LRMAX.GE.LLMIN ) THEN
            CALL REPLAC(ELB,QB,K,*100)
            CALL COUPLD(ELB,QB,K,ELC,NC,*200)
            IF ( NC.GT.0 ) THEN
               WRITE (6,99002)
99002          FORMAT (//'   FOR VIRTUAL SET, GENERATE ',
     &                 ' CONFIGURATION AND COUPLINGS FOR D-REPLACEMENT'
     &                 //)
               SSTR = STR(:N)//ELR(1)//'(2)'
               CALL DEL(SSTR)
               REPL(NR) = SSTR
               CALL PRINT_COUPLINGS(ELC,QB,K,NC,4)
            END IF
         END IF
 100  END DO

      RETURN
 200  RETURN 1
      END SUBROUTINE VPAIR
