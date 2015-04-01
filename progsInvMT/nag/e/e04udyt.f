      SUBROUTINE E04UDY(NDICT,DICTRY,ALPHA,KEY,ENTRY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-724 (DEC 1989).
C
C     Description and usage:
C
C       Performs dictionary lookups.  A pointer is returned if a
C     match is found between the input key and the corresponding
C     initial characters of one of the elements of the dictionary.
C     If a 'synonym' has been provided for an entry, the search is
C     continued until a match to a primary dictionary entry is found.
C     Cases of no match, or multiple matches, are also provided for.
C
C     Dictionary entries must be left-justified, and may be alphabetized
C     for faster searches.  Secondary entries, if any, are composed of
C     two words separated by one or more characters such as blank, tab,
C     comma, colon, or equal sign which are treated as non-significant
C     by E04UDW.  The first entry of each such pair serves as a synonym
C     for the second, more fundamental keyword.
C
C       The ordered search stops after the section of the dictionary
C     having the same first letters as the key has been checked, or
C     after a specified number of entries have been examined.  A special
C     dictionary entry, the currency symbol '$', will also terminate the
C     search.  This will speed things up if an appropriate dictionary
C     length parameter cannot be determined.  Both types of search are
C     sequential.  See 'Notes' below for some suggestions if efficiency
C     is an issue.
C
C
C     Parameters:
C
C     Name    Dimension  Type  I/O/S  Description
C     NDICT               I    I     Number of dictionary entries to be
C                                    examined.
C     DICTRY  NDICT       C    I     Array of dictionary entries,
C                                    left-justified in their fields.
C                                    May be alphabetized for efficiency,
C                                    in which case ALPHA should be
C                                    .TRUE.  Entries with synonyms are
C                                    of the form
C                                    'ENTRY : SYNONYM', where 'SYNONYM'
C                                    is a more fundamental entry in the
C                                    same dictionary.  NOTE: Don't build
C                                    'circular' dictionaries.
C     ALPHA               L    I     Indicates whether the dictionary
C                                    is in alphabetical order, in which
C                                    case the search can be terminated
C                                    sooner.
C     KEY                 C    I/O   String to be compared against the
C                                    dictionary.  Abbreviations are OK
C                                    if they correspond to a unique
C                                    entry in the dictionary.  KEY is
C                                    replaced on termination by its most
C                                    fundamental equivalent dictionary
C                                    entry (uppercase, left-justified)
C                                    if a match was found.
C     ENTRY               I      O   Dictionary pointer.  If .GT. 0, it
C                                    indicates which entry matched KEY.
C                                    In case of trouble, a negative
C                                    value means that a UNIQUE match
C                                    was not found - the absolute value
C                                    of ENTRY points to the second
C                                    dictionary entry that matched KEY.
C                                    Zero means that NO match could be
C                                    found.  ENTRY always refers to the
C                                    last search performed -
C                                    in searching a chain of synonyms,
C                                    a non-positive value will be
C                                    returned if there is any break,
C                                    even if the original input key
C                                    was found.
C
C
C     External references:
C
C     Name    Description
C     E04UDW  Finds first and last significant characters.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C     Notes:
C
C     (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C     (2)  We have assumed that the dictionary is not too big.  If
C         many searches are to be done or if the dictionary has more
C         than a dozen or so entries, it may be advantageous to build
C         an index array of pointers to the beginning of the section
C         of the dictionary containing each letter, then pass in the
C         portion of the dictionary beginning with DICTRY (INDEX).
C         (This won't generally work for dictionaries with synonyms.)
C         For very large problems, a completely different approach may
C         be advisable, e.g. a binary search for ordered dictionaries.
C
C     (3)  E04UDY is case sensitive.  In most applications it will be
C         necessary to use an uppercase dictionary, and to convert the
C         input key to uppercase before calling E04UDY.  Companion
C         routines E04UDV and PAIRS, available from the author, already
C         take care of this.
C
C     (4)  The key need not be left-justified.  Any leading (or
C         trailing) characters which are 'non-significant' to E04UDW
C         will be ignored.  These include blanks, horizontal tabs,
C         commas, colons, and equal signs.  See E04UDW for details.
C
C     (5)  The ASCII collating sequence for character data is assumed.
C         (N.B. This means the numerals precede the alphabet, unlike
C         common practice.)  This should not cause trouble on EBCDIC
C         machines if DICTRY just contains alphabetic keywords.
C         Otherwise it may be necessary to use the FORTRAN lexical
C         library routines to force use of the ASCII sequence.
C
C     (6)  Parameter NUMSIG sets a limit on the length of significant
C         dictionary entries.  Special applications may require that
C         this be increased.  (It is 16 in the present version.)
C
C     (7)  No protection against 'circular' dictionaries is provided:
C         don't claim that A is B, and that B is A.  All synonym chains
C         must terminate.  Other potential errors not checked for
C         include duplicate or mis-ordered entries.
C
C     (8)  The handling of ambiguities introduces some ambiguity:
C
C            ALPHA = .TRUE.  A potential problem, when one entry
C                            looks like an abbreviation for another
C                            (eg. does 'A' match 'A' or 'AB') was
C                            resolved by dropping out of the search
C                            immediately when an 'exact' match is found.
C
C            ALPHA = .FALSE. The programmer must ensure that the above
C                            situation does not arise: each dictionary
C                            entry must be recognizable, at least when
C                            specified to full length.  Otherwise, the
C                            result of a search will depend on the
C                            order of entries.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C     24 Feb. 1984  RAK/DAS  Initial design and coding.
C     25 Feb. 1984    RAK    Combined the two searches by suitable
C                            choice of terminator FLAG.
C     28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                            longer update KEY.
C     29 Mar. 1984    RAK    Put back replacement of KEY by its
C                            corresponding entry.
C     21 June 1984    RAK    Corrected bug in error handling for cases
C                            where no match was found.
C     23 Apr. 1985    RAK    Introduced test for exact matches, which
C                            permits use of dictionary entries which
C                            would appear to be ambiguous (for ordered
C                            case).  Return -I to point to the entry
C                            which appeared ambiguous (had been -1).
C                            Repaired loop termination - had to use
C                            equal length strings or risk quitting too
C                            soon when one entry is an abbreviation
C                            for another.  Eliminated HIT, reduced
C                            NUMSIG to 16.
C     28 May 1986     MPH    Changed to test if an apparently
C                            ambiguous keyword is in fact unique.
C                            If ALPHA is .FALSE. FLAG is now '$'.
C                            Local CHARACTER scalars declared as
C                            CHARACTER*16.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER         BLANK, CURLY
      INTEGER           NUMSIG
      PARAMETER         (BLANK=' ',CURLY='$',NUMSIG=16)
C     .. Scalar Arguments ..
      INTEGER           ENTRY, NDICT
      LOGICAL           ALPHA
      CHARACTER*(*)     KEY
C     .. Array Arguments ..
      CHARACTER*(*)     DICTRY(NDICT)
C     .. Local Scalars ..
      INTEGER           FIRST, I, IFRST, IFRST1, ILAST, ILAST1, ILEN,
     *                  ILEN1, ILST, IMARK, IMARK1, LAST, LENGTH, MARK
      CHARACTER*16      FLAG, TARGET, TRGT, TRGT1
C     .. External Subroutines ..
      EXTERNAL          E04UDW
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, LEN, LLE
C     .. Executable Statements ..
C
      ENTRY = 0
C
C     Isolate the significant portion of the input key (if any).
C
      FIRST = 1
      LAST = MIN(LEN(KEY),NUMSIG)
      CALL E04UDW(KEY,FIRST,LAST,MARK)
C
      IF (MARK.GT.0) THEN
         TARGET = KEY(FIRST:MARK)
C
C        Look up TARGET in the dictionary.
C
   20    CONTINUE
         LENGTH = MARK - FIRST + 1
C
C           Select search strategy by cunning choice of termination test
C           flag.  The left curly bracket follows all the alphabetic
C           characters in the ASCII collating sequence, but precedes the
C           vertical bar.
C
         IF (ALPHA) THEN
            FLAG = TARGET
         ELSE
            FLAG = CURLY
         END IF
C
C
C           Perform search.
C           ---------------
C
         I = 0
   40    CONTINUE
         I = I + 1
         IF (TARGET(1:LENGTH).EQ.DICTRY(I)(1:LENGTH)) THEN
            IF (ENTRY.EQ.0) THEN
C
C                    First 'hit' - must still guard against ambiguities
C                    by searching until we've gone beyond the key
C                    (ordered dictionary) or until the end-of-dictionary
C                    mark is reached (exhaustive search).
C
               ENTRY = I
C
C                    Special handling if match is exact - terminate
C                    search.  We thus avoid confusion if one dictionary
C                    entry looks like an abbreviation of another.
C                    This fix won't generally work for un-ordered
C                    dictionaries.
C
               FIRST = 1
               LAST = NUMSIG
               CALL E04UDW(DICTRY(ENTRY),FIRST,LAST,MARK)
               IF (MARK.EQ.LENGTH) I = NDICT
            ELSE
C                    If two hits check if they are attempting to
C                    indicate the same dictionary entry.
C
C                    Extract keyword from first match found
C
               ILST = NUMSIG
               IFRST = MARK + 2
               CALL E04UDW(DICTRY(ENTRY),IFRST,ILST,IMARK)
               IF (IMARK.GT.0) THEN
                  TRGT = DICTRY(ENTRY) (IFRST:IMARK)
                  ILEN = IMARK - IFRST + 1
               ELSE
                  TRGT = DICTRY(ENTRY) (FIRST:MARK)
                  ILEN = MARK - FIRST + 1
               END IF
C
C                    Extract keyword from next match found
C
               IFRST = 1
               ILAST = NUMSIG
               CALL E04UDW(DICTRY(I),IFRST,ILAST,IMARK)
               ILAST1 = NUMSIG
               IFRST1 = IMARK + 2
               CALL E04UDW(DICTRY(I),IFRST1,ILAST1,IMARK1)
               IF (IMARK1.GT.0) THEN
                  TRGT1 = DICTRY(I) (IFRST1:IMARK1)
                  ILEN1 = IMARK1 - IFRST1 + 1
               ELSE
                  TRGT1 = DICTRY(I) (IFRST:IMARK)
                  ILEN1 = IMARK - IFRST + 1
               END IF
C
C                    If keywords not identical then ambiguity
C
               IF (TRGT(1:ILEN).NE.TRGT1(1:ILEN1)) THEN
C
C                       Oops - two hits.  Abnormal termination.
C                       ---------------------------------------
C
                  ENTRY = -I
                  RETURN
               END IF
            END IF
         END IF
C
C           Check whether we've gone past the appropriate section of the
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.
C
         IF (LLE(DICTRY(I)(1:LENGTH),FLAG) .AND. I.LT.NDICT) GO TO 40
C
C
C           Check for a synonym.
C           --------------------
C
         IF (ENTRY.GT.0) THEN
C
C              Look for a second entry 'behind' the first entry.  FIRST
C              and MARK were determined above when the hit was detected.
C
            FIRST = MARK + 2
            CALL E04UDW(DICTRY(ENTRY),FIRST,LAST,MARK)
            IF (MARK.GT.0) THEN
C
C                 Re-set target and dictionary pointer, then repeat the
C                 search for the synonym instead of the original key.
C
               TARGET = DICTRY(ENTRY) (FIRST:MARK)
               ENTRY = 0
               GO TO 20
C
            END IF
         END IF
C
      END IF
      IF (ENTRY.GT.0) KEY = DICTRY(ENTRY)
C
C
C     Normal termination.
C     -------------------
C
      RETURN
C
C     End of E04UDY.  (CMLOOK/OPLOOK)
      END
