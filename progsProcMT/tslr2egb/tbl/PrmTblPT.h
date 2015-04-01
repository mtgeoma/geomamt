/* PrmTblPT.h                          1:20 PM 02/04/97                       */
/*                                                                            */
/* Defines the parameter table which controls and defines status of all       */
/* functions of the instrument.                                               */
/*                                                                            */
/* 1996-Dec-31   DJD   Create.                                                */
/* 1997-Aug-07   GB    Modified:        include gps.h                         */
/* 1997-Sep-11   GB    Added:    SNUM, SITE, FILE, JCLK, OCTR, XDOS           */
/* 1997-Oct-01   GB    Changed:  L3NP, L4NP to L3NS, L4NS                     */
/* 1997-Oct-10   GB    Added: TOTL, BADR, SGIN                                */
/* 1997-Oct-23   MJM   Added: MODE                                            */
/* 1997-Oct-29   MJM   Added: BAT1, BAT2, BAT3, TEMP                          */
/* 1998-Feb-02   MJM   Moved parameter table definition here from "gps.h".    */
/* 1998-May-7    MJM   Reduce length of P field in ValPT to 13 bytes.         */


#ifndef PrmTblPT                       /* Skip if this file already included. */
#define PrmTblPT        

#include "gps.h"
#include "amx831sd.h"

typedef enum   {                       /* Defines data type of entries in     */
                                       /* parameter table. Corresponds to     */
                                       /* fields of type ValPT.               */    
                  IntPT,               /* Signed 4 byte integer.              */
                  FltPT,               /* 8 byte IEEE floating point (double) */
                  StrPT,               /* String, null terminated, 0-8 bytes. */
                  UTCPT,               /* UTC date and time, 1 s resolution.  */
                  PosPT,               /* Geographical position. (char *)     */
                  AmxPT                /* AMX time and date.                  */
} TypePT;

typedef union {                        /* Defines a value in the parameter    */
                                       /* table. Corresponds to elements of   */
                                       /* type TypePT.                        */
      signed long       I;             /* Signed 4 byte integer.              */
      double            F;             /* 8 byte IEEE floating point.         */
      char              S[9];          /* String, null terminated, 0-8 bytes. */
      char              P[13];         /* Geographical position,5m resolution.*/
      struct amxtds     TD;            /* AMX timedate.                       */
} ValPT;   


typedef        struct {
      char              Code[5];       /* Ascii code for the parameter.       */
                                       /* Up to 4 characters, null terminated */
                                       /* if less than 4. Case ignored.       */
      unsigned short    Grp;           /* Group of parameters of which this   */
                                       /* parameter is a member. One semaphore*/
                                       /* is shared by all parameters in a    */
                                       /* group.                              */
      long int        Smph;            /* ID of the semaphore which protects  */
                                       /* this parameter.                     */
      TypePT            T;             /* The data type of this parameter.    */
      ValPT             V;             /* The value of this parameter.        */
} TableEntryPT; 

int  cdecl ReadPT( const char *code, TableEntryPT *ptp );
int  cdecl WritePT( const char *code, TableEntryPT *ptp );
int  cdecl SavePT( void );
void cdecl GetPT( void );
void cdecl PrmTblRP( void );

#define PrmFile         "STARTUP.TBL"       /* Startup Parameter File   */
#define BlnkFltPT       0
#define BlnkIntPT       0L
#define BlnkPosPT       0L
#define BlnkStrPT       0L
#define ETX             '\3'         

                
      
#else
#endif

