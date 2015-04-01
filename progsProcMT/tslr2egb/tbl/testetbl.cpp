#include <stdio.h>

typedef struct
{                                      /* AMX Time/Date Structure             */
      unsigned char amtdsec;           /* seconds  (0-59)                     */
      unsigned char amtdmin;           /* minutes  (0-59)                     */
      unsigned char amtdhr;            /* hours    (0-23)                     */
      unsigned char amtddy;            /* day      (1-31)                     */
      unsigned char amtdmn;            /* month    (1-12)                     */
      unsigned char amtdyr;            /* year     (0-99)                     */
      unsigned char amtdow;            /* day of week (Mon=1 to Sun=7)        */
      unsigned char amtdcen;           /* 0 if time/date is incorrect         */
                                       /* century if time/date is correct     */
} amxtds;

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
      amxtds            TD;            /* AMX timedate.                       */
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

int main()
{
  printf("sizeof(amxtds)      =%d\n",sizeof(amxtds));
  printf("sizeof(TypePT)      =%d\n",sizeof(TypePT));
  printf("sizeof(ValPT)       =%d\n",sizeof(ValPT));
  printf("sizeof(TableEntryPT)=%d\n",sizeof(TableEntryPT));

  return 0;
}
