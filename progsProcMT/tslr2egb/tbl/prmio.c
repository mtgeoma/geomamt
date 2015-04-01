/* prmio.c - functions for Parameter Table i/o                                */

/* Mar-11-1997  GB      create                                                */
/* Sep-17-1997  GB      added SavePT()                                        */
/* Sep-25-1997  GB      added SemsPT[]  to resolve the din changing &sem1id   */
/* Jan-12-1998  MJM     Close TblFil in GetPT() if opened.                    */
/*                      Stick dummy FILE name into parameter table in GetPT().*/
/* Feb-04-1998  MJM     Overwrite dummy filename if "startup.tbl" is loaded.  */
/* April-8-1998 MJM     Add RQST parameter.                                   */
/* May-1-1998   MJM     Add LATG and LNGG parameters. Delete PGPS parameter.  */
/* May 19, 1998 MJM     Add VER parameter.                                    */
/*                      Disable interrupts before reserving resource semaphore*/
/*                      to prevent a lower priority task from holding up      */
/*                      execution of a higher priority task that wants to     */
/*                      reserve the same semaphore.                           */
/* June 17, 1998 MJM    Re-write GetPT() to read startup.tbl as individual    */
/*                      parameters instead of as a block. Check for write     */
/*                      protection and write using WritePT().                 */
/* June 29, 1998 MJM    Modify SavePT() so that file not written if FILE      */
/*                      parameter is empty.                                   */
/* July 8, 1998. MJM    Set default HGN to 3 instead of 1.                    */
/* July 14, 1998. MJM   Do not load "startup.tbl" if this is an MTU-CLK.      */
/* August 11, 1998. MJM Add "TXPR" to list of writeable parameters.           */
/*                      Add call to CreateFilename() to prevent overwriting   */
/*                      .tbl file specified in startup.tbl if box is powered  */
/*                      up again after performing acquisition.                */
/* August 28, 1998. MJM Add 9 new writeable parameters.                       */
/* September 3, 1998. MJM Add 9 more parameters from document.                */
/* September 8, 1998. MJM Check for overwrite of existing .tbl file in SavePT.*/
/* September 15, 1998. MJM Put default parameters in for coils.               */
/* September 17, 1998. MJM Change coil names to type StrPT from PosPT.        */
/* October 9, 1998. MJM  Change HAMP value to -0.206 from 0.206 for coil phase*/
/* October 13, 1998. MJM Change E line lengths and azimuths to floating point.*/
/* October 29, 1998. MJM Remove reset of XDOS and CLST before saving parameter*/
/*                       table in SavePT().                                   */
/* November 5, 1998. MJM Add check of "mtucl4" in GetPT().                    */
/* December 9, 1998. MJM Write default coil names into table before reading   */
/*                   "Startup.tbl" to avoid overwriting them.                 */
/*                   Clear error in GetPT() after bad parameter type read.    */
/*                   Call fclose() to close "Startup.tbl" if opened.          */
/* January 29, 1999. MJM Add EXAC, EXDC, EYAC, EYDC to parameter table.       */
/* February 1, 1999. MJM Use AMX semaphore sem0id to protect above parameters */
/*                       and CHEX, CHEY, CHHX, CHHY, CHHZ.                    */
/* February 2, 1999. MJM Use new semaphore sem4id to protect CHnn.            */
/* February 10, 1999.MJM Change TCMB value to 2 from 1 to indicate new type   */
/*                       of comb filter used in DSP code.                     */
/* February 22, 1999. GB modify param list for LRMT                           */
/* April 20, 1999. MJM  Remove use of resource semaphore in ReadPT() and      */
/*                      WritePT() as AMX function re-enables interrupts       */
/*                      permitting lower priority task to hold up a higher    */
/*                      priority one leading to the code 6 (timeout) error in */
/*                      the DMA data transfer as proven by the logic analyzer.*/
/* Nov. 18, 1999. GB    Move xdos in PrmTbl to 4-th place.                    */
/* Nov. 24, 1999. GB    Add NOBU parameter for counting not available buffers */
/* Dec. 08, 1999. GB    Initialize CFMN, CFMX to zero for LR box calibration  */


#include <stdio.h>
#include <string.h>
#include <dos.h>
#include "amx831sd.h"
#include "amx831ec.h"
#include "gps.h"
#include "PrmTblPT.h"
#include "misc.h"
#include "Mode.h"
#include "HseKeep.h"
#include "CrFiName.h"
#include "CalTask.h"

extern AMXID sem0id;
extern AMXID sem1id;                      /* Semaphore ids                    */
extern AMXID sem2id;
extern AMXID sem3id;
extern AMXID sem4id;
extern AMXID GetPTid;

AMXID   * SemsPT[] = {&sem0id, &sem1id, &sem2id, &sem3id, &sem4id};

TableEntryPT      PrmTab[] = {            /* The parameter table.             */

     { "RQST", 3, 3L, IntPT, SETUP     }, /* Request a change in MODE.        */
     { "MODE", 3, 3L, IntPT, 1L        }, /* MTU Mode.Indicates what state the*/
                                          /* MTU is running in.               */
     { "AQST", 3, 3L, IntPT, 0L        }, /* Acquisition status. 0 = Waiting, */
                                          /* 1 = In progress, 2 = Completed.  */
     { "XDOS", 3, 3L, IntPT, BlnkIntPT }, /* Exit to DOS                      */
                                          /* = 0 if don't exit to DOS.        */
                                          /* = 1 if quit mtup and exit to DOS.*/
                                          /* = 4 if quit & pwrdown after SCNT.*/
     { "SNUM", 1, 1L, IntPT, 0L        }, /* Serial Number of MTU unit.       */
     { "VER",  1, 1L, StrPT, 0L        }, /* Software version string.         */
     { "HW",   1, 1L, StrPT, BlnkStrPT }, /* Hardware type.                   */
     { "SITE", 1, 1L, StrPT, BlnkStrPT }, /* Site ID of the measurement.      */
     { "CMPY", 1, 1L, PosPT, BlnkStrPT }, /* Company name.                    */
     { "SRVY", 1, 1L, PosPT, BlnkStrPT }, /* Survey ID.                       */
     { "FILE", 1, 1L, StrPT, BlnkStrPT }, /* File Name for data logging.      */
     { "NUTC", 2, 2L, AmxPT, 0L        }, /* AMX timedate.                    */   
     { "AWIN", 1, 1L, IntPT, 1L        }, /* Current data acquisition window. */
     { "STM1", 1, 1L, AmxPT, 0L        }, /* Window1 start time.              */
     { "ETM1", 1, 1L, AmxPT, 0L        }, /* Window1 end time.                */
     { "STM2", 1, 1L, AmxPT, 0L        }, /* Window2 start time.              */
     { "ETM2", 1, 1L, AmxPT, 0L        }, /* Window2 end time.                */
     { "STM3", 1, 1L, AmxPT, 0L        }, /* Window3 start time.              */
     { "ETM3", 1, 1L, AmxPT, 0L        }, /* Window3 end time.                */
     { "SRW1", 1, 1L, IntPT, 0L        }, /* Sample rate for window1          */
     { "SRW2", 1, 1L, IntPT, 3L        }, /* Sample rate for window2          */
     { "SRW3", 1, 1L, IntPT, 3L        }, /* Sample rate for window3          */
     { "SGIN", 1, 1L, IntPT, BlnkIntPT }, /* Signal Input 0=extrnl 1=intr TSIG*/
                                          /* 2=external signal,but test signal*/
                                          /* enabled for 'H' box coil cal.    */
     { "EGN",  1, 1L, IntPT, 3L        }, /* Gain for "E" channels. 3,12,48   */
     { "HGN",  1, 1L, IntPT, 2L        }, /* Gain for "H" channels. 2, 6, 8   */
     { "LFRQ", 1, 1L, IntPT, 60L       }, /* Line frequency, Hz.              */
     
     { "SRPM", 1, 1L, IntPT, 1440L     }, /* Number of samples per minute     */
     { "STIM", 1, 1L, AmxPT, 0L        }, /* Data acquisition start time.     */
     { "ETIM", 1, 1L, AmxPT, 0L        }, /* Data acquisition end time.       */
                                          
     { "CHEX", 4, 4L, IntPT, 0L        }, /* Front end channel for EX.        */
     { "CHEY", 4, 4L, IntPT, 0L        }, /* Front end channel for EY.        */
     { "CHHX", 4, 4L, IntPT, 0L        }, /* Front end channel for HX.        */
     { "CHHY", 4, 4L, IntPT, 0L        }, /* Front end channel for HY.        */
     { "CHHZ", 4, 4L, IntPT, 0L        }, /* Front end channel for HZ.        */
     { "EXAC", 0, 0L, FltPT, 0L        }, /* Ex AC between N-S in Volts.      */
     { "EXDC", 0, 0L, FltPT, 0L        }, /* Ex DC between N-S in Volts.      */
     { "EYAC", 0, 0L, FltPT, 0L        }, /* Ey AC between E-W in Volts.      */
     { "EYDC", 0, 0L, FltPT, 0L        }, /* Ey DC between E-W in Volts.      */
     { "EAZM", 1, 1L, FltPT, 0L        }, /* Ex sensor azimuth.               */
     { "EXLN", 1, 1L, FltPT, 0L        }, /* Ex dipole length, m.             */
     { "EYLN", 1, 1L, FltPT, 0L        }, /* Ey dipole length, m.             */
     { "HAZM", 1, 1L, FltPT, 0L        }, /* Hx sensor azimuth.               */
     { "HXSN", 1, 1L, StrPT, BlnkStrPT }, /* MAG sensor serial number.        */
     { "BSEX", 1, 1L, IntPT, 0L        }, /* MAG sensor base line setting.    */
     { "BSEY", 1, 1L, IntPT, 0L        }, /* MAG sensor base line setting.    */
     { "BSEZ", 1, 1L, IntPT, 0L        }, /* MAG sensor base line setting.    */
     { "TOTL", 1, 1L, IntPT, BlnkIntPT }, /* Total number of 1 s data records */
     { "BADR", 1, 1L, IntPT, BlnkIntPT }, /* Number of records flagged bad.   */
     { "SATR", 1, 1L, IntPT, BlnkIntPT }, /* Number of saturated records.     */
     { "NOBU", 1, 0L, IntPT, 0L        }, /* Not available buffer count       */
     { "NSAT", 2, 2L, IntPT, BlnkIntPT }, /* Number of GPS satellites         */
     { "LFIX", 2, 2L, AmxPT, 0L        }, /* UTC time of last GPS fix.        */
     { "CLST", 2, 2L, IntPT, BlnkIntPT }, /* Clock Status                     */
                                          /* = 0 if the time is uninitialized,*/
                                          /* = 1 if time based on CPU RTC,    */
                                          /* = 2 if time based on a TCXO which*/
                                          /*     was initialized by GPS,      */
                                          /* = 3 if time based on OCXO which  */
                                          /*     was initialized by GPS,      */
                                          /* = 4 if clock locked to GPS.      */
     { "OCTR", 2, 2L, IntPT, BlnkIntPT }, /* Oscillator control register      */  
                                          /* This is always an integral second*/
     { "TERR", 2, 2L, IntPT, BlnkIntPT }, /* Error in sample time, us.        */
                                          /* + for late, - for early.         */      
                                          /* Normally =0, will be non-zero    */      
                                          /* during recovery from GPS dropout.*/
     { "ELEV", 2, 2L, IntPT, BlnkIntPT }, /* Elevation altitude of site [m].  */
     { "LATG", 2, 2L, PosPT, 0L        }, /* GPS Latitude.                    */
     { "LNGG", 2, 2L, PosPT, 0L        }, /* GPS longitude.                   */                                             
     { "CALS", 3, 3L, IntPT, CAL_NOTPRESENT }, /* Box Calibration status.     */
     { "CCLS", 3, 3L, IntPT, CAL_NOTPRESENT }, /* Coil calibration status.    */
     { "CCLT", 3, 3L, IntPT, 2L        }, /* Coil cal time multiplier.        */
     { "CFMN", 3, 3L, FltPT, 0         }, /* Minimum acceptable corner freq.  */
     { "CFMX", 3, 3L, FltPT, 0         }, /* Maximum acceptable corner freq.  */                                          
     { "CCMN", 3, 3L, FltPT, 0         }, /* Coil cal min acceptable corner F.*/
     { "CCMX", 3, 3L, FltPT, 0         }, /* Coil cal max acceptable corner F.*/
     { "HATT", 3, 3L, FltPT, 0         }, /* Gain of coil att. on interconnect*/
                                          /* board, mV/V.                     */
     { "HNOM", 3, 3L, FltPT, 0         }, /* Coil nominal gain, mV/nT.        */
     { "HAMP", 3, 3L, FltPT, 0         }, /* Coil test waveform amplitude, nT.*/
     { "TCMB", 3, 3L, IntPT, 2L        }, /* Type of comb filter.             */
     { "TALS", 3, 3L, IntPT, 1L        }, /* Type of aliasing filter.         */
     { "BAT1", 3, 3L, IntPT, LOW_BATT+1}, /* Battery 1 voltage in mV.         */
     { "BAT2", 3, 3L, IntPT, 0L        }, /* Battery 2 voltage in mV.         */
     { "BAT3", 3, 3L, IntPT, 0L        }, /* Battery 3 voltage in mV.         */
     { "TEMP", 3, 3L, IntPT, 0L        }, /* Temperature in degrees C.        */
     { "GFPG", 3, 3L, IntPT, -1L       }, /* GPS FPGA programming return code.*/
     { "FFPG", 3, 3L, IntPT, -1L       }, /* Front End board(s) FPGA code.    */
     { "DSP",  3, 3L, IntPT, -1L       }, /* Front End board(s) DSP code.     */
     { "DISK", 3, 3L, IntPT, 32000000L }, /* Disk free space, bytes.          */
     { "\3",   3, 3L, IntPT, 0L        }  /* Marking end of table.            */
   };      

#define NENTRIES        (sizeof PrmTab/sizeof (TableEntryPT))

#include "Flags.c"                     /* Code for saving/restoring CPU flags.*/

/* ReadPT()     Function gets one entry from the parameter table.             */
/* Saves CPU interrupt state, disables interrupts, the restores CPU state.    */

int cdecl ReadPT(const char *code, TableEntryPT *ptp)
{
   int I;

   for(I=0; I < NENTRIES; I++) {
      if(stricmp(code, PrmTab[I].Code) == 0)
         break;                        /* Found entry in PrmTab.              */
   }
   if(I >= NENTRIES) {
      return(CERR);                    /* End of table, Code not found.       */
   }
   else {
      PushFlags();                     /* Save CPU flags.                     */
      ajdi();                          /* Disable interrupts.                 */
      *ptp = PrmTab[I];                /* Get the entry (structure)           */
      PopFlags();                      /* Restore CPU flags.                  */
      return 0;
   }
}
    

/* WritePT()      Function sets one entry in the parameter table.             */
/* Saves CPU interrupt state, disables interrupts, the restores CPU state.    */

int cdecl WritePT(const char *code, TableEntryPT *ptp)
{
   int I;

   for(I=0; I < NENTRIES; I++) {
      if(stricmp(code, PrmTab[I].Code) == 0)
         break;                        /* Found entry in PrmTab.              */
   }
   if(I >= NENTRIES) {
      return(CERR);                    /* End of table, Code not found.       */
   }
   else {
      if(ptp->T != PrmTab[I].T) {
         return(TERR);                 /* Wrong type.                         */
      }
      else {
         PushFlags();                  /* Save CPU flags.                     */
         ajdi();                       /* Disable interrupts.                 */
         PrmTab[I].V = ptp->V;         /* Write new value of entry.           */
         PopFlags();                   /* Restore CPU flags.                  */
         return 0;
      }
   }

}                    
                                                               
        
/* GCodePT      function gets PrmTab[index].Code from the parameter table       */
/*              *code is a pointer to a char array of 5 or larger               */
/*              AMX semaphore is not used by this function                      */

int     cdecl   GCodePT(const int index, char *code)
{
    if (index > NENTRIES-2)     /* the last entry is not used   */
        return(CERR);
    else {
        strcpy(code, PrmTab[index].Code);
        return(0);
    }
}


/* SavePT       function saves the PrmTab in a flashdisk file                   */
/*              the file name will be "FILE.TBL"                                */
int cdecl SavePT(void)
{
    int er,result;
    TableEntryPT pte;
    FILE        * TblFil;
    char *FileName = "               ";     /* name of parameter file */
    struct find_t FileInfo;
    unsigned int Error;
    
    er = ReadPT("FILE", &pte);
    strncpy(FileName, pte.V.S, 8);
    if( FileName[0] == NULL )          /* Do not write file if no filename    */
      return 0;                        /* specified.                          */
    strcat(FileName, ".TBL");

    ajpdrs();                          /* Check to see if file exists.        */
    Error = _dos_findfirst( FileName, _A_NORMAL, &FileInfo );
    ajpdrl();
    
    if( Error == 0 ) {                 /* .TBL file exists. Change filename to*/
      CreateFilename();                /* avoid overwriting.                  */
      er = ReadPT("FILE", &pte);
      strncpy(FileName, pte.V.S, 8);
      strcat(FileName, ".TBL");
    }
    
    ajpdrs();
    TblFil = fopen(FileName, "wb");
    ajpdrl();

    ajpdrs();             
    result = fwrite(PrmTab, sizeof(PrmTab), 1, TblFil);
    fclose(TblFil);                    /* Close the parameter file.           */
    ajpdrl();

    return (result);
}


/* Restart procedure for loading the parameter table.                         */

void cdecl PrmTblRP( void )
{
   ajtrig( GetPTid );
}
   

/* Task to load the parameter table from a startup.tbl file in \data.         */
/* Parameters are loaded one at a time and are checked to determine if they   */
/* are writeable parameters. If not, parameter is skipped. Writing is done    */
/* using WritePT().                                                           */

void cdecl GetPT( void )
{
   FILE *Fp;                           /* Startup.tbl file pointer.           */
   TableEntryPT PT, TempPT;            /* One parameter table item.           */
   int Error = 0;
   int ReadErr;
   int HW_Retries;                     /* Number of times to read the HW      */      
                                       /* parameter before giving up.         */

                                       /* Read the HW parameter. If NULL,     */
                                       /* keep reading until set.             */
   for( HW_Retries = 1; HW_Retries <= 10; HW_Retries++ ) {                                       
      ReadPT( "HW", &PT );                
      if( stricmp("", PT.V.S) != 0 )
         break;
      ajwatm(ajtmcnv(2000L));          /* Wait 2 seconds.                     */
   }

                                       /* If this is an MTU-CLK, there is no  */
                                       /* need to read startup.tbl.           */
   if( stricmp("mtuclk", PT.V.S) == 0 ||
       stricmp("mtucl4", PT.V.S) == 0 )
      return;

   ReadPT( "HXSN", &PT );              /* Place default coil names into table.*/
   strcpy( PT.V.S, "CHX" );
   WritePT( "HXSN", &PT );
   
   Fp = NULL;
   ajpdrs();
   Fp = fopen("\\data\\startup.tbl", "rb");
   ajpdrl();

   if( Fp != NULL ) {

      ajpdrs();
      printf("\nReading startup.tbl...");
      ajpdrl();      

      do {
         ajpdrs();                     /* Read one table entry.               */
         ReadErr = fread( &PT, sizeof( TableEntryPT ), 1, Fp ); 
         ajpdrl();

         if( ReadErr == 0 )
            break;
                                       /* Only permit the following parameters*/
                                       /* to be written at startup.           */
         if( ! strcmp( PT.Code, "EGN"  ) ||
             ! strcmp( PT.Code, "FILE" ) ||
             ! strcmp( PT.Code, "HGN"  ) ||
             ! strcmp( PT.Code, "LFRQ" ) ||
             ! strcmp( PT.Code, "RQST" ) ||
             ! strcmp( PT.Code, "SGIN" ) ||
             ! strcmp( PT.Code, "SITE" ) ||
             ! strcmp( PT.Code, "STM1" ) ||
             ! strcmp( PT.Code, "STM2" ) ||              
             ! strcmp( PT.Code, "STM3" ) ||              
             ! strcmp( PT.Code, "SRW1" ) ||              
             ! strcmp( PT.Code, "SRW2" ) ||              
             ! strcmp( PT.Code, "SRW3" ) ||              
             ! strcmp( PT.Code, "XDOS" ) ||
             ! strcmp( PT.Code, "ETM1" ) ||
             ! strcmp( PT.Code, "ETM2" ) ||                
             ! strcmp( PT.Code, "ETM3" ) ||                
             ! strcmp( PT.Code, "CMPY" ) ||
             ! strcmp( PT.Code, "SRVY" ) ||
             ! strcmp( PT.Code, "EAZM" ) ||
             ! strcmp( PT.Code, "EXLN" ) ||
             ! strcmp( PT.Code, "EYLN" ) ||
             ! strcmp( PT.Code, "HAZM" ) ||
             ! strcmp( PT.Code, "HXSN" )
             ) {

            ReadPT( PT.Code, &TempPT );/* Get current parameter.              */

            switch( TempPT.T ) {       /* Check type and set the new value    */
                                       /* appropriately.                      */
               case IntPT: TempPT.V.I = PT.V.I;
                           break; 

               case FltPT: TempPT.V.F = PT.V.F;
                           break;

               case StrPT: strcpy( TempPT.V.S, PT.V.S );
                           break;

               case UTCPT: break;
               
               case PosPT: strcpy( TempPT.V.P, PT.V.P );
                           break;
               
               case AmxPT: TempPT.V.TD = PT.V.TD;
                           break;

               default: Error = 1;     /* Bad type.                           */
            }                                                      
            if( !Error )               /* Do not write parameter if error.    */           
               WritePT( PT.Code, &TempPT );
            else
               Error = 0;              /* Clear error.                        */
         }
                                       /* Keep reading parameters until end of*/
                                       /* table.                              */
      } while( strcmp( PT.Code, "\3" ) );

   fclose( Fp );
   }

   else {
      ajpdrs();
      printf("\nNo startup.tbl");
      ajpdrl();
   }

   /* Set defaults for parameters which are not integer types.                */
   
   ReadPT( "CFMN", &PT );              /* Set min and max corner frequencies. */
   PT.V.F = 0.0;
   WritePT( "CFMN", &PT );

   ReadPT( "CFMX", &PT );
   PT.V.F = 0.0;
   WritePT( "CFMX", &PT );

   ReadPT( "CCMN", &PT );              /* Set coil min and max corner freqs.  */
   PT.V.F = 0.2;
   WritePT( "CCMN", &PT );   

   ReadPT( "CCMX", &PT );
   PT.V.F = 0.4;
   WritePT( "CCMX", &PT );

   ReadPT( "HATT", &PT );              /* Set coil attenuator gain.           */
   PT.V.F = 1.0;
   WritePT( "HATT", &PT );

   ReadPT( "HNOM", &PT );              /* Set coil nominal gain.              */
   PT.V.F = 4;
   WritePT( "HNOM", &PT );

   ReadPT( "HAMP", &PT );              /* Set coil test waveform amplitude.   */
   PT.V.F = -0.206;
   WritePT( "HAMP", &PT );   

}

       

