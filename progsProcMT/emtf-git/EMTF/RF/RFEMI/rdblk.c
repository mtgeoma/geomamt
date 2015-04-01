#include <stdio.h>
#define NCMAX 100
#define NMAX 512

static FILE *fp ;

/* routine to read in next block of EMI data ... then read "trailing block" */
/* returns: +1 for "OK" read  ... 0 for EOF ... -1 for RESET block" */

long rdblk_ ( int *nch , int *n , short *dat )
{
char c[NCMAX];         /* storage for characters in trailing block */
unsigned char dat1[NMAX*2];    /* temporary storage array for one channel of data */
int i , j , nctrail , ierr;
long i2;

i2 = 2;

for ( i = 0 ; i < *nch ; ++i )  {
   ierr = fread(dat1,1,*n*2,fp);        /* read in data blocks for nch channels*/
/*   printf("after reading channel %d , ierr = %d , ichar = %d %d\n",i,ierr,dat1[0],dat1[1]);*/
   if ( feof(fp) > 0 ) return 0 ;
   ierr = fseek(fp,i2,SEEK_CUR);      /* skip new line characters at end of each channel block */
   for ( j = 0 ; j < *n ; ++j )
     dat[*nch*j+i] = 256*((unsigned short)dat1[2*j]-128)+(unsigned short)dat1[2*j+1];  /* copy data into short integer array : dat(nch,n) */
   }

nctrail = 0;
while ( ( c[nctrail] = fgetc(fp) ) != '\r')      /* read trailing block */
  nctrail = nctrail + 1;
fseek(fp,1,SEEK_CUR);       /* skip line feed; position at start of next block */
for (i = 0 ; i < nctrail ; ++i )
  if(c[i] == 'O' ) return 1 ;
return -1 ;
}

/*    open file for reading with C routine rdblk   */
int file_open_c_ ( int *nskip )
{
char cfile[NCMAX];
int ierr ;
FILE *fid ;

printf(" in file_open_c\n");
fid = fopen( "file_name_temp", "r" ) ;
ierr = fscanf(fid,"%s",cfile) ;
printf("file name =  %s \n",cfile);
printf("ierr = %d \n", ierr);
printf("nskip = %d \n", *nskip);
ierr = fclose(fid) ;
fp = fopen(cfile,"rb") ;
printf("fp = %d \n",fp);
ierr = fseek(fp,*nskip,SEEK_SET) ;
printf("ierr = %d \n",ierr);
return ierr ;
}
