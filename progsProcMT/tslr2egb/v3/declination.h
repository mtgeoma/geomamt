#ifndef DECLINATION_H
#define DECLINATION_H

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>                         /* standard C include statements */
#include <math.h>
#include <string.h>
#include <ctype.h>

using namespace std;

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define IEXT 0
#define FALSE 0
#define TRUE 1                                                 /* constants */
#define RECL 81

#define MAXINBUFF RECL+14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF-2
/** Max to read 2 less than total size (just to be safe) **/

#define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define EXT_COEFF1 (float)0
#define EXT_COEFF2 (float)0
#define EXT_COEFF3 (float)0

int declination(float latitude, float longitude, float elev, int isyear, int ismonth, int isday);

#endif // DECLINATION_H
