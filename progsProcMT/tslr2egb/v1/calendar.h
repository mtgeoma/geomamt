/*
 * convcal : dates conversion utility
 *
 * Copyright (c) 1999 Luc Maisonobe
 *
 *                           All Rights Reserved
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * This programs allows you to convert dates between calendar format
 * and numerical format.

 * The following command will compile the program :
 *  cc -o convcal convcal.c -lm

 */

#define lib 	0
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

#define REFDATE "-4713-01-01T12:00:00"

typedef enum   { FMT_iso,
                 FMT_european,
                 FMT_us,
                 FMT_days,
                 FMT_seconds,
                 FMT_nohint
               } Dates_format;

typedef struct { int value;
                 int digits;
               } Int_token;


/*
 * set of functions to convert julian calendar elements
 * with negative years to julian day
 */
static int neg_julian_non_leap (int year);

static long neg_julian_cal_to_jul(int y, int m, int d);

static int neg_julian_year_estimate(long n);

/*
 * set of functions to convert julian calendar elements
 * with positive years to julian day
 */
static int pos_julian_non_leap(int year);

static long pos_julian_cal_to_jul(int y, int m, int d);

static int pos_julian_year_estimate(long n);

/*
 * set of functions to convert gregorian calendar elements to julian day
 */
static int gregorian_non_leap(int year);

static long gregorian_cal_to_jul(int y, int m, int d);

static int gregorian_year_estimate(long n);

/*
 * convert calendar elements to Julian day
 */
long cal_to_jul(int y, int m, int d);

/*
 * convert julian day to calendar elements
 */
static void jul_to_some_cal(long n,
                            int (*some_non_leap) (int),
                            long (*some_cal_to_jul) (int, int, int),
                            int (*some_year_estimate) (long),
                            int *y, int *m, int *d);

/*
 * convert julian day to calendar elements
 */
void jul_to_cal(long n, int *y, int *m, int *d);

/*
 * convert julian day and hourly elements to julian day
 */
double jul_and_time_to_jul(long jul, int hour, int min, double sec);

/*
 * convert calendar and hourly elements to julian day
 */
double cal_and_time_to_jul(int y, int m, int d,
                           int hour, int min, double sec);

/*
 * convert julian day to calendar and hourly elements
 * rounding_tol allows to say 1999-12-31T23:59:59.501
 * should be rounded to 2000-01-01T00:00:00.000 assuming
 * it is set to 0.5 second. It is wise to set it according
 * to the display accuracy of seconds.
 */
void jul_to_cal_and_time(double jday, double rounding_tol,
                         int *y, int *m, int *d,
                         int *hour, int *min, double *sec);

/*
 * check the existence of given calendar elements
 * this includes either number of day in the month
 * and calendars pecularities (year 0 and October 1582)
 */
static int check_date(int century, int wy,
                      Int_token y, Int_token m, Int_token d,
                      long *jul);

/*
 * lexical analyser for float data (knows about fortran exponent
 * markers, return address of following data)
 */
int parse_float(const char* s, double *value, const char **after);

/*
 * lexical analyser for calendar dates
 * return the number of read elements, or -1 on failure
 */
static int parse_calendar_date(const char* s,
                               Int_token tab [5], double *sec);

/*
 * parse a date given either in calendar or numerical format
 */
int parse_date(const char* s, int century, int wy, Dates_format preferred,
               double *jul, Dates_format *recognized);

int convert_and_write(const char *s,
                      int century, int wy, double reference_date,
                      Dates_format input_format, Dates_format output_format);

int string_equal(const char *c1, const char *c2);

int parse_format(const char *s, Dates_format *f);

/*
 * expand a line buffer
 */

void expand_line_buffer(char **adrBuf, int *ptrSize, char **adrPtr);

/*
 * help message
 */

void usage (FILE *stream, const char *progname);
