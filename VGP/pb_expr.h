/*******************************************************************************************
 *
 *  Filter Expression Parser & Evaluator
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 31, 2016
 *
 ********************************************************************************************/

#ifndef _FILTER_EXPR
#define _FILTER_EXPR

#include "pb_expr.h"

typedef void *Filter;

typedef struct
  { int    len;
    int    well;
    int    beg;
    int    end;
    float  qual;
    float  snr[4];
    int    bc[2];
    int    bqual;
    int    nump;
    int    lmax;     //  current size of seq and arr
    char  *seq;
    char  *arr;
    char  *qvs;
    int    dmax;     //  current size of data
    uint8 *data;     //  data buffer
    char  *header;
  } samRecord;

Filter *parse_filter(char *expr);    //  Not re-entrant

int evaluate_bam_filter(Filter *v, samRecord *s);   //  Is re-entrant

#endif // _FILTER_EXPR
