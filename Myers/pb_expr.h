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

#define HAS_ZM  0x01   // Aux value flags
#define HAS_RQ  0x02
#define HAS_BC  0x04
#define HAS_BQ  0x08
#define HAS_NP  0x10
#define HAS_QS  0x20
#define HAS_QE  0x40

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
    int    defined;  //  Flags indicating which aux field values were set
  } samRecord;

Filter *parse_filter(char *expr, int *need);    //  Not re-entrant

int evaluate_bam_filter(Filter *v, samRecord *s);   //  Is re-entrant

#endif // _FILTER_EXPR
