/*  File: vgptypes_0_1.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 28 08:14 2019 (rd109)
 * Created: Tue Jun 11 23:23:10 2019 (rd109)
 *-------------------------------------------------------------------
 */


/*** primary types and subtypes ***/

typedef enum { SEQ = 1, RMP, ALN, JNS, BRK, LIS, MAX_FILE } FileType ;
static char *fileTypeName[] = { 0, "seq", "rmp", "aln", "jns", "brk", "lis" } ;
typedef enum { IRP = 1, PBR, X10, CTG, RMM, RMS, RMA, SXS, RXR, SXR, MAP, LYO, SCF, MAX_SUB } SubType ;
static char *subTypeName[] = { 0, "irp", "pbr", "10x", "ctg", "rmm", "rms", "rma", "sxs", "rxr", "sxr", "map", "lyo", "scf" } ;
static FileType subPrimary[] = { 0, SEQ, SEQ, SEQ, SEQ, RMP, RMP, RMP, ALN, ALN, ALN, ALN, LIS, LIS } ;

#define MAX_FIELD 6	/* may need to increase this if future linetypes have more elements */

/**************** end of file ****************/
