/*  File: vgpformat_0_1.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: this contains the VGP format specification 0.1
 *   the idea is to include this file into vgprd.c, so it is compiled into the executable
 *   it is a separate file so as to make it easy to change the format and version
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 18 18:20 2019 (rd109)
 * * Jul  7 22:14 2019 (rd109): add DNAcodec for sequence S data
 * * Jul  7 22:13 2019 (rd109): added code to build auxiliary structures, pack etc. last 2 days
 * Created: Sun Feb 24 14:48:21 2019 (rd109)
 *-------------------------------------------------------------------
 */

static LineInfo *vgpDefineLine (FieldType f0, FieldType f1, FieldType f2,
				FieldType f3, FieldType f4, FieldType f5)
{ LineInfo *li = new0 (1, LineInfo) ;
  li->fieldSpec[0] = f0 ; li->fieldSpec[1] = f1 ; li->fieldSpec[2] = f2 ;
  li->fieldSpec[3] = f3 ; li->fieldSpec[4] = f4 ; li->fieldSpec[5] = f5 ;
  return li ;
}

static void defineFormat (VgpFile *vf, FileType fileType)
{
  if (!vf) die ("vgpDefineFormat() called without a file") ;

  vf->major = 1 ; vf->minor = 0 ;

  vf->fileType = fileType ;
  switch (fileType)
    {
    case SEQ:
      vf->objectType = 'S' ;
      vf->groupType = 'g' ;
      vf->lineInfo['g'] = vgpDefineLine (INT, STRING, 0, 0, 0, 0) ; /* group number name */
      vf->lineInfo['S'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ; /* the sequence */
      vf->lineInfo['S']->listCodec = DNAcodec ; vf->lineInfo['S']->isUseListCodec = TRUE ;
      vf->lineInfo['Q'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ; /* qualities ascii 33+q */
      vf->lineInfo['Q']->listCodec = vcCreate() ;
      vf->lineInfo['P'] = vgpDefineLine (0, 0, 0, 0, 0, 0) ; /* start of a pair */
      vf->lineInfo['W'] = vgpDefineLine (INT, INT, INT, REAL, 0, 0) ; /* well pulseStart pulseEnd score */
      vf->lineInfo['N'] = vgpDefineLine (REAL, REAL, REAL, REAL, 0, 0) ; /* SNR in A,C,G,T channels */
      vf->lineInfo['A'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ; /* capped pulse widths 1-4 */
      vf->lineInfo['A']->listCodec = vcCreate() ;
      break ;

    case RMP:
      vf->objectType = 'R' ;
      vf->groupType = 'r' ;
      vf->lineInfo['r'] = vgpDefineLine (INT, STRING_LIST, 0, 0, 0, 0) ; /* number restriction_sites */
      vf->lineInfo['R'] = vgpDefineLine (INT, INT_LIST, 0, 0, 0, 0) ; /* len locations(bp) */
      vf->lineInfo['R']->isIntListDiff = TRUE ;
      vf->lineInfo['E'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* sites in list in r line */
      vf->lineInfo['R']->listCodec = vcCreate() ;
      vf->lineInfo['I'] = vgpDefineLine (REAL_LIST, 0, 0, 0, 0, 0) ; /* intensities at each site */
      vf->lineInfo['I']->listCodec = vcCreate() ;
      vf->lineInfo['N'] = vgpDefineLine (REAL_LIST, 0, 0, 0, 0, 0) ; /* SNR_values at each site */
      vf->lineInfo['N']->listCodec = vcCreate() ;
      vf->lineInfo['O'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* object number in referred sequence file */
      break ;

    case ALN:
      vf->objectType = 'A' ;
      vf->groupType = 'g' ;
      vf->lineInfo['g'] = vgpDefineLine (INT, STRING, 0, 0, 0, 0) ; /* group number name */
      vf->lineInfo['A'] = vgpDefineLine (INT, INT, 0, 0, 0, 0) ; /* object numbers of aligned objects */
      vf->lineInfo['I'] = vgpDefineLine (INT, INT, INT, INT, INT, INT) ; /* as ae alen bs be blen */
      vf->lineInfo['I']->fieldCodec = vcCreate() ;
      vf->lineInfo['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* quality in phred units */
      vf->lineInfo['M'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* number of matching bases */
      vf->lineInfo['D'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* number of differences = substitutions + indel bases */
      vf->lineInfo['C'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ; /* cigar string */
      vf->lineInfo['C']->listCodec = vcCreate() ;
      vf->lineInfo['T'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* tracePoint spacing in a (global) */
      vf->lineInfo['U'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* tracePoints in a */
      vf->lineInfo['U']->listCodec = vcCreate() ;
      vf->lineInfo['U']->isIntListDiff = TRUE ;
      vf->lineInfo['V'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* tracePoints in b */
      vf->lineInfo['V']->listCodec = vcCreate() ;
      vf->lineInfo['V']->isIntListDiff = TRUE ;
      vf->lineInfo['W'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* tracePoint spacings in b */
      vf->lineInfo['W']->listCodec = vcCreate() ;
      vf->lineInfo['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* inter-tracePoint diff counts in b */
      vf->lineInfo['X']->listCodec = vcCreate() ;
      break ;

    case JNS:
      vf->objectType = 'J' ;
      vf->lineInfo['J'] = vgpDefineLine (INT, INT, CHAR, INT, INT, CHAR) ; /* a pos_a [s|e] b pos_b [s|e] */
      vf->lineInfo['G'] = vgpDefineLine (INT, INT, 0, 0, 0, 0) ; /* mean and standard deviation of estimated gap size */
      vf->lineInfo['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* confidence in phred units */
      vf->lineInfo['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* alignment objects supporting join */
      break ;

    case BRK:
      vf->objectType = 'B' ;
      vf->lineInfo['B'] = vgpDefineLine (INT, INT, INT, 0, 0, 0) ; /* object start end */
      vf->lineInfo['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* confidence in phred units */
      vf->lineInfo['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* alignment objects supporting join */
      break ;

    case LIS:
      vf->objectType = 'L' ;
      vf->lineInfo['L'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ; /* object identifiers */
      vf->lineInfo['L']->listCodec = vcCreate() ;
      vf->lineInfo['L']->isIntListDiff = TRUE ;
      vf->lineInfo['N'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ; /* optional name for list */
      vf->lineInfo['S'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* seed sequence for scaffold */
      break ;

    default: break ;
    }
}

/********* end of file ***********/
