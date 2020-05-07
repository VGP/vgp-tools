/*****************************************************************************************
 *
 *  File: VGPformat.h
 *      Description: this contains the VGP format specification.  The idea is to include
 *      this file into VGPlib.c, so it is compiled into the executable.  It is a separate
 *      file so as to make it easy to change the format and version implementation for VGPlib.h
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019-
 *
 * HISTORY:
 * Last edited: May  6 23:29 2020 (rd109)
 *   * Dec 27 09:46 2019 (gene): style edits
 *   * Jul  7 22:14 2019 (rd109): add DNAcodec for sequence S data
 *   * Jul  7 22:13 2019 (rd109): added code to build auxiliary structures, pack etc. last 2 days
 *   Created: Sun Feb 24 14:48:21 2019 (rd109)
 *
 ****************************************************************************************/


static LineInfo *vgpDefineLine (FieldType f0, FieldType f1, FieldType f2,
			        FieldType f3, FieldType f4, FieldType f5)
{ LineInfo *li = new0 (1, LineInfo);
  li->fieldSpec[0] = f0;
  li->fieldSpec[1] = f1;
  li->fieldSpec[2] = f2;
  li->fieldSpec[3] = f3;
  li->fieldSpec[4] = f4;
  li->fieldSpec[5] = f5;
  return (li);
}

static void defineFormat(VgpFile *vf, FileType fileType)
{ LineInfo **info;

  if (vf == NULL)
    die("vgpDefineFormat() called without a file");

  vf->major    = MAJOR_NUMBER;
  vf->minor    = MINOR_NUMBER;
  vf->fileType = fileType;

  info = vf->lineInfo;

  switch(fileType)

  { case SEQ:
      vf->objectType = 'S';
      vf->groupType  = 'g';
      info['g'] = vgpDefineLine (INT, STRING, 0, 0, 0, 0);      // group number name
      info['S'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);        // the sequence
      info['Q'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);        // qualities ascii 33+q
      info['P'] = vgpDefineLine (0, 0, 0, 0, 0, 0);             // start of a pair
      info['W'] = vgpDefineLine (INT, INT, INT, REAL, 0, 0);    // well, pulse start & end, score
      info['N'] = vgpDefineLine (REAL, REAL, REAL, REAL, 0, 0); // SNR in A,C,G,T channels
      info['A'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);        // capped pulse widths 1-4
      info['C'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // count (for kmers)
      info['I'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);        // for the identifier
        info['S']->listCodec = DNAcodec;
          info['S']->isUseListCodec = true;
        info['Q']->listCodec = vcCreate();
        info['A']->listCodec = vcCreate();
        info['W']->fieldCodec = vcCreate();
        info['N']->fieldCodec = vcCreate();
        info['S']->fieldCodec = vcCreate();
        info['Q']->fieldCodec = vcCreate();
        info['C']->fieldCodec = vcCreate();
      break;

    case RMP:
      vf->objectType = 'R';
      vf->groupType  = 'r';
      info['r'] = vgpDefineLine (INT, STRING_LIST, 0, 0, 0, 0); // number restriction_sites
      info['R'] = vgpDefineLine (INT, INT_LIST, 0, 0, 0, 0);    // len, locations (in bp)
      info['E'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // sites in list in r line
      info['I'] = vgpDefineLine (REAL_LIST, 0, 0, 0, 0, 0);     // intensities at each site
      info['N'] = vgpDefineLine (REAL_LIST, 0, 0, 0, 0, 0);     // SNR values at each site
      info['O'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // object # in ref'd sequence file
        info['R']->listCodec = vcCreate();
          info['R']->isIntListDiff = true;
        info['E']->listCodec = vcCreate();
          info['E']->isIntListDiff = true;
        info['I']->listCodec = vcCreate();
        info['N']->listCodec = vcCreate();
      break;

    case ALN:
      vf->objectType = 'A';
      vf->groupType  = 'g';
      info['g'] = vgpDefineLine (INT, STRING, 0, 0, 0, 0);      // group number name
      info['A'] = vgpDefineLine (INT, INT, 0, 0, 0, 0);         // object numbers of aligned objects
      info['I'] = vgpDefineLine (INT, INT, INT, INT, INT, INT); // as ae alen bs be blen
      info['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // quality in phred units
      info['M'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // # of matching bases
      info['D'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // # of diffs = subs + indels
      info['C'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);        // cigar string
      info['T'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);           // trace point spacing
      info['U'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // tracePoints in a
      info['V'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // tracePoints in b
      info['W'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // tracePoint spacings in b
      info['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // inter-tracePoint diff counts in b
        info['I']->fieldCodec = vcCreate();
        info['C']->listCodec = vcCreate();
        info['U']->listCodec = vcCreate();
          info['U']->isIntListDiff = true;
        info['V']->listCodec = vcCreate();
          info['V']->isIntListDiff = true;
        info['W']->listCodec = vcCreate();
          info['W']->isIntListDiff = true;
        info['X']->listCodec = vcCreate();
          info['X']->isIntListDiff = true;
      break;

    case HIT:
      vf->objectType = 'H';
      info['H'] = vgpDefineLine (INT, INT_LIST, 0, 0, 0, 0);    // Hits between objects
      info['O'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // Offsets of queries in a target
      info['P'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);      // Positions of a query in targets
        info['H']->listCodec = vcCreate();
          info['H']->isIntListDiff = true;
        info['O']->listCodec = vcCreate();
          info['O']->isIntListDiff = true;
        info['P']->listCodec = vcCreate();
          info['P']->isIntListDiff = true;
      break;

    case JNS:
      vf->objectType = 'J';
      info['J'] = vgpDefineLine (INT, INT, CHAR, INT, INT, CHAR); // a pos_a [s|e] b pos_b [s|e]
      info['G'] = vgpDefineLine (INT, INT, 0, 0, 0, 0);           // gap size estimate & deviation
      info['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);             // confidence in phred units
      info['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);        // alignments supporting join
      break;

    case BRK:
      vf->objectType = 'B';
      info['B'] = vgpDefineLine (INT, INT, INT, 0, 0, 0);  // object start end
      info['Q'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);      // confidence in phred units
      info['X'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0); // alignment objects supporting join
      break;

    case LIS:
      vf->objectType = 'L';
      info['L'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0); // object identifiers
      info['N'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);   // optional name for list
      info['S'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);      // seed sequence for scaffold
        info['L']->listCodec = vcCreate();
          info['L']->isIntListDiff = true;
      break;

    default:
      break;
  }
}

/********* end of file ***********/
