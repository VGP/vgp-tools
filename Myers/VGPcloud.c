#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <dirent.h>

#include "gene_core.h"
#include "radix.h"
#include "exsort.h"

int    VERBOSE;    //  Verbose mode?
char  *SORT_PATH;  //  Director to do external sort in

#define VALID_THRESH 100

static char *Usage = "[-v] [-P<dir(/tmp)>] [-T<int(4)>] <input:irp>";

static int Header[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

static int Value[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

static int IsN[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

static uint32 zero[16] =
  { 0xfffffffcu, 0xffffff3u, 0xffffffcfu, 0xffffff3fu,
    0xfffffcffu, 0xffff3ffu, 0xffffcfffu, 0xffff3fffu,
    0xfffcffffu, 0xff3ffffu, 0xffcfffffu, 0xff3fffffu,
    0xfcffffffu, 0x3ffffffu, 0xcfffffffu, 0x3fffffffu
  };

static uint32 one[16] =
  { 0x00000001u, 0x0000004u, 0x00000010u, 0x00000040u,
    0x00000100u, 0x0000400u, 0x00001000u, 0x00004000u,
    0x00010000u, 0x0040000u, 0x00100000u, 0x00400000u,
    0x01000000u, 0x4000000u, 0x10000000u, 0x40000000u
  };

static uint8 bit[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

static char DNA[4] = { 'a', 'c', 'g', 't' };

static inline void Check_Entry(int hasfs, int hasrs, int hasfq, int hasrq)
{ if (! hasfs)
    { fprintf(stderr,"%s: Entry does not have forward sequence\n",Prog_Name);
      exit (1);
    }
  if (! hasrs)
    { fprintf(stderr,"%s: Entry does not have reverse sequence\n",Prog_Name);
      exit (1);
    }
  if (! hasfq)
    { fprintf(stderr,"%s: Entry does not have forward QVs\n",Prog_Name);
      exit (1);
    }
  if (! hasrq)
    { fprintf(stderr,"%s: Entry does not have reverse QVs\n",Prog_Name);
      exit (1);
    }
}

  //  For barcodes that are either not valid or have N's or QV's < 10
  //    search for a unique correction, and if found change the barcode
  //    in the forward sequence to said and give QV 0 to all changed bases.

static int find(uint32 code, int l, int h, uint32 *count)
{ int m;

  while (h-l > 10)
    { m = (h-l)/2;
      if (m & 0x1)
        m -= 1;
      if (count[m] > code)
        h = m;
      else
        l = m;
    }
  for (m = l; m < h; m += 2)
    if (count[m] == code)
      break;
  return (m+1);
}

static int Correction(uint8 *forw, uint8 *fqvs, uint8 *valid, uint32 *count, int *lookup)
{ int    nbad, stack[16];
  uint32 bar, var, cor;
  uint32 v3h, v3b;
  int    i, pos, yes;

  nbad = 0;
  bar  = 0;
  for (i = 0; i < 16; i++)
    { if (IsN[forw[i]] || fqvs[i] < '+')
        stack[nbad++] = 2*i;
      bar = (bar << 2) | Value[forw[i]];
    }
  if (nbad > 0)
    { int    j;
      uint32 val;

      yes = 0;
      while (1)
        { for (j = 0; j < nbad; j++)
            { pos = stack[j];
              val = (bar >> 2*pos) & 0x3;
              if (val == 0x3)
                bar = bar & zero[pos];
              else
                { bar += one[pos];
                  break;
                }
            }
          if (j >= nbad)
            break;
          v3h = bar>>3;
          v3b = bit[bar&0x3];
          if (valid[v3h] & v3b)
            { if (yes)
                return (0);
              yes = 1;
              var = bar;
            }
        }
      cor = var;
      for (i = 15; i >= 0; i--)
        { val = var & 0x3;
          if (Value[forw[i]] != (int) val)
            { forw[i] = DNA[val];
              fqvs[i] = 0;
            }
          var >>= 2;
        }
    }
  else if ((valid[bar>>3] & bit[bar&0x7]) == 0)
    { int    chr;
      uint32 chg;
      int    i, k;

      yes = 0;
      for (k = 0; k < 4; k++)
        { chg = k;
          for (i = 0; i < 16; i++)
            { var = (bar & zero[i]) | chg;
              v3h = var>>3;
              v3b = bit[var&0x3];
              if (valid[v3h] & v3b)
                { if (yes)
                    return (0);
                  yes = 1;
                  pos = i;
                  chr = k;
                  cor = var;
                }
              chg <<= 2;
            }
        }
      forw[pos] = DNA[chr];
      fqvs[pos] = 0;
    }
  else
    return (1);

  var = (cor >> 16);
  return (++count[find(cor,lookup[var],lookup[var+1],count)]);
}

static void Compress_QV(int len, uint8 *qvec, int qbits, uint8 *map)
{ uint8 *buf, x;
  int    i, rem;

  buf = qvec;
  rem = 0;
  for (i = 0; i < len; i++)
    { x = map[qvec[i]];
      *buf |= (x << rem);
      rem += qbits;
      if (rem >= 8)
        { rem -= 8;
          if (rem > 0)
            *++buf = (x >> (qbits-rem));
          else
            *++buf = 0;
        }
    }
}

static void Uncompress_QV(int len, uint8 *qvec, int qbits, uint8 qmask, uint8 *inv)
{ uint8 *buf, x, v, bmask;
  int    i, rem;

  bmask = (1 << qbits) - 1;
  buf = qvec + (len*qbits-1)/8;
  rem = (len*qbits-1) % 8 + 1;
  x   = *buf--;
  for (i = len-1; i >= 0; i--)
    { if (qbits > rem)
        { v = (x << (qbits-rem));
          rem += 8-qbits;
          x = *buf--;
          qvec[i] = v | (x >> rem);
        }
      else
        qvec[i] = inv[(x >> (rem-qbits)) & qmask];
    }
}

int main(int argc, char *argv[])
{ FILE *input;
  char *fname;        // Input file and name
  int   NTHREADS;

  char *provenance;     // All previous provenance lines
  int   nprov, mprov;   // # of provenance lines and width of maximum line (excluding !)
  int   sprov, aprov;

  uint8  *forw, *revr;  // Forward and reverse sequence & QV buffers & lengths
  uint8  *fqvs, *rqvs;
  int     flen,  rlen;

  uint32 *count;        // Barcode list and eventually valid barcode list & their counts
  int    *lookup;       // Index base on high-order 2 bytes into list of valid barcodes
  int     npair;        // # of read pairs in file

  uint8   map[256];     // Map from QV char to bit code
  uint8   inv[256];     // Map from bit code to QV char
  int     qbits;        // # of bits to encode QV values
  uint8   qmask;        // mask of lowest qbits ina a byte

  uint8  *valid;        // 4^16 bit vector of good codes
  int     ngood;        // # of good barcodes
  int     ndiff;        // # of distinct good barcodes
  int     nused;        // # of pairs that get clustered (ngood <= nused <= npair)
  int     gmax;         // Count of largest group

  int   reclen;         // Compressed record length
  int   fclen, rclen;   // Length of compressed fields
  int   fqlen, rqlen;


  //  Parse command line options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    DIR   *dirp;

    ARG_INIT("VGPcloud")

    SORT_PATH = "/tmp";
    NTHREADS  = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            if ((dirp = opendir(SORT_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -P: Do external sorts in this directory.\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
 
        exit (1);
      }
  }


  //  Open input file

  { int   i;
    char *pwd;
    char *suffix[3] = { ".irp.gz", ".irp", ".gz" };

    pwd = PathTo(argv[1]);
    OPEN(argv[1],pwd,fname,input,suffix,3)
    free(pwd);
  }


  //  Scan 1: Sort barcodes

  if (VERBOSE)
    { fprintf(stderr,"  Scanning barcodes in file %s\n",argv[1]);
      fflush(stderr);
    }

  { int     maxlen;

    //  Header section: Verify is a .irp, get max sequence length & number of pairs 

    { char code, which, *pptr;
      int  len;

      provenance = NULL;
      nprov = mprov = sprov = aprov = 0;
      while (fscanf(input," %c",&code) == 1)       //  Header lines
        if (Header[(int) code])
          { switch (code)
            { case '@':
                fscanf(input," %c",&which);
                if (which == 'S')
                  fscanf(input," %d",&maxlen);
                else if (which == '!')
                  { fscanf(input," %d",&mprov);
                    mprov = 4*(mprov+11)+1;
                  }
                break;
              case '#':
                fscanf(input," %c",&which);
                if (which == 'P')
                  fscanf(input," %d",&npair);
                break;
              case '+':
                fscanf(input," %c",&which);
                break;
              case '1':
              case '2':
                fscanf(input," %d",&len);
                if (len == 3)
                  { char tname[4];
                    fscanf(input," %s",tname);
                    if (code == '1')
                      { if (strcmp(tname,"seq") == 0)
                          break;
                      }
                    else
                      { if (strcmp(tname,"irp") == 0)
                          break;
                      }
                  }
                if (code == '1')
                  fprintf(stderr,"%s: File primary type is not .seq",Prog_Name);
                else
                  fprintf(stderr,"%s: File secondary type is not .irp",Prog_Name);
                exit (1);
              case '!':
                if (nprov <= 0)
                  { fprintf(stderr,"%s: Provenance count is not given\n",Prog_Name);
                    exit (1);
                  }
                if (mprov <= 0)
                  { fprintf(stderr,"%s: Provenance max string is not given\n",Prog_Name);
                    exit (1);
                  }
                if (provenance == NULL)
                  pptr = provenance = (char *) Malloc(nprov*mprov,"Provenance cache");

                { char *p;
                  int   i, d;

                  p = pptr;
                  for (i = 0; i < 4; i++)
                    { fscanf(input," %d",&d);
                      p += sprintf(p," %d ",d);
                      fscanf(input," %s",p);
                      p += d;
                      sprov += d;
                      if (d > aprov)
                        aprov = d;
                    }
                  *p = '\0';
                  pptr += mprov; 
                }
                break;
              default:
                break;
            }
            fscanf(input,"%*[^\n]\n");
          }
        else
          { ungetc(code,input);
            break;
          }
    }

    //  Allocate sequence buffers and barcode sorting array

    forw = (uint8 *) Malloc(sizeof(uint8)*4*(maxlen+1),"Allocating sequence buffers");
    revr = forw + (maxlen+1);
    fqvs = revr + (maxlen+1);
    rqvs = fqvs + (maxlen+1);

    count = (uint32 *) Malloc(sizeof(uint32)*npair,"Allocating barcode array");


    //  Data segment: build list of all barcodes & determine QV-value bit compression mapping

    { char code;
      int  hasfs, hasrs, hasfq, hasrq;
      int  dobar;
      int  usedqvs[256];
      int  i, n;

      bzero(usedqvs,sizeof(int)*256);

      npair = 0;
      flen  = rlen = -1;
      dobar = 0;
      while (fscanf(input," %c",&code) == 1)       //  For each data line do
        { switch (code)
          { case 'P':
              if (npair != 0)
                Check_Entry(hasfs,hasrs,hasfq,hasrq);
              hasfs = hasrs = hasfq = hasrq = 0;
              dobar = 1;
              break;
            case 'S':
              { int len;

                fscanf(input," %d ",&len);
                if (!hasfs)
                  { hasfs = 1;
                    if (flen < 0)
                      flen = len;
                    if (flen != len)
                      { fprintf(stderr,"%s: Forward reads are not all the same length\n",Prog_Name);
                        exit (1);
                      }
                    fscanf(input," %16c",forw);
                  }
                else if (!hasrs)
                  { hasrs = 1;
                    if (rlen < 0)
                      rlen = len;
                    if (rlen != len)
                      { fprintf(stderr,"%s: Reverse reads are not all the same length\n",Prog_Name);
                        exit (1);
                      }
                  }
                else
                  { fprintf(stderr,"%s: Entry has more than 2 S-lines?\n",Prog_Name);
                    exit (1);
                  }
                fscanf(input,"%*[^\n]\n");
              }
              break;
            case 'Q':
              { int i, len;
  
                fscanf(input," %d ",&len);
                if (!hasfq)
                  { hasfq = 1;
                    if (flen < 0)
                      flen = len;
                    if (flen != len)
                      { fprintf(stderr,"%s: Forward QV vectors are not all the same length\n",
                                       Prog_Name);
                        exit (1);
                      }
                    fread(fqvs,1,len,input);
                    for (i = 0; i < len; i++)
                      usedqvs[(int) fqvs[i]] += 1;
                  }
                else if (!hasrq)
                  { hasrq = 1;
                    if (rlen < 0)
                      rlen = len;
                    if (rlen != len)
                      { fprintf(stderr,"%s: Reverse QV vectors are not all the same length\n",
                                       Prog_Name);
                        exit (1);
                      }
                    fread(rqvs,1,len,input);
                    for (i = 0; i < len; i++)
                      usedqvs[(int) rqvs[i]] += 1;
                  }
                else
                  { fprintf(stderr,"%s: Entry has more than 2 Q-lines?\n",Prog_Name);
                    exit (1);
                  }
                fscanf(input,"%*[^\n]\n");
              }
              break;
          }

          if (hasfs && hasfq && dobar)
            { int    i, good;
              uint32 bar;

              dobar = 0;
              good  = 1;
              bar = 0;

              for (i = 0; i < 16; i++)
                { if (IsN[forw[i]] || fqvs[i] < '+')
                    good = 0;
                  bar = (bar << 2) | Value[forw[i]];
                }

              if (good)
                count[npair++] = bar;
            }
        }
      Check_Entry(hasfs,hasrs,hasfq,hasrq);

      n = 0;                         //  Map non-empy QV histogram entries to bit code
      for (i = 0; i < 256; i++)      //    and determine # of bits to encode largest
        if (usedqvs[i])
          { map[i] = n;
            inv[n] = i;
            n += 1;
          }
      for (qbits = 0; n > 0; qbits++)
        n >>= 1;
      qmask  = (1 << qbits) - 1;
      map[0] = 0;   //  Make sure 0-qvs intro'd by Correction get turned to min val after decompr.
    }

  }  //  End of Scan 1


  //  Sort barcode list,then compress to list of valid codes, and finally
  //    build bit vector 'valid' of good codes

  if (VERBOSE)
    { fprintf(stderr,"  About to sort %d barcodes\n",npair);
      fflush(stderr);
    }

  { int    i, g, barsort[5];
    uint32 bar, lst;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
    for (i = 0; i < 4; i++)
      barsort[i] = i;
#else
    for (i = 0; i < 4; i++)
      barsort[i] = 3-i;
#endif
    barsort[4] = -1;

    Set_Radix_Params(6,VERBOSE);

    Radix_Sort(npair,count,count+npair,barsort);

    count[npair] = count[npair] - 1;

    gmax  = 0;
    ngood = 0;
    ndiff = 0;
    g = 0;
    for (i = 1; i <= npair; i++)
      if (count[i] != count[g])
        { g = i-g;
          if (g > VALID_THRESH)
            { count[ndiff++] = count[g];
              count[ndiff++] = g;
              ngood += g;
              if (g > gmax)
                gmax = g;
            }
          g = i;
        }

    count  = Realloc(count,sizeof(uint32)*ndiff,"Resizing valid codes vector");
    valid  = (uint8 *) Malloc(0x20000000ll,"Bit vectors");
    lookup = (int *) Malloc(sizeof(int)*0x10001,"Code index");
    if (count == NULL || valid == NULL || lookup == NULL)
      exit (1);

    bzero(valid,0x20000000ll);

    lst = (count[0] >> 16)-1;
    for (i = 0; i < ndiff; i += 2)
      { bar = count[i];
        valid[bar>>3] |= bit[bar&0x7];
        bar >>= 16;
        while (lst < bar)
          { lst += 1;
            lookup[lst] = i;
          }
      }
    while (lst < 0x10001u)
      { lst += 1;
        lookup[lst] = i;
      }

    if (VERBOSE)
      { fprintf(stderr,"  There are %d valid codes with %d distinct values.\n",ngood,ndiff>>1);
        fflush(stderr);
      }
  }


  //  Scan 2: Correct barcodes when possible, output compressed pairs for 

  { char  code;
    int   nextQ, nextS;
    int   nline, len, cnt;
    FILE *sfile;

    if (VERBOSE)
      { fprintf(stderr,"  Scan to compressing pairs and repair barcodes when possible.\n");
        fflush(stderr);
      }

    while (fscanf(input," %c",&code) == 1)       //  Skip over header lines
      if (Header[(int) code])
        fscanf(input,"%*[^\n]\n");
      else
        { ungetc(code,input);
          break;
        }

    sfile = fopen(Catenate(SORT_PATH,"/",fname,".sort"),"w");
    if (sfile == NULL)
      { fprintf(stderr,"%s: Cannot create %s.sort in directory %s\n",Prog_Name,fname,SORT_PATH);
        exit (1);
      }

    fclen  = COMPRESSED_LEN(flen);
    rclen  = COMPRESSED_LEN(rlen);
    fqlen  = (flen*qbits-1)/8 + 1; 
    rqlen  = (rlen*qbits-1)/8 + 1; 
    reclen = fclen + rclen + fqlen + rqlen;

    nused = 0;
    nextS = 0;
    nextQ = 0;
    while (fscanf(input," %c",&code) == 1)       //  For each data line do
      { if (code == 'P')
          { nextS = nextQ = 1;
            nline = 0;
          }
        else
          { fscanf(input," %d ",&len);
            if (code == 'S')
              { if (nextS)
                  fread(forw,sizeof(char),len,input);
                else
                  fread(revr,sizeof(char),len,input);
                nextS = 0;
              }
            else
              { if (nextQ)
                  fread(fqvs,sizeof(char),len,input);
                else
                  fread(rqvs,sizeof(char),len,input);
                nextQ = 0;
              }
            nline += 1;
          }
        fscanf(input,"%*[^\n]\n");

        if (nline < 4)
          continue;
         
        cnt = Correction(forw,fqvs,valid,count,lookup);
        if (cnt > 0)
          { if (cnt > gmax)
              gmax = cnt;

            Number_Read((char *) forw);
            Compress_Read(flen,(char *) forw);
            fwrite(forw,sizeof(uint8),fclen,sfile);

            Number_Read((char *) revr);
            Compress_Read(rlen,(char *) revr);
            fwrite(revr,sizeof(uint8),rclen,sfile);

            Compress_QV(fclen,fqvs,qbits,map);
            fwrite(fqvs,sizeof(uint8),fqlen,sfile);

            Compress_QV(rclen,rqvs,qbits,map);
            fwrite(rqvs,sizeof(uint8),fqlen,sfile);

            nused += 1;
          }
      }

    fclose(sfile);

    if (VERBOSE)
      { fprintf(stderr,"  %d pairs had good codes, %d repairable codes, %d dropped.\n",
                       ngood,nused-ngood,npair-nused);
        fflush(stderr);
      }
  }

  //  External sort on barcode (1st 4 bytes of each record)

  if (VERBOSE)
    { fprintf(stderr,"  Performing external sort.\n");
      fflush(stderr);
    }

  Ex_sort(Catenate(SORT_PATH,"/",fname,".sort"),reclen,4,NTHREADS);

  //  Read in sorted array and output to standard out

  if (VERBOSE)
    { fprintf(stderr,"  Final scan to produce cloud grouped pairs on stdout.\n");
      fflush(stderr);
    }

  //  Output header

  { char   date[26];
    time_t seconds;
    int    i, clen;

    ndiff >>= 1;
    seconds = time(NULL);
    ctime_r(&seconds,date);
    date[24] = '\0';
    clen = strlen(argv[1]);

    printf(". 10x 1.0\n");
    printf("# ! %d\n",nprov+1);
    printf("# g %d\n",ndiff);
    printf("# P %d\n",nused);
    printf("# S %d\n",2*nused);
    printf("# Q %d\n",2*nused);

    sprov += clen + 35; 
    printf("+ ! %d\n",sprov);
    printf("+ g %d\n",ndiff*16);
    printf("+ S %d\n",nused*((flen-23)+rlen));
    printf("+ Q %d\n",nused*((flen-23)+rlen));

    if (24 > aprov)
      aprov = 24;
    if (clen > aprov)
      aprov = clen;
    printf("@ ! %d\n",aprov);
    printf("@ g 16\n");
    if (flen-23 > rlen)
      { printf("@ S %d\n",flen-23);
        printf("@ Q %d\n",flen-23);
      }
    else
      { printf("@ S %d\n",rlen);
        printf("@ Q %d\n",flen);
      }
    
    printf("%% g # P %d\n",gmax);
    printf("%% g # P %d\n",gmax);
    printf("%% g # Q %d\n",gmax*2);
    printf("%% g + S %d\n",gmax*((flen-23)+rlen));
    printf("%% g + Q %d\n",gmax*((flen-23)+rlen));

    for (i = 0; i < nprov; i++)
      printf("! %s\n",provenance + i*mprov);
    printf("! 8 VGPcloud 3 1.0 %d %s 24 %s\n",clen,argv[1],date);
  }

  { FILE  *sfile;
    int    gc, yes;

    sfile = fopen(Catenate(SORT_PATH,"/",fname,".sort"),"w");
    if (sfile == NULL)
      { fprintf(stderr,"%s: Cannot create %s.sort in directory %s\n",Prog_Name,fname,SORT_PATH);
        exit (1);
      }

    yes = 0;
    gc  = 0;
    while (fread(forw,sizeof(uint8),fclen,sfile) < (uint32) fclen)
      { if ( *((uint32 *) forw) == count[gc])
          { printf("g %d",count[gc+1]);
            gc += 2;
            yes = 1;
          }

        Uncompress_Read(flen,(char *) forw);

        fread(revr,sizeof(uint8),rclen,sfile);
        Uncompress_Read(rlen,(char *) revr);

        fread(fqvs,sizeof(uint8),fqlen,sfile);
        Uncompress_QV(fclen,fqvs,qbits,qmask,inv);

        fread(rqvs,sizeof(uint8),fqlen,sfile);
        Uncompress_QV(rclen,rqvs,qbits,qmask,inv);

        forw[flen] = revr[rlen] = 4;
        Lower_Read((char *) forw);
        Lower_Read((char *) revr);

        if (yes)
          { yes = 0;
            printf(" 16 %.16s\n",forw);
          }

        printf("P\n");
        printf("S %d %s\n",flen-23,forw+23);
        printf("Q %d %s\n",flen-23,fqvs+23);
        printf("S %d %s\n",rlen,revr);
        printf("Q %d %s\n",rlen,rqvs);
      }

    fclose(sfile);
  }

  //  Tidy up just for good form

  fclose(input);
  free(fname);
  exit (0);
}
