/*
   Fragment simulator: generates simulated shotgun data.
 
   Author:            Gene Myers

   Date last Revised: October 12, 1998

   Desciption:
      Got rid of interactive input mode.  Input is now from a file
        parameter and not standard input.
      Got rid of old 4 parameter DNA sequence convention.
      Made it so number of elements can be expressed as a percent of basis.

   * Make it so that # of reads or # of repeats can be expressed as an X
     of genome.
   * Make it so that you can give a series of ramps for error ramp.
   * Make it so that part of repeats can be generated.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

double drand48();
void   srand48();

void  *malloc();

#define MAX_INT 0x7FFFFFFF

int   Seed;             /* Seed for random number generation        */

FILE *sfile;            /* Unit to read spec from. */

int   no_elem;          /* Don't print out sub-elements */
int   comments;         /* Don't print out any comments */


/* >>>> REPETITIVE ELEMENT DEFINITIONS <<<< */

/* An element is either of "type" BASIS or CNCAT, with or without the
   "GLOBAL" bit set.  If of type BASIS then the "len" and "prob" fields
   give the length range and symbol distribution for the basis sequence.
   "rlist" links together all of the element references in the given
   elements definition (in order).  "valstack" is used during sequence
   construction to maintain a stack of current instances for the given
   element.  "instances" points to a list (linked by the "ilist" field
   of the INSTANCE records) of instances of the given element used in
   building the subject strand.  */

#define BASIS  0
#define CNCAT  1
#define GLOBAL 2

typedef struct { struct refr *rlist;
                 int          type;
                 int          minlen, maxlen;
                 double       probA, probAC, probACG;
                 struct inst *valstack;
                 struct inst *instances;
               } DEFINITION;

/* Each reference is to an element ("name" in [0,25]).  The remaining
  fields record the orientation probability, mutation rate range,
  repetition factor range, and unique regeneration range, respectively.
  maxins = 0 if regeneration is not to take place.  */

typedef struct refr { struct refr *rlist;
                      int          name;
                      double       orient;
                      double       minerr, maxerr;
                      int          minrep, maxrep;
                      double       xminrp, xmaxrp;
                      int          minins, maxins;
                    } REFERENCE;

/* Instances of a given element are linked by their "ilist" field, and
   the "elist" field points to a list of the component instances contained
   in the instance.  The length and unique instance number of the instance
   are given in the fields "length" and "level".  "value" points at the
   string giving the instance sequence (only valid during construction
   save for that of the instance that becomes the subject strand).  "nextval"
   links together the stack of instances used during construction. */

typedef struct inst { struct comp *elist;
                      int          length;
                      int          level;
                      char        *value;
                      struct inst *nextval;
                      struct inst *ilist;
                    } INSTANCE;

/* "refer" and "level" identify the instance constituting this component.
   "length" is the components length (including mutations).  "oprand"
   points at the string value of the instance from which this component
   is derived (valid only during construction).  "forw" gives the orientation
   in which the underlying instance was incorporated, "erate" gives the
   mutation rate, and "errors" the number of mutations introduced.  Finally,
   "begpos" and "endpos" specify the beginning and ending indices within
   the containing instance at which this component was incorporated.  */

typedef struct comp { struct comp *elist;
                      int          length;
                      int          level;
                      struct refr *refer;
                      struct inst *source;
                      char        *oprand;
                      int          forw;
                      double       erate;
                      int          errors;
                      int          begpos, endpos;
                    } COMPONENT;
        
DEFINITION Relement[26];

int   Sequence;  /* Element defining string to be fragmented */
FILE *SeqFile;   /* if Sequence < 0 ==> File containing string to be frag'd */


/* >>>> BASIC UTILITIES <<<< */

#define SEGLEN    50    /* Output sequences, SEGLEN chars to a line */
#define DEL_FREQ .333   /* Deletion frequency for mutations */
#define INS_FREQ .333   /* Insertion frequency for mutations */

char SubChar[8][3] = { {'c', 'g', 't'},  /* What chars to substitute */
                       {'a', 'g', 't'},
                       {'a', 'c', 't'},
                       {'a', 'c', 'g'},
                       {'C', 'G', 'T'},
                       {'A', 'G', 'T'},
                       {'A', 'C', 'T'},
                       {'A', 'C', 'G'} };

int Decode[128];     /* Converts chars a,c,g,t to integers 0,1,2,3;
                        and A,C,G,T to integers 4,5,6,7.             */

char FlipChar[128];  /* Gives the Watson-Crick complement */
char TranChar[128];  /* The identity transform */

  /* Allocate space and exit if none */

char *ckalloc(int size)
{ char *m;
  if ((m = (char *) malloc(size)) == NULL)
    { fprintf(stdout,"Out of Memory (frag %d)\n",size);
      exit (2);
    }
  return (m);
}

  /* Generate a character with a coin toss */

char genchar(double pa, double pac, double pacg)
{ double rnd;

  rnd = drand48();
  if (rnd < pac)
    if (rnd < pa)
      return ('a');
    else
      return ('c');
  else
    if (rnd < pacg)
      return ('g');
    else 
      return ('t');
}

  /* Return a substitute for the input character. */

char substitute(int orig)
{ int i;
  i = 2.9999*drand48();
  return (SubChar[Decode[orig]][i]);
}

 
/* >>>> SUBJECT SEQUENCE CONSTRUCTION <<<< */

/*  Place a string at s which is a mutated copy of the string t of
    length len.  The copy is in the direction given by fwd, mutated
    with err changes, ins of which are insertions, and del of which
    are deletions.     */

char *mutate(char *s, char *t, int len, int fwd, int err, int ins, int del) 
{ double x;
  char *tr;

  len += (del - ins);
  del += ins;
  tr   = TranChar;
  if (!fwd)
    { fwd = -1;
      t  += (len-1);
      tr  = FlipChar;
    }

#ifdef DEBUGS
  printf("mutate = '");
#endif
  while (len > 0 || err > 0)
    if ((x = err*drand48()) < ins)
      { if ((len+1)*drand48() < err)
          { *s++ = genchar(.25,.50,.75);
            del -= 1;
            ins -= 1;
            err -= 1;
#ifdef DEBUGS
            printf("!%c",(s[-1]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[*t];
            t += fwd;
            len -= 1;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
      }
    else if (x < del)
      { if (len*drand48() < err)
          { t += fwd;
            del -= 1;
            err -= 1;
#ifdef DEBUGS
            printf(".%c",(t[-fwd]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[*t];
            t += fwd;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
        len -= 1;
      }
    else
      { if (len*drand48() < err)
          { *s++ = substitute(tr[*t]);
            t += fwd;
            err -= 1;
#ifdef DEBUGS
            printf("%c",(s[-1]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[*t];
            t += fwd;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
        len -= 1;
      }
#ifdef DEBUGS
  printf("'\n");
#endif
  return (s);
}

/* Build a string for version lev of element e.  This is done by
   recursively building strings for each reference in e's definition
   (if that version has not yet been built), then determining the
   orientation, mutation level, and number of instances for each
   reference (so the size of e's string is known), and then building
   the string for e from these instance descriptions.   */

static int level;

int COMPARE(COMPONENT **a, COMPONENT **b)
{ return ((*a)->forw - (*b)->forw); }

void build_trav(DEFINITION *e, int lev)
{ REFERENCE  *r;
  DEFINITION *f;
  INSTANCE   *i, *u, *v;
  COMPONENT  *x, *y, **z;
  char *s;
  int n, l;
  int reps, len, num, ins, avl;
  double rv, xrp;
  INSTANCE uhead;
  int *insters;

  ins = 0;
  for (r = e->rlist; r != NULL; r = r->rlist)
    if (r->maxins > 0) ins += 1;

  if (ins > 0)
    insters = (int *) malloc(sizeof(int)*ins);
  else
    insters = NULL;

  ins = 0;
  u = &uhead;
  for (r = e->rlist; r != NULL; r = r->rlist)
    { f = Relement + r->name;
      if (r->maxins > 0)
        { len = r->minins + drand48()*(r->maxins - r->minins + .9999);
          insters[ins++] = len;
          while (len-- > 0)
            { build_trav(f,l = ++level);
              u->nextval  = f->valstack;
              u = u->nextval;
              f->valstack = u->nextval;
              for (n = 0; n < 26; n++)
                { v = Relement[n].valstack;
                  if (v != NULL && v->level == l)
                    { if (no_elem) free(v->value);
                      Relement[n].valstack = v->nextval;
                    }
                }
            }
        }
      else if (f->valstack == NULL ||
               (f->valstack->level != lev && !(f->type & GLOBAL)))
        { if (f->type & GLOBAL)
            build_trav(f,0);
          else
            build_trav(f,lev);
        }
    }

  i = (INSTANCE *) ckalloc(sizeof(INSTANCE));
  i->level     = lev;
  i->ilist     = e->instances;
  e->instances = i;
  i->nextval   = e->valstack;
  e->valstack  = i;
  if ((e->type & 1) == BASIS)
    i->length = e->minlen + drand48()*(e->maxlen - e->minlen);
#ifdef DEBUG
  printf(" Exec: %c %d\n",'A' + (e-Relement),lev); fflush(stdout);
#endif

  len = 0;
  num = 0;
  ins = 0;
  u = uhead.nextval;
  y = (COMPONENT *) i;
  for (r = e->rlist; r != NULL; r = r->rlist)
    { f = Relement + r->name;
      if (r->maxins > 0)
        n = insters[ins++];
      else
        n = 1;
      while (n-- > 0)
        { if (r->maxins > 0)
            { v = u;
              u = u->nextval;
            }
          else
            v = f->valstack;
          if (r->minrep >= 0)
            reps = r->minrep + drand48()*(r->maxrep - r->minrep + .9999);
          else
            { xrp = r->xminrp + drand48()*(r->xmaxrp - r->xminrp);
              avl = v->length * (1. +
                     (INS_FREQ-DEL_FREQ) * (r->minerr + r->maxerr) / 2.);
              reps = rv = ((int) (i->length * xrp )) / (1.*avl);
              if (rv - reps > drand48())
                reps += 1;
            }
          while (reps-- > 0)
            { x = (COMPONENT *) ckalloc(sizeof(COMPONENT));
              x->level  = v->level;
              x->refer  = r;
              x->source = v;
              x->oprand = v->value;
              x->erate  = r->minerr + drand48()*(r->maxerr - r->minerr);
              x->errors = rv = v->length * x->erate;
              if (rv - x->errors > drand48())
                x->errors += 1;
              x->begpos = rv = x->errors * INS_FREQ;
              if (rv - x->begpos > drand48())
                x->begpos += 1;
              if (INS_FREQ > 1 - 1e-20)
                x->endpos = 0;
              else
                { rv = (x->errors - x->begpos)*DEL_FREQ/(1.-INS_FREQ);
                  x->endpos = rv;
                  if (rv - x->endpos > drand48())
                    x->endpos += 1;
                }
              x->length = v->length + x->begpos - x->endpos;
              y = y->elist  = x;
              len += x->length;
              num += 1;
#ifdef DEBUG
              printf("\t%c.%d: o=%d e=%g(%d,%d,%d) ||= %d\n",
                      x->refer->name+'A',x->level,x->forw,x->erate,x->errors,
                      x->begpos,x->endpos,x->length); fflush(stdout);
#endif
            }
        }
    }
  y->elist = NULL;
  if ((e->type & 1) == CNCAT)
    i->length = len;
  i->value = (char *) ckalloc(i->length+1);

  s = i->value;
  if ((e->type & 1) == CNCAT)
    for (x = i->elist; x != NULL; x = x->elist)
      { n = s - i->value;
        x->forw = (drand48() <= x->refer->orient);
        s = mutate(s,x->oprand,x->length,x->forw,
                     x->errors,x->begpos,x->endpos);
        x->begpos = n;
        x->endpos = s - i->value;
      }
  else
    { if (num > 0)
        { z = (COMPONENT **) ckalloc(sizeof(COMPONENT *)*num);
          n = 0;
          for (x = i->elist; x != NULL; x = x->elist)
            { z[n++] = x;
              x->forw = drand48()*MAX_INT;
            }
          qsort(z,num,sizeof(COMPONENT *),COMPARE);
          y = (COMPONENT *) i;
          for (n = 0; n < num; n++)
            y = y->elist = z[n];
          y->elist = NULL;
          free(z);
        }
   
      len = i->length - len;
      x   = i->elist;
      while (len > 0 || num > 0)
        if (drand48() < num/(len+1.))
          { n = s - i->value;
            x->forw = (drand48() <= x->refer->orient);
            s = mutate(s,x->oprand,x->length,x->forw,
                         x->errors,x->begpos,x->endpos);
            x->begpos = n;
            x->endpos = s - i->value;
            x = x->elist;
            num -= 1;
          }
        else
          { *s++ = genchar(e->probA,e->probAC,e->probACG);
            len -= 1;
          }
    }
  *s = '\0';
#ifdef DEBUGS
  printf("Val = '%s'\n",i->value);
#endif

  if (no_elem)
    for (r = e->rlist; r != NULL; r = r->rlist)
      if (r->maxins > 0)
        { u = uhead.nextval;
          uhead.nextval = u->nextval;
          free(u->value);
        }

  if (insters != NULL) free(insters);
}

void report_trav(INSTANCE *i, int lev, int type, int pos)
{ COMPONENT *x;

  for (x = i->elist; x != NULL; x = x->elist)
    { printf("# %*s ",lev,"");
      if (type & 1)
        printf("= ");
      else
        printf("> ");
      printf("%c.%d",x->refer->name+'A',x->level);
      if (x->forw)
        printf("f");
      else
        printf("r");
      printf(" at %d-%d",pos+x->begpos,pos+x->endpos);
      printf(", mutated %.2f",x->erate);
      printf(".\n");
      report_trav(x->source,lev+2,Relement[x->refer->name].type,pos+x->begpos);
    }
}
      

/* Build the subject dna strand. */

void build_seq()
{ int n, k, len;
  INSTANCE  *i;
  COMPONENT *x;

  for (n = 0; n < 26; n++)
    { Relement[n].valstack  = NULL;
      Relement[n].instances = NULL;
    }
  level = 0;
  build_trav(Relement+Sequence,0);

  i = Relement[Sequence].instances;
  if (comments)
    { printf("#  %c.%d (%d)\n",Sequence+'A',i->level,i->length);
      report_trav(i,0,Relement[Sequence].type,0);
      printf("#\n");
    }

  if (!no_elem && comments)
    { for (n = 0; n < 26; n++)
        for (i = Relement[n].instances; i != NULL; i = i->ilist)
          if (n != Sequence)
            { len = i->length;
              printf("#  %c.%d (%d)\n",n+'A',i->level,len);
              for (k = 0; k < len; k += SEGLEN)
                if (len-k < SEGLEN)
                  printf("#     %.*s\n",len-k,i->value+k);
                else
                  printf("#     %.*s\n",SEGLEN,i->value+k);
            }
      printf("#\n");
    }
}


/* >>>> INPUTING AND PARSING A SIMULATION SPECIFICATION <<<< */

/* Any error terminates execution.  So to make error messages informative,
   we adopt the convention of writing out object descriptions as soon as they
   are correctly input.   */

/* Nextsymbol gets the next non blank character.

int nextsymbol(int c)
{ do
    c = fgetc(sfile);
  while (isspace(c));
  return (c);
}

/* Get a postive natural number from the sfile. */

int getnatural(char *label)
{ int datum;
 
  if (fscanf(sfile," %d",&datum) != 1 || datum <= 0)
    { fprintf(stdout,"\n\t*** %s must be a positive int\n",label);
      exit(1);
    }
  return(datum);
}

/* Get a probability from the sfile. */

double getprobability(char *label)
{ double datum;

  if (fscanf(sfile," %lf",&datum) != 1 || datum < 0. || datum > 1.)
    { fprintf(stdout,"\n\t*** %s must be a real number in [0,1]\n",label);
      exit(1);
    }
  return(datum);
}

/* Get a number from the sfile. */

int getfactor(char *label, double *real)
{ int value;
  double xvalue;
  int c;
  int cnt1, cnt2, expo;
  unsigned int frac;

  c = nextsymbol();
  ungetc(c,sfile);
  value = 0;
  if (fscanf(sfile,"%d",&value) != 1 && c != '.')
    { fprintf(stdout,
              "\n\t*** Expecting an integer or real number\n");
      exit(1);
    }
  c = getc(sfile);
  if (c == '.' || c == 'e' || c == 'E')
    { if (c == '.')
        { c = getc(sfile);
          ungetc(c,sfile);
          if (isdigit(c))
            fscanf(sfile,"%n%u%n",&cnt1,&frac,&cnt2);
          else
            cnt1 = cnt2 = frac = 0;
          if (cnt1 == cnt2 || frac == 0) 
            xvalue = 0.;
          else
            { while (frac % 10 == 0)
                { frac /= 10;
                  cnt1 += 1;
                }
              xvalue = frac;
              while (cnt1 < cnt2)
                { xvalue /= 10.;
                  cnt1 += 1;
                }
	    }
          xvalue += value;
        }
      else
        { ungetc(c,sfile);
          xvalue = value;
        }
      c = getc(sfile);
      if (c == 'e' || c == 'E')
        { if (fscanf(sfile,"%d",&expo) != 1)
            ungetc(c,sfile);
          else
            xvalue *= pow(10.,1.*expo);
        }
      else
        ungetc(c,sfile);
      if (xvalue < 0. || xvalue > 1.)
        { fprintf(stdout,"\n\t*** %s %% must be a real ",label);
          fprintf(stdout,"number in [0,1]\n");
          exit(1);
        }
      *real = xvalue;
      return (-1);
    }
  else
    { ungetc(c,sfile);
      if (value < 0)
        { fprintf(stdout,
              "\n\t*** %s must be a non-negative int\n",label);
          exit(1);
        }
      return (value);
    }
}

/* Get a basis definition for element e from the sfile.  Return the
   first non-blank symbol following the definition is returned.  */

void getbasis(DEFINITION *e)
{ int c;
  int minl, maxl;
  double pa, pac, pacg;

/* Get size range */

  minl = getnatural("(Min) Element length");
  c = nextsymbol();
  if (c == '-')
    { maxl = getnatural("Max Element length");
      c = nextsymbol();
    }
  else
    maxl = minl;
  if (minl > maxl)
    { fprintf(stdout,"\n\t*** Min length (%d) > Max length (%d)\n",minl,maxl);
      exit (1);
    }

/* Get symbol bias probabilities (if specified) */

  if (c == 'p')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      pa = getprobability("Odds of A");
      if ((c = nextsymbol()) != ',') ungetc(c,sfile);
      pac = getprobability("Odds of C");
      if ((c = nextsymbol()) != ',') ungetc(c,sfile);
      pacg = getprobability("Odds of G");
      if ((c = nextsymbol()) != ')') ungetc(c,sfile);
      if (pacg > 1.)
        { fprintf(stdout,"\n\t*** Sum (%.2f) of A|C|G odds > 1\n",pacg);
          exit (1);
        }
    }
  else
    { pa = pac = pacg = .25;
      ungetc(c,sfile);
    }

  e->minlen  = minl;   /* Establish record fields */
  e->maxlen  = maxl;
  e->probA   = pa;
  e->probAC  = pac += pa;
  e->probACG = pacg += pac;
}

/* Build and return another reference record for element e's definition. */

REFERENCE *getreference(DEFINITION *e)
{ REFERENCE *r;
  int c, m, minr, maxr, mini, maxi;
  double o, mine, maxe;
  double xminr, xmaxr;

  c = nextsymbol();    /* Get referenced element */
  if (! isupper(c))
    { fprintf(stdout,"\n\t*** %c is not an element name\n",c);
      exit (1);
    }
  m = c - 'A';
  if (Relement[m].type < 0)
    { fprintf(stdout,"\n\t*** %c referenced before being defined\n",c);
      exit (1);
    }
  c = nextsymbol();

  if (c == 'r')     /* Get orientation (if specified) */
    { o = 0.;
      c = nextsymbol();
    }
  else if (c == 'o')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      o = getprobability("Orientation odds");
      if ((c = nextsymbol()) == ')') c = nextsymbol();
    }
  else
    o = 1.;
                    /* Get mutation rate/range (if specified) */
  if (c == 'm')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      mine = maxe = getprobability("(Min) Variability");
      c = nextsymbol();
      if (c == '.' || isdigit(c) || c == ',')
        { if (c != ',') ungetc(c,sfile);
          maxe = getprobability("Max Variability");
          c = nextsymbol();
        }
      if (c == ')') c = nextsymbol();
    }
  else
    mine = maxe = 0.;
  if (mine > maxe)
    { fprintf(stdout,"\n\t*** Min variability (%lf) > Max variability (%lf)\n",
                     mine,maxe);
      exit (1);
    }

  if (c == 'n')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      minr = getfactor("(Min) occurence",&xminr);
      if (minr < 0)
        { if ((e->type & 1) != BASIS)
            { fprintf(stdout,"\n\t*** %% spec. only allowed in basis def's.\n");
              exit (1);
            }
          xmaxr = xminr;
          c = nextsymbol();
          if (isdigit(c) || c == ',' || c == '.')
            { if (c != ',') ungetc(c,sfile);
              xmaxr = getprobability("Max Occurence %");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
          if (xminr > xmaxr)
            { fprintf(stdout,"\n\t*** Min occurrence %% (%g) > ",xminr);
              fprintf(stdout,"Max occurrence %% (%g)\n",xmaxr);
              exit (1);
            }
        }
      else
        { maxr = minr;
          c = nextsymbol();
          if (isdigit(c) || c == ',')
            { if (c != ',') ungetc(c,sfile);
              maxr = getnatural("Max Occurences");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
          if (minr > maxr)
            { fprintf(stdout,"\n\t*** Min occurrences (%d) > ",minr);
              fprintf(stdout,"Max occurrences (%d)\n",maxr);
              exit (1);
            }
        }
    }
  else
    minr = maxr = 1;
  
  if (c == '!')
    { c = nextsymbol();
      if (c == '(' || isdigit(c))
        { if (c != '(') ungetc(c,sfile);
          if (fscanf(sfile," %d",&mini) != 1 || mini < 0)
            { fprintf(stdout,
                      "\n\t*** (Min) Instances must be a non-negative int\n");
              exit(1);
            }
          maxi = mini;
          c = nextsymbol();
          if (isdigit(c) || c == ',')
            { if (c != ',') ungetc(c,sfile);
              maxi = getnatural("Max Instances");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
        }
      else
        maxi = mini = 1;
    }
  else
    mini = maxi = 0;
  if (mini > maxi)
    { fprintf(stdout,"\n\t*** Min instances (%d) > Max instances (%d)\n",
                     mini,maxi);
      exit (1);
    }

  ungetc(c,sfile);
  
/* Build reference record */

  r = (REFERENCE *) ckalloc(sizeof(REFERENCE));
  r->name   = m;
  r->orient = o;
  r->minerr = mine;
  r->maxerr = maxe;
  r->minrep = minr;
  r->maxrep = maxr;
  r->xminrp = xminr;
  r->xmaxrp = xmaxr;
  r->minins = mini;
  r->maxins = maxi;

/* Print reference information */

  if (comments)
    { if ((e->type & 1) == BASIS && e->rlist == NULL)
        printf("#    Containing:\n");
      printf("#\t %c:",r->name + 'A');
      printf(" O'odds = %.2f,",r->orient);
      printf(" Mut's = %.2f-%.2f,",r->minerr,r->maxerr);
      if (r->minrep >= 0)
        printf(" %d-%d Rep's,",r->minrep,r->maxrep);
      else
        printf(" %g-%g %% of Basis,",100*r->xminrp,100*r->xmaxrp);
      printf(" %d-%d Gen's",r->minins,r->maxins);
      printf("\n");
      fflush(stdout);
    }

  return (r);
}


/* Read and parse all element definitions from the input */

int EditMax(int size, double error)
{ int n, q;

  n = size*error;
  if (n < size*error) n += 1;
  q = n*INS_FREQ;
  if (q < n*INS_FREQ) q += 1;
  size += q - ((int) (n*DEL_FREQ));
  return (size);
}

int EditMin(int size, double error)
{ int n, q;

  n = size*error;
  q = n*DEL_FREQ;
  if (q < n*DEL_FREQ) q += 1;
  size += ((int) (n*INS_FREQ)) - q;
  return (size);
}

void getgrammar()
{ int c;
  DEFINITION *e;
  REFERENCE  *r, *rlink;
  int name, type, mlen, llen;

  while (1)
                          /* Is there another definition? */
    { if (Sequence < 0)
        { if (comments) printf("# Seed = %d\n#\n",Seed); }
      else
        { c = nextsymbol();
          if (feof(sfile)) break;
          ungetc(c,sfile);
          if (isdigit(c)) break;
        }
      if (comments)
        { printf("# Element: ");
          fflush(stdout);
        }
                          /* Get name and basis (if any) */
      name = c = nextsymbol();
      if (! isupper(name))
        { if (c == '\n')
            fprintf(stdout,"\n\t*** Expecting a name\n");
          else
            fprintf(stdout,"\n\t*** %c is not an uppercase letter\n",name);
          exit (1);
        }
      else if (Relement[name-'A'].type >= 0)
        { fprintf(stdout,"\n\t*** %c is already defined\n",name);
          exit (1);
        }
      e = Relement + (name - 'A');
      c = nextsymbol();
      if (c != '=' && c != '~')
        { fprintf(stdout,"\n\t*** Expecting an = or ~ after %c\n",name);
          exit(1);
        }
      if (c == '=')
        type = 0;
      else
        type = GLOBAL;
      c = nextsymbol();
      if (isdigit(c))
        { type = BASIS | type;
          ungetc(c,sfile);
          getbasis(e);
          c = nextsymbol();
        }
      else
        { type = CNCAT | type;
          if (Sequence < 0)
            { fprintf(stdout,"\n\t*** 1st declaration must have a basis\n");
              exit(1);
            }
        }
      ungetc(c,sfile);

      e->type  = type;
      e->rlist = NULL;
      Sequence = name - 'A';

      if (comments)
        { printf("%c",name);
          if (e->type & GLOBAL)
            printf(" (static)");
          printf("\n");
          if ((e->type & 1) == BASIS)
            { printf("#    Length = [%d,%d], ACGT odds = %.2f/%.2f/%.2f/%.2f\n",
                     e->minlen,e->maxlen,e->probA,e->probAC - e->probA,
                     e->probACG - e->probAC, 1.0 - e->probACG);
            }
          else
            printf("#    Concatenation of:\n");
          fflush(stdout);
        }

      rlink = (REFERENCE *) e;    /* Get references */
      mlen  = llen = 0;
      while (1)
        { c = nextsymbol();
          if (c == ';')
            { if ((e->type & 1) == CNCAT && e->rlist == NULL)
                { fprintf(stdout,"\n\t *** Must define at least one element\n");
                  exit (1);
                }
              break;
            }
          else
            { ungetc(c,sfile);
              r = getreference(e);
              c = Relement[r->name].maxlen;
              if (r->minrep >= 0)
                { if (r->maxins)
                    { mlen += EditMax(c,r->maxerr)*r->maxrep*r->maxins;
                      llen += EditMin(Relement[r->name].minlen,r->minerr)
                              *r->minrep*r->minins;
                    }
                  else
                    { mlen += EditMax(c,r->maxerr)*r->maxrep;
                      llen += EditMin(Relement[r->name].minlen,r->minerr)
                              *r->minrep;
                    }
                }
              else
                { int x, m, n, s, l;
                  double a;

                  x = e->minlen*r->xmaxrp;
                  a = (1.+(INS_FREQ-DEL_FREQ)*(r->minerr+r->maxerr)/2.);
                  m = c * a;
                  n = x / m;
                  if (n*m < x) n += 1;
                  s = EditMax(c,r->maxerr)*n;
                  if (n*m == x)
                    l = c-1;
                  else
                    { l = x / (n*a);
                      m = l*a;
                      if (m * n == x) l -= 1;
                    }
                  c = Relement[r->name].minlen;
                  m = c*a;
                  if (l >= c)
                    l = EditMax(l,r->maxerr)*(n+1);
                  else if (m * (n+1) <= e->maxlen*r->xmaxrp)
                    l = (EditMax(c,r->maxerr)*x)/m;
                  else
                    l = s;
                  if (s < l) s = l;
                  if (r->maxins)
                    mlen += s*r->maxins;
                  else
                    mlen += s;
                }
              rlink = rlink->rlist = r;
              r->rlist = NULL;
            }
        }

      if ((e->type & 1) == BASIS)
        { if (e->minlen < mlen)
            { fprintf(stdout,"\n\t*** Min basis length(%d) < ",e->minlen);
              fprintf(stdout,"Possible sum of contained elements");
              fprintf(stdout," (%d)\n",mlen);
              exit (1);
            }
        }
      else
        { e->maxlen = mlen;
          e->minlen = llen;
          mlen = 0;
          for (r = e->rlist; r != NULL; r = r->rlist)
            mlen += r->minrep;
          if (mlen == 0)
            { fprintf(stdout,"\n\t*** Potentially empty concatenation\n");
              exit (1);
            }
        }
    }
}

void getfile()
{ int i, n, top;
  char *buffer;

  top     = 1024;
  buffer  = (char *) ckalloc(sizeof(char)*(top+1));
  SeqFile = NULL;
  if (comments)
    { printf("# Seed = %d\n#\n",Seed);
      fflush(stdout);
    }
  n = nextsymbol();
  n = nextsymbol();
  i = 0;
  buffer[i++] = n;
  n = fgetc(sfile);
  while ( ! isspace(n))
    { if (i >= top)
        { buffer = (char *) realloc(buffer,sizeof(char)*(2*top + 1));
          if (buffer == NULL)
            { fprintf(stdout,"Out of Memory (frag %d)\n",2*top+1);
              exit (2);
            }
          top *= 2;
        }
      buffer[i++] = n;
      n = fgetc(sfile);
    }
  buffer[i] = '\0';
  SeqFile = fopen(buffer,"r");
  if (SeqFile == NULL)
    { fprintf(stdout,"\n\t*** Cannot open file '%s'\n",buffer);
      exit(1);
    }
  if (comments)
    { printf("# Sequence from file: %s\n",buffer);
      fflush(stdout);
    }
}

/*  Read in the input specification, honoring the old format
    if it is present.  */

int getinput()
{ int n;

  for (n = 0; n < 26; n++)
    Relement[n].type = -1;
  Sequence = -1;
  
  if (comments)
    printf("#\n");

  n = nextsymbol();
  if (feof(sfile)) exit (0);
  ungetc(n,sfile);
  if (n == '<')
    getfile();
  else
    getgrammar();

  n = nextsymbol();
  if (feof(sfile))
    { if (comments) fprintf(stdout,"#\n");
      return (1);
    }
  ungetc(n,sfile);

  if (comments)
    { printf("# Fragments:\n");   /* Get fragment parameters */
      printf("#    Number = ");
      fflush(stdout);
    }
  NumFrag = getnatural("Number of fragments");
  if (comments)
    { printf("%d, Length Range = ",NumFrag);
      fflush(stdout);
    }
  MinLen = getnatural("Min fragment length");
  MaxLen = getnatural("Max fragment length");
  if (MinLen > MaxLen)
    { fprintf(stdout,"\n\t*** Min len(%d) > Max len(%d)\n",MinLen,MaxLen);
      exit(1);
    }
  if (comments)
    { printf("[%d,%d], F/R odds = ",MinLen,MaxLen);
      fflush(stdout);
    }
  Orient = getprobability("Orientation Odds");
  if (comments)
    { printf("%.2f/%.2f\n",Orient,1.0-Orient);
      printf("# Edit Characteristics:\n");
      printf("#    Error Ramp = ");
      fflush(stdout);
    }
  BegErr = getprobability("Start error rate");
  if (comments)
    { printf("%.2f->",BegErr);
      fflush(stdout);
    }
  EndErr = getprobability("Final error rate");
  if (comments)
    { printf("%.2f, Ins/Del/Sub Odds = ",EndErr);
      fflush(stdout);
    }
  ProbI   = getprobability("Odds of inserting");
  ProbID  = ProbI + getprobability("Odds of deleting");
  if (ProbID > 1.)
    { fprintf(stdout,"\n\t*** Sum (%.2f) of Indel odds > 1\n",ProbID);
      exit(1);
    }
  if (comments)
    { printf("%.2f/%.2f/%.2f\n",ProbI,ProbID-ProbI,1.-ProbID);
      fflush(stdout);
    }

  n = nextsymbol();
  if (n == '\n')
    n = fgetc(sfile);
  ungetc(n,sfile);

  SingleFract = 1.0;
  if (feof(sfile))
    { if (comments) fprintf(stdout,"#\n");
      return (0);
    }

  if (comments)
    { printf("# Dual-End Inserts:\n");   /* Get fragment parameters */
      printf("#    Single Odds = "); /* Get dual-end parameters */
      fflush(stdout);
    }
  SingleFract = getprobability("Odds of single read");
  if (comments)
    { printf("%.2f\n#    Insert Range = ",SingleFract);
      fflush(stdout);
    }
  MinIns  = getnatural("Min insert length");
  MaxIns  = getnatural("Max insert length");
  if (MinIns > MaxIns)
    { fprintf(stdout,"\n\t*** Min insert len(%d) > Max insert len(%d)\n",
                     MinIns,MaxIns);
      exit(1);
    }
  else if (MinIns < MaxLen)
    { fprintf(stdout,"\n\t*** Min insert len(%d) < Max fragment len(%d)\n",
                     MinIns,MaxLen);
      exit(1);
    }
  if (comments)
    { printf("[%d,%d]\n#    Pairing Error Rate = ",MinIns,MaxIns);
      fflush(stdout);
    }
  PairError  = getprobability("Pairing Error Rate");
  if (comments)
    { printf("%.2f\n",PairError);
      printf("#\n");
    }

  return (0);
}

/* Read_seq gets the first FASTA formated sequence from the given input file.
   It is designed to handle lines and strings of arbitrary length.  It
   allocates memory for the sequence and returns a pointer to it.  */

#define LBUFLEN 512

char *read_seq(FILE *input, int *len)
{ static char linebuf[LBUFLEN];

  char *seqbuf, *newbuf;
  int l, e, bol, top;

  top = 1024;
  seqbuf = (char *) ckalloc(sizeof(char)*top);
  if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
  if (*linebuf == '#')
    { do
        { if (bol && *linebuf == '>')
            break;
          else if (bol)
            { if (*linebuf != '#')
                { fprintf(stderr,"No FASTA header after a comment\n");
                  exit (1);
                }
              else
                /* Put your comment processing code here */;
            }
          l = strlen(linebuf);
          bol = (linebuf[l-1] == '\n');
        }
      while (fgets(linebuf,LBUFLEN,input) != NULL);
    }

  if (*linebuf != '>')
    { fprintf(stderr,"First line must start with an >-sign\n");
      exit (1);
    }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  e   = 0;
  while (fgets(linebuf,LBUFLEN,input) != NULL)
    { if (bol && *linebuf == '>')
        break;
      l = strlen(linebuf);
      if (e + l >= top)
        { top = 1.5*(e+l) + 200;
          newbuf = (char *) ckalloc(sizeof(char)*top);
          seqbuf[e] = '\0';
          strcpy(newbuf,seqbuf);
          free(seqbuf);
          seqbuf = newbuf;
        }
      strcpy(seqbuf+e,linebuf);
      bol = (linebuf[l-1] == '\n');
      e = (e+l) - bol;
    }
  seqbuf[e] = '\0';

  if (e == 0)
    { fprintf(stderr,"Input sequence is empty\n");
      exit (1);
    }

  *len = e;
  return (seqbuf);
}


/* >>>> MAIN <<<< */

int main(int argc, char *argv[])
{ char *dna, *pname;
  int   idx, len, i;
  int   no_seed, illegal;

/* Initialize random number generator with process id. */

  pname = argv[0];

  no_seed  = 1;
  no_elem  = 1;
  comments = 1;
  illegal  = 0;
  while (argc > 1 && argv[1][0] == '-')
    { if (argv[1][1] == 's')
        { if (argc > 2)
            { Seed    = atoi(argv[2]);
              no_seed = 0;
              argc   -= 1;
              argv   += 1;
            }
        }
      else
        { for (i = 1; argv[1][i] != '\0'; i++)
            if (argv[1][i] == 'p')
              no_elem = 0;
            else if (argv[1][i] == 'F')
              comments = 0;
            else
              illegal = 1;
          if (i == 1) illegal = 1; 
        }
      argv += 1;
      argc -= 1;
    }

  if (argc != 2 || illegal)
    { fprintf(stderr,"Error: Usage: %s [-s {seed:%%d}] [-pF] specfile\n",pname);
      exit (1);
    }

  if ((sfile = fopen(argv[1],"r")) == NULL)
    { fprintf(stderr,"Error: Cannot open spec file %s\n",argv[1]);
      exit (1);
    }

  if (no_seed)
    Seed = getpid();
  srand48(Seed);

/* Setup mutation tables */

  TranChar['a'] = 'a';   FlipChar['a'] = 't';   Decode['a'] = 0;
  TranChar['c'] = 'c';   FlipChar['c'] = 'g';   Decode['c'] = 1;
  TranChar['g'] = 'g';   FlipChar['g'] = 'c';   Decode['g'] = 2;
  TranChar['t'] = 't';   FlipChar['t'] = 'a';   Decode['t'] = 3;
  TranChar['A'] = 'A';   FlipChar['A'] = 'T';   Decode['A'] = 4;
  TranChar['C'] = 'C';   FlipChar['C'] = 'G';   Decode['C'] = 5;
  TranChar['G'] = 'G';   FlipChar['G'] = 'C';   Decode['G'] = 6;
  TranChar['T'] = 'T';   FlipChar['T'] = 'A';   Decode['T'] = 7;
 
/* Read in the specification from sfile. */

  getinput();

/* Generate a DNA sequence.  */

  if (Sequence < 0)
    dna = read_seq(SeqFile,&len);
  else
    { build_seq();
      dna = Relement[Sequence].valstack->value;
      len = Relement[Sequence].valstack->length;
    }

/* Sample NumFrag samples from this sequence,
   editing each fragment after sampling,
   and generating appropriate output .
*/

#ifdef DEBUGS
  printf("dna = '%s'\n",dna);
#endif

  { int i;
    if (comments)
      { printf("# DNA Sequence:\n");
        printf("#\n");
      }
    printf("> 0\n");
    for (i = 0; i < len; i += SEGLEN)
      if (len-i < SEGLEN)
        printf("%.*s\n",len-i,dna+i);
      else
        printf("%.*s\n",SEGLEN,dna+i);
  }

  exit(0);
}
