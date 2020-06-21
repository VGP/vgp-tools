/*******************************************************************************************
 *
 *  Filter Expression Parser & Evaluator
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 31, 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#undef PRINT_TREE

#include "gene_core.h"
#include "pb_expr.h"

#define OP_OR  0
#define OP_AND 1
#define OP_NOT 2
#define OP_LT  3
#define OP_LE  4
#define OP_GT  5
#define OP_GE  6
#define OP_NE  7
#define OP_EQ  8
#define OP_INT 9
#define OP_ZM  10
#define OP_LN  11
#define OP_RQ  12
#define OP_BC1 13
#define OP_BC2 14
#define OP_BQ  15
#define OP_NP  16
#define OP_QS  17
#define OP_QE  18

#ifdef PRINT_TREE

static char *Symbol[] =
  { "OR", "AND", "NOT", "LT", "LE", "GT", "GE", "NE", "EQ", "INT",
    "ZM", "LN", "RQ", "BC1", "BC2", "BQ", "NP", "QS", "QE" };

#endif

static char *Error_Messages[] =
  { "Out of memory",
    "Unrecognized token",
    "Expecting closing paren",
    "Expecting comparison operator"
  };

static char      *Scan;
static int        Error;
static int        Contents;

#define ERROR(msg)	\
{ Error = msg;		\
  return (NULL);	\
}

typedef struct _node
  { int           op;
    struct _node *lft;
    struct _node *rgt;
  } Node;

static Node *node(int op, Node *lft, Node *rgt)
{ Node *v;

  v = (Node *) malloc(sizeof(Node));
  if (v == NULL)
    ERROR(0);
  v->op  = op;
  v->lft = lft;
  v->rgt = rgt;
  return (v);
}

static Node *terminal()
{ int   op;
  int64 x;

  switch (*Scan)
  { case 'z':
      if (Scan[1] != 'm')
        ERROR(1);
      op = OP_ZM;
      Contents |= HAS_ZM;
      Scan += 2;
      break;
    case 'l':
      if (Scan[1] != 'n')
        ERROR(1);
      op = OP_LN;
      Scan += 2;
      break;
    case 'r':
      if (Scan[1] != 'q')
        ERROR(1);
      op = OP_RQ;
      Contents |= HAS_RQ;
      Scan += 2;
      break;
    case 'b':
      if (Scan[1] == 'c')
        { if (Scan[2] == '1')
            op = OP_BC1;
          else if (Scan[2] == '2')
            op = OP_BC2;
          else
            ERROR(1);
          Contents |= HAS_BC;
          Scan += 3;
        }
      else if (Scan[1] == 'q')
        { op = OP_BQ;
          Contents |= HAS_BQ;
          Scan += 2;
        }
      else
        ERROR(1);
      break;
    case 'n':
      if (Scan[1] != 'p')
        ERROR(1);
      op = OP_NP;
      Contents |= HAS_NP;
      Scan += 2;
      break;
    case 'q':
      if (Scan[1] == 's')
        { Contents |= HAS_QS;
          op = OP_QS;
        }
      else if (Scan[1] == 'e')
        { Contents |= HAS_QE;
          op = OP_QE;
        }
      else
        ERROR(1);
      Scan += 2;
      break;
    default:
      if (!isdigit(*Scan))
        ERROR(1);
      x = *Scan++-'0';
      while (isdigit(*Scan))
        x = 10*x + (*Scan++ - '0');
      return (node(OP_INT,(Node *) x,NULL));
  }
  return (node(op,NULL,NULL));
}

static Node *or();

static Node *pred()
{ Node *v;

  while (isspace(*Scan))
    Scan += 1;
  if (*Scan == '(')
    { Scan += 1;
      v = or();
      if (v == NULL)
        return (NULL);
      while (isspace(*Scan))
        Scan += 1;
      if (*Scan != ')')
        ERROR(2);
      Scan += 1;
      return (v);
    }

  { Node *w;
    int   op;

    v = terminal();
    if (v == NULL)
      return (NULL);

    while (isspace(*Scan))
      Scan += 1;
    if (*Scan == '<')
      { if (Scan[1] == '=')
          { Scan += 2;
            op = OP_LE;
          }
        else
          { Scan += 1;
            op = OP_LT;
          }
      }
    else if (*Scan == '>')
      { if (Scan[1] == '=')
          { Scan += 2;
            op = OP_GE;
          }
        else
          { Scan += 1;
            op = OP_GT;
          }
      }
    else if (*Scan == '!')
      { if (Scan[1] != '=')
          ERROR(3);
        Scan += 2;
        op = OP_NE;
      }
    else if (*Scan == '=')
      { if (Scan[1] != '=')
          ERROR(3);
        Scan += 2;
        op = OP_EQ;
      }
    else
      ERROR(3);

    while (isspace(*Scan))
      Scan += 1;
    w = terminal();
    if (w == NULL)
      return (NULL);

    return (node(op,v,w));
  }
}

static Node *and()
{ Node *v, *w; 

  v = pred();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan != '&')
        return (v);
      if (Scan[1] != '&')
        ERROR(1);
      Scan += 2;
      w = pred();
      if (w == NULL)
        return (NULL);
      v = node(OP_AND,v,w);
    }
}

static Node *or()
{ Node *v, *w; 

  v = and();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan != '|')
        return (v);
      if (Scan[1] != '|')
        ERROR(1);
      Scan += 2;
      w = and();
      if (w == NULL)
        return (NULL);
      v = node(OP_OR,v,w);
    }
}

#ifdef PRINT_TREE

static void print_tree(Node *v, int level)
{ if (v->op == OP_NOT)
    { printf("%*s%s\n",level,"",Symbol[v->op]);
      print_tree(v->lft,level+2);
    }
  else if (v->op <= OP_EQ)
    { printf("%*s%s\n",level,"",Symbol[v->op]);
      print_tree(v->lft,level+2);
      print_tree(v->rgt,level+2);
    }
  else if (v->op == OP_INT)
    printf("%*s%s %d\n",level,"",Symbol[v->op],(int) (v->lft));
  else
    printf("%*s%s\n",level,"",Symbol[v->op]);
}

#endif

Filter *parse_filter(char *expr, int *need)
{ Node *v;

  Contents = 0;
  Scan  = expr;
  v     = or();
  if (v == NULL)
    { if (Error == 0)
        fprintf(stderr,"%s: Out of memory parsing filter expression\n",Prog_Name);
      else
        { fprintf(stderr,"%s: Filter expression syntax error:\n\n",Prog_Name);
          fprintf(stderr,"    %s\n",expr);
          fprintf(stderr,"%*s^ %s\n",(int) ((Scan-expr)+4),"",Error_Messages[Error]);
        }
    }

  *need = Contents;
  return ((Filter *) v);
}

static int eval_S(Node *v, samRecord *s)
{ switch (v->op)
  { case OP_OR:
      return (eval_S(v->lft,s) || eval_S(v->rgt,s));
    case OP_AND:
      return (eval_S(v->lft,s) && eval_S(v->rgt,s));
    case OP_NOT:
      return ( ! eval_S(v->lft,s));
    case OP_LT:
      return (eval_S(v->lft,s) < eval_S(v->rgt,s));
    case OP_LE:
      return (eval_S(v->lft,s) <= eval_S(v->rgt,s));
    case OP_GT:
      return (eval_S(v->lft,s) > eval_S(v->rgt,s));
    case OP_GE:
      return (eval_S(v->lft,s) >= eval_S(v->rgt,s));
    case OP_NE:
      return (eval_S(v->lft,s) != eval_S(v->rgt,s));
    case OP_EQ:
      return (eval_S(v->lft,s) == eval_S(v->rgt,s));
    case OP_INT:
      return ((int) (int64) (v->lft));
    case OP_ZM:
      return (s->well);
    case OP_LN:
      return (s->len);
    case OP_RQ:
      return ((int) (1000*s->qual));
    case OP_BC1:
      return (s->bc[0]);
    case OP_BC2:
      return (s->bc[1]);
    case OP_BQ:
      return (s->bqual);
    case OP_NP:
      return (s->nump);
    case OP_QS:
      return (s->beg);
    case OP_QE:
      return (s->end);
  }
  return (0);
}

int evaluate_bam_filter(Filter *v, samRecord *s)
{ return (eval_S((Node *) v,s));
}
