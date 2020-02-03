#ifndef RADIX_SORT
#define RADIX_SORT

void Set_Radix_Params(int nthread, int verbose, int *seg_beg, int *seg_len);

void *Radix_Sort(long long len, void *src, void *trg, int *bytes);

#endif // RADIX_SORT
