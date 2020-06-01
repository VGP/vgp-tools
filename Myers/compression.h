/*  Last edited: Apr 16 20:13 2020 (rd109) */
/*******************************************************************************************
 *
 *  Length limited Huffman Compressor/decompressor with special 2-bit compressor for DNA.
 *
 *  Author:   Gene Myers
 *  Date:     June 25, 2019
 *
 ********************************************************************************************/

#ifndef _COMPRESSOR
#define _COMPRESSOR

#include <stdio.h>

  //  To create a compressor, get an initially empty object with vcCreate, then
  //    add a significant corpus of the byte data to be compressed with vcAddToTable,
  //    and finally create a Huffman codec based on this corpus by calling
  //    vcCreateCodec.  The parameter "partial" should be set if not all the data
  //    to be compressed has been scanned.  At this point you have a compressor ready
  //    to operate.  You can destroy/free it with vcDestroy.

OneCodec *vcCreate();
void      vcAddToTable(OneCodec *vc, int len, char *bytes);
void      vcCreateCodec(OneCodec *vc, int partial);
void      vcDestroy(OneCodec *vc);

  //  In the instance of accumulating data over multiple threads, vcAddHistogram, will
  //    add the counts in the table for vh, to the table for vc.

void      vcAddHistogram(OneCodec *vc, OneCodec *vh);

  //  A diagnostic routine: shows you the compression scheme and if the distribution
  //    of the scanned corpus is available, it shows you that too.  Output to file 'to'.

void vcPrint(OneCodec *vc, FILE *to);

  //  You can encode and decode where ibytes/ilen are the input and the output
  //    is placed at obytes and the length of the compressed/decompressed result
  //    is returned as the value of the function.  For vcEncode, ilen is the size
  //    of the uncompressed input in bytes, and the return value is the size of
  //    the compressed output in **bits**.  The converse is true for vcDecode, i.e
  //    ilen is the number of bits in the compressed input, and the return value
  //    is the number of bytes in the uncompressed output.  The routines are endian safe.

int vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes);
int vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes);

  //  Rather than directly reading or writing an encoding of a compressor, the routines
  //    below serialize or deserialize the compressor into/outof a user-supplied buffer.
  //    vcMaxSerialSize gives the maximum size of a serialized compressor so the user
  //    can arrange a buffer of the appropriate size.  vcSerialize serializes vc into
  //    buffer 'out' and returns the # of bytes in the encoding.  vcDeserialize will reverse
  //    the process given a serialization.  The routines are endian-safe.

int       vcMaxSerialSize();
int       vcSerialize(OneCodec *vc, void *out);
OneCodec *vcDeserialize(void *in);

#endif // _COMPRESSOR
