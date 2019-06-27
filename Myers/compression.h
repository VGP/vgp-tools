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

typedef void VGPcompressor;  //  Compressor details are supressed

  //  There is special pre-existing compressor one should use for DNA.
  //  It compresses every base to 2-bits, where any non-ACGT letter is
  //  effectively converted to an A.  Compression is case insensitive,
  //  but decompression always delivers lower-case.

extern VGPcompressor *DNAcompressor;

  //  To create a compressor, get an initially empty object with vcCreate, then
  //    add a significant corpus of the byte data to be compressed with vcAddToTable,
  //    and finally create a Huffman code based on this corpus by calling
  //    vcCreateCodec.  The parameter "partial" should be set if not all the data
  //    to be compressed has been added.  At this point you have a compressor ready
  //    to operate.  You can destroy/free it anytime with vcDestroy.

VGPcompressor *vcCreate();
void           vcAddToTable(VGPcompressor *vc, int len, char *bytes);
void           vcCreateCodec(VGPcompressor *vc, int partial);
void           vcDestroy(VGPcompressor *vc);

  //  A useful diagnostic: shows you the compression scheme and if the distribution
  //    of the scanned corpus is available, it show you that too.

void vcPrint(VGPcompressor *vc);

  //  You can encode and decode where ibytes/ilen are the input and the output
  //    is places at obytes and the length of the compressed/decompressed result
  //    is returned as the value of the function.  For vcEncode, ilen is the size
  //    of the uncompressed input in bytes, and the return value is the size of
  //    the compressed output in **bits**.  The converse is true for vcDecode, i.e
  //    ilen is the number of bits in the compressed input, and the return value
  //    is the number of bytes in the uncompressed output.  The routines are endian safe.

int vcEncode(VGPcompressor *vc, int ilen, char *ibytes, char *obytes);
int vcDecode(VGPcompressor *vc, int ilen, char *ibytes, char *obytes);

  //  Rather than directly reading or writing an encoding of a compressor, the routines
  //    below serialize or deserialize the compressor into/outof a user-supplied buffer.
  //    vcMaxSerialSize gives the maximum size of a serialized compressor so the user
  //    can arrange a buffer of the appropriate size.  vcSerialize serializes vc into
  //    buffer out and returns the # of bytes in the encoding.  vcDeserialize will reverse
  //    the process given a buffer.  The routines are endian-safe.

int            vcMaxSerialSize();
int            vcSerialize(VGPcompressor *vc, void *out);
VGPcompressor *vcDeserialize(void *in);

#endif // _COMPRESSOR
