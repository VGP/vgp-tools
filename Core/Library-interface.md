# The VGP-tools C library interface


The interface is defined in `VGPlib.h` (currently in subdirectory Durbin/).  There are 19 functions and 8 macros, with one primary type `VgpFile` which maintains information about the current line.

As a brief synopsis, the following reads a sequence file, prints out some simple stats, and writes a binary file containing the reverse-complemented sequences.

```
{
	VgpFile *in = vgpFileOpenRead(inFile, "seq", 1);
	if (!in) die("can't open VGP sequence file %s to read", inFile);

	VgpFile *out = vgpFileOpenWriteFrom(outFile, in, FALSE, TRUE, 1);
	if (!out) die("can't open VGP sequence file %s to write", outFile);
	vgpAddProvenance(out,"revcomp","1.0","revcomp inFile outFile",0);
	vgpWriteHeader(out);
	
	int totLen = 0, totCount = 0 ;
	while (vgpReadLine(in))
	  if (in->lineType == 'S')
	    { totLen += vgpLen(in);
	      reverseComplement(vgpString(in), vgpLen(in)); // user-provided, assume acts in place
	      vgpWriteLine(out, 'S', vgpLen(in), vgpString(in));
	    }
	  else if (in->lineType == 'C')
	    totCount += vgpInt(in,0);
	printf("total sequence length %d and counts %d\n", totLen, totCount);
	printf("auto-accumulated length %d should be the same\n", in->lineInfo['S']->accum.total);
	vgpFileClose(in);
	vgpFileClose(out);
}
```
For now, more details are provided in `VGPlib.h`.

