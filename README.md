# idba_ud

This repository is a basic fork of `IDBA` (https://github.com/loneknightpy/idba),
specifically only the source files required to build `IDBA-UD`.

## changes
The source code was run through `cppcheck` (http://cppcheck.sourceforge.net/) to remove unused functions.
`cppcheck --enable=unusedFunction .`

The removed functions were:
```
assembly/assembly_utility.cpp:444:0: style: The function 'AlignReadsLocal' is never used. [unusedFunction]
assembly/assembly_utility.cpp:485:0: style: The function 'AlignReadsMultiple' is never used. [unusedFunction]
assembly/assembly_utility.cpp:42:0: style: The function 'CompareHashAlignerRecord' is never used. [unusedFunction]
assembly/assembly_utility.cpp:1016:0: style: The function 'CorrectContigs' is never used. [unusedFunction]
assembly/assembly_utility.cpp:703:0: style: The function 'CorrectPairedReads' is never used. [unusedFunction]
basic/math.cpp:36:0: style: The function 'GenerateNormalCDFTable' is never used. [unusedFunction]
misc/utils.cpp:42:0: style: The function 'IsExist' is never used. [unusedFunction]
graph/contig_info.cpp:46:0: style: The function 'ReadContigInfo' is never used. [unusedFunction]
misc/utils.cpp:79:0: style: The function 'Replace' is never used. [unusedFunction]
misc/utils.cpp:88:0: style: The function 'ToLower' is never used. [unusedFunction]
graph/contig_info.cpp:56:0: style: The function 'WriteContigInfo' is never used. [unusedFunction]
assembly/assembly_utility.cpp:229:0: style: The function 'WriteKmerFile' is never used. [unusedFunction]
idba_ud.cpp:108:0: style: The function 'contig_info_file' is never used. [unusedFunction]

basic/math.cpp:15:0: style: The function 'NormalRand' is never used. [unusedFunction]
misc/utils.cpp:79:0: style: The function 'Replace' is never used. [unusedFunction]
misc/utils.cpp:88:0: style: The function 'ToLower' is never used. [unusedFunction]
graph/contig_info.cpp:56:0: style: The function 'WriteContigInfo' is never used. [unusedFunction]
assembly/assembly_utility.cpp:79:0: style: The function 'WriteHashAlignerRecords' is never used. [unusedFunction]
```
All source code was run through `clang-format`.
