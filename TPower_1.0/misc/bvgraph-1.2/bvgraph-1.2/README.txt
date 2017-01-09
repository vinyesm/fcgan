=The bvgraph library=

Written by David Gleich

Based on code by Sebastiano Vigna, et al. from the Webgraph Framework.
http://webgraph.dsi.unimi.it/

==To install:==

1.  Unzip bvgraph-1.1.zip
2.  Start matlab
3.  cd bvgraph-1.1
4.  addpath(pwd)

Now, bvgraph will be available.

==To test:==

1.  cd test
2.  test_main

==To compile:==

The package is distributed with a precompiled mex dll for Windows,
a precompiled mexglx file for 32-bit Linux, and a precompiled mexa64
file for 64-bit linux.

The software was compiled with 

the Microsoft Visual C++ 2003 Compiler with Matlab 7.0 on Windows
the Microsoft Visual C++ 2005 Compiler with Matlab 2007a on WinXP 64
gcc-3.3.6 with Matlab R14SP3 on on 32-bit linux
gcc-3.4.6 with Matlab 2006b on 64-bit linux.

The software was tested with

Matlab 2007b on amd64 Ubuntu 7.04
Matlab 2007a on amd64 Ubuntu 7.04
Matlab 2006b on amd64 Ubuntu 7.04 

Matlab R14SP3 on qemu (i686) Ubuntu 6.06.1
Matlab 2007b on i686 Ubunutu 5.10 

Matlab 7.0 on i386 Windows
Matlab 2007a on i386 Windows

Matlab 2007a on x86_64 WinXP64

The bvgraph library will attempt to compile itself.  This operation
will fail unless the "-ansi" flag is removed from the mex script 
on linux.

>> mex -setup
(Select gcc)
>> edit [prefdir '/mexopts.sh']
(Edit->Find and Replace, Find what: "-ansi" Replace with: " ")
(File->Save)
>> cd bvgraph-1.1
>> G = bvgraph('data/wb-cs.stanford');
