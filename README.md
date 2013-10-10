BLASR Library
=============

This repository contains an (incomplete) refactor of core BLASR code that  
compiles into a library.  The idea being that other programs needing access   
to BLASR-like functionality can build it right-in instead of forking out   
to a 'blasr' executable.  This has the advantages of enabling tighter   
integration points within your program flow, eliminates the need for  
HDF libraries to be installed, among others.

This first release contains all the basics to enable a sparse dynamic-  
programming alignment between two sequences.  The plan is to eventually fold  
this back into the main BLASR code base, clean up the interface and document  
it more thoroughly.
