Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/bayesianem/myheader.h>

All headers should include `IMP/bayesianem/bayesianem_config.h` as their
first include and surround all code with `IMPBAYESIANEM_BEGIN_NAMESPACE`
and `IMPBAYESIANEM_END_NAMESPACE` to put it in the IMP::bayesianem namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
