Place the private header files in this directory. They will be
available to your code with

     #include <IMP/bayesianem/internal/myheader.h>

All headers should include `IMP/bayesianem/bayesianem_config.h` as their
first include and surround all code with `IMPBAYESIANEM_BEGIN_INTERNAL_NAMESPACE`
and `IMPBAYESIANEM_END_INTERNAL_NAMESPACE` to put it in the
IMP::bayesianem::internal namespace and manage compiler warnings.
