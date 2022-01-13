#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class src_t { off, simple, matching };
//</listing>
    // off - no source
    // simple -   src_dry_distros: new SD are added;                                                               src_dry_sizes: new SD are added
    // matching - src_dry_distros: find similar SD and increase their multiplicity. Add new SD if match not found; src_dry_sizes: new SD are added
  };
};
