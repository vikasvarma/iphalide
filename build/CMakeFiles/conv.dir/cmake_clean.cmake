file(REMOVE_RECURSE
  "conv.h"
  "conv.o"
  "libconv.a"
  "libconv.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/conv.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
