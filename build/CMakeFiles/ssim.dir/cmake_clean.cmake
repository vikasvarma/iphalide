file(REMOVE_RECURSE
  "libssim.a"
  "libssim.pdb"
  "ssim.h"
  "ssim.o"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/ssim.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
