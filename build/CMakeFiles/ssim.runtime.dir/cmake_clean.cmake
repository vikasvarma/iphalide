file(REMOVE_RECURSE
  "libssim.runtime.a"
  "libssim.runtime.pdb"
  "ssim.runtime.o"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/ssim.runtime.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
