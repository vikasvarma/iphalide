file(REMOVE_RECURSE
  "conv.runtime.o"
  "libconv.runtime.a"
  "libconv.runtime.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/conv.runtime.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
