file(REMOVE_RECURSE
  "conv_halide.runtime.o"
  "libconv_halide.runtime.a"
  "libconv_halide.runtime.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/conv_halide.runtime.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
