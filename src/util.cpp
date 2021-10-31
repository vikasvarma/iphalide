#include "util.h"
#include "Halide.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

extern "C" int argsort_buffer(
    halide_buffer_t *in,
    halide_buffer_t *out
){
    if (in->is_bounds_query()) {
        in->dim[0].min = out->dim[0].min;
        in->dim[0].extent = out->dim[0].extent;

    } else {

        Halide::Runtime::Buffer<float> in_buf(*in);
        Halide::Runtime::Buffer<float> out_buf(*out);

        float *in_start = (float *)in->host;
        float *out_start = (float *)out->host;

        // initialize original index locations:
        std::vector<uint16_t> index(in->dim[0].extent);
        std::iota(index.begin(), index.end(), 0);

        // sort the values:
        std::stable_sort(index.begin(), index.end(), 
            [in_start](size_t a, size_t b) {
                return *(in_start + a) < *(in_start + b);
            }
        );

        // store index and sorted value together:
        for (auto ptr = 0; ptr < in->dim[0].extent; ptr++){
            *(out_start + ptr) = index[ptr];
        }
    }
    return 0;
}

extern "C" class SearchSpace()
{
    public:
        Halide::Runtime::Buffer<uint16_t> blocks;
        Halide::Tuple window;
        int dim;

        void SearchSpace(Halide::Tuple index, int S, Halide::Tuple imsize)
        {
            // Construct the object with specified origin and search space 
            // dimensions:
            dim = S;
            int xmin (max(0,index[0]*S));
            int ymin (max(0,index[1]*S));

            window = {
                xmin, ymin,
                min(imsize[0]-1, xmin + S - 1), min(imsize[1]-1, ymin + S - 1)
            };
        }
};