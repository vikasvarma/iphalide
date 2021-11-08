#include "util.h"
#include "Halide.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Halide;

extern "C" int argsort(
    halide_buffer_t* data,
    halide_buffer_t* order
){
    if (data->is_bounds_query()) {
        for (int d = 0; d < order->dimensions; d++){
            data->dim[d].min = order->dim[d].min;
            data->dim[d].extent = order->dim[d].extent;
        }

    } else {
        Halide::Runtime::Buffer<double> databuf(*data);
        Halide::Runtime::Buffer<double> orderbuf(*order);

        // slice dimensions and strides:
        int W = data->dim[2].extent;
        int H = data->dim[3].extent;

        // Create runtime buffers for index and data:
        Halide::Runtime::Buffer<double> slice(W,H);
        Halide::Runtime::Buffer<double> index(W,H);

        // Sort along the first two data dimensions:
        for(int yo = 0; yo < data->dim[1].extent; yo++)
        {
            for (int xo = 0; xo < data->dim[0].extent; xo++)
            {
                // Fill index:
                index.for_each_element([&,W](int x, int y) {
                    index(x,y) = double(y * W + x);
                });

                // Fill slice data:
                slice.for_each_element([&,W](int x, int y) {
                    slice(x,y) = databuf(xo,yo,x,y);
                });

                // sort the values:
                std::stable_sort(index.begin(), index.end(), 
                    [&](size_t a, size_t b) {
                        return *(slice.begin() + a) < *(slice.begin() + b);
                    }
                );

                // Assign output:
                for (int yi = 0; yi < H; yi++)
                {
                    for (int xi = 0; xi < W; xi++)
                    {
                        orderbuf(xo,yo,xi,yi) = index(xi,yi);
                    }
                }
            }
        }
    }
    return 0;
}