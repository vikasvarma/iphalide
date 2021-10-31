#include "Halide.h"
#include "halide_image_io.h"
#include "halide_benchmark.h"
#include "conv.h"
#include "ssim.h"
#include "dct2d.h"
#include "idct2d.h"
#include "blockdct.h"
#include "searchblocks.h"
#include "argsort.h"
#include <math.h>

using namespace Halide;
using namespace Halide::Tools;

void ssim_(){

    // Load an image:
    Halide::Runtime::Buffer<uint8_t> a = load_image("../data/truffle.png");
    Halide::Runtime::Buffer<uint8_t> b = load_image("../data/truffle2.png");
    Halide::Runtime::Buffer<float> score(1);

    double duration = Halide::Tools::benchmark(2, 5, [&]() {
        ssim(a, b, score);
    });
    printf("SSIM: %f Benchmark: [%gms]\n", score(0), duration * 1e3);

}

void conv_(){

    // Load an image
    float h[3][3] = {{0.1,0,0},{0,1,0},{0,0,0.01}};
    Halide::Runtime::Buffer<uint8_t> image = load_image("../data/gray.png");
    Halide::Runtime::Buffer<uint8_t> output(image.width(), image.height());
    Halide::Runtime::Buffer<float> kernel(h);

    double duration = Halide::Tools::benchmark(2, 5, [&]() {
        conv(image, kernel, output);
    });
    printf("conv: 'out.png' Benchmark: [%gms]\n", duration * 1e3);
}

void dct2d_()
{
    uint8_t I[5][4] = {
        {10,20,30,20},
        {23,34,45,56},
        {32,43,54,65},
        {11,22,33,44},
        {0,0,1,0}
    };

    Halide::Runtime::Buffer<uint8_t> img(I);
    Halide::Runtime::Buffer<double> coeff(img.width(),img.height());
    Halide::Runtime::Buffer<uint8_t> ret(img.width(),img.height());

    double d1 = Halide::Tools::benchmark(2, 5, [&]() {
        dct2d(img, coeff);
    });
    printf("2D DCT Benchmark: [%gms]\n", d1 * 1e3);

    double d2 = Halide::Tools::benchmark(2, 5, [&]() {
        idct2d(coeff, ret);
    });
    printf("2D IDCT Benchmark: [%gms]\n", d2 * 1e3);
}

void bm3d_()
{
    uint8_t I[12][12] = {
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1},
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1},
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1},
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1},
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1},
        {1,2,3,4,5,6,7,8,9,10,11,12},
        {12,11,10,9,8,7,6,5,4,3,2,1}
    };

    Halide::Runtime::Buffer<uint8_t> noisy(I);
    Halide::Runtime::Buffer<double> dct(12,12,4,4);
    Halide::Runtime::Buffer<double> sim(2,2,8,8);
    Halide::Runtime::Buffer<uint8_t> sortstatus(2,2);
    double t0 = Halide::Tools::benchmark(2, 5, [&]() {
        blockdct(noisy, 4, dct);
    });

    double t1 = Halide::Tools::benchmark(2, 5, [&]() {
        searchblocks(dct, 4, 8, sim);
    });

    double t2 = Halide::Tools::benchmark(2, 5, [&]() {
        argsort(sim, sortstatus);
    });

    // Check results:
    double exp;
    for(int sx = 0; sx < 1; sx ++){
        for(int sy = 0; sy < 1; sy ++){
            int xmin = sx * 8;
            int ymin = sy * 8;
            int xmax = xmin + 8;
            int ymax = ymin + 8;
            int xref = xmin + 2;
            int yref = ymin + 2;

            for (int x = xmin; x < xmax; x++){
                for (int y = ymin; y < ymax; y++){
                    exp = 0;
                    for (int bx = 0; bx < 4; bx++){
                        for (int by = 0; by < 4; by++){
                            exp += pow(dct(x,y,bx,by) - dct(xref,yref,bx,by),2);
                        }
                    }
                    exp /= pow(4,2);

                    if (exp != sim(sx,sy,x,y))
                    {
                        printf("Error(%d,%d,%d,%d): act[%f] exp[%f]\n", 
                               sx,sy,x,y, sim(sx,sy,x,y), exp);
                    }
                }    
            }
        }
    }

    printf("LUT Computed in: [%f msec]\n", t0);
    printf(" Search done in: [%f msec]\n", t1);
    printf("Sorting done in: [%f msec]\n", t2);
}

int main(){
    bm3d_();
    return 0;
}

