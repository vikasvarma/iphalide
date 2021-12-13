#include "Halide.h"
#include "halide_image_io.h"
#include "halide_benchmark.h"
#include "bm3d.h"
#include "util.h"
#include <math.h>
#include <ctime>

using namespace Halide;
using namespace Halide::Tools;

/**
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
*/

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

    Halide::Runtime::Buffer<uint8_t> image(I);
    Halide::Runtime::Buffer<double> estimate(4,4,10,4,4);
    Halide::Runtime::Buffer<uint16_t> matches(4,4,36);

    /**
    Rdom rWin; //5x5
    RDom bw; //4x4

    scores(x,y,wx,wy) = sum(dct(b*x,b*y,bw.x, bw.y)-dct(x*b+rWin.x, y*b+rWin.y, bw.x,bw.y),'bw')
    */

    double t0 = Halide::Tools::benchmark(2, 5, [&]() {
        bm3d(image, 4, 4, 6, 10, matches, estimate);
    });

    // Print the first block set:
    for(int y = 0; y < 4; y++){
        for(int x = 0; x < 4; x++){
            printf("Matches[%d,%d]: ",x,y);
            for(int z = 0; z < 36; z++){
                printf("%d ", matches(x,y,z));
            }
            printf("\n");
        }    
    }

    printf("\n-------------------------------------------------------------\n");

    // Print the first block set:
    for(int y = 0; y < 4; y++){
        for(int x = 0; x < 4; x++){
            printf("Estimate[%d,%d]: \n",x,y);
            for(int z = 0; z < 10; z++){
                for(int by = 0; by < 4; by++){
                    printf("\t");
                    for(int bx = 0; bx < 4; bx++){
                        printf("%.2f ", estimate(x,y,z,bx,by));
                    }
                    printf(";");
                }
                printf("\n");
            }
        }    
    }

    /*
    for(int y = 0; y < 12; y++){
        for(int x = 0; x < 12; x++){
            printf("Block Coefficient[%d,%d]: ",x,y);
                for(int by = 0; by < 4; by++){
                    for(int bx = 0; bx < 4; bx++){
                        printf("%f  ", coeff(x,y,bx,by));
                    }
                    printf("; ");
                }
                printf("\n");
        }    
    }
    */

    printf("BM3D Computed in: [%f sec]\n", t0);
}

int main(){
    bm3d_();
    return 0;
}

