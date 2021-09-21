#include "Halide.h"
#include "halide_image_io.h"
#include "halide_benchmark.h"
#include "conv.h"
#include "ssim.h"

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
    printf("Score: %f\n", score(0));
    printf("Convolution benchmark with auto schedule: %gms\n", duration * 1e3);

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
    printf("Convolution benchmark with auto schedule: %gms\n", duration * 1e3);
    save_image(output,"out.png");

}

int main(){
    ssim_();
    return 0;
}

