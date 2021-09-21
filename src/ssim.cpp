#include "Halide.h"
#include <math.h>
#include <stdio.h>

#define RAD 5
#define SIGMA 1.5 

using namespace Halide;
using namespace std;

class ssim : public Generator<ssim>
{
    public:
        Input<Buffer<uint8_t>> A{"A",2};
        Input<Buffer<uint8_t>> B{"B",2};
        Output<Buffer<float>> score{"score"};

        /**
         * @brief SSIM Algorithm
         * 
         */
        void generate()
        {
            // Construct separable gaussian kernel:
            Buffer<float> H(2*RAD+1);
            float sum_ = 0.0;
            for (int n = -RAD; n <= RAD ; n++) {
                H(n+RAD) = exp(-((n*n)/(SIGMA*SIGMA))/2);
                sum_ += H(n+RAD);
            }
            for (int n=0 ; n<2*RAD+1 ; n++) {
                H(n) /= sum_;
            }

            RDom px(A);

            // Cast and clamp images:
            src1 = BoundaryConditions::repeat_edge(A);
            src2 = BoundaryConditions::repeat_edge(B);
            a(x,y) = cast<float>(src1(x,y));
            b(x,y) = cast<float>(src2(x,y));

            // Separable filtering based mean calculations:
            RDom k(H);
            mu_ax(x,y) = sum(H(k) * a(x+RAD-k,y));
            mu_bx(x,y) = sum(H(k) * b(x+RAD-k,y));
            mu_a(x,y)  = sum(H(k) * mu_ax(x,y+RAD-k));
            mu_b(x,y)  = sum(H(k) * mu_bx(x,y+RAD-k));

            // Mean squares:
            mu_aa(x,y) = mu_a(x,y) * mu_a(x,y);
            mu_bb(x,y) = mu_b(x,y) * mu_b(x,y);
            mu_ab(x,y) = mu_a(x,y) * mu_b(x,y);

            // Mean corrected images:
            mc_a(x,y)  = a(x,y) - mu_a(x, y);
            mc_b(x,y)  = b(x,y) - mu_b(x, y);
            mc_aa(x,y) = mc_a(x,y) * mc_a(x,y);
            mc_bb(x,y) = mc_b(x,y) * mc_b(x,y);
            mc_ab(x,y) = mc_a(x,y) * mc_b(x,y);

            // Variance computation:
            sig_ax(x,y)  = sum(H(k) * mc_aa(x+RAD-k,y));
            sig_bx(x,y)  = sum(H(k) * mc_bb(x+RAD-k,y));
            sig_abx(x,y) = sum(H(k) * mc_ab(x+RAD-k,y));
            sig_a(x,y)   = sum(H(k) * sig_ax(x,y+RAD-k));
            sig_b(x,y)   = sum(H(k) * sig_bx(x,y+RAD-k));
            sig_ab(x,y)  = sum(H(k) * sig_abx(x,y+RAD-k));

            score_map(x,y) = min(
                ((2 * mu_ab(x,y) + C1) * (2 * sig_ab(x,y) + C2)) /
                ((mu_aa(x,y)+mu_bb(x,y)+C1) * (sig_a(x,y)+sig_b(x,y)+C2)),
                1
            );
            score(index) = sum(score_map(px.x, px.y))/(A.width()*A.height());
        }

        void schedule()
        {
            if (auto_schedule)
            {   
                A.set_estimates({{0,4096},{0,3072}});
                B.set_estimates({{0,4096},{0,3072}});
                score.set_estimates({{0,0}});
            } else 
            {
                score.compute_root();
            }
        }

    private:
        float C1 = pow(0.01*255,2);
        float C2 = pow(0.03*255,2);

        Var x{"x"}, y{"y"}, index{"index"};

        Func src1, src2;
        Func a,b;
        Func mu_ax, mu_bx, mu_a, mu_b;
        Func mu_aa, mu_bb, mu_ab;
        Func mc_a, mc_b, mc_aa, mc_bb, mc_ab;
        Func sig_ax, sig_bx, sig_abx, sig_a, sig_b, sig_ab;
        Func score_map;
};

HALIDE_REGISTER_GENERATOR(ssim, ssim_generator)