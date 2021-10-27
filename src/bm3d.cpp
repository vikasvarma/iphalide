#include "Halide.h"
#include "util.h"
#include <string.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Halide;
using namespace std;

#define PROFILE = "REGULAR"

struct REGULAR 
{
    struct step1
    {
        int block_size = 8;
        int step_size = 3;
        int search_win = 39;
        int thr = 2500;
        int max_match = 16;
    } STEP1;

    struct step2
    {
        int block_size = 8;
        int step_size = 3;
        int search_win = 39;
        int thr = 400;
        int max_match = 32;
    } STEP2;
};

/**
 * @brief BLOCK MATCHING & 3D DENOISING
 * 
 * TODO
 * 
 */
class bm3d : public Generator<bm3d>
{
    public:
        Input<Buffer<uint8_t>> noisy{"noisy",2};
        Output<Buffer<uint16_t>> blocks{"blocks",2};
        Output<Buffer<float>> scores{"scores",1};

        // Image width and height: `W` `H`
        Expr W, H;

        void generate()
        {
            // Populate state variables for image dimensions:
            W = noisy.width();
            H = noisy.height();

            // Edges should be replicated:
            Func I = BoundaryConditions::repeat_edge(noisy);
            Func dct = block_dct2d(I,Expr(4));

            // Hard-threshold the DCT coefficients:
            dct(x,y,p,q) = select(abs(dct(x,y,p,q)) < 0, 0, dct(x,y,p,q));
            
            // Find the matches for a specific block:
            Func bscores = block_matching(
                dct, Expr(3), Expr(3), 8, 4, 10
            );

            Expr bindx = bscores(x)[0];
            Expr bindy = bscores(x)[1];
            Expr dist  = bscores(x)[2];

            blocks(x,y) = select(y == 0, bindx, bindy);
            scores(x)   = dist;
        }

        void schedule()
        {
            if (auto_schedule)
            {
                noisy.set_estimates({{0,12},{0,12}});
                blocks.set_estimates({{0,64},{0,2}});
                scores.set_estimates({{0,64}});
            } else {
                blocks.compute_root();
                scores.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, z{"z"}, p{"p"}, q{"q"};

        Func block_dct2d(Func I, Expr B)
        {
            /**
             * @brief Construct a LUT of DCT coefficients of all blocks in the 
             * image of specified block size `B` in each spatial dimension 
             * spanning `W` and `H` pixels.
             * 
             * Returns LUT(x,y,p,q) = dct<x,y>(p,q) where <x,y> is the top-left 
             * pixel coordinate of the block and <p,q> is the index of the 
             * (p,q) coefficient in the B-by-B coefficient matrix of the block.
             */
            Func LUT {"LUT"};

            // Generate DCT for each block:
            // TODO - Separate both X,Y coefficient calculation (performance).
            Expr PI (M_PI);
            RDom b(0,B,0,B);
            LUT(x,y,p,q) = Expr(2.0f/B) * sum(
                cos((PI*p*(2*b.x+1))/(Expr(2*B))) *
                cos((PI*q*(2*b.y+1))/(Expr(2*B))) *
                I(x+b.x, y+b.y)
            );

            // Divide first row and col of each coefficient block by sqrt(2):
            RDom br(0,1,0,B);
            RDom bc(0,B,0,1);
            LUT(x,y,br.x,br.y) = LUT(x,y,br.x,br.y)/Expr(sqrt(2));
            LUT(x,y,bc.x,bc.y) = LUT(x,y,bc.x,bc.y)/Expr(sqrt(2));

            return LUT;
        }

        Tuple get_search_window(
            Expr xref,
            Expr yref,
            Expr window_size,
            Expr block_size
        ){
            /** Construct search window reduction domain.
             * @brief Construct a reduction domain that spans the extent of the 
             * search window.
             *      
             *     * * * * * * * * * * xxxx
             *     * * * * * * * * * *    x
             *     * * * * * * * * * *    x
             *     * * * * * @ * * * *    x
             *     * * * * * * * * * *    x
             *     * * * * * * * * * *    x
             *     x                      x
             *     xxxxxxxxxxxxxxxxxxxxxxxx
             * 
             *  * - Top left coordainte of the search block.
             *  @ - Reference block top-left coordaintes.
             *  x - Search extents.
             */

            Expr xmin, ymin, xmax, ymax;
            xmin = xref - floor((window_size-block_size)/Expr(2.0));
            ymin = yref - floor((window_size-block_size)/Expr(2.0));
            xmax = xmin + window_size - 1;
            ymax = ymin + window_size - 1;

            Tuple window = {
                max(0, xmin), max(0, ymin),
                min(W - block_size + 1, xmax), min(H - block_size + 1, ymax)
            };

            return window;
        }

        Func block_matching(
            Func LUT,
            Expr xref,
            Expr yref,
            Expr window_size,
            Expr block_size,
            Expr max_blocks
        ){
            /** Block matching routine for BM3D
             * @brief Search the window space of <xref,yref> block to determine 
             * similar blocks.
             * 
             * Returns the top-left block coordiantes that match with the 
             * reference block.
             */

            Func similarity {"similarity"}, 
                 dist       {"distances"}, 
                 matches    {"matches"},
                 bs         {"bs"};
            
            Tuple limits = get_search_window(xref,yref,window_size,block_size);
            similarity(x,y,z) = 
                select(
                    ((limits[0] + x % window_size) < limits[2]) &&
                    ((limits[1] + x / window_size) < limits[3]),
                    abs(
                        LUT(xref,yref,y,z) - 
                        LUT(
                            cast<int>(limits[0] + x % window_size),
                            cast<int>(limits[1] + x / window_size), 
                            y,z
                        )
                    ),
                    cast<double>(INFINITY)
                );
                    
            // Compute distances:
            RDom b (0,block_size,0,block_size,"search_block");
            dist(x) = sum(pow(similarity(x,b.x,b.y),2)) / 
                      pow(window_size,2);
            dist.trace_stores();

            // sorted_dist calls the external C function to return the sorted 
            // argmin order:
            std::vector<ExternFuncArgument> args;
            args.push_back(dist);
            matches.define_extern("argsort_buffer", args, Float(32), 1);

            bs(x) = Tuple(
                cast<uint16_t>(limits[0] + matches(2*x) % window_size),
                cast<uint16_t>(limits[1] + matches(2*x) / window_size),
                cast<float>(matches(2*x + 1))
            );

            return bs;
        }

}; // END - BM3D

HALIDE_REGISTER_GENERATOR(bm3d, bm3d_generator)