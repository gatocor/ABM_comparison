// Simulate randomly migrating cell
#include "../include/dtypes.cuh"
#include "../include/inits.cuh"
#include "../include/polarity.cuh"
#include "../include/solvers.cuh"
#include "../include/utils.cuh"
#include "../include/vtk.cuh"


const auto r_max = 1;
const auto n_cells = 150;
const auto n_time_steps = 150;
const auto dt = 0.05;
const auto D = 1.;


__device__ Po_cell relu_w_migration(
    Po_cell Xi, Po_cell r, float dist, int i, int j)
{
    Po_cell dF{0};
    if (i == j){
        dF.x = D*cosf(Xi.phi)*sinf(Xi.theta);
        dF.y = D*sinf(Xi.phi)*sinf(Xi.theta);
        dF.z = 0;
    };

    return dF;
}

__global__ void update_polarities(Po_cell* d_X, curandState* d_state)
{
    auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i != 0) return;

    //Generate two random normal random numbers
    float x = curand_normal(&d_state[i]);
    float y = curand_normal(&d_state[i]);

    //Generate directions in the plane
    float3 new_dir;
    new_dir.x = x;
    new_dir.y = y;
    new_dir.z = 0;
    auto new_polarity = pt_to_pol(new_dir);
    d_X[i].theta = new_polarity.theta;
    d_X[i].phi = new_polarity.phi;
}


int main(int argc, const char* argv[])
{
    // Prepare initial state
    Solution<Po_cell, Tile_solver> cells{n_cells};
    relaxed_sphere(0.75, cells);
    cells.h_X[0] = Po_cell{0};
    cells.h_X[0].phi = 0.01;
    cells.copy_to_device();
    curandState* d_state;
    cudaMalloc(&d_state, n_cells * sizeof(curandState));
    auto seed = time(NULL);
    setup_rand_states<<<(n_cells + 128 - 1) / 128, 128>>>(
        n_cells, seed, d_state);

    // Integrate cell positions
    Vtk_output output{"random_walk"};
    for (auto time_step = 0; time_step <= n_time_steps; time_step++) {
        cells.copy_to_host();
        update_polarities<<<(n_cells + 32 - 1) / 32, 32>>>(cells.d_X, d_state);
        cells.take_step<relu_w_migration>(dt);
        output.write_positions(cells);
        output.write_polarity(cells);
    }

    return 0;
}
