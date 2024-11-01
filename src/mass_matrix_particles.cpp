#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

    // Number of particles is q.size() / 3 since each particle has x,y,z coordinates
    int n_particles = q.size() / 3;
    
    // Size of mass matrix is 3n x 3n where n is number of particles
    M.resize(3 * n_particles, 3 * n_particles);
    
    // Reserve space for the diagonal entries
    M.reserve(Eigen::VectorXi::Constant(3 * n_particles, 1));
    
    // Create triplet list for efficient sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * n_particles);  // We know exactly how many elements we'll have
    
    // Fill the diagonal entries
    for(int i = 0; i < n_particles; i++) {
        // Add mass to diagonal entries for x, y, and z coordinates of each particle
        triplets.push_back(Eigen::Triplet<double>(3*i, 3*i, mass));
        triplets.push_back(Eigen::Triplet<double>(3*i+1, 3*i+1, mass));
        triplets.push_back(Eigen::Triplet<double>(3*i+2, 3*i+2, mass));
    }
    
    // Construct the sparse matrix from triplets
    M.setFromTriplets(triplets.begin(), triplets.end());
}
