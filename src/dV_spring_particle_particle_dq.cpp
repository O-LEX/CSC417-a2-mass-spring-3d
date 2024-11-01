#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    
    // Calculate the vector between the two particles
    Eigen::Vector3d dx = q1 - q0;
    
    // Calculate the current length of the spring
    double l = dx.norm();
    
    // Compute the direction vector (normalized dx)
    Eigen::Vector3d direction = dx.normalized();
    
    // Calculate the scalar force magnitude
    // f = k(l - l0) where k is stiffness
    double force_magnitude = stiffness * (l - l0);
    
    // Calculate the force vectors for each particle
    // Force on q0 points in the positive direction
    // Force on q1 points in the negative direction
    Eigen::Vector3d f0 = force_magnitude * direction;
    Eigen::Vector3d f1 = -f0;
    
    // Assemble the 6x1 force vector
    // First 3 components for q0, last 3 components for q1
    f.segment<3>(0) = -f0;  // Negative because f = -dV/dq
    f.segment<3>(3) = -f1;
}