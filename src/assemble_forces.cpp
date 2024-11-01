#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        
    // Initialize force vector to zero
    f.setZero(q.size());
    
    // Number of particles
    int n_particles = q.size() / 3;
    
    // Add gravity forces for each particle
    for(int i = 0; i < n_particles; i++) {
        Eigen::Vector3d f_gravity;
        Eigen::Vector3d g(0, 0, 0); // zero gravity
        dV_gravity_particle_dq(f_gravity, mass, g);
        f.segment<3>(3*i) += f_gravity;
    }
    
    // Add spring forces for each spring
    int n_springs = E.rows();
    for(int i = 0; i < n_springs; i++) {
        // Get indices of particles connected by this spring
        int idx1 = E(i, 0);
        int idx2 = E(i, 1);
        
        // Get current positions of the particles
        Eigen::Vector3d q0 = q.segment<3>(3*idx1);
        Eigen::Vector3d q1 = q.segment<3>(3*idx2);
        
        // Compute spring force
        Eigen::Vector6d f_spring;
        dV_spring_particle_particle_dq(f_spring, q0, q1, l0(i), k);
        
        // Add forces to the global force vector
        // Note: each spring affects two particles
        f.segment<3>(3*idx1) += f_spring.segment<3>(0);
        f.segment<3>(3*idx2) += f_spring.segment<3>(3);
    }
};