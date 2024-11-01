#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    V = 0.0;
    // Calculate the vector between the two particles
    Eigen::Vector3d dx = q1 - q0;
    
    // Calculate the current length of the spring
    double l = dx.norm();
    
    // Calculate the potential energy using V = (k/2)(l-l0)^2
    // where k is the spring stiffness
    // l is the current length
    // l0 is the rest length
    V += 0.5 * stiffness * std::pow(l - l0, 2);
}