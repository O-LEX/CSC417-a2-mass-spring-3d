#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {

    // Step 1: Compute the force and stiffness matrix at the current state
    force(tmp_force, q, qdot); // Compute force and store in tmp_force
    stiffness(tmp_stiffness, q, qdot); // Compute stiffness and store in tmp_stiffness

    // Step 2: Construct the linear system (M - dt^2 * K) * qdot_{t+1} = M * qdot_t + dt * f(q)
    Eigen::SparseMatrixd A = mass - dt * dt * tmp_stiffness;
    Eigen::VectorXd rhs = mass * qdot + dt * tmp_force;

    // Step 3: Solve for the updated velocity, qdot_{t+1}
    Eigen::SimplicialLDLT<Eigen::SparseMatrixd> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed.");
    }

    Eigen::VectorXd qdot_new = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving linear system failed.");
    }

    // Step 4: Update the generalized coordinates, q_{t+1} = q_t + dt * qdot_{t+1}
    q += dt * qdot_new;
    qdot = qdot_new;  // Update qdot for the next iteration

}