#include <T_particle.h>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {

    T = 0.0;
    for(int i = 0; i < qdot.size(); i++) {
        T += 0.5 * mass * qdot(i) * qdot(i);
    }
}