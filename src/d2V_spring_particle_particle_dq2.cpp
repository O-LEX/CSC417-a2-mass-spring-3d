#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    // ばねの両端点の間のベクトルとその長さ
    Eigen::Vector3d d = q1 - q0;
    double current_length = d.norm();

    // エネルギー関数のヘッセ行列計算のための前提チェック
    if (current_length == 0) {
        H.setZero();
        return;
    }
    
    // 単位ベクトルを計算
    Eigen::Vector3d d_hat = d / current_length;

    // ヘッセ行列を構成する要素
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d K = (stiffness / current_length) * ((1 - l0 / current_length) * (I - d_hat * d_hat.transpose()) + d_hat * d_hat.transpose());

    // 6x6のヘッセ行列に反映（スプリングの両端点における力の影響を含む）
    H.block<3, 3>(0, 0) = K;
    H.block<3, 3>(0, 3) = -K;
    H.block<3, 3>(3, 0) = -K;
    H.block<3, 3>(3, 3) = K;
}