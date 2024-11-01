#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    
    // 剛性マトリクスの初期化
    std::vector<Eigen::Triplet<double>> triplets;
    int n = V.rows();  // 頂点の数
    K.resize(3 * n, 3 * n);
    
    // 各スプリングを処理
    for (int e = 0; e < E.rows(); ++e) {
        // スプリングの両端の頂点インデックス
        int i = E(e, 0);
        int j = E(e, 1);
        
        // 頂点 i と j の座標
        Eigen::Vector3d qi = q.segment<3>(3 * i);
        Eigen::Vector3d qj = q.segment<3>(3 * j);
        
        // スプリングの未変形長さ
        double l0_e = l0(e);
        
        // ヘッセ行列（スプリングの剛性行列）
        Eigen::Matrix<double, 6, 6> hessian;
        d2V_spring_particle_particle_dq2(hessian, qi, qj, l0_e, k);
        
        // ヘッセ行列をグローバル剛性マトリクスに加算
        for (int ii = 0; ii < 3; ++ii) {
            for (int jj = 0; jj < 3; ++jj) {
                // 頂点 i の剛性行列要素
                triplets.emplace_back(3 * i + ii, 3 * i + jj, hessian(ii, jj));
                // 頂点 j の剛性行列要素
                triplets.emplace_back(3 * j + ii, 3 * j + jj, hessian(3 + ii, 3 + jj));
                // 頂点 i と j の間のクロス項
                triplets.emplace_back(3 * i + ii, 3 * j + jj, hessian(ii, 3 + jj));
                triplets.emplace_back(3 * j + ii, 3 * i + jj, hessian(3 + ii, jj));
            }
        }
    }
    
    // スパースマトリクス K に triplets を設定
    K.setFromTriplets(triplets.begin(), triplets.end());
        
};