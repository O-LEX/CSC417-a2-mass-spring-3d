#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

    unsigned int num_fixed = indices.size();
    unsigned int num_free = q_size - 3 * num_fixed;
    
    // スパース行列の準備
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(num_free);

    // 固定頂点のインデックスをソートしておく
    std::vector<unsigned int> sorted_indices = indices;
    std::sort(sorted_indices.begin(), sorted_indices.end());

    unsigned int free_row = 0;
    unsigned int current_fixed_index = 0;
    
    // 全頂点に対するループ（固定頂点を除くプロジェクション行列を作成）
    for (unsigned int i = 0; i < q_size / 3; ++i) {
        // 固定頂点に属さない場合に限り行を追加
        if (current_fixed_index < sorted_indices.size() && i == sorted_indices[current_fixed_index]) {
            // インデックスが固定されている場合、スキップ
            ++current_fixed_index;
        } else {
            // 固定されていない場合、自由度として行列に追加
            for (int j = 0; j < 3; ++j) {
                triplets.emplace_back(3 * free_row + j, 3 * i + j, 1.0);
            }
            ++free_row;
        }
    }

    // スパース行列 P を triplets から設定
    P.resize(3 * num_free, q_size);
    P.setFromTriplets(triplets.begin(), triplets.end());
}