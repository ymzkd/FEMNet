#include "RigidLink.h"

Eigen::MatrixXd RigidLink::TransformationMatrix()
{
    size_t linkNum = LinkNum();
    size_t slaveNum = SlaveNum();

    Eigen::MatrixXd TBlock = Eigen::MatrixXd::Identity(linkNum, linkNum);

    // TBlockを縦にslaveNum個並べたブロック行列を作成
    Eigen::MatrixXd result(slaveNum * linkNum, linkNum);

    for (size_t i = 0; i < slaveNum; i++) {
        Eigen::MatrixXd TB = Eigen::MatrixXd::Identity(6, 6);
        double dxi = Slaves[i].Location.x - Master.Location.x;
        double dyi = Slaves[i].Location.y - Master.Location.y;
        double dzi = Slaves[i].Location.z - Master.Location.z;
        TB(0,4) = dzi; TB(0,5) = -dyi;
        TB(1,3) = -dzi; TB(1,5) = dxi;
        TB(2,3) = dyi; TB(2,4) = -dxi;


        // flagsがtrueの成分だけを抽出
        Eigen::MatrixXd extractedBlock(linkNum, linkNum);
        size_t row_idx = 0;
        for (size_t r = 0; r < 6; r++) {
            if (flags[r]) {
                size_t col_idx = 0;
                for (size_t c = 0; c < 6; c++) {
                    if (flags[c]) {
                        extractedBlock(row_idx, col_idx) = TB(r, c);
                        col_idx++;
                    }
                }
                row_idx++;
            }
        }
        result.block(i * linkNum, 0, linkNum, linkNum) = extractedBlock;
    }

    return result;
}