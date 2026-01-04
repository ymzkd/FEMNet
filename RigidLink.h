#ifndef _RIGIDLINK_
#define _RIGIDLINK_

#include "Components.h"

#ifndef SWIG
#include <vector>
#include <algorithm>

#include <Eigen/Dense>
#endif

struct RigidLink : public DOFFlags {
public:
    // 仮想的な節点なのでNodeじゃなくて座標情報だけでも良いかも
    Node *Master;
    std::vector<Node> Slaves;

    RigidLink() : DOFFlags(), Master(nullptr) {}

    RigidLink(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz)
        : DOFFlags(ux, uy, uz, rx, ry, rz), Master(nullptr) {}

    size_t
    SlaveNum()
    {
        return Slaves.size();
    }

    size_t LinkNum(){
        size_t count = 0;
        for (int i = 0; i < 6; i++) {
            if (flags[i]) count++;
        }
        return count;
    }

    Eigen::MatrixXd TransformationMatrix();

};

class RigidLinks{
public:
    std::vector<RigidLink> links;

    RigidLinks() {};
    RigidLinks(std::vector<RigidLink> links) : links(links) {};
    
    int SlaveDOFNum() {
        int count = 0;
        for (RigidLink& link : links) {
            for (Node& slave : link.Slaves) {
                for (size_t i = 0; i < 6; i++)
                    if (link.flags[i])
                        count++;
            }
        }
        return count;
    }

    std::vector<int> SlaveDOFIndices() {
        std::vector<int> indices;
        for (RigidLink& link : links) {
            for (Node& slave : link.Slaves) {
                for (size_t i = 0; i < 6; i++)
                    if (link.flags[i])
                        indices.push_back(slave.id * 6 + i);
            }
        }
        return indices;
    }

    int MasterDOFNum() {
        int count = 0;
        for (RigidLink& link : links) {
            for (size_t i = 0; i < 6; i++)
                if (link.flags[i])
                    count++;
        }
        return count;
    }

    // ???要る？
    std::vector<int> MasterDOFIndices() {
        std::vector<int> indices;
        for (RigidLink& link : links) {
            for (size_t i = 0; i < 6; i++)
                if (link.flags[i])
                    indices.push_back(link.Master->id * 6 + i);
        }
        return indices;
    }

    Eigen::MatrixXd TransformationMatrix() {
        size_t slaveDofNum = SlaveDOFNum();
        size_t masterDofNum = MasterDOFNum();
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(slaveDofNum, masterDofNum);

        size_t row = 0;
        size_t col = 0;
        for (RigidLink& link : links) {
            size_t slaveNum = link.SlaveNum();
            size_t linkNum = link.LinkNum();

            // Tlinkをブロック状に配置 (slaveNum*linkNum行 x linkNum列)
            T.block(row, col, slaveNum * linkNum, linkNum) = link.TransformationMatrix();

            row += slaveNum * linkNum;  // 次のスレーブDOFの開始行
            col += linkNum;              // 次のマスターDOFの開始列
        }
        return T;
    }
};

#endif