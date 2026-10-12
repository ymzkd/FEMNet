// LinePolyLoad: 閉形式・独立した数値積分・荷重組立・後からの設定のテスト。
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Model.h"
#include "LoadComponent.h"

namespace {

int failures = 0;

void Check(bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cout << "  [FAIL] " << message << "\n";
        ++failures;
    }
}

void CheckNear(double actual, double expected, const std::string& message)
{
    Check(std::isfinite(actual) &&
          std::abs(actual - expected) <= 1e-12 * std::max(1.0, std::abs(expected)), message);
}

void CheckVector(const Eigen::Vector3d& actual, const Eigen::Vector3d& expected,
                 const std::string& message)
{
    for (int k = 0; k < 3; ++k)
        CheckNear(actual(k), expected(k), message + " component " + std::to_string(k));
}

Eigen::Vector3d Force(const NodeLoadData& load)
{
    return Eigen::Vector3d(load.Px(), load.Py(), load.Pz());
}

void CheckDistribution(const std::vector<double>& w, const std::vector<double>& params,
                       const Vector& direction, double fi, double fj, const std::string& name)
{
    Node ni(4, 1.0, 2.0, 3.0), nj(7, 11.0, 2.0, 3.0); // length = 10
    LinePolyLoad load(w, params, &ni, &nj, direction);
    const auto loads = load.NodeLoads();
    const double norm = std::hypot(direction.x, direction.y, direction.z);
    const Eigen::Vector3d unit(direction.x / norm, direction.y / norm, direction.z / norm);
    Check(loads.size() == 2, name + " returns two node loads");
    Check(loads[0].id == ni.id && loads[1].id == nj.id, name + " node IDs");
    CheckVector(Force(loads[0]), fi * unit, name + " i force");
    CheckVector(Force(loads[1]), fj * unit, name + " j force");
    for (const auto& nl : loads)
    {
        CheckNear(nl.Mx(), 0.0, name + " Mx");
        CheckNear(nl.My(), 0.0, name + " My");
        CheckNear(nl.Mz(), 0.0, name + " Mz");
    }
}

void TestClosedForms()
{
    CheckDistribution({2, 2}, {0, 1}, Vector(0, 0, -5), 10, 10, "uniform");
    CheckDistribution({0, 6}, {0, 1}, Vector(0, 2, 0), 10, 20, "triangular");
    CheckDistribution({2, 8}, {0.2, 0.6}, Vector(2, -3, 6), 11.2, 8.8, "partial trapezoid");
    CheckDistribution({1, 1}, {0.2, 0.4}, Vector::XAxis(), 1.4, 0.6, "partial uniform");
    CheckDistribution({-2, 2}, {0, 1}, Vector::YAxis(), -10.0 / 3, 10.0 / 3,
                      "zero resultant with nonzero moment");
    CheckDistribution({0, 0, 0}, {0, 0.5, 1}, Vector::ZAxis(), 0, 0, "zero loads");
    CheckDistribution({2, -5}, {0.5, 0.5}, Vector::XAxis(), 0, 0, "zero width");
    CheckDistribution({2, 2}, {0, 1}, Vector(0, 0, -1e-100), 10, 10, "tiny direction");
    CheckDistribution({2, 2}, {0, 1}, Vector(0, 0, -1e100), 10, 10, "large direction");
}

void TestPolygonAndSteps()
{
    Node ni(0, 1, 2, 3), nj(1, 4, 6, 15); // oblique line, length = 13
    const Eigen::Vector3d unit(2.0 / 7, -3.0 / 7, 6.0 / 7);
    const Eigen::Vector3d origin(-5, 8, 2);
    const Eigen::Vector3d ui(-0.3, 0.1, 0.7), uj(0.2, -0.8, 0.4);
    const std::vector<std::vector<double>> weights = {
        {0, 4, 2, 0},
        {2, 2, -4, -4, 0, 0, 3, 3}};
    const std::vector<std::vector<double>> positions = {
        {0.1, 0.3, 0.7, 0.9},
        {0.1, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 0.9}};

    for (size_t c = 0; c < weights.size(); ++c)
    {
        const auto& w = weights[c];
        const auto& params = positions[c];
        LinePolyLoad load(w, params, &ni, &nj, Vector(2, -3, 6));
        const auto loads = load.NodeLoads();
        Eigen::Vector3d fi = Eigen::Vector3d::Zero(), fj = Eigen::Vector3d::Zero();
        Eigen::Vector3d total = Eigen::Vector3d::Zero(), moment = Eigen::Vector3d::Zero();
        double work = 0;

        // Simpson積分: 一次形状関数×一次荷重を、端点と中点で独立に評価する。
        for (size_t k = 0; k + 1 < w.size(); ++k)
        {
            if (params[k + 1] == params[k])
                continue;
            for (int p = 0; p < 3; ++p)
            {
                const double u = p / 2.0;
                const double t = (1 - u) * params[k] + u * params[k + 1];
                const double weight = 13 * (params[k + 1] - params[k]) * (p == 1 ? 4 : 1) / 6;
                const Eigen::Vector3d q = ((1 - u) * w[k] + u * w[k + 1]) * unit;
                const Eigen::Vector3d x = (1 - t) * Eigen::Vector3d(1, 2, 3) +
                                         t * Eigen::Vector3d(4, 6, 15);
                fi += weight * (1 - t) * q;
                fj += weight * t * q;
                total += weight * q;
                moment += weight * (x - origin).cross(q);
                work += weight * q.dot((1 - t) * ui + t * uj);
            }
        }

        CheckVector(Force(loads[0]), fi, "polygon i force");
        CheckVector(Force(loads[1]), fj, "polygon j force");
        CheckVector(Force(loads[0]) + Force(loads[1]), total, "resultant conservation");
        CheckVector((Eigen::Vector3d(1, 2, 3) - origin).cross(Force(loads[0])) +
                    (Eigen::Vector3d(4, 6, 15) - origin).cross(Force(loads[1])),
                    moment, "moment conservation about arbitrary origin");
        CheckNear(Force(loads[0]).dot(ui) + Force(loads[1]).dot(uj), work,
                  "virtual work for linear displacement");

        std::vector<double> reversed_w(w.rbegin(), w.rend());
        std::vector<double> reversed_params;
        for (auto it = params.rbegin(); it != params.rend(); ++it)
            reversed_params.push_back(1 - *it);
        LinePolyLoad reversed(reversed_w, reversed_params, &nj, &ni, Vector(2, -3, 6));
        const auto reverse_loads = reversed.NodeLoads();
        CheckVector(Force(reverse_loads[0]), Force(loads[1]), "reverse endpoints j");
        CheckVector(Force(reverse_loads[1]), Force(loads[0]), "reverse endpoints i");
    }
}

void TestScalingAndChanges()
{
    Node ni(0, 0, 0, 0), nj(1, 10, 0, 0);
    LinePolyLoad load({2, 8}, {0.2, 0.6}, &ni, &nj, Vector(1, 2, 2));
    const auto original = load.NodeLoads();
    for (double factor : {-2.5, 0.0, 3.0})
    {
        const auto copy = std::dynamic_pointer_cast<LinePolyLoad>(load.scaled(factor));
        Check(copy != nullptr && copy.get() != &load, "scaled returns a new LinePolyLoad");
        Check(copy->node_i == &ni && copy->node_j == &nj && copy->params == load.params,
              "scaled preserves nodes and positions");
        CheckVector(copy->direction.toEigen(), load.direction.toEigen(), "scaled preserves direction");
        const auto scaled = copy->NodeLoads();
        CheckVector(Force(scaled[0]), factor * Force(original[0]), "scaled i force");
        CheckVector(Force(scaled[1]), factor * Force(original[1]), "scaled j force");
    }
    Check(load.w == std::vector<double>({2, 8}), "scaled does not modify original weights");

    LinePolyLoad refined({2, 5, 8}, {0.2, 0.4, 0.6}, &ni, &nj, Vector(1, 2, 2));
    CheckVector(Force(refined.NodeLoads()[0]), Force(original[0]), "subdivision i");
    CheckVector(Force(refined.NodeLoads()[1]), Force(original[1]), "subdivision j");

    nj.Location.x = 20;
    CheckVector(Force(load.NodeLoads()[0]), 2 * Force(original[0]), "current node positions");
    load.w = {1, 1};
    load.params = {0, 1};
    load.direction = Vector(0, 0, -4);
    CheckVector(Force(load.NodeLoads()[0]), Eigen::Vector3d(0, 0, -10), "current public fields");
}

void TestAssembly()
{
    FEModel model;
    model.Nodes = {Node(0, 0, 0, 0), Node(1, 5, 0, 0), Node(2, 10, 0, 0)};
    // 要素を必要とせず、指定した2節点だけへ載荷する。
    const std::vector<std::shared_ptr<LoadBase>> loads = {
        std::make_shared<LinePolyLoad>(std::vector<double>{2, 2}, std::vector<double>{0, 1},
                                       &model.Nodes[0], &model.Nodes[2], Vector(2, 0, 0)),
        std::make_shared<NodeLoad>(2, 1, 2, 3)};
    const Eigen::VectorXd assembled = model.AssembleLoadVector(loads);
    Check(assembled.size() == 18, "assembly size");
    for (int k = 0; k < assembled.size(); ++k)
    {
        const double expected = k == 0 ? 10 : k == 12 ? 11 : k == 13 ? 2 : k == 14 ? 3 : 0;
        CheckNear(assembled(k), expected, "assembly DOF " + std::to_string(k));
    }
}

void TestDeferredSetup()
{
    Node ni(0, 0, 0, 0), nj(1, 10, 0, 0);
    LinePolyLoad load({}, {}, nullptr, nullptr, Vector());
    load.node_i = &ni;
    load.node_j = &nj;
    load.w = {2, 2};
    load.params = {0, 1};
    load.direction = Vector(0, -3, 0);
    const auto loads = load.NodeLoads();
    CheckVector(Force(loads[0]), Eigen::Vector3d(0, -10, 0), "deferred setup i");
    CheckVector(Force(loads[1]), Eigen::Vector3d(0, -10, 0), "deferred setup j");
}

} // namespace

int main()
{
    std::cout << "LinePolyLoad test\n";
    TestClosedForms();
    TestPolygonAndSteps();
    TestScalingAndChanges();
    TestAssembly();
    TestDeferredSetup();
    std::cout << (failures == 0 ? "RESULT: PASS" : "RESULT: FAIL")
              << " (" << failures << " failures)\n";
    return failures == 0 ? 0 : 1;
}
