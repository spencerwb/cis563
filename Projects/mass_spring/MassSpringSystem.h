#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

#include <cmath>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f)
    {
        // TODO: evaluate spring force
        int s = segments.size();

        for (int i = 0; i < s; i++) {
          Eigen::Vector2i seg = segments.at(i);
          int a = seg(0);
          int b = seg(1);
          TV n = x.at(a) - x.at(b);
          n /= sqrt(n.dot(n));
          T vRel = (v.at(a) - v.at(b)).dot(n);
          f.at(a) -= youngs_modulus * (1.f / rest_length.at(i) - 1.f) * n;
          f.at(b) += youngs_modulus * (1.f / rest_length.at(i) - 1.f) * n;
        }
    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        // TODO: evaluate damping force
        int s = segments.size();

        for (int i = 0; i < s; i++) {
          Eigen::Vector2i seg = segments.at(i);
          int a = seg(0);
          int b = seg(1);
          TV n = x.at(a) - x.at(b);
          n /= sqrt(n.dot(n));
          T vRel = (v.at(a) - v.at(b)).dot(n);
          f.at(a) -= damping_coeff * vRel * n;
          f.at(b) += damping_coeff * vRel * n;
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
};
