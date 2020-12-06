#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

#include <cmath>
#include <math.h>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<Eigen::Matrix<int,2,1> > springs;
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
        int s = springs.size();

        for (int i = 0; i < s; i++) {
          Eigen::Vector2i seg = springs.at(i);
          int a = seg(0);
          int b = seg(1);
          // computing the normalized direction vector
          // and length of the segment
          TV n = x.at(a) - x.at(b);
          T l = sqrt(n.dot(n));
          n /= l;

          // T test1 = l / rest_length.at(i);
          // T test2 = youngs_modulus;
          // T test4 = youngs_modulus * (l / rest_length.at(i) - 1.f);
          // TV test5 = (l / rest_length.at(i) - 1.f) * n;
          // if (test1 != test1) {
          //   std::cout << i << std::endl;
          //   std::cout << "ERROR: l/l0 is nan" << std::endl;
          //   std::cout << l << " / " << rest_length.at(i) << std::endl;
          //   std::cout <<
          // }
          // if (test2 != test2) {
          //   std::cout << i << std::endl;
          //   std::cout << "ERROR: youngs_modulus is nan" << std::endl;
          // }
          if (n != n) {
            std::cout << i << std::endl;
            std::cout << "ERROR: normal is nan: " << n(0, 0) << ", " << n(1, 0) << ", " << n(2, 0) << std::endl;
            std::cout << "x_a: " << x.at(a)(0, 0) << ", " << x.at(a)(1, 0) << ", " << x.at(a)(2, 0) << std::endl;
            std::cout << "x_b: " << x.at(b)(0, 0) << ", " << x.at(b)(1, 0) << ", " << x.at(b)(2, 0) << std::endl;
          }
          if (l != l) {
            std::cout << i << std::endl;
            std::cout << "ERROR: l is nan: " << l << std::endl;
          }
          if (isinf(l)) {
            // std::cout << a << ", " << b << std::endl;
            // std::cout << "ERROR: l is inf: " << l << std::endl;
            // std::cout << "n causing inf: " << n(0, 0) << ", " << n(1, 0) << ", " << n(2, 0) << std::endl;
            // std::cout << "n dot n: " << n.dot(n) << std::endl;
            // std::cout << "sqrt of n dot n: " << sqrt(n.dot(n)) << std::endl;
            // std::cout << "x_a: " << x.at(a)(0, 0) << ", " << x.at(a)(1, 0) << ", " << x.at(a)(2, 0) << std::endl;
            // std::cout << "x_b: " << x.at(b)(0, 0) << ", " << x.at(b)(1, 0) << ", " << x.at(b)(2, 0) << std::endl;
            l = 0.f;
          }
          if (l == 0.f) {
            // std::cout << i << std::endl;
            // std::cout << "ERROR: l is zero: " << l << std::endl;
          }
          // if (test4 != test4) {
          //   std::cout << i << std::endl;
          //   std::cout << "ERROR: youngs_modulus prod w/ i/i0 is nan" << std::endl;
          // }
          // if (test5 != test5) {
          //   std::cout << i << std::endl;
          //   std::cout << "ERROR: l/l0 prod w/ n is nan" << std::endl;
          // }

          // if (isinf(l))

          // compute the spring force with teh youngs modulus model
          f.at(a) -= youngs_modulus * (l / rest_length.at(i) - 1.f) * n;
          f.at(b) += youngs_modulus * (l / rest_length.at(i) - 1.f) * n;
        }
    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        // TODO: evaluate damping force
        int s = springs.size();

        for (int i = 0; i < s; i++) {
          Eigen::Vector2i seg = springs.at(i);
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
