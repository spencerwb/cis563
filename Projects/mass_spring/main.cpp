#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unordered_set>

#include "SimulationDriver.h"

#include <cmath>
#include <iterator>
#include <unordered_map>

// type definitions
using T = float;
constexpr int dim = 3;
using TV = Eigen::Matrix<T,dim,1>;

void obj(std::string filename, std::vector<TV>& X, std::vector<Eigen::Vector4i>& F) {
  std::ofstream fs;
  fs.open(filename);
  int count = 0;
  for (TV x : X) {
      fs << "v";
      for (int i = 0; i < 3; i++)
          fs << " " << x(i);
      fs << "\n";
      count++;
  }
  for (Eigen::Vector4i f : F) {
    f += Eigen::Vector4i(1, 1, 1, 1);
    fs << "f";
    for (int i = 0; i < 4; i++)
        fs << " " << f(i);
    fs << "\n";
  }
  fs.close();
  return;
}

int main(int argc, char* argv[])
{
    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 0.f;
    T damping_coeff = 0.f;
    T dt = 0;

    // node data
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    // segment data
    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<T> rest_length;

    // face data
    std::vector<Eigen::Vector4i> faces;

    if (argc < 2)
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "0") == 0) // cloth case
    {
        // TODO
        /*
            1. Create node data: position, mass, velocity
            2. Fill segments and rest_length, including struct springs, shearing springs and bending springs.
            3. Choose proper youngs_modulus, damping_coeff and dt.
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
            5. Generate quad mesh for rendering.
        */

    		int xN = 10;
    		int yN = 10;
    		int N = xN * yN;
    		T xW = 1.f;
    		T yW = 1.f;
    		T xWHalf = xW / 2.f;
    		T yWHalf = yW / 2.f;
    		T mN = 100.f / N;

        youngs_modulus = 100.f;
        damping_coeff = 100.f;
        dt = 0.0001f;

        // rest displacement will be the starting length of
        // each segment
        T rLStx = xW / float(xN);
        T rLSty = yW / float(yN);
        T rLSh = sqrt(rLStx * rLStx + rLSty * rLSty);
        T rLBdx = rLStx * 2.f;
        T rLBdy = rLSty * 2.f;

    		for (int j = 0; j < yN; j++) {
    			for (int i = 0; i < xN; i++) {
    				x.push_back(TV(xW/xN * i - xWHalf, yW/yN * j - yWHalf, 0.f));
    				v.push_back(TV(0.f, 0.f, 0.f));
    				m.push_back(mN);
            std::cout << x.back()(0,0) << " " << x.back()(1,0) << " " << x.back()(2,0) << std::endl;
            if ((i == 0 && j == yN - 1) || (i == xN - 1 && j == yN - 1)) {
              node_is_fixed.push_back(true);
            } else {
              node_is_fixed.push_back(false);
            }

            int idx = i + j * xN;

            if (j < yN - 1) {
              segments.push_back(Eigen::Vector2i(idx, idx + xN));
              rest_length.push_back(rLSty);
            }
            if (i < xN - 1) {
              segments.push_back(Eigen::Vector2i(idx, idx + 1));
              rest_length.push_back(rLStx);
            }

            if (i < xN - 1 && j < yN - 1) {
              segments.push_back(Eigen::Vector2i(idx, idx + 1 + xN));
              segments.push_back(Eigen::Vector2i(idx + 1, idx + xN));
              rest_length.push_back(rLSh);
              rest_length.push_back(rLSh);
            }

            if (j < yN - 2) {
              segments.push_back(Eigen::Vector2i(idx, idx + 2 * xN));
              rest_length.push_back(rLBdy);
            }
            if (i < xN - 2) {
              segments.push_back(Eigen::Vector2i(idx, idx + 2));
              rest_length.push_back(rLBdx);
            }

            if (i > 0 && j > 0) {
              // the current point should be the bottom right corner
              // of its face. this vertex will be responsible for
              // constructing the face which will eventually be written to an
              // obj file
              faces.push_back(Eigen::Vector4i(idx, idx - 1, idx - 1 - xN, idx - xN));
            }

            obj("cloth.obj", x, faces);

    			}
    		}


        driver.helper = [&](T t, T dt) {
            // TODO
            int n = driver.ms.m.size();
            for (int i = 0; i < n; i++) {
                if (driver.ms.node_is_fixed.at(i)) {
                    driver.ms.x.at(i) += TV(0.f, 0.f, 1.f * dt);
                }
            }
        };
        driver.test="cloth";
    }
    else if (strcmp(argv[1], "1") == 0) // volumetric bunny case
    {
        // TODO
        /*
            1. Create node data from data/points: The first line indicates the number of points and dimension (which is 3).
            2. Fill segments and rest_length from data/cells: The first line indicates the number of tetrahedra and the number of vertices of each tet (which is 6). Each edge in this tetrahedral mesh will be a segment. Be careful not to create duplicate edges.
            3. Choose proper youngs_modulus, damping_coeff, dt
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        */

        std::ifstream pointsStream("data/points");
        std::string line = "";
        getline(pointsStream, line);

        std::istringstream iss(line);
        std::vector<std::string> strings(std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>());

        int n = std::stoi(strings.at(0));
        int ptsDim = std::stoi(strings.at(1));

        std::cout << n << std::endl << dim << std::endl;

        std::unordered_map<int, TV> points;

        int idx = 0;
        while (getline(pointsStream, line)) {
            // Output the text from the file
            iss = std::istringstream(line);
            strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
            TV pt;
            for (int i = 0; i < ptsDim; i++) {
              pt(i, 0) = std::stof(strings.at(i));
            }
            points[idx] = pt;
            // points.insert(std::make_pair<int, TV>(idx, pt));
            // std::cout << points.at(idx)(0, 0) << " " << points.at(idx)(1, 0) << " " << points.at(idx)(2, 0) << std::endl;
            idx++;
        }

        pointsStream.close();

        youngs_modulus = 1.f;
        damping_coeff = 1.f;
        dt = 0.01f;

        driver.helper = [&](T t, T dt) {
            // TODO
        };
        driver.test="bunny";
    }

    else {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    // simulate

    driver.dt = dt;
    driver.ms.segments = segments;
    driver.ms.m = m;
    driver.ms.v = v;
    driver.ms.x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(120);

    return 0;
}
