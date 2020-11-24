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

void writeObj(std::string filename, std::vector<TV>& X, std::vector<Eigen::Vector4i>& F) {
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

void readBunny(std::vector<TV>& X, std::vector<Eigen::Matrix<int,2,1>>& S, std::vector<T>& RL) {
  std::ifstream pointsStream("data/points");
  std::string line = "";
  getline(pointsStream, line);

  std::istringstream iss(line);
  std::vector<std::string> strings(std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>());

  // contains the number of points in the file
  // int n = std::stoi(strings.at(0));
  // contains the dimension of the points which is 3
  int ptsDim = std::stoi(strings.at(1));

  // according to the documentation for getLine, it is unnecessary
  // to empty the string contents since it is done for you, but this 
  // is just for sanity's sake.
  line = "";
  while (getline(pointsStream, line)) {
    iss = std::istringstream(line);
    strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    TV pt;
    for (int i = 0; i < ptsDim; i++) {
      pt(i, 0) = std::stof(strings.at(i));
    }
    X.push_back(pt);
    line = "";
    // points.insert(std::make_pair<int, TV>(idx, pt));
    // std::cout << pt(0, 0) << " " << pt(1, 0) << " " << pt(2, 0) << std::endl;
  }

  pointsStream.close();

  // read in data from cells to obtain the topology
  // of the model
  std::ifstream cellsStream("data/cells");
  getline(cellsStream, line);
  iss = std::istringstream(line);
  strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
  // int numPolys = std::stoi(strings.at(0));
  int sidesPerPoly = std::stoi(strings.at(1));

  std::unordered_map<long long int, Eigen::Matrix<int,2,1>> uniqueSegments;

  line = "";
  while (getline(cellsStream, line)) {
    iss = std::istringstream(line);
    strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    for (int i = 0; i < sidesPerPoly; i++) {
      for (int j = i + 1; j < sidesPerPoly; j++) {
        long long int segIdx1 = 0;
        long long int segIdx2 = 0;
        int pt1 = std::stoi(strings.at(i));
        int pt2 = std::stoi(strings.at(j));
        segIdx1 = pt1;
        segIdx1 |= ((long long int)pt2 << 32);
        segIdx2 = pt2;
        segIdx2 |= ((long long int)pt1 << 32);
        // check both orderings of the pt indices in the segIdx
        if (!uniqueSegments.count(segIdx1) && !uniqueSegments.count(segIdx2)) {
          uniqueSegments[segIdx1] = Eigen::Matrix<int,2,1>(pt1, pt2);
          S.push_back(Eigen::Matrix<int,2,1>(pt1, pt2));
          RL.push_back(sqrt((X.at(pt1) - X.at(pt2)).transpose() * (X.at(pt1) - X.at(pt2))));
        }
      }
    }
    line = "";
    // for (int i = 0; i < sidesPerPoly && (((long unsigned int)i) < strings.size()); i++) {)
    //   long long int segIdx1 = 0;
    //   long long int segIdx2 = 0;
    //   int pt1 = std::stoi(strings.at(i));
    //   int pt2 = std::stoi(strings.at((i + 1) % sidesPerPoly));
    //   segIdx1 = pt1;
    //   segIdx1 |= ((long long int)pt2 << 32);
    //   segIdx2 = pt2;
    //   segIdx2 |= ((long long int)pt1 << 32);
    //   // check both orderings of the pt indices in the segIdx
    //   if (!uniqueSegments.count(segIdx1) && !uniqueSegments.count(segIdx2)) {
    //     uniqueSegments[segIdx1] = Eigen::Matrix<int,2,1>(pt1, pt2);
    //     S.push_back(Eigen::Matrix<int,2,1>(pt1, pt2));
    //     RL.push_back(sqrt((X.at(pt1) - X.at(pt2)).transpose() * (X.at(pt1) - X.at(pt2))));
    //   }
    // }
  }

  cellsStream.close();
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

    // segment/spring data
    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<Eigen::Matrix<int,2,1>> springs;
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

    		int xN = 20;
    		int yN = 20;
    		int N = xN * yN;
    		T xW = 5.f;
    		T yW = 5.f;
    		T xWHalf = xW / 2.f;
    		T yWHalf = yW / 2.f;
    		T mN = 5.f / N;

        youngs_modulus = 100.f;
        damping_coeff = 10.f;
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
            // std::cout << x.back()(0,0) << " " << x.back()(1,0) << " " << x.back()(2,0) << std::endl;
            if ((i == 0 && j == yN - 1) || (i == xN - 1 && j == yN - 1)) {
              node_is_fixed.push_back(true);
            } else {
              node_is_fixed.push_back(false);
            }

            int idx = i + j * xN;

            // structural springs
            if (j < yN - 1) {
              segments.push_back(Eigen::Vector2i(idx, idx + xN));
              springs.push_back(Eigen::Vector2i(idx, idx + xN));
              rest_length.push_back(rLSty);
            }
            if (i < xN - 1) {
              segments.push_back(Eigen::Vector2i(idx, idx + 1));
              springs.push_back(Eigen::Vector2i(idx, idx + 1));
              rest_length.push_back(rLStx);
            }

            // shear springs
            if (i < xN - 1 && j < yN - 1) {
              springs.push_back(Eigen::Vector2i(idx, idx + 1 + xN));
              springs.push_back(Eigen::Vector2i(idx + 1, idx + xN));
              rest_length.push_back(rLSh);
              rest_length.push_back(rLSh);
            }

            // bend (flexion) springs
            if (j < yN - 2) {
              springs.push_back(Eigen::Vector2i(idx, idx + 2 * xN));
              rest_length.push_back(rLBdy);
            }
            if (i < xN - 2) {
              springs.push_back(Eigen::Vector2i(idx, idx + 2));
              rest_length.push_back(rLBdx);
            }

            if (i > 0 && j > 0) {
              // the current point should be the bottom right corner
              // of its face. this vertex will be responsible for
              // constructing the face which will eventually be written to an
              // obj file
              faces.push_back(Eigen::Vector4i(idx, idx - 1, idx - 1 - xN, idx - xN));
            }
    			}
    		}

        writeObj("cloth.obj", x, faces);


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

    } else if (strcmp(argv[1], "1") == 0) {
        // volumetric bunny case
        // TODO
        /*
            1. Create node data from data/points: The first line indicates the number of points and dimension (which is 3).
            2. Fill segments and rest_length from data/cells: The first line indicates the number of tetrahedra and the number of vertices of each tet (which is 4). Each edge in this tetrahedral mesh will be a segment. Be careful not to create duplicate edges.
            3. Choose proper youngs_modulus, damping_coeff, dt
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        */

        // initialize x, segments, and rest_length using the file information
        readBunny(x, segments, rest_length);

        // parameters
        int n = x.size();
        T mN = 0.5f / n;

        // youngs_modulus = 3.f;
        // damping_coeff = 0.1f;
        youngs_modulus = 70.f;
        damping_coeff = 0.5f;
        // unlike in the cloth simulation,
        // i had to increase the time resolution by a 
        // power of 10 to avoid instability
        // with dt = 0.0001f
        dt = 0.00001f;

        v = std::vector<TV>(n, TV(0.f, 0.f, 0.f));
        m = std::vector<T>(n, mN);
        springs = std::vector<Eigen::Matrix<int,2,1>>(segments);

        node_is_fixed = std::vector<bool>(n, false);
        // ear tips
        node_is_fixed.at(2140) = true;
        node_is_fixed.at(2346) = true;
        node_is_fixed.at(1036) = true;

        driver.helper = [&](T t, T dt) {
            // TODO
            if (t >= 1.f) {
              driver.ms.node_is_fixed.at(1036) = false;
            } else {
              driver.ms.x.at(1036) += TV(0.3f * dt, 0.f, 0.f);
            }


            // int n = driver.ms.m.size();
            // for (int i = 0; i < n; i++) {
            //     if (i == 1036) {
            //       if (t >= 3.f) {
            //         driver.ms.node_is_fixed.at(i) = false;
            //       } else {
            //         driver.ms.x.at(i) += TV(0.3f * dt, 0.f, 0.f);
            //       }
            //     }
            //     if (driver.ms.node_is_fixed.at(i)) {
            //         driver.ms.x.at(i) += TV(0.f, 0.f, 0.3f * dt);
            //     }
            // }
        };
        driver.test="bunny";
    } else {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    // simulate
    driver.dt = dt;
    driver.ms.segments = segments;
    driver.ms.springs = springs;
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
