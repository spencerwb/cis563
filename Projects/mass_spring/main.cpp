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
#include <iomanip>

#include "SimulationDriver.h"

#include <iterator>
#include <unordered_map>

// type definitions
using T = float;
constexpr int dim = 3;
using TV = Eigen::Matrix<T,dim,1>;

void readBunny(std::vector<TV>& X, std::vector<Eigen::Matrix<int,2,1>>& S, std::vector<T> RL) {
  std::ifstream pointsStream("data/points");
  std::string line = "";
  getline(pointsStream, line);

  std::istringstream iss(line);
  std::vector<std::string> strings(std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>());

  // contains the dimension of the points which is 3
  // int n = std::stoi(strings.at(0));
  int ptsDim = std::stoi(strings.at(1));

  while (getline(pointsStream, line)) {
      iss = std::istringstream(line);
      strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
      TV pt;
      for (int i = 0; i < ptsDim; i++) {
        pt(i, 0) = std::stof(strings.at(i));
      }
      X.push_back(pt);
      // points.insert(std::make_pair<int, TV>(idx, pt));
      // std::cout << points.at(idx)(0, 0) << " " << points.at(idx)(1, 0) << " " << points.at(idx)(2, 0) << std::endl;
  }

  pointsStream.close();

  // read in data from cells to obtain the topology
  // of the model
  std::ifstream cellsStream("data/cells");
  getline(pointsStream, line);
  iss = std::istringstream(line);
  strings = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
  // int numFaces = std::stoi(strings.at(0));
  int sidesPerPoly = std::stoi(strings.at(1));

  std::unordered_map<long long int, Eigen::Matrix<int,2,1>> uniqueSegments;

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
    T youngs_modulus = 0;
    T damping_coeff = 0; 
    T dt = T(1.f/24.f);

    // node data
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    // segment data
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> rest_length;

    if (argc < 2) 
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "1") == 0) { // bunny case
        youngs_modulus = 0.1; // TODO: iterate in [0.1, 1, 10, 100, 1000]
        damping_coeff = 2;
        // TODO: 
        /* 
            1. Copy the loading codes from your hw1. Fix two ears (2140, 2346) only, and you don't need helper function here.
            2. Set the initial velocities of non_fixed_nodes to be (10, 0, 0)
            
            The output folder will automatically renamed by bunny_[youngs_modulus], don't worry about overwriting.
        */

        // initialize x, segments, and rest_length using the file information
        readBunny(x, segments, rest_length);

        // parameters
        int n = x.size();
        T mN = 1.f / n;

        v = std::vector<TV>(n, TV(10.f, 0.f, 0.f));
        v.at(2140) = TV(0.f, 0.f, 0.f);
        v.at(2346) = TV(0.f, 0.f, 0.f);
        m = std::vector<T>(n, mN);

        node_is_fixed = std::vector<bool>(n, false);
        node_is_fixed.at(2140) = true;
        node_is_fixed.at(2346) = true;

        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << youngs_modulus;
        driver.test="bunny_"+ss.str();
    }

    else if (strcmp(argv[1], "2") == 0) { //brush case
        driver.gravity.setZero();
        driver.collision_stiffness = 0.1;
        youngs_modulus = 10000;
        damping_coeff = 100; // 0
        int N = 32; // z direction
        int M = 4; // y direction
        int L = 32; // x direction
        int N_points = N*M*L;
        T dx = (T)0.1/(N-1);
        m.resize(N_points);
        x.resize(N_points);
        v.resize(N_points);
        node_is_fixed = std::vector<bool>(N_points, false);
        for(int i=0; i<N; i++){ // z
            for(int j=0; j<M; j++) { // y
                for (int k=0; k<L; k++) { // x
                    int id = i * M * L + j * L + k;
                    m[id] = (T)0.001/N_points;
                    x[id](2) = i*dx;
                    x[id](1) = j*dx;
                    x[id](0) = k*dx;
                    v[id] = TV::Zero();
                    if (k <= 2) node_is_fixed[id] = true;
                    // struct spring
                    if (k > 0) { 
                        segments.push_back(Eigen::Matrix<int,2,1>(id, id-1));
                        rest_length.push_back((x[id]-x[id-1]).norm());
                    }
                    // bending spring
                    if (k > 1) {
                        segments.push_back(Eigen::Matrix<int,2,1>(id, id-2));
                        rest_length.push_back((x[id]-x[id-2]).norm());
                    }
                }
            }
        }

        driver.sphere_radius = 0.04;
        driver.sphere_center = TV(0.07, -0.045, 0.05);

        driver.helper = [&](T t, T dt) {
            if(t < 4) {
                // driver.sphere_center = TV(0.12, t/4. * 0.15 + (1-t/4.) * (-0.06), 0.05);
                for(size_t i = 0; i < driver.ms.x.size(); ++i) {
                    if (driver.ms.node_is_fixed[i]) {
                        driver.ms.target_x[i](1) -= 0.15 / 4 * dt;
                    }
                }
            }
        };

        driver.test="brush";
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
    driver.ms.target_x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(180);

    return 0;
}
