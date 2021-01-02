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
using TM = Eigen::Matrix<T,dim,dim>;

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

void scatterPoints(T h, TV o, TV resolution, T density, int n, TV cO, TV cShape, std::vector<TV>& X) {
    // density per unit length as opposed to per grid cell
    density /= h;
    TV ptsPerDim = cShape * density;
    ptsPerDim(0) = std::floor(ptsPerDim(0));
    ptsPerDim(1) = std::floor(ptsPerDim(1));
    ptsPerDim(2) = std::floor(ptsPerDim(2));

    TV unit;
    unit << 1, 1, 1;
    TV max = o + h * (resolution - unit);

    std::cout << "max: (" << max(0) << ", " << max(1) << ", " << max(2) << ")" << std::endl;

    T ptSpacing = 1.f / density;
    for (int k = 0; k < ptsPerDim(2); k++) {
        for (int j = 0; j < ptsPerDim(1); j++) {
            for (int i = 0; i < ptsPerDim(0); i++) {
                TV pos = TV::Zero();
                pos(0) = cO(0) + i * ptSpacing;
                pos(1) = cO(1) + j * ptSpacing;
                pos(2) = cO(2) + k * ptSpacing;
                
                if (pos(0) < o(0) || pos(1) < o(1) || pos(2) < o(2)) {
                    std::cout << "computed scattered point position exceeds minimum @ ("
                    << pos(0) << ", " << pos(1) << ", " << pos(2) << ")" << std::endl;
                    return;
                } else if (pos(0) > max(0) || pos(1) > max(1) || pos(2) > max(2)) {
                    std::cout << "computed scattered point position exceeds maximum @ ("
                    << pos(0) << ", " << pos(1) << ", " << pos(2) << ")" << std::endl;
                    return;
                }
                X.push_back(pos);
            }
        }
    }

    return;

}

void scatterPoints(T h, TV o, TV resolution, T density, std::vector<TV>& X) {
    TV unit;
    unit << 1, 1, 1;
    TV max = o + h * (resolution - unit);

    TV cShape;
    cShape << 3, 3, 3;          // in grid coords
    // TV cLWH = cShape * h;       // size in world space coords
    // assuming that the corners of the cube are initialized to be grid-aligned
    // 4 grid cells wide containing 5 points
    TV cO;
    cO << 3 * h, 9 * h, 3 * h;
    // number of nodes per cell in a single dimension
    T linearDensity = std::pow(density, 1.f / 3.f);

    // density per unit length as opposed to per grid cell
    density /= h;
    T ptSpacing = 1.f / density;

    // in total there are 512 nodes
    for (int k = 0; k < 4 * linearDensity; k++) {
        for (int j = 0; j < 4 * linearDensity; j++) {
            for (int i = 0; i < 4 * linearDensity; i++) {
                TV pos = TV::Zero();
                pos(0) = cO(0) + i * ptSpacing;
                pos(1) = cO(1) + j * ptSpacing;
                pos(2) = cO(2) + k * ptSpacing;
                
                if (pos(0) < o(0) || pos(1) < o(1) || pos(2) < o(2)) {
                    std::cout << "computed scattered point position exceeds minimum @ ("
                    << pos(0) << ", " << pos(1) << ", " << pos(2) << ")" << std::endl;
                    return;
                } else if (pos(0) > max(0) || pos(1) > max(1) || pos(2) > max(2)) {
                    std::cout << "computed scattered point position exceeds maximum @ ("
                    << pos(0) << ", " << pos(1) << ", " << pos(2) << ")" << std::endl;
                    return;
                }
                X.push_back(pos);
            }
        }
    }


}

int main(int argc, char* argv[])
{
    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 0;
    T damping_coeff = 0; 
    T dt = T(1.f/24.f);
    int frames = 180;

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
        youngs_modulus = 1000; // TODO: iterate in [0.1, 1, 10, 100, 1000]
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

    else if (strcmp(argv[1], "3") == 0) {
        // mu = shear term
        T shearTerm = 0.f;
        // lambda = diltional term that penalizes volume changes
        T dilationalTerm = 0.f;

        // grid parameters
        T h = 2.f;                  // grid spacing
        TV o = TV::Zero();              // grid origin or minimum of the grid
        TV resolution = TV::Zero();     // nodes per axis
        resolution << 11, 11, 11;

        // particle parameters
        // cube origin (minimum of cube)
        TV cO = TV::Zero();
        cO << 8, 15, 8;
        // cube dimensions
        TV cShape = TV::Zero();
        cShape << 4, 4, 4;      // length width and height of cube
        // desired number of points per grid cell
        T density = std::pow(2, dim);
        int n = (resolution(0) - 1) * (resolution(1) - 1) * (resolution(2) - 1) * density;
        x.clear();
        scatterPoints(h, o, resolution, density, n, cO, cShape, x);
        x.clear();
        scatterPoints(h, o, resolution, density, x);
        std::cout << x.size() << std::endl;
        T mp = 1.f / n;
        v = std::vector<TV>(n, TV(0.f, 0.f, 0.f));
        m = std::vector<T>(n, mp);

        // driver initializations specific to mpm
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << shearTerm;
        driver.test="mpm_"+ss.str();

        driver.ms.grid = Grid<T, dim>(h, o, resolution);

        TM I;
        I << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        driver.ms.F = std::vector<TM>(n, I);

        // designed the box to occupy 4 grid cells in each direction
        T linearDensity = std::pow(density, 1.f / 3.f);
        T nPtsC = std::pow(4 * linearDensity, 3);
        // volume per particle
        driver.ms.vol = std::pow(4 * h, 3) / nPtsC;

        driver.ms.shearTerm = shearTerm;
        driver.ms.dilationalTerm = dilationalTerm;

        frames = 120;
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

    driver.run(frames);

    return 0;
}
