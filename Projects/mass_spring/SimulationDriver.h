#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"
#include <functional>


template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MassSpringSystem<T,dim> ms;
    T dt;
    TV gravity;
    std::string test;

    std::function<void(T, T)> helper = [&](T, T) {};

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;
    }

    void run(const int max_frame)
    {
        T accumulate_t = 0;
        mkdir("output/", 0777);
        std::string output_folder = "output/" + test;
        mkdir(output_folder.c_str(), 0777);
        std::string filename = output_folder + "/" + std::to_string(0) + ".poly";
        ms.dumpPoly(filename);
        for(int frame=1; frame<=max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;
            // this is necessary so that you know how many time steps have occurred
            // between each frame. this calculation assumes that the frame rate is
            // 24 frames per second.
            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                // std::cout << "Step " << step << std::endl;
                helper(accumulate_t, dt);
                advanceOneStepExplicitIntegration();
                accumulate_t += dt;
            }
            mkdir("output/", 0777);
            std::string output_folder = "output/" + test;
            mkdir(output_folder.c_str(), 0777);
            std::string filename = output_folder + "/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            // std::cout << std::endl;
        }
    }

    void advanceOneStepExplicitIntegration()
    {
        int n = ms.m.size();
        std::vector<TV> f_spring(n);
        ms.evaluateSpringForces(f_spring);
	      std::vector<TV> f_damping(n);
	      ms.evaluateDampingForces(f_damping);

        // TODO: update position and velocity according to Newton's law.
        for (int i = 0; i < n; i++) {
            if (!ms.node_is_fixed.at(i)) {
                ms.v.at(i) += (dt * (gravity + f_spring.at(i) + f_damping.at(i)) / ms.m.at(i));
                // if the force matrix contains any nan values
                // the inequality expression will return true since NaN != NaN
                // is always a true statement
                // ms.v.at(i) += (dt * (gravity + f_spring.at(i)) / ms.m.at(i));
                ms.x.at(i) += (dt * ms.v.at(i));
                if (x.at(i) != x.at(i))
                  std::cout << "ERROR: point " << i << std::endl;
                if (f_spring.at(i) != f_spring.at(i))
                  std::cout << "ERROR: spring force" << std::endl;
            }
        }
    }
};
