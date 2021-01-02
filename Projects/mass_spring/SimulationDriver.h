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
    T frame_dt;
    TV gravity;

    TV sphere_center;
    T sphere_radius;
    std::string test;

    T collision_stiffness = 0;
    std::vector<bool> has_collision;

    std::function<void(T, T)> helper = [&](T, T) {};

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;

        sphere_center = TV::Zero();
        sphere_radius = 0.0;

        frame_dt = (T)1/24;
    }

    void run(const int max_frame)
    {
        // TODO: if you are not sure whether your computation is correct. You can run the following two functions to check f and df/dx
        //          when you run the simulation, comment out these functions.
        // ms.checkGradient();
        // ms.checkHessian();
        T accumulate_t = 0;
        mkdir("output/", 0777);
        std::string output_folder = "output/" + test;
        mkdir(output_folder.c_str(), 0777);
        ms.dumpPoly(output_folder + "/" + std::to_string(0) + ".poly");
        for(int frame=1; frame<=max_frame; frame++) {
            std::cout << "=============================== Frame " << frame << " ===============================" << std::endl;
            T remain_dt = frame_dt;
            T current_dt = std::min(dt, frame_dt);
            while (remain_dt > 0) {
                if (remain_dt > dt * 2)
                    current_dt = dt;
                else if (remain_dt > dt)
                    current_dt = remain_dt / 2;
                else
                    current_dt = remain_dt;
                advanceOneStepMPM(current_dt);
                accumulate_t += current_dt;
                remain_dt -= current_dt;
                std::cout << "Frame " << frame << " finished " << int(100 - 100 * remain_dt/frame_dt) << "%" << std::endl;
            }
            mkdir("output/", 0777);
            std::string output_folder = "output/" + test;
            mkdir(output_folder.c_str(), 0777);
            std::string filename = output_folder + "/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceOneStepMPM(T dt) {

        int n = ms.x.size();

        // std::cout << "CLEAR GRID" << std::endl;
        ms.grid.clear();

        // particles to grid transfer
        for (int i = 0; i < n; i++) {
            // std::cout << "UPDATE GRID MASS" << std::endl;
            ms.grid.updateGridMass(ms.x.at(i), ms.m.at(i));
            // std::cout << "UPDATE GRID MOMENTUM" << std::endl;
            ms.grid.updateGridMomentum(ms.x.at(i), ms.v.at(i), ms.m.at(i));
        }

        // std::cout << "MOMENTUM TO VELOCITY" << std::endl;
        ms.grid.momentumToVelocity();

        // std::cout << "APPLYING GRAVITY" << std::endl;
        ms.grid.applyGravity(dt);

        // std::cout << "COMPUTING ELASTIC FORCES" << std::endl;
        // for (int i = 0; i < n; i++) {
        //     ms.grid.computeElasticForce(ms.x.at(i), ms.vol, ms.F.at(i), ms.shearTerm, ms.dilationalTerm);
        // }
        // std::cout << "APPLYING FORCES" << std::endl;
        // ms.grid.applyElasticForces(dt);

        std::cout << "GRID TO PARTICLE TRANSFER" << std::endl;
        // grid to particles transfer
        for (int i = 0; i < n; i++) {
            // std::cout << "to v" << i;
            ms.v.at(i) = ms.grid.computeParticleVelocity(ms.x.at(i));
            // std::cout << " done ";
            // if (ms.v.at(i) == TV::Zero()) {
            //     std::cout << "error at particle " << i << std::endl;
            //     return;
            // }
            // std::cout << "to x" << i;
            ms.x.at(i) += ms.v.at(i) * dt;
            // std::cout << " done " << std::endl;
            // TV five;
            // five << 0, -5, 0;
            // ms.x.at(i) += five * dt;
        }

        return;

    }
};
