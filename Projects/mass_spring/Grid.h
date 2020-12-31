#include <vector>
#include <cmath>

template<class T, int dim>
class Grid {
    public:
        using TV = Eigen::Matrix<T,dim,1>;

        T h;            // grid spacing
        TV o;           // grid origin or minimum of the grid
        TV resolution;  // nodes per axis
        TV size;        // size of the grid
        int n;          // number of nodes in the grid

        std::vector<T> m;
        std::vector<TV> v;
        std::vector<TV> f;

        Grid() 
        : h(1.f), o(TV::Zero())
        {
            resolution(0) = 10;
            resolution(1) = 10;
            resolution(2) = 10;
            n = int(resolution(0) * resolution(1) * resolution(2));

            TV unit;
            unit << 1, 1, 1;
            size = h * (resolution - unit);

            m = std::vector<T>(n, 0.f);
            v = std::vector<TV>(n, TV::Zero());
            f = std::vector<TV>(n, TV::Zero());
        }

        Grid(T h, TV o, TV resolution) 
        : h(h), o(o), resolution(resolution)
        {
            n = int(resolution(0) * resolution(1) * resolution(2));

            TV unit;
            unit << 1, 1, 1;
            size = h * (resolution - unit);

            m = std::vector<T>(n, 0.f);
            v = std::vector<TV>(n, TV::Zero());
            f = std::vector<TV>(n, TV::Zero());
        }

        void clear() {
            m = std::vector<T>(n, 0.f);
            v = std::vector<TV>(n, TV::Zero());
            f = std::vector<TV>(n, TV::Zero());
        }

        T BSplineInterp(T x) {
            x = std::abs(x);
            if (0 <= x && x < 0.5f) {
                return 0.75f - x * x;
            } else if (0.5 <= x && x < 1.5f) {
                return 0.5f * (1.5f - x) * (1.5f - x);
            } else if (1.5f <= x) {
                return 0.f;
            }
            return INFINITY;
        }

        T kernel(TV xp, TV xi) {
            return BSplineInterp((xp(0) - xi(0)) / h) * BSplineInterp((xp(1) - xi(1)) / h) * BSplineInterp((xp(2) - xi(2)) / h);
        }

        int gridIdx3DTo1D(TV idx) {
            return idx(0) + idx(1) * resolution(0) + idx(2) * resolution(0) * resolution(1);
        }

        TV gridIdx3DToPos(TV idx) {
            return idx * h - o;
        }

        // input is the particle position and mass
        void updateGridMass(TV xp, T mp) {
            // base index
            TV half;
            half << 0.5, 0.5, 0.5;
            TV Bp = (xp - o - half * h) / h;
            Bp(0) = int(Bp(0));
            Bp(1) = int(Bp(1));
            Bp(2) = int(Bp(2));
            for (int i = Bp(0); i <= Bp(0) + 2; i++) {
                for (int j = Bp(1); i <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); i <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        float w = kernel(xp, gridIdx3DToPos(bp));
                        m.at(gridIdx3DTo1D(bp)) += mp * w;
                    }
                }
            }
        }

        // input is the particle position, velocity, and mass
        void updateGridMomentum(TV xp, TV vp, T mp) {
            // base index
            TV half;
            half << 0.5, 0.5, 0.5;
            TV Bp = (xp - o - half * h) / h;
            Bp(0) = int(Bp(0));
            Bp(1) = int(Bp(1));
            Bp(2) = int(Bp(2));
            for (int i = Bp(0); i <= Bp(0) + 2; i++) {
                for (int j = Bp(1); i <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); i <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        float w = kernel(xp, gridIdx3DToPos(bp));
                        v.at(gridIdx3DTo1D(bp)) += mp * w * vp;
                    }
                }
            }
        }

        void momentumToVelocity() {
            for (int i = 0; i < n; i++) {
                if (m.at(i) != 0) {
                    v.at(i) = v.at(i) / m.at(i);
                } else {
                    v.at(i) = TV::Zero();
                }
            }
        }

        void applyForces(float dt) {
            TV g = TV::Zero();
            g(1) = -9.8f;
            for (int i = 0; i < n; i++)
                v.at(i) += g * dt;
        }

        TV computeParticleVelocity(TV xp) {
            TV vp = TV::Zero();
            // base index
            TV half;
            half << 0.5, 0.5, 0.5;
            TV Bp = (xp - o - half * h) / h;
            Bp(0) = int(Bp(0));
            Bp(1) = int(Bp(1));
            Bp(2) = int(Bp(2));
            for (int i = Bp(0); i <= Bp(0) + 2; i++) {
                for (int j = Bp(1); i <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); i <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        float w = kernel(xp, gridIdx3DToPos(bp));
                        vp += v.at(gridIdx3DTo1D(bp)) * w;
                    }
                }
            }
            return vp;
        }
};