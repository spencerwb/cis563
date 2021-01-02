#include <vector>
#include <cmath>

template<class T, int dim>
class Grid {
    public:
        using TV = Eigen::Matrix<T,dim,1>;
        using TM = Eigen::Matrix<T,dim,dim>;

        T h;            // grid spacing
        TV o;           // grid origin or minimum of the grid
        TV resolution;  // nodes per axis
        TV size;        // size of the grid in world space
        int n;          // number of nodes in the grid
        TV max;

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
            max = o + size;

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
            max = o + size;

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
                for (int j = Bp(1); j <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); k <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        TV bpos = gridIdx3DToPos(bp);
                        if (bpos(0) > max(0) || bpos(1) > max(1) || bpos(2) > max(2)) continue;
                        float w = kernel(xp, bpos);
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
                for (int j = Bp(1); j <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); k <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        TV bpos = gridIdx3DToPos(bp);
                        if (bpos(0) > max(0) || bpos(1) > max(1) || bpos(2) > max(2)) continue;
                        float w = kernel(xp, bpos);
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

        void applyGravity(float dt) {
            TV g = TV::Zero();
            g(1) = -9.8f;
            // for (int i = 0; i < n; i++) {
            //     if ()
            //     v.at(i) += g * dt;
            // }
            for (int i = 0; i < resolution(0); i++) {
                for (int j = 0; j < resolution(1); j++) {
                    for (int k = 0; k < resolution(2); k++) {
                        TV bp = TV();
                        bp << i, j, k;
                        
                        TV bpos = gridIdx3DToPos(bp);
                        // if this grid node is at the bottom of the box, then it
                        // is assumed to be the ground
                        if (bpos(1) == 0) {
                            v.at(gridIdx3DTo1D(bp)) = TV::Zero();
                        } else {
                            v.at(gridIdx3DTo1D(bp)) += g * dt;
                        }
                    }
                }
            }
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
                for (int j = Bp(1); j <= Bp(1) + 2; j++) {
                    for (int k = Bp(2); k <= Bp(2) + 2; k++) {
                        TV bp = TV();
                        bp(0) = i;
                        bp(1) = j;
                        bp(2) = k;
                        
                        TV bpos = gridIdx3DToPos(bp);
                        if (bpos(0) > max(0) || bpos(1) > max(1) || bpos(2) > max(2)) continue;
                        float w = kernel(xp, bpos);
                        vp += v.at(gridIdx3DTo1D(bp)) * w;
                    }
                }
            }

            return vp;
        }

        T BSplineInterpDerivative(T x) {
            T ax = std::abs(x);
            if (0 <= ax && ax < 0.5f) {
                return -2 * x;
            } else if (0.5 <= ax && ax < 1.5f) {
                return x * (1.5 - ax) / ax;
            } else if (1.5f <= x) {
                return 0.f;
            }
            return INFINITY;
        }

        TV kernelGradient(TV xp, TV xi) {
            TV gradW;
            TV xDiffh = (xp - xi) / h;
            gradW(0) = BSplineInterpDerivative(xDiffh(0)) * BSplineInterp(xDiffh(1)) * BSplineInterp(xDiffh(2)) / h;
            gradW(1) = BSplineInterp(xDiffh(0)) * BSplineInterpDerivative(xDiffh(1)) * BSplineInterp(xDiffh(2)) / h;
            gradW(1) = BSplineInterp(xDiffh(0)) * BSplineInterp(xDiffh(1)) * BSplineInterpDerivative(xDiffh(2)) / h;
            return gradW;
        }

        // the particle's position, volume, and current deformation gradient
        void computeElasticForce(TV xp, T vol, TM F, T shearTerm, T dilationalTerm) {
            for (int i = 0; i < resolution(0); i++) {
                for (int j = 0; j < resolution(1); j++) {
                    for (int k = 0; k < resolution(2); k++) {
                        TV bp = TV();
                        bp << i, j, k;
                        
                        TV bpos = gridIdx3DToPos(bp);
                        // STOP if the current node is out of bounds
                        if (bpos(0) > max(0) || bpos(1) > max(1) || bpos(2) > max(2)) continue;
                        
                        int idx = gridIdx3DTo1D(bp);

                        Eigen::JacobiSVD<TM> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
                        TM U = svd.matrixU();
                        TM VT = svd.matrixV().transpose();
                        TM R = U * VT;
                        // auto sv = svd.singularValues();
                        // T s0 = sv(0);
                        // T s1 = sv(1);
                        // T s2 = sv(2);

                        T J = F.determinant();
                        // std::cout << s0 << ", " << s1 << ", " << s2 << std::endl;
                        TM P = 2 * shearTerm * (F - R) + dilationalTerm * (J - 1) * J * F.inverse().transpose();
                        TV gradW = kernelGradient(xp, bpos);
                        f.at(idx) -= vol * P * F.transpose() * gradW;
                    }
                }
            }
        }


        void applyElasticForces(float dt) {
            for (int i = 0; i < this->n; i++) {
                if (i >= 1331) {
                    std::cout << "here at idx " << i << std::endl;
                }
                v.at(i) += (dt * f.at(i) / m.at(i));
                // TV vi = (dt * f.at(i) / m.at(i));
                // if (false) {
                //     vi += TV::Zero();
                // }

                // const T vx = vi(0);
                // // std::cout << "here at idx " << i << std::endl;
                // v.at(i) << vx, 0, 0;
            }

            return;
        }


};