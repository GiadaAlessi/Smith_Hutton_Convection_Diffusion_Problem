#include "Smith_Hutton_Upgraded.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace std::chrono;

// Global Variables
double const L = 2.0;           // Lenght
double const H = 1.0;           // Height
double const W = 1.0;           // Width
const int Nx = 100, Ny = 50;    // Number of control volumes (100x50)

// Vectors Definition
double phi[Ny+2][Nx+2], phi_1[Ny+2][Nx+2], phi_g[Ny+2][Nx+2];   // Transported property
double aP[Ny+2][Nx+2], aE[Ny+2][Nx+2], aW[Ny+2][Nx+2],          // Coefficients
        aN[Ny+2][Nx+2], aS[Ny+2][Nx+2], bP[Ny+2][Nx+2];
double xcv[Ny+2][Nx+2], ycv[Ny+2][Nx+2], xP[Ny+2][Nx+2],        // Mesh
        yP[Ny+2][Nx+2], vP[Ny+2][Nx+2];
double u[Ny+2][Nx+2], v[Ny+2][Nx+2];                            // Velocity components

// Function to generate the mesh and define position of nodes, faces and volume of control volumes
void MeshDefinition(double H, double L, double W, double Dx, double Dy) {
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            if (i == 0 && j == 0)
                xcv[i][j] = -L / 2,
                ycv[i][j] = H,
                xP[i][j] = -L / 2,
                yP[i][j] = H,
                vP[i][j] = 0;
            else if (i == 0 && j == Nx+1)
                xcv[i][j] = L / 2,
                ycv[i][j] = H,
                xP[i][j] = L / 2,
                yP[i][j] = H,
                vP[i][j] = 0;
            else if (i == Ny+1 && j == 0)
                xcv[i][j] = -L / 2,
                ycv[i][j] = 0,
                xP[i][j] = -L / 2,
                yP[i][j] = 0,
                vP[i][j] = 0;
            else if (i == Ny+1 && j == Nx+1)
                xcv[i][j] = L / 2,
                ycv[i][j] = 0,
                xP[i][j] = L / 2,
                yP[i][j] = 0,
                vP[i][j] = 0;
            else if (i == 0)
                xcv[i][j] = -L / 2 + Dx * j,
                ycv[i][j] = H,
                xP[i][j] = (xcv[i][j] + xcv[i][j-1]) / 2,
                yP[i][j] = H,
                vP[i][j] = 0;
            else if (j == 0)
                xcv[i][j] = -L / 2,
                ycv[i][j] = H - Dy * i,
                xP[i][j] = -L / 2,
                yP[i][j] = (ycv[i][j] + ycv[i-1][j]) / 2,
                vP[i][j] = 0;
            else if (j == Nx+1)
                xcv[i][j] = L / 2,
                ycv[i][j] = H - Dy * i,
                xP[i][j] = L / 2,
                yP[i][j] = (ycv[i][j] + ycv[i-1][j]) / 2,
                vP[i][j] = 0;
            else if (i == Ny+1)
                xcv[i][j] = -L / 2 + Dx * j,
                ycv[i][j] = 0,
                xP[i][j] = (xcv[i][j] + xcv[i][j-1]) / 2,
                yP[i][j] = 0,
                vP[i][j] = 0;
        }
    }
    for (int i = 1; i < Ny+1; i++) {
        for (int j = 1; j < Nx+1; j++) {
            xcv[i][j] = -L / 2 + Dx * j;
            ycv[i][j] = H - Dy * i;
            xP[i][j] = (xcv[i][j] + xcv[i][j-1]) / 2;
            yP[i][j] = (ycv[i][j] + ycv[i-1][j]) / 2;
            vP[i][j] = W * fabs(xcv[i][j] - xcv[i][j-1]) * fabs(ycv[i][j] - ycv[i-1][j]);
        }
    }
}

// Function to assign velocity value at each face
void VelocityField() {
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            u[i][j] = 2 * yP[i][j] * (1 - pow(xcv[i][j],2));    // ue of node P, if I want uw I'll write u[i][j-1]
            v[i][j] = -2 * xP[i][j] * (1 - pow(ycv[i][j],2));   // vs of node P, if I want vn I'll write u[i-1][j]
        }
    }
}

// Base Class to implement multiple methods
class InternalNodesCalculator {
public:
    virtual void computeCoefficients(double rho, double dt, double gamma) = 0;  // Pure virtual function
    virtual ~InternalNodesCalculator() = default;                               // Virtual destructor
};

// Central Differencing Scheme (CDS) method
class CDSMethod : public InternalNodesCalculator {
    // Now computeCoefficients() is a method of the class, which can be called with CDSMethod object
public:
    void computeCoefficients(double rho, double dt, double gamma) override {
        for (int i = 1; i < Ny + 1; i++) {
            for (int j = 1; j < Nx + 1; j++) {

                double DSe = (ycv[i - 1][j] - ycv[i][j]) * L;
                double DSw = (ycv[i - 1][j] - ycv[i][j]) * L;
                double DSn = (xcv[i][j] - xcv[i][j - 1]) * L;
                double DSs = (xcv[i][j] - xcv[i][j - 1]) * L;

                double dPE = xP[i][j + 1] - xP[i][j];
                double dPW = xP[i][j] - xP[i][j - 1];
                double dPN = yP[i - 1][j] - yP[i][j];
                double dPS = yP[i][j] - yP[i + 1][j];

                aE[i][j] = DSe / vP[i][j] * (-rho * u[i][j] / 2 + gamma / dPE);
                aW[i][j] = DSw / vP[i][j] * (rho * u[i][j - 1] / 2 + gamma / dPW);
                aN[i][j] = DSn / vP[i][j] * (-rho * v[i - 1][j] / 2 + gamma / dPN);
                aS[i][j] = DSs / vP[i][j] * (rho * v[i][j] / 2 + gamma / dPS);
                aP[i][j] = rho / dt;
                bP[i][j] = rho / dt + DSe / vP[i][j] * (-rho * u[i][j] / 2 - gamma / dPE) +
                                                 DSw / vP[i][j] * (rho * u[i][j - 1] / 2 - gamma / dPW) +
                                                 DSn / vP[i][j] * (-rho * v[i - 1][j] / 2 - gamma / dPN) +
                                                 DSs / vP[i][j] * (rho * v[i][j] / 2 - gamma / dPS);
            }
        }
    }
};

// Upwind Differencing Scheme (UDS) method
class UDSMethod : public InternalNodesCalculator {
    // Now computeCoefficients() is another method of the class, which can be called with UDSMethod object
public:
    void computeCoefficients(double rho, double dt, double gamma) override {
        for (int i = 1; i < Ny + 1; i++) {
            for (int j = 1; j < Nx + 1; j++) {

                double DSe = (ycv[i - 1][j] - ycv[i][j]) * L;
                double DSw = (ycv[i - 1][j] - ycv[i][j]) * L;
                double DSn = (xcv[i][j] - xcv[i][j - 1]) * L;
                double DSs = (xcv[i][j] - xcv[i][j - 1]) * L;

                double dPE = xP[i][j + 1] - xP[i][j];
                double dPW = xP[i][j] - xP[i][j - 1];
                double dPN = yP[i - 1][j] - yP[i][j];
                double dPS = yP[i][j] - yP[i + 1][j];

                aE[i][j] = DSe / vP[i][j] * (max(-rho * u[i][j], 0.0) + gamma / dPE);
                aW[i][j] = DSw / vP[i][j] * (max(rho * u[i][j - 1], 0.0) + gamma / dPW);
                aN[i][j] = DSn / vP[i][j] * (max(-rho * v[i - 1][j], 0.0) + gamma / dPN);
                aS[i][j] = DSs / vP[i][j] * (max(rho * v[i][j], 0.0) + gamma / dPS);
                aP[i][j] = rho / dt;
                bP[i][j] = rho / dt + DSe / vP[i][j] * (min(-rho * u[i][j], 0.0) - gamma / dPE) +
                                                 DSw / vP[i][j] * (min(rho * u[i][j - 1], 0.0) - gamma / dPW) +
                                                 DSn / vP[i][j] * (min(-rho * v[i - 1][j], 0.0) - gamma / dPN) +
                                                 DSs / vP[i][j] * (min(rho * v[i][j], 0.0) - gamma / dPS);
            }
        }
    }
};

// Function to choose the method
InternalNodesCalculator* createCalculator(string method) {
    if (method == "CDS") {
        return new CDSMethod();
    } else if (method == "UDS") {
        return new UDSMethod();
    } else {
        cerr << "Not a valid method!" << std::endl;
        return nullptr;
    }
}

// Function to apply boundary conditions
void BoundaryConditions(double alpha) {
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            if (i == Ny+1 && xP[i][j] < 0)                                  // Inlet Condition
                aE[i][j] = 0,
                aW[i][j] = 0,
                aN[i][j] = 0,
                aS[i][j] = 0,
                aP[i][j] = 1,
                bP[i][j] = 1 + tanh((2 * xP[i][j] + 1) * alpha);
            else if (i == Ny+1 && xP[i][j] > 0)                             // Outlet Condition
                aE[i][j] = 0,
                aW[i][j] = 0,
                aN[i][j] = 1,
                aS[i][j] = 0,
                aP[i][j] = 1,
                bP[i][j] = 0;
            else if (i == 0 || j == 0 || j == Nx+1)                         // Rest of the boundaries
                aE[i][j] = 0,
                aW[i][j] = 0,
                aN[i][j] = 0,
                aS[i][j] = 0,
                aP[i][j] = 1,
                bP[i][j] = 1 - tanh(alpha);
        }
    }
}

void Solver(double maxRes, double &t_count, double dt, double rho, double gamma) {
    // Time Loop
    double res = maxRes + 1; // Condition to enter the loop
    while (res > maxRes) {
        double maxDiff = 0.0;

        for (int i = 0; i < Ny + 2; i++) {
            for (int j = 0; j < Nx + 2; j++) {
                if (i == 0 || i == Ny + 1 || j == 0 || j == Nx + 1) {
                    phi_1[i][j] = (aE[i][j] * phi[i][j + 1] + aW[i][j] * phi[i][j - 1] +
                                   aN[i][j] * phi[i - 1][j] + aS[i][j] * phi[i + 1][j] + bP[i][j]) / aP[i][j];
                } else {
                    phi_1[i][j] = (aE[i][j] * phi[i][j + 1] + aW[i][j] * phi[i][j - 1] +
                                   aN[i][j] * phi[i - 1][j] + aS[i][j] * phi[i + 1][j] + phi[i][j] * bP[i][j]) / aP[i][
                                      j];
                }

                double diff = fabs(phi_1[i][j] - phi[i][j]);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
            }
        }
        t_count += dt;
        res = maxDiff;

        // Update phi
        for (int i = 0; i < Ny + 2; i++) {
            for (int j = 0; j < Nx + 2; j++) {
                phi[i][j] = phi_1[i][j];
            }
        }
    }
    std::cout << "Steady state reached in " << t_count << " seconds" << std::endl;
}

int main() {

    // Physical Data
    double phi_0 = 0.0;                         // Initial condition in all the domain
    double gamma = 1.0;                         // Assumed equal to 1
    double ratio = 10;                          // To be changed with 10, 1.000, 1.000.000
    double rho = gamma * ratio;                 // Uniform in the domain
    double alpha = 10;

    // Numerical Data
    double Dx = L / Nx, Dy = H / Ny;
    double dt = 1e-4;                           // Time step
    double t_count = 0.0;                       // Initial time count
    double maxRes = 1e-6;                       // Convergence tolerance
    string method;                              // To choose the scheme

    // Start Timer
    auto start = high_resolution_clock::now();

    // Compute Mesh
    MeshDefinition(H, L, W, Dx, Dy);

    // Compute Velocity Field
    VelocityField();

    // Initial Map
    for (auto &row: phi) {
        for (auto &elem: row) {
            elem = phi_0;
        }
    }

    //Choose Interpolation Scheme
    cout << "Choose the method ('CDS' or 'UDS'): ";
    cin >> method;

    // Create an instance of the chosen method
    InternalNodesCalculator* calculator = createCalculator(method);
    if (!calculator) {
        cerr << "Error: Discretization not valid!" << endl;
        return 1;
    }

    // Compute Internal Nodes Coefficients
    calculator->computeCoefficients(rho, dt, gamma);

    // Compute Boundary Conditions
    BoundaryConditions(alpha);

    // Loop Until Steady-State is Reached
    Solver(maxRes, t_count, dt, rho, gamma);

    // Clear memory
    delete calculator;

    // File to check matrices
    ofstream TestFile ("Test.txt");
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            TestFile << xP[i][j] << " ";
        }
        TestFile << endl;
    }

    // File to print values of phi at the outlet
    ofstream TransportedProperty ("TransportedProperty.txt");
    TransportedProperty << "X" << " " << "Ratio = " << ratio << endl;
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            if (i == Ny && xcv[i][j] >= 0)
                TransportedProperty << xcv[i][j] << " " << (phi_1[i][j] + phi_1[i][j+1]) / 2 << endl;
        }
    }

    // File to print phi as a function of x
    ofstream PhiFunctionX ("PhiFx.txt");
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            if (i == Ny)
                PhiFunctionX << xcv[i][j] << " " << (phi_1[i][j] + phi_1[i][j+1]) / 2 << endl;
        }
    }

    // File to print data for Gnuplot 1
    ofstream GnuplotData1 ("GnuplotData1.txt");
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
           GnuplotData1 << xP[i][j] << " " << yP[i][j] << " " << phi_1[i][j] << endl;
        }
        GnuplotData1 << "\n";
    }

    // File to print data for Gnuplot 2
    ofstream GnuplotData2 ("GnuplotData2.txt");
    double ref[11][4] = {
        {0.0, 1.989, 2.0000, 2.000},
        {0.1, 1.402, 1.9990, 2.000},
        {0.2, 1.146, 1.9997, 2.000},
        {0.3, 0.946, 1.9850, 1.999},
        {0.4, 0.775, 1.8410, 1.964},
        {0.5, 0.621, 0.9510, 1.000},
        {0.6, 0.480, 0.1540, 0.036},
        {0.7, 0.349, 0.0010, 0.001},
        {0.8, 0.227, 0.0000, 0.000},
        {0.9, 0.111, 0.0000, 0.000},
        {1.0, 0.000, 0.0000, 0.000}
    };
    // Prints automatically the correct column
    int index = 1;
    if (ratio == 1000)
        index = 2;
    else if (ratio == 1000000)
        index = 3;
    for (int i = 0; i < 11; i++) {
        GnuplotData2 << ref[i][0] << " " << ref[i][index] << std::endl;
    }

    // File to print property distribution
    ofstream PhiDistribution ("PhiDistribution.txt");
    for (int i = 0; i < Ny+2; i++) {
        for (int j = 0; j < Nx+2; j++) {
            PhiDistribution << phi_1[i][j] << " ";
        }
        PhiDistribution << endl;
    }

    TestFile.close();
    TransportedProperty.close();
    GnuplotData1.close();
    GnuplotData2.close();
    PhiDistribution.close();
    PhiFunctionX.close();

    // Generate Gnuplot script 1
    ofstream Gnuplot1("MagnitudeMap.plt");
    Gnuplot1 << "set terminal pngcairo size 1000,500 enhanced font 'Arial,12'\n";
    Gnuplot1 << "set output 'MagnitudeMap_Plot.png'\n";
    Gnuplot1 << "set pm3d map\n";
    Gnuplot1 << "set palette defined ("
        << "0 'dark-blue', "
        << "0.2 'blue', "
        << "0.4 'cyan', "
        << "0.7 'orange', "
        << "1 'yellow')\n";
    Gnuplot1 << "set colorbox\n";
    Gnuplot1 << "set xlabel 'X (m)'\n";
    Gnuplot1 << "set ylabel 'Y (m)'\n";
    Gnuplot1 << "set title 'Steady state of Φ with {/Symbol r}/{/Symbol G} = " << ratio << "' font 'Arial,20'\n";
    Gnuplot1 << "set xrange [-1:1]\n";
    Gnuplot1 << "set yrange [0:1]\n";
    Gnuplot1 << "set autoscale\n";
    Gnuplot1 << "set cbrange [*:*]\n";
    Gnuplot1 << "splot 'GnuplotData1.txt' using 1:2:3 with pm3d notitle\n";
    Gnuplot1.close();

    // Automatically run Gnuplot
    system("gnuplot MagnitudeMap.plt");
    cout << "Gnuplot script 1 executed: 'MagnitudeMap_Plot.png' generated." << endl;

    // Generate Gnuplot script 2
    ofstream Gnuplot2("ComparisonPlot.plt");
    Gnuplot2 << "set terminal pngcairo size 600,400\n";
    Gnuplot2 << "set output 'Comparison_Plot.png'\n";
    Gnuplot2 << "set xlabel 'x'\n";
    Gnuplot2 << "set ylabel 'φ'\n";
    Gnuplot2 << "set title 'Φ = f(x) with {/Symbol r}/{/Symbol G} = " << ratio << "' font 'Arial,20'\n";
    Gnuplot2 << "set key top left\n";
    Gnuplot2 << "set xrange [-1:1]\n";
    Gnuplot2 << "set yrange [0:2]\n";
    Gnuplot2 << "plot 'PhiFx.txt' with lines lw 2 lc rgb 'blue' title 'Result', "
             << "'GnuplotData2.txt' with points pt 7 ps 1.2 lc rgb 'red' title 'Reference'\n";
    Gnuplot2.close();

    // Automatically run Gnuplot
    system("gnuplot ComparisonPlot.plt");
    cout << "Gnuplot script 2 executed: 'Comparison_Plot.png' generated." << endl;

    // Stop Timer and Print Total Duration
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Total Execution Time: " << duration.count() << " seconds" << endl;

    return 0;
}