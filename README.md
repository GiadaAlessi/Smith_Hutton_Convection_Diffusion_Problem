# Numerical Resolution of the Smith-Hutton Convection-Diffusion Problem

## Description
This project presents a numerical study of the Smith-Hutton convection-diffusion problem, where the steady-state transport of a scalar property Φ is investigated under a predefined velocity field. The numerical solution is obtained using the finite volume method with both Upwind Difference Scheme (UDS) and Central Difference Scheme (CDS) to compare their performance in terms of stability and accuracy.

Two versions of the code are provided:
- **`Smith_Hutton.cpp`** — Original version of the code.
- **`Smith_Hutton_Upgraded.cpp`** — Enhanced version with improved modularity and flexibility in selecting numerical schemes (CDS or UDS).

## Code Structure
- **`Smith_Hutton.cpp`** — Main C++ code implementing the original convection-diffusion model.
- **`Smith_Hutton_Upgraded.cpp`** — Enhanced version allowing the user to choose between CDS and UDS schemes.
- **`Smith_Hutton.h`** and **`Smith_Hutton_Upgraded.h`** — Header files for function declarations and constants.
- **`CMakeLists.txt`** — CMake configuration for building the project.
- **`Es2_Smith_Hutton_AlessiGiada.pdf`** — Report presenting the results and analysis in detail.

## Results
The detailed results, including solution profiles for different ρ/Γ ratios and comparative analysis of CDS and UDS schemes, are presented in the report: **`Es2_Smith_Hutton_AlessiGiada.pdf`**.

## How to Run
1. Ensure **Gnuplot** is installed on your system.
2. Compile the code using CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ./Esercizio3       # For the original code
   ./Esercizio3.2     # For the upgraded code
   ```
3. To visualize results, run the Gnuplot script:
   ```bash
   gnuplot MagnitudeMap.plt
   gnuplot ComparisonPlot.plt
   ```

## Author
**Giada Alessi**  
Master in Thermal Engineering  
Universitat Politècnica de Catalunya

