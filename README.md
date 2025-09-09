# Monte Carlo Particle Tranport code course project
This project implements my personal Monte Carlo particle transport simulation code as part of a course held at Aalto University.

## Requirements
- **Compiler**: GCC with OpenMP support.
- **Operating System**: Linux (tested on WSL system).
- **Libraries**: Standard C99 libraries.

---

## Features
- **Multi-threading**: Parallelized inner loops using OpenMP.
- **Configurable Input**: Parameters can be set via an input file.
- **Progress Bar**: Displays progress during simulations.

---

## Installing and running the project
1. Clone the repository:
   ```
   git clone https://github.com/miikomp/Monte-Carlo-Course
   cd mc-course
   ```
2. Build the project:
   ```
   make all
   ```
3. Run the calculation:
   ```
   ./moca -[KWARGS] [INPUTFILE]
   ```

---

## Example output
```
   Reading input file "input"...

   [NOTE] Using provided seed 12165467978.
   [WARNING] Buffon's needle parameters provided but mode is not set to 2. Ignoring.
   [NOTE] 5 keyword arguments succesfully parsed.
   DONE.

   Approximating pi uisng the 1/4 circle method with 1000 outer, and 1000000 inner iterations...

   [==================================================] 100%

   Mean value = 3.141619 +- 0.000052 [3.141517, 3.141721]

   FOM = 2.91

   ------------------------
   Runtime: 1.2515s
   ------------------------
```
---

## List of commandline arguments
- ```-omp N``` Where N specifies the number of threads

  
   
   
