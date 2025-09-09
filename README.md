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

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/miikomp/Monte-Carlo-Course
   cd mc-course
2. Build the project:
   ```bash
   make all
3. Run the calculation:
   ```bash
   ./moca -[KWARGS] [INPUTFILE]

---

## List of commandline arguments
1. Set the number of threads to use:
   - **-omp N** - Where N specifies the number of threads

  
   
   
