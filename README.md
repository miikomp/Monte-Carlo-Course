# Monte Carlo Particle Transport code course project
This project implements my personal Monte Carlo particle transport simulation code as part of a course held at Aalto University. The current implementation includes the groundwork and two methods for estimating pi. 

## Requirements
- **Compiler**: GCC with Open-MP support.
- **Operating System**: Linux (tested on WSL system).
- **Libraries**: Standard C99 libraries.

---

## Features
- **Multi-threading**: Parallelisation using Open-MP.
- **Configurable Input**: Parameters are set via an input file.
- **Progress Bar**: Displays progress during simulations.

---

## Installing and running the project
1. Clone the repository:
   ```Docker
   git clone https://github.com/miikomp/Monte-Carlo-Course
   cd mc-course
   ```
2. Build the project:
   ```Docker
   make all
   ```
3. Create an input file (or use the 'input' in the repo)
   ```Docker
   # Example input file for the "moca" Monte Carlo code

   # Set the seed (if not set a random seed is used)
   seed 12165467978

   # Set number of inner and outer iterations, or number of histories in a 
   pop 10000000 10000

   # Set parameters. All are optional
   # Set calculation mode, 1: Quarter circle pi, 2: Buffon's needle, 3: Check, default 0: Transport
   set mode 2

   # Length of needle
   set needle 1.0

   # Distance between lines
   set lines 1.0
   ```
4. Run the calculation:
   ```Docker
   ./moca -[KWARGS] [INPUTFILE]
   ```

---

## List of commandline arguments
- ```-omp N``` Where N specifies the number of threads to use

---

## Example output
```
------------------------
  Processing input
------------------------

Reading input file "input"...

[NOTE] Using provided seed 12165467978.

DONE.

[NOTE] 5 keyword arguments succesfully parsed.

------------------------
  Starting simulation
------------------------

Approximating pi using Buffon's Needle with 10000 outer, and 10000000 inner iterations...

[==================================================] 100%

Mean value = 3.141591 +/- 0.000007 [3.141577, 3.141606]

FOM = 4.69

Kolmogorov-Smirnov test for normality:
          K-S D: 0.005487
  Standard crit: 0.013600: TRUE
Lilliefors crit: 0.008860: TRUE


------------------------
  Runtime: 37.5004s
------------------------
```


  
   
   
