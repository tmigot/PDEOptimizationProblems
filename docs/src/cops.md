# Simulations on COPS 3 test set

We run [DCISolver.jl](https://github.com/JuliaSmoothOptimizers/DCISolver.jl) and [NLPModelsIpopt.jl](https://github.com/JuliaSmoothOptimizers/DCISolver.jl) and [NLPModelsKnitro.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsKnitro.jl) on an infinite-dimensional implementation of test problems from [COPS 3.0](https://www.mcs.anl.gov/~more/cops/) test set.

## DCILDL
|                                         name |   nvar |   ncon |      status | objective |   opt_val |  time (s) |     #f |     #c | dual_feas |      feas |
|----------------------------------------------|--------|--------|-------------|-----------|-----------|-----------|--------|--------|-----------|-----------|
|                   Isometrization of α-pinene |    505 |    500 |    max_time |  0.00e+00 |  1.99e+01 |  8.49e+00 |      2 |      2 |  0.00e+00 |  2.84e+01 |
|                   Isometrization of α-pinene |   1005 |   1000 |    max_time |  0.00e+00 |  1.99e+01 |  4.05e+00 |      2 |      2 |  0.00e+00 |  2.84e+01 |
|                   Isometrization of α-pinene |   2005 |   2000 |    max_time |  0.00e+00 |  1.99e+01 |  4.42e+00 |      2 |      2 |  0.00e+00 |  2.84e+01 |
|                              Journal Bearing |   2600 |      0 |   exception |       Inf | -1.55e-01 |       Inf |      0 |      0 |       Inf |       Inf |
|                              Journal Bearing |   5775 |      0 |   exception |       Inf | -1.55e-01 |       Inf |      0 |      0 |       Inf |       Inf |
|                              Journal Bearing |  10200 |      0 |   exception |       Inf | -1.55e-01 |       Inf |      0 |      0 |       Inf |       Inf |
|                              Catalyst Mixing |    602 |    402 |    max_time |  0.00e+00 | -4.81e-02 |  9.96e+01 |      1 |      1 |  1.41e-01 |  2.50e-01 |
|                              Catalyst Mixing |   1202 |    802 |    max_time |  0.00e+00 | -4.81e-02 |  3.94e+00 |      1 |      1 |  9.98e-02 |  2.50e-01 |
|                              Catalyst Mixing |   2402 |   1602 |    max_time |  0.00e+00 | -4.81e-02 |  1.07e+03 |      2 |      2 |  7.07e-02 |  2.50e-01 |
|                            Flow in a Channel |    800 |    800 |    max_time |  0.00e+00 |  1.00e+00 |  5.14e+00 |      2 |      2 |  0.00e+00 |  2.20e-01 |
|                            Flow in a Channel |   1600 |   1600 |    max_time |  0.00e+00 |  1.00e+00 |  2.57e+00 |      2 |      2 |  0.00e+00 |  2.20e-01 |
|                            Flow in a Channel |   3200 |   3200 |    max_time |  0.00e+00 |  1.00e+00 |  2.78e+00 |      2 |      2 |  0.00e+00 |  2.20e-01 |
|  Transition States for the Dirichlet Problem |      9 |      0 | first_order |  0.00e+00 |  1.94e-06 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     19 |      0 | first_order |  0.00e+00 |  1.71e-02 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     39 |      0 | first_order |  0.00e+00 |  3.29e-02 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|                Catalytic Cracking of Gas Oil |    205 |    202 |    max_time |  8.21e-02 |  5.24e-03 |  6.83e+01 |      1 |      1 |  4.19e-02 |  0.00e+00 |
|                Catalytic Cracking of Gas Oil |    405 |    402 |    max_time |  8.11e-02 |  5.24e-03 |  1.11e+00 |      1 |      1 |  2.97e-02 |  0.00e+00 |
|                Catalytic Cracking of Gas Oil |    805 |    802 |    max_time |  8.11e-02 |  5.24e-03 |  2.21e+00 |      1 |      1 |  2.10e-02 |  0.00e+00 |
|                                  Hang Glider |    698 |    498 |    max_time | -9.50e+02 |  1.25e+03 |  4.51e+03 |      2 |      2 |  4.75e-15 |  8.79e+01 |
|                                  Hang Glider |   1398 |    998 |    max_time | -9.50e+02 |  1.25e+03 |  8.86e+03 |      2 |      2 |  6.76e-15 |  6.44e+01 |
|                                  Hang Glider |   2798 |   1998 |    max_time | -9.50e+02 |  1.25e+03 |  1.77e+04 |      2 |      2 |  5.05e-09 |  4.66e+01 |
|      Transition States for the Henon Problem |      9 |      0 | first_order |  0.00e+00 |  7.22e+00 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     19 |      0 | first_order |  0.00e+00 |  7.52e+01 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     39 |      0 | first_order |  0.00e+00 |  1.26e+02 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |      9 |      0 | first_order |  0.00e+00 |  8.49e+00 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     19 |      0 | first_order |  0.00e+00 |  9.11e+00 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     39 |      0 | first_order |  0.00e+00 |  9.28e+00 |  0.00e+00 |      1 |      1 |  0.00e+00 |  0.00e+00 |
|                     Methanol to Hydrocarbons |    308 |    303 |   exception |       Inf |  9.02e-03 |       Inf |      0 |      0 |       Inf |       Inf |
|                     Methanol to Hydrocarbons |    608 |    603 |   exception |       Inf |  9.02e-03 |       Inf |      0 |      0 |       Inf |       Inf |
|                     Methanol to Hydrocarbons |   1208 |   1203 |   exception |       Inf |  9.02e-03 |       Inf |      0 |      0 |       Inf |       Inf |
|                Minimal Surface with Obstacle |   5002 |   2401 |    max_time |  2.33e+00 |  2.51e+00 |  1.94e+01 |      1 |      1 |  6.85e-02 |  7.34e-03 |
|                Minimal Surface with Obstacle |  11252 |   5476 |    max_time |  2.32e+00 |  2.51e+00 |  1.59e+00 |      2 |      2 |  5.29e-02 |  5.08e-03 |
|                Minimal Surface with Obstacle |  20002 |   9801 |    max_time |  2.33e+00 |  2.51e+00 |  2.33e+00 |      2 |      2 |  4.44e-02 |  3.84e-03 |
|                                    Robot Arm |   1200 |    597 |    max_time |  0.00e+00 |  9.14e+00 |  3.74e+00 |      2 |      2 |  0.00e+00 |  2.01e+04 |
|                                    Robot Arm |   2400 |   1197 |    max_time |  0.00e+00 |  9.14e+00 |  3.87e+00 |      2 |      2 |  0.00e+00 |  4.02e+04 |
|                                    Robot Arm |   4800 |   2397 |    max_time |  0.00e+00 |  9.14e+00 |  4.08e+00 |      2 |      2 |  0.00e+00 |  8.04e+04 |
|                               Goddard Rocket |   2400 |   1600 |    max_time | -1.00e+00 |  1.01e+00 |  4.83e+04 |      2 |      2 |  1.85e-02 |  2.64e-01 |
|                               Goddard Rocket |   4800 |   3200 |    max_time | -1.00e+00 |  1.01e+00 |  5.78e+04 |      2 |      2 |  1.32e-02 |  1.90e-01 |
|                               Goddard Rocket |   9600 |   6400 |    max_time | -1.00e+00 |  1.01e+00 |  1.23e+05 |      2 |      2 |  9.39e-03 |  1.69e-01 |
|                            Particle Steering |    800 |    400 |    max_time |  0.00e+00 |  5.55e-01 |  5.91e+01 |      2 |      2 |  0.00e+00 |  1.35e+01 |
|                            Particle Steering |   1600 |    800 |    max_time |  0.00e+00 |  5.55e-01 |  1.80e+02 |      2 |      2 |  0.00e+00 |  9.56e+00 |
|                            Particle Steering |   3200 |   1600 |    max_time |  0.00e+00 |  5.55e-01 |  6.38e+02 |      2 |      2 |  0.00e+00 |  6.79e+00 |
|                      Elastic-Plastic Torsion |   7802 |   7802 |    max_time |  0.00e+00 | -4.18e-01 |  5.58e+01 |      1 |      1 |  5.72e-02 |  9.93e-03 |
|                      Elastic-Plastic Torsion |  17327 |  17327 |    max_time |  0.00e+00 | -4.18e-01 |  2.22e+00 |      1 |      1 |  3.82e-02 |  6.64e-03 |
|                      Elastic-Plastic Torsion |  30602 |  30602 |    max_time |  0.00e+00 | -4.18e-01 |  2.71e+00 |      1 |      1 |  2.87e-02 |  4.98e-03 |
## knitro
|                                         name |   nvar |   ncon |      status | objective |   opt_val |  time (s) |     #f |     #c | dual_feas |      feas |
|----------------------------------------------|--------|--------|-------------|-----------|-----------|-----------|--------|--------|-----------|-----------|
|                   Isometrization of α-pinene |    505 |    500 |    max_time |  0.00e+00 |  1.99e+01 |  4.33e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                   Isometrization of α-pinene |   1005 |   1000 |    max_time |  0.00e+00 |  1.99e+01 |  1.10e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                   Isometrization of α-pinene |   2005 |   2000 |    max_time |  0.00e+00 |  1.99e+01 |  1.26e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                              Journal Bearing |   2600 |      0 |    max_time |  1.20e+01 | -1.55e-01 |  9.39e+00 |      1 |      0 |  5.88e-01 |  0.00e+00 |
|                              Journal Bearing |   5775 |      0 |    max_time |  5.06e+00 | -1.55e-01 |  5.31e+00 |      2 |      0 |  2.31e-01 |  0.00e+00 |
|                              Journal Bearing |  10200 |      0 |    max_time | -1.05e+00 | -1.55e-01 |  1.09e+00 |      5 |      0 |  3.17e-04 |  0.00e+00 |
|                              Catalyst Mixing |    602 |    402 |    max_time |  0.00e+00 | -4.81e-02 |  2.82e+00 |      2 |      2 |  1.00e-02 |  2.49e-01 |
|                              Catalyst Mixing |   1202 |    802 |    max_time |  0.00e+00 | -4.81e-02 |  2.56e+00 |      2 |      2 |  5.00e-03 |  2.49e-01 |
|                              Catalyst Mixing |   2402 |   1602 |    max_time |  0.00e+00 | -4.81e-02 |  3.60e+00 |      2 |      2 |  2.50e-03 |  2.50e-01 |
|                            Flow in a Channel |    800 |    800 | first_order |  0.00e+00 |  1.00e+00 |  5.12e+01 |      2 |      2 |  0.00e+00 |  1.99e-13 |
|                            Flow in a Channel |   1600 |   1600 | first_order |  0.00e+00 |  1.00e+00 |  1.58e+00 |      2 |      2 |  0.00e+00 |  4.44e-16 |
|                            Flow in a Channel |   3200 |   3200 | first_order |  0.00e+00 |  1.00e+00 |  1.55e+00 |      2 |      2 |  0.00e+00 |  6.66e-16 |
|  Transition States for the Dirichlet Problem |      9 |      0 | first_order |  0.00e+00 |  1.94e-06 |  5.16e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     19 |      0 | first_order |  0.00e+00 |  1.71e-02 |  1.99e-01 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     39 |      0 | first_order |  0.00e+00 |  3.29e-02 |  4.02e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|                Catalytic Cracking of Gas Oil |    205 |    202 |    max_time |  7.94e-02 |  5.24e-03 |  5.73e+01 |      2 |      2 |  7.01e-03 |  4.86e-07 |
|                Catalytic Cracking of Gas Oil |    405 |    402 |    max_time |  7.93e-02 |  5.24e-03 |  1.45e+00 |      2 |      2 |  5.73e-03 |  3.20e-07 |
|                Catalytic Cracking of Gas Oil |    805 |    802 |    max_time |  4.65e-02 |  5.24e-03 |  2.88e+00 |      2 |      2 |  3.14e-03 |  7.89e-05 |
|                                  Hang Glider |    698 |    498 |    max_time | -9.50e+02 |  1.25e+03 |  2.77e+01 |      1 |      1 |  4.44e-16 |  9.73e+00 |
|                                  Hang Glider |   1398 |    998 |    max_time | -9.50e+02 |  1.25e+03 |  2.74e+01 |      1 |      1 |  6.66e-16 |  4.86e+00 |
|                                  Hang Glider |   2798 |   1998 |    max_time | -9.50e+02 |  1.25e+03 |  2.77e+01 |      1 |      1 |  6.66e-16 |  2.43e+00 |
|      Transition States for the Henon Problem |      9 |      0 | first_order |  0.00e+00 |  7.22e+00 |  2.82e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     19 |      0 | first_order |  0.00e+00 |  7.52e+01 |  2.74e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     39 |      0 | first_order |  0.00e+00 |  1.26e+02 |  5.00e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |      9 |      0 | first_order |  0.00e+00 |  8.49e+00 |  2.56e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     19 |      0 | first_order |  0.00e+00 |  9.11e+00 |  4.14e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     39 |      0 | first_order |  0.00e+00 |  9.28e+00 |  3.79e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|                     Methanol to Hydrocarbons |    308 |    303 |   exception |  2.33e-01 |  9.02e-03 |  6.73e-01 |      1 |      1 | 1.80e+308 |  0.00e+00 |
|                     Methanol to Hydrocarbons |    608 |    603 |   exception |  2.33e-01 |  9.02e-03 |  5.03e-01 |      1 |      1 | 1.80e+308 |  0.00e+00 |
|                     Methanol to Hydrocarbons |   1208 |   1203 |   exception |  2.33e-01 |  9.02e-03 |  6.91e-01 |      1 |      1 | 1.80e+308 |  0.00e+00 |
|                Minimal Surface with Obstacle |   5002 |   2401 |    max_time |  5.78e+00 |  2.51e+00 |  1.86e+01 |      3 |      3 |  5.09e-02 |  4.30e-04 |
|                Minimal Surface with Obstacle |  11252 |   5476 |    max_time |  6.40e+00 |  2.51e+00 |  1.70e+00 |      3 |      3 |  3.33e-02 |  1.85e-04 |
|                Minimal Surface with Obstacle |  20002 |   9801 |    max_time |  8.47e+00 |  2.51e+00 |  2.43e+00 |      3 |      3 |  2.65e-02 |  1.23e-04 |
|                                    Robot Arm |   1200 |    597 |    max_time |  0.00e+00 |  9.14e+00 |  2.04e+00 |      2 |      2 |  0.00e+00 |  8.31e+03 |
|                                    Robot Arm |   2400 |   1197 |    max_time |  0.00e+00 |  9.14e+00 |  1.79e+00 |      2 |      2 |  0.00e+00 |  1.66e+04 |
|                                    Robot Arm |   4800 |   2397 |    max_time |  0.00e+00 |  9.14e+00 |  2.33e+00 |      2 |      2 |  0.00e+00 |  3.32e+04 |
|                               Goddard Rocket |   2400 |   1600 |    max_time | -1.00e+00 |  1.01e+00 |  3.13e+00 |      3 |      3 |  2.50e-03 |  4.61e-02 |
|                               Goddard Rocket |   4800 |   3200 |    max_time | -1.00e+00 |  1.01e+00 |  3.16e+00 |      3 |      3 |  1.25e-03 |  2.30e-02 |
|                               Goddard Rocket |   9600 |   6400 |    max_time | -1.00e+00 |  1.01e+00 |  3.36e+00 |      3 |      3 |  6.25e-04 |  1.15e-02 |
|                            Particle Steering |    800 |    400 | first_order |  0.00e+00 |  5.55e-01 |  1.80e+01 |      2 |      2 |  0.00e+00 |  7.62e-12 |
|                            Particle Steering |   1600 |    800 | first_order |  0.00e+00 |  5.55e-01 |  3.54e-01 |      2 |      2 |  0.00e+00 |  1.66e-11 |
|                            Particle Steering |   3200 |   1600 | first_order |  0.00e+00 |  5.55e-01 |  5.49e-01 |      2 |      2 |  0.00e+00 |  2.86e-11 |
|                      Elastic-Plastic Torsion |   7802 |   7802 |    max_time |  0.00e+00 | -4.18e-01 |  1.93e+00 |      2 |      2 |  8.06e-04 |  5.92e-04 |
|                      Elastic-Plastic Torsion |  17327 |  17327 |    max_time |  0.00e+00 | -4.18e-01 |  1.19e+01 |      2 |      2 |  3.60e-04 |  2.64e-04 |
|                      Elastic-Plastic Torsion |  30602 |  30602 |    max_time |  0.00e+00 | -4.18e-01 |  2.72e+01 |      2 |      2 |  2.14e-04 |  1.49e-04 |
## ipopt
|                                         name |   nvar |   ncon |      status | objective |   opt_val |  time (s) |     #f |     #c | dual_feas |      feas |
|----------------------------------------------|--------|--------|-------------|-----------|-----------|-----------|--------|--------|-----------|-----------|
|                   Isometrization of α-pinene |    505 |    500 |    max_time |  0.00e+00 |  1.99e+01 |  1.67e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                   Isometrization of α-pinene |   1005 |   1000 |    max_time |  0.00e+00 |  1.99e+01 |  1.83e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                   Isometrization of α-pinene |   2005 |   2000 |    max_time |  0.00e+00 |  1.99e+01 |  2.02e+00 |      1 |      1 |  0.00e+00 |  5.00e+01 |
|                              Journal Bearing |   2600 |      0 |    max_time |  0.00e+00 | -1.55e-01 |  1.06e+00 |      9 |      0 |  2.19e-05 |  0.00e+00 |
|                              Journal Bearing |   5775 |      0 |    max_time |  0.00e+00 | -1.55e-01 |  1.08e+00 |      4 |      0 |  5.14e-02 |  0.00e+00 |
|                              Journal Bearing |  10200 |      0 |    max_time |  0.00e+00 | -1.55e-01 |  1.07e+00 |      3 |      0 |  1.74e-01 |  0.00e+00 |
|                              Catalyst Mixing |    602 |    402 |    max_time |  0.00e+00 | -4.81e-02 |  2.58e+00 |      1 |      1 |  1.00e-02 |  2.50e-01 |
|                              Catalyst Mixing |   1202 |    802 |    max_time |  0.00e+00 | -4.81e-02 |  2.53e+00 |      1 |      1 |  5.00e-03 |  2.50e-01 |
|                              Catalyst Mixing |   2402 |   1602 |    max_time |  0.00e+00 | -4.81e-02 |  5.93e+00 |      1 |      1 |  2.50e-03 |  2.50e-01 |
|                            Flow in a Channel |    800 |    800 |    max_time |  0.00e+00 |  1.00e+00 |  2.68e+01 |      2 |      2 |  2.00e-02 |  3.32e-08 |
|                            Flow in a Channel |   1600 |   1600 |    max_time |  0.00e+00 |  1.00e+00 |  1.52e+00 |      2 |      2 |  2.00e-02 |  6.64e-08 |
|                            Flow in a Channel |   3200 |   3200 |    max_time |  0.00e+00 |  1.00e+00 |  1.15e+00 |      2 |      2 |  2.00e-02 |  1.33e-07 |
|  Transition States for the Dirichlet Problem |      9 |      0 | first_order |  0.00e+00 |  1.94e-06 |  4.90e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     19 |      0 | first_order |  0.00e+00 |  1.71e-02 |  1.48e-01 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|  Transition States for the Dirichlet Problem |     39 |      0 | first_order |  0.00e+00 |  3.29e-02 |  4.40e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|                Catalytic Cracking of Gas Oil |    205 |    202 |    max_time |  0.00e+00 |  5.24e-03 |  2.90e+00 |      2 |      2 |  1.03e-02 |  5.89e-06 |
|                Catalytic Cracking of Gas Oil |    405 |    402 |    max_time |  0.00e+00 |  5.24e-03 |  1.42e+00 |      2 |      2 |  1.47e-02 |  3.88e-06 |
|                Catalytic Cracking of Gas Oil |    805 |    802 |    max_time |  0.00e+00 |  5.24e-03 |  1.88e+00 |      2 |      2 |  3.79e-02 |  1.84e-06 |
|                                  Hang Glider |    698 |    498 |    max_time |  0.00e+00 |  1.25e+03 |  3.71e+01 |      1 |      1 |  1.00e+00 |  9.73e+00 |
|                                  Hang Glider |   1398 |    998 |    max_time |  0.00e+00 |  1.25e+03 |  3.77e+01 |      1 |      1 |  1.00e+00 |  4.86e+00 |
|                                  Hang Glider |   2798 |   1998 |    max_time |  0.00e+00 |  1.25e+03 |  3.86e+01 |      1 |      1 |  1.00e+00 |  2.43e+00 |
|      Transition States for the Henon Problem |      9 |      0 | first_order |  0.00e+00 |  7.22e+00 |  3.60e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     19 |      0 | first_order |  0.00e+00 |  7.52e+01 |  4.00e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|      Transition States for the Henon Problem |     39 |      0 | first_order |  0.00e+00 |  1.26e+02 |  1.36e-01 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |      9 |      0 | first_order |  0.00e+00 |  8.49e+00 |  3.80e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     19 |      0 | first_order |  0.00e+00 |  9.11e+00 |  4.00e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
| Transition States for the Lane-Emden Problem |     39 |      0 | first_order |  0.00e+00 |  9.28e+00 |  3.90e-02 |      1 |      0 |  0.00e+00 |  0.00e+00 |
|                     Methanol to Hydrocarbons |    308 |    303 |     unknown |  0.00e+00 |  9.02e-03 |  1.67e+01 |      0 |      1 |       Inf |       Inf |
|                     Methanol to Hydrocarbons |    608 |    603 |     unknown |  0.00e+00 |  9.02e-03 |  1.20e+00 |      0 |      1 |       Inf |       Inf |
|                     Methanol to Hydrocarbons |   1208 |   1203 |     unknown |  0.00e+00 |  9.02e-03 |  1.03e+01 |      0 |      1 |       Inf |       Inf |
|                Minimal Surface with Obstacle |   5002 |   2401 |    max_time |  0.00e+00 |  2.51e+00 |  1.29e+00 |      1 |      1 |  1.00e+00 |  3.84e-04 |
|                Minimal Surface with Obstacle |  11252 |   5476 |    max_time |  0.00e+00 |  2.51e+00 |  1.14e+01 |      1 |      1 |  1.00e+00 |  1.74e-04 |
|                Minimal Surface with Obstacle |  20002 |   9801 |    max_time |  0.00e+00 |  2.51e+00 |  4.18e+01 |      1 |      1 |  1.00e+00 |  9.85e-05 |
|                                    Robot Arm |   1200 |    597 |    max_time |  0.00e+00 |  9.14e+00 |  2.07e+00 |      1 |      1 |  0.00e+00 |  1.41e+04 |
|                                    Robot Arm |   2400 |   1197 |    max_time |  0.00e+00 |  9.14e+00 |  2.08e+00 |      1 |      1 |  0.00e+00 |  2.83e+04 |
|                                    Robot Arm |   4800 |   2397 |    max_time |  0.00e+00 |  9.14e+00 |  1.83e+00 |      1 |      1 |  0.00e+00 |  5.66e+04 |
|                               Goddard Rocket |   2400 |   1600 |    max_time |  0.00e+00 |  1.01e+00 |  2.92e+00 |      1 |      1 |  1.19e+00 |  1.00e-02 |
|                               Goddard Rocket |   4800 |   3200 |    max_time |  0.00e+00 |  1.01e+00 |  2.92e+00 |      1 |      1 |  1.00e+00 |  1.00e-02 |
|                               Goddard Rocket |   9600 |   6400 |    max_time |  0.00e+00 |  1.01e+00 |  2.36e+00 |      1 |      1 |  1.00e+00 |  1.00e-02 |
|                            Particle Steering |    800 |    400 |    max_time |  0.00e+00 |  5.55e-01 |  1.03e+01 |      2 |      2 |  0.00e+00 |  7.62e-12 |
|                            Particle Steering |   1600 |    800 |    max_time |  0.00e+00 |  5.55e-01 |  1.06e+00 |      4 |      4 |  0.00e+00 |  6.15e-12 |
|                            Particle Steering |   3200 |   1600 |    max_time |  0.00e+00 |  5.55e-01 |  1.07e+00 |      4 |      4 |  0.00e+00 |  5.94e-12 |
|                      Elastic-Plastic Torsion |   7802 |   7802 |    max_time |  0.00e+00 | -4.18e-01 |  7.97e+01 |      2 |      2 |  2.00e-07 |  1.01e-05 |
|                      Elastic-Plastic Torsion |  17327 |  17327 |    max_time |  0.00e+00 | -4.18e-01 |  3.60e+02 |      2 |      2 |  2.00e-07 |  5.46e-06 |
|                      Elastic-Plastic Torsion |  30602 |  30602 |    max_time |  0.00e+00 | -4.18e-01 |  1.14e+00 |      1 |      1 |  1.00e+00 |  5.05e-05 |
