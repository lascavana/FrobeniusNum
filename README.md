# FrobeniusNum
Calculates the Frobenius number for a given vector a.

## Input
Input file must contain the components $a_i$ of a vector $a$ with $a_i>0$ and $a_i \leq a_{i+1}$ for all $i$.

## Output
The code will print the Frobenius number. It will also output an lp file with the corresponding optimization problem ( $a\cdot x = F$ ).

## Requirements
- [Gurobi](https://www.gurobi.com/)
- [NTL](https://libntl.org/)

## Instructions
Compilation and running instructions can be found in ```run.sh```.


