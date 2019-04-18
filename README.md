# symplectic-arxiv18a
MATLAB codes for the 2018 arXiv paper discussing synthesis of logical Clifford operators for stabilizer codes.

Paper: https://arxiv.org/abs/1803.06987

Code: https://github.com/nrenga/symplectic-arxiv18a

Copyright (C) 2019  Narayanan Rengaswamy

This project is licensed under the terms of the GNU Affero General Public License (AGPL) v3.0. See LICENSE.md for details.

**Scripts**:

*logical_cliff_ops_642.m*: Execute this script to reproduce results published in Appendix II of the arXiv paper.                     

*logical_cliff_ops.m*: Generic script that can be easily modified to get logical Clifford operators for any stabilizer code. This script also gives a circuit for each obtained symplectic solution. It contains three examples: the [[6,4,2]] CSS code, the [[5,1,3]] perfect code and the [[15,7,3]] Hamming CSS code.


**Functions**:

*find_symp_mat.m*: Algorithm 1 in the paper.

*find_all_symp_mat.m*: Algorithm 2 in the paper. 

*qfind_all_symp_mat.m*: Algorithm 2 in the paper specialized for the application of finding logical Clifford operators for stabilizer codes.                       

*find_logical_cliff.m*: Algorithm 3 in the paper.

*find_symplectic.m*: A function to calculate the symplectic matrix corresponding to a Clifford circuit.

*symp_mat_decompose.m*: Uses Trung Can's algorithm to decompose a symplectic matrix into a product of elementary ones.

*find_circuit.m*: Uses symp_mat_decompose.m to find the decomposition and then calculates the Clifford circuit for each elementary transformation.

*find_unitary*.m: A function to calculate the unitary operator corresponding to a given Clifford circuit.

*calc_conjugate.m*: A function to calculate the action of a Clifford circuit on an input Pauli operator under conjugation.

*symp_inn_pdt.m*: A one-line function to compute symplectic inner product between the corresponding rows of two matrices.

*gf2rref.m*: A function to reduce a binary matrix to its reduced row echelon form over GF(2). Developed by "esromneb" as "g2rref.m" and modified by Narayanan Rengaswamy (available on GitHub Gist).

*gf2matinv.m*: Uses the g2rref.m function to calculate the inverse of a binary matrix over GF(2).

*gflineq.m*: MATLAB in-built function (from Communications System Toolbox) to solve a system of linear equations by Gaussian elimination.

*gflineq_all.m*: A function to solve a system of linear equations by Gaussian elimination and determine all the solutions.

*gf2lu.m*: A function to perform LU decomposition on a binary matrix. Uses the algorithm given by Trefethen and Bau in the book "Numerical Linear Algebra".


**Data**:

*logical_cliff_ops_642.mat*: File containing all solutions listed in Appendix II of the arXiv paper. Also includes a circuit for each solution.

*logical_cliff_ops_513.mat*: File containing all solutions for the [[5,1,3]] code. Also includes a circuit for each solution.

# License

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

