# symplectic-arxiv18a
MATLAB codes for the 2018 arXiv paper discussing synthesis of logical Clifford operators for stabilizer codes.


**Scripts**:

*logical_cliff_ops_642.m*: Execute this script to reproduce results published in Appendix II of the arXiv paper.                     

*logical_cliff_ops.m*: Generic script that can be easily modified to get logical Clifford operators for any stabilizer code. This script also gives a circuit for each obtained symplectic solution.


**Functions**:

*find_symp_mat.m*: Algorithm 1 in the paper.

*find_all_symp_mat.m*: Algorithm 2 in the paper. 

*qfind_all_symp_mat.m*: Algorithm 2 in the paper specialized for the application of finding logical Clifford operators for stabilizer codes.                       

*find_logical_cliff.m*: Algorithm 3 in the paper.

*symp_mat_decompose.m*: Uses Trung Can's algorithm to decompose a symplectic matrix into a product of elementary ones.

*find_circuit.m*: Uses symp_mat_decompose.m to find the decomposition and then calculates the circuit for each elementary transformation.

*symp_inn_pdt.m*: A one-line function to compute symplectic inner product between the corresponding rows of two matrices.

*g2rref.m*: A function to reduce a binary matrix to its reduced row echelon form over GF(2). Developed by "esromneb" and modified by Narayanan Rengaswamy (available on GitHub Gist).

*g2matinv.m*: Uses the g2rref.m function to calculate the inverse of a binary matrix over GF(2).

*gflineq.m*: MATLAB in-built function (from Communication Toolbox) to solve a system of linear equations by Gaussian elimination.


**Data**:

*logical_ops_642.mat*: File containing all solutions listed in Appendix II of the arXiv paper. Also includes a circuit for each solution.

*logical_ops_513.mat*: File containing all solutions for the [[5,1,3]] code. Also includes a circuit for each solution.
