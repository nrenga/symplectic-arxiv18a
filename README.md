# symplectic-arxiv18a
MATLAB codes for the 2018 arXiv paper discussing synthesis of logical Clifford operators for stabilizer codes.


**Scripts**:

*log_cliff_ops_642.m*: Execute this script to reproduce results published in Appendix II of the arXiv paper.
                     A function inside this script uses the GF structure (from Communication Toolbox) to invert binary matrix.

*log_cliff_ops.m*: Generic script that can be easily modified to get logical Clifford operators for any stabilizer code.


**Functions**:

*find_symp_mat.m*: Algorithm 1 in the paper.

*find_all_symp_mat.m*: Algorithm 2 in the paper. Uses the GF structure (from Communication Toolbox) to invert binary matrix.

*qfind_all_symp_mat.m*: Algorithm 2 in the paper specialized for the application of finding logical Clifford operators for stabilizer codes. 
                      Uses the GF structure (from Communication Toolbox) to invert binary matrix.

*find_logical_cliff.m*: Algorithm 3 in the paper, but does not obtain the circuit yet; only gets all symplectic solutions.

*symp_inn_pdt.m*: A one-line function to compute symplectic inner product between the corresponding rows of two matrices.

*gflineq.m*: MATLAB in-built function (from Communication Toolbox) to solve a system of linear equations by Gaussian elimination.

*gf2dec.m*: User-defined function (written by Murad Qahwash) to convert a GF matrix into a decimal matrix.


**Data**:

*log_ops_642.mat*: File containing all solutions listed in Appendix II of the arXiv paper.
