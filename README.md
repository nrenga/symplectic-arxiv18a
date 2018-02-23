# symplectic-arxiv18a
MATLAB codes for the 2018 arXiv paper discussing synthesis of logical Clifford operators for stabilizer codes.


**Scripts**:

*log_cliff_ops_642.m*: Execute this script to reproduce results published in Appendix II of the arXiv paper.                     

*log_cliff_ops.m*: Generic script that can be easily modified to get logical Clifford operators for any stabilizer code.


**Functions**:

*find_symp_mat.m*: Algorithm 1 in the paper.

*find_all_symp_mat.m*: Algorithm 2 in the paper. 

*qfind_all_symp_mat.m*: Algorithm 2 in the paper specialized for the application of finding logical Clifford operators for stabilizer codes.                       

*find_logical_cliff.m*: Algorithm 3 in the paper.

*symp_inn_pdt.m*: A one-line function to compute symplectic inner product between the corresponding rows of two matrices.

*gflineq.m*: MATLAB in-built function (from Communication Toolbox) to solve a system of linear equations by Gaussian elimination.


**Data**:

*log_ops_642.mat*: File containing all solutions listed in Appendix II of the arXiv paper.
