# kerdock_design
MATLAB codes for the 2018 paper that uses Kerdock codes and matrices to construct unitary 2-designs and logical unitary 2-designs.

To use the scripts and functions in this folder, clone the parent folder (whose link is given below) and add it to your MATLAB path.

You will also need the MATLAB Communication Systems Toolbox since these resources use Galois field arithmetic (e.g., gfmul, gfdiv).

Code: https://github.com/nrenga/symplectic-arxiv18a

Copyright (C) 2018  Narayanan Rengaswamy

This project is licensed under the terms of the GNU Affero General Public License (AGPL) v3.0. See LICENSE.md for details.

**Scripts**:

*logical_kerdock_design.m*: Use this script to translate a Kerdock design into a logical unitary 2-design for any stabilizer code, and obtain circuits. The script contains an example for the [[6,4,2]] CSS code. 

*circuit_complexity.m*: Use this script to reproduce the figures in the paper that plot gate complexity for the 2-design elements.


**Functions**:

*DGSet.m*: Function to construct the Delsarte-Goethals set DG(m,r).

*kerdock_design.m*: Function to produce a unitary 2-design using Kerdock sets of matrices.

*pauli_group.m*: Function to calculate all Pauli group elements, their circuit representation, unitary matrix and binary vector representations.

*pauli_group_logical.m*: Function to calculate a physical realization for one or several logical Pauli circuit(s), up to multiplication by stabilizers.

*gf2trace.m*: A function to calculate the trace map in a Galois field GF(2^m).


**Data**:

*TG4_no_frob.mat*: File containing the Kerdock unitary 2-design for m = 4 qubits (see kerdock_design.m for details).

*TG4_642_no_frob.mat*: File containing the physical realizations for the logical unitary 2-design on the 4 protected qubits of the [[6,4,2]] CSS code (see logical_kerdock_design.m for details).

*TG_with_frob.mat*: File containing the enlarged Kerdock unitary 2-design for m = 4 qubits, using the Frobenius automorphisms in GF(2^m) (see kerdock_design.m for details).

*TG4_642_with_frob.mat*: File containing the physical realizations for the enlarged logical unitary 2-design on the 4 protected qubits of the [[6,4,2]] CSS code (see logical_kerdock_design.m for details).


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

