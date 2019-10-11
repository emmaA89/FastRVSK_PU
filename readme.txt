This software computes Rational RBF-based partition of unity interpolants.
It makes use of Variably Scaled Kernels (VSKs), to enhance the stability of the 
scheme. Furthermore, thanks to the Deflation Accelerated Conjugate Gradient
(DACG), it turns out to be fast.

Authors:  S. De Marchi, A. Martínez, E. Perracchione, 
                 Universita' di Padova, 
                 Dipartimento di Matematica "Tullio Levi-Civita".

Contacts:  demarchi@math.unipd.it; acalomar@math.unipd.it;
                   emma.perracchione@math.unipd.it

Last modified: 10/10/17.

The script Demo.m provides an example for testing the software
For more details concerning input and output parameters and the usage of 
single functions, see the the corresponding Matlab functions.

The folder RvskPuDacg contains:
1. PU_RBF_VSK_RATIONAL_DAGC. m: this is the main funtion for computing the PU interpolant 
     via VSKs and DACG method.

The folder Distance_Matrix contains:
2. DistanceMatrix_VSK.m: distance matrix via VSKs
3. DistanceMatrixCSRBF.m: distance matrix for CSRBFs 
4. Scale_function.m: scale funtion used by the VSKs

The folder DACG contains:
5. DACG.m: eigenvector via DACG method
6. PowerMethod.m: initial approximation for the DACG scheme
7. ApplyTheta: matrix product for the generalized eigenvalue problem 
8. ApplyGamma: coefficients for the DACG scheme

The folder Data_Structure contains:
9.    IntegerBased_MD_Neighbourhood.m
10. IntegerBased_MD_Structure.m
11. IntegerBased_MD_RangeSearch.m
12. IntegerBased_MD_ContainingQuery.m
13. MakeSDGrid.m: computes muktidimensional grids


*** Remarks ***

- Functions 9--12 construct the partitioning data structure for storing the points 
among the PU patches and are downloadble at http://hdl.handle.net/2318/1559094. See also
R. Cavoretto, A. De Rossi, E. Perracchione, Optimal selection of local approximants in RBF-PU 
interpolation,  to appear on J. Sci. Comput. (2017)].

- Function 13 is provided by the book "Meshfree Approximations Methods with Matlab", 
Gregory E. Fasshauer, World Scientific, 2007.
