# The standard, constrained and parallel multiobjective EGO algorithm

1. **Main_Standard_Multiobjective_EGO.m**: The unconstrained multiobjective EGO algorithm using EIM (expected improvement matrix) criteria, which is significant cheaper to evaluate than the state-of-the-art multiobjective EI criteria. For detailed description about the EIM criteria, please refer to [1].
2. **Main_Constrained_Multiobjective_EGO.m**: The constrained multiobjective EGO algorithm using CEIM (constrained expected improvement matrix) criteria to solve expensive constrained multiobjective problems.
3. **Main_Parallel_Multiobjective_EGO.m**: The parallel multiobjective EGO algorithm using PEIM (Pseudo Expected Improvement Matrix) criteria, which is able to select multiple candidates in each cycle to evaluate in parallel [2].
4. **Main_Parallel_Constrained_Multiobjective_EGO.m**: The parallel constrained multiobjective EGO algorithm using PCEIM (Pseudo Constrained Expected Improvement Matrix) criteria, which is able to select multiple candidates in each cycle to evaluate in parallel [2].
5. The dace toolbox [3] is used for building the Kriging models in the implementations.
6. The non-dominated sorting method by Yi Cao [4] is used to identify the non-dominated fronts from all the design points
7. The hypervolume indicators are calculated using the faster algorithm of [5] Nicola Beume et al. (2009).
8. Both the EIM and PEIM criteria are maximized by DE [6] algorithm.

### Reference

1. D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill Criteria for Expensive Multiobjective Optimization, IEEE Transactions on Evolutionary Computation, 2017, 21 (6): 956-975.
2. D. Zhan, J. Qian, J. Liu, et al. Pseudo Expected Improvement Matrix Criteria for Parallel Expensive Multi-objective Optimization. In Advances in Structural and Multidisciplinary Optimization: Proceedings of the 12th World Congress of Structural and Multidisciplinary Optimization (WCSMO12), Schumacher, A.,Vietor, T.,Fiebig, S., et al., Eds. Springer International Publishing: Cham, 2018; 175-190.
3. S. N. Lophaven, H. B. Nielsen, and J. Sodergaard, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002. Available at: http://www2.imm.dtu.dk/projects/dace/.
4. http://www.mathworks.com/matlabcentral/fileexchange/17251-pareto-front.
5. N. Beume, C. M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold, On the Complexity of Computing the Hypervolume Indicator, IEEE Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
6. K. Price, R. M. Storn, and J. A. Lampinen, Differential evolution: a practical approach to global optimization: Springer Science & Business Media, 2006. http://www.icsi.berkeley.edu/~storn/code.html
