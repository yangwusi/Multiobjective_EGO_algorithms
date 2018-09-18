# The parallel multiobjective EGO algorithm

1. The multiobjective EGO algorithm using EIM (expected improvement matrix) criteria, which is significant cheaper to evaluate than the state-of-the-art multiobjective EI criteria. For detailed description about the EIM criteria, please refer to [1].
2. The parallel multiobjective EGO algorithm using PEIM (Pseudo Expected Improvement Matrix) criteria, which is able to select multiple candidates in each cycle to evaluate in parallel [2].
3. The dace toolbox [3] is used for building the Kriging models in the implementations.
4. The non-dominated sorting method by Yi Cao [4] is used to identify the non-dominated fronts from all the design points
5. The hypervolume indicators are calculated using the faster algorithm of [5] Nicola Beume et al. (2009).
6. You need MATLAB 2016b or newer version and the global optimization toolbox to run the codes.

### Reference

1. D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill Criteria for Expensive Multiobjective Optimization, IEEE Transactions on Evolutionary Computation, 2017, 21 (6): 956-975.
2. D. Zhan ,J. Qian ,J. Liu , et al. Pseudo Expected Improvement Matrix Criteria for Parallel Expensive Multi-objective Optimization. In Advances in Structural and Multidisciplinary Optimization: Proceedings of the 12th World Congress of Structural and Multidisciplinary Optimization (WCSMO12), Schumacher, A.,Vietor, T.,Fiebig, S., et al., Eds. Springer International Publishing: Cham, 2018; 175-190.
3. Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002. Available at: http://www2.imm.dtu.dk/projects/dace/.
4. http://www.mathworks.com/matlabcentral/fileexchange/17251-pareto-front.
5. N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold, On the Complexity of Computing the Hypervolume Indicator, IEEE Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
