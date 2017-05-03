# The_expected_improvement_matrix_criteria

1. The multiobjective EGO algorithm using EIM (expected improvement matrix)-based criteria, which is significant cheaper-to-evaluate than the state-of-the-art multiobjective EI criteria. For detailed description about the EIM criteria, please refer to [1].
2. The dace toolbox [2] is used for building the Kriging models in the implementations.
3. The non-dominated sorting method by Yi Cao [3] is used to identify the non-dominated fronts from all the design points
4. The hypervolume indicators are calculated using the faster algorithm of [4] Nicola Beume et al. (2009).
5. You need MATLAB 2016b or newer version and the global optimization toolbox to run the codes.

### Reference

1. D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill Criteria for Expensive Multiobjective Optimization, IEEE Transactions on Evolutionary Computation, DOI: 10.1109/TEVC.2017.2697503
2. Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002. Available at: http://www2.imm.dtu.dk/projects/dace/.
3. http://www.mathworks.com/matlabcentral/fileexchange/17251-pareto-front.
4. N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold, On the Complexity of Computing the Hypervolume Indicator, IEEE Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
