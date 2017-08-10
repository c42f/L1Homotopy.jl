# L1Homotopy

Implementation of the Homotopy method for solving the L1 minimization problem,

    minimize  J(x) = ‖A*x - y‖₂² + λₑ‖x‖₁

    subject to  xᵢ > 0 for all i

The algorithm itself and the connection to other related algorithms such as
matching persuit and Least Angle Regression (LARS) is described in the paper
"Fast solution to L1-norm minimization problems when the solution may be sparse"
by David Donoho and Yaakov Tsaig.  The formulation above is labelled (D_λ) in
the paper, but there are several other equivalent formulations, including two
constrained formulations labelled (P_1) and (L_q) in the paper.  L_q is popular
for statistical model selection and is called Lasso in the statistical
literature.
