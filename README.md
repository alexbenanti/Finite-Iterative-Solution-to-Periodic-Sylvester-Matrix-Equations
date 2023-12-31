# Finite-Iterative-Solution-to-Periodic-Sylvester-Matrix-Equations

In "Finite Iterative Solutions to Periodic Sylvester Equations" by Lv and Zhang consider the periodic Sylvester matrix equations (PSMEs), which take two forms: 


$$A_jX_j + X_{j+1}B_j = C_j\hbox{ (1)}$$

$$A_jX_{j+1} + X_jB_j = C_j\hbox{ (2)}$$

where $A_j,B_j,C_j\in\mathbb{R}^{n\times n}$ are the given coefficient/restriction matrices and $X_j\in\mathbb{R}^{n\times n}$ are the solutions that need to be solved. We call (1) the forward periodic Sylvester matrix equations and (2) the backward periodic Sylvester matrix equations. 


The periodic factor is the key distinction between PSMEs and regular Sylvester Matrix Equations ($AX + XB = C$). If a PSME is periodic with period $T$, then, we have:
$$A_{j+T} = A_j,\hbox{ }B_{j+T} = B_j, \hbox{ }C_{j+T} = C_j,\hbox{ }X_{j+T} = X_j$$

It is clear that, analytically, this problem is nontrivial to solve. Trivializing this process is the motivation behind finding a numerical solution to the PSMEs. In this repository, we implemented the algorithms presented by Lv and Zhang in a straightforward manner and applied it to some new PSMEs to verify their results. We also show that the forward and backward PSMEs are equivalent via a numerical example. 
