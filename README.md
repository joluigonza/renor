
Instructions to run the code:

1) Navigate to the folder in the code editor of your choice 
2) Enter pkg mode by typing ]
3) Type: activate .  (the space and . are needed)
4) Press enter
5) Now run the file final_script_gen_m_d.jl

Alternatively, simply run the notebook final_notebook_gen_m_d.ipynb

Below is the list of parameters that you need to change at the beginning of the code, 
either in final_script_gen_m_d.jl or in final_notebook_gen_m_d.ipynb, in order to reproduce the different proofs
from the paper (Theorem 1.5, Theorem 1.6, and Theorem 7.4). The values used for the proofs presented in the paper are listed in Appendix D.

pre # Precision of the computations

m  # Renormalization order

d  # Degree of the fixed point

ver # index of fixed point, starting from m=5 there is more that one

rho # Bernstein ellipse radius

rstar  # rstar radius in the contraction proof for the fixed point

K0  # Chebyshev order (the actual order is in fact 2K)


In order to reproduce the 500 digit proof for m=2 (Theorem 7.2), you need to run final_script_m2_long.jl

