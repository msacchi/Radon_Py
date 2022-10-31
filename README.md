# RadonPy

## Parabolic Radon Transform

Time domain Parabolic Radon transform solution via the Conjugate Gradient method with implicit forward and adjoint Radon operators.
The Radon domain coefficients $m(\tau,p)$ are found by minimizing the following quadratic cost function

$$J =\| Lm - d \|_2^2 + \mu \| m \|_2^2$$
​
The seismic gather is given by $d(t,x)$, and $L$ provides the Radon forward operator in implicit form. The scalar $\mu$ is the trade-off parameter of the problem.
Notice that CG requires the operator $L$ (Forward) and $L′$ it adjoint (or transpose operator). Check `radon_forward` and `radon_adjoint` in `radon_lib.py`


<img src="Figure_1.png" alt="isolated" width="500"/>



to run the program use: `python main.py`
