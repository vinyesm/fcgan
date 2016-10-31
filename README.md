# fcgan
Fast column generation for atomic norm regularization

We consider linear regression problems of the form

$$\min_{x}\frac{1}{2}\|Xw-y\|^2+\lambda\gamma_{\mathcal{A}}(w)$$

where $X$ is a design matrix and $\gamma_{\mathcal{A}}$ the Latent Group Lasso(LGL) or the sparse-PCA gauge. 
We also considered the constrained version for LGL

$$\min_{x}\frac{1}{2}\|Xw-y\|^2 \quad \st \quad \Olgl(w)\leq\rho. $$

$$I = \int \rho R^{2} dV$$
