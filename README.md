# MANOVABNPTest

![alt text](images/model.jpg)

where ![formula](https://render.githubusercontent.com/render/math?math=N) is the sample size, ![formula](https://render.githubusercontent.com/render/math?math=J) is the number of treatment groups, <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\boldsymbol{y}_i%20\in%20\mathbb{R}^D" /> is the response for the $i$th observation, <img alt="formula" src="https://render.githubusercontent.com/render/math?math=x_i%20\in%20\mathcal{J}%20:=%20\{0,%20\ldots,%20J\}" /> is the associated group label (<img alt="formula" src="https://render.githubusercontent.com/render/math?math=x_i%20=%200" /> in the control group), <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\text{Ga}(a_0,%20b_0)" /> denotes the Gamma distribution with mean <img alt="formula" src="https://render.githubusercontent.com/render/math?math=a_0%20/%20b_0" />, <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\pi_0:\{0,%201\}^J%20\to%20(0,%201)" /> is the probability mass function described in section 2.1, and <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\text{DP}(\alpha,%20\bar{G})" /> denotes a Dirichlet Process (Ferguson, 1973) with concentration parameter <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\alpha" /> and base probability measure <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\bar{G}" />, in this case, <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\bar{G}%20:=%20\prod_{j%20\in%20\mathcal{J}}%20L_0" />, where <img alt="formula" src="https://render.githubusercontent.com/render/math?math=L_0" /> denotes the probability measure associated with a <img alt="formula" src="https://render.githubusercontent.com/render/math?math=\text{NIW}_D(\boldsymbol{u}_0,%20r_0,%20\nu_0,%20\boldsymbol{S}_0)" /> distribution.


