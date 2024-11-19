---
layout: page
title: method
---
\href{https://academic.oup.com/biomet/article/107/3/609/5820553?login=true}{Grace Yoon}

\href{https://www.nature.com/articles/s41467-023-40503-7}{cscore}

### Introduction

Inference of gene co-expression is pivotal for our understanding of biological functions and gene regulations. For scRNA-seq data, the expression level of a specific gene is measured through the observed UMI (unique molecular identifier) count for this gene, and the sequencing depth of a cell is the sum of UMI counts across all genes. The gene co-expression is measured by correlations of UMI counts.

This tutorial aims to describe the process of inferring gene co-expression with a rank-based estimator based on a latent Gaussian copula model and conduct a simulation study that assesses the performance of estimators. To construct co-expression networks, we first review Chang Su's approach called CS-CORE (cell-type-specific co-expressions)[^fn1]. CS-CORE models the underlying expression levels as latent variables, linked to the observed UMI counts through a Poisson measurement model depending on the underlying expression levels and sequencing depth. Under this model, CS-CORE implements an iteratively re-weighted least squares approach for estimating the true correlations between underlying expression levels. Building upon the CS-CORE framework, we incorporate the semiparametric approach proposed by Grace Yoon[^fn2], which models the latent variables with a Gaussian copula and use a rank-based estimator for correlation estimation.

The observed UMI counts for $p$ genes of cell $i, i=1,\cdots,n$ from $n$ cells is denoted by $(x_{i1},\cdots,x_{ip})$. ${s}_{i}=\mathop{\sum }\nolimits_{j=1}^{p}{x}_{ij}$ denotes the sequencing depth of cell $i$, which is the sum of UMI counts across $p$ genes. $(z_{i1},\cdots,z_{ip})$ denote the underlying expression levels from $p$ genes in cell $i$, defined to be the UMI counts from each gene relative to the sequencing depth $s_i$. Assume that:
$$(z_{i1},...,z_{ip}) \sim F_p({\bf \mu,\Sigma}),x_{ij}|z_{ij}\sim \text{Poisson}(s_iz_{ij})$$\\
where $F_p(\bf \mu,\Sigma)$ is an unknown nonnegative $p$-variate distribution with mean vector $\mu=(\mu_1,\cdots,\mu_p)$, $\sum_{j=1}^p \mu_j = 1$ and covariance matrix ${\bf\Sigma}  = (\sigma_{jj'})_{p \times p}$, the element in row $j$ and column $j'$ being $\sigma_{jj'}$. Here, $x_{ij}$ is the UMI count of gene $j$ in cell $i$, assumed to follow a Poisson measurement model depending on the underlying expression level $z_{ij}$ and sequencing depth $s_i$. This Poisson measurement model explicitly accounts for the sequencing depths so that the gene co-expressions measured by ${\bf\Sigma} = (\sigma_{jj'})_{p \times p}$, the covariance of $(z_{i1},...,z_{ip})$, is not biased by sequencing depths.

CS-CORE uses a moment-based iteratively reweighted least squares approach to estimate the covariance $\Sigma$. Given $\{x_{i1},...,x_{ip}\}_i^n$ and $\{s_i\}_i^n$, under the model in assumptions, $E({x}_{ij})={s}_{i}{\mu }_{j}, {{{{{{{\rm{Var}}}}}}}}({x}_{ij})=E[{({x}_{ij}-{s}_{i}{\mu }_{j})}^{2}]={s}_{i}{\mu }_{j}+{s}_{i}^{2}{\sigma }_{jj}, E[({x}_{ij}-{s}_{i}{\mu }_{j})({x}_{i{j}^{{\prime} }}-{s}_{i}{\mu }_{{j}^{{\prime} }})]={s}_{i}^{2}{\sigma }_{j{j}^{{\prime} }}$, so the following set of regression equations hold: 
$$
\begin{array}{rlr}
&{x}_{ij}={s}_{i}{\mu }_{j}+{\epsilon }_{ij}, &\\ 
&{({x}_{ij}-{s}_{i}{\mu }_{j})}^{2}={s}_{i}{\mu }_{j}+{s}_{i}^{2}{\sigma }_{jj}+{\eta }_{ij}, \\ 
&({x}_{ij}-{s}_{i}{\mu }_{j})({x}_{i{j}^{{\prime} }}-{s}_{i}{\mu }_{{j}^{{\prime} }})={s}_{i}^{2}{\sigma }_{j{j}^{{\prime} }}+{\xi }_{ij{j}^{{\prime} }},
\end{array}
$$
where ${\epsilon }_{ij}, {\eta }_{ij}, {\xi }_{ij{j}^{{\prime} }}$ are independent and mean-zero error variables for all $i, j, j'$.  
Then $\mu_j, \sigma_{jj}, \sigma_{jj'}$ are estimated by weighted least squares:
$$
\begin{aligned}
    \hat \mu_j &= {\min }_{\mu }\mathop{\sum }\nolimits_{i=1}^{n}{w}_{ij}{({x}_{ij}-{s}_{i}\mu )}^{2}\\
    \hat \sigma_{jj} &= {\min }_{\sigma }\mathop{\sum }\nolimits_{i=1}^{n}{h}_{ij}{[{({x}_{ij}-{s}_{i}{\hat{\mu }}_{j})}^{2}-{s}_{i}{\hat{\mu }}_{j}-{s}_{i}^{2}\sigma ]}^{2}\\
    \hat \sigma_{jj'} &= {\min }_{\sigma }\mathop{\sum }\nolimits_{i=1}^{n}{g}_{ij{j}^{{\prime} }}{[({x}_{ij}-{s}_{i}{\hat{\mu }}_{j})({x}_{i{j}^{{\prime} }}-{s}_{i}{\hat{\mu }}_{{j}^{{\prime} }})-{s}_{i}^{2}\sigma ]}^{2}\\
\end{aligned}
$$
where ${w}_{ij}, h_{ij}, g_{ijj'}$ are the weights. In the interative process, the weights are updated: ${w}_{ij}=1/{{{{{{{\rm{Var}}}}}}}}({\epsilon }_{ij})=1/({s}_{i}{\mu }_{j}+{s}_{i}^{2}{\sigma }_{jj})$, ${h}_{ij}={w}_{ij}^{2}$ and ${g}_{ij{j}^{{\prime} }}={w}_{ij}{w}_{i{j}^{{\prime} }}$. Finally we have $\hat\mu_j, {\bf\hat\Sigma}=(\hat\sigma_{jk})_{p \times p}$, leading to the estimation of the correlation matrix $\bf \hat R$, whose element in row $j$ and column $k$ is ${\hat{\rho }}_{jk}={\hat{\sigma }}_{jk}/\sqrt{{\hat{\sigma }}_{jj}{\hat{\sigma }}_{kk}}$. We called the estimator R.cscore.

The proposed model in CS-CORE does not impose distributional assumptions on the underlying expression levels, its distribution being an unknown nonnegative $p$-variate distribution, $(z_{i1},\cdots,z_{ip}) \sim F_p(\mu, {\bf \Sigma})$. We reviewed Grace Yoon's paper and adopted its approach, modeling latent variables with a Gaussian copula. The assumption is that the Gaussian copula can be transformed into the distribution of observed count data through unknown monotonic transformations, and due to the invariance of Kendall's $\tau$, the measure of correlation based on ranks, under monotonic transformations, a bridge function can link the Kendall's $\tau$ of the observed data to the correlation matrix of the latent Gaussian copula. We justify applying this assumption to observed UMI counts, that is, the existence of monotonic functions that transform the Gaussian distribution into the Poisson measurement model.

Assume that a truncated Gaussian copula underlies ($z_{i1},...,z_{ip}$), where $(z_{i1}, \cdots, z
_{ip})$ is the positive truncation of a continuous vector $(u_{i1},\cdots,u_{ip})$, that is, there exist $(u_{i1},\cdots,u_{ip})$, a set of monotonically increasing functions $f = (f_j)_{j=1}^p$ and a vector of constants ${\bf C} = (C_1,\cdots C_p)$, $f_j(u_{ij})$ denoted as $v_{ij}$, such that:
$$
\begin{aligned}
    &z_{ij} = I(u_{ij}>C_j)u_{ij}\\
    &(v_{i1},\cdots,v_{ip}) = (f_1(u_{i1}),\cdots,f_p(u_{ip}))\\
    &(v_{i1},\cdots,v_{ip}) \sim N_p(0,{\bf R}),R_{jj}=1\\
\end{aligned}
$$
Our goal is to estimate the latent correlation matrix ${\bf R}$, though transformations $f = (f_j)_{j=1}^p$ are unknown so latent variables $(v_{i1},\cdots,v_{ip})$ are not directly observable. ${\bf R}$ can be connected to the observed data by two steps: 1. Kendall's $\tau$ of the observed $(z_{i1},\cdots,z_{ip})$ and that of latent Gaussian copula $(v_{i1},\cdots,v_{ip})$ are the same, so it bypasses the unknown transformations $f = (f_j)_{j=1}^p$; 2. ${\bf R}$, the correlation matrix of $(v_{i1},\cdots,v_{ip})$, can be linked to its Kendall's $\tau$ by a bridge function, whose  explicit formulas for truncated data are derived in Grace Yoon's paper. 

First, we demonstrate the invariance of Kendall's $\tau$. Given the observed data of the gene pair $(z_{ij},z_{ik}), i=1,\cdots,n$, the sampled Kendall's $\tau$ is defined as $$\hat\tau_{jk} =\frac{2}{n(n-1)}\sum_{1\le i<i'\le n} sign(z_{ij}-z_{i'j})sign(z_{ik}-z_{i'k})$$, and the population Kendall's $\tau$ is defined as 
$\tau_{jk}=E(\hat\tau_{jk})$.  
Denote $\Delta_j = f_j(C_j), j=1,\cdots,p$, because for any monotonic $f_j()$, $$sign(z_{ij}-z_{i'j}) = sign(I(v_{ij}>\Delta_j)v_{ij}-I(v_{i'j}>\Delta_j)v_{i'j})$$, $\tau_{jk}$ of $(z_{ij},z_{ik})$ is equivalent to $\tau_{jk}$ of $(I(v_{ij}>\Delta_j)v_{ij},I(v_{ik}>\Delta_k)v_{ik})$, depending only on $v_{ij}, v_{ik}, \Delta_j, \Delta_k$. $\tau_{jk}$ is a function of ${\bf R}_{jk}$ independent of $(f_j)_{j=1}^p$, and we denote this function as bridge function $\tau_{jk} = F({\bf R}_{jk};\Delta_j,\Delta_k)$ with parameters $\Delta_j,\Delta_k$.

The bridge function establishes the connection between the latent correlation $\bf R$ and the population Kendall's $\tau$. The explicit form of the bridge function F is derived in \href{https://academic.oup.com/biomet/article/107/3/609/5820553?login=true}{Grace Yoon}. That is, $F({\bf R}_{jk};\Delta_j, \Delta_k) = -2 \Phi_4 (-\Delta_j, -\Delta_k, 0,0; \Sigma_{4a}) + 2 \Phi_4 (-\Delta_j, -\Delta_k, 0,0; \Sigma_{4b})$, with 
$$\begin{aligned}
\Sigma_{4a} & = \begin{pmatrix}
1 \; & \; 0 \; &\; 1/\sqrt{2} \; & \; -{\bf R}_{jk}/\sqrt{2}\\
0\; &\; 1 \; &\; -{\bf R}_{jk}/\sqrt{2}\; & \; 1/\sqrt{2}\\
1/\sqrt{2}\; &\; -{\bf R}_{jk}/\sqrt{2} \; &\; 1\; & \; -{\bf R}_{jk}\\
-{\bf R}_{jk}/\sqrt{2}\; & \; 1/\sqrt{2}\; &\; -{\bf R}_{jk}\; &\; 1
\end{pmatrix}\!,
\\
\Sigma_{4b} & = \begin{pmatrix}
1 \; & \; {\bf R}_{jk}\; &\; 1/\sqrt{2}\; & \;{\bf R}_{jk}/\sqrt{2}\\
{\bf R}_{jk}\; &\; 1\; & \;{\bf R}_{jk}/\sqrt{2}\; &\; 1/\sqrt{2}\\
1/\sqrt{2}\; & \;{\bf R}_{jk}/\sqrt{2}\; & \;1\; &\; {\bf R}_{jk}\\
{\bf R}_{jk}/\sqrt{2}\; & \;1/\sqrt{2}\; & \; {\bf R}_{jk}\; & \;1
\end{pmatrix}\!.
\end{aligned}$$
The derivation of the bridge function involves plugging the expression $sign(x) = 2I(x > 0) - 1$ into the definition of $\tau_{jk}$. This allows us to rewrite it as the expectation of indicator functions, and then the joint probability of multivariate normal distribution. Ultimately, this can be rewritten as cumulative normal distribution functions.

For example, without loss of generality we set j = 1 and k = 2, $(v_{i1},v_{i2}),(v_{i'1},v_{i'2})$ as $(v_1,v_2),(v'_1,v'_2)$. 
$$
\begin{aligned}
    \tau_{jk} = &-2E[I(v_1\leq\Delta_1,v'_2-v_2<0)]+2E[I(v'_1\leq\Delta_1,v'_2-v_2<0)]\\
    &+E[I(v_1>\Delta_1,v'_1>\Delta_1)sign(v_1-v'_1)sign(v_2-v'_2)]
\end{aligned}
$$
The last term can be written as joint probability
$$
\begin{aligned}
    &P(v_1>\Delta_1,v'_1>\Delta_1,v_1-v'_1>0,v_2-v'_2>0)\\
    +&P(v_1>\Delta_1,v'_1>\Delta_1,v_1-v'_1<0,v_2-v'_2<0)\\
    -&P(v_1>\Delta_1,v'_1>\Delta_1,v_1-v'_1>0,v_2-v'_2<0)\\
    -&P(v_1>\Delta_1,v'_1>\Delta_1,v_1-v'_1<0,v_2-v'_2>0)\\
\end{aligned}
$$
And then because $(v_1,v'_1,v_1-v'_1,v_2-v'_2)$ follows a multivariate distribution, we can compute the covariance matrix and rewrite it as CDF $\Phi_4(\cdot, \cdot, \cdot, \cdot;\Sigma)$.  
Because $\Delta_j = f_j(C_j)$ is unknown in practice, we used an estimator. Based on the moment equation ${E} \{I(z_{ij}>0) \} = {\rm pr} (z_{ij}>0) = {\rm pr} \{f_j(u_{ij})>\Delta_j\} = 1-\Phi(\Delta_j)$, we use $\hat\Delta_j=\Phi^{-1}(\sum_{i=1}^n I(z_{ij}=0)/n)$.

The estimation procedure of $R$ given the observed data matrix ${\rm X}\in R^{n \times p}$:  
Step 1: $s_i = \sum_{j=1}^p E(x_{ij})$, estimate $\hat s_i = \sum_{j=1}^p x_{ij}$  
Step 2: $\hat z_{ij} = x_{ij}/\hat s_i$  
Step 3: estimate $\hat \tau_{jk}$ from $\hat Z = (\hat z_{ij})_{n \times p}$  
Step 4: estimate $\hat \Delta_j = \Phi^{-1}(\sum_{i=1}^n I(\hat z_{ij}=0)/n)$  
Step 5: estimate R, $\hat R_{jk} = F^{-1}(\hat \tau_{jk})$ with parameters $\hat\Delta_j,\hat\Delta_k$  
We call the rank-based estimator R.rank.  

Simulation Study  
To validate our approach, we conduct a simulation study that assesses the performance of our estimators under various conditions. We followed the simulation method in Chang Su's paper. We specified the distribution of true expression level $z_{ij}$ to be Gamma distribution, marginally $z_{ij} \sim \text{ Gamma }(\alpha_j,\beta_j), \mu_j = \alpha_j/\beta_j, \sigma_{jj} = \alpha_j/\beta_j^2$, where $\mu_j$ and $\sigma_{jj}$ correspond to the marginal mean and variance in assumption $F_p(\mu, {\bf \Sigma})$. Conditional on $z_{ij}$, we simulated counts $x_{ij} \sim \text{Poisson}(s_iz_{ij})$. We obtain $s_i,\mu_j,\sigma_{jj}$ from the sequencing depths, mean and variance of real data, and then obtain $\alpha_j,\beta_j$ from $\mu_j = \alpha_j/\beta_j, \sigma_{jj} = \alpha_j/\beta_j^2$.  
Next, given a correlation matrix $\bf R^*$, we adopted a Gaussian copula to generate correlated Gamma random variables: $(v_{i1},...,v_{ip}) \sim N(0,\bf R^*)$ and then generate ${z}_{ij}={F}_{j}^{-1}({{\Phi }}({v}_{ij}))$ where $F_j$ is the CDF of $\text{ Gamma }(\alpha_j,\beta_j)$, $\bf R^*$ is estimated from real data by Cscore.  
$(z_{i1},...,z_{ip}) \sim F_p(\mu,\Sigma)$, where $\sigma_{jk} = cov({F}_{j}^{-1}({{\Phi }}({v}_{ij})),{F}_{k}^{-1}({{\Phi }}({v}_{ik})))$  
For details and implementation of simulation in R, see

### References

[^fn1]: Su, C., Xu, Z., Shan, X. et al. Cell-type-specific co-expression inference from single cell RNA-sequencing data. Nat Commun 14, 4846 (2023).

[^fn2]: Grace Yoon, Raymond J Carroll, Irina Gaynanova, Sparse semiparametric canonical correlation analysis for data of mixed types, Biometrika, Volume 107, Issue 3, September 2020, Pages 609–625.
