# cTWAS and INTERFACE

## cTWAS

cTWAS (causal-TWAS) includes the tested genes and all genetic confounders in a Bayesian variable selection model. 


$$
\begin{equation}
y = \sum_{j} \beta_j \tilde{X}_j + \sum_{m} \theta_m G_m + \epsilon
\end{equation}
$$

$\beta_j$: Effect size of gene $j$.
$\theta_m$: Effect size of variant $m$.
$\epsilon \sim \mathcal{N}(0, \sigma^2)$: Normally distributed error term.

Write the general model with K groups of explanatory variables as:
$$
\begin{equation}
y=\sum_{k=1}^K\sum_{j\in M_k}\beta_jX_j+\epsilon 
\end{equation}
$$
where $M_1, M_2$ denotes the group of genes and group of SNPs.

Both gene effects $\beta_j$ and variant effects $\theta_m$ are modeled using spike-and-slab priors to capture the sparsity of causal signals:

$$
\begin{equation}
\begin{split}
\gamma_j &\sim \text{Bernoulli}(\pi_k) \\
\beta_j | \gamma_j = 1 &\sim \mathcal{N}(0, \sigma_k^2) \\
\beta_j | \gamma_j = 0 &\sim \delta_0
\end{split}
\end{equation}
$$

$\gamma_j$: Indicator variable denoting whether gene $j$ is causal.
$\pi_k$: Prior probability of a gene or variant being causal within group $k$.
$\sigma_k^2$: Prior variance of the effect size for causal genes or variants.

The inference has two main steps. In the first step, we estimate the prior parameters $\boldsymbol\theta=\{\pi_k, \sigma_k^2\}$ for the two groups, gene effects and  variants effects. In the second step, we use the estimated $\boldsymbol\theta$, and compute the PIP of each variable.

To estimate the prior parameters (\( \pi_k \) and \( \sigma_k^2 \)) of the spike-and-slab priors by treating the causal indicators (\( \gamma_j \)) and effect sizes as missing data, we use the EM algorithm:

- **E-Step:** This involves computing the posterior inclusion probabilities (\( \alpha_j \)) and the posterior second moments of the effect sizes (\( \tau_j^2 \)).
     
- **M-Step:** Update the prior parameters by maximizing the expected log-likelihood obtained from the E-step. The updates for \( \pi_k \) and \( \sigma_k^2 \) are derived based on these expectations.

$$
\begin{align}
\pi_{k}^{(t+1)}=\frac{1}{|\mathbf{M}_{k}|}\sum_{j\in\mathbf{M}_{k}}\alpha_{j}^{(t)}\\
\sigma_{k}^{2,(t+1)}=\frac{\sum_{j\in\mathbf{M}_{k}}\alpha_{j}^{(t)}\cdot\tau_{j}^{2,(t)}}{\sum_{j\in\mathbf{M}_{k}}\alpha_{j}^{(t)}}
\end{align}
$$

> The new parameter $\pi_{k}^{(t+1)}$ is simply  the average PIP of all variables in the group k and the new $\sigma_{k}^{2,(t+1)}$ is the weighted average of the second moment of the posterior effect sizes.

To manage computational complexity and high correlations among variables, **the model assumes that each genomic block contains at most one causal gene or variant**.

cTWAS extends its model to work with GWAS summary statistics by incorporating marginal associations and the correlation structure among variables (linkage disequilibrium, LD).

We need R, the correlation matrix, of all variables for both gene and SNP. Note that if $\boldsymbol{X}_j$ is an imputed gene, we can write it as $\boldsymbol{X}_j=\boldsymbol{X}_j^s \boldsymbol{w}_j^T$. When $\boldsymbol{X}_i$ is variant, then:

$$
\begin{equation}
\begin{split}
R_{i j}&=\frac{\operatorname{Cov}\left(\boldsymbol{X}_{i}, \boldsymbol{X}_{j}\right)}{\sqrt{\operatorname{Var}\left(\boldsymbol{X}_{i}\right) \operatorname{Var}\left(\boldsymbol{X}_{j}\right)}} \\
&=\frac{\operatorname{Cov}\left(\boldsymbol{X}_{i}, \boldsymbol{X}_{j}^{s} \boldsymbol{w}_{j}^{T}\right)}{\sqrt{\operatorname{Var}\left(\boldsymbol{X}_{i}\right) \operatorname{Var}\left(\boldsymbol{X}_{j}^{s} \boldsymbol{w}_{j}^{T}\right)}}\\
&=\frac{\boldsymbol{R}_{i, j}^{s} \boldsymbol{w}_{j}^{T}}{\sqrt{\boldsymbol{w}_{j} \boldsymbol{R}_{j}^{s} \boldsymbol{w}_{j}^{T}}}
\end{split}
\end{equation}
$$

When both variables i and j are genes:
$$
\begin{equation}
\begin{split}
R_{i j}&=\frac{\operatorname{Cov}\left(\boldsymbol{X}_{i}, \boldsymbol{X}_{j}\right)}{\sqrt{\operatorname{Var}\left(\boldsymbol{X}_{i}\right) \operatorname{Var}\left(\boldsymbol{X}_{j}\right)}} \\
&=\frac{\operatorname{Cov}\left(\boldsymbol{X}_{i}^{s} \boldsymbol{w}_{i}^{T}, \boldsymbol{X}_{j}^{s} \boldsymbol{w}_{j}^{T}\right)}{\sqrt{\operatorname{Var}\left(\boldsymbol{X}_{i}^{s} \boldsymbol{w}_{i}^{T}\right) \operatorname{Var}\left(\boldsymbol{X}_{j}^{s} \boldsymbol{w}_{j}^{T}\right)}}\\
&=\frac{\boldsymbol{w}_{i}\boldsymbol{R}_{i, j}^{s} \boldsymbol{w}_{j}^{T}}{\sqrt{\boldsymbol{w}_{i} \boldsymbol{R}_{i}^{s} \boldsymbol{w}_{i}^{T}}\sqrt{\boldsymbol{w}_{j} \boldsymbol{R}_{j}^{s} \boldsymbol{w}_{j}^{T}}}
\end{split}
\end{equation}
$$

## INTERFACE

INTERFACE (INTEgRative Fine-mapping of Causal gEnes) is a probabilistic fine-mapping method designed to identify putative causal genes (PCGs). 

One of the main contributions is the development of interpretable informative priors for probabilistic causal gene fine mapping, integrating marginal TWAS and colocalization evidence. 

In summary, all the required priors in INTERFACE are data-driven and well-justified.

**Structural Equation Model (SEM):** 

INTERFACE employs a SEM that simultaneously models the relationships between multiple genetic variants, multiple molecular traits (e.g., gene expressions or protein levels), and a complex trait of interest. This SEM accounts for both gene-to-trait effects and direct variant-to-trait effects within a genomic region containing multiple gene candidates.

  \[
  \begin{aligned}
  M_i &= \mu_{M,i} + G_i \beta_{E,i} + e_{M,i}, \quad e_{M,i} \sim \mathcal{N}(0, \sigma^2_{M,i} I), \quad i = 1, \dots, q \\
  Y &= \mu_Y 1 + \sum_{i=1}^q \gamma_i M_i + G \beta_Y + e_Y, \quad e_Y \sim \mathcal{N}(0, \sigma^2_Y I)
  \end{aligned}
  \]

  Where:
  - \( M_i \) represents the molecular trait for gene \( i \).
  - \( G \) is the genotype matrix for genetic variants.
  - \( \beta_{E,i} \) captures the cis-eQTL effects for gene \( i \).
  - \( \gamma_i \) denotes the gene-to-trait effect sizes.
  - \( \beta_Y \) represents direct variant-to-trait effect sizes.
  - \( Y \) is the complex trait of interest.

INTERFACE utilizes a Bayesian framework to perform variable selection, determining which genes and genetic variants have non-zero effects on the trait. It assigns posterior inclusion probabilities (PIPs) to each candidate gene and variant, indicating the probability that they are causal.

INTERFACE leverages SuSiE algorithm for efficient Bayesian variable selection and fine-mapping. SuSiE accommodates multiple causal variants and genes within a region, handling complex linkage disequilibrium (LD) patterns effectively.

The three prior distributions critical for INTERFACE inference are:

- The prior probability that the region of interest has neither a PCG nor a direct-effect SNP: $$\Pr(\boldsymbol{\gamma}=\boldsymbol{0}\mathrm{~and~}\boldsymbol{\beta}_{Y}=\boldsymbol{0})=(1-\hat{p}_g)^p$$

- The prior probability that SNP $k$ has a direct effect on the complex trait: $$\Pr(\beta_{Y,k}\neq0)=\hat{p}_g-\hat{p}_c$$

- The prior probability that gene i is a PCG: $$\Pr(\gamma_i\neq0)=\hat{\pi}\cdot\text{GRCP}_i$$

where:
$\hat{p}_g$: the proportion of genome-wide variants that are GWAS hits.

$\hat{p}_m$: the proportion of genome-wide variants that are QTLs.

$\hat{p}_c$: the proportion of genome-wide variants that are colocalized.

$\hat{\pi}$: the proportion of genome-wide genes that are TWAS genes.

$\text{GRCP}_i$: the gene-level colocalization probability for each candidate gene $i$.

The three quantities $\hat{p}_g$, $\hat{p}_m$ and $\hat{p}_c$ are estimated from the enrichment analysis in colocalization analysis. It employs a multiple imputation  scheme that leverages probabilistic association evidence from the fine-mapping analysis of complex trait  GWAS and molecular QTL mapping.
