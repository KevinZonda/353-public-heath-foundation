# Lecture 3 — The Basic Distributions

*Medical Statistics, Chapter 5 of the textbook. Prepared from the lecture transcript of 1 August 2026 (instructor: JG).*

**Learning outcomes.** By the end of this lecture you should be able to:

- standardise a Normal variable and use standard Normal tables in both directions;
- explain *why* standardisation is necessary (areas under the curve as probabilities);
- define the standard error of the mean and relate it to sampling variation;
- explain how the $t$ distribution arises when $\sigma$ is unknown, and use $t$ tables;
- construct and interpret confidence intervals for a population mean ($\sigma$ known and unknown);
- state the defining properties of the Binomial and Poisson distributions and their Normal approximations.

---

## 1 The Normal distribution

**Definition 1.1.** A continuous random variable $X$ is *Normally distributed* with mean $\mu$ and variance $\sigma^2$, written

$$
X \sim N(\mu, \sigma^2).
$$

Here $\mu$ is a *location parameter* (it shifts the curve along the $x$-axis) and $\sigma$ a *shape parameter* (it determines the spread). Consequently the Normal distribution is not one curve but a *family* of curves, one for each pair $(\mu, \sigma)$.

**Remark 1.2 (why areas matter).** The area under the curve over an interval equals the *relative frequency*, and hence — in the limit — the *probability* that $X$ falls in that interval. Almost every application in this lecture reduces, sooner or later, to computing an area under a curve.

### 1.1 Standardisation

Since the family is infinite, computing areas curve by curve is impractical. We therefore *standardise*:

$$
Z = \frac{X - \mu}{\sigma} \sim N(0, 1).
$$

Think of a traveller converting dollars into renminbi before shopping: standardisation is a change of units — a single, unified scale on which all area questions can be answered from one table.

Areas are integrals; for any interval $[x_1, x_2]$,

$$
P(x_1 \le X \le x_2) = F(x_2) - F(x_1),
$$

where $F(x)$ is the distribution function, i.e. the area from $-\infty$ up to $x$. Two results worth memorising:

$$
P(\mu - \sigma \le X \le \mu + \sigma) = 68.27\%, \qquad P(\mu - 1.96\sigma \le X \le \mu + 1.96\sigma) = 95\%.
$$

**Using the tables.** The table gives $\Phi(z) = P(Z \le z)$ for $z \le 0$ only. For $z > 0$ we exploit symmetry about $0$:

$$
\Phi(z) = 1 - \Phi(-z).
$$

*Example.* $\Phi(-1.21) = 0.1131$, read directly from the table (row $-1.2$, column $0.01$).

### 1.2 Application I: estimating a proportion

**Example 1.3 (serum uric acid).** A random sample of $n = 180$ healthy adult women gives a sample mean $\bar{x}$ and sample standard deviation $s$. We wish to estimate the proportion of *all* healthy women in the region whose uric acid exceeds $409.34$.

Note the notation discipline, which will matter throughout the course:

| population (parameters) | sample (statistics) |
|---|---|
| mean $\mu$ | mean $\bar{x}$ |
| standard deviation $\sigma$ | standard deviation $s$ |

The population parameters are unknown, so $\mu$ is replaced by $\bar{x}$ and $\sigma$ by $s$. Standardising the cut-off:

$$
z = \frac{409.34 - \bar{x}}{s} = 1.96.
$$

The question "what proportion exceeds $409.34$?" has become "what area lies to the right of $1.96$ under the standard Normal curve?", which is $2.5\%$ by symmetry of the central $95\%$ region. A second example of the same kind (the proportion between $230$ and $360$) standardises to the interval $[-1.21, 1.09]$.

**Remark 1.4.** Observe the chain of reasoning: proportion $\to$ relative frequency $\to$ area $\to$ standardised area. Every applied problem on the Normal distribution follows it.

### 1.3 Application II: medical reference ranges

**Definition 1.5.** A *medical reference range* is the fluctuation range of an anatomical, physiological or biochemical measurement in *the majority* of normal people. It is a range, not a constant, because individuals vary.

Its construction proceeds in six steps:

1. **Select a sufficiently large sample of "normal" individuals.** "Normal" is a *relative* notion: a person may serve as normal for potassium provided their condition (e.g. abnormal blood glucose) does not affect potassium; typically $n \geq 100$, or $150$–$300$ when many factors are involved.
2. **Control measurement error** — use a unified, accurate method.
3. **Decide whether to subgroup.** Where the standard differs by sex or age, separate ranges are required (e.g. uric acid: roughly $200$–$400$ in men, $100$–$300$ in women; bone density in adults vs children).
4. **Decide one-sided or two-sided limits**, according to the clinical purpose. Leucocyte count and blood glucose are abnormal when *too high or too low* (two-sided); vital capacity harms only when too low (lower one-sided); urinary lead harms only when too high (upper one-sided).
5. **Choose a suitable percentile limit**, conventionally the central $95\%$ of normal people. The distributions of normal and diseased individuals overlap, so the placement of the limit is a trade-off: moving it right reduces the chance of labelling a healthy person as diseased; moving it left reduces the chance of labelling a diseased person as healthy. Which error matters more depends on the purpose of the study.
6. **Choose the method of computation** according to the distribution of the data (details omitted here).

**Remark 1.6.** Because the two distributions overlap, conclusions drawn from a reference range are never absolute: it is a boundary for *judgement*, not a law of nature.

---

## 2 Sampling variation and the standard error

**Definition 2.1.** A *parameter* is a characteristic computed from the whole population ($\mu$, $\sigma$); a *statistic* is computed from a sample ($\bar{x}$, $s$). *Sampling error* is the discrepancy between a statistic and its parameter, and between statistics from different samples. It is unavoidable — individuals vary — but it follows regularities, and those regularities are the subject of this section.

**Thought experiment.** From the population of heights of 18-year-old girls, $X \sim N(\mu, \sigma^2)$, draw $100$ samples of size $n = 20$. This yields $100$ sample means $\bar{x}_1, \dots, \bar{x}_{100}$. Plotting them reveals three facts:

1. the means differ from one another and from $\mu$ — sampling error made visible;
2. the $100$ means are themselves (approximately) Normally distributed;
3. the mean of the $100$ means equals $\mu$.

The natural question — *how dispersed are these means?* — has a precise answer. The standard deviation of the sample mean is

$$
\sigma_{\bar{x}} = \frac{\sigma}{\sqrt{n}}.
$$

**Definition 2.2.** The standard deviation of a sample statistic is called its *standard error* (SE). Here the statistic is the sample mean, so $\sigma_{\bar{x}}$ is the *standard error of the mean*. A standard error is still a standard deviation in substance; the new name merely records the new object. It measures, in one number, the dispersion of the means, their typical distance from $\mu$, and hence the size of the sampling error.

Since $\sigma$ is an unknown parameter, we replace it by $s$ and obtain the *estimated* standard error:

$$
s_{\bar{x}} = \frac{s}{\sqrt{n}}.
$$

**Remark 2.3.** As $n$ increases, $s_{\bar{x}}$ decreases: larger samples carry information closer to the population, so the sampling error shrinks. This is visible directly from the formula.

---

## 3 The $t$ distribution

### 3.1 Why a new distribution is needed

Were $\sigma$ known, we could standardise the sample mean exactly:

$$
\frac{\bar{X} - \mu}{\sigma / \sqrt{n}} \sim N(0, 1).
$$

But $\sigma$ is unknown, and substituting $s$ gives an *estimate* in the denominator. The standard Normal curve is exacting — symmetry about $0$, variance exactly $1$ — and a statistic built on an estimated denominator no longer satisfies it exactly. We are therefore forced to introduce a new variable,

$$
T = \frac{\bar{X} - \mu}{S / \sqrt{n}} \sim t_{\,n-1},
$$

the $t$ distribution with $n - 1$ degrees of freedom, due to W. S. Gosset (writing under the pen-name "Student"). It is the theoretical basis of interval estimation and hypothesis testing for a population mean, and is especially important for small samples.

**Definition 3.1 (degrees of freedom).** The number of observations that can vary freely given the constraints. If three numbers have mean $5$ (hence sum $15$), two may be chosen freely but the third is then fixed; in general $\nu = n - (\text{number of constraints})$.

### 3.2 Properties of the $t$ curve

1. Single-peaked and symmetric about $0$.
2. Smaller $\nu$ $\Rightarrow$ greater dispersion: a lower peak and fatter tails.
3. As $\nu \to \infty$, $t_\nu \to N(0,1)$: the standard Normal is a limiting special case. (Indeed, as $n \to \infty$, $s \to \sigma$, so $S/\sqrt{n} \to \sigma/\sqrt{n}$ and the $t$ statistic becomes the $z$ statistic.)
4. For fixed $\nu$: the larger $|t|$, the smaller the tail area.
5. For fixed tail area: the larger $\nu$, the smaller the tabulated $t$.

**Using the $t$ table.** Given $\nu$ and a one- or two-sided tail area $\alpha$, the table yields the unique critical value; only $|t|$ is tabulated, negative values being recovered by symmetry.

*Examples.* $\nu = 7$: one-sided $\alpha = 0.05$ gives $t = 1.895$; two-sided $\alpha = 0.05$ (i.e. $0.025$ in each tail) gives $t = 2.365$. As a check of property 3: for $\nu = \infty$ and two-sided $\alpha = 0.05$ the table gives $1.96$, exactly the standard Normal value.

### 3.3 From areas to probabilities

The step that unlocks everything after this lecture is a change of interpretation. The statement "the two-sided tail area is $\alpha$" may be re-read as

$$
P\!\left(-t_{\alpha/2,\,\nu} \le T \le t_{\alpha/2,\,\nu}\right) = 1 - \alpha.
$$

The tail area is no longer geometry; it is the probability that $T$ falls outside, or inside, a stated interval. *E.g.* with $\nu = 38$, $t_{0.05,\,38} = 1.686$ means $P(T > 1.686) = 0.05$.

---

## 4 Interval estimation of a population mean

**Definition 4.1.** *Point estimation* uses the sample statistic directly as the parameter's value; it ignores sampling error and is rarely adequate. *Interval estimation* combines the statistic with its standard error to give an interval which contains the parameter with a stated high probability $1 - \alpha$, the *confidence level* (commonly $0.95$ or $0.99$). The interval is the *confidence interval* (CI), written e.g. a 95% CI; its endpoints are the lower and upper *confidence limits*, and $\alpha$ is the probability that the procedure goes wrong.

**Construction.** Since $T = (\bar{X} - \mu)/(S/\sqrt{n})$ contains $\mu$, the probability statement of §3.3 can be rearranged from a statement about $T$ into a statement about $\mu$. The result, for $\sigma$ **unknown**, is the $t$ method:

$$
\bar{x} \pm t_{\alpha/2,\,n-1} \cdot s_{\bar{x}} = \bar{x} \pm t_{\alpha/2,\,n-1} \cdot \frac{s}{\sqrt{n}}.
$$

When $\sigma$ is **known**, no substitution is needed, the statistic is exactly standard Normal, and we use the Normal method:

$$
\bar{x} \pm z_{\alpha/2} \cdot \sigma_{\bar{x}} = \bar{x} \pm z_{\alpha/2} \cdot \frac{\sigma}{\sqrt{n}}.
$$

**Example 4.2 (birth weights).** For newborn birth weights in a region, $\sigma = 360$ is given (a *population* value), and a hospital sample of $n = 540$ provides $\bar{x}$ and $s$. Since $\sigma$ is known, the Normal method applies, with $\sigma_{\bar{x}} = 360/\sqrt{540}$; the sample $s$ plays no role. Keep population and sample quantities strictly separate — misreading which is given is the classic error here.

**Interpretation.** A 95% CI is a probability statement about the *procedure*: in $100$ repeated samples, each producing its own interval, about $95$ intervals contain $\mu$ and $5$ do not.

**Accuracy versus precision.** The two quality criteria pull against each other:

- *Accuracy* (confidence level $1-\alpha$): raising the level widens the interval and increases our assurance of capturing $\mu$ — a 99% CI is "safer" than a 95% CI. Taken to the limit, a 100% CI is $(-\infty, +\infty)$: perfectly safe and perfectly useless. We work at $95\%$ or $99\%$ precisely because a useful answer must say *where* $\mu$ plausibly lies.
- *Precision* (width of the interval): the narrower the interval, the sharper the estimate.

The two cannot be improved simultaneously through $\alpha$ alone. The one exception: at a *fixed* confidence level, increasing the sample size $n$ narrows the interval (since $s_{\bar{x}} = s/\sqrt{n}$), improving precision at no cost to accuracy.

---

## 5 The Binomial distribution

**Definition 5.1.** A *Bernoulli trial* has three properties: (i) mutually exclusive outcomes (success/failure); (ii) a constant probability $\pi$ of success at each trial; (iii) independence of trials. The number of successes $X$ in $n$ independent Bernoulli trials follows the *Binomial distribution*, $X \sim B(n, \pi)$.

**Example 5.2.** Mice receiving a dose of poison die with probability $\pi = 0.8$; survival and death are mutually exclusive, each mouse faces the same $\pi$, and mice do not influence one another. For $n = 3$ mice, $X \in \{0, 1, 2, 3\}$ counts the deaths, and each value carries a probability given by a term of the binomial expansion:

$$
P(X = k) = \binom{n}{k} \pi^{k} (1 - \pi)^{n - k}, \qquad k = 0, 1, \dots, n.
$$

**Properties.**

- Mean $\mu = n\pi$ (the average count of successes: $3 \times 0.8 = 2.4$ mice); variance $\sigma^2 = n\pi(1 - \pi)$.
- The variance depends on both $n$ and $\pi$: larger $n$ widens the range of possible outcomes; $\pi$ nearer $0.5$ makes outcomes least predictable and maximises the variance, since $n\pi(1-\pi)$ is largest at $\pi = 0.5$.
- Shape: symmetric when $\pi = 0.5$; increasingly skewed as $\pi$ moves away from $0.5$; for fixed $\pi$, increasingly symmetric as $n$ grows.
- *Normal approximation*: for $n$ large and $\pi$ not close to $0$ or $1$, $B(n, \pi)$ is approximately Normal.

## 6 The Poisson distribution

**Definition 6.1.** The *Poisson distribution*, $X \sim P(\lambda)$, models rare events over many opportunities — a small $\pi$ with a large $n$ — such as genetic defects, cancer incidence, or counts of *E. coli* in $1$ ml of water. Its single parameter $\lambda = n\pi$ is the average number of occurrences.

**Properties.**

1. *Mean equals variance*: $E(X) = \mathrm{Var}(X) = \lambda$. This is its signature.
2. *Relation to the Binomial*: as $n \to \infty$ with $\pi$ small, $n\pi(1-\pi) \to n\pi$, so the Binomial variance approaches the Poisson variance. When their means are equal ($\mu = n\pi$), the Binomial variance $\mu(1-\pi)$ is strictly less than the Poisson variance $\mu$.
3. *Normal approximation*: for $\lambda \geq 20$ the Poisson may be approximated by a Normal distribution — an approximation only, since the Poisson remains discrete.
4. *Additivity*: sums of independent Poisson variables are Poisson (five $1$ ml water samples pooled into $5$ ml still follow a Poisson law); this is what makes Normal-approximation arguments convenient.

---

## Summary

One idea runs through the whole lecture: **areas under a curve are probabilities, and every device in this chapter exists to compute them.**

- Standardisation unifies the Normal family into a single tabulated curve.
- The standard error quantifies sampling variation; the unknown $\sigma$, replaced by $s$, forces us off the standard Normal and onto the $t$ distribution.
- Reading tail areas as probabilities turns the $t$ statement into a confidence interval for $\mu$.
- The Binomial and Poisson distributions handle the discrete case, with Normal approximations when $n$ or $\lambda$ is large.

**Looking ahead.** Reinterpreting a tail area of $0.05$ as "the probability of a $T$ value this extreme" is precisely the $P$-value; hypothesis testing (the $t$ test) begins next lecture.

---

*Transcription notes for readers of the Chinese transcript: "标准物" = standard error (标准误); "被努力实验" = Bernoulli trial (伯努利实验); "派" = $\pi$; "1%减阿尔法" = $1 - \alpha$; the "English statistician" is Gosset ("Student").*
