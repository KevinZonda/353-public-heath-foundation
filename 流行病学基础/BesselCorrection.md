# Bessel's Coorection

样本标准差（抽样调查）：
$$
s = \sqrt{\frac{\sum^n_i (x_i - \bar{x})^2}{n-1}}
$$

总体标准差（所有数据）：
$$
s = \sqrt{\frac{\sum^n_i (x_i - \mu)^2}{N}}
$$


## Proof

考虑数据 $X \sim D(\mu, \sigma^2)$。$|X| = n$, 我们期望：

$$
\mathbb{E}[s_X^2] = \mathbb{E}[\sigma^2]
$$

考虑不修正量方差为 $S^2_n$：

$$
S^2 = \frac{1}{n} \sum^n_{i=1} (x_i - \bar{x})^2
$$

$$
\begin{align*}
(x_i - \bar{x})^2 &= [(x_i - \mu)-(\bar{x} - \mu)]^2
\\
&= (x_i - \mu)^2 + (\bar{x}-\mu)^2 - 2(x_i-\mu)(\bar{x} - \mu)

\end{align*}
$$

$$
\begin{align*}
 \sum^n_{i=1} (x_i - \bar{x})^2
&=  \sum^n_{i=1} \left\{ (x_i - \mu)^2 + (\bar{x}-\mu)^2 - 2(x_i-\mu)(\bar{x} - \mu) \right\}
\\
&=  \sum^n_{i=1} (x_i - \mu)^2 +
\sum^n_{i=1} (\bar{x}-\mu)^2 - 2(\bar{x} - \mu)\sum^n_{i=1}(x_i-\mu)
\\
&=  \sum^n_{i=1} (x_i - \mu)^2 +
\sum^n_{i=1} (\bar{x}-\mu)^2 - 2(\bar{x} - \mu)\sum^n_{i=1}(\bar{x}-\mu)
\\
&=  \sum^n_{i=1} (x_i - \mu)^2 +
\sum^n_{i=1} (\bar{x}-\mu)^2 - 2n(\bar{x} - \mu)^2
\\
&=  \sum^n_{i=1} (x_i - \mu)^2 +
n(\bar{x}-\mu)^2 - 2n(\bar{x} - \mu)^2
\\
&=  \sum^n_{i=1} (x_i - \mu)^2
-n(\bar{x} - \mu)^2
\end{align*}
$$

$$
\begin{align*}
S^2 &= \frac{1}{n} \sum^n_{i=1} (x_i - \bar{x})^2
\\
&= \frac{1}{n} \left\{ 
    \sum^n_{i=1} (x_i - \mu)^2-n(\bar{x} - \mu)^2
\right\}
\\
&= 
\frac{1}{n}\sum^n_{i=1} (x_i - \mu)^2-(\bar{x} - \mu)^2
\end{align*}
$$

求期望：

$$
\begin{align*}
\mathbb{E}[S^2]
&=
\frac{1}{n}\sum^n_{i=1} \mathbb{E}[(x_i - \mu)^2]-\mathbb{E}[(\bar{x} - \mu)^2]

\\ &=
\frac{1}{n}\sum^n_{i=1} \mathbb{E}[\sigma^2]-\frac{\sigma^2}{n}

\\ &=
\sigma^2-\frac{\sigma^2}{n}
\\

&= \frac{n-1}{n}\sigma^2

\end{align*}
$$

因此如果不做修正，算出来的方差期望是 $\frac{n-1}{n}\sigma^2$，因此需要乘上 $\frac{n}{n-1}\sigma^2$ 做系数修正。