# Lecture 3 — 基本分布（The Basic Distributions）

*医学统计学，教材第五章。依据 2026 年 8 月 1 日课堂逐字稿整理（主讲：JG）。*

**Learning outcomes.** 学完本讲，你应当能够：

- 定义 probability density function 与 distribution function，并由二者推出区间概率的计算公式；
- 陈述并证明 standardisation 下概率不变这一事实，并以此完成 Normal probability 的完整计算；
- 定义 standard error，陈述 sample mean 的 sampling distribution；
- 解释当 $\sigma$ 未知时 $t$ distribution 如何产生，并由 $t$ 的概率陈述**推导**出 confidence interval 的公式（而非记忆它）；
- 陈述 Binomial 与 Poisson distribution 的定义、性质及其 Normal approximation。

---

## 0 预备：Random Variable、Density 与 Distribution Function

本讲的一切结论都建立在下面三个定义之上。后续每一节只引用它们，不再另设假设。

**Definition 0.1（random variable）.** *Random variable* 是随机试验结果的数值化表示。按其取值方式分为两类：

- *discrete*：取值可枚举（如死亡只数 $0, 1, 2, 3$），其概率由 probability mass function 逐点给出；
- *continuous*：取值充满某个区间（如身高、尿酸值），其概率由下面的 density 以面积形式给出。

本讲 §1–§6 处理 continuous 情形（Normal 与 $t$），§7–§8 处理 discrete 情形（Binomial 与 Poisson）。

**Definition 0.2（probability density function）.** 设 $X$ 为 continuous random variable。若存在函数 $f(x)$ 满足

$$
f(x) \ge 0, \qquad \int_{-\infty}^{+\infty} f(x)\,dx = 1,
$$

则称 $f(x)$ 为 $X$ 的 *probability density function*（pdf）。$X$ 落入任意区间的概率**定义**为 density 曲线在该区间上的面积，即积分：

$$
P(a \le X \le b) = \int_{a}^{b} f(x)\,dx.
$$

这就是"面积即概率"的准确含义：它不是比喻，而是 continuous 情形下概率的**定义**。全曲线下的总面积为 $1$，对应"概率总和为 $1$"。

**Definition 0.3（distribution function）.** $X$ 的 *distribution function*（cumulative distribution function, CDF）定义为

$$
F(x) = P(X \le x) = \int_{-\infty}^{x} f(t)\,dt,
$$

即从 $-\infty$ 到 $x$ 的曲线下面积。

**Proposition 0.4（区间概率公式）.** 对任意 $x_1 < x_2$，

$$
P(x_1 \le X \le x_2) = F(x_2) - F(x_1).
$$

*Proof.* 由积分的区间可加性，

$$
F(x_2) = \int_{-\infty}^{x_2} f(t)\,dt = \int_{-\infty}^{x_1} f(t)\,dt + \int_{x_1}^{x_2} f(t)\,dt = F(x_1) + P(x_1 \le X \le x_2),
$$

移项即得。$\square$

**Remark 0.5.** Proposition 0.4 是本章一切查表操作的合法性来源：表格存储的是 $F$ 的值，任何区间概率都可化为两个 $F$ 值之差。

---

## 1 Normal Distribution

**Definition 1.1.** 若 continuous random variable $X$ 的 pdf 为

$$
f(x) = \frac{1}{\sqrt{2\pi}\,\sigma}\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right), \qquad -\infty < x < +\infty,
$$

则称 $X$ 服从 mean 为 $\mu$、variance 为 $\sigma^2$ 的 Normal distribution，记作

$$
X \sim N(\mu, \sigma^2).
$$

其中 $\mu$ 是 *location parameter*（决定曲线沿 $x$ 轴的平移），$\sigma$ 是 *shape parameter*（决定曲线的离散形态）。因此 Normal distribution 是一*族*（*family*）曲线——每一对 $(\mu, \sigma)$ 对应一条。

### 1.1 Standardisation 及其概率不变性

由于曲线有无穷多条而表格只能有一张，我们需要把所有 Normal 问题化归到同一条曲线上。下面这个定理是本章最重要的工具。

**Theorem 1.2（standardisation）.** 若 $X \sim N(\mu, \sigma^2)$，则

$$
Z = \frac{X - \mu}{\sigma} \sim N(0, 1).
$$

$N(0,1)$ 称为 *standard Normal distribution*；其 CDF 有专用记号 $\Phi(z) = P(Z \le z)$，且已制成表格。

**Proposition 1.3（standardisation 下概率不变）.** 对任意实数 $a$，

$$
P(X \le a) = \Phi\!\left(\frac{a - \mu}{\sigma}\right).
$$

*Proof.* 对不等式两端做同一变换，概率不变：

$$
P(X \le a) = P(X - \mu \le a - \mu) = P\!\left(\frac{X - \mu}{\sigma} \le \frac{a - \mu}{\sigma}\right) = P\!\left(Z \le \frac{a - \mu}{\sigma}\right) = \Phi\!\left(\frac{a-\mu}{\sigma}\right),
$$

其中第三步用了 $\sigma > 0$（除以正数不改变不等号方向），第四步用了 Theorem 1.2。$\square$

**Remark 1.4.** Proposition 1.3 就是把任意 Normal 问题"翻译"成 standard Normal 问题的全部内容：左边是原始变量 $X$ 的陈述，右边是可查表的 $\Phi$ 的陈述。类比旅客把美元兑换成人民币再购物——standardisation 是统一度量衡，而本命题保证兑换前后"购买力"（概率）不变。

**Corollary 1.5（两个常用面积）.** 由 Proposition 1.3 与 Proposition 0.4，

$$
P(\mu - \sigma \le X \le \mu + \sigma) = \Phi(1) - \Phi(-1) = 68.27\%,
$$

$$
P(\mu - 1.96\sigma \le X \le \mu + 1.96\sigma) = \Phi(1.96) - \Phi(-1.96) = 95\%.
$$

### 1.2 Standard Normal 表：它是什么，为什么存在

**Remark 1.6（表格的本质）.** 由 Definition 0.3，$\Phi(z)$ 是一个积分：

$$
\Phi(z) = \int_{-\infty}^{z} \frac{1}{\sqrt{2\pi}}\, e^{-t^2/2}\, dt.
$$

这个积分**没有初等表达式**——$e^{-t^2/2}$ 的原函数不能用有限个初等函数写出，因而 $\Phi$ 的值无法手算。统计学家的解决办法是：预先用数值方法把 $\Phi(z)$ 在加密网格上的值算好，印成**表**（教材附录中的 standard Normal 表）。因此，表格不是新的数学对象，它就是 Definition 0.3 中那个 $F$（此处为 $\Phi$）的**数值化身**：Proposition 0.4 告诉我们"区间概率 = 两个 $F$ 值之差"，而表格提供 $F$ 的具体数值，二者合起来才构成完整的计算方法。今天计算机（如 Excel 的 `NORM.S.DIST`、Python 的 `scipy.stats.norm.cdf`）给出的是同一组数值，查表只是其前计算机时代的手工等价物。

**表的结构与查法.** 表的查取方向是"**给 $z$，查面积**"：

- 行：$z$ 到第一位小数（如 $-1.2$）；
- 列：第二位小数（如 $0.01$）；
- 行列交叉的单元格：$\Phi(z)$ 的值。

*Example.* 求 $\Phi(-1.21)$：在表中行 $-1.2$、列 $0.01$ 处直接读得 $\Phi(-1.21) = 0.1131$。

**为什么只印 $z \le 0$ 的一半.** 由 $N(0,1)$ 关于 $0$ 的 symmetry，$z$ 右侧的面积等于 $-z$ 左侧的面积，故

$$
\Phi(z) = 1 - \Phi(-z),
$$

$z > 0$ 的值可由 $z < 0$ 的值换算，印一半即可。

### 1.3 应用一：估计比例——一次完整的计算

**Example 1.7（serum uric acid）.** 随机抽取 $n = 180$ 名健康成年女性，得 sample mean $\bar{x}$ 与 sample standard deviation $s$。设该地区健康女性的尿酸值 $X$ 服从 Normal distribution，欲估计其中尿酸值超过 $409.34$ 者所占的比例。

先明确符号纪律——population 与 sample 各有一套记号，不可混用：

| population（parameters） | sample（statistics） |
|---|---|
| mean $\mu$ | mean $\bar{x}$ |
| standard deviation $\sigma$ | standard deviation $s$ |

计算分四步，每一步都有前文依据：

**第一步：比例即概率。** "全体女性中超过 $409.34$ 者所占的比例"即 relative frequency，在总体意义下就是 probability：

$$
\text{所求比例} = P(X > 409.34).
$$

**第二步：以 sample 估计 parameters。** $\mu$、$\sigma$ 未知，以 $\bar{x}$、$s$ 代之（这是估计，其误差由 §2 的 sampling variation 理论处理；此处先接受之）。

**第三步：standardise。** 由 Proposition 1.3，

$$
P(X > 409.34) = 1 - P(X \le 409.34) = 1 - \Phi\!\left(\frac{409.34 - \bar{x}}{s}\right).
$$

代入数据，界值 $409.34$ 经 $z$ 变换后对应

$$
z = \frac{409.34 - \bar{x}}{s} = 1.96.
$$

**第四步：查表。** 问题已化为 standard Normal 的右尾面积：

$$
P(X > 409.34) = 1 - \Phi(1.96) = 1 - 0.975 = 2.5\%.
$$

**Remark 1.8.** 请回看这条推理链：proportion $\to$ probability $\to$ standardised probability $\to$ 查表。第三步中界值 $409.34$ 必须和 $X$ *一起做 $z$ 变换*——只有变换后的 $z = 1.96$ 才能与 standard Normal 表对接；缺少这一步，查表就没有依据。第二个同类例子（介于 $230$ 与 $360$ 之间的比例）同法处理，standardise 后对应区间 $[-1.21, 1.09]$，由 Proposition 0.4 化为 $\Phi(1.09) - \Phi(-1.21)$ 求出。

### 1.4 应用二：Medical Reference Range

**Definition 1.9.** *Medical reference range*（医学参考值范围）指*大多数*正常人的解剖、生理、生化指标测量值的波动范围。它是一个范围而不是一个常数，因为人与人之间存在变异。

其制定分六步：

1. **选取足够例数的"正常人"。** 这里的"正常"是*相对*概念：只要所患疾病不影响待测指标（如血糖异常者不影响血钾测定），即可作为该指标的"正常人"；一般 $n \geq 100$，研究因素多者需 $150$–$300$。
2. **控制 measurement error** —— 采用统一、准确的测量方法。
3. **决定是否分组。** 若标准因性别、年龄而异，须分别制定（如 uric acid：男性约 $200$–$400$，女性约 $100$–$300$；成人与儿童的骨密度）。
4. **决定 one-sided 还是 two-sided 界限**，依研究目的而定。白细胞数和血糖过高或过低均属异常（two-sided）；肺活量仅在过低时有害（lower one-sided）；尿铅仅在过高时有害（upper one-sided）。
5. **选择合适的 percentile limit**，通常取正常人的中间 $95\%$。正常人与病人的分布曲线本有重叠，故界限的位置是一种 trade-off：右移减少"把正常人判为病人"的可能，左移减少"把病人判为正常人"的可能。究竟哪种错误更要紧，取决于研究目的。
6. **按资料的分布特征选择计算方法**（此处从略）。

**Remark 1.10.** 由于两类分布存在重叠，由 reference range 得出的结论从来不是绝对的：它是供*判断*的界限，而非自然法则。

---

## 2 Sampling Variation 与 Standard Error

**Definition 2.1.** *Parameter* 是由整个 population 算出的特征值（$\mu$、$\sigma$）；*statistic* 是由 sample 算出的特征值（$\bar{x}$、$s$）。*Sampling error* 指 statistic 与 parameter 之间、以及不同 sample 的 statistics 之间的差异。它不可避免——个体必然有变异——但它有规律可循。

### 2.1 Sample mean 的 Sampling Distribution

**Thought experiment.** 设 18 周岁女生身高总体 $X \sim N(\mu, \sigma^2)$，从中抽取 $100$ 个 sample，每个 $n = 20$，得到 $100$ 个 sample means $\bar{x}_1, \dots, \bar{x}_{100}$。把它们画出来，可见三个事实：

1. 各 means 彼此不同，也与 $\mu$ 不同——sampling error 由此可见；
2. 这 $100$ 个 means 本身（近似）服从 Normal distribution；
3. 这 $100$ 个 means 的 mean 恰等于 $\mu$。

这些观察可以严格化为：

**Theorem 2.2（sample mean 的 sampling distribution）.** 设 $X \sim N(\mu, \sigma^2)$，从中抽取容量为 $n$ 的 random sample，则

$$
\bar{X} \sim N\!\left(\mu, \frac{\sigma^2}{n}\right).
$$

即：sample mean 仍服从 Normal distribution，其 mean 与原总体相同，variance 缩小为原来的 $1/n$。

*Remark（证明从略，要点如下）.* 均值不变由期望的线性性质得到；方差缩小 $n$ 倍由独立随机变量之和的方差性质得到；正态性由"独立 Normal 变量之和仍是 Normal"得到。

**Remark 2.3（central limit theorem）.** 即使总体**不**服从 Normal distribution，只要 $n$ 足够大，$\bar{X}$ 也近似服从 $N(\mu, \sigma^2/n)$。这就是 *central limit theorem*，它是 Theorem 2.2 在非 Normal 总体下的替身，也是大样本统计推断的基石。

### 2.2 Standard Error

由 Theorem 2.2，$\bar{X}$ 的 standard deviation 为

$$
\sigma_{\bar{x}} = \frac{\sigma}{\sqrt{n}}.
$$

**Definition 2.4.** Sample statistic 的 standard deviation 称为 *standard error*（SE）。此处 statistic 是 sample mean，故 $\sigma_{\bar{x}}$ 即 mean 的 standard error。Standard error 在实质上仍是 standard deviation，新名字只是标明新对象：它刻画的是 *statistic*（而非个体观测）的离散程度，从而刻画 sampling error 的大小。

由于 $\sigma$ 是未知 parameter，以 $s$ 代之，得 standard error 的*估计值*：

$$
s_{\bar{x}} = \frac{s}{\sqrt{n}}.
$$

**Remark 2.5.** 随 $n$ 增大，$s_{\bar{x}}$ 减小：sample 越大，sample 信息越逼近 population，sampling error 越小。这一点从公式中可以直接看出，也是 §4 中"增大 $n$ 可提高 precision"的依据。

---

## 3 $t$ Distribution

### 3.1 为什么需要一个新分布

由 Theorem 2.2 与 Theorem 1.2，若 $\sigma$ **已知**，sample mean 可以精确 standardise：

$$
\frac{\bar{X} - \mu}{\sigma / \sqrt{n}} \sim N(0, 1).
$$

但 $\sigma$ 通常未知。以 $s$ 替代后，分母从常数变成*随机变量*（$s$ 随 sample 而变），整个分式的分布随之改变——它不再精确满足 standard Normal 那条"严苛"的曲线（关于 $0$ 对称、variance 恰为 $1$）。新分布由 W. S. Gosset 以笔名 "Student" 给出：

**Definition 3.1（$t$ distribution）.** 设 $X \sim N(\mu, \sigma^2)$，sample 容量为 $n$，则

$$
T = \frac{\bar{X} - \mu}{S / \sqrt{n}} \sim t_{\,n-1},
$$

即 $T$ 服从 degrees of freedom 为 $n - 1$ 的 $t$ distribution。它是 population mean 的 interval estimation 与 hypothesis testing 的理论基础，在小样本情形下尤为重要。

**Definition 3.2（degrees of freedom）.** 在给定约束下能够自由取值的观测个数。若三个数的 mean 为 $5$（故 sum 为 $15$），前两个可自由取值，第三个随之被确定；一般地 $\nu = n - (\text{约束条件的个数})$。

### 3.2 $t$ 曲线的性质

1. 单峰，关于 $0$ 对称。
2. $\nu$ 越小越离散：峰越矮、尾部越翘。
3. 当 $\nu \to \infty$ 时 $t_\nu \to N(0,1)$：standard Normal 是 $t$ distribution 的极限特例。（由 §2.2，$n \to \infty$ 时 $s \to \sigma$，故 $S/\sqrt{n} \to \sigma/\sqrt{n}$，$T$ 统计量随之化为 §3.1 开头的 $z$ 统计量。）
4. 固定 $\nu$：$|t|$ 越大，尾部面积越小。
5. 固定尾部面积：$\nu$ 越大，表中查得的 $t$ 越小。

**使用 $t$ table（注意方向与 Remark 1.6 相反）.** Standard Normal 表的查询方向是"给 $z$，查面积"（$z \to \Phi(z)$）；$t$ 表（教材附录中的 $t$ 界值表）方向相反，是"**给面积，查界值**"：

- 行：degrees of freedom $\nu$（固定曲线形状）；
- 列：尾部面积 $\alpha$（区分 one-sided 与 two-sided）；
- 行列交叉的单元格：临界值 $t_{\alpha,\,\nu}$，即"使 $P(T > t_{\alpha,\,\nu}) = \alpha$ 成立的那个界值"。

之所以反向制表，是因为 $t$ 表的用途不是算概率，而是为 interval estimation 与 hypothesis testing 提供临界值（§4 的推导直接用到它）。表中只列 $|t|$，负值由 symmetry 补出。

*Examples.* $\nu = 7$：one-sided $\alpha = 0.05$ 对应 $t = 1.895$；two-sided $\alpha = 0.05$（即每尾 $0.025$）对应 $t = 2.365$。作为性质 3 的验证：$\nu = \infty$、two-sided $\alpha = 0.05$ 时表中给出 $1.96$，恰为 standard Normal 之值。

### 3.3 从面积到概率

由 Definition 0.2，$t$ 曲线下的面积即概率。于是"two-sided 尾部面积为 $\alpha$"可以改写为一个关于 $T$ 的概率陈述：

$$
P\!\left(-t_{\alpha/2,\,\nu} \le T \le t_{\alpha/2,\,\nu}\right) = 1 - \alpha.
$$

*E.g.* 当 $\nu = 38$ 时，$t_{0.05,\,38} = 1.686$ 意即 $P(T > 1.686) = 0.05$。**这个概率陈述是 §4 全部推导的起点。**

---

## 4 Population Mean 的 Interval Estimation

### 4.1 两种估计

**Definition 4.1.** *Point estimation* 直接用 sample statistic 充当 parameter 之值；它不考虑 sampling error，故很少够用。*Interval estimation* 则给出一个以给定的高概率 $1 - \alpha$ 包含 parameter 的区间。$1 - \alpha$ 称为 *confidence level*（常取 $0.95$ 或 $0.99$），区间称为 *confidence interval*（CI），端点称为 lower/upper *confidence limits*，$\alpha$ 是估计程序犯错的概率。

### 4.2 Confidence Interval 的推导（$\sigma$ 未知：$t$ method）

此处**不给出公式，而是推出公式**。起点是 §3.3 的概率陈述（取 $\nu = n - 1$）：

$$
P\!\left(-t_{\alpha/2,\,n-1} \le T \le t_{\alpha/2,\,n-1}\right) = 1 - \alpha.
$$

**第一步：代入 $T$ 的定义。** 由 Definition 3.1，$T = (\bar{X} - \mu)/(S/\sqrt{n})$，代入得

$$
P\!\left(-t_{\alpha/2,\,n-1} \le \frac{\bar{X} - \mu}{S/\sqrt{n}} \le t_{\alpha/2,\,n-1}\right) = 1 - \alpha.
$$

**第二步：解关于 $\mu$ 的不等式。** 三边同乘 $S/\sqrt{n}$（正数，不等号方向不变），再移项，得

$$
P\!\left(\bar{X} - t_{\alpha/2,\,n-1}\cdot\frac{S}{\sqrt{n}} \le \mu \le \bar{X} + t_{\alpha/2,\,n-1}\cdot\frac{S}{\sqrt{n}}\right) = 1 - \alpha.
$$

**第三步：读出 confidence interval。** 上式表明：随机区间

$$
\bar{x} \pm t_{\alpha/2,\,n-1}\cdot s_{\bar{x}} = \bar{x} \pm t_{\alpha/2,\,n-1}\cdot\frac{s}{\sqrt{n}}
$$

以概率 $1 - \alpha$ 覆盖 $\mu$。这就是 $\sigma$ 未知时 population mean 的 $100(1-\alpha)\%$ confidence interval。注意公式不是假设而是结论：它由 $t$ distribution 的定义和一次不等式变形唯一确定。

### 4.3 $\sigma$ 已知：Normal method

当 $\sigma$ 已知时无需用 $s$ 替代，§3.1 开头的统计量精确服从 standard Normal。用 $z_{\alpha/2}$ 替换 $t_{\alpha/2,\,n-1}$、用 $\sigma_{\bar{x}}$ 替换 $s_{\bar{x}}$，重复 §4.2 的三步推导，得

$$
\bar{x} \pm z_{\alpha/2}\cdot \sigma_{\bar{x}} = \bar{x} \pm z_{\alpha/2}\cdot\frac{\sigma}{\sqrt{n}}.
$$

**Example 4.2（birth weights）.** 某地新生儿出生体重的 $\sigma = 360$ 已知（这是 *population* 之值），某医院 $n = 540$ 例的 sample 给出 $\bar{x}$ 与 $s$。因 $\sigma$ 已知，用 Normal method，$\sigma_{\bar{x}} = 360/\sqrt{540}$；sample 的 $s$ 不参与计算。务必先看清题目给的是 population 还是 sample 的 standard deviation——这是此处最经典的错误。

### 4.4 解读与评价

**解读.** 95% CI 是关于*程序*的概率陈述：重复抽样 $100$ 次、每次各算一个区间，则约 $95$ 个区间包含 $\mu$，$5$ 个不包含。对于已算出的某一个具体区间，$\mu$ 要么在其中、要么不在——"95%"描述的是产生区间的程序，而非该区间的运气。

**Accuracy 与 precision.** 两条评价标准相互牵制：

- *Accuracy*（即 confidence level $1-\alpha$）：提高 confidence level 会使区间变宽（因 $t_{\alpha/2}$ 随 $\alpha$ 减小而增大，见 §3.2 性质 4）、包含 $\mu$ 的把握更大。推至极限，100% CI 是 $(-\infty, +\infty)$：绝对保险，也绝对无用。
- *Precision*（即区间宽度）：区间越窄，估计越精细。

仅通过 $\alpha$ 无法同时改善两者。唯一的例外：在 confidence level *固定*时增大 sample size $n$。由 $s_{\bar{x}} = s/\sqrt{n}$（§2.2），区间宽度 $2 \cdot t_{\alpha/2,\,n-1} \cdot s/\sqrt{n}$ 随 $n$ 增大而缩小，precision 提高而 accuracy 不受损失。

---

## 5 Binomial Distribution

以下两节转入 discrete random variable（Definition 0.1）：概率不再由面积而由逐点的 probability mass function 给出。

**Definition 5.1.** *Bernoulli trial* 有三条性质：(i) 结果互斥（success/failure）；(ii) 每次试验中 success 的 probability $\pi$ 恒定；(iii) 各次试验相互独立。$n$ 次独立 Bernoulli trials 中 success 出现的次数 $X$ 服从 *Binomial distribution*，记作 $X \sim B(n, \pi)$。

**Proposition 5.2（Binomial 的 probability mass function）.**

$$
P(X = k) = \binom{n}{k} \pi^{k} (1 - \pi)^{n - k}, \qquad k = 0, 1, \dots, n.
$$

*理由.* $n$ 次试验中指定某 $k$ 次为 success：由 (ii)(iii)，任一指定序列的概率为 $\pi^k(1-\pi)^{n-k}$；这样的序列共 $\binom{n}{k}$ 个，由 (i) 它们互斥，概率相加即得。

**Example 5.3.** 小白鼠接受某剂量毒物后死亡概率 $\pi = 0.8$；生与死互斥，每只鼠的 $\pi$ 相同，彼此不受影响，三条性质齐备。取 $n = 3$ 只，$X \in \{0,1,2,3\}$ 表示死亡只数，各概率由 Proposition 5.2 给出——它们恰是二项式 $[\pi + (1-\pi)]^3$ 展开的各项，这正是"Binomial"之名的由来。

**Properties.**

- Mean $\mu = n\pi$（success 的平均次数：$3 \times 0.8 = 2.4$ 只）；variance $\sigma^2 = n\pi(1 - \pi)$。
- Variance 与 $n$、$\pi$ 均有关：$n$ 越大，可能结果的范围越广；$\pi$ 越接近 $0.5$ 结果越难预料，variance 越大——因为 $n\pi(1-\pi)$ 在 $\pi = 0.5$ 处取最大值。
- 图形：$\pi = 0.5$ 时对称；$\pi$ 离 $0.5$ 越远越偏斜；固定 $\pi$ 时 $n$ 越大越趋对称。
- *Normal approximation*：当 $n$ 足够大且 $\pi$ 不接近 $0$ 或 $1$ 时，$B(n, \pi)$ 近似 $N(n\pi,\, n\pi(1-\pi))$。

## 6 Poisson Distribution

**Definition 6.1.** *Poisson distribution* 刻画"机会极多而单次概率极小"的稀有事件——$\pi$ 很小而 $n$ 很大——如遗传缺陷、癌症发病、$1$ ml 水中 *E. coli* 的计数。记作 $X \sim P(\lambda)$，唯一参数 $\lambda = n\pi$ 是平均发生次数。其 probability mass function 为

$$
P(X = k) = \frac{\lambda^{k} e^{-\lambda}}{k!}, \qquad k = 0, 1, 2, \dots
$$

**Properties.**

1. *Mean 等于 variance*：$E(X) = \mathrm{Var}(X) = \lambda$。这是它的标志性特征。
2. *与 Binomial 的关系*：当 $n \to \infty$ 且 $\pi$ 很小时，$n\pi(1-\pi) \to n\pi$，Binomial 的 variance 趋近 Poisson 的 variance——Poisson 可视为 Binomial 在此极限下的化身。当二者 mean 相等（$\mu = n\pi$）时，Binomial 的 variance $\mu(1-\pi)$ 严格小于 Poisson 的 variance $\mu$。
3. *Normal approximation*：当 $\lambda \geq 20$ 时可用 Normal 近似——但仅仅是近似，Poisson 本质上仍是离散的。
4. *可加性*：相互独立的 Poisson 变量之和仍是 Poisson（五份 $1$ ml 水样合并为 $5$ ml 后仍服从 Poisson distribution）；这一性质使得向 Normal 问题转化的做法十分便利。

---

## Summary

本讲从一个公理化的起点（Definition 0.2：面积即概率）出发，每一步都有前文支撑：

- **Standardisation**（Theorem 1.2 + Proposition 1.3）把整个 Normal 族统一到一张可查的曲线，且保证概率在变换下不变；
- **Standard error**（Theorem 2.2 + Definition 2.4）量化 sampling variation；$\sigma$ 未知、以 $s$ 替代，这一替代把统计量带离 standard Normal、带上 $t$ distribution（Definition 3.1）；
- $t$ 的概率陈述（§3.3）经一次不等式变形**推出** confidence interval（§4.2），公式无需记忆；
- **Binomial 与 Poisson** 处理 discrete 情形，在 $n$ 或 $\lambda$ 足够大时均可作 Normal approximation。

**Looking ahead.** 把 $0.05$ 的尾部面积重新解读为"出现如此极端的 $T$ 值的概率"，这正是 $P$-value；hypothesis testing（$t$ test）下一讲开始。

---

*转写术语对照（供阅读中文逐字稿者参考）："标准物" = standard error（标准误）；"被努力实验" = Bernoulli trial（伯努利实验）；"派" = $\pi$；"1%减阿尔法" = $1 - \alpha$；"英国统计学家" = Gosset（"Student"）。*
