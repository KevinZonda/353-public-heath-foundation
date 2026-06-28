# ML Adapter

> **机器学习**里的概率论更像是在问：
> "假设数据生成过程是这样的，模型怎么学到分布？"  
> **流行病学**里的统计更像是在问：
> "现实数据已经烂成这样了，我怎么还能尽量不说瞎话？"

## Term Conversion

| Term in PH | Term in ML |
| --- | --- |
| 暴露因素 | Feature |
| 结局 | Label |
| 偏倚  | Dataset Shift |
| RR | effect size |
| OR | odds transform |
| 筛检灵敏度 | Recall |
| 筛检特异度 | TNR/Specificity |
| 阳性预测值 | Precision/PPV |
| 标准化率 | Reweighting |
| 误报率 | FPR |
| ROC | ROC 曲线 |
| AUC | AUC 曲线 |

## Confusion Matrix

![](img/confmat.png)

$$
\begin{align*}
Recall &= \frac{TP}{TP+FN}
\\
Precision &= \frac{TP}{TP+FP}\\
TNR &= \frac{TN}{TN+FP}
\end{align*}
$$

$$
\begin{align*}
Recall &= p(\text{pred} = +\ |\ \text{actual} = +)
\\
Precision &=
p(\text{actual} = + | \ \text{pred} = +)
\\
TNR &= p(\text{pred} = - \ |\ \text{actual} = -)
\\
FPR &= p(\text{pred} = + \ |\ \text{actual} = -) = 1-TNR
\end{align*} \\
$$

## OR（Odds Ratio）

### Odd

$$
odd = \frac{p}{1-p} = \frac{p(+)}{p(-)}
$$

odd（优势），有如下不对称性：
- =1: 势均力敌，发生与不发生概率 1:1
- \>1: 发生概率 高于 不发生
  - \>\>1: 大幅领先，odd增长不是线性的，而是爆炸性的（趋近于无穷大）
- <1: 发生概率 低于 不发生
  - -> 0: 大幅落后，odd衰减不是线性的，而是趋近0

logit（log odd），“势均力敌”变成了 $0$，“大幅领先”变成了正无穷，“大幅落后”变成了负无穷，实现对称

$$
OR = Odds_1 : Odds_2
$$

> 假设 $Odds_1$ 是实验组，$Odds_2$ 是对照组。

- $OR = 1$: 势均力敌（两组没区别）
  - 含义：$Odds_1 = Odds_2$。
  - 直觉：实验组的领先优势和对照组一模一样。这意味着你吃不吃这个药，康复的几率没有任何改变。这个因素（吃药）对结果（康复）毫无影响。
- $OR \gg 1$（远大于 1）时：实验组大幅领先
  - 含义：$Odds_1$ 远大于 $Odds_2$。
  - 直觉：实验组的胜算几率压倒性地高于对照组。
  - 举个例子：假设对照组（不抽烟）的患癌 $Odds_2 = 0.01$（极低落后）；而实验组（重度抽烟）的患癌 $Odds_1 = 0.15$。此时的$$OR = \frac{0.15}{0.01} = 15$$这说明：抽烟组的患癌几率，是不抽烟组的 15 倍！（抽烟属于强烈的危险因素，抽烟组在患癌风险上大幅领先）。
- $OR < 1$（接近 0）：实验组大幅落后（即对照组大幅领先）
  - 含义：$Odds_1$ 远小于 $Odds_2$。
  - 直觉：实验组的发生几率被严重压缩，说明该因素起到了“保护”作用。
  - 举个例子：如果计算出来经常运动的人，其患心脏病的 $OR = 0.3$。这意味着运动组的患病几率只有不运动组的 $30\%$（或者说，不运动组的患病几率是运动组的 $\frac{1}{0.3} \approx 3.33$ 倍）。

## 