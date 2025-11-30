---
layout: post
title: 数值积分器（序）
author: myx
categories:
- 数值算法
- 数值积分器
math: true
date: 2025-11-25 18:30 +0800
---
<!-- 数值积分器介绍 -->
## **简介**
<div class="header">
由于大多数复杂系统（包括动力系统）不存在简单封闭形式解析解，因此<strong>数值积分器</strong> （Numerical Integrators）一直是求解描述系统常微分方程（组）<strong>ODEs</strong> (Ordinary Differential Functions)有初值问题（Initial Value Problems，IVPs）最为行至有效的方式。
</div>

### **公式记法和数值本质**
对于常见一阶ODEs系统：

$$\begin{equation}
  y'=f(t,y) \label{Eq.1st-order-odes}
\end{equation}$$
或者天体力学（其他动力系统中）常见的二阶ODEs系统：

$$\begin{equation}
  y''=f(t,y) \label{Eq.2nd-order-odes}
\end{equation}$$
通常我们将方程\eqref{Eq.1st-order-odes}或者\eqref{Eq.2nd-order-odes}中**真实系统状态量流**记作$$\phi_{t}(y_{0})=y(t)$$，将对应当前$$t_{n}$$时刻积分格式记作$$y_{n}$$，对应$$t_{n+1}$$时刻积分格式记作$$y_{n+1}$$。由于数值积分算法其本质上是**利用离散格式对系统的连续微分变化进行逼近**，故存在重要参数即相邻积分时刻（节点）间对应间隔，即积分步长$$h=t_{n+1}-t_{n}$$。而对应离散格式映射记作$$\Phi_{h}:y_{n}\rightarrow y_{n+1}$$。

此时，我们有衡量$$p$$-阶积分算法局部精度的量**局部截断误差LTE** (Local Truncation Error)：

$$\begin{equation}
  LTE=\tilde{y}(h)=y(t_{n+1})-y_{n+1}=\mathcal{O}(h^{p+1})
\end{equation}$$
从小参数动力系统求解的角度看，数值积分器本质上也是对于真实解的有限介绍截断近似，而构造解是基于通用数值分析方法，所以适用于各种问题。无论是基于积分节点构造有限阶数多项式积分近似、等间距有限阶数牛顿差分构造、哈密顿系统辛算法构造近似相流，都会由于有限阶数截断导致LTE产生，当然后期发展出保辛/时间对称的几何特征积分器（Geometry Numerical Integrators）可以有效抑制LTE的长期累积，不过均有各自所适用或者表现更好的场景。对于问题是否刚性也需要选择合适积分器处理不同时间尺度变化率的问题。因此，正如REBOUND[^rebound]的主要开发者Rein[^rein]曾经说过：“并不存在可以解决所有问题完美的数值积分器，更重要的是根据不同问题特性选择最合适的数值积分器”。

## **积分器分类**
<div class="content-box">
数值积分器是一个庞大的算法家族，因此针对不同性质的ODEs系统前人研究中开发了通用性/针对性的数值架构。在分类学（Taxonomy）上想要对数值积分器进行细化分类也是较为复杂的工作。这里对于积分器自身算法性质进行一些罗列。
</div>

<div class="list-box">
<ul>
  <li><strong>单步法/多步法</strong>：单步法（如Runge-Kutta方法等）通常仅依赖积分区间$[t_{n},t_{n+1}]$内特定节点变化率信息构造近似多项式积分形式，即前向积分不会依赖$y_{n-i}(i=1,2,3\cdots)$的节点信息。相应的，多步法（如Adams方法等）通常需要依赖前向节点。因此一个特点就是多步法在代码实现时<strong>通常需要非局部变量的额外存储空间</strong>。</li>
  <li><strong>变步长/定步长</strong>：通常可否变阶次或变步长并不是一个积分算法分类的依据或者标准。对于常用积分方法，通常可以采用步长减半$(h\rightarrow \frac{h}{2})$检测LTE是否在误差容限（Tolerance，TOL）内来进行局部积分精度是否满足要求的判断；或者利用相同积分器族内的高阶方法嵌入（embedding）获得的LTE提供变阶次或者变步长依据（如Runge-Kutta-Felberg 7[8]方法）；此外，还可以通过特定积分构造的LTE公式分析进行自适应步长调整策略。变步长更多是一种<strong>积分策略</strong>，选择不同的积分策略不仅需要考虑积分器自身数值架构，也需要考虑处理的问题刚性/非刚性，或者动力系统是否存在奇点等。</li>
  <li><strong>显式/隐式</strong>：显示方法在进行积分时，不需要提前获取积分节点应有的准确信息；相应地，隐式方法需要未知积分节点的理论准确信息来构造特定阶数解。显而易见的是，隐式方法通常在使用时需要进行<strong>迭代（Iteration）或者预估校正（PECE）</strong>，因此影响使用效率；但是相应地，隐式方法通常可以额外获得<strong>数值稳定/时间对称</strong>等良好特质。</li>
  <li><strong>保辛/非保辛</strong>：针对哈密顿动力系统（Hamilton）在相空间中保面积的性质，学者研究了特定的保辛（Symplectic）结构积分器并且进行推广。即便是截断到特定阶数的近似哈密顿解，得益于保辛特性对于相空间几何特征的保持，LTE并不会发生长期累积，误差项通常表现为有界震荡（bounded oscillation）。因此动力系统在位形空间中的动力学特征在长期积分下得以保持，不会出现显著“质变”。例如，对于低阶算法而言，二体椭圆轨道会因为LTE累积呈现内/外螺旋特性（Inward/Outward Spiral）；而辛算法可以对轨道特征进行长期保持。</li>
  <li><strong>时间对称/非对称</strong>：通俗讲，时间对称（Time-symmetric/reversible）积分器的积分节点映射流对于步长符号而言完全可逆。即经过离散数值映射$\Phi_{h}:y_{n}\rightarrow y_{n+1}$和其反向映射$\Phi_{-h}:y_{n+1}\rightarrow y_{n}$满足条件$\Phi_{h}=\Phi_{-h}^{-1}$。由于积分映射流完全可逆，因此也可以很大程度在数值上抑制误差累积，获得较好的特性保持效果。</li>
</ul>
</div>

## **主要积分方法罗列**

| Scheme | Symplectic | Time-Symmetric | Force Evaluations | Variable-Step | ODEs Order |
|:-----|:----:|:----:|:----:|:----:|:----:|
| **经典单步法积分器** |
| 显式Euler | ❌ | ❌ | 1 | ✅[^note1] | 1 |
| 隐式Euler | ❌ | ❌ | 1[^note2] | ✅[^note1] | 1 |
| 梯形法则  | ❌ | ✅ | 1[^note2] | ✅[^note1] | 1 |
| 隐式中点法 | ✅ | ✅ | 1[^note2] | ✅[^note1] | 1 |
| Runge-Kutta[p] | ❌ | ❌ | $\geq p$ | ✅[^note1] | 1 |
| Runge-Kutta-Felberg[p/p+1] | ❌ | ❌ | $\ge p$ | ✅ | 1 |
| **经典多步法积分器** |
| Adams-Bashforth[p]| ❌ | ❌ | 1 | ❌ | 1 |
| Adams-Moulton[p] | ❌ | ❌ | 1[^note2] | ❌ | 1 |
| Stormer[p] | ❌ | ❌ | 1 | ❌ | 2 |
| Cowell[p] | ❌ | ❌ | 1[^note2] | ❌ | 2 |
| Krogh-Shampine-Gordon[p] | ❌ | ❌ | 1 | ❌ | 1/2 |
| Summed Adams[p] | ❌ | ❌ | 1 | ❌ | 1 |
| Gauss-Jackson[p] | ❌ | ❌ | 1 | ❌ | 2 |
| Symmetric Multistep[p] | ❌ | ❌ | 1 | ❌ | 2 |
| Shampine-Gordon | ❌ | ❌ | Varies | ❌ | 1 |
| **Gaussian配点法（Collocation）** |
| Gauss-Legendre Implicit RK[p] | ✅ | ✅ | $\frac{p}{2}$[^note2] | ✅[^note1] | 1 |
| Gauss-Radau Spacing[p] | ❌ | ❌ | $\frac{p+1}{2}$[^note2] | ✅[^note1] | 1 |
| IAS15[^IAS15] | ❌ | ❌ | $\frac{p+1}{2}$[^note2] | ✅[^note1] | 1 |
| **经典基于Richardon外推的方法** |
| Bulirsch-Stoer | ❌ | ❌ | Varies | ❌ | 1 |
| **经典构造法（Composition）**|
| 保辛Euler | ✅ | ❌ | 1 | ✅[^note1] | 1 |
| Verlet/Leapfrog | ✅ | ✅ | 1 | ✅[^note1] | 2 |
| SABA[s] | ✅ | ❌ | $\geq s$ | ❌ | 2 |
| **基于Wisdom-Holman架构的（混合）方法** |
| Wisdom-Holman (WH)[^WH] | ✅ | ❌ | 1 | ❌ | 2 |
| RMVS[^RMVS] | ✅ | ❌ | 1 | ❌ | 2 |
| MERCURY[^Mer] | ✅ | ❌ | 1 | ❌ | 2 |
| TRACE[^TRACE] | ❌ | ✅ | 1 | ❌ | 2 |

[^note1]: Via step doubling or embedded methods
[^note2]: Implicit method (iterative solution required)

<div class="content-box">
后续会逐渐根据Table所示主要积分器分类进行逐一介绍，主要基于日常科研使用、文献阅读和个人理解。对于不同积分架构通过原理、代码实现、以及具体算例进行相关呈现。
</div>



## 参考文献和注解
[^rebound]: [REBOUND GitHub Repository](https://github.com/hannorein/rebound)
[^rein]: [REBOUND 主页](https://rebound.hanno-rein.de/)
[^IAS15]: [IAS15 文章](https://academic.oup.com/mnras/article/446/2/1424/2892331?login=true)
[^WH]: [Wisdom-Holman 文章](https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract)
[^RMVS]: [RMVS 文章](https://www.sciencedirect.com/science/article/pii/S0019103584710396)
[^Mer]: [Mercury 文章](https://ui.adsabs.harvard.edu/abs/1999MNRAS.304..793C/abstract)
[^TRACE]: [TRACE 文章](https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.3708L/abstract)


<link rel="stylesheet" href="{{ '/assets/css/mystyle.css' | relative_url }}">