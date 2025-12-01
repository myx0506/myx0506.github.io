---
layout: post
title: 数值积分器：数值微分、定积分和IVP的关系
author: myx
categories:
- 数值算法
- 数值积分器
math: true
date: 2025-11-30 22:44 +0800
---
## **简介**
<div class="header">
无论是<strong>数值微分</strong>、<strong>定积分问题</strong>还是<a href="{{ '/posts/numerical-integrator-intro/' | relative_url }}" class="internal-link">数值积分器（序）</a>中提到的IVP问题，解决问题的本质是利用离散的节点信息完成对于连续函数（未必存在封闭解析形式）的复现和近似，在数值分析的基础上完成一定阶数的近似拟合。
</div>

## **数值微分、定积分问题和IVP问题的关系**
### **数值微分**
通常，对于给定函数形式 
$$\begin{equation}
  y=f(x) \label{Eq.function}
\end{equation}$$

而言，无论其自身存在解析形式还是离散计算形式，在数值计算方法中可以依赖前后间隔长度$h$信息点，采用**前向差分**、**后向差分**以及**中间差分**格式构造特定点$x=x_{0}$处的一阶甚至高阶导数近似$y^{(n)}(x_{0})$。这里数值微分的本质是利用等间隔信息点的泰勒展开（Taylor exansion）级数：
$$\begin{equation}
\begin{aligned}
 f(x_{0}+h)&=\sum_{i=0}^{\infty}\frac{h^{i}}{i!}f^{(i)}(x_{0}) \\
 &=\underbrace{f(x_{0})+hf'(x_{0})+\frac{h^{2}}{2!}f''(x_{0})+\cdots+\frac{h^{n}}{n!}f^{(n)}(x_{0})}_{\text{n terms}}+\mathcal{O}(h^{n+1})
\end{aligned}
\label{Eq.Taylor Expansion}
\end{equation}$$
进行$x=x_{0}$处的特定阶导数差分构造。

- #### **前向差分**

对于最基础的前向差分方法，可以采用导数的定义并且进行递推，得到函数n阶导数的差分形式：
  $$\begin{equation}
  \begin{aligned}
    f'(x_{0})&\approx\frac{f(x_{0}+h)-f(x_{0})}{h} \\
              &=\frac{f(x_{0}+h)-f(x_{0})}{h}-\frac{h}{2!}f''(x_{0})-\frac{h^{2}}{3!}f'''(x_{0})+\cdots \\
              &=\frac{f(x_{0}+h)-f(x_{0})}{h}+\mathcal{O}(h) \\
    f''(x_{0})&\approx\frac{f'(x_{0}+h)-f'(x_{0})}{h} \\
              &\approx\frac{\frac{f(x_{0}+2h)-f(x_{0}+h)}{h}-\frac{f(x_{0}+h)-f(x_{0})}{h}}{h} \\
              &\approx\frac{f(x_{0}+2h)-2f(x_{0}+h)+f(x_{0})}{h^{2}} \\
              &=\frac{f(x_{0}+2h)-2f(x_{0}+h)+f(x_{0})}{h^{2}}-hf'''(x_{0})-\frac{7}{12}h^{2}f^{(4)}(x_{0})+\cdots \\
              &=\frac{f(x_{0}+2h)-2f(x_{0}+h)+f(x_{0})}{h^{2}}+\mathcal{O}(h) \\
    \vdots & \\
    f^{(n)}(x_{0})&\approx\frac{\Delta^{n}f(x_{0})}{h^{n}} \\
                  &\approx\frac{\sum_{k=0}^{n}(-1)^{n-k}\binom{n}{k}f(x_{0}+kh)}{h^{n}} \\
                  &=\frac{\sum_{k=0}^{n}(-1)^{n-k}\binom{n}{k}f(x_{0}+kh)}{h^{n}}+\mathcal{O}(h)
  \end{aligned}
  \label{Eq.Forward_Difference}
  \end{equation}$$
注意到，基于前向差分公式的函数微分近似和真实导数的差别是$\mathcal{O}(h)$。可以看到这里分子的形式就是**n阶前向差分算子$\Delta^{n}$**，其展开形式与二次项系数相关。如果令$f_{k}=f(x_{0}+kh)$，那么差分算子定义为：
  $$\begin{equation}
    \Delta^{n}f_{k}=\Delta^{n-1}f_{k+1}-\Delta^{n-1}f_{k} \label{Eq.Differential}
  \end{equation}$$
从算子的角度看可以推导\eqref{Eq.Differential}定义的前向差分算子具有如\eqref{Eq.Forward_Difference}的形式。定义算子$I$为不变算符$If_{k}=f_{k}$，算子$E$为前向算符$Ef_{k}=f_{k+1}$，那么一阶差分算符可以写成$\Delta=E-I$,且有：
  $$\begin{equation}
  \begin{aligned}
    \Delta^{n}f_{0}&=(E-I)^{n}f_{0} \\
                   &=(-I+E)^{n}f_{0} \\
                   &=\left(\sum_{k=0}^{n}(-1)^{n-k}\binom{n}{k}E^{k}\right)f_{0} \\
                   &=\sum_{k=0}^{n}(-1)^{n-k}\binom{n}{k}f(x_{0}+kh)
  \end{aligned}
  \end{equation}$$
当然，用前向差分获得函数微分近似也可以从根本的Newton插值多项式进行理解。对于$n+1$个插值点$(x_{i},y_{i})$，Newton插值多项式利用逐个点添加以后的斜率修正进行$n$阶递推插值，有一般形式：
  $$\begin{equation}
    \begin{aligned}
      P_{n}(x)&=f[x_{0}] \\
              &+f[x_{0},x_{1}](x-x_{0}) \\
              &+f[x_{0},x_{1},x_{2}](x-x_{0})(x-x_{1}) \\
              &+\cdots \\
              &+f[x_{0},\cdots,x_{n}](x-x_{0})\cdots(x-x_{n-1})
    \end{aligned}
    \label{Eq.Newton Interpolation}
  \end{equation}$$
其中$f[x_{0},\cdots,x_{n}]$是**差商**，其存在递推关系：
  $$\begin{equation}
    \begin{aligned}
      &f[x_{0}]=f(x_{0}) \\
      &f[x_{0},x_{1}]=\frac{f(x_{1})-f(x_{0})}{x_{1}-x_{0}} \\
      &f[x_{0},x_{1},x_{2}]=\frac{f[x_{1},x_{2}]-f[x_{0},x_{1}]}{x_{2}-x_{0}} \\
      &\vdots \\
      &f[x_{0},x_{1},\cdots,x_{k}]=\frac{f[x_{1},x_{2},\cdots,x_{k}]-f[x_{0},x_{1},\cdots,x_{k-1}]}{x_{k}-x_{0}}
    \end{aligned}
    \label{Eq.Difference quotient}
  \end{equation}$$
从多项式插值角度看，本质上，各阶次差商代表了添加信息点$(x_{i},y_{i})$后的高阶信息，因此如\eqref{Eq.Forward_Difference}给出的函数微分形式完全可以由\eqref{Eq.Newton Interpolation}的多项式形式给出。将$P_{n}(x)$看作$n+1$个点插值得到的函数$f(x)$的$n$阶多项式近似，并且让离散数据点$(x_{i},y_{i})$对应前向差分的等距情况，则**差商**存在简单形式：
  $$\begin{equation}
    f[x_{0},x_{1},\cdots,x_{k}]=\frac{\Delta^{k}f(x_{0})}{k!h^{k}} \label{Eq.Equal Difference quotient}
  \end{equation}$$
此时\eqref{Eq.Newton Interpolation}可以改写成为：
  $$\begin{equation}
    P_{n}(x)=f(x_{0})+\frac{\Delta f(x_{0})}{h}(x-x_{0})+\frac{\Delta^{2}f(x_{0})}{2!h^{2}}(x-x_{0})(x-x_{1})+\cdots \label{Eq.Equal Newton Interpolation}
  \end{equation}$$
而函数$f(x)$和依赖$n+1$个点进行插值的多项式误差余项应该有：
  $$\begin{equation}
    f(x)-P_{n}(x)=\frac{\Delta^{k+1}f(x_{0})}{(k+1)!}(x-x_{0})(x-x_{1})\cdots(x-x_{k})=\mathcal{O}(h^{n+1})
    \label{Eq.Error}
  \end{equation}$$
这时结合\eqref{Eq.Taylor Expansion}的$f(x)$在$x=x_{0}$附近泰勒级数展开形式，\eqref{Eq.Forward_Difference}的数值微分关系完全由对应项数$(n+1)$的$P_{n}(x)$对应最高阶次斜率导出，即：
  $$\begin{equation}
  \begin{aligned}
    f'(x_{0})&\approx P_{1}'(x_{0})+\mathcal{O}(h)\approx\frac{\Delta f(x_{0})}{h}+\mathcal{O}(h) \\
    f''(x_{0})&\approx P_{2}''(x_{0})+\mathcal{O}(h)\approx\frac{\Delta^{2} f(x_{0})}{h^{2}}+\mathcal{O}(h) \\
    &\vdots \\
    f^{(n)}(x_{0})&\approx P_{n}^{(n)}(x_{0})+\mathcal{O}(h)\approx\frac{\Delta^{n} f(x_{0})}{h^{n}}+\mathcal{O}(h)
  \end{aligned}
  \label{Eq.Forward_Difference_from_P}
  \end{equation}$$
这里实际上利用了Newton插值多项式和待拟合函数泰勒展开的内禀联系，也解释了对于函数$f(x)$的n阶导数，采用$n+1$个前向差分点构造近似微分为什么只能达到$h$的一阶精度。**从多项式拟合来看本质采用$n+1$个信息点插值获取$n$阶导数信息，截断误差自然截止到目标的一阶精度**。

而从Newton插值多项式出发，自然我们可以想办法构造$y=f(x)$的$n$阶导数的任意阶次精度，**只要添加前向差分的等距点即可**。以函数$f(x)$的一阶微分为例，由\eqref{Eq.Equal Newton Interpolation}不难看出，若利用$P_{n}(x)$的一阶微分近似$f'(x_{0})$，我们完全可以通过提高阶数$n$的方式得到更加精确的结果：
  $$\begin{equation}
    f'(x_{0})\approx P_{n}(x_{0})\approx\frac{\Delta f(x_{0})}{h}-\frac{\Delta^{2}f(x_{0})}{2h}+\frac{\Delta^{3}f(x_{0})}{3h}+\cdots
    \label{Eq.12}
  \end{equation}$$
> **$n=2$**: 
> $$f'(x_{0})\approx P_{2}'(x_{0})\approx\frac{1}{2h}\left(-f_{2}+4f_{1}-3f_{0}\right)+\mathcal{O}(h^{2})$$
{:.prompt-info}

该数值微分任意阶次的推导仍然可以用算符关系得到，出发点仍然是泰勒级数\eqref{Eq.Taylor Expansion}，定义微分算符$D$：
  $$\begin{equation}
    f(x_{0}+h)=\sum_{i=0}^{\infty}\frac{h^{i}}{i!}f^{(i)}(x_{0})=\sum_{i=0}^{\infty}\frac{h^{i}D^{i}}{i!}f(x_{0})=e^{hD}f(x_{0})
  \end{equation}$$
则有关系$$\Delta=e^{hD}-I$$，经过简单推导有：
  $$\begin{equation}
    \begin{aligned}
      hD&=\ln{(\Delta+I)}=\Delta-\frac{\Delta^{2}}{2}+\frac{\Delta^{3}}{3}+\cdots \\
      D&=\frac{1}{h}\left(\Delta-\frac{\Delta^{2}}{2}+\frac{\Delta^{3}}{3}+\cdots \right)
    \end{aligned}
    \label{Eq.14}
  \end{equation}$$
不难看出这和\eqref{Eq.12}给出的形式完全一致。

至此，我们从**多项式插值**和**算符推导**理解了**前向差分**形式下如何推导函数任意阶精度到间隔$h$任意阶次的微分形式。除了上述两个角度，也可以从\eqref{Eq.Taylor Expansion}在不同节点展开联立方程组待定系数的方式看待数值微分问题。如果我们用$n+1$个节点来构造$f'(x_{0})$，那么有以下形式：
  $$\begin{equation}
    f'(x_{0})=\frac{1}{h}\sum_{i=0}^{n}\alpha_{i}f_{i}
  \end{equation}$$
当$n=1$时，可以有方程组：
  $$\begin{equation}
  \left\{
  \begin{aligned}
    \alpha_{0} + \alpha_{1} &= 0 \\
    \alpha_{1} &= 1
  \end{aligned}
  \right.
  \end{equation}$$
当$n=2$时，可以有方程组：
  $$\begin{equation}
  \left\{
  \begin{aligned}
    \alpha_{0}+\alpha_{1}+\alpha_{2}&=0 \\
    \alpha_{1}+2\alpha_{2}&=1 \\ 
    \frac{1}{2}\alpha_{1}+2\alpha_{2}&=0
  \end{aligned}
  \right.
  \end{equation}$$
我们同样可以得到$f'(x_{0})$的三点前向差分公式。利用$n+1$个点，最多可以求解$n+1$个独立方程，自然可以将数值微分的近似精度提高到$\mathcal{O}(h^{n})$。

当$n=k$时，可以有方程组：
  $$\begin{equation}
  \left\{
  \begin{aligned}
    \sum_{i=0}^{k}i^{0}\alpha_{i}&=0 \\
    \sum_{i=0}^{k}i^{1}\alpha_{i}&=1 \\ 
    \sum_{i=0}^{k}i^{2}\alpha_{i}&=0 \\ 
    &\vdots \\
    \sum_{i=0}^{k}i^{k}\alpha_{i}&=0
  \end{aligned}
  \right.
  \label{Eq.18}
  \end{equation}$$

这里的三个思想：
> - 利用$n+1$个点进行$n$阶多项式构造
> - 算符推导
> - 利用$n+1$个离散点级数展开联立多项式待定系数
{:.prompt-tip}
是后续进行定积分和IVP求解问题中数值方法构造的基本思想，包括如何利用$n+1$个点进行更高阶次代数精度构造（Gauss/Lobatto Collocation）的出发点。

> 需要说明的是，如果待估函数本身就是$n$阶多项式形式，那么需要最多不超过$n+1$个点就可以得到多项式的全部信息（即多项式插值本质）。从方程组待定系数角度看，由于高阶导数($>n$)均为0，待定系数独立方程个数不会大于$n+1$个。
{:.prompt-info}

- #### **后向差分**
后向差分和前向差分本质完全一致，只需要将$(0,\cdots,k)$的对应index起点放到$k$即可，仅存在表示形式的差别。同时，由于后向差分的往期信息依赖特性，其也是**Adams-Bashforth**为首的**多步法积分器**的重要基础。

- #### **中间差分**
顾名思义，中间差分的基本数值格式依赖与$x=x_{0}$想对称的前后$m$个节点信息。对于函数$f(x)$而言，其位于$x=x_{0}$处任意阶微分的基本中间差分形式为：
  $$\begin{equation}
  \begin{aligned}
    f'(x_{0})&\approx\frac{f(x_{0}+\frac{h}{2})-f(x_{0}-\frac{h}{2})}{h}+\mathcal{O}(h^{2}) \\
             &\approx\frac{f(x_{0}+h)-f(x_{0}-h)}{2h}+\mathcal{O}(h^{2}) \\
    f''(x_{0})&\approx\frac{f'(x_{0}+\frac{h}{2})-f'(x_{0}-\frac{h}{2})}{h} \\
              &\approx\frac{\frac{f(x_{0}+h)-f(x_{0})}{h}-\frac{f(x_{0})-f(x_{0}-h)}{h}}{h} \\
              &\approx\frac{f(x_{0}+h)-2f(x_{0})+f(x_{0}-h)}{h^{2}} \\
              &=\frac{f(x_{0}+h)-2f(x_{0})+f(x_{0}-h)}{h^{2}}+\mathcal{O}(h^{2}) \\
    \vdots & \\
    f^{(n)}(x_{0})&\approx\frac{\delta^{n}f(x_{0})}{h^{n}} \\
                  &\approx\frac{\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}f\left(x_{0}+(\frac{n}{2}-k)h\right)}{h^{n}} \\
                  &=\frac{\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}f\left(x_{0}+(\frac{n}{2}-k)h\right)}{h^{n}}+\mathcal{O}(h^{2})
  \end{aligned}
  \label{Eq.Central_Difference}
  \end{equation}$$
由于其节点相关$x=x_{0}$的对称特性，从\eqref{Eq.Taylor Expansion}中泰勒级数展开的角度出发，可以看到偶数项由于线性组合消除，因此$n$阶微分的精度自然达到了$\mathcal{O}(h^{2})$，比前向/后向差分更高。同样，\eqref{Eq.Central_Difference}中的中间差分算子$\delta^{n}$有类似前向差分算子的二项式系数特征：
  $$\begin{equation}
  \begin{aligned}
    \delta^{n}f_{0}&=\left(E^{\frac{1}{2}}-E^{-\frac{1}{2}}\right)^{n}f_{0} \\
    &=\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}E^{\frac{n-2k}{2}}f_{0} \\
    &=\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}f\left(x_{0}+(\frac{n}{2}-k)h\right)
  \end{aligned}
  \end{equation}$$
类似在**前向差分**部分的推导，我们可以利用添加信息点的方式来构造高阶多项式逼近任意阶微分的任意阶次精度形式（这里中间差分的Newton插值多项式并不直观，或者如\eqref{Eq.14}利用微分算符关系进行递推。同样以一阶微分为例，存在中间差分格式和微分算符关系：
  $$\begin{equation}
    \delta=e^{\frac{h}{2}D}-e^{-\frac{h}{2}D}=2\sinh{\frac{h}{2}D}
  \end{equation}$$
那么一阶微分形式的算子表达为：
  $$\begin{equation}
    hD=2\sinh^{-1}{\frac{\delta}{2}}=\delta-\frac{1}{24}\delta^{3}+\frac{3}{640}\delta^{5}+\cdots \label{Eq.22}
  \end{equation}$$
注意到，这里的中间差分算符$\delta=f(x_{0}+\frac{h}{2})-f(x_{0}-\frac{h}{2})$是半点格式，我们可以改写整点的基本中间差分算符：$\delta'=f(x_{0}+h)-f(x_{0}-h)$，此时\eqref{Eq.22}变成：
  $$\begin{equation}
    hD=\sinh^{-1}{\frac{\delta'}{2}}=\frac{\delta'}{2}-\frac{1}{48}\delta'^{3}+\frac{3}{1280}\delta'^{5}+\cdots \label{Eq.23}
  \end{equation}$$

> 将当前$x=x_{0}$和前后$m$个对称点总共$n=2m+1$个点看成多项式插值节点<br>
> **$m=1$**: $$f'(x_{0})\approx \frac{\delta'}{2h}\approx \frac{f_{1}-f_{-1}}{2h}+\mathcal{O}(h^{2})$$ <br>
> **$m=2$**: $$f'(x_{0})\approx \frac{1}{h}\left(\frac{\delta'}{2}-\frac{\delta^{3}}{48}\right)\approx \frac{-f_{3}+27f_{1}-27f_{-1}+f_{3}}{48h}+\mathcal{O}(h^{4})$$
{:.prompt-info}

> 注意到，当$m=2$时，由微分算符推导的$f(x)$一阶微分的五点四阶精度公式并不是标准中心差分公式：
> $$f'(x_{0})\approx \frac{-f_{2}+8f_{1}-8f_{-1}+f_{-2}}{12h}+\mathcal{O}(h^{4})$$
{:.prompt-warning}

这个问题我们可以从联立多项式待定系数求解角度出发，基于中间差分的对称性，线性方程组的形式可以写成：
  $$\begin{equation}
    f'(x)\approx\frac{1}{h}\sum_{k=-m}^{m}\alpha_{k}f(x+kh) \label{Eq.24}
  \end{equation}$$
其中，对于一阶微分，上述方程组系数满足性质:
 - $\alpha_{-k}=-\alpha_{k}$
 - $\sum_{k=-m}^{m}\alpha_{k}=0$
 - $\sum_{k=-m}^{m}k\alpha_{k}=1$


此时我们可以利用性质对\eqref{Eq.24}进行改写并且代入泰勒级数进行展开：
  $$\begin{equation}
  \begin{aligned}
    f'(x)&\approx\frac{1}{h}\sum_{k=0}^{m}\alpha_{k}\left[f(x+kh)-f(x-kh)\right] \\
    &\approx2\sum_{n odd}\frac{h^{n-1}}{n!}f^{(n)}(x)\left(\sum_{k=1}^{m}\alpha_{k}k^{n}\right)
  \end{aligned}
  \label{Eq.25}
  \end{equation}$$
对于需要达到的$\mathcal{O}(h^{2m})$阶精度，只要对应$m$个独立线性方程组满足要求：
  $$\begin{equation}
  \left\{
  \begin{aligned}
    \sum_{k=1}^{m}2k\alpha_{k}&=1 \\
    \sum_{k=1}^{m}k^{3}\alpha_{k}&=0 \\ 
    \sum_{k=1}^{m}k^{5}\alpha_{k}&=0 \\ 
    &\vdots \\
    \sum_{k=1}^{m}k^{2m-1}\alpha_{k}&=0 \\ 
  \end{aligned}
  \right.
  \label{Eq.26}
  \end{equation}$$
可以看到**微分算符**和**联立方程**得到的高阶公式表达不同本质原因是**中间算符形式选取不同，不同的节点分布密度会导致高阶泰勒项系数不同，一般而言节点分布越密集，带来的局部误差常数越小**，这也符合我们一般的直觉认知（越近的信息越能准确反映局部变率信息）。

- #### **Richardson外推**
<div class="content-box">
无论在<strong>前向差分</strong>、<strong>后向差分</strong>还是<strong>中间差分</strong>中，对于$n$阶微分特定阶次精度近似中，似乎提高精度的方式只有<strong>减小等间隔步长$h$</strong>。而Richardson外推采取了更加巧妙的方式，其核心思想则是：<strong>利用不同步长下渐进展开的低阶次（精度）结果，通过线性组合（插值）逐渐消去较大的误差主项，从而得到更高精度的估算值</strong>。
</div>

由于离散差分方式的微分近似本质上是截断到一定精度的泰勒展开级数的线性组合，误差均可以写成步长$h$的渐进展式：
  $$\begin{equation}
    T(h)=I+A_{1}h^{p}+A_{2}h^{q}+A_{3}h^{r}+\cdots
  \end{equation}$$
这里$T(h)$是目标量的数值近似值（步长$h$的函数），$I$是目标量准确值，$p<q<r<\cdots$是渐进展式的误差阶数。如果我们改变步长$h$构成序列$(h_{i},T(h_{i}))$，那么我们可以利用线性组合消除误差主项$A_{1}h^{p}$，从而提高近似表达式精度。这里以一阶Richardson外插举例：
  $$\begin{equation}
  \begin{aligned}
    T(h)&=I+A_{1}h^{p}+A_{2}h^{q}+A_{3}h^{r}+\cdots \\
    T\left(\frac{h}{2}\right)&=I+A_{1}\left(\frac{h}{2}\right)^{p}+A_{2}\left(\frac{h}{2}\right)^{q}+A_{3}\left(\frac{h}{2}\right)^{r}+\cdots
  \end{aligned}
  \end{equation}$$
显然，我们可以利用线性组合消除$A_{1}h^{p}$来得到更高精度的新待估量：
  $$\begin{equation}
    I\approx T_{1}\left(\frac{h}{2}\right)=\frac{2^{p}\cdot T\left(\frac{h}{2}\right)-T(h)}{2^{p}-1}
  \end{equation}$$
此时，待估量的误差主项阶数升高到$\mathcal{O}(h^{q})$。

> **常见应用**：<br>
> - Romberg积分法是Richardson外推应用于梯形法的算法实现
> - 显著提高数值微分（差分方法）的近似精度
> - 微分方程数值求解/CFD网格收敛性分析
> - **Bulirsch-Store外推方法构造**
{:.prompt-tip}

Richardson外推完全可以添加不同步长信息点，利用一阶外推序列进行二阶外推递推，实现方法非常简单。

| $i$ | $h_i$ | 列0: $O(h^p)$ | 列1: $O(h^{2p})$ | 列2: $O(h^{3p})$ | 列3: $O(h^{4p})$ | 说明 |
|----------|------------|----------------|------------------|------------------|------------------|------|
| 0 | $h$   | $T_{0,0}$ | | | | 原始计算值 |
| 1 | $h/2$ | $T_{1,0}$ | $T_{1,1}$ | | | 第一次外推 |
| 2 | $h/4$ | $T_{2,0}$ | $T_{2,1}$ | $T_{2,2}$ | | 第二次外推 |
| 3 | $h/8$ | $T_{3,0}$ | $T_{3,1}$ | $T_{3,2}$ | $T_{3,3}$ | 第三次外推 |

其中，通用递推公式可以写成：
  $$\begin{equation}
    T_{i,j}=T_{i,j-1}+\frac{T_{i,j-1}-T_{i-1,j-1}}{(h_{i-j}/h_{i})^{p}-1}
  \end{equation}$$
该形式本质上就是多项式插值中的**Neville**算法。以多项式插值的角度来看Richardson外推，其实就是**以$(h_{i}^{p},T(h_{i}))$为数据序列进行多项式插值，而近似待估量就是多项式在$h=0$时的外推值**。**Neville**三角递推表格形式如下：

| $i$ | $x_i$ | $y_i$ |  $P_{i,0}$ | $P_{i,1}$ | $P_{i,2}$ | $P_{i,3}$ | 说明 |
|----------|------------|----------------|------------------|------------------|------------------|------|
| 0 | $x_{0}$ | $y_{0}$ | $P_{0,0}=y_0$ | | | | 通过点$(x_0,y_0)$的0次多项式 |
| 1 | $x_{1}$ | $y_{1}$ | $P_{1,0}=y_1$ | $P_{1,1}$ | | | 通过$(x_0,y_0),(x_1,y_1)$的1次多项式 |
| 2 | $x_{2}$ | $y_{2}$ | $P_{2,0}=y_2$ | $P_{2,1}$ | $P_{2,2}$ | | 通过$(x_0,y_0),(x_1,y_1),(x_2,y_2)$的2次多项式 |
| 3 | $x_{3}$ | $y_{3}$ | $P_{3,0}=y_3$ | $P_{3,1}$ | $P_{3,2}$ | $P_{3,3}$ | 通过所有4点的3次多项式 |

Neville算法通用递推公式可以写成：
  $$\begin{equation}
  \begin{aligned}
    P_{i,j}(x)&=\frac{(x-x_{i-j})P_{i,j-1}(x)-(x-x_{i})P_{i-1,j-1}(x)}{x_{i}-x_{i-j}} \\
    &=P_{i,j-1}(x)+\frac{P_{i,j-1}(x)-P_{i-1,j-1}(x)}{\frac{x-x_{i-j}}{x-x_{i}}-1}
  \end{aligned}
  \end{equation}$$
当$x=0$时，形式为：
  $$\begin{equation}
  \begin{aligned}
    P_{i,j}(0)&=\frac{x_{i}P_{i-1,j-1}(0)-x_{i-j}P_{i,j-1}(0)}{x_{i}-x_{i-j}} \\
    &=\frac{P_{i-1,j-1}(0)-\frac{x_{i-j}}{x_{i}}P_{i,j-1}(0)}{1-\frac{x_{i-j}}{x_{i}}} \\
    &=P_{i,j-1}(0)+\frac{P_{i,j-1}(0)-P_{i-1,j-1}(0)}{r-1}
  \end{aligned}
  \end{equation}$$
此时$r=(h_{i-j}/h_{i})^{p}$，形式和Richardson外推方法完全一致！

> **Neville算法不仅可以进行多项式递推，也可以进行有理函数形式的递推**
> $$$$
{:.prompt-tip}

<link rel="stylesheet" href="{{ '/assets/css/mystyle.css' | relative_url }}">