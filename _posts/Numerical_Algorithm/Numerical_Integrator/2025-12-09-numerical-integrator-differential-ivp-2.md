---
layout: post
title: 数值积分器：数值微分、定积分和IVP的关系（二）
author: myx
categories:
- 数值算法
- 数值积分器
math: true
date: 2025-12-09 16:46 +0800
---
## **简介**
<div class="header">
无论是<strong>数值微分</strong>、<strong>定积分问题</strong>还是<a href="{{ '/posts/numerical-integrator-intro/' | relative_url }}" class="internal-link">数值积分器（序）</a>中提到的IVP问题，解决问题的本质是利用离散的节点信息完成对于连续函数（未必存在封闭解析形式）的复现和近似，在数值分析的基础上完成一定阶数的近似拟合。
</div>

## **数值微分、定积分问题和IVP问题的关系**
### **定积分问题**
<div class="content-box">
在<a href="{{ '/posts/numerical-integrator-differential-ivp-1/' | relative_url }}" class="internal-link">数值积分器：数值微分、定积分和IVP的关系（一）</a>中，我们已经对数值分析中离散节点近似数值微分以及利用Richardson外推进行升阶有了比较详细的介绍和代码实现。其数值本质是对于没有简单解析形式或者仅有离散数据点的待估函数特定节点，依靠<strong>泰勒展开级数的线性组合</strong>消除误差主项，从而得到一定精度的微分近似量。<br><br>
</div>

在数值分析中，一定区间上的定积分求解问题往往也伴随着**解析积分不可行（无解析原函数、离散数据格式、复杂函数形式**的问题。区别于IVP问题中的一阶ODEs求解函数形式：
$$\begin{equation}
    y'=f(t,y)
    \label{Eq.1}
\end{equation}$$

在定积分问题中，一般存在微分方程并不会如\eqref{Eq.1}中依赖自身变量，而通常有简单形式：
$$\begin{equation}
    y'=f(t)\quad or \quad y'=f(x)
    \label{Eq.2}
\end{equation}$$
本质上，相较于IVP问题，定积分问题会更加简单，因为不涉及IVP问题外推（Propagation）时节点信息的隐式问题。举一个简单的例子，在利用梯形法则进行常微分方程积分过程中。对于IVP问题而言，离散积分格式为：

$$\Phi_{h}:y_{n+1}=y_{n}+h\frac{f(t_{n},y_{n})+f(t_{n+1},y_{n+1})}{2}$$

直观地可以看出$y_{n+1}$是隐式的，通常需要**迭代**步骤才能得到。而对于不依赖自身信息的定积分问题（假设积分区间为$\left[x_{0},x_{0}+h\right]$），积分格式为：

$$\Phi_{h}:\int^{x_{0}+h}_{x_{0}}y\,dx \approx y_{n+1}-y_{n}=h\frac{f(x_{n})+f(x_{n+1})}{2}$$

此时，积分格式的每个节点信息完全可以通过精确计算得到。而得益于此，积分可以脱离小区间约束，不需要类似IVP问题进行积分节点的实时更新。通常，进行区间定积分的数值方法主要有以下几种：

<div class="list-box">
<ul>
  <li>Newton-Cotes</li>
</ul>
</div>

<link rel="stylesheet" href="{{ '/assets/css/mystyle.css' | relative_url }}">