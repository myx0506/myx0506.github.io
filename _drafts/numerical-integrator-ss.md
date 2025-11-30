---
layout: post
title: 数值积分器：数值微分、定积分和IVP问题
author: myx
categories:
- 数值算法
- 数值积分器
math: true
---
## **简介**
<div class="header">
在<a href="{{ '/posts/numerical-integrator-intro/' | relative_url }}" class="internal-link">数值积分器（序）</a>中，对于常微分方程组IVP问题我们常见采用的经典单步法积分器主要来自Runge-Kutta族，其最显著的数值特征就是可以用Butcher's Tableau系数表格进行数值架构的表征。对于RK方法，显式和隐式构造均存在，同样也可以利用辛积分条件构造保辛RK算法；此外，针对二阶ODEs，同样存在直接进行二阶量求解的Runge-Kutta-Nyström方法。
</div>




[^RK]: [Runge-Kutta 方法](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/457E8B0413C29B9D7BBCD7D2A45A23D5/S1446788700027932a.pdf/coefficients_for_the_study_of_rungekutta_integration_processes.pdf)

<link rel="stylesheet" href="{{ '/assets/css/mystyle.css' | relative_url }}">