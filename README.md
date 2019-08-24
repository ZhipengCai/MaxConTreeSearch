# Demo---Consensus-maximization-tree-search-revisited (under construction)

A significantly accelerated tree search method for globally optimal consensus maximization. [(Paper link)](https://arxiv.org/abs/1908.02021)

Published in ICCV 2019 as **oral** presentations.

About
=====

![alt text](https://github.com/ZhipengCai/ZhipengCai.github.io/blob/master/papers/ICCV19.png " ")

Consensus maximization is an effective tool for robust fitting in computer vision. A* Tree Search is one of the most efficient method for globally optimal consensus maximization. In this work, we propose two new techniques that significantly accelerate A* Tree Search, making it capable of handling problems with much larger number of outliers.

The demo is free for non-commercial academic use. Any commercial use is strictly 
prohibited without the authors' consent. Please acknowledge the authors by citing:

```
@article{cai2019consensus,
  title={Consensus Maximization Tree Search Revisited},
  author={Cai, Zhipeng and Chin, Tat-Jun and Koltun, Vladlen},
  journal={arXiv preprint arXiv:1908.02021},
  year={2019}
}
```
in any academic publications that have made use of this package or part of it.

------------------------
Contact
------------------------

Homepage:https://zhipengcai.github.io/

Email: czptc2h@gmail.com

Do not hesitate to contact the authors if you have any question or find any bugs :)


Getting Started
===============

The demo is implemented using MATLAB 2018b and has been tested on Ubuntu 14.04 LTS 64-bit. 

-------------
Run the demo
-------------

1. Clone this repository. 

2. Run "demo.m" in MATLAB.

Please refer to "demo.m" file for detailed code explanations.

-----------------------------------------
List of addressed problems in the demo
-----------------------------------------

Linear problem:

1. Linearized Fundamental matrix estimation (ignoring the rank-2 constraint)

Nonlinear problem:

2. Homography estimation

------------------------
List of included methods
------------------------

Previous A* tree search variants:

1. [A* tree search (Chin et al. CVPR'15)](https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Chin_Efficient_Globally_Optimal_2015_CVPR_paper.pdf)

2. [A* tree search + True Outlier Detection (TOD) for branch pruning (Chin et al. TPAMI'17)](https://ieeexplore.ieee.org/document/7755788)

Variants with our new techniques:

3. A* tree search + Non-Adjacent Path Avoidance (NAPA)

4. A* tree search + NAPA + TOD

5. A* tree search + NAPA + Dimension-Insensitive Branch Pruning (DIBP)

