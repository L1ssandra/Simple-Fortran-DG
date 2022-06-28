这是我在研0学习DG的过程中练习用的程序，比较简单，不过我觉得我现在写得已经很直观了，有需要的话可以自取

几点说明：
1. MATLAB中，contour命令是等高线，而mesh命令是画surface的。所以，如果要改成画等高线图，请把contour的注释取消，然后把mesh注释掉；另外把axis(s)这一句也要注释掉，它是一个三维的画布。如果要看数值解的其他分量，请把Q1 = str2num(fileread('Q2.txt'))这一句中后面的“Q2”改成其他的。
2. Fortran中，DG_2D_New是主程序，RK3DG是子程序，现在默认的是用80×80的网格解Euler方程，初值为等熵涡流。
3. 关于程序的详细说明可参见知乎笔记：https://zhuanlan.zhihu.com/p/533958396

Reference：

[1] Chi-Wang Shu, Discontinuous Galerkin Methods: General Approach and Stability.

[2] Bernardo Cockburn, Chi-Wang Shu, The Runge–Kutta Discontinuous Galerkin Method for Conservation Laws V: Multidimensional Systems.
