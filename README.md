# Numerical Implementation of Differential Equations

As part of my undergraduate studies at the University of Basel, I undertook a project (1 ECTS) in 2019 with the objective of acquiring proficiency in MATLAB-based implementation of differential equations. The sole purpose of this repository is to serve as an archive for the project. It is important to mention that the project is conducted in the German language.

The project is divided into three distinct sections: Part I, Part II, and Part III, with each section having a corresponding "main"-MATLAB file. The source is found in the folder [matlab-source](https://github.com/sabrimeyer/numerical-de/tree/main/matlab-source). Each part studies the behaviour or four different numerical methods for solving ordinanry differential equations; the [explicit Euler method](https://en.wikipedia.org/wiki/Euler_method), [Heun's method](https://en.wikipedia.org/wiki/Heun%27s_method), and two [Runge-Kutta methods](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods).



**Part I:** First, we compare the four aforementioned methods for solving a differential equation of the form
$$\frac{dy}{dt}=f(t,y(t)),\quad y(t_0)=y_0.$$
We are particularly interested in the convergence rates.

**Part II:** Here we study

The file [presentation.pdf](https://github.com/sabrimeyer/numerical-de/blob/main/documentation/presentation.pdf) was used to present and discuss the results of the project. Complementary to the presentation, the file [documentation.pdf](https://github.com/sabrimeyer/numerical-de/blob/main/documentation/documentation.pdf) written in LaTeX is included.
