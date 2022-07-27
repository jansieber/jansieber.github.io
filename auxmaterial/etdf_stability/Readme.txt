# Supplementary material for Generic stabilisability for time-delayed feedback
control

Author: Jan Sieber

This folder contains Matlab code and its output for illustrating examples used
in

Jan Sieber: _Generic stabilisability for time-delayed feedback control_,
Preprint download: [arxiv:1508.05671](http://arxiv.org/abs/1508.05671).

List of Matlab scripts and their outputs:

  * [Hopf example output](html/HopfExample.html) (html): published Matlab output of [HopfExample.m](HopfExample.m) script. This shows the Matlab code with formatted comments and graphs similar to the example section in the paper.
  * [Output](html/ParametricPendulum.html) (html) for demo [ParametricPendulum.m](ParametricPendulum.m) using extended time-delayed feedback (ETDF) control to stabilise unstable oscillations in a **parametrically excited pendulum**. 
List of Matlab functions used in the demonstration scripts:

  * [SpecExpAssign.m](SpecExpAssign.m): computes gains K to assign spectrum of matrix P exp(bK) for nxn matrix P and nx1 vector b. 
  * [etdf_controlled_rhs.m](etdf_controlled_rhs.m): wrapper around right-hand side f (in x=f(x,p,u)) to make it fit into DDE-Biftool's format
  * [etdf_control.m](etdf_control.m): constructs ETDF gains for asymptotic parameters and compute true stability
  * [gain_of_x.m](gain_of_x.m): convert gain as contructed by Brunovsky to true x-dependent gain K(x)
Note that the Matlab functions and scripts rely on DDE-Biftool 3.1, available
here: [sourceforge.net/projects/ddebiftool](https://sourceforge.net/projects/d
debiftool).

