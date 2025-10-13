# rk4mrg
This repository contains the codes used in the paper 'Explicit integrators for nonlocal equations: The case of the Maxey-Riley-Gatignol equation' by Divya Jaganathan, Rama Govindarajan, and Vishal Vasan. See: https://doi.org/10.1090/qam/1693 (on arXiv: https://arxiv.org/abs/2308.09714v1)

Fortran codes (example1_rkxtd.f and example2_rkxtd.f) implement Eq. 4.1 for 2-stage (Section 4.1) and 4-stage (Section 4.2) in the two numerical experiments namely, the oscillating-in-time flow (Section 6.1) and particle in 2D stationary Lamb-Oseen vortex (Section 6.2). In addition, MATLAB codes (example3_main.m and auxiliary .m functions) used in the experiment with the alternative embedding procedure (Section 6.3) too are shared.

Figures 2, 3, 4 correspond to Examples 1, 2, 3 respectively.

For Example 1: To switch between 2-stage and 4-stage schemes, toggle the variable 'nstage' to 2 and 4 respectively. Control the time-step, dt, in inverse powers of 2, by changing the variable 'n2pow'. The variables dt and n2pow are related according to: dt = 1/2^n2pow. The rest of the parameters used in the paper are hard-coded. Running the Fortran code will generate a file "traj(nstage)f_h(n2pow).dat", comprising of two columns, namely time and solution.

For Example 2: Similar instruction as in Example 1. Except, the code generates a file comprising of five columns of information namely, time, x-component slip velocity, y-component slip velocity, x-component position, and y-component position.

For Example 3: Control the time-step, dt, using the variable 'npow'. The variable dt and npow are related via: dt = 1/2^npow. Number of even-indexed Hermite bases is controlled by the variable 'N'. 
