# 2D-TOpt
Research tool developed for conducting topology optimization in 2D. The main purpose of this tool is to generate a image dataset for deep learning purposes.
Additional features includes design space variation using active and passive elements, various boundary and loading conditions. These tool contains four matlab files:

1. initial_model.m :- This file contains a function responsible for generating volume fraction, boundary and loading condition images.
2. TO_main.m:- A object designed to carry out SIMP algorithm based on initial conditions provided by TO.m
3. TO.m:- A object designed to generate, run and save different cases at once.
4. run.m:- RUN the cases generated by TO.m.

Use run.m to generate and run cases. Instructions are embedded in rum.m as comments.