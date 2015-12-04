# DASH_Optimal_Set
ILP problem to optimize the representations set in DASH systems

Files Contents 

settings_model.cpp

This code generates a testing scenario file for the ILP model in model.cpp using the default settings in Section 6 of [1], that is, an instance of values for the input parameters of the ILP model described in [2]. This file is saved in a folder named 'scenarios/' at the same level in the directories tree as 'settings_model.cpp'. The objective of this code is merely to illustrate how to assign values to the input parameters of model.cpp to build an instance of the optimization problem.

NOTE: This code assumes that you will use the scenario file to feed model.cpp code. model.cpp solves the ILP formulation in [2] by using IBM ILOG CPLEX Optimizer [3]. This CPLEX solver makes use of the ILOG Concert Technology, that offers a C++ library of classes and functions enabling to define models for optimization problems and to apply algorithms to those models. Thus, you need to install IBM ILOG CPLEX and to prepare a makefile according to the IBM ILOG CPLEX indications to be able to compile and run this code.


model.cpp

This code computes an optimal set of video representations in an Adaptive Streaming scenario using the ILP model described in [2]. The model assumes that a CDN entity prepare multiple representations of the videos in its catalogue. Each representation is characterized by only 3 parameters: the own video, the encoding rate and the spatial resolution. The set is optimal in terms of maximizing the sum of individual satisfactions of a given population of users, where users are defined by their video requests, bandwidth of internet connections and display sizes. Finally, the CDN entity is constrained by a limited CDN budget to buy bandwidth to ISPs and a limited number of representation that it can prepare. The ILP model is solved using IBM ILOG CPLEX solver. More details about the ILP model are in [2].

The 'model.cpp'code loads from a folder named 'scenarios/' (at the same level in the directories tree as 'model.cpp') a file containing all the input parameters characterizing a given scenario: video catalogue, encoding settings, satisfaction model, users population, .... These parameters are the input of the model. After finding a solution to the optimization problem, 'mode.cpp' code saves the output of the optimization process in the folders named 'optimalValues/', 'status/' and 'representations/'. In 'optimalValues/, it is saved a file that' contains the optimal valued found by the CPLEX solver for the decision variables of the model. In 'status/', it is saved a file that contains the status of solution found by CPLEX solver, for example, if the solution is optimal, the status is the string'optimal'. Finally, in 'representations/', it is saved a file that contains a list of the optimal representations found by the solver. Each row in the file corresponds to a representation, where first, second and third column corresponds to video index, encoding rate (in kbps) and resolution (e.g. 360p or 1080p).

NOTE: This code uses IBM ILOG CPLEX optimizer [3] to solve the ILP formulation in [2]. Thus, you need to install IBM ILOG CPLEX and to prepare a makefile according to the IBM ILOG CPLEX indications to be able to compile and run this code.


References

[1] Laura Toni, Ramon Aparicio-Pardo, Gwendal Simon, Alberto Blanc, Pascal Frossard, "Optimal Set of Video Representations in Adaptive Streaming," in Proc. ACM Multimedia Systems conference, MMSys 2014, Singapore, March 2014

[2] Laura Toni, Ramon Aparicio-Pardo, Karine Pires, Alberto Blanc, Gwendal Simon, Pascal Frossard, “Optimal Selection of Adaptive Streaming Representations,” ACM Transactions on Multimedia Computing, Communications and Applications (ACM TOMM), vol. 11, no. 2s, article 43, February 2015

[3] IBM ILOG CPLEX Optimizer http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
