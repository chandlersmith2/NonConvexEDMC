# NonConvexEDMC
This is a repository associated with the paper "Provable Non-Convex Euclidean Distance Matrix Completion: A Study of Convergence, Robustness, and Incoherence" seen at https://arxiv.org/abs/2508.00091


## Paper Intro

This paper seeks to solve the Euclidean Distance Matrix Completion problem. In this paper, we establish convergence guarantees to a ground truth Gram matrix given partially observed distances for a first-order method on the manifold of rank-$r$ matrices. This methodology uses a dual basis approach, i.e. $D_{{ij}} = \left< X,\mathbf{w}_\alpha \right>$ where $X$ is a Gram matrix, and $D_{ij} = \lVert p_{i}-p_{j} \rVert_{2}^{2}$.

## Algorithm

The algorithm of note is distgeo_M_omega.m, with an example of how it runs seen in synthetic_data_test_script.m.
