!*************************************************************************
!*                                                                       *
!*     A-PWLSVR - Adaptive Piecewise Linear Support Vector Regression    *
!*                                                                       *
!*     by Napsu Karmitsa 2019, Adil Bagirov and Kaisa Joki               *
!*     (last modified 10.01.2020).                                       *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but we ask you to refer the appropriate references       *
!*     given below, if you use it.                                       *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*
!*     spr.f03               - Mainprogram (this file).
!*     constants.f03         - Parameters and constants.
!*     initspr.f03           - Initialization of parameters for A-PWLSVR and DBDC.
!*     functions.f03         - Computation of DC components f_1 and f_2 and their subgradients for the PWLSVR problem.
!*     bundle1.f03           - Bundle of DC component f_1.
!*     bundle2.f03           - Bundle of DC component f_2.
!*     dbdc.f03              - DBDC method.
!*     plqdf1.f              - Quadratic solver by L. Luksan.
!*
!*     Makefile              - makefile.
!*
!*
!*     Calling sequence:
!*
!*       ./spr infile nrecord ntrain nft maxmin [noutcom maxlin eps]
!*
!*     with the following arguments:
!*
!*     infile                - the data set;
!*     nrecord               - number of records in the data set;
!*     ntrain                - size of the training set, ntrain <= nrecord;
!*     nft                   - number of features in the data set;
!*     maxmin                - maximum number of min functions under maximum in the model;
!*     noutcom               - number of output column, optional, nft by default;
!*     maxlin                - number of linear functions under each minimum in the model, optional,
!*                             Default value = 3;
!*     eps                   - stopping tolerance, optional, 10^(-2) by default.
!*
!*     If you want to tune other parameters, modify initspr.f03. 
!*     As default the results are written in
!*        
!*     predictions.txt       - Evaluation criteria (e.g. RMSE, MAE,...) with training and predictions;
!*     reg_results.txt       - Regression results with hyperplanes.
!*     pwl_points.txt        - Points of original data and pwl-model obtained
!*                             for drawing purposis. Works only with nft = 2 and maxlin = 3 (otherwise empty). 
!*
!*
!*     References:
!*       A. Bagirov, S. Taheri, N. Karmitsa, K. Joki, M.M. Mäkelä, 
!*       "Adaptive piecewise linear support vector regression", submitted 2020.
!*
!*
!*     Acknowledgements: 
!*       The research work by: 
!*       Prof. Adil Bagirov and Dr. Sona Taheri was supported by the Australian Government through the Australian Research Council's      !*       Discovery Projects funding scheme (Project No. DP190100580); and  
!*       Dr. Napsu Karmitsa and and Dr. Kaisa Joki was supported by the Academy of Finland (Project No. 289500 and 319274).
