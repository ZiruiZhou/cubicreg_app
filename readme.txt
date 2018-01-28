This zip-file contains the MATLAB codes for the numerical experiments of the paper:

"On the Quadratic Convergence of the Cubic Regularization Method under a Local Error Bound Condition"

by Man Chung Yue, Zirui Zhou, and Anthony Man-cho So. 

For any enquries, please contact "zrzhou01 at gmail.com".



************ User Manual *************


To generate the figures in that paper, run the following commands in MATLAB.

	
1. For phase retrieval, run: "PR_Test(n,sav)";
		
	n: positive integer, the dimension of the target signal;
	sav: 1 or 0. 1 stands for save the data, 0 stands for not save. Default: 0.

	
2. For low-rank matrix recovery: "LRM_Test(n,r,sav):;

	[n,r]: the dimension of the true matrix, n*n with rank=r.
 	sav: 1 or 0. 1 stands for save the data, 0 stands for not save. Default: 0.


