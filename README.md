This is the source code for CSE8803 project.



fSVT is the folder from the paper
# Programs of fast SVT algorithm
---

##1.Main algorithms

rSVD/eigSVD.m ---- computes the singular value decomposition with eigendecomposition in [1]

rSVD/rSVDPI.m ---- fast randomized SVD with power iteration in [1]

rSVD/SVDBKI.m ---- fast randomized SVD with block Krylov-subspace iteration in [1]

SVT/fastSVT_Q.m ---- fast SVT with subspace Q reuse in [1]

SVT/fastSVT_U.m ---- fast SVT with subspace U reuse in [1]

##2.Other algorithms for comparson

rSVD/basicrSVD.m ---- Basic randomized SVD in [2]

rSVD/pcafast.m ---- randomized SVD in [3]

rSVD/rSVDpack.m ---- randomized SVD in [4]

rSVD/cSVD.m ---- Compressed SVD using Gaussian projection matrix in [5]

SVT/SVT.m ---- Singular value thresholding algorithm with 'svds' in Matlab in [6]

SVT/SVTlansvd.m ---- Singular value thresholding algorithm with 'lansvd' in PROPACK[7]

##3.Experiments of the testing

(1) Comparision of randomized SVD algorithms

rSVD/test.m is used for the experiment in comparing the cputime and relative error of basic rSVD, rsvdpack, pcafast, rsvdcs and rSVD-PI. The sequences of calculating the ralative error are commented because of the huge memory it needs. There are 3 nonzeros can be used in the test.

rSVD/test2.m is used for the experiment in comparing the rSVD-PI with rSVD-BKI with the svds a basic line.

rSVD/ml.mat is the matrix cut from 10M movielens dataset[7]. rSVD/sample.m is used to sample some point from the initial sparse matrix.

(2) The SVT for image recovery

SVT/new4.jpg is the image used in this section. SVT/pic_Omega.mat owns Ome (20% nonzeros) and Ome10 (10% nonzeros). 

SVT/testSVT_pic.m is used to test the SVT in image recovery. The comment can be changed to see the differences between SVT algorithm and fastSVT algorithm.

(3) The SVT for rating matrix completion

Firstly you should get the Movielens datasets [8] from there website https://grouplens.org/datasets/movielens/, then read the load data matrix with CSR format with variable 'Origin'(3 columns with CSR format), and the run the SVT/data_init.m to divide the dataset into 80% and 20% or 90% and 10%(controlled by the variable 'percent'). Then run the SVT/testSVT_ml.m. The delta of 20M movielen dataset is 4 compared with 5 of 10M dataset.

(4) The output in SVT/SVT.m and SVT/SVTlansvd.m every iteration is [#iters rank relative_error], and the output in SVT/fastSVT_Q.m and SVT/fastSVT_U.m every iteration is [#iters rank relative_error power_iteration_parameter].


# fSVT-py
fSVT-py is the python version of code translated from Matlab to provide better adaptibility for different enviroments

# fSVT-pictures
pictures under fSVT folders are results by running the code with different pictures