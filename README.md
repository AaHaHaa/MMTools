# MMTools
This is the shared package to simulate pulse propagation in a fiber with GMMNLSE/UPPE with MATLAB.

It solves the pulse propagation with RK4IP if single mode and MPA if multimode. Besides, it is implemented with an adaptive step-size control for both methods. Both passive and gain fibers can be simulated, where gain model includes Gaussian gain and rate-equation gain, for both single mode and multimode. Random mode coupling can be included as well.

For multimode, GPU computations (with Nvidia cuda) is highly recommended. I have written a lot of cuda files to speed up simulations. It is controlled by "sim.gpu_yes=true or false".

For details, please read the Supplement of our paper: (link here; wait to be published).

There is a readme.pdf in the Documentations/ folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics. This code is basically an upgraded and highly-optimized version of our "WiseLabAEP/GMMNLSE-Solver-FINAL" from "https://github.com/WiseLabAEP/GMMNLSE-Solver-FINAL.git," with much more functionalities, which however might overwhelm users and thus require more fiber-optic background. It can run order-of-magnitude faster than our old code due to optimizing with CUDA+shared memory, as well as reducing the usage of for-loops. Although our old one claims to be fast with GPU, its CUDA implementation is not optimized, let alone its CPU implementation with a lot of slow for-loops. Besides, this package includes adaptive step-size control, which improves the performance significantly and allow users to free from worrying the reliability of a simulation. For optimization details, please see the Supplement of our paper mentioned previously.
