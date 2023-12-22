# MMTools
This is the shared package to simulate pulse propagation in a solid-core fiber with GMMNLSE/MM-UPPE with MATLAB.

It solves the pulse propagation with RK4IP if single mode and MPA if multimode. Both scalar and polarized scenarios can be simulated. Besides, it is implemented with an adaptive step-size control for both methods. Both passive and gain fibers can be simulated, where gain model includes Gaussian gain and rate-equation gain, for both single mode and multimode. Random mode coupling can be included as well.

For rate-equation-gain modeling, all pumping schemes are implemented: copumping, counterpumping, co+counter-pumping, as well as with and without ASE.

For multimode, GPU computations (with Nvidia cuda) is highly recommended. I have written a lot of cuda files to speed up simulations. It is controlled by `sim.gpu_yes=true or false`.

For details, please read the supplement of our paper: https://doi.org/10.1364/JOSAB.500586.  
Please don't forget to cite our paper if you find this code useful in your work. I, the young and early-career researcher, need your support. Similarly, if you need help or have questions about the code, please feel free to send me an email.

There is a readme.pdf in the Documentations/ folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics. This code is basically an upgraded and highly-optimized version of our "WiseLabAEP/GMMNLSE-Solver-FINAL" from "https://github.com/WiseLabAEP/GMMNLSE-Solver-FINAL," with much more functionalities, which however might overwhelm users and thus require more fiber-optic background. It can run order-of-magnitude faster than our old code due to optimizing with CUDA+shared memory, as well as reducing the usage of for-loops. Although our old one claims to be fast with GPU, its CUDA implementation is not optimized, let alone its CPU implementation with a lot of slow for-loops. Besides, this package includes adaptive step-size control, which improves the performance significantly and allows users to be free from worrying the reliability of a simulation. For optimization details, please see the supplement of our paper mentioned previously. 

If you have questions, feel free to ask them here or send me an email (email address is in my paper).

Note (12/22/2023):
I've been constantly modifying it, especially the ASE part, these days.
If you find a bug, download the latest one since I've been changing it often. If it still exists, let me know.
