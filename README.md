# Rayleigh Benard Flow Simulation
Rayleigh Benard convection is created when a fluid is heated by a hot plate from below and cooled from a cold plate from above. The hot liquid rises as the cold liquid falls, creating beautiful vortices.

Serial lattice boltzmann implementations in Python (Jupyter Notebooks) and C++, GPU accelerated version in CUDA. Based on Prof. Jonas Latt's Matlab implementation and CAPT Stu Blair's LBM code.

Results: In a trial with ~5000 lattice points, the C++ implementation performed 2.97*10^6 lattice point updates per second whereas the CUDA implementation performed 1.15*10^8 lattice point updates per second. Since more lattice points creates more parallelism, the CUDA implementation has greater lattice point updates per second with larger lattices. In a 2 million lattice point problem the CUDA implementation performed just under 1*10^9 lattice point updates per second.
