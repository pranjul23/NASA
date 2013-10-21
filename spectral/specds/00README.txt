SpecDS -- simple tools for spectral learning of dynamical systems.

The most important function is specds.m, for learning a dynamical
system from data.

An example of usage is in hmmex.m -- it builds a simple HMM, samples
from it, and learns it back from the sampled data.

See the file 01GPL.txt for licensing information.

This is a very simple implementation; most of the complexity comes
from the way we compute covariance matrices in specds.m, which is more
complex than it should be to keep Matlab from running out of memory by
making many copies of the observations.  A bit more complexity could
probably make this computation a lot faster as well, but that's not
implemented yet.

The best documentation of this algorithm is probably the ICML tutorial
slides, available at:

http://www.cs.cmu.edu/~ggordon/spectral-learning/

The following paper contains more information, but some details are
slightly out of date:

Byron Boots, Sajid Siddiqi, and Geoff Gordon.  Closing the
Learning-Planning Loop with Predictive State Representations.
International Journal of Robotics Research (IJRR), 30(7):954-966,
2011.

