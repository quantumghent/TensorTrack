TensorTrack
============

*An open-source tensor network library for MATLAB.*


Package summary
----------------

This is a package which aims to efficiently implement the various elementary tensor algorithms that arise in the context of tensor networks.
This includes various index manipulations, as well as contractions, factorizations, eigen- and linear solvers, ...
Additionally, for tensors which are invariant under general global symmetries, various mechanisms are in place to minimize both memory and cpu usage.

   
.. toctree::
   :caption: Manual
   :maxdepth: 2

   man/tensor
   man/symmetries
   man/algorithms


.. toctree::
   :caption: Tutorials
   :maxdepth: 1

   examples/examples.rst


.. toctree::
   :caption: Library
   :maxdepth: 2

   lib/tensors
   lib/sparse
   lib/mps
   lib/environments
   lib/algorithms
   lib/models
   lib/utility
   lib/caches
