Tensors
=======

.. toctree::
   :maxdepth: 2
   
.. module:: src
   
This section contains the API documentation for the :mod:`.tensors` module.
   

Symmetry sectors
----------------

Type hierarchy
``````````````

.. autoclass::  src.tensors.charges.AbstractCharge
.. autoclass::  src.tensors.charges.FusionStyle
.. autoclass::  src.tensors.charges.BraidingStyle
.. autoclass::  src.tensors.charges.ProductCharge


Concrete charge types
`````````````````````

.. autoclass:: src.tensors.charges.Z1
.. autoclass:: src.tensors.charges.Z2
.. autoclass:: src.tensors.charges.ZN
.. autoclass:: src.tensors.charges.fZ2
.. autoclass:: src.tensors.charges.U1
.. autoclass:: src.tensors.charges.fU1
.. autoclass:: src.tensors.charges.SU2
.. autoclass:: src.tensors.charges.fSU2
.. autoclass:: src.tensors.charges.SUN
.. autoclass:: src.tensors.charges.O2
.. autoclass:: src.tensors.charges.A4


Helper routines
```````````````

.. autoclass:: src.tensors.charges.GtPattern


Fusion trees
------------

.. autoclass:: src.tensors.FusionTree


Spaces
------

Type hierarchy
``````````````

.. autoclass::  src.tensors.spaces.AbstractSpace
.. autoclass::  src.tensors.spaces.SumSpace


Concrete space types
````````````````````

.. autoclass::  src.tensors.spaces.CartesianSpace
.. autoclass::  src.tensors.spaces.ComplexSpace
.. autoclass::  src.tensors.spaces.GradedSpace


Convenience constructor wrappers
````````````````````````````````

.. autofunction::  src.tensors.spaces.Z2Space
.. autofunction::  src.tensors.spaces.fZ2Space
.. autofunction::  src.tensors.spaces.U1Space
.. autofunction::  src.tensors.spaces.SU2Space

Helper classes
``````````````

.. autoclass::  src.tensors.spaces.Arrow


Kernels
-------

.. autoclass:: src.tensors.kernels.AbstractBlock
.. autoclass:: src.tensors.kernels.TrivialBlock
   :no-members:
.. autoclass:: src.tensors.kernels.MatrixBlock
   :no-members:
.. autoclass:: src.tensors.kernels.AbelianBlock
   :no-members:


Tensors
-------

.. autoclass:: src.tensors.AbstractTensor
.. autoclass:: src.tensors.Tensor
