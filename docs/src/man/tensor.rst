Tensor manipulations
====================

Definitions
-----------

We start with the definition of a tensor. According to the `Link Wikipedia page <https://en.wikipedia.org/wiki/Tensor#Definition>` on tensors, there are 3 equivalent ways of describing a tensor object:

* As a multidimensional array.

* As a multilinear map.

* Using tensor products.

From a MATLAB perspective, the former approach is the most natural, as such tensors are easily created with several builtin methods. We thus say a tensor of size ``sz`` can be thought of as the array that is created by calling any of the standard functions ``t = zeros(sz)``, ``t = ones(sz)``, ``t = rand(sz)``, ...
We will use a graphical representation of such a tensor as an object with several legs, where each leg represents a tensor index. By convention we will order the indices from top-left to top-right, going around the object counter-clockwise.

Nevertheless, for our purpose it will often be useful to be able to think of such an abstract object as a map, or thus as a matrix.
For this, we additionally specify the rank of the tensor ``[n1 n2]``, which results in a partition of the indices into ``n1`` left and ``n2`` right indices, which we will call **codomain** and **domain** respectively.
If we then reshape the original array into a matrix. :code:`m = reshape(t, prod(sz(1:n1)), prod(sz((1:n2)+n1)))`.
With this definition however, we still run into issues when we want to perform something like matrix multiplication.
As MATLAB uses column-major ordering of array elements, ``reshape`` effectively groups indices together from left to right, where the left index is changing the fastest.
If we then multiply two tensors, we connect the domain (right indices) of the first tensor with the codomain (left indices) of the second, and note that the counter-clockwise ordering of our indices now causes the domain indices to be grouped bottom-to-top, while the codomain indices are grouped top-to-bottom.
Thus, in order for our matrix multiplication to be consistent, we reverse the order of the domain indices, and define the tensor map as :code:`m = reshape(permute(t, [1:n1 n1+flip(1:n2)]), prod(sz(1:n1)), prod(sz(n1+(1:n2))`.

The mathematical origin of this unfortunate permutation is found when considering the dual of a tensor product of spaces, which is isomorphic to the tensor product of the dual spaces, but reversed.
In other words, we start out with a tensor as follows:

.. math::
    t \in V_1 \otimes V_2 \otimes \dots \otimes V_n.

Then we can reinterpret this as a map, by partitioning some of the indices to the domain:

.. math::
    t^\prime : V_1 \otimes V_2 \otimes \dots \otimes V_{n_1} \leftarrow V_{n_2}^* \otimes V_{n_2-1}^* \otimes \dots \otimes V_{n_1+1}^*.

However, in doing so, we had to reverse their order. 

In summary, we consider a tensor to be an object which can be represented in either one of the previously discussed forms, i.e. as a multidimensional array, or as a matrix.
These definitions are equivalent, and it is always possible to go from one to the other.
The main reason for introducing both is that for some algorithms, one is more convenient than the other, and thus this conversion will prove useful.

Creating tensors
----------------

There are several options provided for creating tensors, which are designed to either mimic the existing MATLAB constructors, or to be conveniently generalisable when dealing with symmetries. Almost each of these methods is a wrapper of some sorts around the general static constructor ``Tensor.new``:

.. automethod:: src.tensors.Tensor.new
    :noindex:

In summary, for tensors without any internal symmetry, the easiest way of constructing tensors is by using the form :code:`Tensor.new(fun, [dims...])`, where optionally it is possible to specify arrows for each of the dimensions with :code:`Tensor.new(fun, [dims...], [arrows...])`.
Additionally, it is possible to make the partition of domain and codomain by supplying a rank for the tensor :code:`Tensor.new(fun, [dims...], [arrows...], 'Rank', [r1 r2])`.
While the latter two options do not alter any of the results in the non symmetric case, adding arrows is often useful from a debugging point of view, as it will then become less likely to accidentaly connect tensors that are not compatible. Additionally, the rank determines how the tensor acts as a linear map, which becomes important for decompositions, eigenvalues, etc.
Finally, the most general form is found by explicitly constructing the vector spaces upon which the tensor will act.
For tensors without symmetries, :class:`ComplexSpace` and :class:`CartesianSpace` represent spaces with or without arrows respectively, while :class:`GradedSpace` is used for tensors which have internal symmetries.
These allow for the final form :code:`Tensor.new(fun, codomain, domain)`.
On top of this, several convenience methods exist for commonly used initializers, such as:

* :code:`Tensor.zeros`
* :code:`Tensor.ones`
* :code:`Tensor.eye`
* :code:`Tensor.rand`
* :code:`Tensor.randn`
* :code:`Tensor.randc`
* :code:`Tensor.randnc`

Often, it can prove useful to create new tensors based on existing ones, either to ensure that they can be added, subtracted, or contracted.
In these cases, :code:`Tensor.similar` will be the method of choice, which acts as a wrapper to avoid manually having to select spaces from the tensors.

.. automethod:: src.tensors.Tensor.similar
    :noindex:


Accessing tensor data
---------------------

In general, :class:`Tensor` objects do not allow indexing or slicing operations, as this is not compatible with their internal structure.
Nevertheless, in accordance with the two ways in which a tensor can be represented, it is possible to access the tensor data in two ways.

First off, when we consider a tensor as a linear operator which maps the domain to the codomain, we can ask for the matrix representation of this map.
For tensors with internal symmetries this will in general be a list of matrices, that represent the block-diagonal part of the tensor, along with the corresponding charge through :code:`matrixblocks(t)`.

Secondly, when we want to think of a tensor more in terms of an N-dimensional array, we can access these as :code:`tensorblocks(t)`.
For tensors with internal symmetries, this will generate a list of all channels that are not explicitly zero by virtue of the symmetry, along with a representation of these channels, which is called a :class:`FusionTree`.

Additional details for the symmetric tensors can be found in the :ref:`Symmetries` section.

As an example, it could prove educational to understand the sizes of the lists and the sizes of the blocks generated by the example code below:

.. code-block:: matlab
    
    >> A = Tensor.zeros([2 2 2], 'Rank', [2 1]);
    >> Ablocks = matrixblocks(A)

    Ablocks =

    1×1 cell array

        {4×2 double}

    >> Atensors = tensorblocks(A)

    Atensors =

    1×1 cell array

        {2×2×2 double}
    
    >> z2space = GradedSpace.new(Z2(0, 1), [1 1], false);
    >> B = Tensor.zeros([z2space z2space], z2space);
    >> [Bblocks, Bcharges] = matrixblocks(B)

    Bblocks =

    1×2 cell array

        {2×1 double}    {2×1 double}


    Bcharges = 

    1×2 Z2:

    logical data:
    0   1

    >> [Btensors, Btrees] = tensorblocks(B)

    Btensors =

    4×1 cell array

        {[0]}
        {[0]}
        {[0]}
        {[0]}


    Btrees = 

    (2, 1) Z2 FusionTree array:

        isdual:
        0  0  0
        charges:
        + + | + | +
        - - | + | +
        - + | - | -
        + - | - | -

In the very same way, in order to write data into a tensor, the same two formats can be used.

First off, :code:`t = fill_matrix(t, blocks)` will take a list of blocks and fill these into the tensor.
This requires the list to be full, and sorted according to the charge, or in other words it has to be of the same shape and order as the output of :code:`matrixblocks`.
If it is only necessary to change some of the blocks, :code:`t = fill_matrix(t, blocks, charges)` additionally passes on an array of charges which specifies which block will be filled.

Similarly, :code:`t = fill_tensor(t, tensors)` will take a list of N-dimensional arrays and fill these into the tensor, in the same order and shape of the output of :code:`tensorblocks`.
If it is required to only change some of the tensors, an array of :class:`FusionTree` s can be passed in as :code:`t = fill_tensor(t, tensors, trees)` to specify which tensors should be changed.

.. code-block:: matlab
    
    >> Ablocks{1} = ones(size(Ablocks{1}));
    >> A = fill_matrix(A, Ablocks);
    >> Atensors{1} = rand(size(Atensors{1}));
    >> A = fill_matrix(A, Atensors);
    
    >> Bblocks = cellfun(@(x) ones(size(x)), Bblocks, 'UniformOutput', false);
    >> B = fill_matrix(B, Bblocks);
    >> Bblocks{1} = Bblocks{1} + 1;
    >> B = fill_matrix(B, Bblocks(1), Bcharges(1));

Additionally, it is also possible to use initializers instead of a list of data.
These initializers should have signature :code:`fun(dims, identifier)`.
For non symmetric tensors, ``identifier`` will be omitted, but for symmetric tensors the matrix case uses charges as ``identifier``, while the tensor case uses fusion trees as ``identifier``.
Again, it is possible to select only some of the blocks through the third argument.

.. code-block:: matlab
    
    >> f = @(dims, identifier) ones(dims);
    >> A = fill_matrix(A, f);
    >> A = fill_tensor(A, f);
    
    >> g = @(dims, charge) qdim(charge) * ones(dims);
    >> B = fill_matrix(B, g);
    >> h = @(dims, tree) qdim(f.coupled) * ones(dims);
    >> B = fill_tensor(B, h);

Finally, we mention that for most tensors, it is possible to generate an N-dimensional array representation, at the cost of losing all information about the symmetries.
This can sometimes be useful as a tool for debugging, and can be accessed through :code:`a = double(t)`.


Index manipulations
-------------------

Once a tensor has been created, it is possible to manipulate the order and partition of the indices through the use of :code:`permute(t, p, r)`.
This methods works similarly to :code:`permute` for arrays, as it requires a permutation vector ``p`` for determining the new order of the tensor legs.
Additionally and optionally, one may specify a rank `r` to determine the partition of the resulting tensor.
In order to only change the partition without permuting indices, :code:`repartition(t, r)` also is implemented.

.. code-block:: matlab
    
    >> A = Tensor.rand([1 2 3])

    A = 

    Rank (3, 0) Tensor:

    1.	CartesianSpace of dimension 1

    2.	CartesianSpace of dimension 2

    3.	CartesianSpace of dimension 3

    >> A2 = repartition(A, [1 2])

    A2 = 

    Rank (1, 2) Tensor:

    1.	CartesianSpace of dimension 1

    2.	CartesianSpace of dimension 2

    3.	CartesianSpace of dimension 3
    
    >> A3 = permute(A, [3 2 1])

    A3 = 

    Rank (3, 0) Tensor:

    1.	CartesianSpace of dimension 3

    2.	CartesianSpace of dimension 2

    3.	CartesianSpace of dimension 1

    >> A3 = permute(A, [3 2 1], [2 1])

    A3 = 

    Rank (2, 1) Tensor:

    1.	CartesianSpace of dimension 3

    2.	CartesianSpace of dimension 2

    3.	CartesianSpace of dimension 1


.. note:: 
    
    While the partition of tensor indices might seem of little importance for tensors without internal structure, it can still have non-trivial consequences.
    This is demonstrated by comparing the ``matrixblocks`` and the ``tensorblocks`` before and after repartitioning.


Contractions 
------------

It is also possible to combine multiple tensors by contracting one or more indices.
This is only possible if the contracted spaces are compatible, i.e. one is the dual space of the other.
In general, all contractions will be performed pairwise, such that contracting ``A`` and ``B`` consists of permuting all uncontracted indices of ``A`` to its codomain, all contracted indices of ``A`` to its domain, and the reverse for ``B``.
Then, contraction is just a composition of linear maps, hence the need for the contracted spaces to be compatible.

The syntax for specifying tensor contractions is based on the ``NCon`` (network contraction) algorithm described here (https://arxiv.org/abs/1402.0939).
The core principle is that contracted indices are indicated by incrementing positive integers, which are then pairwise contracted in ascending order.
Uncontracted indices are specified with negative integers, which are sorted in descending order (ascending absolute value).
It is also possible to specify the rank of the resulting tensor with a name-value argument ``'Rank'``, and use in-place conjugation with the name-value argument ``'Conj'``.

.. autofunction:: src.utility.linalg.contract
    :noindex:


Factorizations 
--------------

The reverse process, of splitting a single tensor into different, usually smaller, tensors with specific properties is done by means of factorizations.
In short, these are all analogies of their matrix counterpart, by performing the factorization on the tensor interpreted as a linear map.

Many of these factorizations use the notion of orthogonality (unitarity when working over complex numbers).
An orthogonal tensor map ``t`` is characterised by the fact that ``t' * t = eye = t * t'``, which can be further relaxed to left- or right- orthogonal by respectively dropping the right- and left- hand side of the previous equation.

.. automethod:: src.tensors.Tensor.leftorth
    :noindex:

.. automethod:: src.tensors.Tensor.rightorth
    :noindex:
    
.. automethod:: src.tensors.Tensor.tsvd
    :noindex:

.. automethod:: src.tensors.Tensor.leftnull
    :noindex:

.. automethod:: src.tensors.Tensor.rightnull
    :noindex:

.. automethod:: src.tensors.Tensor.eig
    :noindex:


Matrix functions
----------------

Finally, several functions that are defined for matrices have a generalisation to tensors, again by interpreting them as linear maps.

.. automethod:: src.tensors.Tensor.expm
    :noindex:

.. automethod:: src.tensors.Tensor.sqrtm
    :noindex:

.. automethod:: src.tensors.Tensor.inv
    :noindex:

.. automethod:: src.tensors.Tensor.mpower
    :noindex:
