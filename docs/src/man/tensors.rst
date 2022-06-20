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
If we then reshape the original array into a matrix. ``m = reshape(t, prod(sz(1:n1)), prod(sz((1:n2)+n1)))``.
With this definition however, we still run into issues when we want to perform something like matrix multiplication.
As MATLAB uses column-major ordering of array elements, ``reshape`` effectively groups indices together from left to right, where the left index is changing the fastest.
If we then multiply two tensors, we connect the domain (right indices) of the first tensor with the codomain (left indices) of the second, and note that the counter-clockwise ordering of our indices now causes the domain indices to be grouped bottom-to-top, while the codomain indices are grouped top-to-bottom.
Thus, in order for our matrix multiplication to be consistent, we reverse the order of the domain indices, and define the tensor map as ``m = reshape(permute(t, [1:n1 n1+flip(1:n2)]), prod(sz(1:n1)), prod(sz(n1+(1:n2))``.

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


Index manipulations
-------------------


Contractions 
------------


Factorizations 
--------------
