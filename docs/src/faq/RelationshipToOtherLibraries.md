# Relationship of ITensor to other tensor libraries

Here we will describe the relationship of ITensor to more traditional Julia Arrays or deep learning libraries like TensorFlow and PyTorch. There are a few things that distinguish ITensor from those approaches:

1. ITensors have dimensions with labels that get passed around, which makes it simple to perform certain operations like contraction, addition, and tensor decompositions with a high level interface, independent of memory layout. This is along the same lines as Julia packages like [NamedDims.jl](https://github.com/invenia/NamedDims.jl) and [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl) and libraries in Python like [xarray](https://xarray.pydata.org/en/stable/index.html), however I would argue that the ITensor approach is a little more sophisticated (the dimensions have more metadata which makes them easier to manipulate for different situations, random ids to help avoid name clashes, etc.). This design was inspired by the needs of tensor network algorithms, where there are many tensor dimensions in the computation (of which many of them are dynamically created during the calculation), but would be helpful for writing other algorithms too.

2. The ITensor type has a dynamic high level interface, where the type itself is mutable and the data can be swapped out. This allows for conveniently allocating the data of an ITensor on the fly "as needed", which makes for a nicer, more flexible interface (like initializing an empty ITensor before a loop, and filling it with the correct data type when the first value is set), at the expense of a small overhead for accessing data in the ITensor. We have found this tradeoff is worth it, since we expect ITensors to be used for medium to large scale calculations where operations on the tensors like contraction, addition, and tensor decomposition dominate the cost of the calculation, and code can be designed with function barriers to speed up operations when data is being accessed repeatedly.

3. Another feature that ITensor has that goes beyond what is available in standard Julia, TensorFlow, and PyTorch is tensors which are [symmetric under a group action](https://arxiv.org/pdf/1008.4774.pdf). The physical interpretation of these tensors are ones that have a conserved quantity (like a quantum state with a conserved number of particles), so that feature is more physics-oriented, but could have applications in other areas like machine learning as well. In practice, these tensors are block sparse, and have extra metadata on the dimensions labeling representations of the group.

4. Based on the features above, the ITensor library provides high level implementations of tensor network algorithms (algebraic operations of very high dimensional tensors, such as addition, multiplication, and finding dominant eigenvectors). In general these algorithms can (and have been) written on top of other libraries like standard Julia Arrays/AD, PyTorch, or TensorFlow, but they might have various downsides (a less convenient interface for dealing with tensor operations, no support for the types of symmetric tensors we often need, limited support for tensors with complex numbers in the case of libraries like PyTorch, though perhaps that has improved since I last checked, etc.).

Although ITensor has primarily focused on quantum physics and quantum computing applications, there is work using ITensor for machine learning applications (so far focused on applications of tensor networks to machine learning, so no neural network calculations yet as far as I know). In general, these different libraries (ITensor, Flux, PyTorch, TensorFlow) are biased towards their specific methods and application areas that they are used for the most: ITensor is more biased towards tensor network calculations and quantum physics/quantum computing applications, based on the available features and interface, while PyTorch and TensorFlow are more biased towards neural network calculations. However, our goal would be to provide more features to ITensor that would make it useful for neural network applications as well, such as better support for slicing operations.