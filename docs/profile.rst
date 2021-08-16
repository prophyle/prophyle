.. _profile:

Estimating abundances
=====================

After classification, abundances can be estimated with the following command:

.. code-block:: bash

	prophyle profile <index_dir> <assignments.sam> <output_prefix>

The similarity matrix for the linear model is expected to be in the index
directory under the name ``sim_mat.npy``, which is the default output of
``prophyle similarity`` and the default path in pre-built indexes; if this
is not the case, the path can be specified using the option ``-s``.

The option ``-a`` (float>=0) can be used to adjust the regularization effect, while
the option ``-l`` (0<=float<=1) modifies the ratio between L1 and L2 regularization
penalties. For more information, see the scikit-learn reference for
`Elastic Net <http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html>`_.
