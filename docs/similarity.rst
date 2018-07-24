.. _similarity:

Build similarity matrix
=======================

In order to estimate abundances from raw assignment, ProPhyle needs a similarity
matrix for the genomes in the index. Although a default matrix is shipped
together with the pre-built indexes, you may need to adjust parameters to match
the sequencing technology used for your experiment. The following command will
create `Snakefiles <https://snakemake.readthedocs.io/en/stable/>`_ for simulating
reads from the reference genomes and classifying them with the index, in order to
assess the assignments distribution:

.. code-block:: bash

	prophyle similarity <index_dir> -o <sim_mat.npy> -l <read_len>

This will create the Snakefiles without executing them. Add the ``-R -j <n_jobs>``
option to run ``snakemake`` automatically.
Since this step requires a considerable amount of computational resources, we
recommend to run it in a cluster environment. You can follow the instruction at
`this page <https://snakemake.readthedocs.io/en/stable/executable.html?highlight=cluster#cluster-execution>`_,
depending on your cluster environment.
After running ``snakemake``, you can run the command again with ``-R`` to complete
the process.
