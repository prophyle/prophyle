#!/usr/bin/env python3

import os
import sys
import pysam
import argparse
import numpy as np

from ete3 import Tree
from sklearn.linear_model import ElasticNetCV


def analyse_sam(sam_fn, tree):

    nodes2leaves = {
        node.name: {
            leaf.name for leaf in node
        } for node in tree.traverse("postorder")
    }
    leaves = [leaf.name for leaf in tree]
    leaves_set = set(leaves)
    ref_size = len(leaves)

    total = 0
    vec_unique = np.zeros(ref_size)
    vec_shared = np.zeros(ref_size)
    vec_pos = {leaf: i for i, leaf in enumerate(leaves)}

    sam_f = pysam.AlignmentFile(sam_fn)
    prev_read_name = ""
    cur_ref = []
    for read in sam_f.fetch(until_eof=True):
        if not read.is_unmapped:
            read_name = read.qname
            read_ref = read.reference_name
            if read_name == prev_read_name:
                cur_ref.append(read_ref)
            else:
                # first read
                if len(cur_ref) == 0:
                    pass
                # unique assignment
                elif len(cur_ref) == 1:
                    ref = cur_ref[0]
                    # leaf (REAL unique assignment)
                    if ref in leaves_set:
                        vec_unique[vec_pos[ref]] += 1
                    # internal node (multiple ass to all its leaves)
                    else:
                        for leaf in nodes2leaves[ref]:
                            vec_shared[vec_pos[leaf]] += 1
                else:
                    for ref in cur_ref:
                        for leaf in nodes2leaves[ref]:
                            vec_shared[vec_pos[leaf]] += 1

                total += 1
                cur_ref = [read_ref]
                prev_read_name = read_name

    # last assignment
    # unique assignment
    if len(cur_ref) == 1:
        total += 1
        ref = cur_ref[0]
        # leaf (REAL unique assignment)
        if ref in leaves_set:
            vec_unique[vec_pos[ref]] += 1
        # internal node (multiple ass to all its leaves)
        else:
            for leaf in nodes2leaves[ref]:
                vec_shared[vec_pos[leaf]] += 1
    elif len(cur_ref) > 0:
        total += 1
        for ref in cur_ref:
            for leaf in nodes2leaves[ref]:
                vec_shared[vec_pos[leaf]] += 1

    sam_f.close()

    return(vec_unique, vec_shared, total)


def main():

    parser = argparse.ArgumentParser(
        description='Estimate abundances from ProPhyle pseudoalignments using a Elastic Net-based GLM'
    )

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nw>',
        help='taxonomic tree (tree.preliminary.nw)'
    )

    parser.add_argument(
        'sam_fn',
        type=str,
        metavar='<pseudo_aln.sam>',
        help='output of prophyle'
    )

    parser.add_argument(
        'sim_mat_fn',
        type=str,
        metavar='<sim_matrix.npy>',
        help='similarity matrix'
    )

    parser.add_argument(
        'out_fn',
        type=str,
        metavar='<output_fn>',
        help='output file'
    )

    args = parser.parse_args()

    tree_fn = args.tree_fn
    sam_fn = args.sam_fn
    sim_mat_fn = args.sim_mat_fn
    out_fn = args.out_fn

    count_fn = '.'.join(args.sam_fn.split('.')[:-1]+['npy'])
    tree = Tree(tree_fn, format=1)
    leaves = [leaf.name for leaf in tree]

    if os.path.isfile(count_fn):
        print("Loading counts from existing npy file", file=sys.stderr)
        ## tsv version
        # map_counts = np.array()
        # with open(count_fn, 'r') as count_f:
        #     for line in count_f:
        #         node, count = map(str.strip, line.split('\t'))
        #         np.append(map_counts, float(count))
        map_counts = np.load(count_fn)

    else:
        unique, shared, total = analyse_sam(sam_fn, tree)
        map_counts = unique + shared
        np.save(count_fn, map_counts)

    assert len(leaves) == len(map_counts)

    temp_sim_mat = np.load(sim_mat_fn)
    assert len(leaves) == len(temp_sim_mat)
    sim_mat = np.zeros((len(map_counts), len(map_counts)))

    # sim_mat[i,j] = fraction (0<f<1) of reads from ref j mapped to ref i
    for i in range(len(map_counts)):
        # sim_mat[i, :] = temp_sim_mat[i, :] / temp_sim_mat[i, i]
        for j in range(len(map_counts)):
            sim_mat[i, j] = temp_sim_mat[j, i] / temp_sim_mat[j, j]

    # groups = [str(i) for i in range(len(map_counts))]
    #
    # fam = Poisson()
    # ind = Independence()
    # gee_model = GEE(map_counts, sim_mat, groups, family=fam, cov_struct=ind)
    # result = gee_model.fit()
    # print(result.summary())

    # alphas_positive_enet, coefs_positive_enet, other = enet_path(S, m, eps=5e-5, l1_ratio=0.1, positive=True, fit_intercept=False)

    # enet = ElasticNet(
    #     alpha=0.8,
    #     l1_ratio=0.3,
    #     fit_intercept=False,
    #     max_iter=10000,
    #     normalize=True,
    #     positive=True,
    #     precompute=True
    # )

    ratios = np.array([.05, .1, .3, .5, .6, .7, .8, .9, .93, .95, .97, .99, 1])
    ## Cross validation version
    enet_cv = ElasticNetCV(cv=3, verbose=True, l1_ratio=ratios, n_jobs=8, selection='random', fit_intercept=False, positive=True)
    enet_cv.fit(sim_mat, map_counts)

    print("Model score: {} (1 is best)".format(enet_cv.score(sim_mat, map_counts)), file=sys.stderr)
    print("Number of Iterations: {}".format(enet_cv.n_iter_))
    print("Best l1 ratio: {}".format(enet_cv.l1_ratio_))
    print("Best alpha: {}".format(enet_cv.alpha_))
    print("(tested grid saved in alphas.npy and ratios.npy)")

    with open('alphas.npy', 'wb') as out_f:
        np.save(out_f, enet_cv.alphas_)
    with open('ratios.npy', 'wb') as out_f:
        np.save(out_f, ratios)

    assert enet_cv.intercept_ == 0.

    with open(out_fn, 'w') as out_f:
        for leaf, ab in zip(leaves, enet_cv.coef_):
            print(leaf, ab, sep='\t', file=out_f)


if __name__ == '__main__':
    main()
