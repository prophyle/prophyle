#ifndef BWTEXK_H
#define BWTEXK_H

#include <stdint.h>
#include "bwt.h"
#include "bwtaln.h"
#include "bwa.h"
#include "klcp.h"
#include "prophyle_utils.h"

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint64_t position;
	int strand;
	int rid;
	int node;
} bwt_position_t;

	exk_opt_t *exk_init_opt();
  bwa_seqio_t *bwa_seq_open(const char *fn);
  bwa_seqio_t *bwa_bam_open(const char *fn, int which);
  void bwa_seq_close(bwa_seqio_t *bs);
  void seq_reverse(int len, ubyte_t *seq, int is_comp);
  bwa_seq_t *bwa_read_seq(bwa_seqio_t *seq, int n_needed, int *n, int mode, int trim_qual);

	void bwa_exk_core(const char *prefix, const char *fn_fa, const exk_opt_t *opt);

	void bwa_cal_sa(bwaidx_t* idx, int n_seqs, bwa_seq_t *seqs, const exk_opt_t *opt, klcp_t* klcp);

#ifdef __cplusplus
}
#endif

#endif
