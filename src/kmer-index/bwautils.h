#ifndef BWAUTILS_H
#define BWAUTILS_H

#include <stdint.h>
#include "bwt.h"
#include "bwtaln.h"
#include "bwa.h"
#include "exkutils.h"

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

bwa_seqio_t *bwa_open_reads_new(int mode, const char *fn_fa);
void bwa_destroy_unused_fields(bwaidx_t* idx);
void bns_destroy_without_names_and_annos(bntseq_t* bns);
void bwa_idx_destroy_without_bns_name_and_anno(bwaidx_t *idx);
bntseq_t *bns_restore_core_partial(const char *ann_filename, const char* amb_filename, const char* pac_filename);
bntseq_t *bns_restore_partial(const char *prefix);
bwaidx_t *bwa_idx_load_partial(const char *hint, int which);

#endif // BWAUTILS_H
