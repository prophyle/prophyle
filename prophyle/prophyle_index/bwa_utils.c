#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "utils.h"
#include "bwa.h"
#include "kstring.h"
#include "prophyle_utils.h"
#include "khash.h"
#include "contig_node_translator.h"

KHASH_MAP_INIT_STR(str, int)

void bwa_destroy_unused_fields(bwaidx_t* idx) {
	int64_t i;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		free(idx->bns->anns[i].name);
		free(idx->bns->anns[i].anno);
	}
	if (idx->pac) {
		free(idx->pac);
	}
}

void bns_destroy_without_names_and_anno(bntseq_t* bns) {
	if (bns == 0) return;
	else {
		if (bns->fp_pac) err_fclose(bns->fp_pac);
		free(bns->ambs);
		free(bns->anns);
		free(bns);
	}
}

void bwa_idx_destroy_without_bns_name_and_anno(bwaidx_t* idx)
{
	if (idx == 0) return;
	if (idx->mem == 0) {
		if (idx->bwt) bwt_destroy(idx->bwt);
		if (idx->bns) bns_destroy_without_names_and_anno(idx->bns);
	} else {
		free(idx->bwt); free(idx->bns->anns); free(idx->bns);
		if (!idx->is_shm) free(idx->mem);
	}
	free(idx);
}

bntseq_t* bns_restore_core_partial(const char* ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[8192];
	FILE* fp;
	const char* fname;
	bntseq_t* bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		if (scanres != 3) goto badread;
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t* p = bns->anns + i;
			char* q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			if (scanres != 2) goto badread;

			add_contig(str, i);

			// read fasta comments
			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);
			if (c == EOF) {
				scanres = EOF;
				goto badread;
			}
			*q = 0;

			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			if (scanres != 3) goto badread;
			p->offset = xx;
		}
		err_fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		if (scanres != 3) goto badread;
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t* p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			if (scanres != 3) goto badread;
			p->offset = xx;
			p->amb = str[0];
		}
		err_fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;

badread:
	if (EOF == scanres) {
		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fname);
}

bntseq_t* bns_restore_partial(const char* prefix)
{
	char ann_filename[1024], amb_filename[1024], pac_filename[1024], alt_filename[1024];
	FILE* fp;
	bntseq_t* bns;
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	bns = bns_restore_core_partial(ann_filename, amb_filename, pac_filename);
	if (bns == 0) return 0;
	if ((fp = fopen(strcat(strcpy(alt_filename, prefix), ".alt"), "r")) != 0) { // read .alt file if present
		fprintf(stderr, "[prophyle_index:%s] .alt file is present, something may work wrong!\n", __func__);
		char str[1024];
		khash_t(str)* h;
		int c, i, absent;
		khint_t k;
		h = kh_init(str);
		for (i = 0; i < bns->n_seqs; ++i) {
			k = kh_put(str, h, bns->anns[i].name, &absent);
			kh_val(h, k) = i;
		}
		i = 0;
		while ((c = fgetc(fp)) != EOF) {
			if (c == '\t' || c == '\n' || c == '\r') {
				str[i] = 0;
				if (str[0] != '@') {
					k = kh_get(str, h, str);
					if (k != kh_end(h))
						bns->anns[kh_val(h, k)].is_alt = 1;
				}
				while (c != '\n' && c != EOF) c = fgetc(fp);
				i = 0;
			} else str[i++] = c; // FIXME: potential segfault here
		}
		kh_destroy(str, h);
		fclose(fp);
	}
	return bns;
}

bwt_t* bwa_idx_load_bwt_with_time(const char* hint, int need_log, FILE* log_file)
{
	char* tmp;
	char* prefix;
	bwt_t* bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[prophyle_index:%s] fail to locate the index files\n", __func__);
		return 0;
	}
	clock_t t = clock();
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt");
	bwt = bwt_restore_bwt(tmp);
	if (need_log) {
		fprintf(log_file, "bwt_loading\t%.2fs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	t = clock();
	strcat(strcpy(tmp, prefix), ".sa");
	bwt_restore_sa(tmp, bwt);
	if (need_log) {
		fprintf(log_file, "sa_loading\t%.2fs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	free(tmp); free(prefix);
	return bwt;
}

bwaidx_t* bwa_idx_load_partial(const char* hint, int which, int need_log, FILE* log_file)
{
	bwaidx_t* idx;
	char* prefix;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[prophyle_index::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = calloc(1, sizeof(bwaidx_t));
	if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt_with_time(hint, need_log, log_file);
	if (which & BWA_IDX_BNS) {
		int i, c;
		clock_t t = clock();
		idx->bns = bns_restore_partial(prefix);
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		if (need_log) {
			fprintf(log_file, "bns_loading\t%.2fs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		}
	}
	free(prefix);
	return idx;
}

bwt_t* bwa_idx_load_bwt_without_sa(const char* hint)
{
	char* tmp;
	char* prefix;
	bwt_t* bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[prophyle_index::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt");
	bwt = bwt_restore_bwt(tmp);
	free(tmp); free(prefix);
	return bwt;
}

void bwt_destroy_without_sa(bwt_t* bwt)
{
	if (bwt == 0) return;
	free(bwt->bwt);
	free(bwt);
}
