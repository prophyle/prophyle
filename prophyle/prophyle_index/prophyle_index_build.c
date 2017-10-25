#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "prophyle_index_build.h"
#include "klcp.h"
#include "prophyle_utils.h"
#include "bwa_utils.h"
#include "utils.h"

typedef struct {
	klcp_t* klcp;
	bwt_t* bwt;
	int kmer_length;
	const char* prefix;
	int sa_intv;
} klcp_data_t;

void* construct_klcp_parallel(void* data) {
	klcp_data_t* klcp_data = (klcp_data_t*)data;
	klcp_data->klcp = construct_klcp(klcp_data->bwt, klcp_data->kmer_length);
	return 0;
}

void* construct_sa_parallel(void* data) {
	klcp_data_t* klcp_data = (klcp_data_t*)data;
	bwt_cal_sa(klcp_data->bwt, klcp_data->sa_intv);
	char* fn = malloc((strlen(klcp_data->prefix) + 10) * sizeof(char));
	strcpy(fn, klcp_data->prefix);
	strcat(fn, ".sa");
	bwt_dump_sa(fn, klcp_data->bwt);
	fprintf(stderr, "[prophyle_index:%s] SA dumped\n", __func__);
	return 0;
}

void build_index(const char *prefix, const prophyle_index_opt_t *opt, int sa_intv) {
	bwt_t *bwt;
	{
		if ((bwt = bwa_idx_load_bwt_without_sa(prefix)) == 0) {
			fprintf(stderr, "[prophyle_index:%s] Couldn't load idx from %s\n", __func__, prefix);
			return;
		}
	}

	klcp_t* klcp;
	if (opt->construct_sa_parallel) {
		klcp_data_t* klcp_data = malloc(sizeof(klcp_data_t));
		klcp_data->bwt = bwt;
		klcp_data->kmer_length = opt->kmer_length;
		klcp_data->prefix = prefix;
		klcp_data->sa_intv = sa_intv;
		pthread_t tid[2];
		int status_klcp = pthread_create(&tid[0], NULL, construct_klcp_parallel, (void*)klcp_data);
		int status_sa = pthread_create(&tid[1], NULL, construct_sa_parallel, (void*)klcp_data);
		xassert(!status_klcp, "[prophyle_index] error while creating thread for klcp parallel construction, try construction separate from sa\n");
		xassert(!status_sa, "[prophyle_index] error while creating thread for sa parallel construction, try construction separate from klcp\n");
		fprintf(stderr, "[prophyle_index] parallel construction of klcp and sa started\n");
		int status_addr_klcp = pthread_join(tid[0], (void**)&status_addr_klcp);
		int status_addr_sa = pthread_join(tid[1], (void**)&status_addr_sa);
		xassert(!status_addr_klcp, "[prophyle_index] error while klcp parallel construction, try construction separate from sa\n");
		xassert(!status_addr_sa, "[prophyle_index] error sa parallel construction, try construction separate from klcp\n");
		klcp = klcp_data->klcp;
	} else {
		klcp = construct_klcp(bwt, opt->kmer_length);
	}
	char* fn = malloc((strlen(prefix) + 10) * sizeof(char));
	strcpy(fn, prefix);
	strcat(fn, ".");
	char* kmer_length_str = malloc(5 * sizeof(char));
	sprintf(kmer_length_str, "%d", opt->kmer_length);
	strcat(fn, kmer_length_str);
	strcat(fn, ".klcp");
	klcp_dump(fn, klcp);
	fprintf(stderr, "[prophyle_index:%s] klcp dumped\n", __func__);
	if (opt->construct_sa_parallel) {
		bwt_destroy(bwt);
	} else {
		bwt_destroy_without_sa(bwt);
	}
}

int debwtupdate(const char* bwt_input_file, const char* bwt_output_file) {
	bwtint_t i, k, n_occ;
	uint32_t *buf;
	bwt_t *bwt = bwt_restore_bwt(bwt_input_file);
	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	// fprintf(stderr, "seq_len: %d n_occ: %d old_size: %d new_size: %d\n",
	// 	bwt->seq_len, n_occ, bwt->bwt_size, bwt->bwt_size - n_occ * sizeof(bwtint_t));
	bwt->bwt_size -= n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	// c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		// fprintf(stderr, "i: %d k: %d\n", i, k);
		if (i % OCC_INTERVAL == 0) {
			// memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[i/16] = bwt->bwt[k++]; // 16 == sizeof(uint32_t)/2
		// ++c[bwt_B00(bwt, i)];
	}
	// fprintf(stderr, "end loop\n");
	free(bwt->bwt);
	bwt->bwt = buf;
	// fprintf(stderr, "begore dump\n");
	bwt_dump_bwt(bwt_output_file, bwt);
	// fprintf(stderr, "after dump\n");
	bwt_destroy(bwt);
	return 0;
}
