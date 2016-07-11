Table of Contents
=================

* [01_1from10_bin_k31](#01_1from10_bin_k31)
  * [01_1from10_bin_k31/1.1_kmer_propagation.log](#01_1from10_bin_k3111_kmer_propagationlog)
  * [01_1from10_bin_k31/1.2_merging_fasta.log](#01_1from10_bin_k3112_merging_fastalog)
  * [01_1from10_bin_k31/2.1_bwa_fa2pac.log](#01_1from10_bin_k3121_bwa_fa2paclog)
  * [01_1from10_bin_k31/2.2_bwa_pac2bwtgen.log](#01_1from10_bin_k3122_bwa_pac2bwtgenlog)
  * [01_1from10_bin_k31/2.3_bwa_bwtupdate.log](#01_1from10_bin_k3123_bwa_bwtupdatelog)
  * [01_1from10_bin_k31/2.4_bwa_bwt2sa.log](#01_1from10_bin_k3124_bwa_bwt2salog)
  * [01_1from10_bin_k31/2.5_klcp.log](#01_1from10_bin_k3125_klcplog)
  * [01_1from10_bin_k31/3.1a_matching_rolling.log](#01_1from10_bin_k3131a_matching_rollinglog)
  * [01_1from10_bin_k31/3.1b_matching_rolling.log](#01_1from10_bin_k3131b_matching_rollinglog)
  * [01_1from10_bin_k31/3.2a_matching_restarted.log](#01_1from10_bin_k3132a_matching_restartedlog)
  * [01_1from10_bin_k31/3.2b_matching_restarted.log](#01_1from10_bin_k3132b_matching_restartedlog)
  * [01_1from10_bin_k31/3.3a_matching_rolling_skipping.log](#01_1from10_bin_k3133a_matching_rolling_skippinglog)
  * [01_1from10_bin_k31/3.3b_matching_rolling_skipping.log](#01_1from10_bin_k3133b_matching_rolling_skippinglog)
  * [01_1from10_bin_k31/3.4a_matching_restarted_skipping.log](#01_1from10_bin_k3134a_matching_restarted_skippinglog)
  * [01_1from10_bin_k31/3.4b_matching_restarted_skipping.log](#01_1from10_bin_k3134b_matching_restarted_skippinglog)
  * [01_1from10_bin_k31/4.1_contigs_stats.log](#01_1from10_bin_k3141_contigs_statslog)
  * [01_1from10_bin_k31/4.2_index_size.log](#01_1from10_bin_k3142_index_sizelog)
* [02_fullBacteria_orig_k31](#02_fullBacteria_orig_k31)
  * [02_fullBacteria_orig_k31/1.1_kmer_propagation.log](#02_fullbacteria_orig_k3111_kmer_propagationlog)
  * [02_fullBacteria_orig_k31/1.2_merging_fasta.log](#02_fullbacteria_orig_k3112_merging_fastalog)
  * [02_fullBacteria_orig_k31/2.1_bwa_fa2pac.log](#02_fullbacteria_orig_k3121_bwa_fa2paclog)
  * [02_fullBacteria_orig_k31/2.2_bwa_pac2bwtgen.log](#02_fullbacteria_orig_k3122_bwa_pac2bwtgenlog)
  * [02_fullBacteria_orig_k31/2.3_bwa_bwtupdate.log](#02_fullbacteria_orig_k3123_bwa_bwtupdatelog)
  * [02_fullBacteria_orig_k31/2.4_bwa_bwt2sa.log](#02_fullbacteria_orig_k3124_bwa_bwt2salog)
  * [02_fullBacteria_orig_k31/2.5_klcp.log](#02_fullbacteria_orig_k3125_klcplog)
  * [02_fullBacteria_orig_k31/3.1a_matching_rolling.log](#02_fullbacteria_orig_k3131a_matching_rollinglog)
  * [02_fullBacteria_orig_k31/3.1b_matching_rolling.log](#02_fullbacteria_orig_k3131b_matching_rollinglog)
  * [02_fullBacteria_orig_k31/3.2a_matching_restarted.log](#02_fullbacteria_orig_k3132a_matching_restartedlog)
  * [02_fullBacteria_orig_k31/3.2b_matching_restarted.log](#02_fullbacteria_orig_k3132b_matching_restartedlog)
  * [02_fullBacteria_orig_k31/3.3a_matching_rolling_skipping.log](#02_fullbacteria_orig_k3133a_matching_rolling_skippinglog)
  * [02_fullBacteria_orig_k31/3.3b_matching_rolling_skipping.log](#02_fullbacteria_orig_k3133b_matching_rolling_skippinglog)
  * [02_fullBacteria_orig_k31/3.4a_matching_restarted_skipping.log](#02_fullbacteria_orig_k3134a_matching_restarted_skippinglog)
  * [02_fullBacteria_orig_k31/3.4b_matching_restarted_skipping.log](#02_fullbacteria_orig_k3134b_matching_restarted_skippinglog)
  * [02_fullBacteria_orig_k31/4.1_contigs_stats.log](#02_fullbacteria_orig_k3141_contigs_statslog)
  * [02_fullBacteria_orig_k31/4.2_index_size.log](#02_fullbacteria_orig_k3142_index_sizelog)

***
## 01_1from10_bin_k31

### 01_1from10_bin_k31/1.1_kmer_propagation.log
* Fri Jul  8 14:45:40 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:37:23
System time:   0:00:56
Elapsed time:  0:02:26
CPU usage:        1568%

Memory peak:      1.12 GB
```

* 0inputs+3923136outputs (0major+34902647minor)pagefaults 0swaps


### 01_1from10_bin_k31/1.2_merging_fasta.log
* Fri Jul  8 14:48:07 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:54
System time:   0:00:02
Elapsed time:  0:00:59
CPU usage:          96%

Memory peak:      0.01 GB
```

* 0inputs+2539240outputs (0major+4569minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.1_bwa_fa2pac.log
* Fri Jul  8 14:49:06 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:12
System time:   0:00:01
Elapsed time:  0:00:23
CPU usage:          62%

Memory peak:      1.36 GB
```

* 88inputs+1836576outputs (0major+552129minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul  8 14:49:29 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen index.fa.pac index.fa.bwt

```
User time:     0:12:47
System time:   0:00:00
Elapsed time:  0:12:47
CPU usage:         100%

Memory peak:      0.61 GB
```

* 0inputs+1042696outputs (0major+12693minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.3_bwa_bwtupdate.log
* Fri Jul  8 15:02:17 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:04
System time:   0:00:01
Elapsed time:  0:00:14
CPU usage:          40%

Memory peak:      1.49 GB
```

* 64inputs+2085384outputs (0major+3149minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.4_bwa_bwt2sa.log
* Fri Jul  8 15:02:32 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:04:06
System time:   0:00:01
Elapsed time:  0:04:12
CPU usage:          98%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+27897minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.5_klcp.log
* Fri Jul  8 15:06:44 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     0:11:31
System time:   0:00:01
Elapsed time:  0:11:36
CPU usage:          99%

Memory peak:      1.49 GB
```

* 32inputs+1042696outputs (0major+47365minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1a_matching_rolling.log
* Fri Jul  8 15:18:21 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:29
System time:   0:00:02
Elapsed time:  0:01:31
CPU usage:          99%

Memory peak:      2.47 GB
```

* 0inputs+873456outputs (0major+368455minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1b_matching_rolling.log
* bwt_loading	0.50s
* sa_loading	0.24s
* bns_loading	4.50s
* klcp_loading	0.24s
* matching_time	85.59s
* reads	1000000
* kmers	70000000
* rpm	701052
* kpm	49073644


### 01_1from10_bin_k31/3.2a_matching_restarted.log
* Fri Jul  8 15:19:52 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:07
System time:   0:00:02
Elapsed time:  0:02:11
CPU usage:          99%

Memory peak:      1.97 GB
```

* 0inputs+873456outputs (0major+274488minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.44s
* sa_loading	0.28s
* bns_loading	4.65s
* matching_time	124.40s
* reads	1000000
* kmers	70000000
* rpm	482310
* kpm	33761722


### 01_1from10_bin_k31/3.3a_matching_rolling_skipping.log
* Fri Jul  8 15:22:04 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:28
System time:   0:00:02
Elapsed time:  0:01:32
CPU usage:          98%

Memory peak:      2.47 GB
```

* 16inputs+872584outputs (0major+271004minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.44s
* sa_loading	0.26s
* bns_loading	4.77s
* klcp_loading	0.20s
* matching_time	84.77s
* reads	1000000
* kmers	70000000
* rpm	707787
* kpm	49545114


### 01_1from10_bin_k31/3.4a_matching_restarted_skipping.log
* Fri Jul  8 15:23:36 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:07
System time:   0:00:01
Elapsed time:  0:02:10
CPU usage:          99%

Memory peak:      1.97 GB
```

* 0inputs+872584outputs (0major+395822minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.50s
* sa_loading	0.23s
* bns_loading	4.60s
* matching_time	123.91s
* reads	1000000
* kmers	70000000
* rpm	484213
* kpm	33894913


### 01_1from10_bin_k31/4.1_contigs_stats.log
* Number of contigs: 8495372
* Total length: 1067712578
* Average length: 125.68167444580413

* Median length: 34.0
* Number of k-mer occurencies: 812851418


### 01_1from10_bin_k31/4.2_index_size.log
* 510M	index.fa.31.bit.klcp
* 388M	index.fa.ann
* 1019M	index.fa.bwt
* 510M	index.fa.pac
* 510M	index.fa.sa

***
## 02_fullBacteria_orig_k31

### 02_fullBacteria_orig_k31/1.1_kmer_propagation.log
* Fri Jul  8 16:00:44 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     4:53:31
System time:   0:08:25
Elapsed time:  0:16:14
CPU usage:        1858%

Memory peak:      9.31 GB
```

* 17825456inputs+34175048outputs (7major+248658426minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/1.2_merging_fasta.log
* Fri Jul  8 16:16:59 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:09:45
System time:   0:00:21
Elapsed time:  0:11:54
CPU usage:          85%

Memory peak:      0.01 GB
```

* 11455984inputs+24270152outputs (21major+43634minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.1_bwa_fa2pac.log
* Fri Jul  8 16:28:53 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:05
System time:   0:00:23
Elapsed time:  0:04:00
CPU usage:          62%

Memory peak:     14.08 GB
```

* 672inputs+18388752outputs (1major+7863420minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul  8 16:32:54 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen index.fa.pac index.fa.bwt

```
User time:     3:15:23
System time:   0:00:21
Elapsed time:  3:16:08
CPU usage:         100%

Memory peak:      5.41 GB
```

* 3197600inputs+9968312outputs (2major+16324858minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.3_bwa_bwtupdate.log
* Fri Jul  8 19:49:02 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:31
Elapsed time:  0:02:47
CPU usage:          44%

Memory peak:     14.26 GB
```

* 840inputs+19936624outputs (3major+3793035minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.4_bwa_bwt2sa.log
* Fri Jul  8 19:51:50 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:42:56
System time:   0:01:02
Elapsed time:  0:44:46
CPU usage:          98%

Memory peak:     14.26 GB
```

* 288inputs+9968312outputs (0major+54060232minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.5_klcp.log
* Fri Jul  8 20:36:36 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     2:00:06
System time:   0:03:13
Elapsed time:  2:04:06
CPU usage:          99%

Memory peak:     14.26 GB
```

* 1768inputs+9968312outputs (4major+302487638minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1a_matching_rolling.log
* Fri Jul  8 22:40:43 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:48
System time:   0:01:12
Elapsed time:  0:05:19
CPU usage:          94%

Memory peak:     22.83 GB
```

* 13699656inputs+3943496outputs (2major+19263970minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	4.05s
* sa_loading	8.30s
* bns_loading	60.30s
* klcp_loading	27.21s
* matching_time	199.65s
* reads	1000000
* kmers	70000000
* rpm	300521
* kpm	21036448


### 02_fullBacteria_orig_k31/3.2a_matching_restarted.log
* Fri Jul  8 22:46:02 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:55
System time:   0:00:34
Elapsed time:  0:11:20
CPU usage:          93%

Memory peak:     18.08 GB
```

* 18506832inputs+3943496outputs (4major+15578537minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	8.88s
* sa_loading	5.58s
* bns_loading	58.28s
* matching_time	555.96s
* reads	1000000
* kmers	70000000
* rpm	107921
* kpm	7554446


### 02_fullBacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Fri Jul  8 22:57:23 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:01
System time:   0:01:29
Elapsed time:  0:06:10
CPU usage:          89%

Memory peak:     22.83 GB
```

* 12680800inputs+3942912outputs (3major+57484780minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	19.79s
* sa_loading	12.34s
* bns_loading	59.31s
* klcp_loading	15.01s
* matching_time	222.54s
* reads	1000000
* kmers	70000000
* rpm	269610
* kpm	18872729


### 02_fullBacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Fri Jul  8 23:03:33 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:45
System time:   0:01:08
Elapsed time:  0:11:22
CPU usage:          96%

Memory peak:     18.08 GB
```

* 9132736inputs+3942912outputs (0major+55222188minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	16.14s
* sa_loading	11.38s
* bns_loading	59.09s
* matching_time	565.96s
* reads	1000000
* kmers	70000000
* rpm	106015
* kpm	7421022


### 02_fullBacteria_orig_k31/4.1_contigs_stats.log
* Number of contigs: 96088810
* Total length: 10207548814
* Average length: 106.23035933112294

* Median length: 34.0
* Number of k-mer occurencies: 7324884514


### 02_fullBacteria_orig_k31/4.2_index_size.log
* 4.8G	index.fa.31.bit.klcp
* 4.1G	index.fa.ann
* 9.6G	index.fa.bwt
* 4.8G	index.fa.pac
* 4.8G	index.fa.sa

