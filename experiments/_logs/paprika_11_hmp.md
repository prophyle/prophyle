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
* [03_fullHMP_k31](#03_fullHMP_k31)
  * [03_fullHMP_k31/1.1_kmer_propagation.log](#03_fullhmp_k3111_kmer_propagationlog)
  * [03_fullHMP_k31/1.2_merging_fasta.log](#03_fullhmp_k3112_merging_fastalog)
  * [03_fullHMP_k31/2.1_bwa_fa2pac.log](#03_fullhmp_k3121_bwa_fa2paclog)
  * [03_fullHMP_k31/2.2_bwa_pac2bwtgen.log](#03_fullhmp_k3122_bwa_pac2bwtgenlog)
  * [03_fullHMP_k31/2.3_bwa_bwtupdate.log](#03_fullhmp_k3123_bwa_bwtupdatelog)
  * [03_fullHMP_k31/2.4_bwa_bwt2sa.log](#03_fullhmp_k3124_bwa_bwt2salog)
  * [03_fullHMP_k31/2.5_klcp.log](#03_fullhmp_k3125_klcplog)
  * [03_fullHMP_k31/3.1a_matching_rolling.log](#03_fullhmp_k3131a_matching_rollinglog)
  * [03_fullHMP_k31/3.1b_matching_rolling.log](#03_fullhmp_k3131b_matching_rollinglog)
  * [03_fullHMP_k31/3.2a_matching_restarted.log](#03_fullhmp_k3132a_matching_restartedlog)
  * [03_fullHMP_k31/3.2b_matching_restarted.log](#03_fullhmp_k3132b_matching_restartedlog)
  * [03_fullHMP_k31/3.3a_matching_rolling_skipping.log](#03_fullhmp_k3133a_matching_rolling_skippinglog)
  * [03_fullHMP_k31/3.3b_matching_rolling_skipping.log](#03_fullhmp_k3133b_matching_rolling_skippinglog)
  * [03_fullHMP_k31/3.4a_matching_restarted_skipping.log](#03_fullhmp_k3134a_matching_restarted_skippinglog)
  * [03_fullHMP_k31/3.4b_matching_restarted_skipping.log](#03_fullhmp_k3134b_matching_restarted_skippinglog)
  * [03_fullHMP_k31/4.1_contigs_stats.log](#03_fullhmp_k3141_contigs_statslog)
  * [03_fullHMP_k31/4.2_index_size.log](#03_fullhmp_k3142_index_sizelog)

***
## 01_1from10_bin_k31

### 01_1from10_bin_k31/1.1_kmer_propagation.log
* Thu Jul 14 11:51:32 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:37:17
System time:   0:00:54
Elapsed time:  0:02:38
CPU usage:        1445%

Memory peak:      1.12 GB
```

* 1992296inputs+5007784outputs (21major+34058317minor)pagefaults 0swaps


### 01_1from10_bin_k31/1.2_merging_fasta.log
* Thu Jul 14 11:54:11 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:54
System time:   0:00:02
Elapsed time:  0:00:59
CPU usage:          95%

Memory peak:      0.01 GB
```

* 0inputs+2539248outputs (0major+4536minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.1_bwa_fa2pac.log
* Thu Jul 14 11:55:10 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:12
System time:   0:00:01
Elapsed time:  0:00:23
CPU usage:          59%

Memory peak:      1.36 GB
```

* 224inputs+1836576outputs (1major+521668minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.2_bwa_pac2bwtgen.log
* Thu Jul 14 11:55:34 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     0:12:52
System time:   0:00:00
Elapsed time:  0:12:53
CPU usage:         100%

Memory peak:      0.79 GB
```

* 0inputs+1042696outputs (0major+21570minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.3_bwa_bwtupdate.log
* Thu Jul 14 12:08:27 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:04
System time:   0:00:01
Elapsed time:  0:00:14
CPU usage:          40%

Memory peak:      1.49 GB
```

* 72inputs+2085384outputs (0major+2691minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.4_bwa_bwt2sa.log
* Thu Jul 14 12:08:42 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:04:09
System time:   0:00:01
Elapsed time:  0:04:15
CPU usage:          98%

Memory peak:      1.49 GB
```

* 32inputs+1042696outputs (0major+36490minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.5_klcp.log
* Thu Jul 14 12:12:57 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     0:11:26
System time:   0:00:01
Elapsed time:  0:11:32
CPU usage:          99%

Memory peak:      1.49 GB
```

* 32inputs+1042696outputs (0major+46970minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1a_matching_rolling.log
* Thu Jul 14 13:45:48 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:25
System time:   0:00:01
Elapsed time:  0:01:28
CPU usage:          99%

Memory peak:      2.77 GB
```

* 397160inputs+221576outputs (0major+466930minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1b_matching_rolling.log
* bwt_loading	0.53s
* sa_loading	0.24s
* bns_loading	5.17s
* klcp_loading	0.23s
* matching_time	81.14s
* reads	1000000
* kmers	70000000
* rpm	739448
* kpm	51761394


### 01_1from10_bin_k31/3.2a_matching_restarted.log
* Thu Jul 14 13:47:17 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:01
System time:   0:00:01
Elapsed time:  0:02:03
CPU usage:         100%

Memory peak:      2.27 GB
```

* 0inputs+221576outputs (0major+339248minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.43s
* sa_loading	0.21s
* bns_loading	5.43s
* matching_time	116.78s
* reads	1000000
* kmers	70000000
* rpm	513797
* kpm	35965789


### 01_1from10_bin_k31/3.3a_matching_rolling_skipping.log
* Thu Jul 14 13:49:20 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:22
System time:   0:00:01
Elapsed time:  0:01:25
CPU usage:         100%

Memory peak:      2.77 GB
```

* 0inputs+221576outputs (0major+432595minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.50s
* sa_loading	0.23s
* bns_loading	5.31s
* klcp_loading	0.22s
* matching_time	78.41s
* reads	1000000
* kmers	70000000
* rpm	765254
* kpm	53567784


### 01_1from10_bin_k31/3.4a_matching_restarted_skipping.log
* Thu Jul 14 13:50:45 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:07
System time:   0:00:01
Elapsed time:  0:02:09
CPU usage:         100%

Memory peak:      2.27 GB
```

* 0inputs+221576outputs (0major+422602minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.45s
* sa_loading	0.21s
* bns_loading	5.42s
* matching_time	123.40s
* reads	1000000
* kmers	70000000
* rpm	486243
* kpm	34037013


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
* Thu Jul 14 12:32:06 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     4:54:50
System time:   0:08:47
Elapsed time:  0:16:51
CPU usage:        1801%

Memory peak:      9.31 GB
```

* 20591904inputs+44018032outputs (4major+255127124minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/1.2_merging_fasta.log
* Thu Jul 14 12:48:58 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:09:51
System time:   0:00:21
Elapsed time:  0:11:59
CPU usage:          85%

Memory peak:      0.01 GB
```

* 9984840inputs+24270160outputs (21major+31370minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.1_bwa_fa2pac.log
* Thu Jul 14 13:00:58 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:05
System time:   0:00:22
Elapsed time:  0:03:57
CPU usage:          62%

Memory peak:     14.08 GB
```

* 560inputs+18388752outputs (1major+7682427minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Thu Jul 14 13:52:55 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:24:52
System time:   0:00:07
Elapsed time:  2:25:07
CPU usage:         100%

Memory peak:      5.58 GB
```

* 1893768inputs+9968312outputs (0major+370064minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.3_bwa_bwtupdate.log
* Thu Jul 14 16:18:03 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:36
Elapsed time:  0:02:52
CPU usage:          46%

Memory peak:     14.26 GB
```

* 536inputs+19936624outputs (0major+4418787minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.4_bwa_bwt2sa.log
* Thu Jul 14 16:20:55 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:50:03
System time:   0:01:45
Elapsed time:  0:52:42
CPU usage:          98%

Memory peak:     14.26 GB
```

* 240inputs+9968312outputs (0major+123254907minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.5_klcp.log
* Thu Jul 14 17:13:38 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     2:24:10
System time:   0:05:33
Elapsed time:  2:30:36
CPU usage:          99%

Memory peak:     14.26 GB
```

* 288inputs+9968312outputs (0major+441120609minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1a_matching_rolling.log
* Thu Jul 14 19:44:14 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:15
System time:   0:05:35
Elapsed time:  0:08:56
CPU usage:          99%

Memory peak:     26.06 GB
```

* 9589104inputs+400864outputs (0major+30271123minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	6.48s
* sa_loading	11.72s
* bns_loading	145.50s
* klcp_loading	197.52s
* matching_time	167.15s
* reads	1000000
* kmers	70000000
* rpm	358952
* kpm	25126673


### 02_fullBacteria_orig_k31/3.2a_matching_restarted.log
* Thu Jul 14 19:53:11 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:26
System time:   0:01:37
Elapsed time:  0:11:23
CPU usage:          97%

Memory peak:     21.30 GB
```

* 12768264inputs+400864outputs (4major+89420755minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	27.59s
* sa_loading	9.49s
* bns_loading	65.45s
* matching_time	559.09s
* reads	1000000
* kmers	70000000
* rpm	107318
* kpm	7512268


### 02_fullBacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Thu Jul 14 20:04:34 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:12
System time:   0:00:47
Elapsed time:  0:05:18
CPU usage:          75%

Memory peak:     26.06 GB
```

* 21341016inputs+400864outputs (8major+8717431minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	10.16s
* sa_loading	6.48s
* bns_loading	66.94s
* klcp_loading	15.79s
* matching_time	138.74s
* reads	1000000
* kmers	70000000
* rpm	432457
* kpm	30271961


### 02_fullBacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Thu Jul 14 20:09:52 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:54
System time:   0:00:22
Elapsed time:  0:11:35
CPU usage:          89%

Memory peak:     21.31 GB
```

* 21871432inputs+400864outputs (3major+2418519minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	12.09s
* sa_loading	4.70s
* bns_loading	60.48s
* matching_time	539.16s
* reads	1000000
* kmers	70000000
* rpm	111284
* kpm	7789904


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

***
## 03_fullHMP_k31

### 03_fullHMP_k31/1.1_kmer_propagation.log
* Thu Jul 14 20:21:32 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     2:51:30
System time:   0:04:44
Elapsed time:  0:11:48
CPU usage:        1492%

Memory peak:      5.28 GB
```

* 12902952inputs+25577768outputs (16major+148523967minor)pagefaults 0swaps


### 03_fullHMP_k31/1.2_merging_fasta.log
* Thu Jul 14 20:33:21 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:04:34
System time:   0:00:11
Elapsed time:  0:05:19
CPU usage:          90%

Memory peak:      0.01 GB
```

* 30640inputs+13586264outputs (2major+17785minor)pagefaults 0swaps


### 03_fullHMP_k31/2.1_bwa_fa2pac.log
* Thu Jul 14 20:38:40 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:00
System time:   0:00:10
Elapsed time:  0:01:58
CPU usage:          60%

Memory peak:      6.15 GB
```

* 1288inputs+8866880outputs (4major+1460529minor)pagefaults 0swaps


### 03_fullHMP_k31/2.2_bwa_pac2bwtgen.log
* Thu Jul 14 20:40:38 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     1:25:13
System time:   0:00:06
Elapsed time:  1:25:19
CPU usage:         100%

Memory peak:      3.48 GB
```

* 0inputs+6060464outputs (0major+2022974minor)pagefaults 0swaps


### 03_fullHMP_k31/2.3_bwa_bwtupdate.log
* Thu Jul 14 22:05:58 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:25
System time:   0:00:20
Elapsed time:  0:01:54
CPU usage:          40%

Memory peak:      8.67 GB
```

* 384inputs+12120928outputs (0major+2621194minor)pagefaults 0swaps


### 03_fullHMP_k31/2.4_bwa_bwt2sa.log
* Thu Jul 14 22:07:52 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:34:18
System time:   0:01:28
Elapsed time:  0:36:18
CPU usage:          99%

Memory peak:      8.67 GB
```

* 248inputs+6060464outputs (1major+83809277minor)pagefaults 0swaps


### 03_fullHMP_k31/2.5_klcp.log
* Thu Jul 14 22:44:11 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     1:03:06
System time:   0:01:34
Elapsed time:  1:05:13
CPU usage:          99%

Memory peak:      8.67 GB
```

* 384inputs+6060464outputs (1major+46871436minor)pagefaults 0swaps


### 03_fullHMP_k31/3.1a_matching_rolling.log
* Thu Jul 14 23:49:24 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:12
System time:   0:00:37
Elapsed time:  0:03:06
CPU usage:          91%

Memory peak:     14.12 GB
```

* 8501264inputs+305736outputs (3major+14707144minor)pagefaults 0swaps


### 03_fullHMP_k31/3.1b_matching_rolling.log
* bwt_loading	7.54s
* sa_loading	7.05s
* bns_loading	23.00s
* klcp_loading	11.20s
* matching_time	120.67s
* reads	1000000
* kmers	70000000
* rpm	497230
* kpm	34806127


### 03_fullHMP_k31/3.2a_matching_restarted.log
* Thu Jul 14 23:52:31 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:58
System time:   0:00:19
Elapsed time:  0:05:19
CPU usage:         100%

Memory peak:     11.23 GB
```

* 49384inputs+305736outputs (2major+11220131minor)pagefaults 0swaps


### 03_fullHMP_k31/3.2b_matching_restarted.log
* bwt_loading	3.08s
* sa_loading	5.66s
* bns_loading	23.00s
* matching_time	286.19s
* reads	1000000
* kmers	70000000
* rpm	209654
* kpm	14675803


### 03_fullHMP_k31/3.3a_matching_rolling_skipping.log
* Thu Jul 14 23:57:50 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:05
System time:   0:00:24
Elapsed time:  0:02:30
CPU usage:         100%

Memory peak:     14.12 GB
```

* 0inputs+305736outputs (0major+9228533minor)pagefaults 0swaps


### 03_fullHMP_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	2.49s
* sa_loading	5.58s
* bns_loading	23.11s
* klcp_loading	7.23s
* matching_time	111.26s
* reads	1000000
* kmers	70000000
* rpm	539267
* kpm	37748662


### 03_fullHMP_k31/3.4a_matching_restarted_skipping.log
* Fri Jul 15 00:00:21 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:05:00
System time:   0:00:19
Elapsed time:  0:05:19
CPU usage:         100%

Memory peak:     11.23 GB
```

* 80inputs+305736outputs (1major+11103713minor)pagefaults 0swaps


### 03_fullHMP_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	2.45s
* sa_loading	5.67s
* bns_loading	23.76s
* matching_time	286.92s
* reads	1000000
* kmers	70000000
* rpm	209121
* kpm	14638452


### 03_fullHMP_k31/4.1_contigs_stats.log
* Number of contigs: 33613728
* Total length: 6205913878
* Average length: 184.6243855486663

* Median length: 34.0
* Number of k-mer occurencies: 5197502038


### 03_fullHMP_k31/4.2_index_size.log
* 2.9G	index.fa.31.bit.klcp
* 1.4G	index.fa.ann
* 5.8G	index.fa.bwt
* 2.9G	index.fa.pac
* 2.9G	index.fa.sa

