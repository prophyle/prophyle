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
* [03_fullBacteria_orig_k31_lowcomp_masked](#03_fullBacteria_orig_k31_lowcomp_masked)
  * [03_fullBacteria_orig_k31_lowcomp_masked/1.1_kmer_propagation.log](#03_fullbacteria_orig_k31_lowcomp_masked11_kmer_propagationlog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/1.2_merging_fasta.log](#03_fullbacteria_orig_k31_lowcomp_masked12_merging_fastalog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/2.1_bwa_fa2pac.log](#03_fullbacteria_orig_k31_lowcomp_masked21_bwa_fa2paclog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/2.2_bwa_pac2bwtgen.log](#03_fullbacteria_orig_k31_lowcomp_masked22_bwa_pac2bwtgenlog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/2.3_bwa_bwtupdate.log](#03_fullbacteria_orig_k31_lowcomp_masked23_bwa_bwtupdatelog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/2.4_bwa_bwt2sa.log](#03_fullbacteria_orig_k31_lowcomp_masked24_bwa_bwt2salog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/2.5_klcp.log](#03_fullbacteria_orig_k31_lowcomp_masked25_klcplog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.1a_matching_rolling.log](#03_fullbacteria_orig_k31_lowcomp_masked31a_matching_rollinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.1b_matching_rolling.log](#03_fullbacteria_orig_k31_lowcomp_masked31b_matching_rollinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.2a_matching_restarted.log](#03_fullbacteria_orig_k31_lowcomp_masked32a_matching_restartedlog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.2b_matching_restarted.log](#03_fullbacteria_orig_k31_lowcomp_masked32b_matching_restartedlog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.3a_matching_rolling_skipping.log](#03_fullbacteria_orig_k31_lowcomp_masked33a_matching_rolling_skippinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.3b_matching_rolling_skipping.log](#03_fullbacteria_orig_k31_lowcomp_masked33b_matching_rolling_skippinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.4a_matching_restarted_skipping.log](#03_fullbacteria_orig_k31_lowcomp_masked34a_matching_restarted_skippinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/3.4b_matching_restarted_skipping.log](#03_fullbacteria_orig_k31_lowcomp_masked34b_matching_restarted_skippinglog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/4.1_contigs_stats.log](#03_fullbacteria_orig_k31_lowcomp_masked41_contigs_statslog)
  * [03_fullBacteria_orig_k31_lowcomp_masked/4.2_index_size.log](#03_fullbacteria_orig_k31_lowcomp_masked42_index_sizelog)

***
## 01_1from10_bin_k31

### 01_1from10_bin_k31/1.1_kmer_propagation.log
* Tue Jul 12 00:10:24 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:37:46
System time:   0:00:56
Elapsed time:  0:02:33
CPU usage:        1516%

Memory peak:      1.12 GB
```

* 743312inputs+5007784outputs (1major+35153638minor)pagefaults 0swaps


### 01_1from10_bin_k31/1.2_merging_fasta.log
* Tue Jul 12 00:12:57 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:55
System time:   0:00:02
Elapsed time:  0:00:59
CPU usage:          98%

Memory peak:      0.01 GB
```

* 16inputs+2539240outputs (0major+7497minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.1_bwa_fa2pac.log
* Tue Jul 12 00:13:57 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:12
System time:   0:00:01
Elapsed time:  0:00:20
CPU usage:          70%

Memory peak:      1.36 GB
```

* 56inputs+1836576outputs (0major+471507minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.2_bwa_pac2bwtgen.log
* Tue Jul 12 00:14:17 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     0:12:50
System time:   0:00:00
Elapsed time:  0:12:51
CPU usage:         100%

Memory peak:      0.79 GB
```

* 0inputs+1042696outputs (0major+21932minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.3_bwa_bwtupdate.log
* Tue Jul 12 00:27:08 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:04
System time:   0:00:01
Elapsed time:  0:00:15
CPU usage:          39%

Memory peak:      1.49 GB
```

* 80inputs+2085384outputs (0major+2689minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.4_bwa_bwt2sa.log
* Tue Jul 12 00:27:24 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:04:08
System time:   0:00:01
Elapsed time:  0:04:14
CPU usage:          98%

Memory peak:      1.49 GB
```

* 32inputs+1042696outputs (0major+28205minor)pagefaults 0swaps


### 01_1from10_bin_k31/2.5_klcp.log
* Tue Jul 12 00:31:38 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     0:11:29
System time:   0:00:01
Elapsed time:  0:11:35
CPU usage:          99%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+47034minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1a_matching_rolling.log
* Tue Jul 12 19:26:55 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:33
System time:   0:00:10
Elapsed time:  0:01:53
CPU usage:          92%

Memory peak:      2.47 GB
```

* 4991376inputs+873456outputs (1major+4278316minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.1b_matching_rolling.log
* bwt_loading	2.56s
* sa_loading	1.77s
* bns_loading	5.43s
* klcp_loading	1.72s
* matching_time	92.15s
* reads	1000000
* kmers	70000000
* rpm	651117
* kpm	45578209


### 01_1from10_bin_k31/3.2a_matching_restarted.log
* Tue Jul 12 19:28:49 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:05
System time:   0:00:03
Elapsed time:  0:02:11
CPU usage:          98%

Memory peak:      1.97 GB
```

* 0inputs+873456outputs (0major+1164512minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.89s
* sa_loading	0.55s
* bns_loading	5.11s
* matching_time	122.50s
* reads	1000000
* kmers	70000000
* rpm	489800
* kpm	34285971


### 01_1from10_bin_k31/3.3a_matching_rolling_skipping.log
* Tue Jul 12 19:31:00 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:27
System time:   0:00:01
Elapsed time:  0:01:30
CPU usage:          98%

Memory peak:      2.47 GB
```

* 0inputs+872584outputs (0major+275896minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.42s
* sa_loading	0.19s
* bns_loading	4.56s
* klcp_loading	0.19s
* matching_time	83.51s
* reads	1000000
* kmers	70000000
* rpm	718516
* kpm	50296106


### 01_1from10_bin_k31/3.4a_matching_restarted_skipping.log
* Tue Jul 12 19:32:31 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:04
System time:   0:00:01
Elapsed time:  0:02:07
CPU usage:          99%

Memory peak:      1.97 GB
```

* 16inputs+872584outputs (0major+298415minor)pagefaults 0swaps


### 01_1from10_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.40s
* sa_loading	0.19s
* bns_loading	4.74s
* matching_time	120.49s
* reads	1000000
* kmers	70000000
* rpm	497976
* kpm	34858299


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
* Tue Jul 12 00:50:47 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     4:52:05
System time:   0:07:35
Elapsed time:  0:16:12
CPU usage:        1849%

Memory peak:      9.31 GB
```

* 11288072inputs+44018032outputs (7major+233868317minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/1.2_merging_fasta.log
* Tue Jul 12 01:06:59 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:09:48
System time:   0:00:22
Elapsed time:  0:11:58
CPU usage:          85%

Memory peak:      0.01 GB
```

* 14234168inputs+24270192outputs (21major+22878minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.1_bwa_fa2pac.log
* Tue Jul 12 01:18:58 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:06
System time:   0:00:26
Elapsed time:  0:03:40
CPU usage:          69%

Memory peak:     14.08 GB
```

* 1288inputs+18388752outputs (4major+8170067minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Tue Jul 12 01:22:38 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:24:41
System time:   0:00:06
Elapsed time:  2:24:49
CPU usage:         100%

Memory peak:      5.58 GB
```

* 166088inputs+9968312outputs (0major+365927minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.3_bwa_bwtupdate.log
* Tue Jul 12 03:47:27 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:34
Elapsed time:  0:02:34
CPU usage:          50%

Memory peak:     14.26 GB
```

* 432inputs+19936624outputs (0major+4623366minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.4_bwa_bwt2sa.log
* Tue Jul 12 03:50:02 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     1:01:19
System time:   0:03:01
Elapsed time:  1:05:03
CPU usage:          99%

Memory peak:     14.26 GB
```

* 240inputs+9968312outputs (0major+210317194minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/2.5_klcp.log
* Tue Jul 12 04:55:05 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     2:07:38
System time:   0:04:23
Elapsed time:  2:12:40
CPU usage:         100%

Memory peak:     14.26 GB
```

* 1248inputs+9968312outputs (2major+295534353minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1a_matching_rolling.log
* Tue Jul 12 19:34:39 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:05
System time:   0:01:24
Elapsed time:  0:07:14
CPU usage:          76%

Memory peak:     22.83 GB
```

* 48293728inputs+3943496outputs (0major+23761767minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	28.53s
* sa_loading	16.00s
* bns_loading	62.57s
* klcp_loading	16.41s
* matching_time	204.01s
* reads	1000000
* kmers	70000000
* rpm	294109
* kpm	20587619


### 02_fullBacteria_orig_k31/3.2a_matching_restarted.log
* Tue Jul 12 19:41:53 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:41
System time:   0:00:51
Elapsed time:  0:11:47
CPU usage:          89%

Memory peak:     18.08 GB
```

* 37691824inputs+3943496outputs (0major+19193354minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	14.04s
* sa_loading	14.51s
* bns_loading	60.85s
* matching_time	541.79s
* reads	1000000
* kmers	70000000
* rpm	110743
* kpm	7752024


### 02_fullBacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Tue Jul 12 19:53:41 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:51
System time:   0:04:30
Elapsed time:  0:09:41
CPU usage:          86%

Memory peak:     22.83 GB
```

* 47305960inputs+3942912outputs (1major+25292199minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	11.45s
* sa_loading	183.78s
* bns_loading	93.66s
* klcp_loading	15.61s
* matching_time	195.46s
* reads	1000000
* kmers	70000000
* rpm	306961
* kpm	21487242


### 02_fullBacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Tue Jul 12 20:03:23 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:44
System time:   0:00:54
Elapsed time:  0:11:38
CPU usage:          91%

Memory peak:     18.08 GB
```

* 31081232inputs+3942912outputs (0major+32465928minor)pagefaults 0swaps


### 02_fullBacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	11.96s
* sa_loading	13.38s
* bns_loading	59.59s
* matching_time	551.74s
* reads	1000000
* kmers	70000000
* rpm	108746
* kpm	7612232


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
## 03_fullBacteria_orig_k31_lowcomp_masked

### 03_fullBacteria_orig_k31_lowcomp_masked/1.1_kmer_propagation.log
* Tue Jul 12 07:40:48 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     5:02:39
System time:   0:09:03
Elapsed time:  0:19:07
CPU usage:        1630%

Memory peak:      9.22 GB
```

* 25637320inputs+44984216outputs (675major+248045896minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/1.2_merging_fasta.log
* Tue Jul 12 07:59:56 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:10:57
System time:   0:00:23
Elapsed time:  0:13:13
CPU usage:          86%

Memory peak:      0.01 GB
```

* 14475664inputs+25364072outputs (6major+72741minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/2.1_bwa_fa2pac.log
* Tue Jul 12 08:13:09 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:17
System time:   0:00:34
Elapsed time:  0:04:17
CPU usage:          67%

Memory peak:     15.71 GB
```

* 2822792inputs+19995136outputs (3major+9046713minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/2.2_bwa_pac2bwtgen.log
* Tue Jul 12 08:17:27 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:29:54
System time:   0:00:07
Elapsed time:  2:30:10
CPU usage:         100%

Memory peak:      5.70 GB
```

* 2194312inputs+10185336outputs (3major+374271minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/2.3_bwa_bwtupdate.log
* Tue Jul 12 10:47:37 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:30
Elapsed time:  0:02:34
CPU usage:          48%

Memory peak:     14.57 GB
```

* 456inputs+20370672outputs (0major+3381363minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/2.4_bwa_bwt2sa.log
* Tue Jul 12 10:50:11 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:54:19
System time:   0:02:47
Elapsed time:  0:57:49
CPU usage:          99%

Memory peak:     14.57 GB
```

* 393296inputs+10185336outputs (1major+146521481minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/2.5_klcp.log
* Tue Jul 12 11:48:01 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     2:10:37
System time:   0:04:29
Elapsed time:  2:15:57
CPU usage:          99%

Memory peak:     14.57 GB
```

* 1176inputs+10185336outputs (2major+264031250minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/3.1a_matching_rolling.log
* Tue Jul 12 20:15:02 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:27
System time:   0:01:10
Elapsed time:  0:07:17
CPU usage:          77%

Memory peak:     23.85 GB
```

* 50211000inputs+3927312outputs (0major+16602406minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/3.1b_matching_rolling.log
* bwt_loading	12.22s
* sa_loading	13.45s
* bns_loading	68.97s
* klcp_loading	15.51s
* matching_time	225.69s
* reads	1000000
* kmers	70000000
* rpm	265855
* kpm	18609866


### 03_fullBacteria_orig_k31_lowcomp_masked/3.2a_matching_restarted.log
* Tue Jul 12 20:22:20 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 -v index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:47
System time:   0:00:55
Elapsed time:  0:11:09
CPU usage:          96%

Memory peak:     18.99 GB
```

* 16883512inputs+3927312outputs (6major+33110977minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/3.2b_matching_restarted.log
* bwt_loading	6.96s
* sa_loading	14.39s
* bns_loading	71.52s
* matching_time	547.81s
* reads	1000000
* kmers	70000000
* rpm	109527
* kpm	7666883


### 03_fullBacteria_orig_k31_lowcomp_masked/3.3a_matching_rolling_skipping.log
* Tue Jul 12 20:33:29 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:10
System time:   0:01:03
Elapsed time:  0:06:03
CPU usage:          86%

Memory peak:     23.85 GB
```

* 16191096inputs+3925608outputs (5major+34558314minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/3.3b_matching_rolling_skipping.log
* bwt_loading	11.55s
* sa_loading	10.10s
* bns_loading	68.37s
* klcp_loading	10.55s
* matching_time	211.30s
* reads	1000000
* kmers	70000000
* rpm	283962
* kpm	19877331


### 03_fullBacteria_orig_k31_lowcomp_masked/3.4a_matching_restarted_skipping.log
* Tue Jul 12 20:39:33 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -v -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:10:38
System time:   0:01:06
Elapsed time:  0:12:50
CPU usage:          91%

Memory peak:     18.99 GB
```

* 26574440inputs+3925608outputs (3major+40363720minor)pagefaults 0swaps


### 03_fullBacteria_orig_k31_lowcomp_masked/3.4b_matching_restarted_skipping.log
* bwt_loading	15.47s
* sa_loading	15.13s
* bns_loading	68.47s
* matching_time	603.63s
* reads	1000000
* kmers	70000000
* rpm	99399
* kpm	6957959


### 03_fullBacteria_orig_k31_lowcomp_masked/4.1_contigs_stats.log
* Number of contigs: 111846835
* Total length: 10429783705
* Average length: 93.25059314373983

* Median length: 35
* Number of k-mer occurencies: 7074378655


### 03_fullBacteria_orig_k31_lowcomp_masked/4.2_index_size.log
* 4.9G	index.fa.31.bit.klcp
* 4.7G	index.fa.ann
* 9.8G	index.fa.bwt
* 4.9G	index.fa.pac
* 4.9G	index.fa.sa

