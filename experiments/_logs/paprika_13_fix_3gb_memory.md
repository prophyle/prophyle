Table of Contents
=================

* [01_1from10_bacteria_bin_k31](#01_1from10_bacteria_bin_k31)
  * [01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log](#01_1from10_bacteria_bin_k3111_kmer_propagationlog)
  * [01_1from10_bacteria_bin_k31/1.2_merging_fasta.log](#01_1from10_bacteria_bin_k3112_merging_fastalog)
  * [01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log](#01_1from10_bacteria_bin_k3121_bwa_fa2paclog)
  * [01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log](#01_1from10_bacteria_bin_k3122_bwa_pac2bwtgenlog)
  * [01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log](#01_1from10_bacteria_bin_k3123_bwa_bwtupdatelog)
  * [01_1from10_bacteria_bin_k31/2.4_bwa_bwt2sa.log](#01_1from10_bacteria_bin_k3124_bwa_bwt2salog)
  * [01_1from10_bacteria_bin_k31/2.5_klcp.log](#01_1from10_bacteria_bin_k3125_klcplog)
  * [01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log](#01_1from10_bacteria_bin_k3131a_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log](#01_1from10_bacteria_bin_k3131b_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log](#01_1from10_bacteria_bin_k3132a_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log](#01_1from10_bacteria_bin_k3132b_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/3.3a_matching_rolling_skipping.log](#01_1from10_bacteria_bin_k3133a_matching_rolling_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.3b_matching_rolling_skipping.log](#01_1from10_bacteria_bin_k3133b_matching_rolling_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.4a_matching_restarted_skipping.log](#01_1from10_bacteria_bin_k3134a_matching_restarted_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.4b_matching_restarted_skipping.log](#01_1from10_bacteria_bin_k3134b_matching_restarted_skippinglog)
  * [01_1from10_bacteria_bin_k31/4.1_read_assignment.log](#01_1from10_bacteria_bin_k3141_read_assignmentlog)
  * [01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log](#01_1from10_bacteria_bin_k3142_read_assignment_simlcalog)
  * [01_1from10_bacteria_bin_k31/5.1_contigs_stats.log](#01_1from10_bacteria_bin_k3151_contigs_statslog)
  * [01_1from10_bacteria_bin_k31/5.2_index_size.log](#01_1from10_bacteria_bin_k3152_index_sizelog)
* [02_bacteria_orig_k31](#02_bacteria_orig_k31)
  * [02_bacteria_orig_k31/1.1_kmer_propagation.log](#02_bacteria_orig_k3111_kmer_propagationlog)
  * [02_bacteria_orig_k31/1.2_merging_fasta.log](#02_bacteria_orig_k3112_merging_fastalog)
  * [02_bacteria_orig_k31/2.1_bwa_fa2pac.log](#02_bacteria_orig_k3121_bwa_fa2paclog)
  * [02_bacteria_orig_k31/2.2_bwa_pac2bwtgen.log](#02_bacteria_orig_k3122_bwa_pac2bwtgenlog)
  * [02_bacteria_orig_k31/2.3_bwa_bwtupdate.log](#02_bacteria_orig_k3123_bwa_bwtupdatelog)
  * [02_bacteria_orig_k31/2.4_bwa_bwt2sa.log](#02_bacteria_orig_k3124_bwa_bwt2salog)
  * [02_bacteria_orig_k31/2.5_klcp.log](#02_bacteria_orig_k3125_klcplog)
  * [02_bacteria_orig_k31/3.1a_matching_rolling.log](#02_bacteria_orig_k3131a_matching_rollinglog)
  * [02_bacteria_orig_k31/3.1b_matching_rolling.log](#02_bacteria_orig_k3131b_matching_rollinglog)
  * [02_bacteria_orig_k31/3.2a_matching_restarted.log](#02_bacteria_orig_k3132a_matching_restartedlog)
  * [02_bacteria_orig_k31/3.2b_matching_restarted.log](#02_bacteria_orig_k3132b_matching_restartedlog)
  * [02_bacteria_orig_k31/3.3a_matching_rolling_skipping.log](#02_bacteria_orig_k3133a_matching_rolling_skippinglog)
  * [02_bacteria_orig_k31/3.3b_matching_rolling_skipping.log](#02_bacteria_orig_k3133b_matching_rolling_skippinglog)
  * [02_bacteria_orig_k31/3.4a_matching_restarted_skipping.log](#02_bacteria_orig_k3134a_matching_restarted_skippinglog)
  * [02_bacteria_orig_k31/3.4b_matching_restarted_skipping.log](#02_bacteria_orig_k3134b_matching_restarted_skippinglog)
  * [02_bacteria_orig_k31/4.1_read_assignment.log](#02_bacteria_orig_k3141_read_assignmentlog)
  * [02_bacteria_orig_k31/4.2_read_assignment_simlca.log](#02_bacteria_orig_k3142_read_assignment_simlcalog)
  * [02_bacteria_orig_k31/5.1_contigs_stats.log](#02_bacteria_orig_k3151_contigs_statslog)
  * [02_bacteria_orig_k31/5.2_index_size.log](#02_bacteria_orig_k3152_index_sizelog)
* [03_hmp_orig_k31](#03_hmp_orig_k31)
  * [03_hmp_orig_k31/1.1_kmer_propagation.log](#03_hmp_orig_k3111_kmer_propagationlog)
  * [03_hmp_orig_k31/1.2_merging_fasta.log](#03_hmp_orig_k3112_merging_fastalog)
  * [03_hmp_orig_k31/2.1_bwa_fa2pac.log](#03_hmp_orig_k3121_bwa_fa2paclog)
  * [03_hmp_orig_k31/2.2_bwa_pac2bwtgen.log](#03_hmp_orig_k3122_bwa_pac2bwtgenlog)
  * [03_hmp_orig_k31/2.3_bwa_bwtupdate.log](#03_hmp_orig_k3123_bwa_bwtupdatelog)
  * [03_hmp_orig_k31/2.4_bwa_bwt2sa.log](#03_hmp_orig_k3124_bwa_bwt2salog)
  * [03_hmp_orig_k31/2.5_klcp.log](#03_hmp_orig_k3125_klcplog)
  * [03_hmp_orig_k31/3.1a_matching_rolling.log](#03_hmp_orig_k3131a_matching_rollinglog)
  * [03_hmp_orig_k31/3.1b_matching_rolling.log](#03_hmp_orig_k3131b_matching_rollinglog)
  * [03_hmp_orig_k31/3.2a_matching_restarted.log](#03_hmp_orig_k3132a_matching_restartedlog)
  * [03_hmp_orig_k31/3.2b_matching_restarted.log](#03_hmp_orig_k3132b_matching_restartedlog)
  * [03_hmp_orig_k31/3.3a_matching_rolling_skipping.log](#03_hmp_orig_k3133a_matching_rolling_skippinglog)
  * [03_hmp_orig_k31/3.3b_matching_rolling_skipping.log](#03_hmp_orig_k3133b_matching_rolling_skippinglog)
  * [03_hmp_orig_k31/3.4a_matching_restarted_skipping.log](#03_hmp_orig_k3134a_matching_restarted_skippinglog)
  * [03_hmp_orig_k31/3.4b_matching_restarted_skipping.log](#03_hmp_orig_k3134b_matching_restarted_skippinglog)
  * [03_hmp_orig_k31/4.1_read_assignment.log](#03_hmp_orig_k3141_read_assignmentlog)
  * [03_hmp_orig_k31/4.2_read_assignment_simlca.log](#03_hmp_orig_k3142_read_assignment_simlcalog)
  * [03_hmp_orig_k31/5.1_contigs_stats.log](#03_hmp_orig_k3151_contigs_statslog)
  * [03_hmp_orig_k31/5.2_index_size.log](#03_hmp_orig_k3152_index_sizelog)

***
## 01_1from10_bacteria_bin_k31

### 01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log
* Fri Jul 22 19:35:45 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:38:29
System time:   0:01:27
Elapsed time:  0:02:54
CPU usage:        1377%

Memory peak:      1.12 GB
```

* 0inputs+5007784outputs (0major+41134683minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/1.2_merging_fasta.log
* Fri Jul 22 19:38:39 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:58
System time:   0:00:01
Elapsed time:  0:01:01
CPU usage:          98%

Memory peak:      0.01 GB
```

* 0inputs+2443112outputs (0major+4543minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log
* Fri Jul 22 19:39:41 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:12
System time:   0:00:02
Elapsed time:  0:00:22
CPU usage:          64%

Memory peak:      1.32 GB
```

* 64inputs+1741464outputs (0major+555955minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul 22 19:40:03 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     0:13:19
System time:   0:00:00
Elapsed time:  0:13:19
CPU usage:         100%

Memory peak:      0.79 GB
```

* 0inputs+1042696outputs (0major+71377minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log
* Fri Jul 22 19:53:23 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:04
System time:   0:00:02
Elapsed time:  0:00:15
CPU usage:          45%

Memory peak:      1.49 GB
```

* 64inputs+2085384outputs (0major+334022minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.4_bwa_bwt2sa.log
* Fri Jul 22 19:53:39 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:04:17
System time:   0:00:03
Elapsed time:  0:04:26
CPU usage:          98%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+2757392minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.5_klcp.log
* Fri Jul 22 19:58:06 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     0:11:39
System time:   0:00:04
Elapsed time:  0:11:56
CPU usage:          98%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+2932504minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log
* Fri Jul 29 12:27:46 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:23
System time:   0:00:06
Elapsed time:  0:01:41
CPU usage:          88%

Memory peak:      2.50 GB
```

* 4828072inputs+173160outputs (0major+1194426minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log
* bwt_loading	1.23s
* sa_loading	1.46s
* bns_loading	5.66s
* klcp_loading	1.44s
* matching_time	79.35s
* reads	1000000
* kmers	70000000
* rpm	756135
* kpm	52929462


### 01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log
* Fri Jul 29 12:29:27 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:01
System time:   0:00:02
Elapsed time:  0:02:11
CPU usage:          94%

Memory peak:      2.00 GB
```

* 2452432inputs+173160outputs (0major+361795minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log
* bwt_loading	1.39s
* sa_loading	0.45s
* bns_loading	4.69s
* matching_time	116.97s
* reads	1000000
* kmers	70000000
* rpm	512970
* kpm	35907900


### 01_1from10_bacteria_bin_k31/3.3a_matching_rolling_skipping.log
* Fri Jul 29 12:31:39 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:24
System time:   0:00:01
Elapsed time:  0:01:26
CPU usage:         100%

Memory peak:      2.50 GB
```

* 0inputs+173160outputs (0major+290436minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.45s
* sa_loading	0.22s
* bns_loading	4.61s
* klcp_loading	0.23s
* matching_time	80.30s
* reads	1000000
* kmers	70000000
* rpm	747159
* kpm	52301101


### 01_1from10_bacteria_bin_k31/3.4a_matching_restarted_skipping.log
* Fri Jul 29 12:33:05 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:02
System time:   0:00:01
Elapsed time:  0:02:03
CPU usage:         100%

Memory peak:      2.00 GB
```

* 0inputs+173160outputs (0major+327602minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.45s
* sa_loading	0.22s
* bns_loading	4.58s
* matching_time	118.36s
* reads	1000000
* kmers	70000000
* rpm	506937
* kpm	35485613


### 01_1from10_bacteria_bin_k31/4.1_read_assignment.log
* Fri Jul 29 12:35:09 CEST 2016
* jobs: 1
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a

```
User time:     0:00:29
System time:   0:00:00
Elapsed time:  0:00:30
CPU usage:          97%

Memory peak:      0.05 GB
```

* 19272inputs+0outputs (52major+20249minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log
* Fri Jul 29 12:35:39 CEST 2016
* jobs: 1
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a -t

```
User time:     0:00:44
System time:   0:00:00
Elapsed time:  0:00:44
CPU usage:         100%

Memory peak:      0.05 GB
```

* 0inputs+0outputs (0major+22426minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/5.1_contigs_stats.log
* Number of contigs: 8495372
* Total length: 1067712578
* Average length: 125.68167444580413

* Median length: 34.0
* Number of k-mer occurencies: 812851418


### 01_1from10_bacteria_bin_k31/5.2_index_size.log
* 510M	index.fa.31.bit.klcp
* 342M	index.fa.ann
* 1019M	index.fa.bwt
* 510M	index.fa.pac
* 510M	index.fa.sa

***
## 02_bacteria_orig_k31

### 02_bacteria_orig_k31/1.1_kmer_propagation.log
* Fri Jul 22 20:17:09 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     4:59:08
System time:   0:16:03
Elapsed time:  0:22:41
CPU usage:        1389%

Memory peak:      9.31 GB
```

* 32611512inputs+43545296outputs (46major+376635908minor)pagefaults 0swaps


### 02_bacteria_orig_k31/1.2_merging_fasta.log
* Fri Jul 22 20:39:50 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:12:08
System time:   0:00:27
Elapsed time:  0:17:37
CPU usage:          71%

Memory peak:      0.01 GB
```

* 16204336inputs+23879936outputs (21major+27143minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.1_bwa_fa2pac.log
* Fri Jul 22 20:57:29 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:09
System time:   0:00:40
Elapsed time:  0:04:43
CPU usage:          60%

Memory peak:     13.94 GB
```

* 1688208inputs+18043024outputs (5major+9165965minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul 22 21:02:12 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:42:13
System time:   0:00:07
Elapsed time:  2:42:32
CPU usage:         100%

Memory peak:      5.54 GB
```

* 2145672inputs+9887296outputs (2major+422671minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.3_bwa_bwtupdate.log
* Fri Jul 22 23:44:45 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:28
Elapsed time:  0:02:43
CPU usage:          44%

Memory peak:     14.15 GB
```

* 360inputs+19774584outputs (0major+3352956minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.4_bwa_bwt2sa.log
* Fri Jul 22 23:47:28 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:50:12
System time:   0:01:44
Elapsed time:  0:52:49
CPU usage:          98%

Memory peak:     14.15 GB
```

* 840inputs+9887296outputs (4major+125206730minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.5_klcp.log
* Sat Jul 23 00:40:17 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     1:59:46
System time:   0:03:38
Elapsed time:  2:04:34
CPU usage:          99%

Memory peak:     14.15 GB
```

* 5372496inputs+9887296outputs (3major+187924523minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1a_matching_rolling.log
* Fri Jul 29 12:36:24 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:17
System time:   0:01:13
Elapsed time:  0:05:38
CPU usage:          80%

Memory peak:     23.00 GB
```

* 30041464inputs+360824outputs (0major+20751283minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	25.31s
* sa_loading	16.57s
* bns_loading	61.22s
* klcp_loading	12.39s
* matching_time	153.91s
* reads	1000000
* kmers	70000000
* rpm	389834
* kpm	27288369


### 02_bacteria_orig_k31/3.2a_matching_restarted.log
* Fri Jul 29 12:42:02 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:08:48
System time:   0:01:09
Elapsed time:  0:10:59
CPU usage:          91%

Memory peak:     18.29 GB
```

* 27085616inputs+360824outputs (1major+36165728minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	25.85s
* sa_loading	15.78s
* bns_loading	60.14s
* matching_time	494.63s
* reads	1000000
* kmers	70000000
* rpm	121302
* kpm	8491135


### 02_bacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Fri Jul 29 12:53:01 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:19
System time:   0:00:56
Elapsed time:  0:05:01
CPU usage:          85%

Memory peak:     23.00 GB
```

* 12100360inputs+360824outputs (2major+16821263minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	18.52s
* sa_loading	11.45s
* bns_loading	59.00s
* klcp_loading	10.27s
* matching_time	154.92s
* reads	1000000
* kmers	70000000
* rpm	387294
* kpm	27110570


### 02_bacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Fri Jul 29 12:58:03 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:08:47
System time:   0:01:05
Elapsed time:  0:10:17
CPU usage:          96%

Memory peak:     18.29 GB
```

* 7971728inputs+360824outputs (1major+45899156minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	21.36s
* sa_loading	11.48s
* bns_loading	58.85s
* matching_time	499.28s
* reads	1000000
* kmers	70000000
* rpm	120173
* kpm	8412088


### 02_bacteria_orig_k31/4.1_read_assignment.log
* Fri Jul 29 13:08:20 CEST 2016
* jobs: 1
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:02:07
System time:   0:00:00
Elapsed time:  0:02:08
CPU usage:          99%

Memory peak:      0.07 GB
```

* 407840inputs+0outputs (180major+122733minor)pagefaults 0swaps


### 02_bacteria_orig_k31/4.2_read_assignment_simlca.log
* Fri Jul 29 13:10:29 CEST 2016
* jobs: 1
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:04:39
System time:   0:00:00
Elapsed time:  0:04:40
CPU usage:         100%

Memory peak:      0.07 GB
```

* 0inputs+0outputs (0major+64637minor)pagefaults 0swaps


### 02_bacteria_orig_k31/5.1_contigs_stats.log
* Number of contigs: 95208424
* Total length: 10124583056
* Average length: 106.34125249253154

* Median length: 34.0
* Number of k-mer occurencies: 7268330336


### 02_bacteria_orig_k31/5.2_index_size.log
* 4.8G	index.fa.31.bit.klcp
* 3.9G	index.fa.ann
* 9.5G	index.fa.bwt
* 4.8G	index.fa.pac
* 4.8G	index.fa.sa

***
## 03_hmp_orig_k31

### 03_hmp_orig_k31/1.1_kmer_propagation.log
* Sat Jul 23 03:15:37 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     2:53:08
System time:   0:04:24
Elapsed time:  0:12:04
CPU usage:        1470%

Memory peak:      5.81 GB
```

* 11552456inputs+26263224outputs (9major+147070940minor)pagefaults 0swaps


### 03_hmp_orig_k31/1.2_merging_fasta.log
* Sat Jul 23 03:27:42 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:05:27
System time:   0:00:12
Elapsed time:  0:05:57
CPU usage:          95%

Memory peak:      0.01 GB
```

* 56inputs+14285504outputs (2major+25802minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.1_bwa_fa2pac.log
* Sat Jul 23 03:33:39 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:07
System time:   0:00:11
Elapsed time:  0:02:18
CPU usage:          57%

Memory peak:      7.06 GB
```

* 1104inputs+9814544outputs (3major+2687492minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.2_bwa_pac2bwtgen.log
* Sat Jul 23 03:35:57 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     1:31:47
System time:   0:00:03
Elapsed time:  1:31:51
CPU usage:         100%

Memory peak:      3.56 GB
```

* 0inputs+6201448outputs (0major+223255minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.3_bwa_bwtupdate.log
* Sat Jul 23 05:07:49 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:25
System time:   0:00:17
Elapsed time:  0:01:52
CPU usage:          38%

Memory peak:      8.87 GB
```

* 272inputs+12402896outputs (0major+2441202minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.4_bwa_bwt2sa.log
* Sat Jul 23 05:09:42 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:32:54
System time:   0:01:01
Elapsed time:  0:34:30
CPU usage:          98%

Memory peak:      8.87 GB
```

* 168inputs+6201448outputs (0major+68580216minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.5_klcp.log
* Sat Jul 23 05:44:12 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     1:10:18
System time:   0:01:17
Elapsed time:  1:12:10
CPU usage:          99%

Memory peak:      8.87 GB
```

* 304inputs+6201448outputs (1major+106897669minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1a_matching_rolling.log
* Fri Jul 29 13:15:10 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:34
System time:   0:00:52
Elapsed time:  0:04:39
CPU usage:          74%

Memory peak:     13.75 GB
```

* 28127488inputs+258864outputs (0major+24883301minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1b_matching_rolling.log
* bwt_loading	15.51s
* sa_loading	10.32s
* bns_loading	27.05s
* klcp_loading	10.56s
* matching_time	142.40s
* reads	1000000
* kmers	70000000
* rpm	421349
* kpm	29494412


### 03_hmp_orig_k31/3.2a_matching_restarted.log
* Fri Jul 29 13:19:50 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:51
System time:   0:00:11
Elapsed time:  0:05:35
CPU usage:          90%

Memory peak:     10.80 GB
```

* 8493968inputs+258864outputs (1major+418263minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.2b_matching_restarted.log
* bwt_loading	7.20s
* sa_loading	2.14s
* bns_loading	23.78s
* matching_time	269.49s
* reads	1000000
* kmers	70000000
* rpm	222644
* kpm	15585046


### 03_hmp_orig_k31/3.3a_matching_rolling_skipping.log
* Fri Jul 29 13:25:25 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:32
System time:   0:00:39
Elapsed time:  0:03:11
CPU usage:         100%

Memory peak:     13.75 GB
```

* 0inputs+258864outputs (0major+24877865minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	9.35s
* sa_loading	6.66s
* bns_loading	26.68s
* klcp_loading	6.31s
* matching_time	141.68s
* reads	1000000
* kmers	70000000
* rpm	423504
* kpm	29645266


### 03_hmp_orig_k31/3.4a_matching_restarted_skipping.log
* Fri Jul 29 13:28:37 CEST 2016
* jobs: 1
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:50
System time:   0:00:06
Elapsed time:  0:04:56
CPU usage:         100%

Memory peak:     10.80 GB
```

* 0inputs+258864outputs (0major+462039minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	2.50s
* sa_loading	1.31s
* bns_loading	23.59s
* matching_time	268.55s
* reads	1000000
* kmers	70000000
* rpm	223419
* kpm	15639315


### 03_hmp_orig_k31/4.1_read_assignment.log
* Fri Jul 29 13:33:34 CEST 2016
* jobs: 1
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a

```
User time:     0:01:15
System time:   0:00:00
Elapsed time:  0:01:17
CPU usage:          99%

Memory peak:      0.09 GB
```

* 270632inputs+0outputs (10major+180898minor)pagefaults 0swaps


### 03_hmp_orig_k31/4.2_read_assignment_simlca.log
* Fri Jul 29 13:34:51 CEST 2016
* jobs: 1
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a -t

```
User time:     0:02:30
System time:   0:00:00
Elapsed time:  0:02:31
CPU usage:         100%

Memory peak:      0.09 GB
```

* 0inputs+0outputs (0major+137040minor)pagefaults 0swaps


### 03_hmp_orig_k31/5.1_contigs_stats.log
* Number of contigs: 42285628
* Total length: 6350278678
* Average length: 150.17581571686722

* Median length: 34.0
* Number of k-mer occurencies: 5081709838


### 03_hmp_orig_k31/5.2_index_size.log
* 3.0G	index.fa.31.bit.klcp
* 1.8G	index.fa.ann
* 6.0G	index.fa.bwt
* 3.0G	index.fa.pac
* 3.0G	index.fa.sa

