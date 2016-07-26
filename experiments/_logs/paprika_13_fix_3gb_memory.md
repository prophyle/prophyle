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
* Tue Jul 26 12:06:30 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:23
System time:   0:00:03
Elapsed time:  0:01:39
CPU usage:          87%

Memory peak:      2.50 GB
```

* 4869568inputs+173160outputs (0major+240306minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log
* bwt_loading	1.23s
* sa_loading	0.61s
* bns_loading	5.14s
* klcp_loading	0.60s
* matching_time	78.73s
* reads	1000000
* kmers	70000000
* rpm	762140
* kpm	53349826


### 01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log
* Tue Jul 26 12:08:09 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:07
System time:   0:00:02
Elapsed time:  0:02:10
CPU usage:         100%

Memory peak:      2.00 GB
```

* 0inputs+173160outputs (0major+508694minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.47s
* sa_loading	0.21s
* bns_loading	5.10s
* matching_time	124.65s
* reads	1000000
* kmers	70000000
* rpm	481330
* kpm	33693086


### 01_1from10_bacteria_bin_k31/3.3a_matching_rolling_skipping.log
* Tue Jul 26 12:10:20 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:25
System time:   0:00:03
Elapsed time:  0:01:29
CPU usage:         100%

Memory peak:      2.50 GB
```

* 0inputs+173160outputs (0major+1249828minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.49s
* sa_loading	0.24s
* bns_loading	5.31s
* klcp_loading	0.91s
* matching_time	81.86s
* reads	1000000
* kmers	70000000
* rpm	732985
* kpm	51308965


### 01_1from10_bacteria_bin_k31/3.4a_matching_restarted_skipping.log
* Tue Jul 26 12:11:49 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:00
System time:   0:00:01
Elapsed time:  0:02:02
CPU usage:         100%

Memory peak:      2.00 GB
```

* 0inputs+173160outputs (0major+546494minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.46s
* sa_loading	0.21s
* bns_loading	5.09s
* matching_time	116.08s
* reads	1000000
* kmers	70000000
* rpm	516900
* kpm	36183018


### 01_1from10_bacteria_bin_k31/4.1_read_assignment.log
* Tue Jul 26 12:08:09 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a

```
User time:     0:00:28
System time:   0:00:00
Elapsed time:  0:00:29
CPU usage:         100%

Memory peak:      0.06 GB
```

* 488inputs+0outputs (0major+23891minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log
* Tue Jul 26 12:08:09 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a -t

```
User time:     0:00:45
System time:   0:00:00
Elapsed time:  0:00:45
CPU usage:         100%

Memory peak:      0.06 GB
```

* 0inputs+0outputs (0major+30651minor)pagefaults 0swaps


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
* Tue Jul 26 12:13:52 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:23
System time:   0:01:35
Elapsed time:  0:06:40
CPU usage:          75%

Memory peak:     23.00 GB
```

* 46029368inputs+360824outputs (0major+48674134minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	30.04s
* sa_loading	16.09s
* bns_loading	64.80s
* klcp_loading	16.11s
* matching_time	169.24s
* reads	1000000
* kmers	70000000
* rpm	354522
* kpm	24816563


### 02_bacteria_orig_k31/3.2a_matching_restarted.log
* Tue Jul 26 12:20:32 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:08:55
System time:   0:01:19
Elapsed time:  0:11:53
CPU usage:          86%

Memory peak:     18.29 GB
```

* 28716656inputs+360824outputs (0major+79908497minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	17.88s
* sa_loading	11.60s
* bns_loading	63.68s
* matching_time	520.55s
* reads	1000000
* kmers	70000000
* rpm	115263
* kpm	8068432


### 02_bacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Tue Jul 26 12:32:26 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:18
System time:   0:01:13
Elapsed time:  0:05:41
CPU usage:          80%

Memory peak:     23.00 GB
```

* 23697992inputs+360824outputs (0major+34765810minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	21.54s
* sa_loading	8.39s
* bns_loading	62.53s
* klcp_loading	17.77s
* matching_time	159.50s
* reads	1000000
* kmers	70000000
* rpm	376169
* kpm	26331848


### 02_bacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Tue Jul 26 12:38:08 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:07
System time:   0:01:23
Elapsed time:  0:10:55
CPU usage:          96%

Memory peak:     18.29 GB
```

* 7015528inputs+360824outputs (2major+86710194minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	19.40s
* sa_loading	10.28s
* bns_loading	63.58s
* matching_time	535.98s
* reads	1000000
* kmers	70000000
* rpm	111945
* kpm	7836141


### 02_bacteria_orig_k31/4.1_read_assignment.log
* Tue Jul 26 12:20:32 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:02:06
System time:   0:00:01
Elapsed time:  0:02:19
CPU usage:          92%

Memory peak:      0.07 GB
```

* 104808inputs+0outputs (214major+88863minor)pagefaults 0swaps


### 02_bacteria_orig_k31/4.2_read_assignment_simlca.log
* Tue Jul 26 12:20:32 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:04:47
System time:   0:00:00
Elapsed time:  0:04:57
CPU usage:          97%

Memory peak:      0.07 GB
```

* 31032inputs+0outputs (191major+51308minor)pagefaults 0swaps


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
* Tue Jul 26 12:49:04 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:33
System time:   0:00:46
Elapsed time:  0:04:33
CPU usage:          73%

Memory peak:     13.75 GB
```

* 28418912inputs+258864outputs (0major+18635182minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1b_matching_rolling.log
* bwt_loading	13.03s
* sa_loading	10.00s
* bns_loading	28.39s
* klcp_loading	10.24s
* matching_time	136.66s
* reads	1000000
* kmers	70000000
* rpm	439051
* kpm	30733580


### 03_hmp_orig_k31/3.2a_matching_restarted.log
* Tue Jul 26 12:53:37 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:05:36
System time:   0:00:20
Elapsed time:  0:06:40
CPU usage:          89%

Memory peak:     10.80 GB
```

* 9029976inputs+258864outputs (1major+10541896minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.2b_matching_restarted.log
* bwt_loading	6.08s
* sa_loading	1.60s
* bns_loading	29.73s
* matching_time	317.73s
* reads	1000000
* kmers	70000000
* rpm	188842
* kpm	13218935


### 03_hmp_orig_k31/3.3a_matching_rolling_skipping.log
* Tue Jul 26 13:00:17 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:30
System time:   0:00:41
Elapsed time:  0:03:28
CPU usage:          92%

Memory peak:     13.75 GB
```

* 6224872inputs+258864outputs (0major+23192518minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	6.90s
* sa_loading	6.60s
* bns_loading	28.19s
* klcp_loading	11.32s
* matching_time	137.88s
* reads	1000000
* kmers	70000000
* rpm	435171
* kpm	30461961


### 03_hmp_orig_k31/3.4a_matching_restarted_skipping.log
* Tue Jul 26 13:03:46 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:53
System time:   0:00:29
Elapsed time:  0:05:23
CPU usage:         100%

Memory peak:     10.80 GB
```

* 56inputs+258864outputs (0major+19892953minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	7.29s
* sa_loading	7.64s
* bns_loading	27.92s
* matching_time	279.55s
* reads	1000000
* kmers	70000000
* rpm	214633
* kpm	15024304


### 03_hmp_orig_k31/4.1_read_assignment.log
* Tue Jul 26 12:53:37 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a

```
User time:     0:01:13
System time:   0:00:00
Elapsed time:  0:01:24
CPU usage:          88%

Memory peak:      0.09 GB
```

* 35336inputs+0outputs (165major+123213minor)pagefaults 0swaps


### 03_hmp_orig_k31/4.2_read_assignment_simlca.log
* Tue Jul 26 12:53:37 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a -t

```
User time:     0:02:43
System time:   0:00:00
Elapsed time:  0:02:53
CPU usage:          95%

Memory peak:      0.09 GB
```

* 36960inputs+0outputs (165major+215060minor)pagefaults 0swaps


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

