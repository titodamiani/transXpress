# Kallisto strandedness test on four transXpress runs

Date of the test: 2026-11-03. Runs: piper06, piper09, piper25, piper26.

## Question

Are the RNA-Seq reads of each run stranded, and in which direction? Does the answer match the `strand_specific` value used for the assembly?

| Run | `strand_specific` in `config.yaml` (checked in each file) |
|---|---|
| piper06 | `--ss rf` |
| piper09 | `""` (unstranded) |
| piper25 | `""` (unstranded) |
| piper26 | `""` (unstranded) |

## Methods

Everything below ran on Aurum. Login node work was limited to reading files and `sbatch`.

### Inputs

- Run folders: `/uochb/project/pluskal_group/pipernet/data/interim/transxpress/<species>`. The transcriptome is `results/summary/transcriptome.fasta`. The libraries (name, R1, R2) are in `samples.txt`: 9 libraries for piper06, piper09 and piper25, 12 for piper26.
- Reads: the raw read files named in `samples.txt` (`/uochb/project/pluskal_group/pipernet/data/raw/rnaseq/...`), not the trimmed reads of the pipeline.
- Pipeline master at the time: `4ec0863`. The test did not use the clone or its code.

### Tools

- Image: `/uochb/soft/users/titodamiani/images/transxpress/2e6428dc0cf7788423f25731215fcbde.simg`, sha256 `71fb09f89add02c08dffd67d51a0c97642f6835d0472c51406c8ad943a183bf8` (computed in job 3277141).
- kallisto 0.48.0 at `/conda-envs/52a57574e85d29e6e4d2485c7d666486/bin/kallisto`. This is the version pinned in `workflow/envs/trinity_utils.yaml` (`kallisto=0.48.0=h15996b6_2`).
- seqkit 2.14.0 at `/conda-envs/39f4a4a140ec9685febf193bc07b3778/bin/seqkit` (version printed only). gzip 1.12, mawk 1.3.4 (from the image).
- No package was installed.

### Scripts and jobs

Two scripts ran the test: `run.sh` (Slurm side, starts the image with `apptainer exec --bind /scratch,/uochb/project/pluskal_group/pipernet,<test folder> <image> bash inner.sh <species>`) and `inner.sh` (runs inside the image). Work files lived in `$SCRATCH` of the node (NVMe, not in RAM), which Slurm deletes at job end. The scripts and the scratch folder were deleted after the report. The commands are in the steps below, and the numbers of every library are in the Results tables.

| Run | Job id | Node | sbatch options | Elapsed |
|---|---|---|---|---|
| piper06 | 3277141 | d021 | `--exclude=d001 --constraint=gen-d --time=4:00:00 --mem=32G --cpus-per-task=8` | 00:10:41 |
| piper09 | 3277142 | d021 | same | 00:10:12 |
| piper25 | 3277143 | d021 | same | 00:14:44 |
| piper26 | 3277144 | d021 | same, `--mem=48G` | 00:24:35 |

All four ended COMPLETED. Peak memory (`sacct` MaxRSS of the batch step): 10.8, 9.2, 14.4 and 16.1 GB. Job 3277134 (first try) stopped at the index step because `kallisto index` of 0.48.0 has no `-t` option. I cancelled it, and its output is not used.

Command line of one job:

```bash
sbatch --parsable --exclude=d001 --constraint=gen-d --time=4:00:00 --mem=32G --cpus-per-task=8 -J ks_piper06 \
  -o /uochb/people/titodamiani/projects/kallisto-strand-test/ks_piper06.%j.out \
  /uochb/people/titodamiani/projects/kallisto-strand-test/run.sh piper06
```

### Steps in `inner.sh`

1. `kallisto index -i $SCRATCH/<species>/idx <run>/results/summary/transcriptome.fasta` (default k-mer 31).
2. For each library: take the first 1,000,000 read pairs of R1 and R2 with `gzip -cd <file> | head -n 4000000`. There is no random seed. The same files give the same reads. The script checks that R1 and R2 have 4,000,000 lines each.
3. Run `kallisto quant -i idx -t 8 [flag] -o <dir> r1.fq r2.fq` three times, with no flag, with `--rf-stranded`, and with `--fr-stranded`.
4. Read `n_processed` and `n_pseudoaligned` from `run_info.json` of each quant. Compute the pseudoaligned fraction (`n_pseudoaligned / n_processed`).
5. `paste` the three `abundance.tsv` files (same transcript order, same index) and compute the strand numbers below with awk.

### Numbers in the tables

| Column | Meaning |
|---|---|
| Unstranded, `--rf-stranded`, `--fr-stranded` | Pseudoaligned fraction of the 1,000,000 pairs. |
| rf share | Sum of `est_counts` of the rf run divided by the sum of the rf and fr runs, over all transcripts. |
| Transcripts | Number of transcripts with `est_counts` of 10 or more in the unstranded run. |
| Share with skew 0.9 or more | Among those transcripts, the share where max(rf, fr) / (rf + fr) is 0.9 or more (transcripts with rf + fr = 0 are not counted as skewed). |
| Read-weighted skew | Sum of max(rf, fr) divided by sum of (rf + fr), over those transcripts. |

The cut-offs 10 and 0.9 are my choice. Why the skew numbers exist: an unstranded assembly gives each contig a random direction, so the global rf share is near 0.5 for stranded and for unstranded reads. The per-transcript skew still separates the two cases. A stranded library puts nearly all reads of a transcript in one direction.

### What the sources say

- kallisto manual (pachterlab.github.io/kallisto/manual): `--fr-stranded` processes only fragments where the first read pseudoaligns to the forward strand of a transcript. `--rf-stranded` is the same for the reverse strand. This supports the three settings and the idea that a wrong setting loses reads.
- The test did not use the pipeline `kallisto` rule. It called `kallisto quant` directly with each flag.

## Results

Each table is one run. Fractions are of 1,000,000 read pairs.

### piper06

| Library | Pairs | Unstranded | `--rf-stranded` | `--fr-stranded` | rf share | Transcripts (unstr. count 10 or more) | Share with skew 0.9 or more | Read-weighted skew |
|---|---|---|---|---|---|---|---|---|
| P06_leaf_rep1 | 1000000 | 0.9105 | 0.8991 | 0.0145 | 0.9842 | 18502 | 0.9784 | 0.9909 |
| P06_leaf_rep2 | 1000000 | 0.9118 | 0.9014 | 0.0132 | 0.9855 | 17903 | 0.9817 | 0.9925 |
| P06_leaf_rep3 | 1000000 | 0.9115 | 0.9010 | 0.0134 | 0.9854 | 17698 | 0.9802 | 0.9917 |
| P06_stem_rep1 | 1000000 | 0.9166 | 0.9048 | 0.0155 | 0.9831 | 18099 | 0.9822 | 0.9920 |
| P06_stem_rep2 | 1000000 | 0.9137 | 0.9013 | 0.0170 | 0.9815 | 18743 | 0.9819 | 0.9932 |
| P06_stem_rep3 | 1000000 | 0.9179 | 0.8994 | 0.0231 | 0.9749 | 18452 | 0.9747 | 0.9877 |
| P06_root_rep1 | 1000000 | 0.9113 | 0.8772 | 0.0409 | 0.9554 | 18362 | 0.9780 | 0.9913 |
| P06_root_rep2 | 1000000 | 0.9112 | 0.8812 | 0.0363 | 0.9604 | 18597 | 0.9816 | 0.9922 |
| P06_root_rep3 | 1000000 | 0.9104 | 0.8778 | 0.0403 | 0.9561 | 18664 | 0.9810 | 0.9932 |

Range over 9 libraries: p_rf 0.88-0.90, p_fr 0.01-0.04, rf share 0.955-0.986, skew share 0.975-0.982, weighted skew 0.988-0.993.

### piper09

| Library | Pairs | Unstranded | `--rf-stranded` | `--fr-stranded` | rf share | Transcripts (unstr. count 10 or more) | Share with skew 0.9 or more | Read-weighted skew |
|---|---|---|---|---|---|---|---|---|
| P09_leaf_rep1 | 1000000 | 0.9052 | 0.7313 | 0.7323 | 0.4996 | 14848 | 0.0445 | 0.5661 |
| P09_leaf_rep2 | 1000000 | 0.8988 | 0.7277 | 0.7286 | 0.4997 | 15006 | 0.0522 | 0.5647 |
| P09_leaf_rep3 | 1000000 | 0.9149 | 0.7349 | 0.7359 | 0.4997 | 15204 | 0.0464 | 0.5641 |
| P09_stem_rep1 | 1000000 | 0.9039 | 0.7058 | 0.7073 | 0.4995 | 16115 | 0.0437 | 0.5686 |
| P09_stem_rep2 | 1000000 | 0.9334 | 0.7338 | 0.7358 | 0.4993 | 15527 | 0.0408 | 0.5664 |
| P09_stem_rep3 | 1000000 | 0.9203 | 0.7259 | 0.7270 | 0.4996 | 15414 | 0.0465 | 0.5672 |
| P09_root_rep1 | 1000000 | 0.9411 | 0.7285 | 0.7310 | 0.4991 | 16216 | 0.0459 | 0.5672 |
| P09_root_rep2 | 1000000 | 0.9406 | 0.7506 | 0.7424 | 0.5028 | 16036 | 0.0539 | 0.5857 |
| P09_root_rep3 | 1000000 | 0.9297 | 0.7203 | 0.7228 | 0.4992 | 15387 | 0.0420 | 0.5646 |

Range over 9 libraries: p_rf 0.71-0.75, p_fr 0.71-0.74, rf share 0.499-0.503, skew share 0.041-0.054, weighted skew 0.564-0.586.

### piper25

| Library | Pairs | Unstranded | `--rf-stranded` | `--fr-stranded` | rf share | Transcripts (unstr. count 10 or more) | Share with skew 0.9 or more | Read-weighted skew |
|---|---|---|---|---|---|---|---|---|
| P25_leaf_rep1 | 1000000 | 0.9212 | 0.7557 | 0.7568 | 0.4996 | 17873 | 0.0639 | 0.5751 |
| P25_leaf_rep2 | 1000000 | 0.9244 | 0.7578 | 0.7587 | 0.4997 | 17510 | 0.0655 | 0.5765 |
| P25_leaf_rep3 | 1000000 | 0.9165 | 0.7549 | 0.7555 | 0.4998 | 17663 | 0.0679 | 0.5758 |
| P25_stem_rep1 | 1000000 | 0.9264 | 0.7571 | 0.7591 | 0.4993 | 17110 | 0.0566 | 0.5723 |
| P25_stem_rep2 | 1000000 | 0.9180 | 0.6898 | 0.6905 | 0.4997 | 21066 | 0.0448 | 0.5810 |
| P25_stem_rep3 | 1000000 | 0.9336 | 0.7668 | 0.7681 | 0.4996 | 16604 | 0.0526 | 0.5704 |
| P25_root_rep1 | 1000000 | 0.9270 | 0.7615 | 0.7614 | 0.5000 | 17592 | 0.0575 | 0.5714 |
| P25_root_rep2 | 1000000 | 0.9323 | 0.7600 | 0.7606 | 0.4998 | 17822 | 0.0571 | 0.5739 |
| P25_root_rep3 | 1000000 | 0.9189 | 0.7411 | 0.7431 | 0.4993 | 16526 | 0.0536 | 0.5705 |

Range over 9 libraries: p_rf 0.69-0.77, p_fr 0.69-0.77, rf share 0.499-0.500, skew share 0.045-0.068, weighted skew 0.570-0.581.

### piper26

| Library | Pairs | Unstranded | `--rf-stranded` | `--fr-stranded` | rf share | Transcripts (unstr. count 10 or more) | Share with skew 0.9 or more | Read-weighted skew |
|---|---|---|---|---|---|---|---|---|
| P26_leaf_rep1 | 1000000 | 0.8009 | 0.6540 | 0.6546 | 0.4997 | 19797 | 0.0835 | 0.5956 |
| P26_leaf_rep2 | 1000000 | 0.7959 | 0.6295 | 0.6324 | 0.4988 | 21329 | 0.0847 | 0.6044 |
| P26_leaf_rep3 | 1000000 | 0.7965 | 0.6319 | 0.6337 | 0.4993 | 20925 | 0.0840 | 0.6011 |
| P26_stem_rep1 | 1000000 | 0.7913 | 0.6340 | 0.6352 | 0.4995 | 19542 | 0.0717 | 0.5913 |
| P26_stem_rep2 | 1000000 | 0.7911 | 0.6390 | 0.6408 | 0.4993 | 18822 | 0.0714 | 0.5870 |
| P26_stem_rep3 | 1000000 | 0.8113 | 0.6574 | 0.6581 | 0.4997 | 19986 | 0.0724 | 0.5921 |
| P26_root_rep1 | 1000000 | 0.7976 | 0.6388 | 0.6401 | 0.4995 | 19285 | 0.0762 | 0.5905 |
| P26_root_rep2 | 1000000 | 0.7930 | 0.6321 | 0.6339 | 0.4993 | 19808 | 0.0710 | 0.5930 |
| P26_root_rep3 | 1000000 | 0.7997 | 0.6309 | 0.6319 | 0.4996 | 20499 | 0.0701 | 0.5938 |
| P26_youfr_rep1 | 1000000 | 0.8244 | 0.6623 | 0.6573 | 0.5019 | 19378 | 0.0650 | 0.5935 |
| P26_youfr_rep2 | 1000000 | 0.8203 | 0.6532 | 0.6481 | 0.5020 | 19743 | 0.0643 | 0.5907 |
| P26_youfr_rep3 | 1000000 | 0.8247 | 0.6573 | 0.6535 | 0.5015 | 19780 | 0.0618 | 0.5909 |

Range over 12 libraries: p_rf 0.63-0.66, p_fr 0.63-0.66, rf share 0.499-0.502, skew share 0.062-0.085, weighted skew 0.587-0.604.


## Conclusion per run

| Run | Reads | Direction | `strand_specific` used | Match |
|---|---|---|---|---|
| piper06 | stranded | rf (first read on the reverse strand of the transcript) | `--ss rf` | Yes |
| piper09 | unstranded | none | `""` | Yes |
| piper25 | unstranded | none | `""` | Yes |
| piper26 | unstranded | none | `""` | Yes |

What I saw:
- piper06, all 9 libraries: `--rf-stranded` keeps 0.88 to 0.90 of the pairs (the unstranded run keeps 0.91). `--fr-stranded` keeps 0.01 to 0.04. rf share is 0.955 to 0.986. 0.975 to 0.982 of the transcripts are skewed 0.9 or more.
- piper09, piper25, piper26, all 30 libraries: rf and fr keep the same fraction (differences of 0.01 or less), rf share is 0.499 to 0.503, and 0.04 to 0.085 of the transcripts are skewed 0.9 or more.

What I think:
- piper06 serves as the positive control: the same method gives a clear stranded signal there. The root libraries (rf share 0.955 to 0.960) have more reads in the wrong direction than leaf and stem (0.975 to 0.986). I do not know why.
- For the three unstranded runs, the skew share of 0.04 to 0.085 and the read-weighted skew of 0.56 to 0.60 are the values I saw with no strand signal. I ran no simulated unstranded control, so I cannot give an expected value from first principles.
- In the unstranded runs the rf and fr runs each keep 0.63 to 0.77 of the pairs, more than half of the unstranded fraction. I think the assemblies hold contigs in both directions for many transcripts, so a pair finds a consistent contig in both modes. I did not test this.

## Limits and open questions

- The test uses the first 1,000,000 pairs of each library, not a random sample. Reads at the start of a file may differ from the rest. I do not know if they do.
- The reads are raw. The pipeline uses trimmed reads. I do not expect trimming to change strandedness, but I did not test it.
- The skew cut-offs (count 10, skew 0.9) are my choice. Other values change the percentages, but I did not test how much.
- The test has no check of the library kit. I do not know which kit made the reads. The sequencing facility could say.
- The test shows the direction of the first read relative to the transcript. It says rf for piper06. I did not check the rnaSPAdes manual to confirm that `--ss rf` means the same direction.
- The piper06 assembly used `--ss rf`, so its contigs have a set direction. Its rf share is a check of the assembly direction and of the reads together.
- All libraries of one species gave the same answer. I did not test other species.
- `aurum.md` lists `--ss rf` for the prepared runs except piper25 and piper26. This test did not cover them. Piper09 is unstranded, so other species may be unstranded too. I do not know. A test like this one on each of them would tell.
