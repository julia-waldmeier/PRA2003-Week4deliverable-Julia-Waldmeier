# PRA2003 — Bacterial Population Analysis (Biology Track)

Julia Waldmeier

## 1. Overview

This project analyzes simulated experimental data representing bacterial
populations across millions of independent "events" (runs). Each event
contains a variable number of bacteria, each tagged with an ID code
identifying its strain (see table below). The goal is to determine, for
each of 12 bacterial strains, the **average number of that strain observed
per event**, together with its **statistical uncertainty**, using the full
5,000,000-event dataset.

## 2. Data

- 10 data files (`output-Set1.txt` ... `output-Set10.txt`), ~500,000 events
  each, ~800MB per file — available from the course surfdrive, not included
  in this repository due to size.
- Each file is plain text, structured as repeating blocks:
  ```
  <eventNumber> <nParticlesInEvent>
  <px> <py> <pz> <bacterialID>       (one line per bacterium, repeated
  ...                                 nParticlesInEvent times)
  ```
- Some events contain 0 bacteria (e.g. a header line like `77 0`). These
  are excluded from the event count used to calculate averages, since
  nothing was observed in them.

### Bacterial ID → strain mapping

| ID    | Strain                                        |
| ----- | --------------------------------------------- |
| 211   | *E. coli* WT (wild type)                      |
| -211  | *E. coli* mutant                              |
| 321   | *Bacillus subtilis* WT                        |
| -321  | *Bacillus subtilis* mutant                    |
| 2212  | *Pseudomonas aeruginosa* WT                   |
| -2212 | *Pseudomonas aeruginosa* antibiotic-resistant |
| 3122  | *Streptococcus pneumoniae*                    |
| -3122 | Capsule-deficient *S. pneumoniae*             |
| 3312  | *Mycobacterium tuberculosis*                  |
| -3312 | Drug-resistant *M. tuberculosis*              |
| 3334  | *Salmonella enterica*                         |
| -3334 | *Salmonella* mutant                           |

## 3. Method

**Step 1 — per-file analysis** (`Week3deliverable.R` for Set 1,
`output2.R`...`output10.R` for Sets 2–10): each script reads one
~800MB file in memory-safe chunks (never loading the whole file at once),
counts how many times each of the 12 strain IDs appears, and computes for
that single file:
- `average per event = total count / number of valid events`
- `uncertainty = sqrt(total count) / number of valid events` (Poisson
  statistics — see week 3)

**Step 2 — combining the 10 sub-samples** (`Week4deliverable.R`): the 10
per-file results are combined into one overall measurement per strain,
using the **sub-sampling method**:
- **Central value** — inverse-variance weighted average across the 10
  sub-samples: <br>
  `⟨x⟩ = Σ wᵢ⟨xᵢ⟩ / Σ wᵢ`, where `wᵢ = 1 / σᵢ²`
- **Uncertainty** — the standard deviation of the 10 sub-sample averages
  themselves (the spread of independent measurements of the same
  quantity), rather than each sub-sample's own Poisson uncertainty.

The 10 files disagree with each other more than the simple per-file
Poisson formula predicts (see section 6 for the actual numbers). In other
words, bacteria within an event aren't showing up purely independently of
each other — something is making the counts swing more from file to file
than random chance alone would cause. That's why the uncertainty here
comes from directly measuring how much the 10 files disagree, rather than
trusting the Poisson formula's guess.

## 4. Main result

Combined across all 10 sub-samples (**4,617,993 valid events**, out of
5,000,000 total — the remainder had 0 bacteria recorded):

| ID    | Strain                                        | Average per event | Uncertainty |
| ----- | --------------------------------------------- | ------------------| ------------|
| 211   | *E. coli* WT                                  |          19.94948 |   ± 0.03274 |
| -211  | *E. coli* mutant                              |          19.91716 |   ± 0.03187 |
| 321   | *Bacillus subtilis* WT                        |           2.50914 |   ± 0.00477 |
| -321  | *Bacillus subtilis* mutant                    |           2.50345 |   ± 0.00550 |
| 2212  | *Pseudomonas aeruginosa* WT                   |           1.20803 |   ± 0.00190 |
| -2212 | *Pseudomonas aeruginosa* antibiotic-resistant |           1.18416 |   ± 0.00241 |
| 3122  | *Streptococcus pneumoniae*                    |           0.27659 |   ± 0.00107 |
| -3122 | Capsule-deficient *S. pneumoniae*             |           0.27169 |   ± 0.00098 |
| 3312  | *Mycobacterium tuberculosis*                  |           0.03944 |   ± 0.00028 |
| -3312 | Drug-resistant *M. tuberculosis*              |           0.03900 |   ± 0.00040 |
| 3334  | *Salmonella enterica*                         |           0.00119 |   ± 0.00004 |
| -3334 | *Salmonella* mutant                           |           0.00115 |   ± 0.00005 |

## 5. Is there an asymmetry between WT and mutant/resistant strains?

For each pair, the difference is computed **within each of the 10 files
first** (`ΔNᵢ = averageᵢ(X) − averageᵢ(−X)` for file i), rather than
subtracting the two already-combined final averages. The 10 resulting
differences are then treated the same way as any other sub-sampling
result: the mean of the 10 differences is the reported difference, and
the standard deviation of the 10 differences is the uncertainty on it
(same sub-sampling convention as section 4 — no further division).

This matters because WT and mutant/resistant counts are not independent:
a file with more bacteria overall tends to have more of *both* the normal
and mutant form together. Subtracting within each file first cancels out
that shared file-to-file noise before combining, which is more accurate
than combining the two strains' already-finished uncertainties as if they
were unrelated (`sqrt(σ_WT² + σ_mutant²)` — a formula that assumes no
correlation between them, which does not hold here).

The significance is the difference divided by its uncertainty, expressed
as a number of standard deviations ("sigma") — the same scale used in the
lecture (slide 11: 3σ ≈ 1-in-740 chance of being a statistical fluctuation,
5σ ≈ discovery-level significance).

| Strain pair                                         | Difference | Uncertainty | Significance             |
| --------------------------------------------------- | ---------: | ----------: | ------------------------ |
| *E. coli* (WT vs mutant)                            |    0.03230 |     0.00452 | **7.15σ — significant**  |
| *Bacillus subtilis* (WT vs mutant)                  |    0.00569 |     0.00329 | 1.73σ — not significant  |
| *Pseudomonas aeruginosa* (WT vs resistant)          |    0.02387 |     0.00237 | **10.07σ — significant** |
| *Streptococcus pneumoniae* (normal vs capsule-def.) |    0.00490 |     0.00058 | **8.46σ — significant**  |
| *Mycobacterium tuberculosis* (WT vs drug-resistant) |    0.00044 |     0.00049 | 0.90σ — not significant  |
| *Salmonella* (WT vs mutant)                         |    0.00004 |     0.00007 | 0.53σ — not significant  |

**Conclusion:** three of the six strain pairs show a statistically
significant difference between the normal and mutant/resistant form:
*E. coli* (7.15σ), *Pseudomonas aeruginosa* (10.07σ), and *Streptococcus
pneumoniae* (8.46σ), in each case the WT/normal form is more abundant
than the mutant/resistant form. The other three pairs (*Bacillus
subtilis*, *Mycobacterium tuberculosis*, *Salmonella*) show no significant
asymmetry; their averages agree within uncertainty.

## 6. Consistency check: do the 10 sub-samples agree with their own precision?

Each individual file reports its own Poisson uncertainty (`sqrt(count) /
nEvents`), a prediction of how much that file's answer should wobble due
to pure random counting noise alone. If that prediction is accurate, the
spread across the 10 independent files (the number used as the final
uncertainty in section 4) should come out close to that per-file Poisson
uncertainty. The ratio of the two tells us whether the naive Poisson
formula correctly describes the real spread, or underestimates it.

| Strain                                        | Spread across 10 files | Individual file's Poisson unc. |     Ratio |
| --------------------------------------------- | ---------------------: | -----------------------------: | --------: |
| *E. coli* WT                                  |                0.03274 |                        0.00657 | **4.98x** |
| *E. coli* mutant                              |                0.03187 |                        0.00657 | **4.85x** |
| *Bacillus subtilis* WT                        |                0.00477 |                        0.00233 |     2.05x |
| *Bacillus subtilis* mutant                    |                0.00550 |                        0.00233 |     2.36x |
| *Pseudomonas aeruginosa* WT                   |                0.00190 |                        0.00162 |     1.17x |
| *Pseudomonas aeruginosa* antibiotic-resistant |                0.00241 |                        0.00160 |     1.51x |
| *Streptococcus pneumoniae*                    |                0.00107 |                        0.00077 |     1.39x |
| Capsule-deficient *S. pneumoniae*             |                0.00098 |                        0.00077 |     1.28x |
| *Mycobacterium tuberculosis*                  |                0.00028 |                        0.00029 |     0.97x |
| Drug-resistant *M. tuberculosis*              |                0.00040 |                        0.00029 |     1.39x |
| *Salmonella enterica*                         |                0.00004 |                        0.00005 |     0.83x |
| *Salmonella* mutant                           |                0.00005 |                        0.00005 |     1.04x |

**Conclusion:** the ratio is close to 1× (consistent with pure Poisson
noise) for the rarest strains, but climbs as high as ~5× for the most
abundant strain (E. coli). This is a systematic, abundance-dependent
pattern, not random noise: for high-count strains, the individual files'
Poisson uncertainty clearly *understates* the real uncertainty, since it
assumes fully independent random counting and ignores real extra
correlation/clustering in how bacteria appear within an event. This is
exactly why the sub-sampling method (the empirical spread across 10
independent files) is used as the final reported uncertainty in section 4,
rather than trusting any single file's own Poisson estimate.

## 7. Asymmetry as a function of momentum

**Method** (`momentumAsymmetry.R`): every occurrence of one of the 12
target strains, across all 10 files (the full ~5M-event sample), has its
momentum computed (`p = sqrt(px² + py² + pz²)`) and is sorted into one of
5 momentum bins (boundaries computed from Set 1's quintiles, then reused
for every file so all 10 are binned consistently). Within each bin, the
**sub-sampling method** is used, matching sections 4-6 of this report
rather than a plain Poisson-counting formula: each of the 10 files gets
its own asymmetry for that bin (`A_file = (count_WT - count_mutant) /
(count_WT + count_mutant)`), and the reported asymmetry is the mean of
the 10 per-file values, with the standard deviation of those 10 values as
the uncertainty. (A file with zero particles of a strain in a bin is
excluded from that bin's mean/SD, since 0/0 is undefined - this did not
happen for any strain/bin here, all 30 rows used all 10 files.)

Significance = |asymmetry| / uncertainty, tested against zero (i.e. "is
there any asymmetry at all in this momentum range"), using the same 3σ
threshold as the rest of this report.

| Strain                     | Momentum range |   Count WT | Count mutant | Asymmetry          |    Sigma | Significant? |
| -------------------------- | -------------- | ---------: | -----------: | ------------------ | -------: | ------------ |
| E. coli                    | < 0.65         | 18,387,463 |   18,384,806 | 0.00029 ± 0.00082  |     0.35 |              |
| E. coli                    | 0.65 - 1.32    | 17,722,346 |   17,698,495 | 0.00067 ± 0.00052  |     1.29 |              |
| E. coli                    | 1.32 - 2.68    | 19,507,546 |   19,482,867 | 0.00066 ± 0.00041  |     1.61 |              |
| E. coli                    | 2.68 - 5.87    | 21,361,848 |   21,314,609 | 0.00109 ± 0.00032  | **3.41** | **yes**      |
| E. coli                    | > 5.87         | 15,147,485 |   15,096,765 | 0.00293 ± 0.00401  |     0.73 |              |
| Bacillus subtilis          | < 0.65         |  1,238,914 |    1,237,780 | 0.00093 ± 0.00202  |     0.46 |              |
| Bacillus subtilis          | 0.65 - 1.32    |  2,029,612 |    2,031,222 | -0.00022 ± 0.00188 |     0.12 |              |
| Bacillus subtilis          | 1.32 - 2.68    |  2,605,834 |    2,605,381 | 0.00005 ± 0.00131  |     0.04 |              |
| Bacillus subtilis          | 2.68 - 5.87    |  3,017,075 |    3,006,537 | 0.00164 ± 0.00086  |     1.91 |              |
| Bacillus subtilis          | > 5.87         |  2,695,792 |    2,680,026 | 0.00151 ± 0.00483  |     0.31 |              |
| Pseudomonas aeruginosa     | < 0.65         |    337,024 |      331,332 | 0.00801 ± 0.00416  |     1.93 |              |
| Pseudomonas aeruginosa     | 0.65 - 1.32    |    778,588 |      766,503 | 0.00811 ± 0.00333  |     2.44 |              |
| Pseudomonas aeruginosa     | 1.32 - 2.68    |  1,249,066 |    1,229,314 | 0.00803 ± 0.00231  | **3.48** | **yes**      |
| Pseudomonas aeruginosa     | 2.68 - 5.87    |  1,580,395 |    1,549,824 | 0.00982 ± 0.00295  | **3.33** | **yes**      |
| Pseudomonas aeruginosa     | > 5.87         |  1,633,620 |    1,591,474 | 0.01674 ± 0.01181  |     1.42 |              |
| Streptococcus pneumoniae   | < 0.65         |     63,000 |       61,838 | 0.00419 ± 0.01863  |     0.22 |              |
| Streptococcus pneumoniae   | 0.65 - 1.32    |    155,835 |      153,540 | 0.00748 ± 0.00614  |     1.22 |              |
| Streptococcus pneumoniae   | 1.32 - 2.68    |    278,197 |      274,302 | 0.00705 ± 0.00337  |     2.09 |              |
| Streptococcus pneumoniae   | 2.68 - 5.87    |    372,039 |      365,814 | 0.00816 ± 0.00454  |     1.80 |              |
| Streptococcus pneumoniae   | > 5.87         |    408,259 |      399,196 | 0.00979 ± 0.00503  |     1.95 |              |
| Mycobacterium tuberculosis | < 0.65         |      6,932 |        6,891 | 0.00261 ± 0.03556  |     0.07 |              |
| Mycobacterium tuberculosis | 0.65 - 1.32    |     18,756 |       18,621 | 0.00405 ± 0.01339  |     0.30 |              |
| Mycobacterium tuberculosis | 1.32 - 2.68    |     37,778 |       37,410 | 0.00466 ± 0.01005  |     0.46 |              |
| Mycobacterium tuberculosis | 2.68 - 5.87    |     54,399 |       53,968 | 0.00445 ± 0.01029  |     0.43 |              |
| Mycobacterium tuberculosis | > 5.87         |     64,274 |       63,214 | 0.01199 ± 0.01549  |     0.77 |              |
| Salmonella                 | < 0.65         |        164 |          146 | 0.05829 ± 0.18672  |     0.31 |              |
| Salmonella                 | 0.65 - 1.32    |        466 |          428 | 0.03333 ± 0.09047  |     0.37 |              |
| Salmonella                 | 1.32 - 2.68    |      1,037 |          979 | 0.02705 ± 0.07961  |     0.34 |              |
| Salmonella                 | 2.68 - 5.87    |      1,661 |        1,628 | 0.01485 ± 0.05539  |     0.27 |              |
| Salmonella                 | > 5.87         |      2,154 |        2,137 | 0.02252 ± 0.07785  |     0.29 |              |

**Does the asymmetry itself change with momentum?** The table above tests
each bin against zero. A separate question is whether asymmetry *changes*
across momentum - tested here the same way WT vs. mutant was compared in
section 5: highest momentum bin minus lowest momentum bin, using each
strain's own combined uncertainty.

| Strain                     | Highest bin − lowest bin | Significance            |
| -------------------------- | -----------------------: | ----------------------- |
| E. coli                    |          0.0026 ± 0.0041 | 0.65σ — not significant |
| Bacillus subtilis          |          0.0006 ± 0.0052 | 0.11σ — not significant |
| Pseudomonas aeruginosa     |          0.0087 ± 0.0125 | 0.70σ — not significant |
| Streptococcus pneumoniae   |          0.0056 ± 0.0193 | 0.29σ — not significant |
| Mycobacterium tuberculosis |          0.0094 ± 0.0388 | 0.24σ — not significant |
| Salmonella                 |         -0.0358 ± 0.2023 | 0.18σ — not significant |

**Conclusion:** out of 30 strain/momentum-bin combinations, only 3 show a
statistically significant asymmetry *against zero*: *E. coli* in its
highest-momentum bin (3.41σ), and *Pseudomonas aeruginosa* in its two
highest-momentum bins (3.48σ and 3.33σ). But testing whether asymmetry
*changes* across momentum (highest bin vs. lowest bin, table above) shows
no significant trend for any strain - every strain's momentum dependence
is consistent with being flat. These two findings don't contradict each
other: asymmetry is real in a few specific momentum ranges, but it isn't
detectably *growing or shrinking* with momentum across the full range.

**Why the highest-momentum bin is so much noisier for the common strains**
(E. coli, Bacillus subtilis, Pseudomonas): one file is a major outlier
there. In E. coli's highest bin, Sets 1-9 each have roughly 1.68 million
WT particles, but **Set 10 has only 8,032** - over 200x fewer. This isn't
a data quality problem (Set 10's event and file structure are otherwise
identical to every other file - see section 2); Set 10's particles simply
have systematically lower momentum on average across the board (consistent
with its data lines being ~5% shorter than Set 1's, i.e. smaller numbers
throughout), so very few of them clear the high-momentum threshold that
was calibrated from Set 1. That one outlier file dominates the spread
across the 10 sub-samples for this specific bin, inflating its
uncertainty well beyond the other bins'.


## 8. How to reproduce

1. Download the 10 data files from surfdrive and place them in this
   folder (not committed to the repo — too large).
2. Run `analyzeAllSubsamples.R` once — analyzes all 10 files (one shared
   function, looped over each filename) and writes `subsample_results.csv`.
3. Run `Week4deliverable.R` — reads `subsample_results.csv`, combines the
   10 sub-samples into the final result above, and writes
   `final_results.csv`. Runs instantly, since it doesn't re-read the raw
   data files.
4. Run `asymmetryTest.R` — reads `subsample_results.csv`, runs the
   paired-difference asymmetry test for the 6 WT/mutant pairs, and writes
   `asymmetry_results.csv`. Also runs instantly.
5. Run `momentumAsymmetry.R` — re-reads all 10 raw data files (needs
   per-particle momentum, which none of the counting-only scripts above
   track), uses the sub-sampling method per momentum bin, and writes
   `momentum_asymmetry_results.csv` plus a cached per-file breakdown
   (`momentum_perfile_results.csv`). Takes about 20-25 minutes.

## 9. Files in this repository

| File                             | Purpose                                                                                               |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| `analyzeAllSubsamples.R`         | Analyzes all 10 data files, writes `subsample_results.csv`                                            |
| `Week4deliverable.R`             | Combines the 10 sub-samples into the final result                                                     |
| `asymmetryTest.R`                | Runs the paired-difference asymmetry test for the 6 pairs                                             |
| `momentumAsymmetry.R`            | Runs the momentum-binned asymmetry test (re-reads raw data)                                           |
| `subsample_results.csv`          | Per-file results (10 files x 12 strains, 120 rows)                                                    |
| `final_results.csv`              | Final combined results (12 strains)                                                                   |
| `asymmetry_results.csv`          | Asymmetry test results (6 pairs)                                                                      |
| `momentum_asymmetry_results.csv` | Momentum-binned asymmetry results (6 pairs x 5 bins)                                                  |
| `momentum_perfile_results.csv`   | Cached per-file momentum bin counts (avoids re-reading raw data if the momentum method changes again) |
| `README_JW.md`                   | This file                                                                                             |
