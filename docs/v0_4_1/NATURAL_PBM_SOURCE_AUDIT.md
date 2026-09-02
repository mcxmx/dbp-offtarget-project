# v0.4.1 Natural PBM Source Audit

Audit date: 2026-09-02

## Candidate Sources

| database | assay_type | number_of_proteins | protein_sequence_availability | kmer_length | score_type | replicate_availability | raw_probe_availability | species | family_diversity | downloadability | license | recommended_use |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| UniPROBE | universal PBM | publication-specific; v0.4.1 downloaded 112 profile files | factor names in profile paths; construct sequences are not generally embedded in contig8mer files | 8 | contiguous 8-mer E-score plus median/z-score where available | replicate IDs inferable for some entries | deBruijn probe files available but license-protected | mixed by publication | broad across UniPROBE publications | direct publication zip URLs work | academic research use license for files containing PBM probe sequences; E-score profiles are publicly downloadable from UniPROBE | PRIMARY SOURCE for v0.4.1 natural PBM pilot/training benchmark |
| CIS-BP 3.10 | PBM/PWM compiled specificity catalog | large database | TF information and protein metadata available in database downloads | mostly 8-mer E-scores/PWMs | E-score, Z-score, intensity, motif/PWM | varies | available in bulk archives | broad | broad | full E-score archive is ~1.19 GB; species POST archive was tested but large download was interrupted in this environment | CIS-BP database terms/citation required | SECONDARY SOURCE for future sequence metadata and larger natural PBM benchmark |
| DREAM/PBM challenge-style datasets | PBM benchmark/challenge | challenge-specific | varies | probe or k-mer level depending release | challenge-specific | varies | varies | mainly natural TFs | moderate | not selected in this pass | dataset-specific | Future validation source after harmonization audit |

## v0.4.1 Selection

PRIMARY SOURCE: UniPROBE publication-level contiguous 8-mer E-score files.

Selected UniPROBE dataset codes: NAR10, NBT06, DEV12, PNAS12, Path10, MBE14, PO10, NAR11, PNAS08, GB11, LIU18A

Rationale: the files are directly downloadable, have explicit contiguous 8-mer E-score tables, and each profile contains reverse-complement-paired 8-mer rows suitable for per-protein ranking. This round does not claim construct-level protein sequence recovery from UniPROBE contig8mer files.

## Important Limitation

UniPROBE contig8mer profile files do not reliably provide experimental construct sequences. v0.4.1 therefore treats protein sequence recovery as a separate provenance task. Reference full-length sequences, when later added from UniProt, must be flagged as `sequence_match_to_assay=false` unless construct sequence evidence is available.
