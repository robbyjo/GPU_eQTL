# Synthetic BCF dictionary fixtures

`sparse-dictionary.vcf` contains three invented variants and four invented samples.
`sparse-dictionary.bcf.base64` is the BGZF BCF 2.2 encoding produced by bcftools
1.23/htslib 1.23 with `bcftools view --no-version -Ob`. Its shared tag indices
include 0, 8, 27, 31 and 53, with AF declared before lower-index tags. LowQual is
both an INFO flag and a FILTER under the same dictionary index. The AF description
contains quoted `IDX=999` text that is not a header attribute.

`sparse-contigs.bcf.base64` additionally gives chromosomes 1 and 2 the sparse
contig indices 3 and 7. Genotypes, dosage, phasing, missing values, INFO and FILTER
values are otherwise identical. No cohort data or production identifiers are used.

Run `scripts/generate-sparse-bcf-fixtures.sh NEW_DIRECTORY` to regenerate the
native fixtures and a version/checksum receipt. Compression bytes can vary with
the native library version; decoded values must match the supplied VCF oracle.
Tests mutate only headers of these tiny synthetic records to exercise invalid
indices and verify that gaps are never treated as declared fields.

The mapping follows HTSlib 1.23's `bcf_hdr_set_idx` and
`bcf_hdr_register_hrec` in [vcf.c](https://github.com/samtools/htslib/blob/1.23/vcf.c):
first occurrences register immediately, absent IDX uses the current dictionary
length (largest preceding index plus one), contigs have a separate dictionary,
and INFO/FILTER/FORMAT share IDs. The Java adapter deliberately rejects
inconsistent repeated-ID indices instead of silently keeping the first mapping.
The reader uses a checked bridge to the pinned HTSJDK 5.0.0 lookup fields because
its dictionary initializers are private and its lookup methods are final. These
tests must continue to pass when upgrading that dependency.
