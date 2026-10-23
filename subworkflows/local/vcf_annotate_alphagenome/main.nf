//
// Run AlphaGenome to annotate VCF files
//

include { ALPHAGENOME_SCORE } from '../../../modules/local/alphagenome_score'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix'

workflow VCF_ANNOTATE_ALPHAGENOME {
    take:
    ch_vcf          // channel: [ val(meta), path(vcf) ]
    api_key_file       // channel: [ path(api_key_file) ]
    val_organism  // string:  human, mouse
    val_sequence_length // string:  sequence length (e.g., '1Mb')
    val_max_variants // string: maximum number of variants to consider (default: 3)')
    val_max_scorers // string: maximum number of scorers to consider (default: 1)')
    val_scorers // string: comma-separated list of scorer groups to use (e.g., 'atac,dnase,rna_seq'). If not set, val_max_scorers is used.

    main:
    ch_versions = Channel.empty()

// Run AlphaGenome scoring

    ALPHAGENOME_SCORE(ch_vcf, api_key_file, val_organism, val_sequence_length, val_max_variants, val_max_scorers, val_scorers)

// Compress and index the output VCF file using bgzip and tabix
    TABIX_BGZIPTABIX(ALPHAGENOME_SCORE.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ALPHAGENOME_SCORE.out.versions)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    vcf_tbi   = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
    summary_scores   = ALPHAGENOME_SCORE.out.summary_scores // channel: [ val(meta), path(summary_scores.tsv) ]
    tidy_scores   = ALPHAGENOME_SCORE.out.tidy_scores // channel: [ val(meta), path(tidy_scores.tsv) ]
    versions  = ch_versions // channel: [ path(versions.yml) ]
}
