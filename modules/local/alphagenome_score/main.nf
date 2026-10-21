process ALPHAGENOME_SCORE {

    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b4aad86558f56aa4eda7b2cc8d89c79ee9459c4760230025e38647e701eefb4b/data'
        : 'community.wave.seqera.io/library/alphagenome-score:1e780ad6ad243bba'}"

    input:
    tuple val(meta), path(vcf)
    path(api_key_file)
    val organism
    val sequence_length
    val max_variants
    val max_scorers
    val scorers

    output:
    tuple val(meta), path("*.alphagenome.vcf"), emit: vcf
    tuple val(meta), path("*.alphagenome_summary_scores.tsv"), optional: true, emit: summary_scores
    tuple val(meta),path("*.alphagenome_tidy_scores.tsv"), optional: true, emit: tidy_scores
    path "versions.yml", emit: versions


script:
"""
python3 ${projectDir}/bin/run_alphagenome.py \\
    --vcf ${vcf} \\
    --api-key-file ${api_key_file} \\
    --organism ${organism} \\
    --sequence-length ${sequence_length} \\
    --output ${meta.id}.alphagenome.vcf \\
    --summary-scores ${meta.id}.alphagenome_summary_scores.tsv \\
    --tidy-scores ${meta.id}.alphagenome_tidy_scores.tsv \\
    --max-variants ${max_variants} \\
    --max-scorers ${max_scorers} \\
    --scorers "${scorers}"
    

echo "Files after AlphaGenome script:"
ls -lah

test -s ${meta.id}.alphagenome.vcf || {
    echo "ERROR: ${meta.id}.alphagenome.vcf was not created or is empty" >&2
    exit 1
}

test -s ${meta.id}.alphagenome_summary_scores.tsv || {
    echo "ERROR: ${meta.id}.alphagenome_summary_scores.tsv was not created or is empty" >&2
    exit 1
}

test -f ${meta.id}.alphagenome_tidy_scores.tsv || {
    echo "ERROR: ${meta.id}.alphagenome_tidy_scores.tsv was not created" >&2
    exit 1
}

cat > versions.yml <<-END_VERSIONS
"${task.process}":
    python: "\$(python3 -c 'import platform; print(platform.python_version())')"
    alphagenome: "\$(python3 -c 'import importlib.metadata as metadata; print(metadata.version("alphagenome"))')"
    cyvcf2: "\$(python3 -c 'import importlib.metadata as metadata; print(metadata.version("cyvcf2"))')"
    pandas: "\$(python3 -c 'import importlib.metadata as metadata; print(metadata.version("pandas"))')"
END_VERSIONS
"""

}