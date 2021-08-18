import pandas as pd

SAMPLES = glob.wildcards("DATA/RAW/{sample}_1.fq.gz")
print(SAMPLES)

# Pull reference sample name
with open('refseq_accessions.txt','r') as ref:
    tmp_name = pd.read_csv(ref, sep = '\t', header = None)
    name = str(tmp_name[tmp_name[1] == 'reference'][0].tolist()[0])

rule all:
    input:
        expand('DATA/CLEAN/{sample}_qual_R1.fq.gz', sample=SAMPLES), expand('DATA/CLEAN/{sample}_qual_R2.fq.gz', sample=SAMPLES),
        'DATA/REF/contigs.1.bt2',
        expand('DATA/SAM/{sample}_syn.sam', sample=SAMPLES),
        expand('DATA/MAP_BAM/{sample}_mapped.bam', sample=SAMPLES),
        expand('DATA/SYN_READS/{sample}_onlySyn1.fq',sample=SAMPLES),
        expand('DATA/SYN_READS_NORM/{sample}_norm_onlySyn1.fq',sample=SAMPLES), expand('DATA/SYN_READS_NORM/{sample}_norm_onlySyn2.fq',sample=SAMPLES),
        expand('DATA/SPADES_SYN/{sample}_spades_contigs.fa', sample=SAMPLES),
        'DATA/SYN_PARSNP/parsnp.tree',
        expand('DATA/UNMAP_BAM/{sample}_unmapped.bam', sample=SAMPLES),
        expand('DATA/NOTSYN_READS/{sample}_notSyn1.fq', sample=SAMPLES), expand('DATA/NOTSYN_READS/{sample}_notSyn2.fq', sample=SAMPLES),
        expand('DATA/METAPHLAN/{sample}_metaphlan_profile.txt', sample=SAMPLES), expand('DATA/METAPHLAN/{sample}_metaphlan.bowtie2.bz2', sample=SAMPLES),
        'DATA/METAPHLAN/merged_abundance_table_species.txt', 'DATA/METAPHLAN/merged_abundance_summary.csv',
        'DATA/METAPHLAN/metaphlan_full_taxonomy.csv',
        'DATA/METAPHLAN/abundance_heatmap_species.png',
        expand('DATA/HUMANN/{sample}_pathcoverage.tsv', sample=SAMPLES), expand('DATA/HUMANN/{sample}_pathabundance.tsv', sample=SAMPLES), expand('DATA/HUMANN/{sample}_genefamilies.tsv', sample=SAMPLES),
        'DATA/HUMANN/merged_genefamilies.tsv',
        'DATA/METAPHLAN/metaphlan_merged_abundance_table.txt', 'DATA/METAPHLAN/metaphlan_full_taxonomy.csv',
        'DATA/HUMANN/merged_pathabundance.tsv', 'DATA/HUMANN/merged_genefamilies_cpm.tsv', 'DATA/HUMANN/merged_pathcoverage.tsv'

rule clean:
    input:
        R1='DATA/RAW/{sample}_1.fq.gz',
        R2='DATA/RAW/{sample}_2.fq.gz'
    conda:'env1.yml'
    log: "LOGS/{sample}_clean.log"
    output:
        R1='DATA/CLEAN/{sample}_qual_R1.fq.gz',
        R2='DATA/CLEAN/{sample}_qual_R2.fq.gz'
    shell:
        'mkdir -p DATA/CLEAN && '
        'bbduk.sh overwrite=true in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq=25 2> {log}'

rule remove_syn_ind:
    input:'DATA/REF/'+name+'_ref_genome.fa'
    conda:'env1.yml'
    params:'DATA/REF/'+name+'_index'
    log: "LOGS/remove_syn_ind.log"
    output:'DATA/REF/contigs.1.bt2'
    shell:
        'REF_STRAIN=$(grep "reference" refseq_accessions.txt | cut -f1) && '
        'mkdir -p {params} && bowtie2-build {input} {params}/contigs 2>{log} && '
         'cp {params}/contigs.1.bt2 {output}'

rule syn_map:
    input:
        R1='DATA/CLEAN/{sample}_qual_R1.fq.gz',
        R2='DATA/CLEAN/{sample}_qual_R2.fq.gz',
        ind='DATA/REF/contigs.1.bt2'
    params:'REF/'+name+'_index'
    log: "LOGS/{sample}_syn_map.log"
    conda:'env1.yml'
    output:'DATA/SAM/{sample}_syn.sam'
    shell:
        'touch {input.ind} && mkdir -p DATA/SAM && bowtie2 --threads 16 -x {params}/contigs -1 {input.R1} -2 {input.R2} -S {output} 2> {log}'

# UPDATE: go straight to bam files at some point -> bowtie2 (...option) | samtools view -o out.bam

rule only_syn_mapped:
    input:'DATA/SAM/{sample}_syn.sam'
    conda: 'env1.yml'
    log: "LOGS/{sample}_only_syn_mapped.log"
    output:'DATA/MAP_BAM/{sample}_mapped.bam'
    shell:
        'mkdir -p DATA/MAP_BAM && (samtools view -F 4 -b {input} | samtools sort -n - > {output}) 2> {log}'

rule only_syn_subset:
    input:'DATA/MAP_BAM/{sample}_mapped.bam'
    conda:'env1.yml'
    log: "LOGS/{sample}_only_syn_subset.log"
    output:
        R1='DATA/SYN_READS/{sample}_onlySyn1.fq',
        R2='DATA/SYN_READS/{sample}_onlySyn2.fq'
    shell:
        'mkdir -p DATA/SYN_READS && bedtools bamtofastq -i {input} -fq {output.R1} -fq2 {output.R2} 2> {log}'

rule only_syn_norm:
    input:
        R1='DATA/SYN_READS/{sample}_onlySyn1.fq',
        R2='DATA/SYN_READS/{sample}_onlySyn2.fq'
    conda:'env1.yml'
    log: "LOGS/{sample}_only_syn_norm.log"
    output:
        R1='DATA/SYN_READS_NORM/{sample}_norm_onlySyn1.fq',
        R2='DATA/SYN_READS_NORM/{sample}_norm_onlySyn2.fq'
    shell:
        'mkdir -p DATA/SYN_READS_NORM && '
        'bbnorm.sh overwrite=true bits=16 prefilter=t threads=16 -Xmx99g \
        in={input.R1} in2={input.R2} \
        out={output.R1} out2={output.R2} target=60 min=5 2> {log}'

rule only_syn_assembly:
    input:
        R1='DATA/SYN_READS_NORM/{sample}_norm_onlySyn1.fq',
        R2='DATA/SYN_READS_NORM/{sample}_norm_onlySyn2.fq'
    conda:'env1.yml'
    log: "LOGS/{sample}_only_syn_assembly.log"
    output:'DATA/SPADES_SYN/{sample}_spades_contigs.fa'
    shell:
        'mkdir -p DATA/SPADES_SYN && spades.py --phred-offset 33 -t 16 -1 {input.R1} -2 {input.R2} -o DATA/SPADES_SYN/{wildcards.sample} 2> {log} && '
        'cp DATA/SPADES_SYN/{wildcards.sample}/contigs.fasta {output}'

rule parsnp:
    input:
        assemblies=expand('DATA/SPADES_SYN/{sample}_spades_contigs.fa',sample=SAMPLES),
        ref_LA31='DATA/REF/'+name+'_ref_genome.fa'
    conda:'env2.yml'
    log: "LOGS/Parsnp.log"
    output:'DATA/SYN_PARSNP/parsnp.tree'
    shell:
        'mkdir -p DATA/SYN_PARSNP && cp {input.assemblies} DATA/SYN_PARSNP/ && parsnp -r {input.ref_LA31} -d DATA/SYN_PARSNP -p 16 -o tmp_parsnp_out 2> {log} &&'
        'cp tmp_parsnp_out/*  DATA/SYN_PARSNP/ && '
        'rm -r tmp_parsnp_out && '
        'touch {output}'

rule remove_syn_bam:
    input:'DATA/SAM/{sample}_syn.sam'
    conda:'env1.yml'
    log: "LOGS/{sample}_remove_syn_bam.log"
    output:'DATA/UNMAP_BAM/{sample}_unmapped.bam'
    shell:
        'mkdir -p DATA/UNMAP_BAM && '
        '(samtools view -f 4 -b {input} | samtools sort -n -@ 16 -m 3G - > {output}) 2> {log}'

rule remove_syn_subset:
    input:
        bam='DATA/UNMAP_BAM/{sample}_unmapped.bam'
    conda:'env1.yml'
    log: "LOGS/{sample}_remove_syn_subset.log"
    output:
        R1='DATA/NOTSYN_READS/{sample}_notSyn1.fq',
        R2='DATA/NOTSYN_READS/{sample}_notSyn2.fq'
    shell:
        'mkdir -p DATA/NOTSYN_READS && '
        'bedtools bamtofastq -i {input.bam} -fq {output.R1} -fq2 {output.R2} 2> {log}'

rule metaphlan_profile:
    input:
        R1='DATA/CLEAN/{sample}_qual_R1.fq.gz',
        R2='DATA/CLEAN/{sample}_qual_R2.fq.gz'
    conda:'env4.yml'
    log: "LOGS/{sample}_metaphlan_profile.log"
    output:
        metaP='DATA/METAPHLAN/{sample}_metaphlan_profile.txt',
        bw2='DATA/METAPHLAN/{sample}_metaphlan.bowtie2.bz2'
    shell:
        'mkdir -p DATA/METAPHLAN && '
        'metaphlan {input.R1},{input.R2} --bowtie2db DBs/metaphlan_db/ --bowtie2out {output.bw2} --nproc 16 --input_type fastq --bt2_ps very-sensitive-local --add_viruses --unknown_estimation -o {output.metaP} 2> {log}'

# -x metaphlan_db/mpa_v30_CHOCOPhlAn_201901

rule metaphlan_merge:
    input:expand('DATA/METAPHLAN/{sample}_metaphlan_profile.txt',sample=SAMPLES)
    conda:'env4.yml'
    log: "LOGS/Metaphlan_merge.log"
    output:'DATA/METAPHLAN/metaphlan_merged_abundance_table.txt'
    shell:
        'merge_metaphlan_tables.py {input} > {output} 2> {log}'

rule metaphlan_convert_table:
    input:'DATA/METAPHLAN/metaphlan_merged_abundance_table.txt'
    conda:'env0.yml'
    log: "LOGS/Metaphlan_convert_table.log"
    output:
        spp_tab='DATA/METAPHLAN/merged_abundance_table_species.txt',
        sum_tab='DATA/METAPHLAN/merged_abundance_summary.csv'
    shell:
        "grep -E 's__|clade' {input} | sed 's/^.*s__//g' | cut -f1,3-14 | sed -e 's/clade_name/strain/g' | sed 's/_metaphlan_profile//g' > tmp.txt && "
        #"grep -E 'UNKNOWN' {input} | cut -f1,3-14 >> tmp.txt && "
        "mv tmp.txt {output.spp_tab} && "
        "Rscript SCRIPTS/quick_metaphlan_summary.r {output.spp_tab} 2> {log} && "
        "touch {output.sum_tab}"

rule metaphlan_full_tax:
    input:'DATA/METAPHLAN/metaphlan_merged_abundance_table.txt'
    output:'DATA/METAPHLAN/metaphlan_full_taxonomy.csv'
    log: "LOGS/Metaphlan_full_tax.log"
    shell:
        "grep 's__' {input} | cut -f1 | sed 's/|.__/,/g' | sed 's/k__//' > {output} 2> {log}"

rule metaphlan_hclust:
    input:'DATA/METAPHLAN/merged_abundance_table_species.txt'
    conda:'env4.yml'
    log: "LOGS/Metaphlan_hclust.log"
    output:'DATA/METAPHLAN/abundance_heatmap_species.png'
    shell:
        'hclust2.py -i {input} -o {output} --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300 2> {log}'

rule humann_profile:
    input:
        R1='DATA/CLEAN/{sample}_qual_R1.fq.gz',
        R2='DATA/CLEAN/{sample}_qual_R2.fq.gz'
    conda:'env4.yml'
    log: "LOGS/{sample}_humann_profile.log"
    output:
        pathcov='DATA/HUMANN/{sample}_pathcoverage.tsv',
        pathabund='DATA/HUMANN/{sample}_pathabundance.tsv',
        families='DATA/HUMANN/{sample}_genefamilies.tsv'
    shell:
        'mkdir -p DATA/HUMANN && cat {input.R1} {input.R2} > merge_{wildcards.sample}_notSyn.fq && '
        'humann --nucleotide-database DBs/humann_db --output-basename {wildcards.sample} --input merge_{wildcards.sample}_notSyn.fq --output DATA/HUMANN 2> {log} && '
        'rm merge_{wildcards.sample}_notSyn.fq && '
        'touch {output.pathcov} && touch {output.pathabund} && touch {output.families}'

rule humann_merge_cpm:
    input:expand('DATA/HUMANN/{sample}_genefamilies.tsv',sample=SAMPLES)
    conda:'env4.yml'
    log: "LOGS/Humann_merge_cpm.log"
    output:'DATA/HUMANN/merged_genefamilies.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input HUMANN --output {output} --file_name _genefamilies 2> {log}'

rule humann_norm_merged_table:
    input:'DATA/HUMANN/merged_genefamilies.tsv'
    conda:'env4.yml'
    log: "LOGS/Humann_norm_merged_table.log"
    output:'DATA/HUMANN/merged_genefamilies_cpm.tsv'
    shell:
        'humann_renorm_table --input {input} --output {output} 2> {log}'

rule humann_merge_coverage:
    input:expand('DATA/HUMANN/{sample}_pathcoverage.tsv',sample=SAMPLES)
    conda:'env4.yml'
    log: "LOGS/Humann_merge_coverage.log"
    output:'HUMANN/merged_pathcoverage.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input humann_out --output {output} --file_name pathcoverage 2> {log}'

rule humann_merge_abundance:
    input:expand('DATA/HUMANN/{sample}_pathabundance.tsv',sample=SAMPLES)
    conda:'env4.yml'
    log: "LOGS/Humann_merge_abundance.log"
    output:'DATA/HUMANN/merged_pathabundance.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input HUMANN --output {output} --file_name pathabundance 2> {log}'
