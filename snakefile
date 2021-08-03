SAMPLES=['LA31','LA3','LA21','LA20','LA27','LA29','LA101','LA103','LA117','LA126','LA127']
SCRIPTS='../Scripts'
CONDA=''

rule all:
    input:
        expand('UNMAP_BAM/{sample}_unmapped.bam',sample=SAMPLES),
        expand('MAP_BAM/{sample}_mapped.bam',sample=SAMPLES),
        expand('SPADES_SYN/{sample}_spades_contigs.fa',sample=SAMPLES),
        #expand('NOTSYN_NORM_READS/{sample}_norm_notSyn1.fq',sample=SAMPLES),
        #expand('NOTSYN_NORM_READS/{sample}_norm_notSyn2.fq',sample=SAMPLES),
        expand('SYN_READS_NORM/{sample}_norm_onlySyn1.fq',sample=SAMPLES),
        expand('SYN_READS_NORM/{sample}_norm_onlySyn2.fq',sample=SAMPLES),
        'SYN_PARSNP/parsnp.tree',
        #'NOTSYN_ASSEMBLIES/norm_noSyn_assembly_megahit.fa',
        'NOTSYN_ASSEMBLIES/noSyn_assembly_metaspades.fa',
        expand('BAM/{sample}_metaspades.bam',sample=SAMPLES),
        expand('METAPHLAN/metaphlan_merged_abundance_table.txt'),
        'METAPHLAN/metaphlan_merged_abundance_table.txt','METAPHLAN/metaphlan_full_taxonomy.csv',
        'HUMANN/merged_pathabundance.tsv','HUMANN/merged_genefamilies_cpm.tsv','HUMANN/merged_pathcoverage.tsv',
        'BINS/METABAT2/bin.1.fa','BINS/CONCOCT/0.fa'

rule clean:
    input:
        R1='RAW/{sample}_1.fq.gz',
        R2='RAW/{sample}_2.fq.gz'
    conda:'env1.yml'
    output:
        R1='CLEAN/{sample}_qual_R1.fq.gz',
        R2='CLEAN/{sample}_qual_R2.fq.gz'
    shell:
        'mkdir -p CLEAN &&'
        'bbduk.sh overwrite=true in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq=25'

rule remove_syn_ind:
    input:'REF/LA31_falcon_pacbio_simple.fa'
    conda:'env1.yml'
    params:'REF/LA31_index'
    output:'REF/contigs.1.bt2'
    shell:
        'mkdir -p {params} && bowtie2-build {input} {params}/contigs && cp {params}/contigs.1.bt2 {output}'

rule syn_map:
    input:
        R1='CLEAN/{sample}_qual_R1.fq.gz',
        R2='CLEAN/{sample}_qual_R2.fq.gz',
        ind='REF/contigs.1.bt2'
    params:'REF/LA31_index'
    conda:'env1.yml'
    output:'SAM/{sample}_syn.sam'
    shell:
        'touch {input.ind} && mkdir -p SAM && bowtie2 --threads 16 -x {params}/contigs -1 {input.R1} -2 {input.R2} -S {output}'

# UPDATE: go straight to bam files at some point -> bowtie2 (...option) | samtools view -o out.bam

rule only_syn_mapped:
    input:'SAM/{sample}_syn.sam'
    conda: 'env1.yml'
    output:'MAP_BAM/{sample}_mapped.bam'
    shell:
        'mkdir -p MAP_BAM && samtools view -F 4 -b {input} | samtools sort -n - > {output}'

rule only_syn_subset:
    input:'MAP_BAM/{sample}_mapped.bam'
    conda:'env1.yml'
    output:
        R1='SYN_READS/{sample}_onlySyn1.fq',
        R2='SYN_READS/{sample}_onlySyn2.fq'
    shell:
        'mkdir -p SYN_READS && bedtools bamtofastq -i {input} -fq {output.R1} -fq2 {output.R2}'

rule only_syn_norm:
    input:
        R1='SYN_READS/{sample}_onlySyn1.fq',
        R2='SYN_READS/{sample}_onlySyn2.fq'
    conda:'env1.yml'
    output:
        R1='SYN_READS_NORM/{sample}_norm_onlySyn1.fq',
        R2='SYN_READS_NORM/{sample}_norm_onlySyn2.fq'
    shell:
        'mkdir -p SYN_READS_NORM &&'
        'bbnorm.sh overwrite=true bits=16 prefilter=t threads=16 -Xmx99g \
        in={input.R1} in2={input.R2} \
        out={output.R1} out2={output.R2} target=60 min=5'

rule only_syn_assembly:
    input:
        R1='SYN_READS_NORM/{sample}_norm_onlySyn1.fq',
        R2='SYN_READS_NORM/{sample}_norm_onlySyn2.fq'
    conda:'env1.yml'
    params:
        dir='SPADES_SYN'
    output:'SPADES_SYN/{sample}_spades_contigs.fa'
    shell:
        'mkdir -p {params.dir} && spades.py --phred-offset 33 -t 16 -1 {input.R1} -2 {input.R2} -o {params.dir}/{wildcards.sample} && cp {params.dir}/{wildcards.sample}/contigs.fasta {output}'

rule parsnp:
    input:
        assemblies=expand('SPADES_SYN/{sample}_spades_contigs.fa',sample=SAMPLES),
        ref_LA31='REF/LA31_falcon_pacbio_simple.fa'
    conda:'env2.yml'
    params:
        dir='SYN_PARSNP'
    output:'SYN_PARSNP/parsnp.tree'
    shell:
        'mkdir -p {params.dir} && cp {input.assemblies} ./{params.dir}/ && parsnp -r {input.ref_LA31} -d {params.dir} -p 16 -o tmp_parsnp_out &&'
        'cp tmp_parsnp_out/* {params.dir}/ &&'
        'rm -r tmp_parsnp_out &&'
        'touch {output}'

rule remove_syn_bam:
    input:'SAM/{sample}_syn.sam'
    conda:'env1.yml'
    output:'UNMAP_BAM/{sample}_unmapped.bam'
    shell:
        'mkdir -p UNMAP_BAM &&'
        'samtools view -f 4 -b {input} | samtools sort -n -@ 16 -m 3G - > {output}'

rule remove_syn_subset:
    input:
        bam='UNMAP_BAM{sample}_unmapped.bam'
    conda:'env1.yml'
    output:
        R1='NOTSYN_READS/{sample}_notSyn1.fq',
        R2='NOTSYN_READS/{sample}_notSyn2.fq'
    shell:
        'mkdir -p NOTSYN_READS &&'
        'bedtools bamtofastq -i {input.bam} -fq {output.R1} -fq2 {output.R2}'

# need to install metaphlan dbs outside of SNAKEMAKE
## /scratch/joshuakl/syn_final/.snakemake/conda/0eea8dbb/bin/metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_db

rule metaphlan_profile:
    input:
        R1='CLEAN/{sample}_qual_R1.fq.gz',
        R2='CLEAN/{sample}_qual_R2.fq.gz'
    conda:'env6.yml'
    output:
        metaP='METAPHLAN/{sample}_metaphlan_profile.txt',
        bw2='METAPHLAN/{sample}_metaphlan.bowtie2.bz2'
    shell:
        'mkdir -p METAPHLAN && '
        'metaphlan {input.R1},{input.R2} --bowtie2db metaphlan_db/ --bowtie2out {output.bw2} --nproc 16 --input_type fastq --bt2_ps very-sensitive-local --add_viruses -o {output.metaP}'

# -x metaphlan_db/mpa_v30_CHOCOPhlAn_201901

rule metaphlan_merge:
    input:expand('METAPHLAN/{sample}_metaphlan_profile.txt',sample=SAMPLES)
    conda:'env6.yml'
    output:'METAPHLAN/metaphlan_merged_abundance_table.txt'
    shell:
        'merge_metaphlan_tables.py {input} > {output}'

rule metaphlan_convert_table:
    input:'METAPHLAN/metaphlan_merged_abundance_table.txt'
    params:
        script=SCRIPTS
    conda:'env0.yml'
    output:
        spp_tab='METAPHLAN/merged_abundance_table_species.txt',
        sum_tab='METAPHLAN/merged_abundance_summary.csv'
    shell:
        "grep -E 's__|clade' {input} | sed 's/^.*s__//g' | cut -f1,3-14 | sed -e 's/clade_name/strain/g' | sed 's/_metaphlan_profile//g' > tmp.txt && "
        #"grep -E 'UNKNOWN' {input} | cut -f1,3-14 >> tmp.txt && "
        "mv tmp.txt {output.spp_tab} && "
        "Rscript {params}/quick_metaphlan_summary.r && "
        "touch {output.sum_tab}"

rule metaphlan_full_tax:
    input:'METAPHLAN/metaphlan_merged_abundance_table.txt'
    output:'METAPHLAN/metaphlan_full_taxonomy.csv'
    shell:
        "grep 's__' {input} | cut -f1 | sed 's/|.__/,/g' | sed 's/k__//' > {output}"

rule metaphlan_hclust:
    input:'METAPHLAN/merged_abundance_table_species.txt'
    conda:'env6.yml'
    output:'METAPHLAN/abundance_heatmap_species.png'
    shell:
        'hclust2.py -i {input} -o {output} --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300'

rule humann_profile:
    input:
        R1='CLEAN/{sample}_qual_R1.fq.gz',
        R2='CLEAN/{sample}_qual_R2.fq.gz'
    conda:'env6.yml'
    params:
        dir='HUMANN'
    output:
        pathcov='HUMANN/{sample}_pathcoverage.tsv',
        pathabund='HUMANN/{sample}_pathabundance.tsv',
        families='HUMANN/{sample}_genefamilies.tsv'
    shell:
        'mkdir -p {params.dir} && cat {input.R1} {input.R2} > merge_{wildcards.sample}_notSyn.fq && '
        'humann --nucleotide-database humann_db --output-basename {wildcards.sample} --input merge_{wildcards.sample}_notSyn.fq --output {params.dir} && '
        'rm merge_{wildcards.sample}_notSyn.fq && '
        'touch {output.pathcov} && touch {output.pathabund} && touch {output.families}'

rule humann_merge_cpm:
    input:expand('HUMANN/{sample}_genefamilies.tsv',sample=SAMPLES)
    conda:'env6.yml'
    output:'HUMANN/merged_genefamilies.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input HUMANN --output {output} --file_name _genefamilies'

rule humann_norm_merged_table:
    input:'HUMANN/merged_genefamilies.tsv'
    conda:'env6.yml'
    output:'HUMANN/merged_genefamilies_cpm.tsv'
    shell:
        'humann_renorm_table --input {input} --output {output}'

rule humann_merge_coverage:
    input:expand('HUMANN/{sample}_pathcoverage.tsv',sample=SAMPLES)
    conda:'env6.yml'
    output:'HUMANN/merged_pathcoverage.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input humann_out --output {output} --file_name pathcoverage'

rule humann_merge_abundance:
    input:expand('HUMANN/{sample}_pathabundance.tsv',sample=SAMPLES)
    conda:'env6.yml'
    output:'HUMANN/merged_pathabundance.tsv'
    shell:
        'touch {input} && '
        'humann_join_table --input HUMANN --output {output} --file_name pathabundance'

#rule norm_not_syn:
#    input:
#        R1='NOTSYN_READS/{sample}_notSyn1.fq',
#        R2='NOTSYN_READS/{sample}_notSyn2.fq'
#    conda:'env1.yml'
#    output:
#        R1='NOTSYN_NORM_READS/{sample}_norm_notSyn1.fq',
#        R2='NOTSYN_NORM_READS/{sample}_norm_notSyn2.fq'
#    shell:
#        'mkdir -p NOTSYN_NORM_READS &&'
#        'bbnorm.sh overwrite=true bits=16 prefilter=t threads=16 -Xmx175g \
#        in={input.R1} in2={input.R2} \
#        out={output.R1} out2={output.R2} target=60 min=5'

#rule combine_norm:
#    input:
#        R1=expand('NOTSYN_NORM_READS/{sample}_norm_notSyn1.fq',sample=SAMPLES),
#        R2=expand('NOTSYN_NORM_READS/{sample}_norm_notSyn2.fq',sample=SAMPLES)
#    output:
#        R1='NOTSYN_NORM_READS/combined_notSyn_R1.fq',
#        R2='NOTSYN_NORM_READS/combined_notSyn_R2.fq'
#    shell:
#        'cat {input.R1} > {output.R1} && cat {input.R2} > {output.R2}'

#rule combine_norm:
#    input:
#        R1=expand('NOTSYN_READS/{sample}_notSyn1.fq',sample=SAMPLES),
#        R2=expand('NOTSYN_READS/{sample}_notSyn2.fq',sample=SAMPLES)
#    output:
#        R1='NOTSYN_READS/combined_notSyn_R1.fq',
#        R2='NOTSYN_READS/combined_notSyn_R2.fq'
#    shell:
#        'cat {input.R1} > {output.R1} && cat {input.R2} > {output.R2}'

#rule not_syn_metaG_megahit:
#    input:
#        R1='NOTSYN_READS/combined_notSyn_R1.fq',
#        R2='NOTSYN_READS/combined_notSyn_R2.fq'
#    conda:'env1.yml'
#    output:'NOTSYN_ASSEMBLIES/noSyn_assembly_megahit.fa'
#    shell:
#        'mkdir -p NOTSYN_ASSEMBLIES'
#        'megahit --mem-flag 2 -1 {input.R1} -2 {input.R2} -o NOTSYN_ASSEMBLIES/MEGAHIT && cp NOTSYN_ASSEMBLIES/MEGAHIT/final.contigs.fa {output} &&'
#        'mkdir -p BINS'

rule combine_not_norm:
    input:
        R1=expand('NOTSYN_READS/{sample}_notSyn1.fq',sample=SAMPLES),
        R2=expand('NOTSYN_READS/{sample}_notSyn2.fq',sample=SAMPLES)
    output:
        R1='combined_notSyn_R1.fq',
        R2='combined_notSyn_R2.fq'
    shell:
        'cat {input.R1} > {output.R1} && cat {input.R2} > {output.R2}'

rule not_syn_metaG_metaspades:
    input:
        R1='combined_notSyn_R1.fq',
        R2='combined_notSyn_R2.fq'
    conda:'env1.yml'
    output:'NOTSYN_ASSEMBLIES/noSyn_assembly_metaspades.fa'
    shell:
        'mkdir -p NOTSYN_ASSEMBLIES && '
        'spades.py --meta --plasmid --phred-offset 33 -1 {input.R1} -2 {input.R2} -o NOTSYN_ASSEMBLIES/METASPADES && cp NOTSYN_ASSEMBLIES/METASPADES/contigs.fasta {output} && '
        'mkdir -p BINS'

# already ran this on the normalized/megahit output, trying on the metaspades output now
rule map_metaG_index:
    input:'NOTSYN_ASSEMBLIES/noSyn_assembly_metaspades.fa'
    conda:'env1.yml'
    params:'NOTSYN_ASSEMBLIES/MS_IND/contigs_meta'
    output:'NOTSYN_ASSEMBLIES/MS_IND/contigs_meta.1.bt2'
    shell:
        'mkdir -p NOTSYN_ASSEMBLIES/MS_IND &&'
        'bowtie2-build {input} {params} && touch {output}'

rule metaG_map:
    input:
        R1='NOTSYN_READS/{sample}_notSyn1.fq',
        R2='NOTSYN_READS/{sample}_notSyn2.fq',
        rdy='NOTSYN_ASSEMBLIES/MS_IND/contigs_meta.1.bt2'
    conda:'env1.yml'
    params:'NOTSYN_ASSEMBLIES/MS_IND/contigs_meta'
    output:
        sam='SAM/{sample}_metaspades.sam'
    shell:
        'mkdir -p SAM && touch {input.rdy} && bowtie2 --threads 16 --very-sensitive-local -x {params} -1 {input.R1} -2 {input.R2} -S {output.sam}'

rule metaG_bam_fix:
    input:'SAM/{sample}_metaspades.sam'
    conda:'env1.yml'
    output:'BAM/{sample}_metaspades.bam'
    shell:
        'mkdir -p BAM && samtools view -b {input} | samtools sort -@ 16 -m 2G - > {output}'

rule metaG_bin_metabat2:
    input:
        bam=expand('BAM/{sample}_metaspades.bam',sample=SAMPLES),
        assembly='NOTSYN_ASSEMBLIES/noSyn_assembly_metaspades.fa'
    conda:'env7.yml'
    output:'BINS/METABAT2/bin.1.fa'
    shell:
        'jgi_summarize_bam_contig_depths --outputDepth depth.txt {input.bam} &&'
        'metabat2 -i {input.assembly} -a depth.txt -o BINS/METABAT/ &&'
        'touch {output}'

rule metaG_bin_index_bams:
    input:'BAM/{sample}_metaspades.bam'
    conda:'env7.yml'
    output:'BAM/{sample}_metaspades.bam.bai'
    shell:
        'samtools index -@ 16 {input} && touch {output}'

rule metaG_bin_concoct:
    input:
        bam=expand('BAM/{sample}_metaspades.bam',sample=SAMPLES),
        bai=expand('BAM/{sample}_metaspades.bam.bai',sample=SAMPLES),
        assembly='NOTSYN_ASSEMBLIES/noSyn_assembly_metaspades.fa'
    conda:'env7.yml'
    output:'BINS/CONCOCT/0.fa'
    shell:
        'mkdir -p BINS/CONCOCT &&'
        'mkdir -p INT/CONCOCT &&'
        'touch {input.bam} &&'
        'touch {input.bai} &&'
        'cut_up_fasta.py {input.assembly} -c 10000 -o 0 --merge_last -b INT/CONCOCT/contigs_10K.bed > INT/CONCOCT/contigs_10K.fa &&'
        'concoct_coverage_table.py INT/CONCOCT/contigs_10K.bed {input.BAM} > INT/CONCOCT/coverage_table.tsv &&'
        'concoct --composition_file INT/CONCOCT/contigs_10K.fa --coverage_file INT/CONCOCT/coverage_table.tsv -b INT/CONCOCT/ &&'
        'merge_cutup_clustering.py INT/CONCOCT/clustering_gt1000.csv > INT/CONCOCT/clustering_merged.csv &&'
        'extract_fasta_bins.py {input.assembly} INT/CONCOCT/clustering_merged.csv --output_path BINS/CONCOCT/'
        'touch {output}'
