#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/UMInator
========================================================================================
 MaestSi/UMInator analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/UMInator
----------------------------------------------------------------------------------------
*/
def helpMessage() {
        log.info"""
    Usage:
    nextflow -c UMInator.conf run UMInator.nf --fastq_files = "/path/to/files*.fastq" --scripts_dir = "/path/to/scripts_dir" --results_dir = "/path/to/results_dir" -profile docker

    Mandatory argument:
    -profile                                                              Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the UMInator.conf file
    --fastq_files                                                         Path to fastq files, use wildcards to select multiple samples
    --results_dir                                                         Path to a folder where to store results
    --scripts_dir                                                         Directory containing all scripts
    --FW_adapter                                                          Forward adapter sequence
    --RV_adapter                                                          Reverse adapter sequence
    --FW_primer                                                           Forward primer sequence
    --RV_primer                                                           Reverse primer sequence
    --minQ                                                                min Q value for reads filtering
    --minLen                                                              min read length for reads filtering
    --maxLen                                                              max read length for reads filtering
    --searchLen                                                           Amount of bases at the beginning and end of each read to search for UMIs
    --tolCutadaptErr                                                      Cutadapt maximum allowed error rate [0, 1]
    --minLenOvlp                                                          Min overlap between read and adapter
    --UMILen                                                              UMI length (before merging UMI1 and UMI2)
    --UMIPattern                                                          UMI structure (after merging UMI1 and UMI2) in the form of a regex of the type: [nucl.]{cardinality}
    --UMIClustID                                                          UMI clustering identity
    --maxDiff                                                             BWA aln maximum number of differences 
    --max_NM_mean                                                         Maximum tolerated mean mapping difference between UMI and query reads
    --max_NM_sd                                                           Maximum tolerated sd of mapping difference between UMI and query reads
    --min_UMI_freq                                                        Minimum number of reads assigned to UMI for generating a consensus sequence
    --target_reads_consensus                                              Maximum number of reads used for consensus calling
    --target_reads_polishing                                              Maximum number of reads used for consensus polishing
    --fast_consensus_flag                                                 Set fast_consensus_flag = 1 for obtaining draft consensus with VSEARCH, instead of  MAFFT + EMBOSS cons
    --fast_polishing_flag                                                 Set fast_polishing_flag = 1 for polishing with Racon, instead of Racon + Medaka
    --plurality                                                           MAFFT plurality value: minimum fraction of aligned reads supporting a basis for including it in the preliminary consensus
    --fast_alignment_flag                                                 Set fast_alignment_flag=1 if you want to perform fast multiple sequence alignment; otherwise set fast_alignment_flag=0
    --medaka_model                                                        Medaka model for consensus polishing
    --maxF                                                                Maximum forks
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of sample names and conditions
Channel
    .fromPath(params.fastq_files)
    .map {tuple( it.name.split('\\.')[0], it )}
    .set{inputFiles_readsFiltering}

// readsFiltering
process readsFiltering {
  input:
    tuple val(sample), val(fastq)
  output:
    val sample
  script:
  if(params.readsFiltering)
  """
    mkdir -p ${params.results_dir}
    mkdir -p ${params.results_dir}/readsFiltering
    mkdir -p ${params.results_dir}/readsFiltering/${sample}

    #filter reads
    cat ${fastq} | NanoFilt -q ${params.minQ} --length ${params.minLen} --maxlength ${params.maxLen} > ${params.results_dir}/readsFiltering/${sample}/${sample}_filtered.fastq
  """
  else
  """
    cp ${fastq} ${params.results_dir}/readsFiltering/${sample}/${sample}_filtered.fastq
  """
}
    
// candidateUMIsExtraction
process candidateUMIsExtraction {
  input:
    val sample
  output:
    val sample
  script:
  if(params.candidateUMIsExtraction)
  """
    mkdir -p ${params.results_dir}
    mkdir -p ${params.results_dir}/candidateUMIsExtraction
    mkdir -p ${params.results_dir}/candidateUMIsExtraction/${sample}

    #obtain reverse complement of primers and adapters sequences
    FW_primer_R=\$(echo -e \">tmp\\n\" ${params.FW_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_primer_R=\$(echo -e \">tmp\\n\" ${params.RV_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    FW_adapter_R=\$(echo -e \">tmp\\n\" ${params.FW_adapter} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_adapter_R=\$(echo -e \">tmp\\n\" ${params.RV_adapter} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')

    #double UMI design
    #obtain first and last bases of reads
    READS_START=${params.results_dir}/candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
    READS_END=${params.results_dir}/candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fastq
    seqtk trimfq -L ${params.searchLen} ${params.results_dir}/readsFiltering/${sample}/${sample}_filtered.fastq > \$READS_START
    seqtk seq -r ${params.results_dir}/readsFiltering/${sample}/${sample}_filtered.fastq | seqtk trimfq -L ${params.searchLen} - | seqtk seq -r - > \$READS_END
    
    #search candidate UMIs with exact length between adapters and primers
    cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m ${params.UMILen} -M ${params.UMILen} \
    --discard-untrimmed --match-read-wildcards \
    -g ${params.FW_adapter}...${params.FW_primer} -g ${params.RV_adapter}...${params.RV_primer} \
    -G \$RV_primer_R...\$RV_adapter_R -G \$FW_primer_R...\$FW_adapter_R \
    -o ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_db_tmp1.fastq \
    -p ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_db_tmp1.fastq \
    \$READS_START \$READS_END

    #collapse the 5' and 3' end UMIs of each read
    paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_db_tmp1.fastq ) \
    <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_db_tmp1.fastq ) | cut -d " " -f1 \
    > ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_db_tmp1.fasta

    #search candidate UMIs with approximate length between adapters and primers
    cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m ${params.UMILen} -l ${params.UMILen} \
    --discard-untrimmed --match-read-wildcards \
    -g ${params.FW_adapter} -g ${params.RV_adapter} \
    -G \$RV_primer_R -G \$FW_primer_R \
    -o ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq \
    -p ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq \
    \$READS_START \$READS_END

    #collapse the 5' and 3' end UMIs of each read
    paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq ) \
    <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq ) |  \
    cut -d " " -f1 > ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fastq
      
    #convert UMI candidates to fasta
    seqtk seq -A ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fastq > ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fasta
  """
  else
  """
    echo "Skipped"
  """
}

//candidateUMIsFiltering
process candidateUMIsFiltering {
  input:
    val(sample)
  output:
    val(sample)
  script:
  if(params.candidateUMIsFiltering)
  """
    mkdir -p ${params.results_dir}/candidateUMIsFiltering
    mkdir -p ${params.results_dir}/candidateUMIsFiltering/${sample}
    
    #filter candidate UMIs for pattern
    cat ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_db_tmp1.fasta | grep -B1 -E ${params.UMIPattern} | sed \'/^--\$/d\' \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp2.fasta

    #evaluate the total UMI length (UMI1 + UMI2)
    totUMILen=\$(echo 2*${params.UMILen} | bc)
    
    #dereplicate high-quality UMIs
    vsearch --derep_fulllength ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp2.fasta \
    --minseqlength \$totUMILen --sizeout --relabel umi --strand both \
    --output ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp3.fasta

    #cluster dereplicated high-quality UMIs to obtain a database of high-quality UMIs
    vsearch --cluster_smallmem ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp3.fasta \
    --usersort --id ${params.UMIClustID} --iddef 1 --strand both --clusterout_sort --minseqlength \$totUMILen \
    --sizein --sizeout --consout ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    --minwordmatches 0

    #index the database of tmp high-quality UMIs
    bwa index ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta

    #align candidate UMIs to tmp high-quality UMIs
    bwa aln \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fasta \
    -n ${params.maxDiff} \
    -t ${task.cpus} \
    -N > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sai

    bwa samse \
    -n 10000000 \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sai \
    ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fasta | \
    samtools view -F 4 - \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sam

    #retain only UMIs from the database that are supported by at least 3 reads
    cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sam \
    | cut -f3 | sort | uniq -c | sort -nr | awk \'BEGIN {FS=\" \"} { if ( \$1 > 2 ) print "\t"\$2"\t"}\' \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3.txt

    #extract from fasta only UMIs supported by at least 3 reads
    cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3.txt \
    | cut -f2 > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3_noTab.txt

    seqtk subseq ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3_noTab.txt \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta

    #split the database of high-quality UMIs in part 1 and part 2
    #trim UMILen from right
    seqtk trimfq -e ${params.UMILen} ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta

    seqtk seq -r ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta \
    | sed \'/^>/ s/\$/_rc/\' > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1_rv.fasta

    #trim UMILen from left
    seqtk trimfq -b ${params.UMILen} ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta

    seqtk seq -r ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta \
    | sed \'/^>/ s/\$/_rc/\' > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2_rv.fasta
    
    #concatenate UMI1 and reverse complement of UMI2
    cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2_rv.fasta \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.fasta

    #concatenate UMI2 and reverse complement of UMI1
    cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1_rv.fasta \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.fasta

    #align high-quality UMIs (part 1 and part 2) to the terminal portion of reads
    READS_START_FQ=${params.results_dir}/candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
    READS_START_FA=${params.results_dir}/candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fasta
    READS_END_FQ=${params.results_dir}/candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fastq
    READS_END_FA=${params.results_dir}/candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fasta
    
    #convert fastq to fasta
    seqtk seq -A \$READS_START_FQ > \$READS_START_FA
    seqtk seq -A \$READS_END_FQ > \$READS_END_FA

    #index db
    bwa index \$READS_START_FA
    bwa index \$READS_END_FA

    maxDiffSingle=\$(echo ${params.maxDiff}/2 | bc)
    
    #map UMI_db_p1 to reads start
    bwa aln \
    \$READS_START_FA \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.fasta \
    -n \$maxDiffSingle \
    -t ${task.cpus} \
    -N > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.sai

    bwa samse \
    -n 10000000 \
    \$READS_START_FA \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.sai \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.fasta | \
    samtools view -F 20 - \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.sam

    #map UMI_db_p2 to reads end
    bwa aln \
    \$READS_END_FA \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.fasta \
    -n \$maxDiffSingle \
    -t ${task.cpus} \
    -N > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.sai

    bwa samse \
    -n 10000000 \
    \$READS_END_FA \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.sai \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.fasta | \
    samtools view -F 20 - \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.sam

    #filter alignments
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Filter_UMIs.R alignment_file_1=${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p1.sam alignment_file_2=${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_p2.sam map_file=${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_read_map.txt max_NM_mean=${params.max_NM_mean} max_NM_sd=${params.max_NM_sd}
  """
  else
  """
    echo "Skipped"
  """
}

//readsUMIsAssignment
process readsUMIsAssignment {
  maxForks params.maxF
  input:
    val sample
  output:
    tuple val(sample), env(UMI_all)
  script:
  if(params.readsUMIsAssignment)
  """
    mkdir -p ${params.results_dir}/readsUMIsAssignment
    mkdir -p ${params.results_dir}/readsUMIsAssignment/${sample}

    #read the current fastq file and, if a read matches a single UMI, assign it to it
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Bin_reads.R fastq_file=${params.results_dir}/readsFiltering/${sample}/${sample}_filtered.fastq map_file=${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_read_map.txt outdir=${params.results_dir}/readsUMIsAssignment/${sample}

    #extract all UMIs
    UMI_all_tmp=\$(cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_read_map.txt | cut -f2 | sed \'s/^centroid=//\' | sed \'s/;seqs=.*//\' | sort | uniq)

    #add sample name before UMI ID
    sn=${sample}
    UMI_all=\$(echo \$UMI_all_tmp | sed \"s/umi/\$sn\"\\|umi\"/g\")
  """
  else
  """
    echo "Skipped"
  """
}

//draftConsensusCalling
process draftConsensusCalling {
  maxForks params.maxF
  input:
    tuple val(sample), val(UMI)
  output:
    tuple val(sample), val(UMI)
  script:
  if(params.draftConsensusCalling)
  """
    mkdir -p ${params.results_dir}/draftConsensusCalling
    mkdir -p ${params.results_dir}/draftConsensusCalling/${sample}

    #concatenate files with the same UMI obtained from different reads chunks
    chunks_binned_files_curr_UMI=\$(find ${params.results_dir}/readsUMIsAssignment/${sample}/ | grep ${UMI}_chunk);
    cat \$chunks_binned_files_curr_UMI > ${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq
    rm \$chunks_binned_files_curr_UMI

    #convert fastq to fasta
    seqtk seq -A ${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq > ${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fasta

    #obtain draft consensus sequence
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Obtain_draft_consensus.R fastq_file=${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq TRC=${params.target_reads_consensus} PLUR=${params.plurality} num_threads=${task.cpus} fast_alignment_flag=${params.fast_alignment_flag} fast_consensus_flag=${params.fast_consensus_flag} min_UMI_freq=${params.min_UMI_freq} 
  """
  else
  """
    echo "Skipped"
  """
}

//QC
process QC {
  input:
    val sample
  output:
  script:
  if(params.QC)
  """
    mkdir -p ${params.results_dir}/QC
    mkdir -p ${params.results_dir}/QC/${sample}
    #concatenate files with the same UMI obtained from different reads chunks
    chunks_unbinned_files=\$(find ${params.results_dir}/readsUMIsAssignment/${sample}/ -name \"*unbinned_chunk*\");
    unbinned_reads_files=${params.results_dir}/readsUMIsAssignment/${sample}/unbinned.fastq
    binned_reads_files=${params.results_dir}/readsUMIsAssignment/${sample}/binned.fastq
    if [[ -f "\$chunks_unbinned_files" ]]; then cat \$chunks_unbinned_files > \$unbinned_reads_files; rm \$chunks_unbinned_files; fi
    fastq_files_binned=\$(find ${params.results_dir}/readsUMIsAssignment/${sample} -name \"umi*\\.fastq\")
    fastq_files=\$(find ${params.results_dir}/readsUMIsAssignment/${sample} -name \"*\\.fastq\")
    if [[ ! -z "\$fastq_files_binned" ]]; then
      #cat \$fastq_files_binned > \$binned_reads_files
      for f in \$fastq_files_binned; do
        cat \$f >> \$binned_reads_files;
      done
    fi
    
    #do QC plot for unbinned reads
    if [[ -f "\$unbinned_reads_files" ]]; then
      NanoPlot -t ${task.cpus} --fastq \$unbinned_reads_files -o ${params.results_dir}/QC/${sample}/QC_unbinned_reads
    fi

    #do QC plot for binned reads
    if [[ -f "\$binned_reads_files" ]]; then
      NanoPlot -t ${task.cpus} --fastq \$binned_reads_files -o ${params.results_dir}/QC/${sample}/QC_binned_reads 
      rm \$binned_reads_files
    fi
    
    #produce tsv files with read-UMI assignment stats
    for f in \$fastq_files; do
      reads_names=\$(seqtk seq -A \$f | grep \"^>\" | sed \'s/>//\' | paste -sd ",")
      num_reads=\$(seqtk seq -A \$f | grep \"^>\" | sed \'s/>//\' | wc -l)
      echo -e \"\$(basename \$f)\t\$num_reads\t\$reads_names\" >> ${params.results_dir}/QC/${sample}/${sample}_UMI_stats_tmp.tsv
    done
    
    cat ${params.results_dir}/QC/${sample}/${sample}_UMI_stats_tmp.tsv | sort -k2,2nr > ${params.results_dir}/QC/${sample}/${sample}_UMI_stats.tsv

    rm ${params.results_dir}/QC/${sample}/${sample}_UMI_stats_tmp.tsv
  """
  else
  """
    echo "Skipped"
  """
}

//consensusPolishing
process consensusPolishing {
  maxForks params.maxF
  input:
    tuple val(sample), val(UMI)
  output:
    tuple val(sample), val(UMI)
  script:
  if(params.consensusPolishing)
  """
    mkdir -p ${params.results_dir}/consensusPolishing
    mkdir -p ${params.results_dir}/consensusPolishing/${sample}
    
    #polish consensus sequence with racon and medaka
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Polish_consensus.R draft_consensus=${params.results_dir}/draftConsensusCalling/${sample}/${UMI}/${UMI}_draft_consensus.fasta fastq_file=${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq TRP=${params.target_reads_polishing}  num_threads=${task.cpus} fast_polishing_flag=${params.fast_polishing_flag} medaka_model=${params.medaka_model}
  """
  else
  """
    cp ${params.results_dir}/draftConsensusCalling/${sample}/${UMI}/${UMI}_draft_consensus.fasta ${params.results_dir}/consensusPolishing/${sample}/${UMI}/${UMI}_polished_consensus.fasta
    
    echo "Skipped"
  """
}

//primers trimming
process primersTrimming {
  input:
    val(sample)
  output:
  script:
  if(params.primersTrimming)
  """
    mkdir -p ${params.results_dir}/primersTrimming
    mkdir -p ${params.results_dir}/primersTrimming/${sample}

    #find consensus sequences for all UMIs of one sample and concatenate them
    polished_consensus=\$(find ${params.results_dir}/consensusPolishing/${sample} | grep "_polished_consensus.fasta")
    cat \$polished_consensus > ${params.results_dir}/primersTrimming/${sample}/${sample}_consensus_polished.fasta

    #obtain reverse complement of primers and adapters sequences
    FW_primer_R=\$(echo -e \">tmp\\n\" ${params.FW_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_primer_R=\$(echo -e \">tmp\\n\" ${params.RV_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
  
    #trim PCR primers and external sequence
    if [[ -f "${params.results_dir}/primersTrimming/${sample}/${sample}_consensus_polished.fasta" ]]; then
      cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} \
      --discard-untrimmed --match-read-wildcards \
      -g ${params.FW_primer} -g \$RV_primer_R  \
      -a ${params.RV_primer} -a \$FW_primer_R \
      -o ${params.results_dir}/primersTrimming/${sample}/${sample}_consensus_polished_primersTrimmed.fasta \
      ${params.results_dir}/primersTrimming/${sample}/${sample}_consensus_polished.fasta
    fi
  """
  else
  """
    mkdir -p ${params.results_dir}/primersTrimming
    mkdir -p ${params.results_dir}/primersTrimming/${sample}

    #find consensus sequences for all UMIs of one sample and concatenate them
    polished_consensus=\$(find ${params.results_dir}/consensusPolishing/${sample} | grep "_polished_consensus.fasta")
    cat \$polished_consensus > ${params.results_dir}/primersTrimming/${sample}/${sample}_consensus_polished.fasta
  """
}

workflow {
  //filter reads
  readsFiltering(inputFiles_readsFiltering)

  //extract candidate UMIs
  candidateUMIsExtraction(readsFiltering.out)
  
  //filter candidate UMIs to buid db of high-confidence UMIs
  candidateUMIsFiltering(candidateUMIsExtraction.out)
  
  //assign reads to UMIs
  readsUMIsAssignment(candidateUMIsFiltering.out)

  //analyze binned reads with the same UMI in parallel: obtain Sample-UMI tuples for each reads chunk
  readsUMIsAssignment.out
  .groupTuple(by:0)
  .map { it -> it[1]}
  .set{UMIs_tmp}

  //obtain unique couples of sampleName-UMI for each sample
  UMIs_tmp
  .flatten()
  .distinct()
  .splitCsv( sep: ' ')
  .distinct()
  .flatten()
  .distinct()
  .splitCsv( sep: '|')
  .set{SN_UMIs}

  //perform draft consensus calling
  draftConsensusCalling(SN_UMIs)

  //group the output of consensus sequences by sample
  draftConsensusCalling.out
  .groupTuple(by:0)
  .map{ it -> it[0]}
  .set{samplesQC}

  //run QC
  QC(samplesQC)
  
  //perform consensus polishing
  consensusPolishing(draftConsensusCalling.out)
  
  //perform primers trimming
  primersTrimming(consensusPolishing.out.collect().map{ it -> it[0]})
}
