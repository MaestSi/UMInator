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
    --UMI design                                                          "double" or "single", depending on whether reads have UMIs at both ends or not
    --FW_adapter                                                          Forward adapter sequence
    --RV_adapter                                                          Reverse adapter sequence
    --FW_primer                                                           Forward primer sequence
    --RV_primer                                                           Reverse primer sequence
    --searchLen                                                           Amount of bases at the beginning and end of each read to search for UMIs
    --tolCutadaptErr                                                      Cutadapt maximum allowed error rate [0, 1]
    --minLenOvlp                                                          Min overlap between read and adapter
    --UMILen                                                              UMI length (before merging UMI1 and UMI2 in case of double UMI design)
    --UMILenTol                                                           Tolerated candidate UMI discrepancy in length
    --UMIPattern                                                          UMI structure (after merging UMI1 and UMI2, in case of double UMI design) in the form of a regex of the type: [nucl.]{cardinality}
    --UMIClustID                                                          UMI clustering identity
    --seedLen                                                             BWA seed length
    --readsChunkSize                                                      Number of lines in each fastq split file (should be multiple of 4)
    --min_UMI_freq                                                        Minimum number of reads assigned to UMI for generating a consensus sequence
    --target_reads_consensus                                              Maximum number of reads used for consensus calling
    --target_reads_polishing                                              Maximum number of reads used for consensus polishing
    --plurality                                                           MAFFT plurality value: minimum fraction of aligned reads supporting a basis for including it in the preliminary consensus
    --fast_alignment_flag                                                 Set fast_alignment_flag=1 if you want to perform fast multiple sequence alignment; otherwise set fast_alignment_flag=0
    --medaka_model                                                        Medaka model for consensus polishing
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
    .set{inputFiles_candidateUMIsExtraction}
    
// candidateUMIsExtraction
process candidateUMIsExtraction {
  input:
    tuple val(sample), val(fastq)
  output:
    tuple val(sample), val(fastq)
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
    if [[ ${params.UMIDesign} == "double" ]]; then
      #obtain first and last bases of reads
      READS_START=${params.results_dir}/candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
      READS_END=${params.results_dir}/candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fastq
      seqtk trimfq -L ${params.searchLen} ${fastq} > \$READS_START
      seqtk seq -r ${fastq} | seqtk trimfq -L ${params.searchLen} - | seqtk seq -r - > \$READS_END
    
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

      #evaluate candidate UMI min and max length
      UMIMinLen=\$(echo ${params.UMILen} - ${params.UMILenTol} | bc)
      UMIMaxLen=\$(echo ${params.UMILen} + ${params.UMILenTol} | bc)

      #search candidate UMIs with approximate length between adapters and primers
      cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m \$UMIMinLen -M \$UMIMaxLen \
      --discard-untrimmed --match-read-wildcards \
      -g ${params.FW_adapter}...${params.FW_primer} -g ${params.RV_adapter}...${params.RV_primer} \
      -G \$RV_primer_R...\$RV_adapter_R -G \$FW_primer_R...\$FW_adapter_R \
      -o ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq \
      -p ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq \
      \$READS_START \$READS_END

      #collapse the 5' and 3' end UMIs of each read
      paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq ) \
      <( sed -n '1~4s/^@/>/p;2~4p' ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq ) |  \
      cut -d " " -f1 > ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fastq
      
    #single UMI design
    else
      #obtain first bases at both ends of reads
      READS_START=${params.results_dir}/candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
      seqtk trimfq -L ${params.searchLen} ${fastq} > \$READS_START
      seqtk seq -r ${fastq} | seqtk trimfq -L ${params.searchLen} | sed s\'/>/>RC_/\' - >> \$READS_START
    
      #search candidate UMIs with exact length between adapters and primers
      cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m ${params.UMILen} -M ${params.UMILen} \
      --discard-untrimmed --match-read-wildcards \
      -g ${params.FW_adapter}...${params.FW_primer} \
      -o ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_db_tmp1.fastq \
      \$READS_START

      #evaluate candidate UMI min and max length
      UMIMinLen=\$(echo ${params.UMILen} - ${params.UMILenTol} | bc)
      UMIMaxLen=\$(echo ${params.UMILen} + ${params.UMILenTol} | bc)

      #search candidate UMIs with approximate length between adapters and primers
      cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m \$UMIMinLen -M \$UMIMaxLen \
      --discard-untrimmed --match-read-wildcards \
      -g ${params.FW_adapter}...${params.FW_primer} \
      -o ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fastq \
      \$READS_START 
    fi
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
    tuple val(sample), val(fastq)
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
    if [[ ${params.UMIDesign} == "double" ]]; then
      totUMILen=\$(echo 2*${params.UMILen} | bc)
    else
      totUMILen=${params.UMILen}
    fi

    #dereplicate high-quality UMIs
    vsearch --derep_fulllength ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp2.fasta \
    --minseqlength \$totUMILen --sizeout --relabel umi --strand both \
    --output ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp2_derep.fasta

    #cluster dereplicated high-quality UMIs to obtain a database of high-quality UMIs
    vsearch --cluster_smallmem ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db_tmp2_derep.fasta \
    --usersort --id ${params.UMIClustID} --iddef 1 --strand both --clusterout_sort --minseqlength \$totUMILen \
    --sizein --sizeout --consout ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta \
    --minwordmatches 0

    #index the database of high-quality UMIs
    bwa index ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta

    #align candidate UMIs to high-quality UMIs
    bwa aln \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta \
    ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fasta \
    -l ${params.seedLen} \
    -t ${task.cpus} \
    -N > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map.sai

    bwa samse \
    -n 10000000 \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_db.fasta \
    ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map.sai \
    ${params.results_dir}/candidateUMIsExtraction/${sample}/UMI_candidates.fasta | \
    samtools view -F 4 - \
    > ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map.sam
  """
  else
  """
    echo "Skipped"
  """
}

//readsUMIsAssignment
process readsUMIsAssignment {
  input:
    val sample
    each path('readsChunk.fastq')
  output:
    tuple val(sample), env(UMI_all)
  script:
  if(params.readsUMIsAssignment)
  """
    mkdir -p ${params.results_dir}/readsUMIsAssignment
    mkdir -p ${params.results_dir}/readsUMIsAssignment/${sample}
    
    readsCurrChunk=\$(basename \$(realpath readsChunk.fastq))

    #read the current chunk of reads and, if a read matches to a single UMI, assign it to it
    cat readsChunk.fastq | while read -r l; do
      read_id=\$(echo \$l | sed \'s/^@//\' | cut -d \' \' -f1);
      UMI_curr_chunk=\$(cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map.sam | cut -f1,3 | grep -P "\$read_id\t" | cut -f2 | sed \'s/^centroid=//\' | sed \'s/;seqs=.*//\');
      Matches=\$(echo \$UMI_curr_chunk | grep \"umi\" -o | wc -l);
      if [ -z \"\$UMI_curr_chunk\" ] || [ \"\$Matches\" -gt 1 ]; then UMI_curr_chunk=unbinned; fi;
        echo @\$read_id >> ${params.results_dir}/readsUMIsAssignment/${sample}/\$UMI_curr_chunk\"_chunk_\"\$readsCurrChunk
        read -r l; printf \'%s\n\' \$l >> ${params.results_dir}/readsUMIsAssignment/${sample}/\$UMI_curr_chunk\"_chunk_\"\$readsCurrChunk
        read -r l; printf \'%s\n\' \$l >> ${params.results_dir}/readsUMIsAssignment/${sample}/\$UMI_curr_chunk\"_chunk_\"\$readsCurrChunk
        read -r l; printf \'%s\n\' \$l >> ${params.results_dir}/readsUMIsAssignment/${sample}/\$UMI_curr_chunk\"_chunk_\"\$readsCurrChunk
      done
      #extract all UMIs
      UMI_all_tmp=\$(cat ${params.results_dir}/candidateUMIsFiltering/${sample}/UMI_candidates_map.sam | cut -f3 | sed \'s/^centroid=//\' | sed \'s/;seqs=.*//\' | sort | uniq)

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
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Obtain_draft_consensus.R fastq_file=${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq TRC=${params.target_reads_consensus} PLUR=${params.plurality} num_threads=${task.cpus} fast_alignment_flag=${params.fast_alignment_flag} min_UMI_freq=${params.min_UMI_freq}
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
    val sample
  script:
  if(params.QC)
  """
    mkdir -p ${params.results_dir}/QC
    mkdir -p ${params.results_dir}/QC/${sample}
    mkdir -p ${params.results_dir}/QC/${sample}/QC_all_reads
    mkdir -p ${params.results_dir}/QC/${sample}/QC_binned_reads

    #concatenate files with the same UMI obtained from different reads chunks
    chunks_unbinned_files=\$(find ${params.results_dir}/readsUMIsAssignment/${sample}/ | grep unbinned_chunk);
    cat \$chunks_unbinned_files > ${params.results_dir}/readsUMIsAssignment/${sample}/unbinned.fastq
    rm \$chunks_unbinned_files
    
    #do QC plot for all reads
    fastq_files=\$(find ${params.results_dir}/readsUMIsAssignment/${sample} | grep \"\\.fastq\")
    NanoPlot -t ${task.cpus} --fastq \$fastq_files -o ${params.results_dir}/QC/${sample}/QC_all_reads

    #do QC plot for all reads
    fastq_files_binned=\$(find ${params.results_dir}/readsUMIsAssignment/${sample} | grep \"umi.*\\.fastq\")
    NanoPlot -t ${task.cpus} --fastq \$fastq_files_binned -o ${params.results_dir}/QC/${sample}/QC_binned_reads 
    
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
    /opt/conda/envs/UMInator_env/bin/Rscript ${params.scripts_dir}/Polish_consensus.R draft_consensus=${params.results_dir}/draftConsensusCalling/${sample}/${UMI}/${UMI}_draft_consensus.fasta fastq_file=${params.results_dir}/readsUMIsAssignment/${sample}/${UMI}.fastq TRP=${params.target_reads_polishing}  num_threads=${task.cpus} medaka_model=${params.medaka_model}
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
    val(sample)
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
    echo "Skipped"
  """
}

workflow {
  //extract candidate UMIs
  candidateUMIsExtraction(inputFiles_candidateUMIsExtraction)
  
  //filter candidate UMIs to buid db of high-confidence UMIs
  candidateUMIsFiltering(candidateUMIsExtraction.out)
  
  //split reads into chunks
  Channel
  .fromPath(params.fastq_files)
  .splitFastq( by: params.readsChunkSize, file:true )
  .set{readsChunk}
  
  //assign chunks of reads to UMIs
  readsUMIsAssignment(candidateUMIsFiltering.out, readsChunk)

  //analyze binned reads with the same UMI in parallel: obtain Sample-UMI tuples for each reads chunk
  readsUMIsAssignment.out
  .groupTuple(by:0)
  .multiMap { it ->
    SN: it[0]
    UMI: it[1]}
  .set{UMIs_tmp}

  //obtain unique couples of sampleName-UMI for each sample
  UMIs_tmp.UMI
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
