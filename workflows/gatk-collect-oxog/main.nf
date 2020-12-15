#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
name = 'gatk-collect-oxog'
version = '1.0.0'

/*
==============================================================================================
                        ICGC-ARGO Data Quality Assurance Workflows - Collect OxoG Workflow
==============================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/data-qa-wfs/workflows/gatk-collect-oxog
#### Authors
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              song study ID
--analysis_id                           song sequencing_experiment analysis ID
--ref_genome_fa                         reference genome '.fa' file, other secondary files are expected to be under the same folder
--song_url                              song server URL
--score_url                             score server URL
--api_token                             song/score API Token

General Parameters (with defaults):
--cpus                                  cpus given to all process containers (default 1)
--mem                                   memory (GB) given to all process containers (default 1)

Download Parameters (object):
--download
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for download process (defaults to main song_url param)
    score_url                           score url for download process (defaults to main score_url param)
    api_token                           song/score API token for download process (defaults to main api_token param)
    song_cpus
    song_mem
    score_cpus
    score_mem
    score_transport_mem                 TODO: Description
}

gatkCollectOxogMetrics (object):
--gatkCollectOxogMetrics
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for align container, defaults to cpus parameter
    mem                                 memory (GB) for align container, defaults to mem parameter
    oxog_scatter                        number of parallel tasks for scattering OxoG metrics collection
}

Upload Parameters (object):
--upload
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for upload process (defaults to main song_url param)
    score_url                           score url for upload process (defaults to main score_url param)
    api_token                           song/score API token for upload process (defaults to main api_token param)
    song_cpus                           cpus for song container, defaults to cpus parameter
    song_mem                            memory (GB) for song container, defaults to mem parameter
    score_cpus                          cpus for score container, defaults to cpus parameter
    score_mem                           memory (GB) for score container, defaults to mem parameter
    score_transport_mem                 memory (GB) for score_transport, defaults to mem parameter
    extract_cpus                        cpus for extract container, defaults to cpus parameter
    extract_mem                         memory (GB) extract score container, defaults to mem parameter
}

*/

params.study_id = ""
params.analysis_id = ""
params.ref_genome_fa = ""
params.cleanup = true
params.cpus = 1
params.mem = 1
params.tempdir = "NO_DIR"
params.publish_dir = ""
params.analysis_metadata = "NO_FILE"
params.sequencing_files = []
params.song_url = ""
params.score_url = ""
params.api_token = ""
params.download = [:]
params.payloadGenDnaSeqQc = [:]
params.uploadQc = [:]
params.gatkCollectOxogMetrics = [:]


download_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

payloadGenDnaSeqQc_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenDnaSeqQc ?: [:])
]

uploadQc_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.uploadQc ?: [:])
]

gatkCollectOxogMetrics_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'oxog_scatter': 8,  // default, may be overwritten by params file
    *:(params.gatkCollectOxogMetrics ?: [:])
]


// Include all modules and pass params
include { songScoreDownload as dnld } from './song-score-utils/song-score-download' params(download_params)
include { payloadGenDnaSeqQc as pGenDnaSeqQc } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-dna-seq-qc.0.5.3.0/tools/payload-gen-dna-seq-qc/payload-gen-dna-seq-qc.nf" params(payloadGenDnaSeqQc_params)
include { gatkSplitIntervals as splitItvls; getSecondaryFiles as getSIIdx } from "./modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.4.1-1.0/tools/gatk-split-intervals/gatk-split-intervals"
include { gatkCollectOxogMetrics as oxog; getOxogSecondaryFiles; gatherOxogMetrics as gatherOM } from "./modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-collect-oxog-metrics.4.1.8.0-2.0/tools/gatk-collect-oxog-metrics/gatk-collect-oxog-metrics" params(gatkCollectOxogMetrics_params)
include { songScoreUpload as upQc } from './song-score-utils/song-score-upload' params(uploadQc_params)
include { cleanupWorkdir as cleanup } from './modules/raw.githubusercontent.com/icgc-argo/nextflow-data-processing-utility-tools/2.3.0/process/cleanup-workdir'


workflow CollectOxog {
    take:
        study_id
        analysis_id
        ref_genome_fa
        analysis_metadata
        sequencing_files

    main:
        // detect local mode or not
        local_mode = false
        if (analysis_metadata != "NO_FILE" && sequencing_files.size() > 0){
            local_mode = true
            if (params.cleanup == true) {
                log.info "Will not perform 'cleanup' step when running in local mode."
            }
            analysis_metadata = file(analysis_metadata)
            sequencing_files = Channel.fromPath(sequencing_files)
        } else {
            // download files and metadata from song/score (analysis type: sequencing_experiment)
            dnld(study_id, analysis_id)
            analysis_metadata = dnld.out.song_analysis
            sequencing_files = dnld.out.files
        }

        // prepare oxog_scatter intervals
        splitItvls(gatkCollectOxogMetrics_params.oxog_scatter, file(ref_genome_fa),
            Channel.fromPath(getSIIdx(ref_genome_fa), checkIfExists: true).collect(), file('NO_FILE'))

        // perform gatkCollectOxogMetrics in parallel tasks
        oxog(sequencing_files.flatten().first(), sequencing_files.flatten().last(),
            file(ref_genome_fa),
            Channel.fromPath(getOxogSecondaryFiles(ref_genome_fa), checkIfExists: true).collect(),
            splitItvls.out.interval_files.flatten(), true)  // run after alignedSeqQC

        // gatherOxogMatrics
        gatherOM(oxog.out.oxog_metrics.collect())

        // prepare song payload for qc metrics
        pGenDnaSeqQc(analysis_metadata,
                gatherOM.out.oxog_metrics.collect(),
            name, version)

        // upload aligned file and metadata to song/score
        if (!local_mode) {
            upQc(study_id, pGenDnaSeqQc.out.payload, pGenDnaSeqQc.out.qc_files.collect())
        }

        if (params.cleanup && !local_mode) {
            cleanup(
                sequencing_files.concat(oxog.out).collect(),
                upQc.out.analysis_id.collect())  // wait until upAln and upQc is done
        }

    emit:
        qc_metrics_payload = pGenDnaSeqQc.out.payload
        qc_metrics_files = pGenDnaSeqQc.out.qc_files
}


workflow {
    CollectOxog(
        params.study_id,
        params.analysis_id,
        params.ref_genome_fa,
        params.analysis_metadata,
        params.sequencing_files
    )
}
