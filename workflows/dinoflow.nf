/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDinoflow.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
ch_star_index = Channel.fromPath(params.star_index).map{ fn -> [ [id: fn.baseName], fn ]}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { EXTRACT_BARCODE_TXT } from '../modules/local/extract_barcode_txt/main'
include { STARSOLO } from '../modules/nf-core/star/starsolo/main'


/*
HELPER FUNCTIONS

*/

def extract_csv(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, samplesheet_line_count = 0;
        while ((line = reader.readLine()) != null) {samplesheet_line_count++}
        if (samplesheet_line_count < 2) {
            error("Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines.")
        }
    }
    // Additional check of sample sheet:
    // 1. If params.step == "mapping", then each row should specify a lane and the same combination of pool, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different pools.
    def pool_lane_combinations = []
    def sample2pool = [:]

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
            if ( !row.lane ) {  // This also handles the case where the lane is left as an empty string
                error('The sample sheet should specify a lane for pool "' + row.pool.toString() + '".')
            }
            def pool_lane = [row.pool.toString(), row.lane.toString()]
            if (pool_lane in pool_lane_combinations) {
                error('The pool-lane combination "' + row.pool.toString() + '", and "' + row.lane.toString() + '" is present multiple times in the sample sheet.')
            } else {
                pool_lane_combinations.add(pool_lane)
            }
        }
    sample_count_all = 0
    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by pool and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            if (!(row.anno && row.fastq_1 && row.fastq_2 && row.lane)) {
                error("Missing field in csv file header. The csv file must have fields named 'pool', 'meta', 'fastq_1', 'fastq_2' and 'lane'.")
            }
            else if (row.pool.contains(" ")) {
                error("Invalid value in csv file. Values for 'pool' and 'sample' can not contain space.")
            }
            [ [ row.pool.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both pool and sample are mandatory
        // Several sample can belong to the same pool
        // Sample should be unique for the pool
        if (row.pool) meta.pool = row.pool.toString()

        // mapping with fastq
        meta.id         = "${row.pool}-${row.lane}".toString()
        meta.anno = file(row.anno, checkIfExists: true)
        def fastq_1     = file(row.fastq_1, checkIfExists: true)
        def fastq_2     = file(row.fastq_2, checkIfExists: true)

        meta.num_lanes  = num_lanes.toInteger()
        meta.data_type  = 'fastq'

        meta.size       = 1 // default number of splitted fastq

        return [ meta, [ fastq_1, fastq_2 ] ]
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DINOFLOW {

    ch_versions = Channel.empty()

    input = extract_csv(ch_input)

    EXTRACT_BARCODE_TXT( input.map { meta, reads -> meta.anno})

    //EXTRACT_BARCODE_TXT.out.barcode.view()

    reads_barcode = input.combine(EXTRACT_BARCODE_TXT.out.barcode)
	.map { meta, reads, barcode -> [ meta + [whitelist: barcode, barcode_len: params.starsolo_barcode_len, umi_len: params.starsolo_umi_len, umi_start: params.starsolo_umi_start, cb_len: params.starsolo_cb_len, cb_start: params.starsolo_cb_start], params.starsolo_algorithm, reads ] }

    reads_barcode.view()
	ch_star_index.view()
    STARSOLO( reads_barcode, ch_star_index )
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // INPUT_CHECK (
    //     ch_input
    // )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    // workflow_summary    = WorkflowDinoflow.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowDinoflow.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
