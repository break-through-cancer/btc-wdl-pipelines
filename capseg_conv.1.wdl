### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Shahab Sarmashghi, ssarmash@broadinstitute.org
### Date last updated: July 30, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute

version 1.0
#this is a workflow converts iconicc segfile output to capseg format
workflow make_capseg {

    input {
        File capseg_rscript
        File processed_counts
        File segfile
        String participant_id
    }

    #INSERT CALLS HERE
    call make_capseg {
        input:
            capseg_rscript = capseg_rscript,
            processed_counts = processed_counts,
            segfile = segfile,
            participant_id = participant_id
    }

    output{
        File capseg_file = make_capseg.seg_file
    }
}


task make_capseg {
    input{
        File capseg_rscript
        File processed_counts
        File segfile
        String participant_id

        String memory = "10 GB"
        Int timeMinutes = 1 + ceil(size(processed_counts, "G"))
        String r_dockerImage = "wchukwu/r-docker:latest"
    }

    command <<<
        Rscript ~{capseg_rscript} --segfile ~{segfile} --processed_cts ~{processed_counts} --participant_id ~{participant_id}
    >>>

    output {
        File seg_file = "~{participant_id}.capseg.txt"
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:r_dockerImage
    }
}