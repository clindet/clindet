rule fastp_trim:
    input:
        unpack(get_rna_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    conda:
        config['softwares']['fastp']['conda']
    shell:
        """{config[softwares][fastp][call]} --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 16 -Q -c -L \
            -h {output.html} -j {output.json}\
            -o {output.R1} -O {output.R2}
        """


# rule fastp_tumor_sample:
#     input:
#         unpack(get_tumor_fastq)
#     output:
#         R1="{project}/{genome_version}/results/trimmed/{sample}-T_R1.fastq.gz",
#         R2="{project}/{genome_version}/results/trimmed/{sample}-T_R2.fastq.gz"
#     params:
#         adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
#         extra="--merge"
#     threads: 2
#     conda:
#         config['softwares']['fastp']['conda']
#     shell:
#         """{config[softwares][fastp][call]} --thread {threads} \
#             -i {input.R1} -I {input.R2} \
#             -w 8 -Q -c -L \
#             -o {output.R1} -O {output.R2}
#         """

