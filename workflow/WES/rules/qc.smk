rule fastp_normal_sample:
    input:
        unpack(get_normal_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}-NC_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}-NC_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-NC-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-NC-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    conda:
        config['softwares']['fastp']['conda']
    shell:
        """ fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 8 -Q -c -L \
            -h {output.html} -j {output.json}\
            -o {output.R1} -O {output.R2}
        """

rule fastp_tumor_sample:
    input:
        unpack(get_tumor_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}-T_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}-T_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-T-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-T-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    conda:
        config['softwares']['fastp']['conda']
    shell:
        """fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 8 -Q -c -L \
            -h {output.html} -j {output.json} \
            -o {output.R1} -O {output.R2}
        """

