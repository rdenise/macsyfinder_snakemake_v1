##########################################################################


rule macsyfinder:
    input:
        install=os.path.join(
                OUTPUT_FOLDER,
                "log",
                "installation.done"
            ),
        fasta=os.path.join(FASTA_FOLDER, f"{{replicon}}.{EXT_FILE}"),
    output:
        results=directory(
            os.path.join(
                OUTPUT_FOLDER,
                "REPLICONS",
                "{replicon}",
            ),
        ),
        report=os.path.join(
            OUTPUT_FOLDER,
            "REPLICONS",
            "{replicon}",
        ),
        out=os.path.join(
            OUTPUT_FOLDER,
            "REPLICONS",
            "{replicon}",
            "macsyfinder.out",
        ),
        summary=os.path.join(
            OUTPUT_FOLDER,
            "REPLICONS",
            "{replicon}",
        ),
        tab=os.path.join(
            OUTPUT_FOLDER,
            "REPLICONS",
            "{replicon}",
        ),        
    params:
        macsydata_def=config["macsydata_def"],
        macsydata_profile=config["macsydata_profile"],
        db_type=config["db_type"],
        replicon_topology=config["replicon_topology"],
        models=lambda wildcards: dict_models[wildcards.replicon] if wildcards.replicon in dict_models else config["models"],
        macsyfinder=config["macsyfinder"],
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "macsyfinder",
            "{replicon}.log",
        ),
    threads: 1
    conda:
        "../envs/macsyfinder.yaml"
    shell:
        """
        if [[ -e {params.folder} ]]; then
            rm -rf {params.folder}
        fi

        if [[ -s {input.fasta:q} ]] ; then
            python2.7 {params.macsyfinder} -d {params.def_folder} -p {params.profile} --profile-suffix .hmm \
            --sequence-db {input.fasta:q} -o {output.results:q} --db-type {params.db_type} --replicon-topology {params.replicon_topology} \
            -w {threads} --models {params.models} &> {log:q}

            touch tmp_{params.replicon}
            rm tmp_{params.replicon}
        else
            echo {input.fasta}
            mkdir -p {output.results:q}
        fi

        touch {output.report}
        touch {output.out}
        touch {output.summary}
        touch {output.tab}
        
        rm {input.seq}.[ip]*       
        """


##########################################################################

rule macsymerge_out:
    input:
        all_out=sorted(
            expand(
                os.path.join(
                    OUTPUT_FOLDER, 
                    "REPLICONS",
                    "{replicon}", 
                    "macsyfinder.out"
                ), 
                replicon=REPLICON_NAME
            )
        )
    output:
        out=os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.out"
        ) 
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "macsymerge",
            "macsymerge.out.log",
        ),       
    threads: 1
    conda: "../envs/pandas.yaml"
    script: "../scripts/merge_out.py"


##########################################################################


rule macsymerge_report:
    input:
        all_out=sorted(
            expand(
                os.path.join(
                    OUTPUT_FOLDER, 
                    "REPLICONS",
                    "{replicon}", 
                    "macsyfinder.out"
                ), 
                replicon=REPLICON_NAME
            )
        )
    output:
        out=os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.report"
        )  
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "macsymerge",
            "macsymerge.report.log",
        ), 
    threads: 1
    conda: "../envs/pandas.yaml"          
    script: "../scripts/merge_report.py"


##########################################################################


rule merge_summary:
    input:
        all_out=sorted(
            expand(
                os.path.join(
                    OUTPUT_FOLDER, 
                    "REPLICONS",
                    "{replicon}", 
                    "macsyfinder.summary"
                ), 
                replicon=REPLICON_NAME
            )
        )    
    output:
        out=os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.summary"
        )  
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "macsymerge",
            "macsymerge.summary.log",
        ),   
    threads: 1
    conda: "../envs/pandas.yaml"          
    script: "../scripts/merge_summary.py"


##########################################################################


rule merge_tab:
    input:
        all_out=sorted(
            expand(
                os.path.join(
                    OUTPUT_FOLDER, 
                    "REPLICONS",
                    "{replicon}", 
                    "macsyfinder.tab"
                ), 
                replicon=REPLICON_NAME
            )
        )      
    output:
        out=os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.tab"
        )  
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "macsymerge",
            "macsymerge.tab.log",
        ),   
    threads: 1
    conda: "../envs/pandas.yaml"          
    script: "../scripts/merge_tab.py"


##########################################################################

