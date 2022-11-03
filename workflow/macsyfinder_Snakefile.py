#snakemake --snakefile macsyfinder_Snakefile.py --config output="/Users/rdenise/Documents/PostDoc/PLP/macsyfinder/pdx/" systems="pdx" def="/Users/rdenise/Documents/PostDoc/PLP/macsyfinder/pdx/DEF_model/" profiles="/Users/rdenise/Documents/PostDoc/PLP/macsyfinder/pdx/HMM_profile/" -j 8 --restart-time 10

import os

##########################################################################
##########################################################################
# GLOBAL VARIABLES

# HOME = "/ufrc/lagard/remi.denise/"
# GEM = "/ufrc/lagard/remi.denise/Gembase/"

HOME = "/Users/rdenise/Documents"
GEM = "/Volumes/Seagate_5TB/Gembase"


# DATABASE = "/pasteur/projets/policy01/gem/gembases/Microbial_B_1116/Genomes/Proteins"
# DATABASE = os.path.join(GEM, "gembases","Microbial_C_0218", "Genomes", "Proteins")
DATABASE = os.path.join(GEM,"Microbial_D_0419", "Genomes", "Proteins")

ALL_FILE, = glob_wildcards(os.path.join(DATABASE , "{replicon}.prt"))
# ALL_FILE = glob.glob(os.path.join(DATABASE , "*.prt"))
# ALL_FILE = [os.path.basename(i)[:-4] for i in ALL_FILE]

##########################################################################
# NEEDED BEFORE THE RULE OF INTERREST

# MACSYFINDER
FOLDER_MACSYFINDER_ANALYSIS = config["output"]

FOLDER_MACSYFINDER_REPLICON = os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "ALL_REPLICON")


##########################################################################
##########################################################################

rule all:
    input:
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.report"),
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.out"),
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.summary"),
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.tab"),
    params:
        error_db = os.path.join(DATABASE, "formatdb.err")
    shell :
        """
        if [[ -e {params.error_db} ]]; then
            rm {params.error_db}
        fi
        """

##########################################################################
##########################################################################

######
# ETAPE 1 Utilisation de macsyfinder
#####

PROFIL_FOLDER = config['profiles']

# MACSYFINDER = os.path.join(HOME, "script", "python2", "macsyfinder", "macsyfinder")
MACSYFINDER = '/Users/rdenise/Desktop/macsyfinder/macsyfinder/macsyfinder'

##########################################################################
# On lance macsyfinder

rule macsyfinder:
    input :
        seq = os.path.join(DATABASE, "{replicon}.prt")
    output :
        report = os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.report"),
        out = os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.out"),
        tab = os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.tab"),
        summary = os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.summary"),
    params :
        system_studied = config['systems'],
        folder = os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}"),
        def_folder = config['def'],
        macsyfinder = MACSYFINDER,
        profile = PROFIL_FOLDER,
        replicon = "{replicon}"
    threads :
        1
    shell :
        """
        if [[ -e {params.folder} ]]; then
            rm -rf {params.folder}
        fi

        
        if [[ -s {input.seq} ]] ; then
            python2.7 {params.macsyfinder} --sequence-db {input.seq} --db-type gembase --replicon-topology circular -d {params.def_folder} -p {params.profile} --profile-suffix .hmm -o {params.folder} -w {threads} {params.system_studied} > tmp_{params.replicon}
            touch tmp_{params.replicon}
            rm tmp_{params.replicon}
        else
            mkdir -p {params.folder}
        fi

        touch {output.report}
        touch {output.out}
        touch {output.summary}
        touch {output.tab}
        
        rm {input.seq}.[ip]*

        """
##########################################################################
##########################################################################   

########
# Etape 2 : Extraction info
########

##########################################################################

rule merge_out:
    input:
        all_out = sorted(expand(os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.out"), replicon=ALL_FILE))
    output :
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.out")     
    run :
        with open(str(output), "wt") as w_file :
            with open(str(input.all_out[0]), "rt") as r_file :
                for line in r_file.read() :
                    if line.startswith("*") :
                        break
                    else :
                        w_file.write(line)
            
            Replicon_nohit = "--- Replicons with no hits: ---\n\n"

            for file_out in input.all_out :
                nohit_begin = True
                hit_begin = False

                with open(file_out, "rt") as r_file :
                    for line in r_file :

                        if line.startswith("*") and hit_begin == False :
                            hit_begin = True
                            w_file.write(line)
                        elif "Replicons" in line :
                            hit_begin = False
                            nohit_begin = False
                        elif hit_begin :
                            w_file.write(line)

                    if nohit_begin :
                        Replicon_nohit += file_out.split("/")[-2]+"\n"

            w_file.write(Replicon_nohit)

##########################################################################
##########################################################################   

rule merge_report:
    input:
        all_out = sorted(expand(os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.report"), replicon=ALL_FILE))
    output :
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.report")    
    run :
        import pandas as pd

        header = "Hit_Id\tReplicon_name\tPosition\tSequence_length\tGene\tReference_system\tPredicted_system\tSystem_Id\tSystem_status\tGene_status\ti-evalue\tScore\tProfile_coverage\tSequence_coverage\tBegin_match\tEnd_match"

        with open(str(output), "wt") as w_file :
            w_file.write(header+"\n")

            for file_out in input.all_out : 
                pd.read_table(file_out, names=header.split(), comment="#", skiprows=1).to_csv(w_file, index=False, header=False, sep="\t")

##########################################################################
##########################################################################   

rule merge_summary:
    input:
        all_out = sorted(expand(os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.summary"), replicon=ALL_FILE))
    output :
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.summary")       
    run :
        import pandas as pd

        header = "Replicon_name\tSystem_Id\tReference_system\tSystem_status\tNb_loci\tNb_Ref_mandatory\tNb_Ref_accessory\tNb_Ref_Genes_detected_NR\tNb_Genes_with_match\tSystem_length\tNb_Mandatory_NR\tNb_Accessory_NR\tNb_missing_mandatory\tNb_missing_accessory\tList_missing_mandatory\tList_missing_accessory\tLoci_positions\tOccur_Mandatory\tOccur_Accessory\tOccur_Forbidden"

        with open(str(output), "wt") as w_file :
            w_file.write(header+"\n")

            for file_out in input.all_out : 
                pd.read_table(file_out, names=header.split(), comment="#", skiprows=1).to_csv(w_file, index=False, header=False, sep="\t")

##########################################################################
##########################################################################   

rule merge_tab:
    input:
        all_out = sorted(expand(os.path.join(FOLDER_MACSYFINDER_REPLICON, "{replicon}", "macsyfinder.tab"), replicon=ALL_FILE))
    output :
        os.path.join(FOLDER_MACSYFINDER_ANALYSIS, "macsyfinder.tab")        
    run :
        import pandas as pd

        for file_tab in input.all_out :
            if os.stat(file_tab).st_size != 0 :
                with open(file_tab, "rt") as r_file :
                    header = r_file.readline().rstrip().replace("#", "")
                    break

        with open(str(output), "wt") as w_file :
            w_file.write(header+"\n")

            for file_out in input.all_out : 
                pd.read_table(file_out, names=header.split(), comment="#", skiprows=1).to_csv(w_file, index=False, header=False, sep="\t")
