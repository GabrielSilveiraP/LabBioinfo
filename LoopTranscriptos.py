import subprocess
import glob
import re
import os

# 1. Definindo os diretorios
# 'working_dir' is where your FASTQ files are and where SAM will be generated.
working_dir = "/resultados/Analises_Gabriel/Testes1/fiocruz/Fatbody"
bwa_exe = "/resultados/bin_software/bwa-mem2-v2.2.1/bwa-mem2"
# genome de referencia
ref_genome = "/resultados/Analises_Gabriel/Testes1/fiocruz/ref/TrangeliSC58"
#--- Parte 1: Achar os arquivos
print("Scanning for files...")
#Procurando os arqvs r1 e identificando as amostras
search_pattern = os.path.join(working_dir, "*_Unmapped.out_R1.trim.fastq.gz")
# os. -> aqui ele ta só fazendo o caminho do glob
files = glob.glob(search_pattern)
# glob. é pra organizar e aqui ta sendo usado para procurar semelhanças dentro dos arqvs com aql final
samples = sorted({
    int(os.path.basename(f).split("_")[0])
    for f in glob.glob(os.path.join(working_dir, "*_Unmapped.out_R1.trim.fastq.gz"))
})
#--- Parte 2: parte do processamento
for sample in samples:
    print(f"\n>>> Processing Sample {sample} <<<")

    # Define Filenames (Absolute paths for Host execution)
    r1_file = os.path.join(working_dir, f"{sample}_Unmapped.out_R1.trim.fastq.gz")
    r2_file = os.path.join(working_dir, f"{sample}_Unmapped.out_R2.trim.fastq.gz")

    # We use simple names for the Docker mounting logic, but full paths for BWA
    sam_filename = f"transcriptosMapeados.{sample}.sam"
    bam_filename = f"transcriptosMapeados.{sample}.bam"
    bam_filtrado = f"transcriptoFiltrado.{sample}.bam"

    sam_full_path = os.path.join(working_dir, sam_filename)
    # Essa variavel n é muito necessaria, mas já que tamo indeando tudo deixa essa aqui tbm
    # ---------------------------------------------------------
    # STEP 1: BWA
    # ---------------------------------------------------------
    # Lembrar que esse bwa ta em python, ent algumas coisas vão ficar diferentes
    cmd_bwa = [
        bwa_exe, "mem", "-t", "15", ref_genome, r1_file, r2_file, "-a", "-o", sam_full_path
    ]

    print(f"1. Running BWA for sample {sample}...")
    subprocess.run(cmd_bwa, check=True)

    # ---------------------------------------------------------
    # STEP 2: Samtools pt.1
    # ---------------------------------------------------------
    samtoolsSamtoBam_cmd = (
        f"docker run -i --rm "
        f"-v {working_dir}:/data "
        f"-w /data "
        f"staphb/samtools:1.21 "
        f"samtools view -b -@ 15 {sam_filename} -o {bam_filename} "
    )
    # ---------------------------------------------------------
    # STEP 2: Samtools pt.2
    # ---------------------------------------------------------
    # "working_dir é o diretorio montado dentro do docker

    samtools_cmd = (
        f"docker run -i --rm "
        f"-v {working_dir}:/data "
        f"-w /data "
        f"staphb/samtools:1.21 "
        f"samtools view -f 0x2 -@ 15 {bam_filename} -o {bam_filtrado} "
    )
    print(f"2. Running Samtools (Docker) for sample {sample}...")
    subprocess.run(samtoolsSamtoBam_cmd, shell=True, check=True)
    subprocess.run(samtools_cmd, shell=True, check=True)

    # ---------------------------------------------------------
    # STEP 3: Removedor de .sam
    # ---------------------------------------------------------
    print(f"3. Removing intermediate SAM file...")
    subprocess.run(["rm", sam_full_path], check=True)

print("\nAll processing finished.")