import os
import re
import subprocess
import pandas as pd
from glob import glob

def run_command(cmd):
    """Выполняет shell команду и проверяет результат"""
    print(f"Выполняю: {cmd}")
    result = subprocess.run(cmd, shell=True, check=True)
    if result.returncode != 0:
        raise RuntimeError(f"Ошибка выполнения команды: {cmd}")
    return result

def get_sample_names(input_dir):
    """Автоматически определяет имена образцов по FASTQ файлам"""
    fastq_files = glob(f"{input_dir}/*.fastq.gz")
    if not fastq_files:
        raise ValueError(f"Не найдены FASTQ файлы в директории {input_dir}")

    samples = set()
    for f in fastq_files:
        fname = os.path.basename(f).replace('.fastq.gz', '')
        sample = re.sub(r'_R[12](_\d+)?$', '', fname)
        samples.add(sample)
    
    return sorted(samples)

def process_bed_file(file_path):
    """Заменяет имена хромосом в BED файле"""
    chrom_mapping = {
        'CM000663.2': 'chr1',
        'CM000664.2': 'chr2',
        'CM000665.2': 'chr3',
        'CM000666.2': 'chr4',
        'CM000667.2': 'chr5',
        'CM000668.2': 'chr6',
        'CM000669.2': 'chr7',
        'CM000670.2': 'chr8',
        'CM000671.2': 'chr9',
        'CM000672.2': 'chr10',
        'CM000673.2': 'chr11',
        'CM000674.2': 'chr12',
        'CM000675.2': 'chr13',
        'CM000676.2': 'chr14',
        'CM000677.2': 'chr15',
        'CM000678.2': 'chr16',
        'CM000679.2': 'chr17',
        'CM000680.2': 'chr18',
        'CM000681.2': 'chr19',
        'CM000682.2': 'chr20',
        'CM000683.2': 'chr21',
        'CM000684.2': 'chr22',
        'CM000685.2': 'chrX',
        'CM000686.2': 'chrY',
        'J01415.2': 'chrM'
    }
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    for old_name, new_name in chrom_mapping.items():
        content = content.replace(old_name, new_name)
    
    with open(file_path, 'w') as f:
        f.write(content)

def quality_control(input_dir, sample_name):
    """Выполняет контроль качества"""
    print(f"\n=== Контроль качества для образца {sample_name} ===")
    os.makedirs("fqc", exist_ok=True)
    run_command(f"fastqc -o fqc {input_dir}/{sample_name}*.fastq.gz")
    
    os.makedirs("mqc", exist_ok=True)
    run_command("multiqc -o mqc fqc")

def trim_reads(input_dir, sample_name):
    """Подрезает чтения"""
    print(f"\n=== Подрезка чтений для образца {sample_name} ===")
    os.makedirs("trimmed", exist_ok=True)
    r1 = f"{input_dir}/{sample_name}_R1*.fastq.gz"
    r2 = f"{input_dir}/{sample_name}_R2*.fastq.gz"
    
    cmd = (
        f"fastp --detect_adapter_for_pe --overrepresentation_analysis --correction "
        f"--cut_right --thread 2 --html trimmed/{sample_name}.fastp.html "
        f"--json trimmed/{sample_name}.fastp.json -i {r1} -I {r2} "
        f"-o trimmed/{sample_name}_R1_trimmed.fastq.gz -O trimmed/{sample_name}_R2_trimmed.fastq.gz"
    )
    run_command(cmd)
    
    # QC после подрезки
    os.makedirs("fqc_trim", exist_ok=True)
    run_command(f"fastqc -o fqc_trim trimmed/{sample_name}_R*_trimmed.fastq.gz")
    os.makedirs("mqc_trim", exist_ok=True)
    run_command("multiqc -o mqc_trim fqc_trim")

def align_reads(sample_name, reference_dir):
    """Выравнивает чтения на референсный геном"""
    print(f"\n=== Выравнивание для образца {sample_name} ===")
    
    # Находим референсный файл
    ref_files = glob(f"{reference_dir}/*.fna.gz") or glob(f"{reference_dir}/*.fa.gz")
    if not ref_files:
        raise FileNotFoundError(f"Не найден референсный файл в {reference_dir}")
    
    reference_gz = ref_files[0]
    reference_path = os.path.join("reference", os.path.basename(reference_gz).replace('.gz', ''))
    
    # Создаем директорию reference если её нет
    os.makedirs("reference", exist_ok=True)
    
    # Распаковываем если нужно
    if not os.path.exists(reference_path):
        print(f"Распаковываю референсный геном в {reference_path}...")
        run_command(f"gunzip -k -c {reference_gz} > {reference_path}")
    
    # Индексирование референса
    if not os.path.exists(f"{reference_path}.bwt"):
        run_command(f"bwa index {reference_path}")
    
    if not os.path.exists(f"{reference_path}.fai"):
        run_command(f"samtools faidx {reference_path}")
    
    if not os.path.exists(reference_path.replace('.fna', '.dict')):
        run_command(f"gatk CreateSequenceDictionary -R {reference_path}")
    
    # Выравнивание
    os.makedirs("bam", exist_ok=True)
    r1 = f"trimmed/{sample_name}_R1_trimmed.fastq.gz"
    r2 = f"trimmed/{sample_name}_R2_trimmed.fastq.gz"
    
    cmd = (
        f"bwa mem -t 4 {reference_path} {r1} {r2} | "
        f"samtools view -b | samtools sort -o bam/{sample_name}.bam"
    )
    run_command(cmd)
    
    # Добавление read groups
    os.makedirs("bam_rg", exist_ok=True)
    run_command(
        f"gatk AddOrReplaceReadGroups -I bam/{sample_name}.bam "
        f"-O bam_rg/{sample_name}_rg.bam --RGID {sample_name} --RGLB lib1 "
        f"--RGPL ILLUMINA --RGPU unit1 --RGSM {sample_name}"
    )
    
    # Индексирование BAM
    run_command(f"samtools index bam_rg/{sample_name}_rg.bam")

def variant_calling(sample_name, reference_dir):
    """Вызов вариантов"""
    print(f"\n=== Вызов вариантов для образца {sample_name} ===")
    os.makedirs("res_vcf", exist_ok=True)
    
    # Находим распакованный референс
    ref_files = glob("reference/*.fna") or glob("reference/*.fa")
    if not ref_files:
        raise FileNotFoundError("Не найден распакованный референсный файл")
    
    reference_path = ref_files[0]
    run_command(
        f"gatk HaplotypeCaller --java-options '-Xmx4g' -R {reference_path} "
        f"-I bam_rg/{sample_name}_rg.bam -O res_vcf/{sample_name}_rg.vcf"
    )

def coverage_analysis(sample_name):
    """Анализ покрытия"""
    print(f"\n=== Анализ покрытия для образца {sample_name} ===")
    os.makedirs("bed", exist_ok=True)
    run_command(
        f"bedtools genomecov -bg -ibam bam_rg/{sample_name}_rg.bam > bed/{sample_name}_rg.bed"
    )
    
    # Обработка BED файла
    bed_file = f"bed/{sample_name}_rg.bed"
    process_bed_file(bed_file)
    
    # Фильтрация только chr
    os.makedirs("bed_clean", exist_ok=True)
    with open(bed_file, 'r') as infile, open(f"bed_clean/{sample_name}_rg.bed", 'w') as outfile:
        for line in infile:
            if line.startswith('chr'):
                outfile.write(line)

def intersect_with_target(sample_name, target_bed):
    """Пересечение с целевыми регионами"""
    print(f"\n=== Пересечение с целевыми регионами для образца {sample_name} ===")
    os.makedirs("intersect", exist_ok=True)
    os.makedirs("target", exist_ok=True)
    
    # Копируем целевые регионы (если нужно)
    target_name = os.path.basename(target_bed)
    if not os.path.exists(f"target/{target_name}"):
        run_command(f"cp {target_bed} target/")
    
    run_command(
        f"bedtools intersect -a bed_clean/{sample_name}_rg.bed "
        f"-b target/{target_name} > intersect/{sample_name}_rg.bed"
    )

def calculate_coverage_stats(sample_name):
    """Вычисляет статистики покрытия"""
    print(f"\n=== Расчет статистик покрытия для образца {sample_name} ===")
    os.makedirs("mean_less10_all", exist_ok=True)
    
    bed_file = f"intersect/{sample_name}_rg.bed"
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"Файл {bed_file} не найден")
    
    file = pd.read_csv(bed_file, sep='\t', header=None, quoting=3)
    
    mean_all = file[3].mean()
    less10_all = (file[3] < 10).mean()
    total_sum_reads = file[3].sum()
    
    results = pd.DataFrame([{
        'sample': sample_name,
        'mean_all': mean_all,
        'less10_all': less10_all,
        'total_sum_reads': total_sum_reads
    }])
    
    # Сохраняем индивидуальные результаты
    results.to_csv(f"mean_less10_all/{sample_name}_stats.csv", index=False)
    results.to_csv(f"mean_less10_all/{sample_name}_stats.txt", sep='\t', index=False)
    
    # Объединяем все результаты в один файл
    all_results = []
    for stats_file in glob("mean_less10_all/*_stats.csv"):
        all_results.append(pd.read_csv(stats_file))
    
    if all_results:
        pd.concat(all_results).to_csv("mean_less10_all/all_samples_stats.csv", index=False)
        pd.concat(all_results).to_csv("mean_less10_all/all_samples_stats.txt", sep='\t', index=False)

def process_sample(sample_name, input_dir, reference_dir, target_bed):
    """Обрабатывает один образец"""
    try:
        quality_control(input_dir, sample_name)
        trim_reads(input_dir, sample_name)
        align_reads(sample_name, reference_dir)
        variant_calling(sample_name, reference_dir)
        coverage_analysis(sample_name)
        intersect_with_target(sample_name, target_bed)
        calculate_coverage_stats(sample_name)
        
        print(f"\n=== Обработка образца {sample_name} завершена успешно ===")
    except Exception as e:
        print(f"\n!!! Ошибка при обработке образца {sample_name}: {e}")

def main():
    # Параметры
    input_dir = "data"
    reference_dir = "reference"
    target_bed = "target/LM_v4__roi.bed"
    
    try:
        # Получаем список образцов автоматически
        sample_names = get_sample_names(input_dir)
        if not sample_names:
            raise ValueError("Не найдены образцы для обработки")
        
        print(f"Найдены образцы: {', '.join(sample_names)}")
        
        # Обрабатываем каждый образец
        for sample_name in sample_names:
            process_sample(sample_name, input_dir, reference_dir, target_bed)
        
        print("\n=== Весь анализ успешно завершен ===")
    except Exception as e:
        print(f"\n!!! Критическая ошибка: {e}")
        print("Анализ прерван")

if __name__ == "__main__":
    main()