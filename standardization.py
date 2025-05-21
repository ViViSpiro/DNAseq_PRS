import os
import pandas as pd
import gzip

def split_alleles(row):
    """Функция для разделения строк с несколькими аллелями"""
    if pd.isna(row['effect_allele']):
        return [row]
    
    alleles = str(row['effect_allele']).split("/")
    rows = []
    for allele in alleles:
        if allele.strip():
            new_row = row.copy()
            new_row['effect_allele'] = allele.strip()
            rows.append(new_row)
    return rows

# Пути к папкам
input_dirs = ["pgs_19", "pgs_38"]
output_dirs = ["pgs_19/standardized_19", "pgs_38/standardized_38"]

# Создаем выходные директории
for output_dir in output_dirs:
    os.makedirs(output_dir, exist_ok=True)

# Обрабатываем файлы в каждой папке
for input_dir, output_dir in zip(input_dirs, output_dirs):
    file_names = [f for f in os.listdir(input_dir) if f.endswith('.txt.gz')]
    
    for file_name in file_names:
        input_path = os.path.join(input_dir, file_name)
        output_path = os.path.join(output_dir, file_name.replace('.gz', ''))
        
        try:
            # Определяем заголовки
            with gzip.open(input_path, 'rt') as f:
                header_line = None
                for line in f:
                    if not line.startswith('###'):
                        header_line = line.strip()
                        break
            
            if not header_line:
                print(f"Файл {file_name} не содержит заголовков. Пропускаем.")
                continue
            
            # Читаем файл по чанкам
            chunk_size = 10000
            chunks = pd.read_csv(input_path, 
                               compression='gzip',
                               sep='\t',
                               comment='#',
                               on_bad_lines='warn',
                               chunksize=chunk_size,
                               low_memory=False)
            
            processed_chunks = []
            for chunk in chunks:
                if chunk.empty:
                    continue
                
                # ОРИГИНАЛЬНАЯ ОБРАБОТКА СТОЛБЦОВ (без изменений)
                if 'chr_position' not in chunk.columns and 'hm_pos' in chunk.columns:
                    chunk.rename(columns={'hm_pos': 'chr_position'}, inplace=True)
                
                if 'rsID' not in chunk.columns:
                    if 'chr_position' in chunk.columns:
                        chunk['rsID'] = chunk['chr_position']
                    else:
                        continue    
                
                if 'other_allele' not in chunk.columns:
                    if 'effect_allele' in chunk.columns:
                        chunk['other_allele'] = chunk['effect_allele']

                        if 'hm_inferOtherAllele' in chunk.columns:
                            if not chunk['hm_inferOtherAllele'].isnull().all():
                                chunk['effect_allele'] = chunk['hm_inferOtherAllele']
                            elif 'variant_description' in chunk.columns:
                                def extract_alleles(desc):
                                    if pd.isna(desc):
                                        return None, None
                                    parts = desc.split(":")
                                    if len(parts) >= 4:
                                        return parts[-1], parts[-2]
                                    return None, None
                                
                                chunk[['effect_allele_new', 'other_allele_new']] = chunk['variant_description'].apply(extract_alleles).apply(pd.Series)
                                chunk['effect_allele'] = chunk['effect_allele_new']
                                chunk['other_allele'] = chunk['other_allele_new']
                                chunk = chunk.drop(['effect_allele_new', 'other_allele_new'], axis=1)
                        elif 'variant_description' in chunk.columns:
                            def extract_alleles(desc):
                                if pd.isna(desc):
                                    return None, None
                                parts = desc.split(":")
                                if len(parts) >= 4:
                                    return parts[-1], parts[-2]
                                return None, None
                            
                            chunk[['effect_allele_new', 'other_allele_new']] = chunk['variant_description'].apply(extract_alleles).apply(pd.Series)
                            chunk['effect_allele'] = chunk['effect_allele_new']
                            chunk['other_allele'] = chunk['other_allele_new']
                            chunk = chunk.drop(['effect_allele_new', 'other_allele_new'], axis=1)

                if 'chr_name' not in chunk.columns and 'hm_chr' in chunk.columns:
                    chunk['chr_name'] = chunk['hm_chr']
                
                # Разделение аллелей
                split_rows = []
                for _, row in chunk.iterrows():
                    if 'effect_allele' in chunk.columns and "/" in str(row.get('effect_allele', '')):
                        split_rows.extend(split_alleles(row))
                    else:
                        split_rows.append(row)
                
                if split_rows:
                    processed_chunk = pd.DataFrame(split_rows)
                    
                    # Сохраняем только нужные столбцы
                    required_columns = ['rsID', 'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight']
                    available_columns = [col for col in required_columns if col in processed_chunk.columns]
                    if available_columns:
                        processed_chunks.append(processed_chunk[available_columns])
            
            # Сохраняем результат
            if processed_chunks:
                final_data = pd.concat(processed_chunks, ignore_index=True)
                final_data.to_csv(output_path, sep='\t', index=False)
                print(f"Успешно обработан файл: {input_path} -> {output_path}")
            else:
                print(f"Файл {file_name} не содержит данных для сохранения.")
                
        except Exception as e:
            print(f"Ошибка при обработке файла {file_name}: {str(e)}")
            continue