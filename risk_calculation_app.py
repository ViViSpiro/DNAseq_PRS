import os
import tempfile
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import cyvcf2
from shiny import App, ui, render, reactive
import warnings
warnings.filterwarnings('ignore')

# Установка параметров для matplotlib
rcParams.update({
    'figure.autolayout': True,
    'font.family': 'DejaVu Sans',  # Шрифт поддерживающий кириллицу
    'axes.titlesize': 14,          # Размер заголовка
    'axes.labelsize': 12,          # Размер подписей осей
    'xtick.labelsize': 10,         # Размер подписей по оси X
    'ytick.labelsize': 10          # Размер подписей по оси Y
})

# Установка переменной окружения для Qt
os.environ['QT_QPA_PLATFORM'] = 'xcb'

def load_pgs_dict(file_path):
    """Загрузка словаря PGS с разделением на русские и английские названия"""
    pgs_dict = {}
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if ':' in line:
                    pgs_id, trait = line.split(':', 1)
                    # Разделяем русское и английское название
                    trait = trait.strip().strip('",')
                    if '(' in trait and ')' in trait:
                        rus_part = trait.split('(')[0].strip()
                        eng_part = trait.split('(')[1].split(')')[0].strip()
                        pgs_dict[pgs_id.strip('" ')] = {
                            'rus': rus_part,
                            'eng': eng_part,
                            'full': f"{rus_part} ({eng_part})"
                        }
                    else:
                        pgs_dict[pgs_id.strip('" ')] = {
                            'rus': trait,
                            'eng': '',
                            'full': trait
                        }
    except Exception as e:
        print(f"Ошибка загрузки словаря PGS: {str(e)}")
    return pgs_dict

def load_chr_names(file_path):
    """Загрузка соответствия названий хромосом"""
    chr_map = {}
    try:
        with open(file_path, 'r') as f:
            next(f)  # Пропускаем заголовок
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chr_num, genbank, refseq = parts[:3]
                    chr_map[genbank] = chr_num
                    chr_map[refseq] = chr_num
    except Exception as e:
        print(f"Ошибка загрузки соответствия хромосом: {str(e)}")
    return chr_map

# Конфигурация
CONFIG = {
    'pgs_dict': "dict.txt",
    'pgs_19_dir': "pgs_19/standardized_19",
    'pgs_38_dir': "pgs_38/standardized_38",
    'chr_names': "chr_names.txt"
}

# Загрузка данных
PGS_LIST = load_pgs_dict(CONFIG['pgs_dict'])
CHR_NAME_MAP = load_chr_names(CONFIG['chr_names'])

# Создаем список для выбора PGS с английскими названиями
pgs_choices = {k: f"{k} - {v['eng']}" if v['eng'] else k for k, v in PGS_LIST.items()}

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_select("genome_build", "Сборка генома", 
                       {"37": "GRCh37/hg19", "38": "GRCh38/hg38"}),
        ui.input_file("vcf_file", "Загрузить VCF файл",
                     accept=[".vcf", ".vcf.gz"],
                     multiple=False),
        ui.input_selectize("pgs_scores", "Выберите PGS шкалы",
                         pgs_choices,
                         multiple=True),
        ui.input_action_button("calculate", "Рассчитать риски"),
        width=350,
    ),
    ui.card(
        ui.output_plot("risk_plot", height="auto"),
        ui.output_text("status_text"),
        ui.output_table("results_table"),
        full_screen=True,
    ),
    title="Калькулятор геномных рисков"
)

def server(input, output, session):
    @reactive.Calc
    def get_pgs_data():
        """Загрузка стандартизированных PGS данных"""
        genome_build = input.genome_build()
        selected_pgs = input.pgs_scores()
        
        if not selected_pgs:
            return None
        
        # Извлекаем только ID из выбранных значений (удаляем английские названия)
        selected_ids = [pgs.split(' - ')[0] for pgs in selected_pgs]
        
        pgs_data = {}
        base_dir = CONFIG['pgs_19_dir'] if genome_build == "37" else CONFIG['pgs_38_dir']
        
        for pgs_id in selected_ids:
            filename = f"{pgs_id}_hmPOS_GRCh{genome_build}.txt"
            filepath = os.path.join(base_dir, filename)
            
            try:
                df = pd.read_csv(
                    filepath, 
                    sep='\t',
                    usecols=['rsID', 'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight'],
                    dtype={
                        'rsID': str,
                        'chr_name': str,
                        'chr_position': int,
                        'effect_allele': str,
                        'other_allele': str,
                        'effect_weight': float
                    }
                )
                
                # Очистка данных
                df = df.rename(columns={
                    'chr_name': 'chr',
                    'chr_position': 'pos',
                    'effect_weight': 'weight'
                })
                
                # Нормализация хромосом
                df['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False)
                df['chr'] = df['chr'].str.replace(r'\.0$', '', regex=True)
                
                # Добавление метаданных
                pgs_info = PGS_LIST.get(pgs_id, {"rus": "Неизвестный признак", "eng": ""})
                df['trait'] = pgs_info['rus']
                df['trait_full'] = pgs_info['full']
                df['pgs_id'] = pgs_id
                
                pgs_data[pgs_id] = df
                
            except Exception as e:
                print(f"Ошибка обработки {pgs_id}: {str(e)}")
                continue
        
        return pgs_data if pgs_data else None

    @reactive.Calc
    def process_vcf():
        """Обработка VCF файла"""
        if not input.vcf_file():
            return pd.DataFrame()
            
        file_info = input.vcf_file()[0]
        
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.vcf') as tmp:
                with open(file_info['datapath'], 'rb') as src:
                    tmp.write(src.read())
                tmp_path = tmp.name
            
            if not os.path.exists(tmp_path) or os.path.getsize(tmp_path) == 0:
                print("Ошибка: Неверный временный файл")
                return pd.DataFrame()
            
            vcf_reader = cyvcf2.VCF(tmp_path)
            variants = []
            
            for variant in vcf_reader:
                try:
                    original_chrom = str(variant.CHROM)
                    normalized_chrom = CHR_NAME_MAP.get(original_chrom, original_chrom)
                    normalized_chrom = normalized_chrom.replace('chr', '')
                    
                    # Получаем генотип (0/0, 0/1, 1/1)
                    gt = variant.genotypes[0][:2] if variant.genotypes else (0, 0)
                    
                    variants.append({
                        'chr': normalized_chrom,
                        'pos': int(variant.POS),
                        'ref': str(variant.REF),
                        'alt': str(variant.ALT[0]) if variant.ALT else str(variant.REF),
                        'gt': gt,
                        'dosage': sum(gt)  # Сумма аллелей (0, 1 или 2)
                    })
                except Exception as e:
                    print(f"Ошибка варианта: {str(e)}")
                    continue
            
            return pd.DataFrame(variants) if variants else pd.DataFrame()
            
        except Exception as e:
            print(f"Ошибка обработки VCF: {str(e)}")
            return pd.DataFrame()
            
        finally:
            if 'tmp_path' in locals() and os.path.exists(tmp_path):
                os.unlink(tmp_path)
            if 'vcf_reader' in locals():
                vcf_reader.close()
    
    @reactive.Calc
    def calculate_risks():
        """Расчет полигенных рисков"""
        pgs_data = get_pgs_data()
        vcf_data = process_vcf()
        
        if pgs_data is None or not isinstance(pgs_data, dict) or not pgs_data:
            return None
        
        results = []
        
        for pgs_id, pgs_df in pgs_data.items():
            # Создаем запись даже если нет VCF данных
            if vcf_data is None or not isinstance(vcf_data, pd.DataFrame) or vcf_data.empty:
                results.append({
                    'PGS ID': pgs_id,
                    'Признак': pgs_df['trait'].iloc[0],
                    'Оценка риска': 0,
                    'Совпадений': 0,
                    'Всего SNP': len(pgs_df)
                })
                continue
                
            # Объединение по позиции и аллелям
            merged = pd.merge(
                vcf_data,
                pgs_df,
                left_on=['chr', 'pos', 'ref', 'alt'],
                right_on=['chr', 'pos', 'other_allele', 'effect_allele'],
                how='inner'
            )
            
            # Расчет score (даже если нет совпадений)
            score = 0
            matched = 0
            
            if not merged.empty:
                # Расчет score
                merged['effect'] = merged['dosage'] * merged['weight']
                score = merged['effect'].sum()
                matched = len(merged)        
            
            
            results.append({
                'PGS ID': pgs_id,
                'Признак': pgs_df['trait'].iloc[0],
                'Оценка риска': score,
                'Совпадений': matched,
                'Всего SNP': len(pgs_df)                
            })
        
        return pd.DataFrame(results) if results else None

    @output
    @render.plot
    @reactive.event(input.calculate)
    def risk_plot():
        """Визуализация результатов с динамическими размерами шрифтов"""
        results_df = calculate_risks()
        
        if results_df is None or results_df.empty:
            fig, ax = plt.subplots(figsize=(10, 2))
            ax.text(0.5, 0.5, 'Нет данных для отображения', 
                ha='center', va='center', fontsize=12)
            ax.axis('off')
            return fig
        
        # Сортируем по оценке риска
        results_df = results_df.sort_values('Оценка риска', ascending=True)
        num_pgs = len(results_df)
        
        # Динамические параметры в зависимости от количества PGS
        if num_pgs <= 5:
            base_height = 6
            font_sizes = {
                'title': 14,
                'axis_label': 12,
                'tick_labels': 9,
                'bar_labels': 9,
                'y_labels': 10
            }
        elif num_pgs <= 10:
            base_height = 8
            font_sizes = {
                'title': 13,
                'axis_label': 11,
                'tick_labels': 8,
                'bar_labels': 8,
                'y_labels': 9
            }
        elif num_pgs <= 20:
            base_height = 10
            font_sizes = {
                'title': 11,
                'axis_label': 9,
                'tick_labels': 6,
                'bar_labels': 6,
                'y_labels': 7
            }
        else:
            base_height = 12
            font_sizes = {
                'title': 12,
                'axis_label': 10,
                'tick_labels': 7,
                'bar_labels': 7,
                'y_labels': 8
            }
        
        # Динамическая высота графика
        height_per_bar = 0.6
        fig_height = max(base_height, num_pgs * height_per_bar)
        fig_width = 12
        
        # Создаем фигуру
        fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
        ax = fig.add_subplot()
        
        # Цвета столбцов
        colors = ['#1f77b4' if x >= 0 else '#d62728' for x in results_df['Оценка риска']]
        
        # Рисуем горизонтальные столбцы
        y_pos = np.arange(num_pgs)
        bar_height = 0.5
        bars = ax.barh(y_pos, results_df['Оценка риска'], height=bar_height, color=colors)
        
        # Настройка осей и заголовков с динамическими размерами шрифтов
        ax.set_xlabel('Полигенная оценка риска', fontsize=font_sizes['axis_label'], labelpad=10)
        ax.set_title('Результаты анализа PGS', fontsize=font_sizes['title'], pad=20)
        ax.grid(axis='x', linestyle='--', alpha=0.6)
        
        # Подписи оси Y (русские названия) с динамическим размером
        ax.set_yticks(y_pos)
        ax.set_yticklabels(results_df['Признак'], fontsize=font_sizes['y_labels'])
        
        # Подписи значений на столбцах с динамическим размером
        for i, (bar, score) in enumerate(zip(bars, results_df['Оценка риска'])):
            width = bar.get_width()
            x_pos = width + 0.05*abs(width) if width >=0 else width - 0.05*abs(width)
            ha = 'left' if width >=0 else 'right'
            color = 'black'
            
            ax.text(x_pos, i, f"{score:.4f}", 
                    ha=ha, va='center', color=color, fontsize=font_sizes['bar_labels'])
        
        # Размер подписей на оси X
        ax.tick_params(axis='x', labelsize=font_sizes['tick_labels'])
        
        # Ручная настройка отступов
        left_margin = 0.3 + 0.01*max([len(x) for x in results_df['Признак']])
        fig.subplots_adjust(
            left=left_margin,
            right=0.9,
            top=0.9 - 0.01*num_pgs,
            bottom=0.15
        )
        
        # Автомасштабирование с запасом
        x_min = results_df['Оценка риска'].min()
        x_max = results_df['Оценка риска'].max()
        x_pad = max(abs(x_min), abs(x_max)) * 0.25
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        
        return fig

    @output
    @render.table
    def results_table():
        """Таблица результатов"""
        results_df = calculate_risks()
        if results_df is None or results_df.empty:
            return None
        
        # Форматируем числа для отображения
        display_df = results_df.copy()
        display_df['Оценка риска'] = display_df['Оценка риска'].apply(lambda x: f"{x:.4f}")
        return display_df

    @output
    @render.text
    def status_text():
        """Статус выполнения"""
        if not input.vcf_file():
            return "Загрузите VCF файл для анализа"
            
        vcf_status = process_vcf()
        if vcf_status is None or vcf_status.empty:
            return "Ошибка: Не удалось обработать VCF файл или файл пуст"
            
        pgs_status = get_pgs_data()
        if not pgs_status:
            return "Ошибка: Нет данных PGS для анализа"
            
        results = calculate_risks()
        if results is None:
            return "Ошибка при расчете рисков"
            
        total_snps = sum(len(df) for df in pgs_status.values())
        matched_snps = sum(results['Совпадений'])
        
        
        return (
            f"Анализ завершен\n"
            f"Обработано PGS: {len(pgs_status)}\n"
            f"Совпадений: {matched_snps} из {total_snps} SNP\n"            
        )

app = App(app_ui, server)