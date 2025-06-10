import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import bootstrap, mannwhitneyu
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings("ignore")

# Настройки стиля
sns.set_theme(style="whitegrid", palette="husl")
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 12,
    'figure.figsize': (12, 8),
    'figure.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans'
})

# Создаем папки для результатов
os.makedirs("results/plots", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

def load_patient_data():
    """Загрузка данных пациентов"""
    try:
        df1 = pd.read_csv("ngs_pgs_intersections/pgs_scores_results.txt", sep="\t")
        df2 = pd.read_csv("ngs_pgs_intersections_cd/pgs_scores_results_cd_63.txt", sep="\t")
        patients = pd.concat([df1, df2])
        patients['Group'] = 'Patient'
        return patients
    except Exception as e:
        print(f"Ошибка загрузки данных пациентов: {str(e)}")
        raise

def load_healthy_data(pgs_list):
    """Загрузка данных здоровых доноров"""
    healthy_data = []
    
    for pgs in tqdm(pgs_list, desc="Loading healthy data"):
        paths = [
            f"selected_data/raw/{pgs}_hmPOS_GRCh37_res/plink.profile",
            f"selected_data_cd/raw_cd/{pgs}_hmPOS_GRCh37_res/plink.profile"
        ]
        
        for path in paths:
            if os.path.exists(path):
                try:
                    df = pd.read_csv(path, sep='\s+')
                    df['PGS'] = pgs
                    healthy_data.append(df[['PGS', 'SCORE']].rename(columns={'SCORE': 'Score'}))
                    break
                except Exception as e:
                    print(f"Error loading {path}: {str(e)}")
    
    if not healthy_data:
        raise ValueError("No healthy data found for any PGS")
    
    return pd.concat(healthy_data)

def calculate_bootstrap_ci(data, percentiles=[25, 50, 75], n_resamples=1000):
    """Расчет доверительных интервалов методом бутстрепа"""
    def percentile_func(x):
        return np.percentile(x, percentiles)
    
    res = bootstrap((data,), percentile_func, n_resamples=n_resamples)
    point_estimate = np.percentile(data, percentiles)
    
    if hasattr(res, 'confidence_interval'):
        ci_low = res.confidence_interval.low
        ci_high = res.confidence_interval.high
    else:
        ci_low = np.percentile(res.bootstrap_distribution, 2.5, axis=0)
        ci_high = np.percentile(res.bootstrap_distribution, 97.5, axis=0)
    
    return point_estimate, ci_low, ci_high

def create_final_dataset():
    """Создание итогового датасета с расчетом перцентилей"""
    print("Loading patient data...")
    patients = load_patient_data()
    unique_pgs = patients['PGS'].unique()
    
    print("Loading healthy data...")
    healthy = load_healthy_data(unique_pgs)
    
    print("Calculating statistics...")
    results = []
    
    for pgs in tqdm(unique_pgs, desc="Processing PGS"):
        healthy_sub = healthy[healthy['PGS'] == pgs]['Score'].values
        patient_sub = patients[patients['PGS'] == pgs]['Score'].values
        
        if len(healthy_sub) < 10 or len(patient_sub) < 1:
            continue
        
        # Вычисляем перцентили 25%, 50%, 75%
        point_est, ci_low, ci_high = calculate_bootstrap_ci(healthy_sub)
        
        try:
            stat, pval = mannwhitneyu(healthy_sub, patient_sub, alternative='two-sided')
        except Exception as e:
            print(f"Error in MWU test for {pgs}: {str(e)}")
            stat, pval = np.nan, np.nan
        
        results.append({
            'PGS': pgs,
            'p25': point_est[0], 'p25_ci_lower': ci_low[0], 'p25_ci_upper': ci_high[0],
            'p50': point_est[1], 'p50_ci_lower': ci_low[1], 'p50_ci_upper': ci_high[1],
            'p75': point_est[2], 'p75_ci_lower': ci_low[2], 'p75_ci_upper': ci_high[2],
            'mann_whitney_u': stat,
            'p_value': pval,
            'n_patients': len(patient_sub),
            'n_healthy': len(healthy_sub)
        })
    
    stats_df = pd.DataFrame(results)
    final_df = pd.merge(patients, stats_df, on='PGS', how='left')
    
    # Классификация риска по перцентилям
    conditions = [
        (final_df['Score'] < final_df['p25']),
        (final_df['Score'] >= final_df['p25']) & (final_df['Score'] <= final_df['p75']),
        (final_df['Score'] > final_df['p75'])
    ]
    choices = ['Низкий', 'Средний', 'Высокий']
    final_df['Risk'] = np.select(conditions, choices, default='Не определен')
    
    # Сохраняем результаты
    final_df.to_csv('results/tables/final_dataset_with_stats.txt', sep='\t', index=False, header=True)
    stats_df.to_csv('results/tables/statistical_summary.txt', sep='\t', index=False, header=True)
    
    return final_df, healthy, stats_df

def create_percentile_scatter_plot(final_df, healthy_df, stats_df):
    """Визуализация с перцентилями здоровых по X и значениями PGS пациентов по Y"""
    print("\nСоздание перцентильного scatter plot...")
    
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'figure.figsize': (12, 8),
        'figure.autolayout': True  # Автоподгонка layout
    })
    
    for _, row in tqdm(stats_df.iterrows(), total=len(stats_df), desc="Обработка PGS"):
        pgs = row['PGS']
        healthy = healthy_df[healthy_df['PGS'] == pgs]['Score'].dropna()
        patients = final_df[final_df['PGS'] == pgs]
        
        if len(healthy) < 10 or len(patients) < 1:
            continue
            
        # Создаем временный DataFrame для визуализации
        plot_df = patients.copy()
        
        # Рассчитываем перцентили для каждого пациента относительно здоровых
        plot_df['Pop_Percentile'] = plot_df['Score'].apply(
            lambda x: np.mean(healthy <= x) * 100
        )
        
        # Добавляем небольшой случайный разброс по X (не более ±2%)
        plot_df['Pop_Percentile_jittered'] = plot_df['Pop_Percentile'] + \
            np.random.uniform(-2, 2, size=len(plot_df))
        
        # Создаем категории риска
        plot_df['Risk_Category'] = pd.cut(
            plot_df['Pop_Percentile'],
            bins=[0, 25, 75, 100],
            labels=['Низкий', 'Средний', 'Высокий'],
            include_lowest=True
        )
        
        # Получаем фактические значения пациентов
        patient_scores = plot_df['Score'].values
        patient_min = np.min(patient_scores)
        patient_max = np.max(patient_scores)
        
        # Рассчитываем границы с запасом
        x_min = -5  # Запас слева по X
        x_max = 105  # Запас справа по X
        
        # Для оси Y используем 10% запас сверху и снизу
        y_padding = (patient_max - patient_min) * 0.1 if patient_max != patient_min else 1
        y_min = patient_min - y_padding
        y_max = patient_max + y_padding
        
        # Создаем фигуру с увеличенными границами
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 1. Scatter plot с реальными значениями PGS пациентов
        for risk, color in zip(['Низкий', 'Средний', 'Высокий'], ['green', 'orange', 'red']):
            subset = plot_df[plot_df['Risk_Category'] == risk]
            if len(subset) > 0:
                ax.scatter(
                    x=subset['Pop_Percentile_jittered'],
                    y=subset['Score'],
                    color=color,
                    s=100,
                    alpha=0.7,
                    label=f"{risk} (n={len(subset)})"
                )
        
        # 2. Устанавливаем жесткие границы с запасом
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        
        # 3. Добавляем линии перцентилей
        ax.axvline(25, color='blue', linestyle='--', alpha=0.7, label='25-й перцентиль')
        ax.axvline(50, color='black', linestyle='-', alpha=0.7, label='Медиана (50-й перцентиль)')
        ax.axvline(75, color='blue', linestyle='--', alpha=0.7, label='75-й перцентиль')
        
        # 4. Добавляем горизонтальную линию медианы здоровых
        healthy_median = np.median(healthy)
        ax.axhline(y=healthy_median, color='gray', linestyle=':', alpha=0.5, 
                  label=f'Медиана здоровых: {healthy_median:.6f}')
        
        # 5. Настройки графика
        title = f"PGS: {pgs}\n"
        title += f"Значения пациентов: {patient_min:.7f} - {patient_max:.7f}\n"
        title += f"(p-value: {row['p_value']:.3e})"
        
        ax.set_title(title)
        ax.set_xlabel("Перцентиль в здоровой популяции (%)")
        ax.set_ylabel("Значение PGS пациента")
        
        # Форматирование оси Y
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.5f'))
        
        # Настройка расположения элементов
        ax.legend(title='Категории риска', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Сохраняем с увеличенными границами
        plt.savefig(
            f"results/plots/{pgs}_percentile_scatter.png", 
            dpi=300, 
            bbox_inches='tight',
            pad_inches=0.5  # Дополнительный запас при сохранении
        )
        plt.close()
        
        print(f"Создан график: results/plots/{pgs}_percentile_scatter.png")

if __name__ == "__main__":
    print("=== PGS Risk Analysis Pipeline ===")
    try:
        final_df, healthy_df, stats_df = create_final_dataset()
        create_percentile_scatter_plot(final_df, healthy_df, stats_df)
        
        print("\n=== Analysis Complete ===")
        print(f"Results saved to:")
        print(f"- Final dataset: results/tables/final_dataset_with_stats.txt")
        print(f"- Statistical summary: results/tables/statistical_summary.txt")
        
        print("\nTop significant PGS:")
        print(stats_df.sort_values('p_value').head(10)[['PGS', 'p_value', 'n_patients']])
        
    except Exception as e:
        print(f"\n!!! Critical Error: {str(e)}")