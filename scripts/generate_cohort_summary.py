#!/usr/bin/env python3
"""
Gera relatório sumário de todo o cohort
Overview de todas as amostras analisadas
"""

import sqlite3
import pandas as pd
from pathlib import Path
from datetime import datetime
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Parâmetros do Snakemake
db_path = snakemake.input.db
output_html = snakemake.output.cohort_html
output_pdf = snakemake.output.cohort_pdf if hasattr(snakemake.output, 'cohort_pdf') else None

def query_db(db_path, query, params=()):
    """Query no database"""
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn, params=params)
    conn.close()
    return df

def plot_to_html(fig):
    """Converte plot Plotly para HTML"""
    return fig.to_html(full_html=False, include_plotlyjs='cdn')

# Carrega dados do cohort
samples_df = query_db(db_path, "SELECT * FROM v_sample_summary")
n_samples = len(samples_df)

# Estatísticas gerais
purecn_df = query_db(db_path, "SELECT * FROM purecn_results")
mrd_df = query_db(db_path, "SELECT * FROM mrd_results")
classification_df = query_db(db_path, "SELECT * FROM leukemia_classification")
drivers_df = query_db(db_path, "SELECT * FROM driver_mutations")
fusions_df = query_db(db_path, "SELECT * FROM gene_fusions")

# ========== GRÁFICOS ==========

# 1. Distribuição de pureza
fig_purity = px.histogram(
    purecn_df, x='purity',
    nbins=20,
    labels={'purity': 'Pureza Tumoral', 'count': 'Número de Amostras'},
    title='Distribuição de Pureza Tumoral no Cohort'
)
fig_purity.update_traces(marker_color='#667eea')
fig_purity.update_layout(xaxis_tickformat='.0%')

# 2. Status MRD
mrd_counts = mrd_df['mrd_status'].value_counts().reset_index()
mrd_counts.columns = ['status', 'count']
colors_mrd = {'NEGATIVE': '#28a745', 'LOW': '#ffc107', 'MODERATE': '#fd7e14', 'HIGH': '#dc3545'}
fig_mrd = px.pie(
    mrd_counts, names='status', values='count',
    title='Status MRD - Distribuição no Cohort',
    color='status',
    color_discrete_map=colors_mrd
)

# 3. Tipos de leucemia
if not classification_df.empty:
    leuk_counts = classification_df['leukemia_type'].value_counts().reset_index()
    leuk_counts.columns = ['type', 'count']
    fig_leukemia = px.bar(
        leuk_counts, x='type', y='count',
        labels={'type': 'Tipo de Leucemia', 'count': 'Número de Amostras'},
        title='Distribuição de Tipos de Leucemia'
    )
    fig_leukemia.update_traces(marker_color='#764ba2')
else:
    fig_leukemia = None

# 4. Categorias de risco
if not classification_df.empty:
    risk_counts = classification_df['risk_category'].value_counts().reset_index()
    risk_counts.columns = ['risk', 'count']
    colors_risk = {'favorable': '#28a745', 'intermediate': '#ffc107', 'adverse': '#dc3545'}
    fig_risk = px.pie(
        risk_counts, names='risk', values='count',
        title='Estratificação de Risco',
        color='risk',
        color_discrete_map=colors_risk
    )
else:
    fig_risk = None

# 5. Genes driver mais frequentes
if not drivers_df.empty:
    driver_counts = drivers_df[drivers_df['impact'].isin(['HIGH', 'MODERATE'])].groupby('gene').size().reset_index(name='count')
    driver_counts = driver_counts.sort_values('count', ascending=False).head(15)
    fig_drivers = px.bar(
        driver_counts, x='count', y='gene', orientation='h',
        labels={'count': 'Número de Amostras', 'gene': 'Gene'},
        title='Top 15 Genes Driver Mais Frequentes'
    )
    fig_drivers.update_traces(marker_color='#e74c3c')
else:
    fig_drivers = None

# 6. Fusões patognomônicas
if not fusions_df.empty:
    pathog_fusions = fusions_df[fusions_df['pathognomonic'] == 1].copy()
    if not pathog_fusions.empty:
        pathog_fusions['fusion'] = pathog_fusions['gene1'] + '-' + pathog_fusions['gene2']
        fusion_counts = pathog_fusions.groupby('fusion').size().reset_index(name='count')
        fusion_counts = fusion_counts.sort_values('count', ascending=False)
        fig_fusions = px.bar(
            fusion_counts, x='count', y='fusion', orientation='h',
            labels={'count': 'Número de Amostras', 'fusion': 'Fusão'},
            title='Fusões Patognomônicas Detectadas'
        )
        fig_fusions.update_traces(marker_color='#9b59b6')
    else:
        fig_fusions = None
else:
    fig_fusions = None

# 7. Pureza vs Ploidy
if not purecn_df.empty:
    fig_purity_ploidy = px.scatter(
        purecn_df, x='purity', y='ploidy',
        labels={'purity': 'Pureza Tumoral', 'ploidy': 'Ploidy'},
        title='Pureza vs Ploidy',
        hover_data=['sample_id']
    )
    fig_purity_ploidy.update_traces(marker=dict(size=10, color='#667eea', line=dict(width=1, color='white')))
    fig_purity_ploidy.update_layout(xaxis_tickformat='.0%')
else:
    fig_purity_ploidy = None

# 8. Timeline de análises
samples_timeline = samples_df.copy()
samples_timeline['analysis_date'] = pd.to_datetime(samples_timeline['analysis_date'])
samples_timeline = samples_timeline.sort_values('analysis_date')
samples_timeline['cumulative'] = range(1, len(samples_timeline) + 1)
fig_timeline = px.line(
    samples_timeline, x='analysis_date', y='cumulative',
    labels={'analysis_date': 'Data', 'cumulative': 'Amostras Acumuladas'},
    title='Timeline de Análises'
)
fig_timeline.update_traces(line_color='#667eea', line_width=3)

# ========== GERA HTML ==========

html_content = f"""
<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VKaryo Cohort Summary</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            padding: 20px;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 50px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 3em;
            margin-bottom: 15px;
        }}
        .header .subtitle {{
            font-size: 1.3em;
            opacity: 0.95;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            padding: 40px;
            background: #f8f9fa;
        }}
        .stat-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
            border-left: 4px solid #667eea;
        }}
        .stat-card .label {{
            font-size: 0.9em;
            color: #666;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 10px;
        }}
        .stat-card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #667eea;
        }}
        .section {{
            padding: 40px;
            border-bottom: 1px solid #e0e0e0;
        }}
        .section h2 {{
            color: #667eea;
            margin-bottom: 30px;
            font-size: 2em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 15px;
        }}
        .plot-container {{
            margin: 30px 0;
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .grid-2 {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }}
        .footer {{
            text-align: center;
            padding: 30px;
            background: #2c3e50;
            color: white;
            font-size: 0.95em;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        table th {{
            background: #667eea;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
        }}
        table tr:hover {{
            background: #f8f9fa;
        }}
        .badge {{
            display: inline-block;
            padding: 5px 12px;
            border-radius: 15px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        .badge.favorable {{ background: #28a745; color: white; }}
        .badge.intermediate {{ background: #ffc107; color: #333; }}
        .badge.adverse {{ background: #dc3545; color: white; }}
    </style>
</head>
<body>
    <div class="container">
        <!-- HEADER -->
        <div class="header">
            <h1>🧬 VKaryo Cohort Summary</h1>
            <div class="subtitle">Relatório Consolidado do Cohort</div>
            <div style="margin-top: 25px; font-size: 1.1em;">
                <strong>Total de Amostras:</strong> {n_samples}<br>
                <strong>Gerado em:</strong> {datetime.now().strftime('%d/%m/%Y às %H:%M')}
            </div>
        </div>

        <!-- ESTATÍSTICAS GERAIS -->
        <div class="stats-grid">
            <div class="stat-card">
                <div class="label">Total Amostras</div>
                <div class="value">{n_samples}</div>
            </div>
            <div class="stat-card">
                <div class="label">Pureza Média</div>
                <div class="value">{purecn_df['purity'].mean()*100:.1f}%</div>
            </div>
            <div class="stat-card">
                <div class="label">MRD Negativo</div>
                <div class="value">{len(mrd_df[mrd_df['mrd_status'] == 'NEGATIVE'])}</div>
            </div>
            <div class="stat-card">
                <div class="label">Alto Risco</div>
                <div class="value">{len(classification_df[classification_df['risk_category'] == 'adverse'])}</div>
            </div>
            <div class="stat-card">
                <div class="label">Fusões Detectadas</div>
                <div class="value">{len(fusions_df)}</div>
            </div>
            <div class="stat-card">
                <div class="label">Drivers Únicos</div>
                <div class="value">{drivers_df['gene'].nunique()}</div>
            </div>
        </div>

        <!-- TIMELINE -->
        <div class="section">
            <h2>📅 Timeline de Análises</h2>
            <div class="plot-container">
                {plot_to_html(fig_timeline)}
            </div>
        </div>

        <!-- PUREZA E PLOIDY -->
        <div class="section">
            <h2>🔬 Análise de Pureza e Ploidy</h2>
            <div class="grid-2">
                <div class="plot-container">
                    {plot_to_html(fig_purity)}
                </div>
                <div class="plot-container">
                    {plot_to_html(fig_purity_ploidy) if fig_purity_ploidy else '<p>Dados não disponíveis</p>'}
                </div>
            </div>
        </div>

        <!-- MRD -->
        <div class="section">
            <h2>📉 Status de Doença Residual Mínima</h2>
            <div class="plot-container">
                {plot_to_html(fig_mrd)}
            </div>
        </div>

        <!-- CLASSIFICAÇÃO -->
        <div class="section">
            <h2>🎯 Classificação e Estratificação</h2>
            <div class="grid-2">
                <div class="plot-container">
                    {plot_to_html(fig_leukemia) if fig_leukemia else '<p>Dados de classificação não disponíveis</p>'}
                </div>
                <div class="plot-container">
                    {plot_to_html(fig_risk) if fig_risk else '<p>Dados de risco não disponíveis</p>'}
                </div>
            </div>
        </div>

        <!-- DRIVER MUTATIONS -->
        <div class="section">
            <h2>🎯 Landscape de Driver Mutations</h2>
            <div class="plot-container">
                {plot_to_html(fig_drivers) if fig_drivers else '<p>Nenhuma driver mutation detectada</p>'}
            </div>
        </div>

        <!-- FUSÕES -->
        {f'''
        <div class="section">
            <h2>🔗 Fusões Patognomônicas</h2>
            <div class="plot-container">
                {plot_to_html(fig_fusions)}
            </div>
        </div>
        ''' if fig_fusions else ''}

        <!-- TABELA DE AMOSTRAS -->
        <div class="section">
            <h2>📊 Detalhes das Amostras</h2>
            <table>
                <thead>
                    <tr>
                        <th>Amostra</th>
                        <th>Data</th>
                        <th>Pureza</th>
                        <th>Ploidy</th>
                        <th>MRD</th>
                        <th>Tipo</th>
                        <th>Risco</th>
                        <th>Drivers</th>
                        <th>Fusões</th>
                    </tr>
                </thead>
                <tbody>
"""

for _, row in samples_df.iterrows():
    risk_badge = f'<span class="badge {row.get("risk_category", "").lower()}">{row.get("risk_category", "N/A")}</span>'
    html_content += f"""
                    <tr>
                        <td><strong>{row['sample_id']}</strong></td>
                        <td>{row.get('analysis_date', 'N/A')}</td>
                        <td>{row.get('purity', 0)*100:.1f}%</td>
                        <td>{row.get('ploidy', 0):.2f}</td>
                        <td>{row.get('mrd_status', 'N/A')}</td>
                        <td>{row.get('leukemia_type', 'N/A')}</td>
                        <td>{risk_badge}</td>
                        <td>{row.get('n_driver_genes', 0)}</td>
                        <td>{row.get('n_fusions', 0)}</td>
                    </tr>
"""

html_content += """
                </tbody>
            </table>
        </div>

        <!-- FOOTER -->
        <div class="footer">
            <h3>🧬 VKaryo - Virtual Karyotyping for Leukemia Analysis</h3>
            <p style="margin-top: 15px;">
                Sistema completo de cariotipagem virtual para análise de leucemias<br>
                Integrando WES + RNA-Seq para diagnóstico molecular de precisão
            </p>
            <p style="margin-top: 20px; font-size: 0.9em; opacity: 0.8;">
                Desenvolvido com Snakemake | Database: SQLite | Dashboard: Streamlit
            </p>
        </div>
    </div>
</body>
</html>
"""

# Salva HTML
with open(output_html, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"✓ Cohort summary HTML generated: {output_html}")

# Gera PDF se solicitado
if output_pdf:
    try:
        import weasyprint
        weasyprint.HTML(string=html_content).write_pdf(output_pdf)
        print(f"✓ Cohort summary PDF generated: {output_pdf}")
    except ImportError:
        print("⚠ WeasyPrint not available, skipping PDF generation")
    except Exception as e:
        print(f"⚠ Could not generate PDF: {e}")

print(f"\n✓ Cohort analysis complete!")
print(f"  - Total samples: {n_samples}")
print(f"  - Average purity: {purecn_df['purity'].mean()*100:.1f}%")
print(f"  - MRD negative: {len(mrd_df[mrd_df['mrd_status'] == 'NEGATIVE'])}/{n_samples}")
