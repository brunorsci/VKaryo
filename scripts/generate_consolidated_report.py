#!/usr/bin/env python3
"""
Gera relatório consolidado HTML/PDF com todos os resultados VKaryo
Integra: PureCN, OptiType, VarDict MRD, Mutect2, CNVs, Fusões, Classificação
"""

import sqlite3
import pandas as pd
from pathlib import Path
from datetime import datetime
import base64
import json

# Parâmetros do Snakemake
sample_id = snakemake.wildcards.sample
db_path = snakemake.input.db
output_html = snakemake.output.html
output_pdf = snakemake.output.pdf if hasattr(snakemake.output, 'pdf') else None

# Plots (se disponíveis)
plots = {}
if hasattr(snakemake.input, 'purecn_plot'):
    plots['purecn'] = snakemake.input.purecn_plot
if hasattr(snakemake.input, 'mrd_plot'):
    plots['mrd'] = snakemake.input.mrd_plot
if hasattr(snakemake.input, 'aneuploidy_plot'):
    plots['aneuploidy'] = snakemake.input.aneuploidy_plot

def img_to_base64(img_path):
    """Converte imagem para base64 para embed no HTML"""
    if not Path(img_path).exists():
        return None
    with open(img_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def query_db(db_path, query, params=()):
    """Query no database"""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute(query, params)
    rows = cursor.fetchall()
    conn.close()
    return rows

def get_sample_overview(db_path, sample_id):
    """Busca overview da amostra"""
    query = "SELECT * FROM v_sample_summary WHERE sample_id = ?"
    rows = query_db(db_path, query, (sample_id,))
    return dict(rows[0]) if rows else {}

def get_purecn_data(db_path, sample_id):
    """Busca dados PureCN"""
    query = "SELECT * FROM purecn_results WHERE sample_id = ?"
    rows = query_db(db_path, query, (sample_id,))
    return dict(rows[0]) if rows else {}

def get_hla_data(db_path, sample_id):
    """Busca dados HLA"""
    query = "SELECT * FROM optitype_results WHERE sample_id = ?"
    rows = query_db(db_path, query, (sample_id,))
    return dict(rows[0]) if rows else {}

def get_mrd_data(db_path, sample_id):
    """Busca dados MRD"""
    query = "SELECT * FROM mrd_results WHERE sample_id = ?"
    rows = query_db(db_path, query, (sample_id,))
    mrd_summary = dict(rows[0]) if rows else {}

    # Variantes MRD
    query = "SELECT * FROM mrd_variants WHERE sample_id = ? ORDER BY vaf DESC"
    variants = [dict(row) for row in query_db(db_path, query, (sample_id,))]

    return mrd_summary, variants

def get_driver_mutations(db_path, sample_id):
    """Busca driver mutations"""
    query = """
        SELECT * FROM driver_mutations
        WHERE sample_id = ?
        ORDER BY
            CASE impact WHEN 'HIGH' THEN 1 WHEN 'MODERATE' THEN 2 ELSE 3 END,
            gene
    """
    return [dict(row) for row in query_db(db_path, query, (sample_id,))]

def get_cnvs(db_path, sample_id):
    """Busca CNVs"""
    query = """
        SELECT * FROM cnv_calls
        WHERE sample_id = ? AND cn_status != 'neutral'
        ORDER BY ABS(log2_ratio) DESC
        LIMIT 50
    """
    return [dict(row) for row in query_db(db_path, query, (sample_id,))]

def get_fusions(db_path, sample_id):
    """Busca fusões"""
    query = "SELECT * FROM gene_fusions WHERE sample_id = ? ORDER BY confidence DESC"
    return [dict(row) for row in query_db(db_path, query, (sample_id,))]

def get_classification(db_path, sample_id):
    """Busca classificação"""
    query = "SELECT * FROM leukemia_classification WHERE sample_id = ?"
    rows = query_db(db_path, query, (sample_id,))
    return dict(rows[0]) if rows else {}

def generate_html_report(data, plots):
    """Gera relatório HTML"""

    html = f"""
<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VKaryo Report - {data['sample_id']}</title>
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
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        .header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        .section {{
            padding: 30px 40px;
            border-bottom: 1px solid #e0e0e0;
        }}
        .section:last-child {{
            border-bottom: none;
        }}
        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        .section h3 {{
            color: #555;
            margin: 20px 0 10px 0;
            font-size: 1.3em;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .card .label {{
            font-size: 0.9em;
            color: #666;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        .card .value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #333;
            margin-top: 5px;
        }}
        .badge {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: bold;
            text-transform: uppercase;
        }}
        .badge.high {{ background: #dc3545; color: white; }}
        .badge.moderate {{ background: #fd7e14; color: white; }}
        .badge.low {{ background: #ffc107; color: #333; }}
        .badge.negative {{ background: #28a745; color: white; }}
        .badge.favorable {{ background: #28a745; color: white; }}
        .badge.intermediate {{ background: #ffc107; color: #333; }}
        .badge.adverse {{ background: #dc3545; color: white; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
        }}
        table th {{
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        table td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ddd;
        }}
        table tr:hover {{
            background: #f8f9fa;
        }}
        .plot {{
            margin: 20px 0;
            text-align: center;
        }}
        .plot img {{
            max-width: 100%;
            border: 1px solid #ddd;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .alert {{
            padding: 15px;
            margin: 20px 0;
            border-radius: 8px;
            border-left: 4px solid;
        }}
        .alert.info {{
            background: #e3f2fd;
            border-color: #2196f3;
            color: #0d47a1;
        }}
        .alert.warning {{
            background: #fff3cd;
            border-color: #ffc107;
            color: #856404;
        }}
        .alert.danger {{
            background: #f8d7da;
            border-color: #dc3545;
            color: #721c24;
        }}
        .alert.success {{
            background: #d4edda;
            border-color: #28a745;
            color: #155724;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            background: #f8f9fa;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- HEADER -->
        <div class="header">
            <h1>🧬 VKaryo - Cariotipagem Virtual</h1>
            <div class="subtitle">Relatório Consolidado de Análise</div>
            <div style="margin-top: 20px;">
                <strong>Amostra:</strong> {data['sample_id']}<br>
                <strong>Data de Análise:</strong> {data.get('overview', {}).get('analysis_date', datetime.now().strftime('%Y-%m-%d'))}
            </div>
        </div>

        <!-- OVERVIEW -->
        <div class="section">
            <h2>📊 Visão Geral</h2>
            <div class="grid">
                <div class="card">
                    <div class="label">Pureza Tumoral</div>
                    <div class="value">{data.get('purecn', {}).get('purity', 0)*100:.1f}%</div>
                </div>
                <div class="card">
                    <div class="label">Ploidy</div>
                    <div class="value">{data.get('purecn', {}).get('ploidy', 0):.2f}</div>
                </div>
                <div class="card">
                    <div class="label">MRD Status</div>
                    <div class="value">
                        <span class="badge {data.get('mrd_summary', {}).get('mrd_status', 'unknown').lower()}">
                            {data.get('mrd_summary', {}).get('mrd_status', 'N/A')}
                        </span>
                    </div>
                </div>
                <div class="card">
                    <div class="label">Driver Mutations</div>
                    <div class="value">{len(data.get('drivers', []))}</div>
                </div>
            </div>
        </div>

        <!-- CLASSIFICAÇÃO -->
        {generate_classification_section(data.get('classification', {}))}

        <!-- PURECN -->
        {generate_purecn_section(data.get('purecn', {}), plots.get('purecn'))}

        <!-- HLA TYPING -->
        {generate_hla_section(data.get('hla', {}))}

        <!-- MRD -->
        {generate_mrd_section(data.get('mrd_summary', {}), data.get('mrd_variants', []), plots.get('mrd'))}

        <!-- DRIVER MUTATIONS -->
        {generate_drivers_section(data.get('drivers', []))}

        <!-- CNVs -->
        {generate_cnv_section(data.get('cnvs', []), plots.get('aneuploidy'))}

        <!-- FUSÕES -->
        {generate_fusions_section(data.get('fusions', []))}

        <!-- FOOTER -->
        <div class="footer">
            <p><strong>VKaryo</strong> - Virtual Karyotyping for Leukemia Analysis</p>
            <p>Gerado automaticamente em {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
"""
    return html

def generate_classification_section(classification):
    """Gera seção de classificação"""
    if not classification:
        return ""

    risk_badge = classification.get('risk_category', 'unknown').lower()
    return f"""
        <div class="section">
            <h2>🎯 Classificação e Estratificação de Risco</h2>
            <div class="alert info">
                <strong>Tipo de Leucemia:</strong> {classification.get('leukemia_type', 'N/A')}<br>
                <strong>Subtipo:</strong> {classification.get('subtype', 'N/A')}<br>
                <strong>Categoria de Risco:</strong>
                <span class="badge {risk_badge}">{classification.get('risk_category', 'N/A')}</span>
            </div>
            <h3>Terapia Recomendada</h3>
            <p>{classification.get('recommended_therapy', 'Consultar protocolo institucional')}</p>
        </div>
    """

def generate_purecn_section(purecn, plot_path):
    """Gera seção PureCN"""
    plot_html = ""
    if plot_path and Path(plot_path).exists():
        img_b64 = img_to_base64(plot_path)
        if img_b64:
            plot_html = f'<div class="plot"><img src="data:image/png;base64,{img_b64}"></div>'

    return f"""
        <div class="section">
            <h2>🔬 Pureza Tumoral e Ploidy (PureCN)</h2>
            <div class="grid">
                <div class="card">
                    <div class="label">Pureza</div>
                    <div class="value">{purecn.get('purity', 0)*100:.1f}%</div>
                </div>
                <div class="card">
                    <div class="label">Ploidy</div>
                    <div class="value">{purecn.get('ploidy', 0):.2f}</div>
                </div>
                <div class="card">
                    <div class="label">Amplificações</div>
                    <div class="value">{purecn.get('n_amplifications', 0)}</div>
                </div>
                <div class="card">
                    <div class="label">Deleções</div>
                    <div class="value">{purecn.get('n_deletions', 0)}</div>
                </div>
                <div class="card">
                    <div class="label">LOH Events</div>
                    <div class="value">{purecn.get('n_loh', 0)}</div>
                </div>
                <div class="card">
                    <div class="label">Contaminação</div>
                    <div class="value">{purecn.get('contamination', 0):.4f}</div>
                </div>
            </div>
            {plot_html}
        </div>
    """

def generate_hla_section(hla):
    """Gera seção HLA"""
    if not hla:
        return ""

    return f"""
        <div class="section">
            <h2>🧬 HLA Typing (OptiType)</h2>
            <div class="alert info">
                <strong>Relevância:</strong> Essencial para matching de transplante alogênico e avaliação de risco de GvHD
            </div>
            <table>
                <tr>
                    <th>Locus</th>
                    <th>Alelo 1</th>
                    <th>Alelo 2</th>
                </tr>
                <tr>
                    <td><strong>HLA-A</strong></td>
                    <td>{hla.get('hla_a1', 'N/A')}</td>
                    <td>{hla.get('hla_a2', 'N/A')}</td>
                </tr>
                <tr>
                    <td><strong>HLA-B</strong></td>
                    <td>{hla.get('hla_b1', 'N/A')}</td>
                    <td>{hla.get('hla_b2', 'N/A')}</td>
                </tr>
                <tr>
                    <td><strong>HLA-C</strong></td>
                    <td>{hla.get('hla_c1', 'N/A')}</td>
                    <td>{hla.get('hla_c2', 'N/A')}</td>
                </tr>
            </table>
            <p><strong>Reads analisados:</strong> {hla.get('reads', 'N/A'):,} | <strong>Objective:</strong> {hla.get('objective', 0):.2f}</p>
        </div>
    """

def generate_mrd_section(mrd_summary, mrd_variants, plot_path):
    """Gera seção MRD"""
    status = mrd_summary.get('mrd_status', 'UNKNOWN')
    status_class = status.lower()

    alert_class = {
        'NEGATIVE': 'success',
        'LOW': 'info',
        'MODERATE': 'warning',
        'HIGH': 'danger'
    }.get(status, 'info')

    plot_html = ""
    if plot_path and Path(plot_path).exists():
        img_b64 = img_to_base64(plot_path)
        if img_b64:
            plot_html = f'<div class="plot"><img src="data:image/png;base64,{img_b64}"></div>'

    variants_html = ""
    if mrd_variants:
        variants_rows = ""
        for var in mrd_variants[:10]:  # Top 10
            variants_rows += f"""
            <tr>
                <td>{var.get('gene', 'N/A')}</td>
                <td>{var.get('chromosome', '')}:{var.get('position', '')}</td>
                <td>{var.get('ref', '')}&gt;{var.get('alt', '')}</td>
                <td>{var.get('vaf', 0):.3f}%</td>
                <td>{var.get('cosmic_id', 'N/A')}</td>
            </tr>
            """
        variants_html = f"""
            <h3>Variantes MRD Detectadas (Top 10)</h3>
            <table>
                <tr>
                    <th>Gene</th>
                    <th>Posição</th>
                    <th>Alteração</th>
                    <th>VAF</th>
                    <th>COSMIC</th>
                </tr>
                {variants_rows}
            </table>
        """

    return f"""
        <div class="section">
            <h2>📉 Doença Residual Mínima (MRD)</h2>
            <div class="alert {alert_class}">
                <strong>Status MRD:</strong> <span class="badge {status_class}">{status}</span><br>
                <strong>VAF Máxima Detectada:</strong> {mrd_summary.get('max_vaf', 0):.3f}%<br>
                <strong>Variantes MRD Relevantes:</strong> {mrd_summary.get('mrd_variants', 0)} / {mrd_summary.get('total_variants', 0)} totais
            </div>
            {plot_html}
            {variants_html}
        </div>
    """

def generate_drivers_section(drivers):
    """Gera seção de driver mutations"""
    if not drivers:
        return """
        <div class="section">
            <h2>🎯 Driver Mutations</h2>
            <div class="alert success">
                Nenhuma driver mutation detectada em genes conhecidos de leucemia.
            </div>
        </div>
        """

    rows = ""
    for driver in drivers[:20]:  # Top 20
        impact_class = driver.get('impact', 'low').lower()
        rows += f"""
        <tr>
            <td>{driver.get('gene', 'N/A')}</td>
            <td>{driver.get('chromosome', '')}:{driver.get('position', '')}</td>
            <td>{driver.get('ref', '')}&gt;{driver.get('alt', '')}</td>
            <td><span class="badge {impact_class}">{driver.get('impact', 'N/A')}</span></td>
            <td>{driver.get('exonic_func', 'N/A')}</td>
            <td>{driver.get('aa_change', 'N/A')}</td>
        </tr>
        """

    return f"""
        <div class="section">
            <h2>🎯 Driver Mutations (Mutect2)</h2>
            <div class="alert warning">
                <strong>{len(drivers)} driver mutations</strong> detectadas em genes recorrentes de leucemia
            </div>
            <table>
                <tr>
                    <th>Gene</th>
                    <th>Posição</th>
                    <th>Alteração</th>
                    <th>Impacto</th>
                    <th>Função</th>
                    <th>Mudança AA</th>
                </tr>
                {rows}
            </table>
        </div>
    """

def generate_cnv_section(cnvs, plot_path):
    """Gera seção CNVs"""
    plot_html = ""
    if plot_path and Path(plot_path).exists():
        img_b64 = img_to_base64(plot_path)
        if img_b64:
            plot_html = f'<div class="plot"><img src="data:image/png;base64,{img_b64}"></div>'

    if not cnvs:
        return f"""
        <div class="section">
            <h2>📊 Alterações no Número de Cópias (CNVs)</h2>
            <div class="alert info">
                Nenhuma alteração significativa detectada.
            </div>
            {plot_html}
        </div>
        """

    rows = ""
    for cnv in cnvs[:30]:  # Top 30
        status_badge = 'high' if cnv.get('cn_status') == 'gain' else 'moderate'
        validated = "✓" if cnv.get('validated') else ""
        rows += f"""
        <tr>
            <td>{cnv.get('chromosome', 'N/A')}</td>
            <td>{cnv.get('start', 0):,} - {cnv.get('end', 0):,}</td>
            <td>{cnv.get('gene', 'N/A')}</td>
            <td>{cnv.get('log2_ratio', 0):.2f}</td>
            <td><span class="badge {status_badge}">{cnv.get('cn_status', 'N/A')}</span></td>
            <td>{validated}</td>
        </tr>
        """

    return f"""
        <div class="section">
            <h2>📊 Alterações no Número de Cópias (CNVs)</h2>
            {plot_html}
            <h3>CNVs Detectados (Top 30)</h3>
            <table>
                <tr>
                    <th>Cromossomo</th>
                    <th>Região</th>
                    <th>Gene</th>
                    <th>Log2 Ratio</th>
                    <th>Status</th>
                    <th>Validado</th>
                </tr>
                {rows}
            </table>
        </div>
    """

def generate_fusions_section(fusions):
    """Gera seção de fusões"""
    if not fusions:
        return ""

    rows = ""
    for fusion in fusions:
        conf_badge = fusion.get('confidence', 'low').lower()
        pathog = "⭐" if fusion.get('pathognomonic') else ""
        rows += f"""
        <tr>
            <td>{pathog} {fusion.get('gene1', '')}-{fusion.get('gene2', '')}</td>
            <td>{fusion.get('breakpoint1', 'N/A')}</td>
            <td>{fusion.get('breakpoint2', 'N/A')}</td>
            <td><span class="badge {conf_badge}">{fusion.get('confidence', 'N/A')}</span></td>
            <td>{fusion.get('prognosis', 'N/A')}</td>
            <td>{fusion.get('therapy', 'N/A')}</td>
        </tr>
        """

    return f"""
        <div class="section">
            <h2>🔗 Fusões de Genes</h2>
            <table>
                <tr>
                    <th>Fusão</th>
                    <th>Breakpoint 1</th>
                    <th>Breakpoint 2</th>
                    <th>Confiança</th>
                    <th>Prognóstico</th>
                    <th>Terapia</th>
                </tr>
                {rows}
            </table>
        </div>
    """

# Coleta dados do database
overview = get_sample_overview(db_path, sample_id)
purecn = get_purecn_data(db_path, sample_id)
hla = get_hla_data(db_path, sample_id)
mrd_summary, mrd_variants = get_mrd_data(db_path, sample_id)
drivers = get_driver_mutations(db_path, sample_id)
cnvs = get_cnvs(db_path, sample_id)
fusions = get_fusions(db_path, sample_id)
classification = get_classification(db_path, sample_id)

# Consolida dados
data = {
    'sample_id': sample_id,
    'overview': overview,
    'purecn': purecn,
    'hla': hla,
    'mrd_summary': mrd_summary,
    'mrd_variants': mrd_variants,
    'drivers': drivers,
    'cnvs': cnvs,
    'fusions': fusions,
    'classification': classification
}

# Gera HTML
html_content = generate_html_report(data, plots)

# Salva HTML
with open(output_html, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"✓ HTML report generated: {output_html}")

# Gera PDF se solicitado
if output_pdf:
    try:
        import weasyprint
        weasyprint.HTML(string=html_content).write_pdf(output_pdf)
        print(f"✓ PDF report generated: {output_pdf}")
    except ImportError:
        print("⚠ WeasyPrint not available, skipping PDF generation")
    except Exception as e:
        print(f"⚠ Could not generate PDF: {e}")
