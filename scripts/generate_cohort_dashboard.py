#!/usr/bin/env python3
"""
Generate cohort dashboard with interactive visualizations
VKaryo - Virtual Karyotyping Workflow
"""

import os
import json
import glob
from datetime import datetime
from collections import defaultdict, Counter


def load_sample_summaries(summary_dir):
    """Load all sample JSON summaries"""
    summaries = []

    json_files = glob.glob(os.path.join(summary_dir, "*_summary.json"))

    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                summaries.append(data)
        except Exception as e:
            print(f"Warning: Could not load {json_file}: {e}")

    return summaries


def aggregate_statistics(summaries):
    """Calculate cohort-wide statistics"""
    stats = {
        'total_samples': len(summaries),
        'qc_pass': 0,
        'qc_fail': 0,
        'total_snvs': 0,
        'total_fusions': 0,
        'total_cnvs': 0,
        'avg_coverage': 0,
        'avg_purity': 0,
        'subtype_distribution': Counter(),
        'driver_mutation_frequency': Counter(),
        'fusion_frequency': Counter()
    }

    coverage_sum = 0
    purity_sum = 0
    purity_count = 0

    for sample in summaries:
        # QC status
        qc_status = sample.get('qc', {}).get('qc_status', 'UNKNOWN')
        if qc_status == 'PASS':
            stats['qc_pass'] += 1
        elif qc_status == 'FAIL':
            stats['qc_fail'] += 1

        # Alterations
        sample_stats = sample.get('statistics', {})
        stats['total_snvs'] += sample_stats.get('total_snvs', 0)
        stats['total_fusions'] += sample_stats.get('total_fusions', 0)
        stats['total_cnvs'] += sample_stats.get('total_cnvs', 0)

        # Coverage
        try:
            coverage = float(sample.get('qc', {}).get('mean_coverage', 0))
            coverage_sum += coverage
        except:
            pass

        # Purity
        try:
            purity = float(sample.get('qc', {}).get('tumor_purity', 0))
            if purity > 0:
                purity_sum += purity
                purity_count += 1
        except:
            pass

        # Classification
        subtype = sample.get('classification', {}).get('predicted_subtype', 'Unknown')
        stats['subtype_distribution'][subtype] += 1

    # Calculate averages
    if stats['total_samples'] > 0:
        stats['avg_coverage'] = coverage_sum / stats['total_samples']

    if purity_count > 0:
        stats['avg_purity'] = purity_sum / purity_count

    return stats


def generate_dashboard_html(summaries, stats, output_file):
    """Generate interactive HTML dashboard"""

    # Prepare data for JavaScript
    samples_data = []
    for sample in summaries:
        samples_data.append({
            'name': sample.get('sample', 'Unknown'),
            'qc': sample.get('qc', {}).get('qc_status', 'UNKNOWN'),
            'snvs': sample.get('statistics', {}).get('total_snvs', 0),
            'fusions': sample.get('statistics', {}).get('total_fusions', 0),
            'cnvs': sample.get('statistics', {}).get('total_cnvs', 0),
            'subtype': sample.get('classification', {}).get('predicted_subtype', 'Unknown'),
            'coverage': sample.get('qc', {}).get('mean_coverage', 'NA'),
            'purity': sample.get('qc', {}).get('tumor_purity', 'NA')
        })

    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VKaryo Cohort Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-2.18.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: #f5f5f5;
            padding: 20px;
        }}

        .container {{
            max-width: 1600px;
            margin: 0 auto;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .stat-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }}

        .stat-card h3 {{
            color: #667eea;
            font-size: 1.1em;
            margin-bottom: 10px;
        }}

        .stat-card .number {{
            font-size: 2.5em;
            font-weight: bold;
            color: #333;
        }}

        .chart-container {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}

        .chart-container h2 {{
            color: #667eea;
            margin-bottom: 20px;
            border-left: 4px solid #667eea;
            padding-left: 15px;
        }}

        table {{
            width: 100%;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        table th {{
            background: #667eea;
            color: white;
            padding: 15px;
            text-align: left;
        }}

        table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
        }}

        table tr:hover {{
            background: #f8f8f8;
        }}

        .badge {{
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
        }}

        .badge-pass {{
            background: #4caf50;
            color: white;
        }}

        .badge-fail {{
            background: #f44336;
            color: white;
        }}

        .filter-bar {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .filter-bar input {{
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 5px;
            width: 300px;
            font-size: 1em;
        }}

        .footer {{
            text-align: center;
            padding: 30px;
            color: #666;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>📊 VKaryo Cohort Dashboard</h1>
            <p>Análise de {stats['total_samples']} amostras | Gerado em {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        </div>

        <div class="stats-grid">
            <div class="stat-card">
                <h3>Total de Amostras</h3>
                <div class="number">{stats['total_samples']}</div>
            </div>
            <div class="stat-card">
                <h3>QC PASS</h3>
                <div class="number" style="color: #4caf50;">{stats['qc_pass']}</div>
            </div>
            <div class="stat-card">
                <h3>SNVs Detectados</h3>
                <div class="number">{stats['total_snvs']}</div>
            </div>
            <div class="stat-card">
                <h3>Fusões Detectadas</h3>
                <div class="number">{stats['total_fusions']}</div>
            </div>
            <div class="stat-card">
                <h3>CNVs Detectados</h3>
                <div class="number">{stats['total_cnvs']}</div>
            </div>
            <div class="stat-card">
                <h3>Cobertura Média</h3>
                <div class="number">{stats['avg_coverage']:.0f}X</div>
            </div>
        </div>

        <div class="chart-container">
            <h2>Distribuição de Subtipos</h2>
            <div id="subtype-chart"></div>
        </div>

        <div class="chart-container">
            <h2>Status de QC por Amostra</h2>
            <div id="qc-chart"></div>
        </div>

        <div class="chart-container">
            <h2>Alterações por Amostra</h2>
            <div id="alterations-chart"></div>
        </div>

        <div class="filter-bar">
            <input type="text" id="search-input" placeholder="Filtrar amostras..." onkeyup="filterTable()">
        </div>

        <table id="samples-table">
            <thead>
                <tr>
                    <th>Amostra</th>
                    <th>QC Status</th>
                    <th>SNVs</th>
                    <th>Fusões</th>
                    <th>CNVs</th>
                    <th>Subtipo</th>
                    <th>Cobertura</th>
                    <th>Pureza</th>
                    <th>Ação</th>
                </tr>
            </thead>
            <tbody>
"""

    # Add table rows
    for sample in samples_data:
        qc_badge = f'<span class="badge badge-pass">PASS</span>' if sample['qc'] == 'PASS' else f'<span class="badge badge-fail">{sample["qc"]}</span>'

        html += f"""
                <tr>
                    <td><strong>{sample['name']}</strong></td>
                    <td>{qc_badge}</td>
                    <td>{sample['snvs']}</td>
                    <td>{sample['fusions']}</td>
                    <td>{sample['cnvs']}</td>
                    <td>{sample['subtype']}</td>
                    <td>{sample['coverage']}X</td>
                    <td>{sample['purity']}%</td>
                    <td><a href="{sample['name']}_vkaryo_report.html" target="_blank">Ver Relatório</a></td>
                </tr>
"""

    # JavaScript data
    js_samples_data = json.dumps(samples_data)
    js_subtype_dist = json.dumps(dict(stats['subtype_distribution']))

    html += f"""
            </tbody>
        </table>

        <div class="footer">
            <p><strong>VKaryo</strong> - Virtual Karyotyping Workflow</p>
            <p>Gerado em {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        </div>
    </div>

    <script>
        const samplesData = {js_samples_data};
        const subtypeDistribution = {js_subtype_dist};

        // Subtype distribution pie chart
        const subtypeData = [{{
            values: Object.values(subtypeDistribution),
            labels: Object.keys(subtypeDistribution),
            type: 'pie',
            marker: {{
                colors: ['#667eea', '#764ba2', '#f093fb', '#4facfe', '#43e97b', '#fa709a']
            }}
        }}];

        const subtypeLayout = {{
            height: 400,
            margin: {{ t: 0, b: 0, l: 0, r: 0 }}
        }};

        Plotly.newPlot('subtype-chart', subtypeData, subtypeLayout);

        // QC status bar chart
        const qcCounts = samplesData.reduce((acc, s) => {{
            acc[s.qc] = (acc[s.qc] || 0) + 1;
            return acc;
        }}, {{}});

        const qcData = [{{
            x: Object.keys(qcCounts),
            y: Object.values(qcCounts),
            type: 'bar',
            marker: {{
                color: Object.keys(qcCounts).map(k => k === 'PASS' ? '#4caf50' : '#f44336')
            }}
        }}];

        const qcLayout = {{
            height: 400,
            xaxis: {{ title: 'Status' }},
            yaxis: {{ title: 'Número de Amostras' }},
            margin: {{ t: 20 }}
        }};

        Plotly.newPlot('qc-chart', qcData, qcLayout);

        // Alterations scatter plot
        const alterationsData = [{{
            x: samplesData.map(s => s.snvs),
            y: samplesData.map(s => s.fusions),
            mode: 'markers',
            type: 'scatter',
            text: samplesData.map(s => s.name),
            marker: {{
                size: samplesData.map(s => s.cnvs + 5),
                color: samplesData.map(s => s.qc === 'PASS' ? '#4caf50' : '#f44336')
            }}
        }}];

        const alterationsLayout = {{
            height: 400,
            xaxis: {{ title: 'SNVs' }},
            yaxis: {{ title: 'Fusões' }},
            hovermode: 'closest',
            margin: {{ t: 20 }}
        }};

        Plotly.newPlot('alterations-chart', alterationsData, alterationsLayout);

        // Table filter function
        function filterTable() {{
            const input = document.getElementById('search-input');
            const filter = input.value.toUpperCase();
            const table = document.getElementById('samples-table');
            const tr = table.getElementsByTagName('tr');

            for (let i = 1; i < tr.length; i++) {{
                const td = tr[i].getElementsByTagName('td')[0];
                if (td) {{
                    const txtValue = td.textContent || td.innerText;
                    if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                        tr[i].style.display = '';
                    }} else {{
                        tr[i].style.display = 'none';
                    }}
                }}
            }}
        }}
    </script>
</body>
</html>
"""

    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html)

    print(f"✓ Cohort dashboard saved: {output_file}")


def main():
    # For Snakemake integration
    if 'snakemake' in globals():
        summary_files = snakemake.input.summaries
        output_file = snakemake.output.dashboard

        # Load summaries
        summaries = []
        for json_file in summary_files:
            with open(json_file, 'r') as f:
                summaries.append(json.load(f))

    else:
        # Standalone usage
        import sys
        if len(sys.argv) < 3:
            print("Usage: python generate_cohort_dashboard.py <summary_dir> <output_file>")
            sys.exit(1)

        summary_dir = sys.argv[1]
        output_file = sys.argv[2]

        summaries = load_sample_summaries(summary_dir)

    # Aggregate statistics
    stats = aggregate_statistics(summaries)

    # Generate dashboard
    generate_dashboard_html(summaries, stats, output_file)

    print(f"Dashboard generated successfully!")
    print(f"Total samples: {stats['total_samples']}")
    print(f"QC PASS: {stats['qc_pass']}")
    print(f"Total SNVs: {stats['total_snvs']}")
    print(f"Total Fusions: {stats['total_fusions']}")


if __name__ == "__main__":
    main()
