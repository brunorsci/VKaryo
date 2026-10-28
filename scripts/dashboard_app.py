#!/usr/bin/env python3
"""
Dashboard interativo VKaryo usando Streamlit
Visualização dinâmica de resultados de cariotipagem virtual
"""

import streamlit as st
import sqlite3
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import json

# Configuração da página
st.set_page_config(
    page_title="VKaryo Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS customizado
st.markdown("""
<style>
    .main {background-color: #f5f5f5;}
    .stMetric {background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);}
    h1 {color: #667eea;}
    h2 {color: #764ba2;}
    .stAlert {border-radius: 8px;}
</style>
""", unsafe_allow_html=True)

@st.cache_resource
def get_db_connection(db_path):
    """Conecta ao database"""
    return sqlite3.connect(db_path, check_same_thread=False)

@st.cache_data
def load_samples(_conn):
    """Carrega lista de amostras"""
    query = """
        SELECT
            s.*,
            p.purity,
            p.ploidy,
            m.mrd_status,
            l.leukemia_type,
            l.risk_category
        FROM samples s
        LEFT JOIN purecn_results p ON s.sample_id = p.sample_id
        LEFT JOIN mrd_results m ON s.sample_id = m.sample_id
        LEFT JOIN leukemia_classification l ON s.sample_id = l.sample_id
        ORDER BY s.analysis_date DESC
    """
    return pd.read_sql_query(query, _conn)

@st.cache_data
def load_sample_details(_conn, sample_id):
    """Carrega detalhes de uma amostra"""
    details = {}

    # PureCN
    query = "SELECT * FROM purecn_results WHERE sample_id = ?"
    details['purecn'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # HLA
    query = "SELECT * FROM optitype_results WHERE sample_id = ?"
    details['hla'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # MRD
    query = "SELECT * FROM mrd_results WHERE sample_id = ?"
    details['mrd'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    query = "SELECT * FROM mrd_variants WHERE sample_id = ? ORDER BY vaf DESC"
    details['mrd_variants'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # Drivers
    query = "SELECT * FROM driver_mutations WHERE sample_id = ? ORDER BY vaf DESC"
    details['drivers'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # CNVs
    query = "SELECT * FROM cnv_calls WHERE sample_id = ? AND cn_status != 'neutral'"
    details['cnvs'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # Fusões
    query = "SELECT * FROM gene_fusions WHERE sample_id = ?"
    details['fusions'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    # Classificação
    query = "SELECT * FROM leukemia_classification WHERE sample_id = ?"
    details['classification'] = pd.read_sql_query(query, _conn, params=(sample_id,))

    return details

def plot_purity_distribution(conn):
    """Gráfico de distribuição de pureza"""
    df = pd.read_sql_query("SELECT sample_id, purity FROM purecn_results", conn)
    if df.empty:
        return None

    fig = px.histogram(df, x='purity', nbins=20,
                      labels={'purity': 'Pureza Tumoral', 'count': 'Número de Amostras'},
                      title='Distribuição de Pureza Tumoral')
    fig.update_traces(marker_color='#667eea')
    return fig

def plot_mrd_status(conn):
    """Gráfico de status MRD"""
    df = pd.read_sql_query("SELECT mrd_status, COUNT(*) as count FROM mrd_results GROUP BY mrd_status", conn)
    if df.empty:
        return None

    colors = {'NEGATIVE': '#28a745', 'LOW': '#ffc107', 'MODERATE': '#fd7e14', 'HIGH': '#dc3545'}
    fig = px.pie(df, names='mrd_status', values='count',
                title='Status MRD por Amostra',
                color='mrd_status',
                color_discrete_map=colors)
    return fig

def plot_driver_genes(conn):
    """Gráfico de genes driver mais frequentes"""
    df = pd.read_sql_query("""
        SELECT gene, COUNT(*) as count
        FROM driver_mutations
        WHERE impact IN ('HIGH', 'MODERATE')
        GROUP BY gene
        ORDER BY count DESC
        LIMIT 15
    """, conn)
    if df.empty:
        return None

    fig = px.bar(df, x='count', y='gene', orientation='h',
                labels={'count': 'Número de Amostras', 'gene': 'Gene'},
                title='Genes Driver Mais Frequentes')
    fig.update_traces(marker_color='#764ba2')
    return fig

def plot_cnv_distribution(conn, sample_id):
    """Gráfico de distribuição de CNVs por cromossomo"""
    df = pd.read_sql_query(f"""
        SELECT chromosome, cn_status, COUNT(*) as count
        FROM cnv_calls
        WHERE sample_id = '{sample_id}'
        GROUP BY chromosome, cn_status
    """, conn)
    if df.empty:
        return None

    # Ordena cromossomos
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df['chromosome'] = pd.Categorical(df['chromosome'], categories=chr_order, ordered=True)
    df = df.sort_values('chromosome')

    colors = {'gain': '#dc3545', 'loss': '#007bff', 'neutral': '#6c757d'}
    fig = px.bar(df, x='chromosome', y='count', color='cn_status',
                labels={'count': 'Número de CNVs', 'chromosome': 'Cromossomo', 'cn_status': 'Status'},
                title=f'Distribuição de CNVs - {sample_id}',
                color_discrete_map=colors,
                barmode='stack')
    return fig

def plot_vaf_distribution(conn, sample_id):
    """Gráfico de distribuição de VAF das variantes"""
    df = pd.read_sql_query(f"""
        SELECT gene, vaf
        FROM driver_mutations
        WHERE sample_id = '{sample_id}'
        ORDER BY vaf DESC
        LIMIT 20
    """, conn)
    if df.empty:
        return None

    fig = px.scatter(df, x='gene', y='vaf', size='vaf',
                    labels={'vaf': 'VAF (%)', 'gene': 'Gene'},
                    title=f'VAF das Driver Mutations - {sample_id}',
                    color='vaf',
                    color_continuous_scale='Reds')
    fig.update_layout(xaxis_tickangle=-45)
    return fig

# ==================== MAIN APP ====================

st.title("🧬 VKaryo Dashboard")
st.markdown("### Dashboard Interativo de Cariotipagem Virtual para Leucemias")

# Sidebar
st.sidebar.title("📊 Navegação")

# Busca database
db_path = st.sidebar.text_input(
    "Caminho do Database",
    value="output_vkaryo/database/vkaryo_results.db",
    help="Caminho para o arquivo SQLite do database"
)

if not Path(db_path).exists():
    st.error(f"❌ Database não encontrado: {db_path}")
    st.info("Por favor, execute o workflow VKaryo primeiro para gerar o database.")
    st.stop()

# Conecta ao database
conn = get_db_connection(db_path)

# Carrega amostras
samples_df = load_samples(conn)

if samples_df.empty:
    st.warning("⚠️ Nenhuma amostra encontrada no database.")
    st.stop()

# Seleção de view
view_mode = st.sidebar.radio(
    "Modo de Visualização",
    ["📋 Visão Geral", "🔍 Detalhes da Amostra", "📊 Análises Comparativas"]
)

# ==================== VISÃO GERAL ====================

if view_mode == "📋 Visão Geral":
    st.header("📋 Visão Geral - Todas as Amostras")

    # Métricas gerais
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total de Amostras", len(samples_df))

    with col2:
        mrd_negative = len(samples_df[samples_df['mrd_status'] == 'NEGATIVE'])
        st.metric("MRD Negativo", mrd_negative)

    with col3:
        avg_purity = samples_df['purity'].mean() * 100 if 'purity' in samples_df else 0
        st.metric("Pureza Média", f"{avg_purity:.1f}%")

    with col4:
        high_risk = len(samples_df[samples_df['risk_category'] == 'adverse'])
        st.metric("Alto Risco", high_risk)

    st.markdown("---")

    # Gráficos overview
    col1, col2 = st.columns(2)

    with col1:
        fig = plot_purity_distribution(conn)
        if fig:
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        fig = plot_mrd_status(conn)
        if fig:
            st.plotly_chart(fig, use_container_width=True)

    # Genes driver
    fig = plot_driver_genes(conn)
    if fig:
        st.plotly_chart(fig, use_container_width=True)

    # Tabela de amostras
    st.subheader("📊 Lista de Amostras")

    # Formatando colunas
    display_df = samples_df.copy()
    if 'purity' in display_df.columns:
        display_df['purity'] = (display_df['purity'] * 100).round(1).astype(str) + '%'
    if 'ploidy' in display_df.columns:
        display_df['ploidy'] = display_df['ploidy'].round(2)

    st.dataframe(
        display_df[['sample_id', 'analysis_date', 'purity', 'ploidy', 'mrd_status', 'leukemia_type', 'risk_category']],
        use_container_width=True,
        hide_index=True
    )

# ==================== DETALHES DA AMOSTRA ====================

elif view_mode == "🔍 Detalhes da Amostra":
    st.header("🔍 Detalhes da Amostra")

    # Seleção de amostra
    selected_sample = st.sidebar.selectbox(
        "Selecionar Amostra",
        samples_df['sample_id'].tolist()
    )

    # Carrega detalhes
    details = load_sample_details(conn, selected_sample)

    # Header da amostra
    sample_info = samples_df[samples_df['sample_id'] == selected_sample].iloc[0]

    col1, col2, col3 = st.columns(3)
    with col1:
        st.info(f"**Amostra:** {selected_sample}")
    with col2:
        if pd.notna(sample_info.get('leukemia_type')):
            st.info(f"**Tipo:** {sample_info['leukemia_type']}")
    with col3:
        if pd.notna(sample_info.get('risk_category')):
            risk_colors = {'favorable': '🟢', 'intermediate': '🟡', 'adverse': '🔴'}
            icon = risk_colors.get(sample_info['risk_category'], '⚪')
            st.info(f"**Risco:** {icon} {sample_info['risk_category']}")

    st.markdown("---")

    # Tabs
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📊 Overview", "🔬 PureCN", "📉 MRD", "🎯 Drivers", "📈 CNVs", "🔗 Fusões"
    ])

    # TAB 1: Overview
    with tab1:
        col1, col2, col3, col4 = st.columns(4)

        if not details['purecn'].empty:
            with col1:
                purity = details['purecn'].iloc[0]['purity'] * 100
                st.metric("Pureza Tumoral", f"{purity:.1f}%")
            with col2:
                ploidy = details['purecn'].iloc[0]['ploidy']
                st.metric("Ploidy", f"{ploidy:.2f}")

        if not details['mrd'].empty:
            with col3:
                mrd_status = details['mrd'].iloc[0]['mrd_status']
                st.metric("MRD Status", mrd_status)

        with col4:
            n_drivers = len(details['drivers'])
            st.metric("Driver Mutations", n_drivers)

        # Classificação
        if not details['classification'].empty:
            st.subheader("🎯 Classificação")
            class_info = details['classification'].iloc[0]

            col1, col2 = st.columns(2)
            with col1:
                st.success(f"**Tipo de Leucemia:** {class_info['leukemia_type']}")
                st.info(f"**Subtipo:** {class_info['subtype']}")
            with col2:
                st.warning(f"**Risco:** {class_info['risk_category']}")
                st.info(f"**Confiança:** {class_info['confidence']:.2%}")

            st.markdown(f"**Terapia Recomendada:** {class_info['recommended_therapy']}")

    # TAB 2: PureCN
    with tab2:
        if not details['purecn'].empty:
            purecn = details['purecn'].iloc[0]

            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Pureza", f"{purecn['purity']*100:.1f}%")
                st.metric("Sexo", purecn['sex'])
            with col2:
                st.metric("Ploidy", f"{purecn['ploidy']:.2f}")
                st.metric("Contaminação", f"{purecn['contamination']:.4f}")
            with col3:
                st.metric("Amplificações", int(purecn['n_amplifications']))
                st.metric("Deleções", int(purecn['n_deletions']))

            st.metric("LOH Events", int(purecn['n_loh']))
        else:
            st.info("Dados PureCN não disponíveis")

    # TAB 3: MRD
    with tab3:
        if not details['mrd'].empty:
            mrd = details['mrd'].iloc[0]

            status_colors = {
                'NEGATIVE': 'success',
                'LOW': 'info',
                'MODERATE': 'warning',
                'HIGH': 'error'
            }
            status_type = status_colors.get(mrd['mrd_status'], 'info')

            getattr(st, status_type)(f"**Status MRD:** {mrd['mrd_status']}")

            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Variantes", int(mrd['total_variants']))
            with col2:
                st.metric("Variantes MRD", int(mrd['mrd_variants']))
            with col3:
                st.metric("VAF Máxima", f"{mrd['max_vaf']:.3f}%")

            # Tabela de variantes MRD
            if not details['mrd_variants'].empty:
                st.subheader("Variantes MRD Detectadas")
                st.dataframe(
                    details['mrd_variants'][['gene', 'chromosome', 'position', 'ref', 'alt', 'vaf']],
                    use_container_width=True,
                    hide_index=True
                )
        else:
            st.info("Dados MRD não disponíveis")

    # TAB 4: Drivers
    with tab4:
        if not details['drivers'].empty:
            st.subheader(f"{len(details['drivers'])} Driver Mutations Detectadas")

            # Gráfico VAF
            fig = plot_vaf_distribution(conn, selected_sample)
            if fig:
                st.plotly_chart(fig, use_container_width=True)

            # Tabela
            st.dataframe(
                details['drivers'][['gene', 'chromosome', 'position', 'ref', 'alt', 'impact', 'exonic_func', 'aa_change']],
                use_container_width=True,
                hide_index=True
            )
        else:
            st.success("✅ Nenhuma driver mutation detectada")

    # TAB 5: CNVs
    with tab5:
        if not details['cnvs'].empty:
            st.subheader(f"{len(details['cnvs'])} CNVs Detectados")

            # Gráfico de distribuição
            fig = plot_cnv_distribution(conn, selected_sample)
            if fig:
                st.plotly_chart(fig, use_container_width=True)

            # Tabela
            st.dataframe(
                details['cnvs'][['chromosome', 'start', 'end', 'gene', 'log2_ratio', 'cn_status', 'validated']],
                use_container_width=True,
                hide_index=True
            )
        else:
            st.info("Nenhum CNV significativo detectado")

    # TAB 6: Fusões
    with tab6:
        if not details['fusions'].empty:
            st.subheader(f"{len(details['fusions'])} Fusões Detectadas")

            for _, fusion in details['fusions'].iterrows():
                with st.expander(f"🔗 {fusion['gene1']}-{fusion['gene2']}", expanded=True):
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.write(f"**Confiança:** {fusion['confidence']}")
                    with col2:
                        st.write(f"**Prognóstico:** {fusion['prognosis']}")
                    with col3:
                        pathog = "⭐ Sim" if fusion['pathognomonic'] else "Não"
                        st.write(f"**Patognomônica:** {pathog}")

                    if pd.notna(fusion['therapy']):
                        st.info(f"💊 **Terapia:** {fusion['therapy']}")
        else:
            st.info("Nenhuma fusão detectada")

# ==================== ANÁLISES COMPARATIVAS ====================

elif view_mode == "📊 Análises Comparativas":
    st.header("📊 Análises Comparativas")

    st.subheader("Comparação Entre Amostras")

    # Seleção de amostras para comparar
    selected_samples = st.multiselect(
        "Selecionar Amostras para Comparar",
        samples_df['sample_id'].tolist(),
        default=samples_df['sample_id'].tolist()[:5]
    )

    if len(selected_samples) < 2:
        st.warning("Selecione pelo menos 2 amostras para comparação")
    else:
        # Comparação de pureza/ploidy
        comp_data = samples_df[samples_df['sample_id'].isin(selected_samples)]

        col1, col2 = st.columns(2)

        with col1:
            if 'purity' in comp_data.columns:
                fig = px.bar(comp_data, x='sample_id', y='purity',
                           labels={'purity': 'Pureza Tumoral', 'sample_id': 'Amostra'},
                           title='Comparação de Pureza Tumoral')
                fig.update_traces(marker_color='#667eea')
                st.plotly_chart(fig, use_container_width=True)

        with col2:
            if 'ploidy' in comp_data.columns:
                fig = px.bar(comp_data, x='sample_id', y='ploidy',
                           labels={'ploidy': 'Ploidy', 'sample_id': 'Amostra'},
                           title='Comparação de Ploidy')
                fig.update_traces(marker_color='#764ba2')
                st.plotly_chart(fig, use_container_width=True)

        # Comparação de drivers
        st.subheader("Genes Driver por Amostra")
        driver_comp = pd.read_sql_query(f"""
            SELECT sample_id, gene, COUNT(*) as count
            FROM driver_mutations
            WHERE sample_id IN ({','.join(['?']*len(selected_samples))})
            GROUP BY sample_id, gene
        """, conn, params=selected_samples)

        if not driver_comp.empty:
            # Pivot para heatmap
            pivot_df = driver_comp.pivot(index='gene', columns='sample_id', values='count').fillna(0)

            fig = px.imshow(pivot_df,
                          labels=dict(x="Amostra", y="Gene", color="Mutações"),
                          title="Heatmap de Driver Mutations",
                          color_continuous_scale="Reds")
            st.plotly_chart(fig, use_container_width=True)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p><strong>VKaryo Dashboard</strong> - Virtual Karyotyping for Leukemia Analysis</p>
    <p>Powered by Streamlit | Database: SQLite</p>
</div>
""", unsafe_allow_html=True)
