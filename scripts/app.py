import os
import json
import dash
from dash import dcc, html
import plotly.graph_objects as go

# Inicializa o app Dash
app = dash.Dash(__name__)

# Função para criar o gráfico de sequências nucleotídicas
def create_sequence_plot(fusion_name, fusion_data):
    fig = go.Figure()

    # Adiciona sequência dos genes envolvidos
    left_gene_seq = fusion_data['left']['reference']
    right_gene_seq = fusion_data['right']['reference']
    
    fig.add_trace(go.Scatter(
        x=list(range(len(left_gene_seq))),
        y=[1] * len(left_gene_seq),
        mode='lines+markers',
        name=f"Gene Esquerdo: {fusion_data['left']['gene_name']}",
        text=left_gene_seq,
        textposition='top center'
    ))

    fig.add_trace(go.Scatter(
        x=list(range(len(right_gene_seq))),
        y=[2] * len(right_gene_seq),
        mode='lines+markers',
        name=f"Gene Direito: {fusion_data['right']['gene_name']}",
        text=right_gene_seq,
        textposition='top center'
    ))

    fig.update_layout(
        title=f"Sequências Nucleotídicas da Fusão: {fusion_name}",
        xaxis_title="Posição",
        yaxis_title="",
        showlegend=True,
        height=300,
        margin=dict(l=20, r=20, t=40, b=20)
    )

    return fig

# Função para gerar layout HTML
def create_layout(fusions):
    layout = []
    for fusion_name, fusion_data in fusions.items():
        layout.append(html.Div([
            html.H2(f"Fusão Genética: {fusion_name}"),
            dcc.Graph(figure=create_sequence_plot(fusion_name, fusion_data))
        ], style={'margin-bottom': '20px'}))
    return html.Div(layout, style={'padding': '20px'})

# Carrega os dados
def load_json_data(directory):
    fusions = {}
    if not os.path.exists(directory):
        raise FileNotFoundError(f"O diretório '{directory}' não foi encontrado.")
    
    for filename in os.listdir(directory):
        if filename.endswith("_fusions.json"):
            with open(os.path.join(directory, filename), 'r') as file:
                data = json.load(file)
                fusion_name = list(data.keys())[0]  # Assume que a fusão é a primeira chave
                fusions[fusion_name] = data[fusion_name]
    return fusions

# Atualize o caminho do diretório com o caminho correto
data_directory = '/mnt/c/Users/bruno/Desktop/genefuse/resultados/'  # Atualize com o caminho correto

try:
    fusions = load_json_data(data_directory)
except FileNotFoundError as e:
    print(e)
    fusions = {}

# Cria o layout com os dados carregados
app.layout = create_layout(fusions)

if __name__ == '__main__':
    app.run_server(debug=True)
