import os
from bs4 import BeautifulSoup

def extract_fusion_info(html_path):
    """Extrai informações sobre o número de fusões e os nomes dos genes dos arquivos HTML."""
    with open(html_path, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')

    # Número de fusões encontradas
    fusion_info = soup.find('div', id='menu')
    num_fusions = fusion_info.find('p').text.strip().split()[1]

    # Nomes dos genes
    protein_table = soup.find('table', class_='protein_table')
    genes = [td.text.strip() for td in protein_table.find_all('td')]
    exons_info = soup.find('div', class_='fusion_head').text.strip()
    exons_info = exons_info.replace(',', '.')

    return num_fusions, genes, exons_info

def generate_index(directory):
    """Gera o arquivo index.html com links para cada HTML e informações adicionais."""
    html_files = [f for f in os.listdir(directory) if f.endswith('_genefuse_report.html')]
    html_files.sort()  # Ordena os arquivos por nome

    index_content = '''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Fusões Genéticas</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                margin: 0;
                background-color: #f4f4f4;
                color: #333;
            }
            header {
                background-color: #343a40;
                color: #fff;
                padding: 20px;
                text-align: center;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }
            header h1 {
                margin: 0;
                font-size: 2.5em;
                font-weight: 700;
            }
            .container {
                width: 80%;
                margin: 20px auto;
            }
            .card {
                background: #fff;
                border-radius: 8px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                margin-bottom: 20px;
                padding: 20px;
                overflow: hidden;
                display: flex;
                flex-direction: column;
                align-items: flex-start;
            }
            .card h2 {
                margin: 0;
                font-size: 1.5em;
                color: #007bff;
                font-weight: 600;
            }
            .card p {
                margin: 5px 0;
                font-size: 1em;
                color: #555;
            }
            .card a {
                text-decoration: none;
                color: #007bff;
            }
            .card a:hover {
                text-decoration: underline;
            }
            .gene-info {
                margin-top: 10px;
                padding: 10px;
                background: #e9ecef;
                border-radius: 5px;
                font-size: 0.9em;
            }
            .gene-info strong {
                display: block;
                margin-bottom: 5px;
            }
            footer {
                background-color: #343a40;
                color: #fff;
                padding: 10px;
                text-align: center;
                position: fixed;
                width: 100%;
                bottom: 0;
                box-shadow: 0 -2px 6px rgba(0,0,0,0.1);
            }
            footer p {
                margin: 0;
            }
        </style>
    </head>
    <body>
        <header>
            <h1>VKaryo: genefuse</h1>
            <h2>Fusões Gênicas</h2>

        </header>
        <div class="container">
    '''

    for file_name in html_files:
        file_path = os.path.join(directory, file_name)
        num_fusions, genes, exons_info = extract_fusion_info(file_path)
        base_name = file_name.replace('_genefuse_report.html', '')
        index_content += f'''
        <div class="card">
            <h2><a href="{file_name}">{base_name}</a></h2>
            <p>Fusões encontradas: {num_fusions}</p>
            <div class="gene-info">
                <strong>Genes envolvidos:</strong><br>
                {genes[0]}<br>
                {genes[1]}<br><br>
                <strong>Exons:</strong><br>
                {exons_info}<br>
            </div>
        </div><br><br>
        '''

    index_content += '''
        </div>
        <footer>
            <p>&copy; 2024 VKaryo: genefuse</p>
        </footer>
    </body>
    </html>
    '''

    index_path = os.path.join(directory, 'index.html')
    with open(index_path, 'w') as file:
        file.write(index_content)

    print(f'Arquivo {index_path} criado com sucesso.')

# Use o diretório onde os arquivos HTML estão localizados
generate_index('resultados')