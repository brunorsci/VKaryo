import os
import webbrowser

# Diretório onde os arquivos HTML estão localizados
input_directory = '/mnt/c/Users/bruno/Desktop/genefuse/'

# Caminho para o arquivo HTML mesclado
output_file = '/mnt/c/Users/bruno/Desktop/genefuse/merged_output.html'

# Função para ler o conteúdo de um arquivo HTML
def read_html(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.read()

# Obtém todos os arquivos HTML no diretório
html_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.html')]

# Inicializa o conteúdo mesclado
merged_content = ""

# Mescla o conteúdo dos arquivos HTML
for file in html_files:
    file_content = read_html(file)
    # Adiciona uma quebra de linha entre os arquivos para separá-los
    merged_content += file_content + "\n\n"

# Salva o conteúdo mesclado em um novo arquivo HTML
with open(output_file, 'w', encoding='utf-8') as file:
    file.write(merged_content)

print(f"Arquivo HTML mesclado salvo em: {output_file}")

# Abre o arquivo HTML no navegador padrão
webbrowser.open(f'file://{os.path.abspath(output_file)}')
