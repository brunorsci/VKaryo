#!/usr/bin/env python3
"""
Consolida os resultados de UM programa (Arriba/Manta/DELLY/CNVkit/Mutect2)
de todas as amostras num único TSV, filtrado por um painel de genes de
interesse. Chamado por uma regra de verdade do workflow (não script solto),
para funcionar rodando o workflow sozinho em outro projeto.

Portabilidade: só `--sample-mapping` e `--driver-genes` são obrigatórios
(ambos genéricos, já usados em outras etapas do pipeline). `--extra-gene-panel`
e `--ground-truth` são específicos de uma coorte com verdade conhecida e
ficam OPCIONAIS -- sem eles, o filtro de genes usa só driver_genes_leukemia.txt
e a coluna "entity" fica vazia (não há como inferir diagnóstico esperado sem
uma referência de verdade conhecida, então isso não é gerado automaticamente
para outras coortes).

NOTA sobre a coluna 'gene' do CNVkit: segmentos de CNV grandes podem
sobrepor centenas de genes (lista separada por vírgula) -- isso e a coluna
'matched_genes' (filtrada pelo painel) ficam por ÚLTIMO em cada tabela, pra
não esconder as colunas importantes (log2/cn/posição/etc) à direita.
"""
import argparse
import csv
import glob
import os
import re
import sys

# PAUSADO (2026-10-20, decisão do usuário): uma classificação em camada de
# confiança (ex.: "confiável"/"razoável") a partir da tag de FILTER foi
# desenhada e chegou a ser testada, mas o usuário pediu pra NÃO escrever
# nenhum rótulo inventado num arquivo de resultado -- nem em português, nem
# sem critério formalmente validado. Por enquanto a saída expõe só a tag de
# FILTER NATIVA da ferramenta (vocabulário original, em inglês, sem
# reclassificação), pra o bioinformata revisar diretamente. A decisão de como
# (e se) resumir "passou ou não nos critérios de filtragem" fica pra depois,
# com critério formal definido -- ver PLANO_correcao_artigo_e_banco_dados.md
# seção 6 pro rascunho que tínhamos desenhado, caso sirva de ponto de partida.

# Func.refGene (ANNOVAR) mantidos na consolidação de SNV -- filtro por REGIÃO
# FUNCIONAL, não por nome de gene. Troca o antigo filtro por painel curado de
# genes (que escondia SOCS1 -- não está em driver_genes_leukemia.txt nem no
# painel extra, apesar de documentada): assim
# nenhuma variante some só por não estar numa lista de genes pré-definida.
# Deixa de fora intronic/intergenic/upstream/downstream/UTR (ruído de
# profundidade muito alto em WES tumor-only sem normal pareado -- ver
# PLANO_correcao_artigo_e_banco_dados.md).
FUNC_REGIAO_MANTIDA = {
    'exonic', 'splicing', 'exonic;splicing',
    'ncRNA_exonic', 'ncRNA_splicing', 'ncRNA_exonic;splicing',
}


def load_gene_filter(driver_genes_file, extra_panel_file=None):
    """Retorna (GENES: set de nomes em maiúsculo, GENE_COORDS: lista de
    (chrom, start, end, name) pra filtro por sobreposição de coordenada)."""
    driver_genes = set()
    gene_coords = []
    with open(driver_genes_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            name = parts[0].strip().upper()
            driver_genes.add(name)
            if len(parts) >= 4:
                try:
                    gene_coords.append((parts[1], int(parts[2]), int(parts[3]), name))
                except ValueError:
                    pass

    genes = set(driver_genes)
    if extra_panel_file and os.path.exists(extra_panel_file):
        with open(extra_panel_file) as f:
            for line in f:
                if line.startswith('>'):
                    name = line[1:].split(',')[0]
                    name = re.sub(r'_(NM|NR)_.*', '', name)
                    genes.add(name.strip().upper())

    return genes, gene_coords


def load_sample_map(sample_mapping_file, ground_truth_file=None):
    """sid (dna ou rna) -> (cell_line, entity, dna_id, rna_id). entity fica
    vazio a menos que --ground-truth seja dado (TSV com:
    RNA_sample, DNA_sample, Cell_line, Entity, ...)."""
    sample_map = {}
    with open(sample_mapping_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            dna, rna, cell = row['DNA_sample'], row['RNA_sample'], row.get('Cell_line', '')
            sample_map[dna] = [cell, '', dna, rna]
            sample_map[rna] = [cell, '', dna, rna]

    if ground_truth_file and os.path.exists(ground_truth_file):
        with open(ground_truth_file, newline='') as f:
            r = csv.reader(f, delimiter='\t')
            next(r, None)
            for row in r:
                if len(row) < 4:
                    continue
                rna, dna, cell, entity = row[0], row[1], row[2], row[3]
                for sid in (dna, rna):
                    if sid in sample_map:
                        sample_map[sid][1] = entity
                        if cell:
                            sample_map[sid][0] = cell

    return {k: tuple(v) for k, v in sample_map.items()}


def make_gene_hit(genes):
    def gene_hit(genestr):
        if not genestr or genestr == '-':
            return None
        toks = re.split(r'[;,/]', genestr)
        hits = [t.strip() for t in toks if t.strip().upper() in genes]
        return ';'.join(hits) if hits else None
    return gene_hit


def make_overlap_genes(gene_coords):
    def overlap_genes(chrom, start, end):
        return [name for c, s, e, name in gene_coords if c == chrom and start <= e and end >= s]
    return overlap_genes


def consolidate_arriba(output_dir, sample_map, gene_hit, out_path):
    """Fusões do Arriba.

    Filtro de painel de genes (gene_hit) deixou de ser critério de EXCLUSÃO
    -- vira coluna informativa (painel_driver_genes), pelo mesmo motivo já
    aplicado ao Mutect2 e ao AnnotSV: 14 das 100 amostras tinham de 8 a 82
    fusões no *_fusions.tsv bruto e ZERO linhas na tabela consolidada, só
    porque nenhum dos dois parceiros batia com a lista de genes curada
    (entre elas DAUDI, NU-DHL-1 e RI-1, justamente casos discutidos no
    artigo). A confiança atribuída pelo próprio Arriba já vem na coluna
    nativa `confidence` -- é ela que o analista deve usar na etapa final de
    filtragem, não uma exclusão silenciosa aqui."""
    n = 0
    with open(out_path, 'w', newline='') as out:
        w = None
        for path in sorted(glob.glob(os.path.join(output_dir, 'arriba', '*_fusions.tsv'))):
            rna_id = os.path.basename(path).replace('_fusions.tsv', '')
            cell, entity, dna_id, _ = sample_map.get(rna_id, ('', '', '', rna_id))
            with open(path) as f:
                header = f.readline().lstrip('#').rstrip('\n').split('\t')
                if w is None:
                    w = csv.writer(out, delimiter='\t', lineterminator='\n')
                    w.writerow(['sample_rna_id', 'sample_dna_id', 'cell_line', 'entity']
                               + header + ['painel_driver_genes'])
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 2:
                        continue
                    em_painel = 'sim' if (gene_hit(fields[0]) or gene_hit(fields[1])) else 'nao'
                    w.writerow([rna_id, dna_id, cell, entity] + fields + [em_painel])
                    n += 1
    return n


def consolidate_annotsv(output_dir, sample_map, gene_hit, program, out_path):
    """Manta ou DELLY, ambos anotados via AnnotSV com o mesmo formato de
    colunas (Gene_name, Annotation_mode, FILTER).

    Filtro de painel de genes (gene_hit) deixou de ser critério de EXCLUSÃO
    -- vira coluna informativa (painel_driver_genes), igual já feito pro
    Mutect2. Motivo: K-562/DELLY tinha 195 eventos PASS no VCF bruto e ZERO
    linhas na tabela consolidada, só porque nenhum batia com a lista de
    genes curada -- ver PLANO_correcao_artigo_e_banco_dados.md. Mantém o
    filtro Annotation_mode=='full' (não é filtro de qualidade, é só a
    granularidade de anotação do próprio AnnotSV -- evita linha duplicada
    'split' por breakpoint).

    A coluna FILTER (já presente no header nativo do AnnotSV) é mantida
    como veio da ferramenta -- sem reclassificação em camada de confiança
    (pausado, ver nota no topo do arquivo)."""
    n = 0
    pattern = os.path.join(output_dir, 'annotsv', program, f'*_annotsv_{program}.tsv')
    with open(out_path, 'w', newline='') as out:
        w = None
        for path in sorted(glob.glob(pattern)):
            dna_id = os.path.basename(path).replace(f'_annotsv_{program}.tsv', '')
            cell, entity, _, rna_id = sample_map.get(dna_id, ('', '', dna_id, ''))
            with open(path) as f:
                header = f.readline().rstrip('\n').split('\t')
                try:
                    gi = header.index('Gene_name')
                    ai = header.index('Annotation_mode')
                except ValueError:
                    continue
                if w is None:
                    w = csv.writer(out, delimiter='\t', lineterminator='\n')
                    w.writerow(['sample_dna_id', 'sample_rna_id', 'cell_line', 'entity']
                                + header + ['painel_driver_genes'])
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) <= gi or fields[ai] != 'full':
                        continue
                    em_painel = 'sim' if gene_hit(fields[gi]) else 'nao'
                    w.writerow([dna_id, rna_id, cell, entity] + fields + [em_painel])
                    n += 1
    return n


def consolidate_cnvkit(output_dir, sample_map, overlap_genes, out_path, source='call'):
    """Coluna 'gene' (CNVkit, pode ter centenas de genes por segmento) e
    'matched_genes' (filtrada pelo painel) ficam por último de propósito.

    source='call' (padrão, .call.cns) -- segmentação processada, funde
    regiões adjacentes de mesmo número de cópias; mais limpa, mas dilui
    picos focais (ver PLANO_correcao_artigo_e_banco_dados.md e
    evidencias_artigo_vkaryo/correcoes_manuscrito_mapeadas.md, item O2 --
    as amplificações e deleções de alto nível vêm do .cns bruto, não
    desta tabela). source='raw' (.cns) -- preserva resolução focal, usar
    pra material suplementar que sustenta números citados individualmente
    no texto."""
    # NOTA: o painel de genes deixou de ser critario de EXCLUSAO aqui, pelo
    # mesmo motivo dos caminhos do Mutect2, do AnnotSV e do Arriba. O
    # `if not genes: continue` que existia antes descartava todo segmento
    # que nao tocasse gene do painel curado -- 11.485 de 15.073 segmentos
    # brutos (80%) e 5.830 de 8.915 processados (70%) nunca chegavam a
    # tabela. Isso tornava impossivel conferir qualquer alteracao de numero
    # de copias fora do painel, e a descricao do material suplementar
    # ("segmentacao bruta") passava a impressao de completude que o arquivo
    # nao tinha. Agora todos os segmentos sao escritos e a pertinencia ao
    # painel vira a coluna informativa `painel_driver_genes`.
    pattern = '*.call.cns' if source == 'call' else '*.cns'
    exclude_suffix = None if source == 'call' else '.call.cns'
    n = 0
    with open(out_path, 'w', newline='') as out:
        w = None
        paths = sorted(glob.glob(os.path.join(output_dir, 'cnvkit', pattern)))
        if exclude_suffix:
            paths = [p for p in paths if not p.endswith(exclude_suffix) and not p.endswith('.bintest.cns')]
        for path in paths:
            dna_id = os.path.basename(path).split('_markduplicates_recal')[0]
            cell, entity, _, rna_id = sample_map.get(dna_id, ('', '', dna_id, ''))
            with open(path) as f:
                header = f.readline().rstrip('\n').split('\t')
                gene_idx = header.index('gene') if 'gene' in header else None
                if w is None:
                    w = csv.writer(out, delimiter='\t', lineterminator='\n')
                    other_header = [h for i, h in enumerate(header) if i != gene_idx]
                    w.writerow(['sample_dna_id', 'sample_rna_id', 'cell_line', 'entity']
                                + other_header
                                + ['painel_driver_genes', 'matched_genes', 'gene'])
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 3:
                        continue
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    genes = overlap_genes(chrom, start, end)
                    gene_col_value = fields[gene_idx] if gene_idx is not None else '-'
                    other_fields = [f for i, f in enumerate(fields) if i != gene_idx]
                    w.writerow([dna_id, rna_id, cell, entity] + other_fields
                                + ['sim' if genes else 'nao',
                                   ';'.join(genes) if genes else '-',
                                   gene_col_value])
                    n += 1
    return n


def consolidate_mutect2(output_dir, sample_map, gene_hit, out_path):
    keep_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                 'ExonicFunc.refGene', 'AAChange.refGene', 'CLNSIG', 'cosmic101',
                 'gnomad41_genome_AF', 'avsnp154']
    n = 0
    # Padrão ESPECÍFICO '_mutect2.hg38_multianno.txt': o diretório annovar/
    # também contém '{sample}_funcregion.hg38_multianno.txt' (saída do
    # estágio 1 da anotação em 2 estágios, rule mutect2_funcregion_scan),
    # que cobre o genoma inteiro e serve só pra decidir quais posições são
    # exônicas. Um glob '*.hg38_multianno.txt' pega os dois e duplica tudo.
    pattern = os.path.join(output_dir, 'mutect2', 'annovar', '*_mutect2.hg38_multianno.txt')
    with open(out_path, 'w', newline='') as out:
        w = None
        for path in sorted(glob.glob(pattern)):
            dna_id = os.path.basename(path).split('_mutect2')[0]
            cell, entity, _, rna_id = sample_map.get(dna_id, ('', '', dna_id, ''))
            with open(path) as f:
                header = f.readline().rstrip('\n').split('\t')
                idx = {h: i for i, h in enumerate(header)}
                keep_idx = [idx[c] for c in keep_cols if c in idx]
                keep_names = [c for c in keep_cols if c in idx]
                gi = idx.get('Gene.refGene')
                fi = idx.get('Func.refGene')
                filt_i = idx.get('Otherinfo10')
                fmt_i = idx.get('Otherinfo12')
                samp_i = idx.get('Otherinfo13')
                if gi is None or fi is None:
                    continue
                if w is None:
                    w = csv.writer(out, delimiter='\t', lineterminator='\n')
                    w.writerow(['sample_dna_id', 'sample_rna_id', 'cell_line', 'entity']
                                + keep_names + ['FILTER', 'painel_driver_genes', 'AF'])
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) <= gi or len(fields) <= fi:
                        continue
                    # Região funcional (exônica/splicing) -- não nome de gene, ver
                    # FUNC_REGIAO_MANTIDA acima. Substitui o antigo gene_hit() como
                    # filtro de INCLUSÃO; o painel de genes agora é só uma coluna
                    # informativa (painel_driver_genes), nunca exclui a linha.
                    if fields[fi] not in FUNC_REGIAO_MANTIDA:
                        continue
                    # FILTER fica exatamente como o Mutect2 gerou (ex.: 'PASS',
                    # 'clustered_events;haplotype') -- sem reclassificação em
                    # camada de confiança (pausado, ver nota no topo do arquivo).
                    filter_value = fields[filt_i] if filt_i is not None and len(fields) > filt_i else ''
                    em_painel = 'sim' if gene_hit(fields[gi]) else 'nao'
                    af = ''
                    if fmt_i is not None and samp_i is not None and len(fields) > samp_i:
                        fk = fields[fmt_i].split(':')
                        sv = fields[samp_i].split(':')
                        if 'AF' in fk and len(sv) == len(fk):
                            af = sv[fk.index('AF')]
                    row = ([dna_id, rna_id, cell, entity] + [fields[i] for i in keep_idx]
                           + [filter_value, em_painel, af])
                    w.writerow(row)
                    n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--program', required=True, choices=['arriba', 'manta', 'delly', 'cnvkit', 'mutect2'])
    ap.add_argument('--output-dir', required=True, help='output_vkaryo/ do projeto')
    ap.add_argument('--sample-mapping', required=True, help='sample_mapping.tsv (DNA_sample/RNA_sample/Cell_line)')
    ap.add_argument('--driver-genes', required=True, help='driver_genes_leukemia.txt (ou equivalente)')
    ap.add_argument('--extra-gene-panel', default=None, help='opcional: painel extra estilo FASTA-header (ex. painel curado da coorte)')
    ap.add_argument('--ground-truth', default=None, help='opcional: tabela de verdade conhecida (entity)')
    ap.add_argument('--cnvkit-source', choices=['call', 'raw'], default='call',
                     help='só --program cnvkit: "call" (.call.cns, padrão) ou "raw" (.cns, preserva picos focais)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    genes, gene_coords = load_gene_filter(args.driver_genes, args.extra_gene_panel)
    sample_map = load_sample_map(args.sample_mapping, args.ground_truth)
    gene_hit = make_gene_hit(genes)
    overlap_genes = make_overlap_genes(gene_coords)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.program == 'arriba':
        n = consolidate_arriba(args.output_dir, sample_map, gene_hit, args.out)
    elif args.program in ('manta', 'delly'):
        n = consolidate_annotsv(args.output_dir, sample_map, gene_hit, args.program, args.out)
    elif args.program == 'cnvkit':
        n = consolidate_cnvkit(args.output_dir, sample_map, overlap_genes, args.out, source=args.cnvkit_source)
    elif args.program == 'mutect2':
        n = consolidate_mutect2(args.output_dir, sample_map, gene_hit, args.out)
    else:
        print(f'Programa desconhecido: {args.program}', file=sys.stderr)
        sys.exit(1)

    # O painel NAO e mais criterio de exclusao -- vira a coluna
    # painel_driver_genes. A mensagem antiga dizia "linhas filtradas", o que
    # era enganoso depois dessa mudanca.
    print(f'{os.path.basename(args.out)}: {n} linhas escritas '
          f'(coluna painel_driver_genes marcada contra {len(genes)} genes)', file=sys.stderr)


if __name__ == '__main__':
    main()
