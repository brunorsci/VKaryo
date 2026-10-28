#!/usr/bin/env python3
"""
Estratificação automática de leucemias/linfomas a partir dos achados do VKaryo.

Arquitetura de 3 camadas:

  Camada 1 - Risco ELN 2022 (LMA): classificação categórica exatamente como
             publicada (Döhner et al., Blood 2022) -- não é mais uma soma de
             pontos aproximada. 100% baseada em achados genéticos/citogenéticos.

  Camada 2 - Entidade WHO-HAEM5 / ICC 2022: atribuída SOMENTE quando uma
             alteração genética definidora de entidade é detectada (fusão,
             mutação ou translocação reconhecida pela classificação). Quando
             nada bate, o relatório diz isso explicitamente em vez de
             adivinhar -- a classificação por diferenciação morfológica
             exige dados que o VKaryo não tem (morfologia/imunofenotipagem).

  Camada 3 - Achados adicionais (camada própria): qualquer variante em gene
             driver recorrente de leucemia (driver_genes_leukemia.txt) que
             não foi consumida pelas camadas 1-2, reportada com uma nota de
             significância clínica geral (não uma reclassificação formal).

Uso:
    python classify_leukemia.py \
        --fusions output_vkaryo/arriba_annotated/sample_fusions_annotated.tsv \
        --cnvs output_vkaryo/cnv_validated/sample_cnv_validated.tsv \
        --snvs output_vkaryo/mutect2/annovar/sample_mutect2.hg38_multianno.txt \
        --aneuploidy-summary output_vkaryo/plots/sample_aneuploidy_summary.txt \
        --sample SAMPLE_ID \
        --output output_vkaryo/classification/sample_classification.json \
        --report output_vkaryo/classification/sample_classification.md
"""

import pandas as pd
import json
import argparse
import logging
import re
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def _dedupe_protein_changes(aachange_refgene):
    """AAChange.refGene do ANNOVAR lista uma entrada POR TRANSCRITO (às vezes
    10+ para o mesmo gene, cada uma com sua própria numeração de aminoácido --
    ex. TP53 c.G743A é p.R248Q no transcrito canônico NM_000546, mas p.R116Q/
    p.R209Q/p.R89Q em transcritos alternativos do mesmo evento). O antigo
    corte ingênuo em 120 caracteres cortava a string ANTES de chegar no
    transcrito canônico, mostrando ao usuário a numeração "errada"/menos
    reconhecida de um hotspot bem conhecido (ex.: KASUMI-1, TP53 R248Q --
    hotspot clássico, aparecia como "R116Q").
    Em vez de truncar por tamanho, deduplica pela notação p.XXXN e mostra só
    as variações de proteína distintas (bem mais curto e não corta no meio)."""
    s = str(aachange_refgene or '').strip()
    if not s or s == '.' or s.lower() == 'nan':
        return '.'
    # Ordena por número do transcrito RefSeq (NM_XXXXXX) ascendente antes de
    # deduplicar -- heurística: números baixos (ex. NM_000546 do TP53) tendem
    # a ser a designação canônica/MANE Select mais antiga e reconhecida,
    # numeros altos (ex. NM_001126112) costumam ser isoformas alternativas
    # menos citadas na literatura. Sem isso, a ordem de aparição na string
    # do ANNOVAR é arbitrária e pode esconder a notação canônica.
    entries = []
    for entry in s.split(','):
        acc_match = re.search(r'NM_(\d+)', entry)
        acc_num = int(acc_match.group(1)) if acc_match else 10**9
        m = re.search(r'p\.[A-Za-z0-9*]+', entry)
        label = m.group(0) if m else entry.split(':')[-1]
        entries.append((acc_num, label))
    entries.sort(key=lambda x: x[0])
    seen = []
    for _, label in entries:
        if label not in seen:
            seen.append(label)
        if len(seen) >= 3:
            break
    return ','.join(seen) if seen else s[:80]


# ============================================================
# CAMADA 1 -- ELN 2022 (risco em LMA), Döhner et al. Blood 2022;140:1345-1377
# Regras categóricas, não aditivas: um único marcador "hard adverse" já
# define risco adverso independentemente de qualquer outro achado.
# ============================================================

ELN2022_FAVORABLE_FUSIONS = {'RUNX1--RUNX1T1', 'CBFB--MYH11'}
ELN2022_ADVERSE_HARD_FUSIONS = {
    'DEK--NUP214', 'BCR--ABL1',  # BCR-ABL1 em contexto de novo AML
}
# KMT2A-rearranged é adverso, EXCETO especificamente t(9;11)/KMT2A::MLLT3
ELN2022_KMT2A_FAVORABLE_EXCEPTION = {'KMT2A--MLLT3', 'MLLT3--KMT2A'}

ELN2022_ADVERSE_HARD_GENES = {'TP53'}
ELN2022_ADVERSE_SOFT_GENES = {
    'ASXL1', 'BCOR', 'EZH2', 'RUNX1', 'SF3B1', 'SRSF2', 'STAG2', 'U2AF1', 'ZRSR2'
}


def eln2022_risk(fusion_names, mutated_genes, cebpa_bzip_inframe, npm1_mut,
                  flt3_itd, complex_karyotype, monosomal_karyotype,
                  adverse_cytogenetics):
    """
    Retorna (grupo, fatores) segundo os critérios categóricos ELN 2022.
    fusion_names: set de nomes de fusão normalizados (GENE1--GENE2)
    mutated_genes: set de genes com variante driver detectada
    adverse_cytogenetics: set de descrições de citogenética adversa já
        identificadas por outras camadas (ex.: 'MECOM rearrangement',
        '-5/del(5q)', '-7', 'abn(17p)') -- opcional, vazio se não avaliado.
    """
    factors = []

    # --- ADVERSO (hard) -- vence qualquer outro achado ---
    hard_adverse = False
    if complex_karyotype:
        hard_adverse = True
        factors.append('Cariótipo complexo')
    if monosomal_karyotype:
        hard_adverse = True
        factors.append('Cariótipo monossômico')
    if 'TP53' in mutated_genes:
        hard_adverse = True
        factors.append('TP53 mutado')
    for f in fusion_names:
        # NOTA: 'KMT2A' in f / 'MECOM' in f sobre a string inteira
        # "GENE1--GENE2" dá falso positivo em qualquer gene que contenha
        # essas letras como prefixo (ex. MLLT10 não é KMT2A/MLL, mas contém
        # "MLL"; mesma classe de bug em annotate_arriba_leukemia.py -- ver
        # parceiro inesperado). Por isso os dois genes da fusão são
        # checados individualmente, por igualdade.
        genes_in_fusion = f.split('--') if '--' in f else [f]
        if f in ELN2022_ADVERSE_HARD_FUSIONS:
            hard_adverse = True
            factors.append(f'Fusão adversa: {f}')
        if ('KMT2A' in genes_in_fusion) and f not in ELN2022_KMT2A_FAVORABLE_EXCEPTION:
            hard_adverse = True
            factors.append(f'Rearranjo de KMT2A (exceto t(9;11)): {f}')
        if 'MECOM' in genes_in_fusion:
            hard_adverse = True
            factors.append(f'Rearranjo de MECOM: {f}')
    for c in adverse_cytogenetics:
        hard_adverse = True
        factors.append(f'Citogenética adversa: {c}')

    if hard_adverse:
        return 'Adverso', factors

    # --- FAVORÁVEL ---
    favorable = False
    for f in fusion_names:
        if f in ELN2022_FAVORABLE_FUSIONS:
            favorable = True
            factors.append(f'Fusão CBF favorável: {f}')
    if npm1_mut and not flt3_itd:
        favorable = True
        factors.append('NPM1 mutado sem FLT3-ITD')
    if cebpa_bzip_inframe:
        favorable = True
        factors.append('CEBPA com mutação bZIP em quadro de leitura')

    if favorable:
        # marcadores "soft adverse" perdem peso quando um achado favorável
        # está presente, por definição do ELN 2022
        return 'Favorável', factors

    # --- ADVERSO (soft) ---
    soft_hits = ELN2022_ADVERSE_SOFT_GENES & mutated_genes
    if soft_hits:
        factors.extend(f'Mutação adversa: {g}' for g in sorted(soft_hits))
        return 'Adverso', factors

    # --- INTERMEDIÁRIO ---
    if npm1_mut and flt3_itd:
        factors.append('NPM1 mutado com FLT3-ITD')
        return 'Intermediário', factors
    if flt3_itd:
        factors.append('FLT3-ITD sem NPM1 mutado (sem marcador adverso)')
        return 'Intermediário', factors
    for f in fusion_names:
        if f in ELN2022_KMT2A_FAVORABLE_EXCEPTION:
            factors.append(f't(9;11)/KMT2A::MLLT3: {f}')
            return 'Intermediário', factors

    factors.append('Sem alteração classificada como favorável ou adversa pelo ELN 2022')
    return 'Intermediário', factors


# ============================================================
# CAMADA 1b -- Risco genômico de LLC e IPSS-R de SMD.
#
# Usa a mesma lógica de detecção de deleção focal por coordenada (não
# contagem de segmentos do CNVkit) usada no resto do módulo, com o mesmo
# aviso de qualidade do cariótipo quando aplicável (ver
# `_detect_karyotype_complexity`).
#
# Assim como o ELN 2022 (Camada 1), calculadas para toda amostra com
# achados compatíveis, independentemente de o VKaryo "saber" que a amostra
# é LLC ou SMD -- a aplicabilidade clínica depende de confirmação
# morfológica/imunofenotípica, explicitada no relatório.
# ============================================================

# coordenadas hg38 aproximadas dos loci focais usados (além de TP53/ATM, já
# cobertos por driver_genes_leukemia.txt): 13q14 (região mínima deletada em
# LLC, ~DLEU2/MIR15A/MIR16-1), 5q31-33 (SMD, região próxima a EGR1) e 20q.
CLL_MDS_FOCAL_COORDS = {
    'del17p': ('chr17', 1000000, 22200000),      # braço curto inteiro (região de TP53, FISH clássico)
    'del11q': ('chr11', 60000000, 135000000),    # braço longo a partir do centrômero (região de ATM)
    'del13q14': ('chr13', 49000000, 51000000),   # região mínima deletada (DLEU2/MIR15A/MIR16-1)
    'del5q': ('chr5', 135000000, 150000000),
    'del20q': ('chr20', 30000000, 64000000),
}


def cll_risk_stratification(mutated_genes, focal_dels, karyotype_quality_flag):
    """
    Estratificação de risco genômico em LLC, baseada nas alterações
    citogenéticas/moleculares clássicas (Döhner et al., NEJM 2000;
    incorporação de TP53/NOTCH1/SF3B1/BIRC3 conforme guidelines iwCLL/ERIC
    mais recentes). Não substitui o status mutacional de IGHV (não
    determinável a partir de WES) nem estadiamento clínico (Rai/Binet).
    """
    features = {
        'del17p': 'del17p' in focal_dels, 'tp53_mut': 'TP53' in mutated_genes,
        'del11q': 'del11q' in focal_dels,
        'del13q': 'del13q14' in focal_dels,
        'trisomy12': 'trisomy12' in focal_dels,
        'notch1_mut': 'NOTCH1' in mutated_genes, 'sf3b1_mut': 'SF3B1' in mutated_genes,
        'atm_mut': 'ATM' in mutated_genes, 'birc3_mut': 'BIRC3' in mutated_genes,
    }

    factors = []
    if features['del17p'] or features['tp53_mut']:
        group = 'Risco muito alto'
        factors.append('del(17p) e/ou TP53 mutado')
    elif features['del11q'] or features['notch1_mut'] or features['birc3_mut']:
        group = 'Risco alto'
        if features['del11q']:
            factors.append('del(11q) (região de ATM)')
        if features['notch1_mut']:
            factors.append('NOTCH1 mutado')
        if features['birc3_mut']:
            factors.append('BIRC3 mutado')
    elif features['trisomy12'] or features['atm_mut']:
        group = 'Risco intermediário'
        if features['trisomy12']:
            factors.append('Trissomia do cromossomo 12')
        if features['atm_mut']:
            factors.append('ATM mutado (sem del(11q) detectada)')
    elif features['del13q']:
        group = 'Risco baixo'
        factors.append('del(13q) isolada')
    else:
        group = 'Indeterminado'
        factors.append('Nenhum marcador citogenético/molecular clássico de LLC detectado')

    if features['sf3b1_mut']:
        factors.append('SF3B1 mutado (adverso, independente do grupo citogenético)')

    note = ('Aplicável apenas se o diagnóstico for LLC/leucemia prolinfocítica -- não considera status '
            'mutacional de IGHV (não determinável a partir de WES) nem estadiamento clínico Rai/Binet.')
    if karyotype_quality_flag:
        note += f' ATENÇÃO: {karyotype_quality_flag}'
    return {'risk_group': group, 'factors': factors, 'features': features, 'applicable_note': note}


def ipss_r_lite(mutated_genes, focal_dels, chrom_calls, complex_karyotype, karyotype_quality_flag):
    """
    Aproximação do componente citogenético + mutacional do IPSS-R/IPSS-M
    para SMD (Greenberg et al., Blood 2012; Bernard et al., NEJM Evid 2022).
    Versão "lite": não inclui contagem de blastos nem citopenias (não
    disponíveis a partir de WES/RNA-seq), portanto reporta apenas o
    componente genômico, não o escore IPSS-R completo.
    """
    score = 0
    factors = []
    cytogenetic_group = 'Bom'

    has_del5q = 'del5q' in focal_dels
    has_del20q = 'del20q' in focal_dels
    # NOTA: monossomia do 7 vem de chrom_calls (nível de cromossomo inteiro)
    # -- precisa passar pela mesma blindagem de qualidade que protege
    # focal_dels/complex_karyotype, senão um artefato sistêmico do CNVkit
    # (muitos cromossomos "loss" ao mesmo tempo) é lido como monossomia real.
    has_monosomy7 = (not karyotype_quality_flag) and chrom_calls.get('chr7') == 'loss'

    if complex_karyotype:
        cytogenetic_group = 'Muito ruim'
        score += 4
        factors.append('Cariótipo complexo (≥3 alterações)')
    elif has_monosomy7:
        cytogenetic_group = 'Ruim'
        score += 3
        factors.append('Monossomia do cromossomo 7')
    elif has_del5q:
        cytogenetic_group = 'Bom'
        factors.append('del(5q)')
    elif has_del20q:
        cytogenetic_group = 'Bom'
        factors.append('del(20q)')

    adverse_mutations = {'TP53', 'ASXL1', 'RUNX1', 'EZH2'} & mutated_genes
    if 'TP53' in mutated_genes:
        score += 2
        factors.append('TP53 mutado')
    if len(adverse_mutations) > 1:
        score += 1
        factors.append(f"Múltiplas mutações adversas: {', '.join(sorted(adverse_mutations))}")

    if score <= 1.5:
        category = 'Muito baixo/Baixo (componente genômico)'
    elif score <= 3:
        category = 'Baixo/Intermediário (componente genômico)'
    elif score <= 4.5:
        category = 'Intermediário/Alto (componente genômico)'
    else:
        category = 'Alto/Muito alto (componente genômico)'

    note = ('Componente citogenético/mutacional apenas -- o IPSS-R/IPSS-M completo exige contagem de '
            'blastos medulares e citopenias periféricas, não disponíveis a partir de sequenciamento. '
            'Aplicável apenas se o diagnóstico for SMD.')
    if karyotype_quality_flag:
        note += f' ATENÇÃO: {karyotype_quality_flag}'
    return {'score_genomic_component': score, 'category': category, 'cytogenetic_group': cytogenetic_group,
            'factors': factors, 'applicable_note': note}


# ============================================================
# CAMADA 2 -- Entidade WHO-HAEM5 / ICC 2022, apenas quando uma alteração
# genética definidora de entidade é detectada. Fusões já vêm resolvidas de
# annotate_arriba_leukemia.py (annot_classification começando com "WHO:").
# Aqui cobrimos as entidades definidas por variantes de nucleotídeo único ou
# CNV que não passam pelo Arriba (não dependem de RNA-seq).
# ============================================================

WHO_SNV_DEFINING = {
    # gene -> (entidade WHO-HAEM5, nota ICC se divergente)
    'NPM1': ('LMA com mutação de NPM1', None),
    'CEBPA': ('LMA com mutação bZIP de CEBPA (em quadro de leitura)', None),
}

WHO_CNV_DEFINING = {
    # driver gene / região -> entidade
    'SOCS1': ('Deleção bialélica de SOCS1 (achado recorrente em PMBL/linfoma de Hodgkin, não formalmente '
              'definidor de entidade WHO-HAEM5, reportado como achado relevante)', None),
}


def assign_who_entity(fusion_rows, mutated_gene_details, cnv_hits):
    """
    Retorna dict {entity, source, icc_note, evidence} ou None se nenhuma
    alteração definidora de entidade foi encontrada.
    """
    # Prioridade 1: fusão com entidade WHO-HAEM5 reconhecida (já resolvida
    # por annotate_arriba_leukemia.py)
    defining_fusions = [f for f in fusion_rows if str(f.get('annot_classification', '')).startswith('WHO:')]
    if defining_fusions:
        f = defining_fusions[0]
        return {
            'entity': f.get('annot_leukemia_type', f.get('annot_classification')),
            'who_classification': f.get('annot_classification'),
            'source': f"Fusão: {f.get('fusion')}",
            'confidence': f.get('confidence', 'unknown'),
            'evidence': [f"{f.get('fusion')} ({f.get('confidence')}) -- {f.get('annot_classification')}"],
            'additional_fusions_same_sample': [d.get('fusion') for d in defining_fusions[1:]],
        }

    # Prioridade 2: mutação definidora (NPM1, CEBPA bZIP)
    for gene, (entity, icc) in WHO_SNV_DEFINING.items():
        if gene in mutated_gene_details:
            detail = mutated_gene_details[gene]
            return {
                'entity': entity,
                'who_classification': f'WHO: {entity}',
                'icc_note': icc,
                'source': f'SNV/indel: {gene}',
                'confidence': 'high' if detail.get('impact') == 'HIGH' else 'medium',
                'evidence': [f"{gene} {detail.get('mutation', '')}"],
            }

    return None


# ============================================================
# CAMADA 3 -- achados adicionais em genes driver (curadoria própria)
# ============================================================

GENE_CLINICAL_NOTES = {
    'FLT3': 'Alvo terapêutico (inibidores de FLT3); prognóstico depende de ITD vs TKD e razão alélica.',
    'IDH1': 'Alvo terapêutico (ivosidenibe).', 'IDH2': 'Alvo terapêutico (enasidenibe).',
    'DNMT3A': 'Mutação de idade avançada / CHIP; cooperante comum em LMA.',
    'TET2': 'Regulador epigenético; recorrente em neoplasias mieloides (LMA/SMD/LMMC).',
    'ASXL1': 'Marcador de risco adverso (ELN 2022); recorrente em SMD/LMMC.',
    'RUNX1': 'Marcador de risco adverso (ELN 2022) quando mutado (distinto da fusão RUNX1::RUNX1T1).',
    'KIT': 'Mutação frequentemente concomitante com RUNX1::RUNX1T1; impacto prognóstico contextual.',
    'JAK2': 'Driver de neoplasia mieloproliferativa (V617F é o hotspot clássico).',
    'CALR': 'Driver de neoplasia mieloproliferativa (trombocitemia essencial/mielofibrose).',
    'MPL': 'Driver de neoplasia mieloproliferativa.',
    'SF3B1': 'Mutação de splicing; associada a SMD com sideroblastos em anel.',
    'SRSF2': 'Mutação de splicing; recorrente em LMMC.',
    'U2AF1': 'Mutação de splicing; marcador de risco adverso (ELN 2022).',
    'ZRSR2': 'Mutação de splicing; marcador de risco adverso (ELN 2022).',
    'EZH2': 'Regulador epigenético (domínio SET); mutações de ganho de função recorrentes em linfoma folicular/DLBCL GC.',
    'NRAS': 'Via de sinalização RAS; driver comum em neoplasias mieloides.',
    'KRAS': 'Via de sinalização RAS; driver comum em neoplasias mieloides.',
    'PTPN11': 'Via de sinalização RAS; recorrente em LMMJ e AML.',
    'CBL': 'Via de sinalização RAS; recorrente em LMMC/LMMJ.',
    'PHF6': 'Recorrente em T-ALL e AML.',
    'BCOR': 'Marcador de risco adverso (ELN 2022); regulador transcricional.',
    'BCORL1': 'Recorrente em neoplasias mieloides.',
    'STAG2': 'Coesina; marcador de risco adverso (ELN 2022).',
    'IKZF1': 'Deleções/mutações associadas a B-ALL de alto risco (incl. Ph-like).',
    'PAX5': 'Recorrente em B-ALL; pode definir subtipo B-ALL com alterações de PAX5.',
    'NOTCH1': 'Driver mais comum em T-ALL; também recorrente em CLL de alto risco.',
    'ATM': 'Via de reparo de DNA; deleção/mutação associada a CLL de risco intermediário-alto (del11q).',
    'BIRC3': 'Recorrente em CLL, associado a resistência e pior prognóstico.',
    'BTK': 'Alvo terapêutico direto (ibrutinibe); mutações associadas a resistência adquirida.',
    'MYD88': 'Hotspot L265P clássico em macroglobulinemia de Waldenström/DLBCL ABC.',
    'PTEN': 'Supressor tumoral; perda recorrente em neoplasias linfoides T.',
    'CDKN2A': 'Supressor tumoral (p16); deleção recorrente em B-ALL e outras neoplasias.',
    'WT1': 'Recorrente em LMA, associado a pior resposta em alguns contextos.',
    'GATA2': 'Predisposição familiar a LMA/SMD quando germinativo; rearranjos somáticos (3q26) são adversos.',
}


def build_layer3_findings(mutated_gene_details, used_genes):
    findings = []
    for gene, detail in mutated_gene_details.items():
        if gene in used_genes:
            continue
        findings.append({
            'gene': gene,
            'mutation': detail.get('mutation', 'N/A'),
            'impact': detail.get('impact', 'MEDIUM'),
            'clinical_note': GENE_CLINICAL_NOTES.get(
                gene, 'Gene driver recorrente de leucemia (driver_genes_leukemia.txt); '
                      'sem nota curada específica -- revisar literatura/COSMIC.'),
        })
    return findings


# ============================================================
# Classe principal -- carregamento de dados e orquestração das 3 camadas
# (funções de ploidia / CLL / IPSS-R reaproveitadas da versão anterior)
# ============================================================

class LeukemiaClassifier:
    def __init__(self, sample_id, driver_genes_file=None):
        self.sample_id = sample_id
        self.fusion_rows = []
        self.mutated_gene_details = {}
        self.aneuploidy_events = []
        self.cnvs_df = None
        self.complex_karyotype = False
        self.monosomal_karyotype = False
        self.karyotype_quality_flag = None
        self.karyotype_arm_level_unreliable = False
        self.chrom_calls = {}
        self.focal_dels = set()
        # Painel autoritativo de genes-driver (driver_genes_leukemia.txt, 75
        # genes, WHO 2022) -- fonte unica de verdade para "e um gene driver?"
        # em vez da chave de GENE_CLINICAL_NOTES (que so cobre os genes com
        # nota curada, um subconjunto menor e nao a lista completa do painel).
        self.driver_genes = set()
        if driver_genes_file:
            with open(driver_genes_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    self.driver_genes.add(line.split('\t')[0])
        if not self.driver_genes:
            self.driver_genes = self._driver_gene_names()

    # -------- loaders --------
    def load_fusions(self, fusion_file):
        logger.info(f"Carregando fusões: {fusion_file}")
        try:
            fusions = pd.read_csv(fusion_file, sep='\t')
        except Exception as e:
            logger.warning(f"  Erro carregando fusões: {e}")
            return None

        defining = fusions[fusions.get('annot_classification', pd.Series(dtype=str)).astype(str).str.startswith('WHO:', na=False)]
        logger.info(f"  Total fusões (pós-filtro confidence do Arriba): {len(fusions)}")
        logger.info(f"  Definidoras de entidade WHO-HAEM5: {len(defining)}")
        self.fusion_rows = defining.to_dict('records')
        return fusions

    def load_cnvs(self, cnv_file):
        logger.info(f"Carregando CNVs: {cnv_file}")
        try:
            cnvs = pd.read_csv(cnv_file, sep='\t')
        except Exception as e:
            logger.warning(f"  Erro carregando CNVs: {e}")
            return None
        logger.info(f"  Total segmentos: {len(cnvs)}")
        self.cnvs_df = cnvs
        self._detect_karyotype_complexity(cnvs)
        self._detect_focal_deletions(cnvs)
        return cnvs

    def _detect_focal_deletions(self, cnvs):
        """Sobreposição de coordenada (mesma estratégia usada em toda a
        consolidação de resultados) para del17p/del11q/del13q14/del5q/del20q
        -- usados por cll_risk_stratification/ipss_r_lite.
        Trissomia 12 e monossomia 7 reaproveitam self.chrom_calls (nível de
        cromossomo inteiro, já calculado por _detect_karyotype_complexity)."""
        if cnvs is None or len(cnvs) == 0 or 'chrom' not in cnvs.columns or 'cnv_type' not in cnvs.columns:
            return
        # Mesma salvaguarda de qualidade usada para cariótipo complexo/monossômico
        # no ELN2022 (ver _detect_karyotype_complexity): se o CNVkit está com desvio
        # sistemático de calibração, TODAS as janelas focais (del17p/del11q/etc.)
        # herdam o mesmo falso sinal -- não populamos focal_dels nesse caso, para
        # a camada 1b não ficar vulnerável ao mesmo artefato já blindado na
        # camada ELN2022.
        if self.karyotype_quality_flag:
            logger.warning(
                "  Deleções/ganhos focais (LLC/SMD): ignorados -- "
                f"{self.karyotype_quality_flag}"
            )
            return
        for name, (chrom, start, end) in CLL_MDS_FOCAL_COORDS.items():
            sub = cnvs[(cnvs['chrom'] == chrom) & (cnvs['start'] < end) & (cnvs['end'] > start)]
            if len(sub) == 0:
                continue
            region_len = end - start
            lost_len = sum(max(0, min(r['end'], end) - max(r['start'], start))
                            for _, r in sub[sub['cnv_type'] == 'loss'].iterrows())
            if region_len > 0 and lost_len / region_len > 0.5:
                self.focal_dels.add(name)
        if self.chrom_calls.get('chr12') == 'gain':
            self.focal_dels.add('trisomy12')
        logger.info(f"  Deleções/ganhos focais (LLC/SMD): {sorted(self.focal_dels) or 'nenhum'}")

    # ExonicFunc.refGene values worth reporting. Intronic/intergenic/UTR/
    # synonymous variants are excluded, as are those failing the gnomAD
    # population-frequency filter, since GERMLINE COMMON VARIANTS ARE NOT DRIVERS.
    _EXONIC_FUNCS_OF_INTEREST = {
        'nonsynonymous SNV', 'frameshift deletion', 'frameshift insertion',
        'stopgain', 'stoploss', 'nonframeshift deletion', 'nonframeshift insertion',
        'startloss',
    }
    _GNOMAD_AF_MAX = 0.01

    def load_snvs(self, snv_file):
        """Carrega e filtra o multianno do ANNOVAR/Mutect2 diretamente.

        NOTA: lê o multianno do ANNOVAR diretamente (PASS + função exônica
        relevante + AF populacional <1%), não
        merged_variants/{sample}_driver_variants_final.tsv (saída de
        merge_mutect2_vardict -> prioritize_driver_variants), que opera sobre
        o VCF bruto do Mutect2/VarDict sem filtro de frequência populacional
        -- sem esse filtro, dezenas de genes driver aparecem "mutados" na
        mesma amostra, a maioria com VAF alto, consistente com polimorfismo
        germinativo comum, não driver somático.
        """
        logger.info(f"Carregando SNVs/indels anotados (ANNOVAR): {snv_file}")
        try:
            df = pd.read_csv(snv_file, sep='\t', low_memory=False)
        except Exception as e:
            logger.warning(f"  Erro carregando SNVs: {e}")
            return None

        gnomad_col = 'gnomad41_genome_AF' if 'gnomad41_genome_AF' in df.columns else None
        filt_col = 'Otherinfo10' if 'Otherinfo10' in df.columns else None  # FILTER
        fmt_col = 'Otherinfo12' if 'Otherinfo12' in df.columns else None  # FORMAT
        sample_col = 'Otherinfo13' if 'Otherinfo13' in df.columns else None  # sample GT field

        mask = pd.Series(True, index=df.index)
        if filt_col:
            mask &= (df[filt_col] == 'PASS')
        if 'ExonicFunc.refGene' in df.columns:
            mask &= df['ExonicFunc.refGene'].isin(self._EXONIC_FUNCS_OF_INTEREST)
        if gnomad_col:
            af = pd.to_numeric(df[gnomad_col], errors='coerce').fillna(0)
            mask &= (af < self._GNOMAD_AF_MAX)

        candidates = df[mask]
        logger.info(f"  Total variantes no arquivo: {len(df)}")
        logger.info(f"  PASS + exônica relevante + AF populacional <1%: {len(candidates)}")

        for _, row in candidates.iterrows():
            gene = str(row.get('Gene.refGene', '')).split(';')[0]
            if not gene or gene == 'nan' or gene not in self.driver_genes:
                continue
            clnsig = str(row.get('CLNSIG', '.'))
            af = 'NA'
            if fmt_col and sample_col and pd.notna(row.get(fmt_col)) and pd.notna(row.get(sample_col)):
                fk = str(row[fmt_col]).split(':')
                sv = str(row[sample_col]).split(':')
                if 'AF' in fk and len(sv) == len(fk):
                    af = sv[fk.index('AF')]
            aachg = _dedupe_protein_changes(row.get('AAChange.refGene', ''))
            mutation = f"{row.get('Chr','')}:{row.get('Start','')} {aachg} (VAF={af}, ClinVar={clnsig})"
            impact = 'HIGH' if 'athogenic' in clnsig else self._assess_mutation_impact_annovar(gene, row)

            prev = self.mutated_gene_details.get(gene)
            rank = {'HIGH': 2, 'MEDIUM': 1}.get(impact, 0)
            if prev is None or rank > {'HIGH': 2, 'MEDIUM': 1}.get(prev.get('impact'), 0):
                self.mutated_gene_details[gene] = {'mutation': mutation, 'impact': impact}
        logger.info(f"  Genes do painel driver ({len(self.driver_genes)}) com variante retida: {len(self.mutated_gene_details)}")
        return df

    _DRIVER_GENE_NAMES_CACHE = None

    @classmethod
    def _driver_gene_names(cls):
        return set(GENE_CLINICAL_NOTES.keys()) | ELN2022_ADVERSE_SOFT_GENES | ELN2022_ADVERSE_HARD_GENES | {'NPM1', 'CEBPA'}

    @staticmethod
    def _assess_mutation_impact_annovar(gene, row):
        exfunc = str(row.get('ExonicFunc.refGene', '')).lower()
        aachg = str(row.get('AAChange.refGene', '')).upper()
        high_impact_markers = {
            'FLT3': ['ins', 'itd'], 'NPM1': ['frameshift'], 'CEBPA': ['frameshift', 'stopgain'],
            'TP53': ['*'], 'RUNX1': ['frameshift', 'stopgain'], 'ASXL1': ['frameshift', 'stopgain'],
        }
        if gene == 'FLT3' and ('insertion' in exfunc or 'ITD' in aachg):
            return 'HIGH'
        if 'frameshift' in exfunc or 'stopgain' in exfunc:
            return 'HIGH' if gene in high_impact_markers else 'MEDIUM'
        return 'MEDIUM'

    def load_aneuploidy_summary(self, summary_file):
        if not summary_file:
            return None
        logger.info(f"Carregando aneuploidias: {summary_file}")
        try:
            with open(summary_file) as f:
                lines = f.readlines()
            aneuploidy = {}
            for line in lines:
                if ':' in line:
                    key, val = line.split(':', 1)
                    aneuploidy[key.strip()] = val.strip()
            return aneuploidy
        except Exception as e:
            logger.warning(f"  Erro carregando aneuploidias: {e}")
            return None

    @staticmethod
    def _assess_mutation_impact(gene, row):
        high_impact_markers = {
            'FLT3': ['itd', 'insertion'], 'NPM1': ['frameshift', 'insertion'],
            'CEBPA': ['frameshift', 'nonsense'], 'TP53': ['missense', 'frameshift', 'nonsense', 'stopgain'],
            'RUNX1': ['frameshift', 'nonsense'], 'ASXL1': ['frameshift', 'nonsense'],
        }
        blob = f"{row.get('Alt','')} {row.get('Priority','')}".lower()
        for pattern in high_impact_markers.get(gene, []):
            if pattern in blob:
                return 'HIGH'
        return 'MEDIUM'

    def _detect_karyotype_complexity(self, cnvs):
        """Conta alterações cromossômicas distintas (não segmentos individuais
        do CNVkit, que podem ser dezenas por cromossomo por ruído de
        segmentação -- ver nota de qualidade do CNVkit na tese/Seção 6.6).
        Cariótipo complexo (ELN2022/OMS): >=3 anormalidades cromossômicas
        NÃO relacionadas. Aproximamos usando 1 chamada por cromossomo/braço
        com alteração predominante, não por segmento."""
        if cnvs is None or len(cnvs) == 0 or 'chrom' not in cnvs.columns:
            return
        chrom_calls = {}
        for chrom in [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']:
            sub = cnvs[cnvs['chrom'] == chrom]
            if len(sub) == 0 or 'cnv_type' not in sub.columns:
                continue
            # pondera por tamanho do segmento, não por contagem de segmentos
            sub = sub.copy()
            sub['length'] = (sub['end'] - sub['start']).clip(lower=0)
            gained = sub.loc[sub['cnv_type'] == 'gain', 'length'].sum()
            lost = sub.loc[sub['cnv_type'] == 'loss', 'length'].sum()
            chrom_len = sub['length'].sum()
            if chrom_len == 0:
                continue
            if lost / chrom_len > 0.5:
                chrom_calls[chrom] = 'loss'
            elif gained / chrom_len > 0.5:
                chrom_calls[chrom] = 'gain'
        n_alterations = len(chrom_calls)
        n_loss = sum(1 for t in chrom_calls.values() if t == 'loss')
        n_gain = sum(1 for t in chrom_calls.values() if t == 'gain')
        self.aneuploidy_events = [f"{c}:{t}" for c, t in chrom_calls.items()]
        self.chrom_calls = chrom_calls  # reused por cll_risk_stratification/ipss_r_lite (trissomia12, -7, del20q)

        # Salvaguarda de qualidade #1: uma referência do CNVkit mal calibrada
        # produz um desvio de log2 GLOBAL, fazendo quase todo o genoma
        # parecer "perdido" (ou "ganho") de uma vez -- isso NÃO é um
        # cariótipo complexo real (que tem uma MISTURA de eventos em
        # cromossomos distintos), e sim um artefato técnico (ver Seção 6.6
        # sobre a correção de GC/RepeatMasker do CNVkit). Se >=18 dos 24
        # cromossomos (75%) apontam na MESMA direção, tratamos como
        # inconclusivo -- isso também suspende a detecção de deleções focais
        # (LLC/IPSS-R), que herdariam o mesmo sinal corrompido (ver
        # _detect_focal_deletions).
        if n_alterations > 0 and max(n_loss, n_gain) / 24 >= 0.75:
            self.karyotype_quality_flag = (
                f"INCONCLUSIVO: {max(n_loss, n_gain)}/24 cromossomos na mesma direção "
                f"({'perda' if n_loss > n_gain else 'ganho'}) -- padrão consistente com desvio "
                "sistemático de calibração do CNVkit (ver Seção 6.6), não com cariótipo complexo "
                "real."
            )
            logger.warning(f"  {self.karyotype_quality_flag}")
        else:
            self.karyotype_quality_flag = None

        # Salvaguarda de qualidade #2: mesmo após corrigir a calibração de GC
        # (Seção 6.6), a agregação por braço cromossômico (>50% do
        # comprimento em ganho/perda) não discrimina caso simples de caso
        # complexo com a referência atual -- sem um painel de normais real
        # (impossível sem amostra normal pareada
        # normal), o CNVkit "flat"/sintético tem ruído residual suficiente
        # para inflar essa métrica em praticamente qualquer amostra (visível,
        # por exemplo, em linhagens com cariótipo simples bem documentado na
        # literatura, como KASUMI-1/SKNO-1, t(8;21)/RUNX1::RUNX1T1). As
        # deleções FOCAIS (del17p/del11q/CDKN2A/SOCS1/etc., Camada 1b) não
        # têm esse problema -- por isso só o cariótipo complexo/monossômico
        # (nível de braço cromossômico) é suspenso aqui, não a detecção
        # focal. self.chrom_calls continua populado e reportado (Camada 3,
        # informativo), só não alimenta mais a Camada 1 (ELN2022 Adverso).
        self.karyotype_arm_level_unreliable = True
        self.complex_karyotype = False
        self.monosomal_karyotype = False
        logger.info(f"  Alterações cromossômicas (>50% do cromossomo, por comprimento): {n_alterations} {chrom_calls} "
                    "-- NÃO usado para cariótipo complexo/monossômico do ELN2022 (ver nota de qualidade #2)")

    # -------- orchestration --------
    def classify(self):
        mutated_genes = set(self.mutated_gene_details.keys())
        fusion_names = set()
        for f in self.fusion_rows:
            fusion_names.add(str(f.get('fusion', '')).upper())

        npm1_mut = 'NPM1' in mutated_genes
        cebpa_bzip = 'CEBPA' in mutated_genes
        flt3_itd = self.mutated_gene_details.get('FLT3', {}).get('impact') == 'HIGH'

        # Camada 1: ELN 2022
        eln_group, eln_factors = eln2022_risk(
            fusion_names, mutated_genes, cebpa_bzip, npm1_mut, flt3_itd,
            self.complex_karyotype, self.monosomal_karyotype, adverse_cytogenetics=set(),
        )

        # Camada 1b: risco genômico de LLC e componente citogenético/mutacional
        # do IPSS-R (SMD).
        cll_risk = cll_risk_stratification(mutated_genes, self.focal_dels, self.karyotype_quality_flag)
        mds_risk = ipss_r_lite(mutated_genes, self.focal_dels, self.chrom_calls,
                                self.complex_karyotype, self.karyotype_quality_flag)

        # Camada 2: entidade WHO-HAEM5/ICC
        entity = assign_who_entity(self.fusion_rows, self.mutated_gene_details, self.cnvs_df)

        # Camada 3: achados adicionais
        used_genes = {'NPM1', 'CEBPA'} if entity else set()
        layer3 = build_layer3_findings(self.mutated_gene_details, used_genes)

        return {
            'sample_id': self.sample_id,
            'layer1_eln2022_risk': {
                'applicable_note': 'Estratificação ELN 2022 é validada para LMA; reportada para toda amostra '
                                    'com achados genéticos compatíveis -- a aplicabilidade clínica depende de '
                                    'confirmação morfológica/imunofenotípica de LMA.',
                'risk_group': eln_group,
                'factors': eln_factors,
            },
            'layer1b_cll_risk': cll_risk,
            'layer1b_mds_ipssr_lite': mds_risk,
            'layer2_who_icc_entity': entity if entity else {
                'entity': None,
                'note': 'Nenhuma alteração genética definidora de entidade (WHO-HAEM5/ICC 2022) foi detectada '
                        'nos dados de sequenciamento. Isso NÃO exclui diagnóstico -- classificação morfológica '
                        'e imunofenotípica (fora do escopo do VKaryo) pode ainda definir a entidade.',
            },
            'layer3_additional_findings': layer3,
            'karyotype_summary': {
                'complex_karyotype': self.complex_karyotype,
                'monosomal_karyotype': self.monosomal_karyotype,
                'chromosome_level_events': self.aneuploidy_events,
                'quality_flag': self.karyotype_quality_flag,
                'arm_level_unreliable': self.karyotype_arm_level_unreliable,
                'arm_level_note': (
                    'Eventos em nível de braço cromossômico (ganho/perda) listados acima são '
                    'informativos apenas -- não alimentam o veredito de cariótipo complexo/'
                    'monossômico do ELN 2022. Motivo: sem um painel de amostras normais pareadas '
                    'para calibração da referência do CNVkit, a magnitude de ganhos/perdas em nível '
                    'de braço cromossômico inteiro não atinge confiabilidade suficiente para uso '
                    'automatizado nessa decisão de risco. Deleções/ganhos focais em genes '
                    'específicos (Camada 1b) não têm essa limitação.'
                ) if self.karyotype_arm_level_unreliable else None,
            },
            'additional_fusions_not_who_defining': [
                f.get('fusion') for f in [] # placeholder kept for backward-compat key presence
            ],
        }

    def generate_recommendations(self, result):
        recs = []
        entity = result['layer2_who_icc_entity']
        if entity.get('entity'):
            for f in self.fusion_rows[:1]:
                therapy = f.get('annot_therapy')
                if therapy and therapy != 'N/A':
                    recs.append({'category': 'Terapia direcionada à entidade',
                                 'recommendation': therapy,
                                 'indication': f"Entidade WHO-HAEM5: {entity['entity']}"})
        for gene in ('FLT3', 'IDH1', 'IDH2', 'BTK', 'MYD88'):
            if gene in self.mutated_gene_details:
                recs.append({'category': 'Alvo terapêutico',
                              'recommendation': GENE_CLINICAL_NOTES.get(gene, ''),
                              'indication': f'Mutação em {gene}'})
        risk = result['layer1_eln2022_risk']['risk_group']
        if risk == 'Favorável':
            recs.append({'category': 'Risco ELN 2022', 'recommendation':
                          'Indução padrão + consolidação com altas doses de citarabina',
                          'indication': 'Grupo de risco favorável'})
        elif risk == 'Adverso':
            recs.append({'category': 'Risco ELN 2022', 'recommendation':
                          'Considerar transplante alogênico em primeira remissão',
                          'indication': 'Grupo de risco adverso'})
        else:
            recs.append({'category': 'Risco ELN 2022', 'recommendation':
                          'Protocolo padrão; considerar transplante em segunda remissão ou conforme MRD',
                          'indication': 'Grupo de risco intermediário'})
        if not recs:
            recs.append({'category': 'Padrão', 'recommendation':
                          'Seguir protocolo institucional; nenhuma alteração acionável específica detectada',
                          'indication': 'N/A'})
        return recs


def render_markdown_report(result, recommendations, sample_id):
    lines = [f"# Relatório de Estratificação Genômica — Amostra {sample_id}",
             "", "_Gerado automaticamente pelo módulo de estratificação do VKaryo "
             "(ELN 2022 + WHO-HAEM5/ICC 2022 + achados adicionais curados). "
             "Requer revisão por hematologista/bioinformata qualificado antes de uso clínico._", ""]

    lines.append("## Camada 1 — Risco ELN 2022")
    l1 = result['layer1_eln2022_risk']
    lines.append(f"**Grupo de risco: {l1['risk_group']}**")
    lines.append("")
    if l1['factors']:
        for f in l1['factors']:
            lines.append(f"- {f}")
    else:
        lines.append("- Nenhum fator classificador encontrado.")
    lines.append("")
    lines.append(f"_{l1['applicable_note']}_")
    lines.append("")

    l1b_cll = result['layer1b_cll_risk']
    lines.append("## Camada 1b — Risco genômico de LLC")
    lines.append(f"**Grupo: {l1b_cll['risk_group']}**")
    lines.append("")
    for f in l1b_cll['factors']:
        lines.append(f"- {f}")
    lines.append("")
    lines.append(f"_{l1b_cll['applicable_note']}_")
    lines.append("")

    l1b_mds = result['layer1b_mds_ipssr_lite']
    lines.append("## Camada 1b — IPSS-R (SMD), componente genômico")
    lines.append(f"**Categoria: {l1b_mds['category']}** (grupo citogenético: {l1b_mds['cytogenetic_group']})")
    lines.append("")
    for f in l1b_mds['factors']:
        lines.append(f"- {f}")
    if not l1b_mds['factors']:
        lines.append("- Nenhum marcador citogenético/molecular do IPSS-R detectado.")
    lines.append("")
    lines.append(f"_{l1b_mds['applicable_note']}_")
    lines.append("")

    lines.append("## Camada 2 — Entidade WHO-HAEM5 / ICC 2022")
    l2 = result['layer2_who_icc_entity']
    if l2.get('entity'):
        lines.append(f"**Entidade: {l2['entity']}**")
        lines.append("")
        lines.append(f"- Classificação: {l2.get('who_classification', 'N/A')}")
        lines.append(f"- Fonte da evidência: {l2.get('source', 'N/A')}")
        for e in l2.get('evidence', []):
            lines.append(f"- Evidência: {e}")
        if l2.get('icc_note'):
            lines.append(f"- Nota ICC: {l2['icc_note']}")
        if l2.get('additional_fusions_same_sample'):
            lines.append(f"- Outras fusões definidoras na mesma amostra: {', '.join(l2['additional_fusions_same_sample'])} "
                          "(possível biologia duplo/triplo-hit — revisar manualmente)")
    else:
        lines.append("**Nenhuma entidade atribuída.**")
        lines.append("")
        lines.append(f"_{l2.get('note', '')}_")
    lines.append("")

    lines.append("## Camada 3 — Achados adicionais em genes driver")
    l3 = result['layer3_additional_findings']
    if l3:
        lines.append("| Gene | Alteração | Impacto | Nota |")
        lines.append("|---|---|---|---|")
        for f in l3:
            lines.append(f"| {f['gene']} | {f['mutation']} | {f['impact']} | {f['clinical_note']} |")
    else:
        lines.append("Nenhum achado adicional em genes driver além das camadas 1-2.")
    lines.append("")

    ks = result['karyotype_summary']
    lines.append("## Cariótipo virtual (resumo)")
    if ks.get('quality_flag'):
        lines.append(f"⚠️ **{ks['quality_flag']}**")
        lines.append("")
    lines.append(f"- Cariótipo complexo (≥3 alterações cromossômicas não relacionadas, por >50% do comprimento "
                  f"do cromossomo): **{'Sim' if ks['complex_karyotype'] else 'Não'}**"
                  + (" _(nível de braço não usado como critério -- ver nota abaixo)_" if ks.get('arm_level_unreliable') else ""))
    lines.append(f"- Cariótipo monossômico: **{'Sim' if ks['monosomal_karyotype'] else 'Não'}**")
    if ks['chromosome_level_events']:
        lines.append(f"- Eventos em nível cromossômico (informativo, ver nota): {', '.join(ks['chromosome_level_events'])}")
    lines.append("")
    if ks.get('arm_level_note'):
        lines.append(f"⚠️ _{ks['arm_level_note']}_")
        lines.append("")
    lines.append("_Nota de qualidade: deleções/ganhos FOCAIS (gene específico, Camada 1b) usam correção de "
                  "conteúdo GC e são consideradas confiáveis. Eventos em nível de braço cromossômico "
                  "inteiro seguem informativos apenas (ver nota acima)._")
    lines.append("")

    lines.append("## Recomendações")
    for r in recommendations:
        lines.append(f"- **[{r['category']}]** {r['recommendation']} _(indicação: {r['indication']})_")
    lines.append("")
    lines.append("---")
    lines.append("_Este relatório é gerado computacionalmente e não substitui a interpretação de um "
                 "bioinformata/hematologista qualificado (ver limitações discutidas na tese do VKaryo)._")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description='Estratificação de leucemia/linfoma (ELN2022 + WHO-HAEM5/ICC + achados adicionais)')
    parser.add_argument('--fusions', help='Fusões anotadas TSV (arriba_annotated)')
    parser.add_argument('--cnvs', help='CNVs validados TSV (cnv_validated)')
    parser.add_argument('--snvs', help='ANNOVAR multianno TSV do Mutect2 (mutect2/annovar/*_mutect2.hg38_multianno.txt)')
    parser.add_argument('--driver-genes', help='Lista de genes driver (reference/driver_genes_leukemia.txt); '
                         'se omitido, usa como fallback as chaves de GENE_CLINICAL_NOTES')
    parser.add_argument('--aneuploidy-summary', help='Sumário de aneuploidias TXT')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Arquivo JSON de saída')
    parser.add_argument('--report', help='Arquivo Markdown de saída (relatório legível)')
    args = parser.parse_args()

    classifier = LeukemiaClassifier(args.sample, driver_genes_file=args.driver_genes)
    if args.fusions:
        classifier.load_fusions(args.fusions)
    if args.cnvs:
        classifier.load_cnvs(args.cnvs)
    if args.snvs:
        classifier.load_snvs(args.snvs)
    if args.aneuploidy_summary:
        classifier.load_aneuploidy_summary(args.aneuploidy_summary)

    result = classifier.classify()
    recommendations = classifier.generate_recommendations(result)
    result['recommendations'] = recommendations

    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    logger.info(f"✓ Classificação salva: {args.output}")

    if args.report:
        md = render_markdown_report(result, recommendations, args.sample)
        with open(args.report, 'w') as f:
            f.write(md)
        logger.info(f"✓ Relatório salvo: {args.report}")

    logger.info("\n" + "=" * 60)
    logger.info("SUMÁRIO")
    logger.info("=" * 60)
    logger.info(f"Risco ELN 2022: {result['layer1_eln2022_risk']['risk_group']}")
    logger.info(f"Entidade WHO-HAEM5/ICC: {result['layer2_who_icc_entity'].get('entity') or 'Não atribuída'}")
    logger.info(f"Achados adicionais (camada 3): {len(result['layer3_additional_findings'])}")


if __name__ == '__main__':
    main()
