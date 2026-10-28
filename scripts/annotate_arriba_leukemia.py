#!/usr/bin/env python3
"""
Anota fusões do Arriba com significância clínica em leucemias
Adiciona informações sobre prognóstico, terapia, e classificação

Uso:
    python annotate_arriba_leukemia.py \
        --input output_vkaryo/arriba/ERR123_fusions.tsv \
        --output output_vkaryo/arriba_annotated/ERR123_fusions_annotated.tsv
"""

import pandas as pd
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


# Database de fusões recorrentes em leucemias
LEUKEMIA_FUSIONS = {
    'BCR--ABL1': {
        'alias': ['ABL1--BCR'],
        'leukemia_type': 'CML / Ph+ ALL',
        'classification': 'WHO: BCR::ABL1-positive',
        'frequency': 'CML: 95%, ALL: 25% adultos, 3-5% pediátrico',
        'prognosis': 'Bom com TKI',
        'risk_group': 'Intermediário (melhora com TKI)',
        'therapy': 'Imatinib, Dasatinib, Nilotinib, Ponatinib',
        'monitoring': 'BCR-ABL1 quantitativo a cada 3 meses',
        'clinical_significance': 'PATOGNOMÔNICA - Diagnóstica e alvo terapêutico',
        'references': 'PMID: 15659725, 16258524'
    },

    'PML--RARA': {
        'alias': ['RARA--PML'],
        'leukemia_type': 'APL (Leucemia Promielocítica Aguda - M3)',
        'classification': 'WHO: APL with PML::RARA',
        'frequency': '95-98% das APL',
        'prognosis': 'Excelente (>90% cura)',
        'risk_group': 'Favorável',
        'therapy': 'ATRA + Arsênico trióxido (ATO)',
        'monitoring': 'PML-RARA por RT-PCR',
        'clinical_significance': 'PATOGNOMÔNICA - Emergência oncológica',
        'references': 'PMID: 23512909, 30180821'
    },

    'RUNX1--RUNX1T1': {
        'alias': ['RUNX1T1--RUNX1', 'AML1--ETO', 'ETO--AML1'],
        'leukemia_type': 'AML com t(8;21)(q22;q22)',
        'classification': 'WHO: AML with RUNX1::RUNX1T1',
        'frequency': '5-10% AML adultos, 12-15% pediátrico',
        'prognosis': 'Favorável',
        'risk_group': 'Favorável (ELN 2022)',
        'therapy': 'Quimioterapia intensiva + altas doses citarabina',
        'monitoring': 'RUNX1-RUNX1T1 por RT-PCR',
        'clinical_significance': 'PATOGNOMÔNICA - Favorável prognóstico',
        'references': 'PMID: 27895058, 35797463'
    },

    'CBFB--MYH11': {
        'alias': ['MYH11--CBFB'],
        'leukemia_type': 'AML com inv(16)(p13.1q22) ou t(16;16)',
        'classification': 'WHO: AML with CBFB::MYH11',
        'frequency': '5-8% AML',
        'prognosis': 'Favorável',
        'risk_group': 'Favorável (ELN 2022)',
        'therapy': 'Quimioterapia intensiva + altas doses citarabina',
        'monitoring': 'CBFB-MYH11 por RT-PCR',
        'clinical_significance': 'PATOGNOMÔNICA - Core-binding factor AML',
        'references': 'PMID: 27895058'
    },

    'KMT2A--MLLT3': {
        'alias': ['MLLT3--KMT2A', 'MLL--AF9', 'AF9--MLL'],
        'leukemia_type': 'AML com t(9;11)(p21.3;q23.3)',
        'classification': 'WHO: AML with KMT2A rearrangement',
        'frequency': '2-3% AML',
        'prognosis': 'Intermediário',
        'risk_group': 'Intermediário (ELN 2022)',
        'therapy': 'Quimioterapia intensiva ± transplante',
        'monitoring': 'KMT2A-MLLT3 por RT-PCR ou FISH',
        'clinical_significance': 'MLL-rearranged leukemia',
        'references': 'PMID: 27895058'
    },

    'DEK--NUP214': {
        'alias': ['NUP214--DEK'],
        'leukemia_type': 'AML com t(6;9)(p22.3;q34.1)',
        'classification': 'WHO: AML with DEK::NUP214',
        'frequency': '<1% AML',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso (ELN 2022)',
        'therapy': 'Quimioterapia intensiva + transplante alogênico',
        'monitoring': 'Morfologia + citogenética',
        'clinical_significance': 'Prognóstico desfavorável',
        'references': 'PMID: 27895058'
    },

    'ETV6--RUNX1': {
        'alias': ['RUNX1--ETV6', 'TEL--AML1', 'AML1--TEL'],
        'leukemia_type': 'B-ALL pediátrica',
        'classification': 'WHO: B-ALL with ETV6::RUNX1',
        'frequency': '25% ALL pediátrica',
        'prognosis': 'Favorável',
        'risk_group': 'Favorável',
        'therapy': 'Protocolo ALL padrão',
        'monitoring': 'MRD por flow cytometry ou PCR',
        'clinical_significance': 'Prognóstico favorável em ALL infantil',
        'references': 'PMID: 27895058'
    },

    'TCF3--PBX1': {
        'alias': ['PBX1--TCF3', 'E2A--PBX1', 'PBX1--E2A'],
        'leukemia_type': 'B-ALL com t(1;19)(q23;p13.3)',
        'classification': 'WHO: B-ALL with TCF3::PBX1',
        'frequency': '5-6% ALL',
        'prognosis': 'Intermediário (melhorou com terapia moderna)',
        'risk_group': 'Intermediário',
        'therapy': 'Protocolo ALL padrão',
        'monitoring': 'MRD',
        'clinical_significance': 'Resposta boa com terapia moderna',
        'references': 'PMID: 27895058'
    },

    'BCR--ABL2': {
        'alias': ['ABL2--BCR'],
        'leukemia_type': 'B-ALL Ph-like',
        'classification': 'WHO: B-ALL BCR::ABL1-like',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'TKI (Dasatinib) + quimioterapia',
        'monitoring': 'MRD',
        'clinical_significance': 'Ph-like ALL - considerar TKI',
        'references': 'PMID: 27895058'
    },

    'NUP98--NSD1': {
        'alias': ['NSD1--NUP98'],
        'leukemia_type': 'AML',
        'classification': 'WHO: AML with NUP98 rearrangement',
        'frequency': 'Raro (<1%)',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'Quimioterapia intensiva + transplante',
        'monitoring': 'Morfologia + MRD',
        'clinical_significance': 'PATOGNOMÔNICA - Associada a FLT3-ITD, prognóstico ruim',
        'references': 'PMID: 29439119'
    },

    'NPM1--ALK': {
        'alias': ['ALK--NPM1'],
        'leukemia_type': 'Linfoma Anaplásico de Grandes Células (ALCL), ALK+',
        'classification': 'WHO: ALCL, ALK-positive',
        'frequency': '~80% dos ALCL sistêmicos',
        'prognosis': 'Favorável (melhor que ALCL ALK-negativo)',
        'risk_group': 'Favorável',
        'therapy': 'CHOP + considerar inibidor de ALK (crizotinibe)',
        'monitoring': 'NPM1-ALK por RT-PCR/FISH',
        'clinical_significance': 'PATOGNOMÔNICA - Define ALCL ALK-positive',
        'references': 'PMID: 30990397'
    },
    'CCND1--IGH': {
        'alias': ['IGH--CCND1'],
        'leukemia_type': 'Linfoma de Células do Manto (LCM)',
        'classification': 'WHO: Mantle cell lymphoma with CCND1::IGH',
        'frequency': '~95% dos LCM',
        'prognosis': 'Variável (depende de Ki-67, TP53, blastoide)',
        'risk_group': 'Intermediário (MIPI)',
        'therapy': 'Quimioimunoterapia ± inibidor de BTK',
        'monitoring': 'MRD por PCR (IGH-CCND1)',
        'clinical_significance': 'PATOGNOMÔNICA - Define linfoma de células do manto',
        'references': 'WHO-HAEM5 2022'
    },
    'MYC--IGH': {
        'alias': ['IGH--MYC'],
        'leukemia_type': 'Linfoma de Burkitt (se rearranjo único) ou componente de linfoma B de alto grau',
        'classification': 'WHO: Burkitt lymphoma (MYC::IGH) / High-grade B-cell lymphoma with MYC and BCL2/BCL6 rearrangements se concomitante',
        'frequency': '100% Burkitt; variável em DLBCL',
        'prognosis': 'Agressivo, potencialmente curável (Burkitt)',
        'risk_group': 'Adverso se duplo/triplo-hit (concomitante BCL2/BCL6)',
        'therapy': 'Quimioterapia intensiva tipo Burkitt (ex: DA-EPOCH-R, CODOX-M/IVAC)',
        'monitoring': 'MRD, imagem',
        'clinical_significance': 'Rearranjo de MYC - checar concomitância com BCL2/BCL6 (duplo/triplo-hit)',
        'references': 'WHO-HAEM5 2022; ICC 2022'
    },
    'BCL2--IGH': {
        'alias': ['IGH--BCL2'],
        'leukemia_type': 'Linfoma Folicular / componente de DLBCL t(14;18)',
        'classification': 'WHO: Follicular lymphoma with BCL2::IGH (ou DLBCL GC se contexto difuso)',
        'frequency': '~85% dos linfomas foliculares',
        'prognosis': 'Variável',
        'risk_group': 'Adverso se concomitante com rearranjo de MYC (duplo-hit)',
        'therapy': 'Depende do contexto histológico (indolente vs. agressivo)',
        'monitoring': 'MRD, imagem',
        'clinical_significance': 'Rearranjo de BCL2 - checar concomitância com MYC (duplo-hit) e BCL6',
        'references': 'WHO-HAEM5 2022'
    },
    'BCL6--IGH': {
        'alias': ['IGH--BCL6'],
        'leukemia_type': 'Componente de DLBCL / linfoma B de alto grau',
        'classification': 'WHO: rearranjo de BCL6 (contribui para "double/triple expressor" ou high-grade B-cell lymphoma NOS)',
        'frequency': 'Variável em DLBCL',
        'prognosis': 'Variável',
        'risk_group': 'Adverso se concomitante com MYC e BCL2',
        'therapy': 'Depende do contexto',
        'monitoring': 'MRD, imagem',
        'clinical_significance': 'Rearranjo de BCL6 - checar concomitância com MYC/BCL2',
        'references': 'WHO-HAEM5 2022'
    },
    'FGFR1OP2--FGFR1': {
        'alias': ['FGFR1--FGFR1OP2'],
        'leukemia_type': 'Neoplasia mieloide/linfoide com rearranjo de FGFR1 (8p11 myeloproliferative syndrome)',
        'classification': 'WHO: Myeloid/lymphoid neoplasm with FGFR1 rearrangement',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável sem transplante',
        'risk_group': 'Adverso',
        'therapy': 'Inibidor de FGFR1 (ex: pemigatinibe) + transplante alogênico',
        'monitoring': 'FISH para rearranjo de FGFR1',
        'clinical_significance': 'PATOGNOMÔNICA - Define categoria WHO de neoplasia com rearranjo de FGFR1',
        'references': 'WHO-HAEM5 2022'
    },
    'SET--NUP214': {
        'alias': ['NUP214--SET'],
        'leukemia_type': 'AML (frequentemente M7) ou T-ALL',
        'classification': 'AML/T-ALL with SET::NUP214',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'Quimioterapia intensiva + considerar transplante',
        'monitoring': 'RT-PCR',
        'clinical_significance': 'Associada a superexpressão de HOXA/MEIS1',
        'references': 'PMID: 24124338'
    },
    'FIP1L1--PDGFRA': {
        'alias': ['PDGFRA--FIP1L1'],
        'leukemia_type': 'Leucemia Eosinofílica Crônica / AML eosinofílica',
        'classification': 'WHO: Myeloid/lymphoid neoplasm with PDGFRA rearrangement',
        'frequency': 'Principal causa de eosinofilia clonal',
        'prognosis': 'Excelente com terapia-alvo',
        'risk_group': 'Favorável com imatinibe',
        'therapy': 'Imatinibe (extremamente sensível, doses baixas)',
        'monitoring': 'FIP1L1-PDGFRA por RT-PCR (a deleção CHIC2 impede detecção por FISH convencional)',
        'clinical_significance': 'PATOGNOMÔNICA - Define categoria WHO de neoplasia com rearranjo de PDGFRA',
        'references': 'WHO-HAEM5 2022'
    },
    'RBM6--CSF1R': {
        'alias': ['CSF1R--RBM6'],
        'leukemia_type': 'AML megacariocítica com rearranjo de CSF1R',
        'classification': 'Rearranjo de CSF1R (receptor tirosina quinase)',
        'frequency': 'Raro',
        'prognosis': 'Desconhecido / potencialmente sensível a inibidor de CSF1R',
        'risk_group': 'Indeterminado',
        'therapy': 'Considerar inibidor de CSF1R + quimioterapia padrão',
        'monitoring': 'RT-PCR',
        'clinical_significance': 'Rearranjo raro de receptor tirosina quinase',
        'references': 'Literatura de caso'
    },
    'KMT2A--AFF1': {
        'alias': ['AFF1--KMT2A', 'MLL--AF4', 'AF4--MLL'],
        'leukemia_type': 'Pre-B-ALL com t(4;11)(q21;q23)',
        'classification': 'WHO: B-ALL with KMT2A::AFF1',
        'frequency': 'Predomina em ALL de lactentes',
        'prognosis': 'Desfavorável, especialmente em lactentes',
        'risk_group': 'Adverso',
        'therapy': 'Protocolo intensivo específico para KMT2A-r ± transplante',
        'monitoring': 'MRD',
        'clinical_significance': 'PATOGNOMÔNICA - Rearranjo de KMT2A de pior prognóstico em ALL',
        'references': 'WHO-HAEM5 2022'
    },
    'KMT2A--MLLT1': {
        'alias': ['MLLT1--KMT2A', 'MLL--ENL', 'ENL--MLL'],
        'leukemia_type': 'Pre-B-ALL/AML com t(11;19)(q23;p13)',
        'classification': 'WHO: leukemia with KMT2A::MLLT1',
        'frequency': 'Raro',
        'prognosis': 'Intermediário a desfavorável',
        'risk_group': 'Intermediário/Adverso',
        'therapy': 'Protocolo intensivo específico para KMT2A-r',
        'monitoring': 'MRD',
        'clinical_significance': 'Rearranjo de KMT2A',
        'references': 'WHO-HAEM5 2022'
    },
    'MLLT10--PICALM': {
        'alias': ['PICALM--MLLT10'],
        'leukemia_type': 'AML/T-ALL com t(10;11)(p13;q14)',
        'classification': 'AML/T-ALL with PICALM::MLLT10',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'Quimioterapia intensiva + considerar transplante',
        'monitoring': 'RT-PCR',
        'clinical_significance': 'Rearranjo recorrente em AML monocítica e T-ALL',
        'references': 'PMID: 26837932'
    },
    'CBFA2T3--ABHD12': {
        'alias': ['ABHD12--CBFA2T3'],
        'leukemia_type': 'AML eritroide (M6)',
        'classification': 'Rearranjo de CBFA2T3 em AML eritroide',
        'frequency': 'Raro',
        'prognosis': 'Desconhecido',
        'risk_group': 'Indeterminado',
        'therapy': 'Protocolo padrão AML',
        'monitoring': 'RT-PCR',
        'clinical_significance': 'Achado recorrente em AML eritroide',
        'references': 'Literatura de caso'
    },
    'STIL--TAL1': {
        'alias': ['TAL1--STIL', 'SIL--TAL1', 'TAL1--SIL'],
        'leukemia_type': 'T-ALL',
        'classification': 'WHO: T-ALL with STIL::TAL1 (del(1)(p32))',
        'frequency': '~25% dos T-ALL pediátricos',
        'prognosis': 'Variável',
        'risk_group': 'Intermediário',
        'therapy': 'Protocolo T-ALL padrão',
        'monitoring': 'MRD',
        'clinical_significance': 'PATOGNOMÔNICA - Driver oncogênico bem estabelecido em T-ALL',
        'references': 'WHO-HAEM5 2022'
    },
    'IGH--SPIB': {
        'alias': ['SPIB--IGH'],
        'leukemia_type': 'DLBCL ABC',
        'classification': 'Rearranjo de SPIB em DLBCL ABC',
        'frequency': 'Raro',
        'prognosis': 'Variável',
        'risk_group': 'Intermediário',
        'therapy': 'R-CHOP ou protocolo DLBCL',
        'monitoring': 'Imagem, LDH',
        'clinical_significance': 'Rearranjo recorrente em DLBCL do subtipo ABC',
        'references': 'Literatura de caso'
    },
    'TBL1XR1--TP63': {
        'alias': ['TP63--TBL1XR1'],
        'leukemia_type': 'Linfoma/leucemia de células T maduras',
        'classification': 'Rearranjo de TP63 em neoplasia de células T maduras',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'Protocolo específico para neoplasia de células T maduras',
        'monitoring': 'Morfologia + imunofenotipagem',
        'clinical_significance': 'Rearranjo recorrente de TP63 em neoplasias T maduras agressivas',
        'references': 'PMID: 25281473'
    },
    'CASP8--ERBB4': {
        'alias': ['ERBB4--CASP8'],
        'leukemia_type': 'Linfoma/leucemia de células T maduras',
        'classification': 'Rearranjo recorrente em neoplasia de células T madura',
        'frequency': 'Raro',
        'prognosis': 'Desconhecido',
        'risk_group': 'Indeterminado',
        'therapy': 'Protocolo específico para neoplasia de células T maduras',
        'monitoring': 'Morfologia + imunofenotipagem',
        'clinical_significance': 'Achado recorrente em neoplasias T maduras',
        'references': 'Literatura de caso'
    },
    'CIITA--C15ORF65': {
        'alias': ['C15ORF65--CIITA', 'CIITA--PIERCE2', 'PIERCE2--CIITA'],
        'leukemia_type': 'Linfoma de Hodgkin / linfoma B mediastinal primário',
        'classification': 'Rearranjo de CIITA (evasão imune por perda de MHC-II)',
        'frequency': '~15% do linfoma de Hodgkin clássico, ~38% do PMBL',
        'prognosis': 'Variável',
        'risk_group': 'Intermediário',
        'therapy': 'Protocolo padrão + considerar checkpoint inibidor (contexto de evasão imune)',
        'monitoring': 'Imagem (PET-CT)',
        'clinical_significance': 'Mecanismo de evasão imune por downregulação de MHC classe II',
        'references': 'PMID: 27540013'
    },
    'FOXK2--TP63': {
        'alias': ['TP63--FOXK2'],
        'leukemia_type': 'Linfoma cutâneo de células T / neoplasia de células T madura',
        'classification': 'Rearranjo de TP63 em neoplasia de células T madura',
        'frequency': 'Raro',
        'prognosis': 'Desfavorável',
        'risk_group': 'Adverso',
        'therapy': 'Protocolo específico para neoplasia de células T maduras',
        'monitoring': 'Morfologia + imunofenotipagem',
        'clinical_significance': 'Rearranjo recorrente de TP63 em CTCL/neoplasias T maduras',
        'references': 'PMID: 25281473'
    },
    'R3HCC1L--HTRA1': {
        'alias': ['HTRA1--R3HCC1L'],
        'leukemia_type': 'B-NHL: CLL/PLL',
        'classification': 'Rearranjo recorrente em CLL/PLL',
        'frequency': 'Raro',
        'prognosis': 'Desconhecido',
        'risk_group': 'Indeterminado',
        'therapy': 'Protocolo padrão CLL/PLL (ver estratificação por del17p/TP53/IGHV)',
        'monitoring': 'MRD',
        'clinical_significance': 'Achado recorrente, significado clínico não plenamente estabelecido',
        'references': 'Literatura de caso'
    },
}

# Genes cujo parceiro de fusão é o lócus de imunoglobulina (IGH/IGK/IGL) ou um
# gene "proxy" adjacente frequentemente usado por Arriba em vez do
# nome literal do segmento de Ig (ex: CASC11/PVT1 ficam ao lado de MYC,
# MYEOV/ELK2AP ao lado de CCND1). Mapeia gene-alvo -> chave em LEUKEMIA_FUSIONS.
IG_PARTNER_PROXY = {
    'MYC': 'MYC--IGH', 'CASC11': 'MYC--IGH', 'PVT1': 'MYC--IGH',
    'BCL2': 'BCL2--IGH', 'KDSR': 'BCL2--IGH',
    'BCL6': 'BCL6--IGH',
    'CCND1': 'CCND1--IGH', 'MYEOV': 'CCND1--IGH',
}
IG_SEGMENT_PREFIXES = ('IGH', 'IGK', 'IGL')


def normalize_fusion_name(gene1, gene2):
    """Normaliza nome da fusão para match no database"""
    # Remove isoformas/transcritos (ex: BCR(1) -> BCR)
    gene1_clean = gene1.split('(')[0].strip().upper()
    gene2_clean = gene2.split('(')[0].strip().upper()

    fusion_standard = f"{gene1_clean}--{gene2_clean}"
    fusion_reverse = f"{gene2_clean}--{gene1_clean}"

    return fusion_standard, fusion_reverse


def annotate_fusion(gene1, gene2):
    """Anota fusão com informações clínicas"""

    fusion_std, fusion_rev = normalize_fusion_name(gene1, gene2)

    # Busca no database
    if fusion_std in LEUKEMIA_FUSIONS:
        return LEUKEMIA_FUSIONS[fusion_std]

    if fusion_rev in LEUKEMIA_FUSIONS:
        return LEUKEMIA_FUSIONS[fusion_rev]

    # Busca por alias
    for fusion_name, info in LEUKEMIA_FUSIONS.items():
        if 'alias' in info:
            if fusion_std in info['alias'] or fusion_rev in info['alias']:
                return info

    # Parceiro de fusão é o lócus de IGH/IGK/IGL (ou um proxy adjacente do
    # Ensembl como CASC11/PVT1/MYEOV) -- comum em translocações de linfoma B
    # onde Arriba nomeia o segmento de Ig (IGHJ4, IGHM, IGHG1...) em vez do
    # gene MYC/BCL2/BCL6/CCND1 propriamente dito.
    gene1_clean, gene2_clean = fusion_std.split('--')
    for a, b in [(gene1_clean, gene2_clean), (gene2_clean, gene1_clean)]:
        if a in IG_PARTNER_PROXY and (b.startswith(IG_SEGMENT_PREFIXES) or b in
                                       ('IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6',
                                        'IGHD', 'IGHM', 'IGHE', 'IGHA1', 'IGHA2',
                                        'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHGP',
                                        'IGHEP1', 'IGHJ3P', 'IGHV1-69D', 'ELK2AP')):
            return LEUKEMIA_FUSIONS[IG_PARTNER_PROXY[a]]

    # Fusões envolvendo KMT2A (MLL) - sempre relevante.
    # NOTA: comparação exata do nome do gene, não substring ('MLL' in
    # fusion_std) -- substring dá falso positivo em qualquer gene começando
    # com "MLL" (ex. MLLT10, parceiro de fusão distinto, não o próprio
    # KMT2A/MLL1; um parceiro inesperado em entidade nao mieloide,
    # não uma leucemia KMT2A-rearranged). Mesma lição aplicada ao filtro
    # CDKN2A/CDKN2AIP.
    if gene1_clean in ('KMT2A', 'MLL') or gene2_clean in ('KMT2A', 'MLL'):
        return {
            'leukemia_type': 'AML/ALL with KMT2A rearrangement',
            'classification': 'WHO: KMT2A-rearranged',
            'frequency': 'Variável',
            'prognosis': 'Geralmente desfavorável',
            'risk_group': 'Intermediário a Adverso',
            'therapy': 'Quimioterapia intensiva + considerar transplante',
            'monitoring': 'KMT2A rearrangement por FISH',
            'clinical_significance': 'MLL-rearranged - prognóstico variável dependendo do parceiro',
            'references': 'PMID: 27895058'
        }

    # Fusões envolvendo genes de Ph-like ALL
    ph_like_genes = ['ABL1', 'ABL2', 'PDGFRA', 'PDGFRB', 'FGFR1', 'JAK2', 'EPOR', 'CRLF2']
    genes_in_fusion = [fusion_std.split('--')[0], fusion_std.split('--')[1]]

    if any(g in genes_in_fusion for g in ph_like_genes):
        return {
            'leukemia_type': 'Possível B-ALL Ph-like',
            'classification': 'Investigar perfil Ph-like',
            'frequency': 'Variável',
            'prognosis': 'Geralmente desfavorável sem TKI',
            'risk_group': 'Adverso',
            'therapy': 'Considerar TKI dirigido + quimioterapia',
            'monitoring': 'MRD',
            'clinical_significance': 'Ph-like signature - pode se beneficiar de TKI',
            'references': 'PMID: 27895058'
        }

    # Desconhecida/rara
    return {
        'leukemia_type': 'Não classificada / Rara',
        'classification': 'Requer investigação adicional',
        'frequency': 'Rara ou não reportada',
        'prognosis': 'Desconhecido',
        'risk_group': 'Indeterminado',
        'therapy': 'Protocolo padrão baseado em tipo celular',
        'monitoring': 'MRD',
        'clinical_significance': 'FUSÃO NOVA/RARA - Investigar literatura',
        'references': 'N/A'
    }


def process_arriba_output(input_file, output_file):
    """
    Processa output do Arriba e adiciona anotações clínicas
    """
    logger.info(f"Processando: {input_file}")

    # Lê arquivo Arriba
    df = pd.read_csv(input_file, sep='\t')

    # Remove '#' do header se existir
    df.columns = df.columns.str.replace('^#', '', regex=True)

    logger.info(f"  Total fusões: {len(df)}")

    # Filtra baixa confiança
    if 'confidence' in df.columns:
        df_filtered = df[df['confidence'] != 'low'].copy()
        logger.info(f"  Após filtrar 'low confidence': {len(df_filtered)}")
    else:
        df_filtered = df.copy()

    # Adiciona anotações
    annotations = []
    for _, row in df_filtered.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']

        annot = annotate_fusion(gene1, gene2)
        annotations.append(annot)

    # Adiciona colunas de anotação
    for key in ['leukemia_type', 'classification', 'frequency', 'prognosis',
                'risk_group', 'therapy', 'monitoring', 'clinical_significance', 'references']:
        df_filtered[f'annot_{key}'] = [a.get(key, 'N/A') for a in annotations]

    # Cria coluna fusion padronizada
    df_filtered['fusion'] = df_filtered.apply(
        lambda row: f"{row['gene1'].split('(')[0].upper()}--{row['gene2'].split('(')[0].upper()}",
        axis=1
    )

    # Reordena colunas para ter anotações no início
    annot_cols = [c for c in df_filtered.columns if c.startswith('annot_')]
    other_cols = [c for c in df_filtered.columns if not c.startswith('annot_')]

    # Ordem: fusion, type, confidence, anotações, resto
    priority_cols = ['fusion', 'gene1', 'gene2', 'type', 'confidence']
    priority_cols = [c for c in priority_cols if c in df_filtered.columns]

    remaining_cols = [c for c in other_cols if c not in priority_cols]

    df_filtered = df_filtered[priority_cols + annot_cols + remaining_cols]

    # Salva
    df_filtered.to_csv(output_file, sep='\t', index=False)

    logger.info(f"✓ Anotações salvas: {output_file}")

    # Estatísticas
    logger.info("\n=== RESUMO DAS FUSÕES ===")

    # NOTA: usa annot_classification (começa com "WHO:" em toda entrada com
    # entidade reconhecida), não a string 'PATOGNOMÔNICA' em
    # annot_clinical_significance -- esse texto só foi digitado manualmente
    # em poucas entradas do database (BCR-ABL1, PML-RARA, RUNX1-RUNX1T1,
    # CBFB-MYH11), deixando de fora ETV6-RUNX1, TCF3-PBX1 e outras fusões
    # igualmente definidoras de entidade WHO-HAEM5. "WHO:" é o sinal mais
    # completo e confiável para "fusão definidora de entidade" (usado também
    # por classify_leukemia.py).
    defining = df_filtered[df_filtered['annot_classification'].str.startswith('WHO:', na=False)]
    logger.info(f"Fusões definidoras de entidade WHO-HAEM5: {len(defining)}")

    for _, row in defining.iterrows():
        logger.info(f"  • {row['fusion']}: {row['annot_leukemia_type']} ({row['annot_prognosis']})")

    known = df_filtered[~df_filtered['annot_leukemia_type'].isin(['Não classificada / Rara', 'Possível B-ALL Ph-like'])]
    logger.info(f"\nFusões recorrentes conhecidas: {len(known)}")

    novel = df_filtered[df_filtered['annot_leukemia_type'] == 'Não classificada / Rara']
    logger.info(f"Fusões novas/raras: {len(novel)}")

    if len(novel) > 0:
        logger.info("  Fusões para investigar:")
        for _, row in novel.iterrows():
            logger.info(f"    • {row['fusion']} ({row.get('reading_frame', 'N/A')} frame)")

    return df_filtered


def main():
    parser = argparse.ArgumentParser(
        description='Anota fusões do Arriba com significância clínica em leucemias'
    )
    parser.add_argument('--input', required=True, help='Arquivo Arriba TSV de entrada')
    parser.add_argument('--output', required=True, help='Arquivo TSV anotado de saída')

    args = parser.parse_args()

    df_annotated = process_arriba_output(args.input, args.output)

    logger.info("\n" + "="*60)
    logger.info("ANOTAÇÃO COMPLETA!")
    logger.info("="*60)
    logger.info(f"\nOutput: {args.output}")
    logger.info("\nColunas adicionadas:")
    logger.info("  • annot_leukemia_type: Tipo de leucemia associado")
    logger.info("  • annot_classification: Classificação WHO")
    logger.info("  • annot_prognosis: Prognóstico")
    logger.info("  • annot_risk_group: Grupo de risco (ELN 2022)")
    logger.info("  • annot_therapy: Terapia recomendada")
    logger.info("  • annot_monitoring: Monitoramento sugerido")
    logger.info("  • annot_clinical_significance: Significância clínica")


if __name__ == '__main__':
    main()
