#!/usr/bin/env python3
"""
Popula banco de dados VKaryo com resultados de análise
Consolida todos os outputs em um database SQLite centralizado
"""

import sqlite3
import json
import pandas as pd
import gzip
from pathlib import Path
from datetime import datetime
import re

# Parâmetros do Snakemake
sample_id = snakemake.wildcards.sample
db_path = snakemake.output.db
config = snakemake.config

# Inputs (todos os resultados)
inputs = snakemake.input

def init_database(db_path):
    """Inicializa database com schema"""
    conn = sqlite3.connect(db_path)

    # Carrega schema
    schema_path = Path(__file__).parent / "database_schema.sql"
    with open(schema_path, 'r') as f:
        schema = f.read()

    conn.executescript(schema)
    conn.commit()
    return conn

def insert_sample(conn, sample_id, sample_type="DNA"):
    """Insere ou atualiza registro de amostra"""
    cursor = conn.cursor()
    cursor.execute("""
        INSERT OR REPLACE INTO samples (sample_id, sample_type, analysis_date)
        VALUES (?, ?, ?)
    """, (sample_id, sample_type, datetime.now()))
    conn.commit()

def insert_purecn_results(conn, sample_id, purecn_summary):
    """Insere resultados PureCN"""
    if not Path(purecn_summary).exists():
        return

    # Parse sumário PureCN
    data = {}
    with open(purecn_summary, 'r') as f:
        content = f.read()

        # Extrai valores
        purity_match = re.search(r'Purity: ([\d.]+)%', content)
        ploidy_match = re.search(r'Ploidy: ([\d.]+)', content)
        sex_match = re.search(r'Sex: (\w+)', content)
        contam_match = re.search(r'Contamination: ([\d.]+)', content)
        amp_match = re.search(r'Amplifications.*?:\s*(\d+)', content)
        del_match = re.search(r'Deletions.*?:\s*(\d+)', content)
        loh_match = re.search(r'LOH events:\s*(\d+)', content)
        ll_match = re.search(r'Log-likelihood: ([-\d.]+)', content)

        data = {
            'purity': float(purity_match.group(1)) / 100 if purity_match else None,
            'ploidy': float(ploidy_match.group(1)) if ploidy_match else None,
            'sex': sex_match.group(1) if sex_match else None,
            'contamination': float(contam_match.group(1)) if contam_match else None,
            'n_amplifications': int(amp_match.group(1)) if amp_match else 0,
            'n_deletions': int(del_match.group(1)) if del_match else 0,
            'n_loh': int(loh_match.group(1)) if loh_match else 0,
            'log_likelihood': float(ll_match.group(1)) if ll_match else None
        }

    cursor = conn.cursor()
    cursor.execute("""
        INSERT OR REPLACE INTO purecn_results
        (sample_id, purity, ploidy, sex, contamination, n_amplifications, n_deletions, n_loh, log_likelihood)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (sample_id, data['purity'], data['ploidy'], data['sex'], data['contamination'],
          data['n_amplifications'], data['n_deletions'], data['n_loh'], data['log_likelihood']))
    conn.commit()

def insert_optitype_results(conn, sample_id, hla_typing):
    """Insere resultados OptiType"""
    if not Path(hla_typing).exists():
        return

    with open(hla_typing, 'r') as f:
        content = f.read()

        # Extrai alelos HLA
        hla_a = re.search(r'HLA-A:\s*([^,]+),\s*([^\n]+)', content)
        hla_b = re.search(r'HLA-B:\s*([^,]+),\s*([^\n]+)', content)
        hla_c = re.search(r'HLA-C:\s*([^,]+),\s*([^\n]+)', content)
        reads = re.search(r'Reads:\s*(\d+)', content)
        objective = re.search(r'Objective:\s*([\d.]+)', content)

        cursor = conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO optitype_results
            (sample_id, hla_a1, hla_a2, hla_b1, hla_b2, hla_c1, hla_c2, reads, objective)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (sample_id,
              hla_a.group(1).strip() if hla_a else None,
              hla_a.group(2).strip() if hla_a else None,
              hla_b.group(1).strip() if hla_b else None,
              hla_b.group(2).strip() if hla_b else None,
              hla_c.group(1).strip() if hla_c else None,
              hla_c.group(2).strip() if hla_c else None,
              int(reads.group(1)) if reads else None,
              float(objective.group(1)) if objective else None))
        conn.commit()

def insert_mrd_results(conn, sample_id, mrd_summary, mrd_tracking):
    """Insere resultados MRD"""
    if not Path(mrd_summary).exists():
        return

    with open(mrd_summary, 'r') as f:
        content = f.read()

        total_var = re.search(r'Total variants called:\s*(\d+)', content)
        mrd_var = re.search(r'MRD-relevant variants:\s*(\d+)', content)
        status = re.search(r'MRD Status:\s*(\w+)', content)
        max_vaf = re.search(r'Highest VAF detected:\s*([\d.]+)%', content)

        cursor = conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO mrd_results
            (sample_id, total_variants, mrd_variants, mrd_status, max_vaf)
            VALUES (?, ?, ?, ?, ?)
        """, (sample_id,
              int(total_var.group(1)) if total_var else 0,
              int(mrd_var.group(1)) if mrd_var else 0,
              status.group(1) if status else "UNKNOWN",
              float(max_vaf.group(1)) if max_vaf else 0.0))
        conn.commit()

    # Insere variantes MRD individuais
    if Path(mrd_tracking).exists():
        try:
            df = pd.read_csv(mrd_tracking, sep='\t')
            for _, row in df.iterrows():
                cursor.execute("""
                    INSERT INTO mrd_variants
                    (sample_id, gene, chromosome, position, ref, alt, vaf, cosmic_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (sample_id, row.get('Gene'), row.get('Position', '').split(':')[0] if ':' in str(row.get('Position', '')) else None,
                      int(row.get('Position', '').split(':')[1]) if ':' in str(row.get('Position', '')) else None,
                      row.get('Alleles', '').split('>')[0] if '>' in str(row.get('Alleles', '')) else None,
                      row.get('Alleles', '').split('>')[1] if '>' in str(row.get('Alleles', '')) else None,
                      row.get('VAF (%)'), row.get('COSMIC')))
            conn.commit()
        except Exception as e:
            print(f"Warning: Could not parse MRD tracking: {e}")

def insert_driver_mutations(conn, sample_id, driver_file):
    """Insere driver mutations do Mutect2"""
    if not Path(driver_file).exists():
        return

    try:
        df = pd.read_csv(driver_file, sep='\t')
        cursor = conn.cursor()

        for _, row in df.iterrows():
            cursor.execute("""
                INSERT INTO driver_mutations
                (sample_id, gene, chromosome, position, ref, alt, impact, exonic_func, aa_change, cosmic_id, clinvar)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (sample_id, row.get('Gene'), row.get('Chr'), row.get('Start'),
                  row.get('Ref'), row.get('Alt'), row.get('Impact'),
                  row.get('ExonicFunc'), row.get('AAChange'),
                  row.get('COSMIC'), row.get('ClinVar')))
        conn.commit()
    except Exception as e:
        print(f"Warning: Could not parse driver mutations: {e}")

def insert_cnv_calls(conn, sample_id, cnv_file, validated_file=None):
    """Insere CNVs do CNVkit"""
    if not Path(cnv_file).exists():
        return

    try:
        df = pd.read_csv(cnv_file, sep='\t')

        # Carrega validação se disponível
        validated = set()
        if validated_file and Path(validated_file).exists():
            val_df = pd.read_csv(validated_file, sep='\t')
            if 'validation_status' in val_df.columns:
                validated = set(val_df[val_df['validation_status'] == 'validated'].index)

        cursor = conn.cursor()
        for idx, row in df.iterrows():
            # Determina status do CN
            log2 = row.get('log2', 0)
            if log2 > 0.3:
                cn_status = 'gain'
            elif log2 < -0.3:
                cn_status = 'loss'
            else:
                cn_status = 'neutral'

            cursor.execute("""
                INSERT INTO cnv_calls
                (sample_id, chromosome, start, end, gene, log2_ratio, depth, weight, cn_status, validated)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (sample_id, row.get('chromosome'), row.get('start'), row.get('end'),
                  row.get('gene'), log2, row.get('depth'), row.get('weight'),
                  cn_status, idx in validated))
        conn.commit()
    except Exception as e:
        print(f"Warning: Could not parse CNVs: {e}")

def insert_gene_fusions(conn, sample_id, fusion_file):
    """Insere fusões de genes do Arriba (anotado)"""
    if not Path(fusion_file).exists():
        return

    try:
        df = pd.read_csv(fusion_file, sep='\t')
        cursor = conn.cursor()

        for _, row in df.iterrows():
            cursor.execute("""
                INSERT INTO gene_fusions
                (sample_id, gene1, gene2, breakpoint1, breakpoint2, type, split_reads, discordant_pairs,
                 confidence, pathognomonic, prognosis, therapy)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (sample_id, row.get('#gene1'), row.get('gene2'),
                  row.get('breakpoint1'), row.get('breakpoint2'),
                  row.get('type'), row.get('split_reads1', 0) + row.get('split_reads2', 0),
                  row.get('discordant_mates'), row.get('confidence'),
                  row.get('pathognomonic', False), row.get('prognosis'),
                  row.get('therapy')))
        conn.commit()
    except Exception as e:
        print(f"Warning: Could not parse fusions: {e}")

def insert_classification(conn, sample_id, classification_file):
    """Insere classificação de leucemia"""
    if not Path(classification_file).exists():
        return

    try:
        with open(classification_file, 'r') as f:
            data = json.load(f)

        cursor = conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO leukemia_classification
            (sample_id, leukemia_type, subtype, risk_category, eln_2022_risk, confidence, key_features, recommended_therapy)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (sample_id, data.get('leukemia_type'), data.get('subtype'),
              data.get('risk_category'), data.get('eln_2022_risk'),
              data.get('confidence'), json.dumps(data.get('key_features', {})),
              data.get('recommended_therapy')))
        conn.commit()
    except Exception as e:
        print(f"Warning: Could not parse classification: {e}")

# Inicializa database
conn = init_database(db_path)

# Insere amostra
insert_sample(conn, sample_id)

# Insere resultados de cada análise
try:
    if hasattr(inputs, 'purecn_summary'):
        insert_purecn_results(conn, sample_id, inputs.purecn_summary)
except: pass

try:
    if hasattr(inputs, 'hla_typing'):
        insert_optitype_results(conn, sample_id, inputs.hla_typing)
except: pass

try:
    if hasattr(inputs, 'mrd_summary') and hasattr(inputs, 'mrd_tracking'):
        insert_mrd_results(conn, sample_id, inputs.mrd_summary, inputs.mrd_tracking)
except: pass

try:
    if hasattr(inputs, 'driver_mutations'):
        insert_driver_mutations(conn, sample_id, inputs.driver_mutations)
except: pass

try:
    validated = inputs.cnv_validated if hasattr(inputs, 'cnv_validated') else None
    if hasattr(inputs, 'cnv_calls'):
        insert_cnv_calls(conn, sample_id, inputs.cnv_calls, validated)
except: pass

try:
    if hasattr(inputs, 'fusions_annotated'):
        insert_gene_fusions(conn, sample_id, inputs.fusions_annotated)
except: pass

try:
    if hasattr(inputs, 'classification'):
        insert_classification(conn, sample_id, inputs.classification)
except: pass

conn.close()

print(f"✓ Database populated for sample {sample_id}")
print(f"✓ Database location: {db_path}")
