#!/usr/bin/env python3
"""
Consolida databases individuais de amostras em um único database central
"""

import sqlite3
import shutil
from pathlib import Path

# Parâmetros do Snakemake
input_dbs = snakemake.input.dbs
output_db = snakemake.output.merged_db

# Cria database central com schema
schema_path = Path(__file__).parent / "database_schema.sql"

# Inicializa database central
conn_central = sqlite3.connect(output_db)
with open(schema_path, 'r') as f:
    schema = f.read()
conn_central.executescript(schema)
conn_central.commit()

# Tabelas para copiar
TABLES = [
    'samples',
    'quality_metrics',
    'purecn_results',
    'optitype_results',
    'mrd_results',
    'mrd_variants',
    'driver_mutations',
    'cnv_calls',
    'gene_fusions',
    'snv_calls',
    'aneuploidy_events',
    'leukemia_classification',
    'therapy_recommendations'
]

# Merge cada database individual
for db_path in input_dbs:
    if not Path(db_path).exists():
        print(f"Warning: {db_path} não encontrado, pulando...")
        continue

    print(f"Merging {db_path}...")

    # Conecta ao database individual
    conn_individual = sqlite3.connect(db_path)
    cursor_individual = conn_individual.cursor()

    # Copia dados de cada tabela
    for table in TABLES:
        try:
            # Verifica se tabela existe
            cursor_individual.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table}'")
            if not cursor_individual.fetchone():
                continue

            # Pega schema da tabela
            cursor_individual.execute(f"PRAGMA table_info({table})")
            columns_info = cursor_individual.fetchall()
            columns = [col[1] for col in columns_info]

            # Seleciona todos os dados
            cursor_individual.execute(f"SELECT * FROM {table}")
            rows = cursor_individual.fetchall()

            if rows:
                # Insere no database central
                placeholders = ','.join(['?'] * len(columns))
                insert_query = f"INSERT OR REPLACE INTO {table} VALUES ({placeholders})"

                cursor_central = conn_central.cursor()
                cursor_central.executemany(insert_query, rows)
                conn_central.commit()

                print(f"  ✓ {table}: {len(rows)} registros")

        except Exception as e:
            print(f"  ⚠ Erro ao copiar {table}: {e}")

    conn_individual.close()

conn_central.close()

print(f"\n✓ Database central criado: {output_db}")
print(f"✓ Total de amostras consolidadas: {len(input_dbs)}")
