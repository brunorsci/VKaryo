-- VKaryo Database Schema
-- Banco de dados centralizado para todos os resultados de cariotipagem virtual

-- Tabela de amostras
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,
    analysis_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    sample_type TEXT,  -- DNA, RNA, ou BOTH
    disease_type TEXT,  -- AML, ALL, CML, CLL, etc.
    notes TEXT
);

-- Tabela de métricas de qualidade
CREATE TABLE IF NOT EXISTS quality_metrics (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    total_reads BIGINT,
    mapped_reads BIGINT,
    duplicate_rate REAL,
    mean_coverage REAL,
    coverage_10x REAL,
    coverage_30x REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela PureCN (pureza e ploidy)
CREATE TABLE IF NOT EXISTS purecn_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    purity REAL,
    ploidy REAL,
    sex TEXT,
    contamination REAL,
    n_amplifications INTEGER,
    n_deletions INTEGER,
    n_loh INTEGER,
    log_likelihood REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela OptiType (HLA typing)
CREATE TABLE IF NOT EXISTS optitype_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    hla_a1 TEXT,
    hla_a2 TEXT,
    hla_b1 TEXT,
    hla_b2 TEXT,
    hla_c1 TEXT,
    hla_c2 TEXT,
    reads INTEGER,
    objective REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela MRD (VarDict)
CREATE TABLE IF NOT EXISTS mrd_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    total_variants INTEGER,
    mrd_variants INTEGER,
    mrd_status TEXT,  -- NEGATIVE, LOW, MODERATE, HIGH
    max_vaf REAL,
    analysis_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de variantes MRD
CREATE TABLE IF NOT EXISTS mrd_variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    gene TEXT,
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    vaf REAL,
    depth INTEGER,
    cosmic_id TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela Mutect2 driver mutations
CREATE TABLE IF NOT EXISTS driver_mutations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    gene TEXT,
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    vaf REAL,
    impact TEXT,  -- HIGH, MODERATE, LOW
    exonic_func TEXT,
    aa_change TEXT,
    cosmic_id TEXT,
    clinvar TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de CNVs (CNVkit)
CREATE TABLE IF NOT EXISTS cnv_calls (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    chromosome TEXT,
    start INTEGER,
    end INTEGER,
    gene TEXT,
    log2_ratio REAL,
    depth REAL,
    weight REAL,
    cn INTEGER,  -- Copy number
    cn_status TEXT,  -- gain, loss, neutral
    validated BOOLEAN DEFAULT 0,  -- Validado com RNA
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de fusões (Arriba)
CREATE TABLE IF NOT EXISTS gene_fusions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    gene1 TEXT,
    gene2 TEXT,
    breakpoint1 TEXT,
    breakpoint2 TEXT,
    type TEXT,
    split_reads INTEGER,
    discordant_pairs INTEGER,
    confidence TEXT,  -- high, medium, low
    pathognomonic BOOLEAN,
    prognosis TEXT,  -- favorable, unfavorable, intermediate
    therapy TEXT,  -- Terapia sugerida
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de SNVs (Strelka/Mutect2)
CREATE TABLE IF NOT EXISTS snv_calls (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    caller TEXT,  -- strelka ou mutect2
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    vaf REAL,
    depth INTEGER,
    quality REAL,
    gene TEXT,
    exonic_func TEXT,
    aa_change TEXT,
    cosmic_id TEXT,
    clinvar TEXT,
    gnomad_af REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de aneuploidias
CREATE TABLE IF NOT EXISTS aneuploidy_events (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    chromosome TEXT,
    event_type TEXT,  -- trisomy, monosomy, tetrasomy, partial_gain, partial_loss
    mean_log2 REAL,
    size_mb REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de classificação de leucemia
CREATE TABLE IF NOT EXISTS leukemia_classification (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    leukemia_type TEXT,  -- AML, ALL, CML, CLL, MDS
    subtype TEXT,
    risk_category TEXT,  -- favorable, intermediate, adverse
    eln_2022_risk TEXT,  -- Para AML
    confidence REAL,
    key_features TEXT,  -- JSON com features principais
    recommended_therapy TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Tabela de terapias recomendadas
CREATE TABLE IF NOT EXISTS therapy_recommendations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    therapy_class TEXT,  -- chemotherapy, targeted, immunotherapy
    drug_name TEXT,
    indication TEXT,
    evidence_level TEXT,  -- FDA_approved, clinical_trial, preclinical
    based_on TEXT,  -- Qual alteração justifica (ex: FLT3-ITD)
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Índices para performance
CREATE INDEX IF NOT EXISTS idx_mrd_variants_sample ON mrd_variants(sample_id);
CREATE INDEX IF NOT EXISTS idx_mrd_variants_gene ON mrd_variants(gene);
CREATE INDEX IF NOT EXISTS idx_driver_mutations_sample ON driver_mutations(sample_id);
CREATE INDEX IF NOT EXISTS idx_driver_mutations_gene ON driver_mutations(gene);
CREATE INDEX IF NOT EXISTS idx_cnv_sample ON cnv_calls(sample_id);
CREATE INDEX IF NOT EXISTS idx_cnv_gene ON cnv_calls(gene);
CREATE INDEX IF NOT EXISTS idx_fusions_sample ON gene_fusions(sample_id);
CREATE INDEX IF NOT EXISTS idx_snv_sample ON snv_calls(sample_id);
CREATE INDEX IF NOT EXISTS idx_snv_gene ON snv_calls(gene);

-- Views úteis
CREATE VIEW IF NOT EXISTS v_sample_summary AS
SELECT
    s.sample_id,
    s.analysis_date,
    p.purity,
    p.ploidy,
    m.mrd_status,
    m.max_vaf as mrd_max_vaf,
    l.leukemia_type,
    l.risk_category,
    COUNT(DISTINCT d.gene) as n_driver_genes,
    COUNT(DISTINCT f.id) as n_fusions,
    COUNT(DISTINCT c.id) as n_cnvs
FROM samples s
LEFT JOIN purecn_results p ON s.sample_id = p.sample_id
LEFT JOIN mrd_results m ON s.sample_id = m.sample_id
LEFT JOIN leukemia_classification l ON s.sample_id = l.sample_id
LEFT JOIN driver_mutations d ON s.sample_id = d.sample_id
LEFT JOIN gene_fusions f ON s.sample_id = f.sample_id
LEFT JOIN cnv_calls c ON s.sample_id = c.sample_id
GROUP BY s.sample_id;

CREATE VIEW IF NOT EXISTS v_actionable_mutations AS
SELECT
    d.sample_id,
    d.gene,
    d.aa_change,
    d.vaf,
    d.impact,
    t.drug_name,
    t.indication,
    t.evidence_level
FROM driver_mutations d
LEFT JOIN therapy_recommendations t ON d.sample_id = t.sample_id
    AND d.gene = t.based_on
WHERE d.impact IN ('HIGH', 'MODERATE')
ORDER BY d.vaf DESC;
