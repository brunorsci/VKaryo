### VKARYO_DB ###

library(DBI)
library(RPostgres)

library(DBI)
library(RPostgres)

# 1️⃣ Conectar como superusuário (por exemplo: postgres)
#con_admin <- dbConnect(
#  RPostgres::Postgres(),
#  dbname = "postgres",
#  host = "10.23.5.102",
#  port = 5432,
#  user = "postgres",        # admin
#  password = "burloga1"
#)

# 2️⃣ Criar usuário e banco
#dbExecute(con_admin, "CREATE USER vkaryo WITH PASSWORD 'darwin58';")
#dbExecute(con_admin, "CREATE DATABASE \"VKaryo_DB\" OWNER vkaryo;")
#dbExecute(con_admin, "GRANT ALL PRIVILEGES ON DATABASE \"VKaryo_DB\" TO vkaryo;")

#dbDisconnect(con_admin)

library(DBI)
library(RPostgres)

con_vkaryo <- dbConnect(
  RPostgres::Postgres(),
  dbname = "VKaryo_DB",
  host = "10.23.5.102",
  port = 5432,
  user = "vkaryo",
  password = "darwin58"
)

# ARRIBA TABLE
dbExecute(con_vkaryo, "
CREATE TABLE IF NOT EXISTS arriba_fusions (
  samples_id TEXT,
  fusion TEXT,
  type TEXT,
  confidence TEXT,
  gene1 TEXT,
  gene2 TEXT,
  gene_id1 TEXT,
  gene_id2 TEXT,
  breakpoint1 TEXT,
  breakpoint2 TEXT,
  reading_frame TEXT,
  coverage1 INTEGER,
  coverage2 INTEGER,
  site1 TEXT,
  site2 TEXT,
  direction1 TEXT,
  direction2 TEXT,
  retained_protein_domains TEXT,
  fusion_transcript TEXT,
  peptide_sequence TEXT
);
")
# subindo os dados para o db
dbWriteTable(con_vkaryo, "arriba_fusions", Arriba_clean, append = TRUE, row.names = FALSE)

## tabela annotsv delly#####

#create table
dbExecute(con_vkaryo, "
CREATE TABLE IF NOT EXISTS delly_sv (
  id TEXT,
  sample TEXT,
  hgvs_notation TEXT,
  svnote TEXT,
  chrom TEXT,
  start_pos BIGINT,
  end_pos BIGINT,
  length BIGINT,
  type TEXT,
  cyto_band TEXT,
  gene_name TEXT,
  gene_count INTEGER,
  frameshift TEXT,
  exon_count INTEGER,
  gen_cc_disease TEXT,
  gen_cc_classification TEXT,
  omim_phenotype TEXT,
  omim_inheritance TEXT,
  p_gain_phen TEXT,
  p_loss_phen TEXT,
  p_snvindel_phen TEXT
);
")

# filing table
dbWriteTable(con_vkaryo, "delly_sv", Annotsv_clean, append = TRUE, row.names = FALSE)


dbGetQuery(con_vkaryo, "SELECT * FROM delly_sv LIMIT 10;")


