#!/usr/bin/env python3
"""Resolve a entidade WHO-HAEM5/ICC combinando a alteracao definidora com a linhagem medida.

Varias alteracoes definidoras ocorrem em mais de uma entidade, e a entidade so
fica determinada quando se sabe de que celula se trata. A t(11;14) define linfoma
de celulas do manto em celula B madura e mieloma com t(11;14) em celula
plasmocitaria; um rearranjo de MYC junto a um locus de imunoglobulina define
Burkitt em celula B madura e mieloma com rearranjo de MYC em plasmocito; um
rearranjo de KMT2A define LMA em celula mieloide e LLA-B em precursor B.

Atribuir a entidade so pela fusao produz exatamente os erros observados na coorte
observados na pratica: mielomas reportados como linfoma de celulas do manto ou como Burkitt, e
linfoma B maduro reportado como "LMA/LLA com KMT2A".

A linhagem vem do arquivo de referencia linhagem_celular_por_amostra.tsv, que
registra o imunofenotipo de cada amostra -- em uso clinico, o resultado da
citometria de fluxo, disponivel de rotina em sangue periferico e aspirado de
medula e especifico do clone por ser lido com gate nos blastos. Nao se infere a
linhagem a partir da expressao porque a amostra de paciente e uma mistura
celular: os marcadores refletiriam a composicao do aspirado, e nao o clone.

Quando a linhagem e desconhecida, a entidade e mantida como estava e recebe
ressalva, em vez de ser descartada.
"""

import argparse
import csv
import io
import re

# padrao da alteracao -> {familia ou painel de linhagem: entidade resolvida}
REGRAS = [
    (re.compile(r"\b(CCND1|MYEOV)\b", re.I), {
        "plasmocitaria": "Mieloma múltiplo com t(11;14)(q13;q32)/IGH::CCND1",
        "B_maduro": "Linfoma de Células do Manto (LCM)",
    }),
    (re.compile(r"\b(MYC|PVT1)\b", re.I), {
        "plasmocitaria": "Mieloma múltiplo com rearranjo de MYC",
        "B_maduro": "Linfoma de Burkitt (se rearranjo único) ou componente de linfoma B de alto grau",
    }),
    (re.compile(r"\bKMT2A\b", re.I), {
        "mieloide": "LMA com rearranjo de KMT2A",
        "B_precursor": "LLA-B com rearranjo de KMT2A",
        "B_maduro": "Neoplasia de células B maduras com rearranjo de KMT2A (apresentação atípica)",
        "T": "LLA-T com rearranjo de KMT2A",
    }),
]

# Entidades cuja linhagem implicada permite checar coerencia sem reescrever a atribuicao.
LINHAGEM_IMPLICADA = {
    "mieloide": re.compile(r"\b(LMA|AML|APL|promielo|eosinof|mieloide)\b", re.I),
    "linfoide_B": re.compile(r"\b(B-ALL|LLA-B|Pre-B|Burkitt|manto|DLBCL|folicular|mieloma|linfoma B)\b", re.I),
    "linfoide_T": re.compile(r"\b(T-ALL|LLA-T|ALCL|anaplásico|células T)\b", re.I),
}

FAMILIA = {"B_maduro": "linfoide_B", "B_precursor": "linfoide_B", "plasmocitaria": "linfoide_B", "T": "linfoide_T",
           "NK": "NK", "mieloide": "mieloide", "eritroide": "mieloide",
           "megacariocitica": "mieloide"}

COMPATIVEL = {"mieloide": {"mieloide"}, "linfoide_B": {"linfoide_B"},
              "linfoide_T": {"linfoide_T", "NK"}, "NK": {"NK", "linfoide_T"}}


# Uma atribuicao so e refinada pela linhagem quando ela propria e ambigua, isto e,
# quando nomeia mais de uma linhagem ou mais de uma entidade alternativa. Uma
# atribuicao que ja nomeia a translocacao especifica -- "LMA com t(9;11)" ou
# "LLA-B com t(4;11)" -- e mais informativa do que a forma generica por gene e
# nunca e substituida por ela; nesses casos resta apenas conferir a coerencia.
AMBIGUA = re.compile(r"/|\bou\b|\bse rearranjo", re.I)


def e_ambigua(entidade):
    """Verdadeiro quando a entidade oferece alternativas em vez de uma atribuicao unica."""
    sem_citogenetica = re.sub(r"\([^)]*\)|t\(\d+;\d+\)|inv\(\d+\)", "", entidade)
    return bool(AMBIGUA.search(sem_citogenetica))


def _preservar_citogenetica(alvo, original):
    """Mantem a translocacao nomeada na atribuicao original, que e mais informativa."""
    achado = re.search(r"(t\(\d+;\d+\)\([^)]*\)|inv\(\d+\)\([^)]*\))", original)
    if achado and achado.group(1) not in alvo:
        return re.sub(r"rearranjo de (\w+)", f"\\g<0> {achado.group(1)}", alvo, count=1)
    return alvo


def resolver(entidade, evidencia, painel, familia):
    """Devolve (entidade, ressalva). Nao descarta atribuicao: corrige ou sinaliza."""
    if not entidade or entidade in ("—", "-", ""):
        return entidade, ""
    if familia in ("indeterminada", "ausente", None):
        return entidade, ("linhagem celular indeterminada pela expressão de marcadores; "
                          "entidade mantida como atribuída pela alteração, sem verificação de coerência")

    for padrao, por_linhagem in REGRAS:
        if not padrao.search(evidencia or ""):
            continue
        alvo = por_linhagem.get(painel) or por_linhagem.get(familia)
        if not alvo:
            continue
        if alvo == entidade:
            return entidade, ""
        # Refina quando a atribuicao e ambigua, e corrige quando ela e justamente
        # a alternativa da mesma regra para outra linhagem -- caso do mieloma
        # reportado como linfoma de celulas do manto, em que as duas entidades
        # sao da mesma familia e so o estagio de maturacao as separa.
        alternativa_de_outra_linhagem = entidade in por_linhagem.values()
        if e_ambigua(entidade) or alternativa_de_outra_linhagem:
            return _preservar_citogenetica(alvo, entidade), (
                f"entidade resolvida pela linhagem medida ({painel}); "
                f"a alteração também ocorre em: {entidade}")
        return entidade, ""

    for esperada, padrao in LINHAGEM_IMPLICADA.items():
        if padrao.search(entidade):
            if familia not in COMPATIVEL.get(esperada, {esperada}):
                return entidade, (f"atribuição incoerente com a linhagem medida: a entidade implica "
                                  f"{esperada} e a expressão de marcadores indica {familia}")
            break
    return entidade, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tabela", required=True, help="CSV de classificacao (Tabela S3)")
    ap.add_argument("--linhagens", required=True,
                    help="reference/linhagem_celular_por_amostra.tsv")
    ap.add_argument("--saida", required=True)
    ap.add_argument("--col-entidade", default="Entidade_inferida_WHO_HAEM5_ICC")
    ap.add_argument("--col-evidencia", default="Evidencia_que_definiu_a_entidade")
    ap.add_argument("--col-ressalva", default="Ressalva_da_atribuicao")
    ap.add_argument("--col-linhagem", default="Linhagem_celular")
    args = ap.parse_args()

    linhagens = {}
    rna = {}
    for r in csv.DictReader(io.open(args.linhagens, encoding="utf-8"), delimiter="\t"):
        linhagens[r["amostra_RNA"]] = r
        rna[r["linhagem"]] = r["amostra_RNA"]

    linhas = list(csv.DictReader(io.open(args.tabela, encoding="utf-8")))
    campos = list(linhas[0].keys())
    if "Linhagem_celular_imunofenotipo" not in campos:
        campos.insert(campos.index(args.col_entidade), "Linhagem_celular_imunofenotipo")

    mudadas = 0
    sinalizadas = 0
    for r in linhas:
        info = linhagens.get(rna.get(r[args.col_linhagem], ""), {})
        painel = info.get("estagio_de_linhagem", "-")
        familia = info.get("familia_de_linhagem", "ausente")
        r["Linhagem_celular_imunofenotipo"] = familia if familia != "ausente" else "-"
        antes = r[args.col_entidade]
        nova, ressalva = resolver(antes, r.get(args.col_evidencia, ""), painel, familia)
        if nova != antes:
            mudadas += 1
            print(f"  {r[args.col_linhagem]:<13} {antes}\n  {'':<13} -> {nova}")
        if ressalva:
            sinalizadas += 1
            atual = (r.get(args.col_ressalva) or "").strip()
            if atual in ("—", "-", ""):
                r[args.col_ressalva] = ressalva
            elif ressalva not in atual:
                r[args.col_ressalva] = f"{atual}; {ressalva}"
        r[args.col_entidade] = nova

    with io.open(args.saida, "w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, fieldnames=campos, lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    print(f"\nentidades corrigidas: {mudadas}   com ressalva: {sinalizadas}")


if __name__ == "__main__":
    main()
