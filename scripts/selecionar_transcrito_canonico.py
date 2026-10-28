#!/usr/bin/env python3
"""Escolhe o transcrito canonico entre os reportados pelo ANNOVAR em AAChange.refGene.

O ANNOVAR lista todos os transcritos do gene e as ferramentas a jusante costumam
tomar o primeiro da lista, que nao e o canonico: para TP53, o primeiro e o
NM_001126115 (261 aa), cuja numeracao desloca os hotspots em 132 residuos e faz
a troca M133K aparecer como perda de inicio (p.M1?).

A escolha aqui e o transcrito de maior CDS, com desempate pelo menor numero de
acesso RefSeq -- criterio que favorece sistematicamente os acessos originais
(NM_000546) sobre as isoformas descritas depois (NM_001126112).

O comprimento de CDS vem do proprio GTF do refGene ja usado pelo workflow, de
modo que a regra nao depende de arquivo externo nem de curadoria manual.
"""

import argparse
import gzip
import re
from collections import defaultdict

ACESSO = re.compile(r"(?:NM|NR|XM|XR)_(\d+)")


def comprimentos_de_cds(caminho_gtf):
    """Soma as bases de CDS de cada transcrito do GTF do refGene."""
    por_transcrito = defaultdict(int)
    abrir = gzip.open if caminho_gtf.endswith(".gz") else open
    with abrir(caminho_gtf, "rt") as gtf:
        for linha in gtf:
            if linha.startswith("#"):
                continue
            campos = linha.rstrip("\n").split("\t")
            if len(campos) < 9 or campos[2] != "CDS":
                continue
            achado = re.search(r'transcript_id "([^"]+)"', campos[8])
            if not achado:
                continue
            por_transcrito[achado.group(1)] += int(campos[4]) - int(campos[3]) + 1
    return por_transcrito


def numero_de_acesso(transcrito):
    achado = ACESSO.search(transcrito)
    return int(achado.group(1)) if achado else float("inf")


def escolher(aachange, comprimentos):
    """Devolve a anotacao do transcrito canonico, ou a primeira se nao houver base para escolher."""
    if not aachange or aachange in (".", ""):
        return aachange
    candidatos = [p for p in aachange.split(",") if p.count(":") >= 2]
    if not candidatos:
        return aachange

    def chave(parte):
        transcrito = parte.split(":")[1]
        # CDS maior primeiro; entre iguais, o menor numero de acesso.
        return (-comprimentos.get(transcrito, 0), numero_de_acesso(transcrito))

    return min(candidatos, key=chave)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="GTF do refGene usado pelo workflow")
    ap.add_argument("--aachange", required=True,
                    help="conteudo do campo AAChange.refGene")
    args = ap.parse_args()
    print(escolher(args.aachange, comprimentos_de_cds(args.gtf)))


if __name__ == "__main__":
    main()
