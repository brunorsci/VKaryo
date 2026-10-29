"""Figura 3 do Artigo 1: desenho da comparação com os cariótipos do DSMZ (A) e
recuperação das alterações declaradas contra a taxa de bloco espúrio em
cromossomos sem alteração declarada (B).

Entrada: Tabela D10 (Tabelas_Suplementares_DSMZ_PT.xlsx), reproduzida abaixo.
Saída: TIFF (LZW, 600 dpi) e PNG, 17,4 cm de largura, em português e inglês.
Uso: python figura3_validacao_dsmz.py [pasta_de_saida]
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
# Tabela D10: (declaradas, recuperadas, controle, espúrios)
D10 = [('diplo', 256, 189, 1796, 170), ('tri', 282, 185, 194, 50), ('tetra', 194, 119, 446, 58),
       ('spec', 85, 67, 886, 68), ('total', 732, 493, 2436, 278), ('plena', 732, 135, 2436, 26)]
DIST = [26, 16, 35, 15]  # linhagens: todas, 80-99%, 50-79%, <50% (Tabela S14)
NAVY, GRAY, LIGHT = '#1B2A4A', '#9AA3AF', '#EEF2F7'

TXT = {
 'pt': dict(
   A='A', B='B',
   ref1='Alterações-chave do LL-100\n(RT-PCR, FISH, sequenciamento)', ref2='Cariótipos em bandas do DSMZ\n(bandeamento G, 100 linhagens)',
   ev1='Avaliação I\nsensibilidade\n70 alterações', steps=['Notações de dosagem\n(732 no escopo;\n175 excluídas)',
   'Perfil de cobertura\njanelas de 0,25 Mb', 'Maior bloco no\nsentido esperado\n(±0,20; ≥3 Mb)',
   'Veredito:\nconcordante,\nrefinamento ou\nnão reproduzida'],
   ctrl='Controle: mesma regra em\ncromossomos sem alteração\ndeclarada (2.436 avaliações)', ev2='Avaliação II',
   vk='Saídas do VKaryo\n(fusões, SV,\nSNV, cobertura)',
   cats=['Próximo de diploide', 'Triploide', 'Tetraploide', 'Cariótipo completamente\nespecificado', 'Total', 'Total, só concordância\nplena'],
   C='C', ct='Linhagens segundo a proporção de alterações declaradas recuperadas (n = 92)', cl=['Todas', '80–99%', '50–79%', 'Menos de 50%'],
   l1='Alterações declaradas recuperadas', l2='Bloco espúrio no controle', xl='Proporção (%)', ratio='×'),
 'en': dict(
   A='A', B='B',
   ref1='LL-100 key alterations\n(RT-PCR, FISH, sequencing)', ref2='DSMZ banded karyotypes\n(G-banding, 100 cell lines)',
   ev1='Evaluation I\nsensitivity\n70 alterations', steps=['Dosage notations\n(732 in scope;\n175 excluded)',
   'Coverage profile\n0.25 Mb windows', 'Largest block in\nthe expected direction\n(±0.20; ≥3 Mb)',
   'Verdict:\nconcordant,\nrefinement or\nnot reproduced'],
   ctrl='Control: same rule on\nchromosomes without declared\nalteration (2,436 assessments)', ev2='Evaluation II',
   vk='VKaryo outputs\n(fusions, SV,\nSNV, coverage)',
   cats=['Near-diploid', 'Triploid', 'Tetraploid', 'Fully specified\nkaryotype', 'Total', 'Total, full\nconcordance only'],
   C='C', ct='Cell lines by proportion of declared alterations recovered (n = 92)', cl=['All', '80–99%', '50–79%', 'Less than 50%'],
   l1='Declared alterations recovered', l2='Spurious block in control', xl='Proportion (%)', ratio='×'),
}


def box(ax, x, y, w, h, t, fc=LIGHT, ec=NAVY, tc='#1F2328', fs=7, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.012', fc=fc, ec=ec, lw=0.8))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs, color=tc, fontweight='bold' if bold else 'normal', linespacing=1.15)


def arr(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='-|>', color='#4B5563', lw=0.8, mutation_scale=7))


def make(lang):
    T = TXT[lang]
    dec = ',' if lang == 'pt' else '.'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7})
    fig = plt.figure(figsize=(17.4 / 2.54, 18.5 / 2.54))
    # ---------- A: esquema
    ax = fig.add_axes([0.0, 0.61, 1.0, 0.39]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.text(0.01, 0.97, T['A'], fontsize=10, fontweight='bold', va='top')
    box(ax, 0.02, 0.42, 0.16, 0.22, T['vk'], fc='white')
    box(ax, 0.24, 0.76, 0.27, 0.18, T['ref1'])
    box(ax, 0.58, 0.76, 0.18, 0.18, T['ev1'], fc=NAVY, tc='white', bold=True)
    box(ax, 0.24, 0.44, 0.27, 0.18, T['ref2'])
    ax.text(0.625, 0.53, T['ev2'], ha='center', va='center', fontsize=7.5, fontweight='bold', color=NAVY)
    box(ax, 0.74, 0.42, 0.25, 0.22, T['ctrl'], fc='white', fs=6.2)
    arr(ax, 0.18, 0.57, 0.24, 0.85); arr(ax, 0.18, 0.52, 0.24, 0.53)
    arr(ax, 0.51, 0.85, 0.58, 0.85)
    for i, st in enumerate(T['steps']):
        x = 0.24 + i * 0.19
        box(ax, x, 0.04, 0.17, 0.24, st, fs=6)
        if i:
            arr(ax, x - 0.02, 0.16, x, 0.16)
    arr(ax, 0.325, 0.44, 0.325, 0.28)
    arr(ax, 0.895, 0.42, 0.895, 0.28)
    # ---------- B: barras
    bx = fig.add_axes([0.25, 0.28, 0.70, 0.30])
    bx.text(-0.33, 1.08, T['B'], transform=bx.transAxes, fontsize=10, fontweight='bold', va='top')
    n = len(D10); y = range(n)
    rec = [100 * r[2] / r[1] for r in D10]; fp = [100 * r[4] / r[3] for r in D10]
    h = 0.36
    bx.barh([i - h / 2 for i in y], rec, h, color=NAVY, label=T['l1'])
    bx.barh([i + h / 2 for i in y], fp, h, color=GRAY, label=T['l2'])
    for i in y:
        bx.text(rec[i] + 1, i - h / 2, f'{rec[i]:.1f}%'.replace('.', dec), va='center', fontsize=6)
        bx.text(fp[i] + 1, i + h / 2, f'{fp[i]:.1f}%'.replace('.', dec), va='center', fontsize=6, color='#4B5563')
        bx.text(101, i, f'{rec[i] / fp[i]:.1f}{T["ratio"]}'.replace('.', dec), va='center', ha='left', fontsize=6.5, fontweight='bold', color=NAVY)
    bx.set_yticks(list(y)); bx.set_yticklabels(T['cats'], fontsize=6.5); bx.invert_yaxis()
    bx.set_xlim(0, 100); bx.set_xlabel(T['xl'], fontsize=7)
    bx.axhline(4.5, color='#D1D5DB', lw=0.6, ls='--')
    for s in ('top', 'right'):
        bx.spines[s].set_visible(False)
    bx.tick_params(labelsize=6.5)
    bx.legend(loc='lower center', bbox_to_anchor=(0.45, -0.36), ncol=2, frameon=False, fontsize=6.5)
    # ---------- C: distribuição por linhagem (Tabela S14)
    cx = fig.add_axes([0.25, 0.06, 0.70, 0.06])
    cx.text(-0.33, 1.9, T['C'], transform=cx.transAxes, fontsize=10, fontweight='bold', va='top')
    cx.text(0, 1.35, T['ct'], transform=cx.transAxes, fontsize=7, va='bottom')
    left = 0
    cols = [NAVY, '#4A5D82', '#8FA3C4', '#D6DEEA']
    for (lab, n), col in zip(zip(T['cl'], DIST), cols):
        cx.barh(0, n, left=left, color=col, edgecolor='white', lw=0.6)
        cx.text(left + n / 2, 0, f'{n}', ha='center', va='center', fontsize=7, color='white' if col in (NAVY, '#4A5D82') else '#1F2328', fontweight='bold')
        cx.text(left + n / 2, -0.75, lab, ha='center', va='top', fontsize=6)
        left += n
    cx.set_xlim(0, 92); cx.set_ylim(-0.5, 0.5); cx.axis('off')
    name = 'Figura3_validacao_DSMZ_PT' if lang == 'pt' else 'Figure3_DSMZ_validation_EN'
    fig.savefig(f'{OUT}/{name}.png', dpi=600, facecolor='white')
    plt.close(fig)
    im = Image.open(f'{OUT}/{name}.png').convert('RGB')
    im.save(f'{OUT}/{name}.tif', compression='tiff_lzw', dpi=(600, 600))
    print(name, im.size)


for lg in ('pt', 'en'):
    make(lg)
