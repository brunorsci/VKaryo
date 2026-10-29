"""Taxa de bloco espurio: a mesma regra aplicada a bracos que o cariotipo NAO
declara alterados. Sem esse controle, a taxa de recuperacao nao e interpretavel.

Classificador congelado, sem mudanca de parametro.
"""
import collections, io, json, os, re, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confrontar_cariotipo_bandas as C

OUT=os.environ.get('VKARYO_OUTPUT','results/output_vkaryo')
DSMZ=os.environ.get('VKARYO_KARYOTYPES','.')
REF=os.environ.get('VKARYO_REFERENCE','reference')
AUTOS=[str(i) for i in range(1,23)]
PARC=os.environ.get('VKARYO_TMP','.')+'/fp_parcial.json'
CARGA=34.0

def carga():
    try: return os.getloadavg()[0]
    except OSError: return 0.0

def main():
    kar={}
    with io.open(f'{DSMZ}/dsmz_karyotypes.tsv',encoding='utf-8') as fh:
        ch=fh.readline().rstrip('\n').split('\t')
        for l in fh:
            v=l.rstrip('\n').split('\t')
            if v and v[0]: kar[v[0]]=dict(zip(ch,v))
    dna={}
    with io.open(f'{REF}/ll100_sample_info.tsv',encoding='utf-8') as fh:
        next(fh)
        for l in fh:
            c=l.rstrip('\n').split('\t')
            if len(c)>2: dna[c[2]]=c[1]

    feito=json.load(io.open(PARC,encoding='utf-8')) if os.path.exists(PARC) else {}
    nomes=[n for n in sorted(kar) if n not in feito]
    print(f'a processar: {len(nomes)}', flush=True)

    for i,nome in enumerate(nomes,1):
        while carga()>CARGA:
            print(f'  [carga {carga():.1f}, pausa]',flush=True); time.sleep(90)
        acc=dna.get(nome)
        cnr=f'{OUT}/cnvkit/{acc}_markduplicates_recal.cnr' if acc else None
        if not cnr or not os.path.exists(cnr):
            feito[nome]={'erro':1}; continue
        texto=kar[nome].get('karyotype','').split('Molec')[0]
        # cromossomos tocados por qualquer notacao: saem do controle
        tocados=set()
        for m in re.finditer(r'(?:[+-]|del\(|dup\(|add\(|der\(|inv\(|i\(|t\()([0-9]{1,2}|X|Y)',texto):
            tocados.add(m.group(1))
        for m in re.finditer(r't\(([0-9XY;]+)\)',texto):
            for x in m.group(1).split(';'):
                if x: tocados.add(x)
        perf,_=C.perfil(cnr)
        res={'concordante':0,'refino':0,'avaliados':0}
        for cr in AUTOS:
            if cr in tocados: continue
            for sent in ('perda','ganho'):
                # mesma combinacao de bracos usada nos eventos declarados
                melhor='nao_reproduzido'
                for br in ('p','q'):
                    cls,_,_=C.classificar(perf,cr,br,sent)
                    if cls=='concordante': melhor='concordante'; break
                    if cls=='concordante_com_refino' and melhor=='nao_reproduzido': melhor='concordante_com_refino'
                res['avaliados']+=1
                if melhor=='concordante': res['concordante']+=1
                elif melhor=='concordante_com_refino': res['refino']+=1
        # sentido oposto nos eventos declarados
        feito[nome]={'controle':res,'ploidia':kar[nome].get('ploidy',''),
                     'modal_low':kar[nome].get('modal_low',''),'modal_high':kar[nome].get('modal_high',''),
                     'cariotipo':kar[nome].get('karyotype','')}
        io.open(PARC,'w').write(json.dumps(feito))
        print(f'  {i:3d}/{len(nomes)} {nome:<14} controle {res["concordante"]}+{res["refino"]}/{res["avaliados"]} [carga {carga():.1f}]',flush=True)
    print('FIM',flush=True)

if __name__=='__main__': main()
