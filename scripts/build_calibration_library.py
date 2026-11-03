"""Build versioned offline ion assets from broad, formula-bearing ion tables.
Usage: python build_library.py POSITIVE.tsv NEGATIVE.tsv POS_CHECKS.tsv NEG_CHECKS.tsv OUTPUT_DIR
Sources are read only. No network queries or manuscript paths are needed at runtime.
"""
import csv,sys,re,json,hashlib,pathlib,datetime,collections
ATOMS={'C':12.,'H':1.00782503207,'N':14.0030740048,'O':15.99491461956,'P':30.97376163,'S':31.97207100,'Fe':55.9349375,'Na':22.98976928,'K':38.96370649,'Cl':34.96885268}
E=0.00054857990946
SHIFTS={'[M+H]+':{'H':1},'[M+Na]+':{'Na':1},'[M+K]+':{'K':1},'[M+NH4]+':{'N':1,'H':4},'[M+H-H2O]+':{'H':-1,'O':-1},'[M]+':{},'[M-H]-':{'H':-1},'[M+Cl]-':{'Cl':1},'[M+HCOO]-':{'C':1,'H':1,'O':2},'[M-CH3]-':{'C':-1,'H':-3},'[M-H-H2O]-':{'H':-3,'O':-1}}
def formula(s):
 parts=re.findall(r'([A-Z][a-z]?)(\d*)',s)
 if ''.join(a+b for a,b in parts)!=s: raise ValueError(s)
 out=collections.Counter()
 for a,b in parts: out[a]+=int(b or 1)
 return out
def canonical(d): return ''.join(a+(str(d[a]) if d[a]!=1 else '') for a in sorted(d,key=lambda a:(a!='C',a!='H',a)) if d[a])
def build(src,checks,polarity,out):
 rows=list(csv.DictReader(open(src),delimiter='\t')); ver=list(csv.DictReader(open(checks),delimiter='\t'))
 result=[]
 for r in rows:
  ion=formula(r['formula']); ad=r['adduct']; charge=1 if polarity=='positive' else -1
  ion.update(SHIFTS[ad]); assert all(v>=0 for v in ion.values())
  computed=sum(ATOMS[a]*n for a,n in ion.items())-charge*E
  # Preserve the documented instrument lock reference, with its difference explicit.
  override=r['name']=='Leu-Enk (lock mass)'
  theo=float(r['mz']) if override else computed
  assert override or abs(computed-float(r['mz']))<0.000002,(r['name'],computed,r['mz'])
  direct=[v for v in ver if v['name']==r['name']]
  ev=direct or [v for v in ver if v['cls']==r['class'] and v['check']=='class rule']
  result.append(dict(name=r['name'],adduct=ad,formula=r['formula'],ion_formula=canonical(ion),charge=charge,theo=format(theo,'.15g'),class_name=r['class'],prior=r['prior'],reference_family=r['name'],source_ids='; '.join(sorted({v[k] for v in ev for k in ('db1','db2') if v[k]})),source_scope='compound record' if direct else 'computed class rule; database examples',mass_basis='instrument lock reference' if override else 'ion formula minus charge times electron mass',library_version='2026.10.27.1'))
 # The class-generating source already contains lipids and small metabolites; include all entries.
 for r in result: r['reference_id']='ion_'+hashlib.sha256((r['ion_formula']+'|'+str(r['charge'])).encode()).hexdigest()[:16]
 fields=list(result[0]);path=out/(polarity+'.tsv')
 with path.open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(result)
 return dict(rows=len(result),classes=sorted({r['class_name'] for r in result}),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(pathlib.Path(src).read_bytes()).hexdigest(),verification_sha256=hashlib.sha256(pathlib.Path(checks).read_bytes()).hexdigest())
if __name__=='__main__':
 pos,neg,pv,nv,dest=sys.argv[1:];out=pathlib.Path(dest);out.mkdir(parents=True,exist_ok=True)
 manifest={'version':'2026.10.27.1','description':'Formula-derived lipid and small-metabolite ion library. Database example support is distinguished from direct records. Singly charged ions only. Not an exhaustive chemical database.','positive':build(pos,pv,'positive',out),'negative':build(neg,nv,'negative',out)}
 # A common formula-bearing compound index links both polarity libraries and their allowed adducts.
 compounds={}
 for pol in ('positive','negative'):
  for row in csv.DictReader((out/(pol+'.tsv')).open(),delimiter='\t'):
   key=(row['name'],row['formula'])
   if key not in compounds: compounds[key]={'name':key[0],'formula':key[1],'class_name':row['class_name'],'positive_adducts':set(),'negative_adducts':set()}
   compounds[key][pol+'_adducts'].add(row['adduct'])
 with (out/'compounds.tsv').open('w') as f:
  writer=csv.DictWriter(f,['name','formula','class_name','positive_adducts','negative_adducts'],delimiter='\t');writer.writeheader()
  for row in compounds.values():
   row['positive_adducts']=';'.join(sorted(row['positive_adducts']));row['negative_adducts']=';'.join(sorted(row['negative_adducts']));writer.writerow(row)
 manifest['compounds']=len(compounds)
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({k:manifest[k]['rows'] for k in ('positive','negative')}))
