# 62_gse12276_chars.py - export all GSE12276 characteristics
import gzip
sample_ids=[]
char_lines=[]
with gzip.open(r'D:\BCBM_Project\data\raw\GSE12276_series_matrix.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            sample_ids=[v.strip().strip('"') for v in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            char_lines.append([v.strip().strip('"') for v in line.split('\t')[1:]])
        elif line.startswith('!Sample_title'):
            titles=[v.strip().strip('"') for v in line.split('\t')[1:]]

print(f"Total samples: {len(sample_ids)}")
print(f"Characteristic lines: {len(char_lines)}")
# build dict
data={sid:{'title':titles[i] if i<len(titles) else ''} for i,sid in enumerate(sample_ids)}
for vals in char_lines:
    for i,sid in enumerate(sample_ids):
        if i<len(vals) and ':' in vals[i]:
            k,v=vals[i].split(':',1)
            data[sid][k.strip()]=v.strip()

# print all unique keys
all_keys=set()
for sid,d in data.items(): all_keys.update(d.keys())
print("\nAll characteristic keys:")
for k in sorted(all_keys): print(f"  {k}")

# print first 5 samples full
print("\n=== First 5 samples ===")
for sid in sample_ids[:5]:
    print(f"\n{sid} ({data[sid].get('title','')}):")
    for k,v in data[sid].items():
        if k!='title': print(f"  {k}: {v}")

# check for brain relapse / metastasis site
print("\n=== Brain-related fields ===")
for k in sorted(all_keys):
    if any(w in k.lower() for w in ['brain','metast','relapse','site','surviv','event','status','vital','follow']):
        vals=set(data[sid].get(k,'') for sid in sample_ids)
        print(f"  {k}: {sorted(vals)[:10]}")

# save all
import csv
with open(r'D:\BCBM_Project\results\tables\GSE12276_all_characteristics.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['sample']+sorted(all_keys))
    w.writeheader()
    for sid in sample_ids:
        row={'sample':sid}; row.update(data[sid]); w.writerow(row)
print("\nSaved GSE12276_all_characteristics.csv")
