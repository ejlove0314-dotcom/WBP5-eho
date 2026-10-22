import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

df = pd.read_csv('C:/Temp/reinvent_wbp5/samples.csv')
results = []

for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol is None:
        continue
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    arom_rings = Descriptors.NumAromaticRings(mol)
    lipinski_pass = (mw < 500 and logp < 5 and hbd <= 5 and hba <= 10)
    wbp5_match = (arom_rings >= 1 and (hbd >= 1 or hba >= 2))
    if lipinski_pass and wbp5_match:
        results.append({
            'SMILES': row['SMILES'], 'MW': round(mw,1), 'LogP': round(logp,2),
            'HBD': hbd, 'HBA': hba, 'RotBonds': rotb, 'AromaticRings': arom_rings,
            'NLL': row['NLL']
        })

filtered = pd.DataFrame(results)
filtered = filtered.sort_values('NLL')
filtered.to_csv('C:/Temp/reinvent_wbp5/filtered_candidates.csv', index=False)
print(f'Original: {len(df)} molecules')
print(f'After filtering: {len(filtered)} molecules')
print(f'\nTop 15 candidates:')
print(filtered.head(15).to_string(index=False))