"""
Candidate_1 재도킹 — 피리딘형 (REINVENT4 원본, MW 347.44)

배경
----
all_candidates.smi 의 Candidate_1 문자열 'Nc2ccnc2' 는 kekulize 불가한
손상된 SMILES 였다. 2026-08-07 세션에서 이를 피롤형('Nc2cc[nH]c2', MW 335.43)
으로 복원했으나, REINVENT4 자체 리포트의 Rank 1 물성(MW 347, TPSA 88.16)은
피리딘형(MW 347.44, TPSA 88.16)과 일치한다. 따라서 피리딘형으로 되돌린다.

이 스크립트는
  1) 피리딘형 Candidate_1 의 3D 좌표를 생성하고 (RDKit ETKDGv3 + MMFF94)
  2) Meeko 로 pdbqt 를 만든 뒤
  3) P_0 포켓(ccp4 부피 격자에서 중심 산출)에서 5 시드로 도킹한다.

나머지 9 종은 이미 P_0 에서 도킹이 끝났으므로 다시 돌리지 않는다.

실행
----
  cd C:\\Users\\user\\Desktop\\wbp5
  python redock_candidate1_pyridine.py
"""

import os
import re
import struct
import subprocess
import statistics

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# ── 대상 ────────────────────────────────────────────────────────────
NAME  = "Candidate_1"
SMILES = "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccncc2)cc1"   # pyridin-4-yl
EXPECT_MW = 347.44

# ── 경로 (redock_P0_correct.py 와 동일) ─────────────────────────────
VINA       = r'"C:\Program Files (x86)\PyRx\vina.exe"'
RECEPTOR   = os.path.join(os.environ["USERPROFILE"],
                          r".mgltools\PyRx\Macromolecules"
                          r"\WBP5_AlphaFold_model1\WBP5_AlphaFold_model1.pdbqt")
LIGAND_DIR = r"C:\Users\user\Desktop\wbp5\ligands_pdbqt"
OUTPUT_DIR = r"C:\Users\user\Desktop\wbp5\docking_results_P0"
CCP4       = (r"C:\Users\user\Desktop\wbp5"
              r"\wbp5alphafoldpdb0f462a25-d694-42fa-a0c1-f397f6b16f2f_P_0_gpsAll.ccp4")

SEEDS = [42, 1042, 2042, 3042, 4042]
BOX = (25.0, 25.0, 25.0)
EXHAUSTIVENESS = 8
NUM_MODES = 9
EMBED_SEED = 42          # 기존 리간드 준비와 동일


def pocket_center_from_ccp4(path):
    d = open(path, "rb").read()
    nc, nr, ns, _m   = struct.unpack("<4i", d[0:16])
    ncs, nrs, nss    = struct.unpack("<3i", d[16:28])
    nx, ny, nz       = struct.unpack("<3i", d[28:40])
    cell             = struct.unpack("<6f", d[40:64])
    mapc, mapr, maps = struct.unpack("<3i", d[64:76])
    nsym             = struct.unpack("<i",  d[92:96])[0]
    vox  = np.array(cell[:3]) / np.array([nx, ny, nz])
    vals = np.frombuffer(d, dtype="<f4", count=nc * nr * ns,
                         offset=1024 + nsym).reshape(ns, nr, nc)
    occ = np.argwhere(vals > 0)
    s_m, r_m, c_m = occ.mean(0)
    ax = {mapc: 0, mapr: 1, maps: 2}
    g = np.zeros(3)
    g[ax[1]] = c_m + ncs
    g[ax[2]] = r_m + nrs
    g[ax[3]] = s_m + nss
    return g * vox, len(occ) * float(np.prod(vox))


def build_ligand():
    m = Chem.MolFromSmiles(SMILES)
    if m is None:
        raise SystemExit("[중단] SMILES 파싱 실패")
    mw = Descriptors.MolWt(m)
    print(f"  분자식 {rdMolDescriptors.CalcMolFormula(m)}  MW={mw:.2f}  "
          f"heavy={m.GetNumHeavyAtoms()}")
    if abs(mw - EXPECT_MW) > 0.05:
        raise SystemExit(f"[중단] MW가 {EXPECT_MW} 와 다릅니다 — SMILES 확인 필요")

    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = EMBED_SEED
    p.maxIterations = 2000
    p.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mh, p) == -1:
        p.useRandomCoords = True
        AllChem.EmbedMolecule(mh, p)
        print("  (임베딩: randomCoords 사용)")
    ff = "MMFF94" if AllChem.MMFFHasAllMoleculeParams(mh) else "UFF"
    (AllChem.MMFFOptimizeMolecule if ff == "MMFF94"
     else AllChem.UFFOptimizeMolecule)(mh, maxIters=1000)
    print(f"  최소화: {ff}")

    mh.SetProp("_Name", NAME)
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation().prepare(mh)[0])
    if not ok:
        raise SystemExit(f"[중단] pdbqt 변환 실패: {err}")
    path = os.path.join(LIGAND_DIR, f"{NAME}.pdbqt")
    if os.path.exists(path):
        os.replace(path, path + ".pyrrole_backup")
        print(f"  기존 pdbqt 백업 → {NAME}.pdbqt.pyrrole_backup")
    open(path, "w").write(txt)
    tors = [l for l in txt.splitlines() if l.startswith("TORSDOF")][0].split()[1]
    print(f"  pdbqt 생성 완료 (TORSDOF {tors})")
    return path


def parse_affinity(log):
    for line in log.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def main():
    for label, path in [("vina", VINA.strip('"')), ("수용체", RECEPTOR),
                        ("리간드 폴더", LIGAND_DIR), ("ccp4", CCP4)]:
        if not os.path.exists(path):
            raise SystemExit(f"[중단] {label} 경로 없음:\n  {path}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 64)
    print("STEP 1  리간드 생성 (pyridin-4-yl)")
    print("=" * 64)
    lig = build_ligand()

    center, volume = pocket_center_from_ccp4(CCP4)
    print("\n" + "=" * 64)
    print("STEP 2  P_0 그리드")
    print("=" * 64)
    print(f"  center X={center[0]:.3f} Y={center[1]:.3f} Z={center[2]:.3f}")
    print(f"  volume {volume:.2f} A^3   <- 230.11 이면 정상")

    print("\n" + "=" * 64)
    print("STEP 3  도킹 (5 시드)")
    print("=" * 64)
    energies = {}
    for seed in SEEDS:
        out = os.path.join(OUTPUT_DIR, f"{NAME}_seed{seed}_out.pdbqt")
        cmd = (f'{VINA} --receptor "{RECEPTOR}" --ligand "{lig}" '
               f'--center_x {center[0]:.3f} --center_y {center[1]:.3f} '
               f'--center_z {center[2]:.3f} '
               f'--size_x {BOX[0]} --size_y {BOX[1]} --size_z {BOX[2]} '
               f'--exhaustiveness {EXHAUSTIVENESS} --num_modes {NUM_MODES} '
               f'--seed {seed} --out "{out}"')
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        e = parse_affinity(proc.stdout)
        if e is None:
            print(f"  seed {seed}: 실패\n{proc.stderr[:300]}")
            continue
        energies[seed] = e
        print(f"  seed {seed:>5}: {e:6.1f} kcal/mol")

    vals = list(energies.values())
    sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
    print("\n" + "=" * 64)
    print(f"  {NAME}  mean = {statistics.mean(vals):.3f}  SD = {sd:.3f}  "
          f"(n={len(vals)})")
    print(f"  참고) 피롤형 P_0 결과: -5.340 +- 0.089")
    print("=" * 64)

    csv = os.path.join(OUTPUT_DIR, "candidate1_pyridine_P0.csv")
    with open(csv, "w", encoding="utf-8") as fh:
        fh.write("ligand,n_replicates,mean_kcal_per_mol,SD,min,max,"
                 + ",".join(f"seed_{s}" for s in SEEDS) + "\n")
        fh.write(f"{NAME},{len(vals)},{statistics.mean(vals):.3f},{sd:.3f},"
                 f"{min(vals)},{max(vals)},"
                 + ",".join(str(energies.get(s, "")) for s in SEEDS) + "\n")
    print(f"\n완료 → {csv}")


if __name__ == "__main__":
    main()
