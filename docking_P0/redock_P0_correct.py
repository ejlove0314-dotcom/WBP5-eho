"""
WBP5 재도킹 — DoGSiteScorer P_0 포켓 (진짜 P_0, Drug Score 0.837)

배경
----
기존 도킹은 그리드 중심을 잔기 64-74 의 원자 좌표 평균 (2.402, -3.010, 1.561) 로
잡았으나, DoGSiteScorer 작업 0f462a25 의 출력을 검증한 결과 그 위치는
포켓 P_3 (Drug Score 0.281) 이며, Drug Score 0.837 을 받은 P_0 는
잔기 38-58 에 위치하고 두 부위는 32.75 A 떨어져 있다.

이 스크립트는 그리드 중심을 P_0 포켓의 '부피 격자(ccp4)'에서 직접 계산한다.
좌표를 손으로 입력하지 않으므로 전사 오류가 발생하지 않는다.

전제
----
  - 리간드 pdbqt 10 종이 이미 준비되어 있을 것 (기존 도킹과 동일 파일 재사용)
  - 수용체 pdbqt 가 PyRx 로 생성되어 있을 것
  - DoGSiteScorer 출력 P_0_gpsAll.ccp4 가 있을 것

실행
----
  python redock_P0_correct.py
"""

import os
import re
import struct
import subprocess
import statistics

import numpy as np

# ── 경로 (기존 prep_and_dock_all10_v2.py 와 동일하게 맞출 것) ─────────────
VINA      = r'"C:\Program Files (x86)\PyRx\vina.exe"'
RECEPTOR  = os.path.join(os.environ["USERPROFILE"],
                         r".mgltools\PyRx\Macromolecules"
                         r"\WBP5_AlphaFold_model1\WBP5_AlphaFold_model1.pdbqt")
LIGAND_DIR = r"C:\Users\user\Desktop\wbp5\ligands_pdbqt"
OUTPUT_DIR = r"C:\Users\user\Desktop\wbp5\docking_results_P0"
CCP4       = (r"C:\Users\user\Desktop\wbp5"
              r"\wbp5alphafoldpdb0f462a25-d694-42fa-a0c1-f397f6b16f2f_P_0_gpsAll.ccp4")

SEEDS        = [42, 1042, 2042, 3042, 4042]   # 기존과 동일
BOX          = (25.0, 25.0, 25.0)             # 기존과 동일 (P_0 실크기 22.8x15.0x13.5 를 포함)
EXHAUSTIVENESS = 8                            # 기존과 동일
NUM_MODES      = 9                            # 기존과 동일


def pocket_center_from_ccp4(path):
    """CCP4 맵에서 밀도 > 0 인 복셀의 중심을 직교좌표로 반환한다."""
    d = open(path, "rb").read()
    nc, nr, ns, _mode      = struct.unpack("<4i", d[0:16])
    ncs, nrs, nss          = struct.unpack("<3i", d[16:28])
    nx, ny, nz             = struct.unpack("<3i", d[28:40])
    cell                   = struct.unpack("<6f", d[40:64])
    mapc, mapr, maps       = struct.unpack("<3i", d[64:76])
    nsym                   = struct.unpack("<i",  d[92:96])[0]

    vox  = np.array(cell[:3]) / np.array([nx, ny, nz])
    vals = np.frombuffer(d, dtype="<f4", count=nc * nr * ns,
                         offset=1024 + nsym).reshape(ns, nr, nc)
    occ = np.argwhere(vals > 0)
    if not len(occ):
        raise SystemExit("[중단] ccp4 에서 포켓 복셀을 찾지 못했습니다.")

    s_m, r_m, c_m = occ.mean(0)
    axis = {mapc: 0, mapr: 1, maps: 2}          # 1=X, 2=Y, 3=Z
    g = np.zeros(3)
    g[axis[1]] = c_m + ncs
    g[axis[2]] = r_m + nrs
    g[axis[3]] = s_m + nss

    extent = np.zeros(3)
    span = occ.max(0) - occ.min(0) + 1
    extent[axis[1]], extent[axis[2]], extent[axis[3]] = span[2], span[1], span[0]

    volume = len(occ) * float(np.prod(vox))
    return g * vox, extent * vox, volume


def parse_affinity(log_text):
    """Vina 로그에서 mode 1 의 결합에너지를 뽑는다."""
    for line in log_text.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def main():
    for label, path in [("vina", VINA.strip('"')), ("수용체", RECEPTOR),
                        ("리간드 폴더", LIGAND_DIR), ("ccp4", CCP4)]:
        if not os.path.exists(path):
            raise SystemExit(f"[중단] {label} 경로를 찾을 수 없습니다:\n  {path}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    center, extent, volume = pocket_center_from_ccp4(CCP4)
    print("=" * 68)
    print("Grid box from DoGSiteScorer pocket P_0 (drugScore 0.837)")
    print("=" * 68)
    print(f"  center  : X={center[0]:.3f}  Y={center[1]:.3f}  Z={center[2]:.3f}")
    print(f"  extent  : {np.round(extent, 1)} A")
    print(f"  volume  : {volume:.2f} A^3   <- desc.txt 의 229.62 와 일치해야 정상")
    print(f"  box     : {BOX[0]} x {BOX[1]} x {BOX[2]} A")
    print(f"  (참고) 기존 그리드 (2.402, -3.010, 1.561) 로부터 "
          f"{np.linalg.norm(center - np.array([2.402, -3.010, 1.561])):.2f} A")
    print()

    ligands = sorted(f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdbqt"))
    if not ligands:
        raise SystemExit(f"[중단] {LIGAND_DIR} 에 pdbqt 파일이 없습니다.")
    print(f"리간드 {len(ligands)} 종 x 시드 {len(SEEDS)} = {len(ligands)*len(SEEDS)} 회 실행\n")

    results = {}
    for lig in ligands:
        name = os.path.splitext(lig)[0]
        energies = {}
        for seed in SEEDS:
            out = os.path.join(OUTPUT_DIR, f"{name}_seed{seed}_out.pdbqt")
            cmd = (f'{VINA} --receptor "{RECEPTOR}" '
                   f'--ligand "{os.path.join(LIGAND_DIR, lig)}" '
                   f'--center_x {center[0]:.3f} --center_y {center[1]:.3f} '
                   f'--center_z {center[2]:.3f} '
                   f'--size_x {BOX[0]} --size_y {BOX[1]} --size_z {BOX[2]} '
                   f'--exhaustiveness {EXHAUSTIVENESS} --num_modes {NUM_MODES} '
                   f'--seed {seed} --out "{out}"')
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            e = parse_affinity(proc.stdout)
            if e is None:
                print(f"  [경고] {name} seed={seed} 실패\n{proc.stderr[:300]}")
                continue
            energies[seed] = e
        if energies:
            results[name] = energies
            vals = list(energies.values())
            sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
            print(f"  {name:<18} mean={statistics.mean(vals):7.3f}  SD={sd:5.3f}  "
                  f"{['%.1f' % v for v in vals]}")

    csv_path = os.path.join(OUTPUT_DIR, "docking_replicates_P0_summary.csv")
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("ligand,n_replicates,mean_kcal_per_mol,SD,min,max,"
                 + ",".join(f"seed_{s}" for s in SEEDS) + "\n")
        for name, en in results.items():
            vals = [en.get(s) for s in SEEDS]
            got = [v for v in vals if v is not None]
            sd = statistics.stdev(got) if len(got) > 1 else 0.0
            fh.write(f"{name},{len(got)},{statistics.mean(got):.3f},{sd:.3f},"
                     f"{min(got)},{max(got)},"
                     + ",".join("" if v is None else f"{v}" for v in vals) + "\n")

    print(f"\n완료 → {csv_path}")
    print("이 CSV 를 올려주시면 통계 재계산과 본문 반영을 진행합니다.")


if __name__ == "__main__":
    main()
