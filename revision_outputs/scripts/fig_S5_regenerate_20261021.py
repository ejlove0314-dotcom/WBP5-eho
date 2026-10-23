#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S5 재생성 (YAP/TAZ) — matched-null 10,000 반영 (2026-10-21, v2)
===================================================================
모든 점수를 전장 매트릭스에서 직접 계산(z-module)하므로 별도 점수 파일 불필요.
  WBP5~YAP-target, YAP1 mRNA, p-EMT, 그리고 p-EMT 통제 부분상관 재계산
  Panel B: 10,000 matched-null (발현구간 매칭 22-gene) — perm p ≈ 0.0002

필요 파일(바탕화면 WBP5 revision):
  GSE103322_HNSCC_all_data.txt (또는 .txt.gz)   전장 매트릭스

실행:
  python fig_S5_regenerate_20261021.py --matrix GSE103322_HNSCC_all_data.txt --outdir figures_revision

산출: figures_revision/Fig_S5_updated_300dpi.png
      figures_revision/YAP_matched_null_distribution.csv
화면에 rho_obs≈0.173 / perm_p≈0.0002 / partial≈0.040 이 찍히면 §3.3과 일치.
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE_ALIASES=["WBP5","TCEAL9"]
YAP22=["MYOF","AMOTL2","LATS2","CTGF","CYR61","ANKRD1","ASAP1","AXL","F3","IGFBP3",
       "CRIM1","FJX1","FOXF2","GADD45A","CCDC80","NT5E","DOCK5","PTPN14","ARHGEF17",
       "NUAK2","TGFB2","RBMS3"]
# Puram p-EMT 상위 30 (P0_signature_gene_mapping.csv 순서)
PURAM30=["SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC",
         "MMP2","EMP3","INHBA","LAMB3","VIM","SEMA3C","PRKCDBP","ANXA5","DHRS7","ITGB1",
         "ACTN1","CXCR7","ITGB6","IGFBP7","THBS1","PTHLH","TNFRSF6B","PDLIM7","CAV1","DKK3"]
SEED=42; NPERM=10000; NBIN=25

def _comp(p):
    with open(p,"rb") as f: return "gzip" if f.read(2)==b"\x1f\x8b" else None
def parse_patient(name):
    name=str(name).strip().strip("'\"")
    if re.match(r"^HNSCC[_]?combo",name,re.I): return None
    m=re.match(r"^HN(?:SCC)?[_]?(\d+)",name,re.I); return f"HN{int(m.group(1))}" if m else None

def load_matrix_malignant(path):
    if not os.path.exists(path):
        raise SystemExit(f"매트릭스 없음: {path}\n폴더: {sorted(os.listdir('.'))}")
    df=pd.read_csv(path,sep="\t",index_col=0,low_memory=False,compression=_comp(path))
    df.index=[str(i).strip().strip("'\"").strip() for i in df.index]
    df.columns=[str(c).strip().strip("'\"").strip() for c in df.columns]
    nk=lambda s:" ".join(str(s).replace('"','').replace("'","").lower().split())
    crow=next((i for i in df.index if nk(i)=="classified as cancer cell"),None)
    if crow is None: raise SystemExit("'classified as cancer cell' 행 없음")
    is_cancer=pd.to_numeric(df.loc[crow],errors="coerce").fillna(0).astype(int)
    mkeys={"processed by maxima enzyme","maxima enriched","lymph node",
           "classified as cancer cell","classified as non-cancer cells","non-cancer cell type"}
    expr=df.drop(index=[i for i in df.index if nk(i) in mkeys]).apply(pd.to_numeric,errors="coerce")
    pats=pd.Series([parse_patient(c) for c in expr.columns],index=expr.columns)
    keep=(is_cancer.reindex(expr.columns)==1)&(pats.notna())
    mal=expr.columns[keep.values]
    sys.stderr.write(f"악성세포 {len(mal):,} / 유전자 {expr.shape[0]:,}\n")
    return expr[mal]

def zmodule(expr_sub, genes):
    g=[x for x in genes if x in expr_sub.index]
    sub=expr_sub.loc[g].astype(float)
    z=sub.sub(sub.mean(axis=1),axis=0).div(sub.std(axis=1).replace(0,np.nan),axis=0)
    return z.mean(axis=0,skipna=True).values, g

def partial_spearman(x,y,z):
    from scipy.stats import rankdata
    rx,ry,rz=rankdata(x),rankdata(y),rankdata(z)
    res=lambda t,u:(t-np.c_[np.ones_like(u),u]@np.linalg.lstsq(np.c_[np.ones_like(u),u],t,rcond=None)[0])
    return np.corrcoef(res(rx,rz),res(ry,rz))[0,1]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--matrix",required=True)
    ap.add_argument("--scores",default=None)   # 미사용(호환용)
    ap.add_argument("--yap-map",default=None)
    ap.add_argument("--outdir",default="figures_revision")
    a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True); rng=np.random.default_rng(SEED)

    expr=load_matrix_malignant(a.matrix)
    gene=next((g for g in GENE_ALIASES if g in expr.index),None)
    if gene is None: raise SystemExit("WBP5/TCEAL9 없음")
    w=expr.loc[gene].astype(float).values

    yap=YAP22
    if a.yap_map and os.path.exists(a.yap_map):
        mp=pd.read_csv(a.yap_map)
        c0=[c for c in mp.columns if c.lower() in ("list_symbol","symbol","gene")]
        if c0: yap=[str(x) for x in mp[c0[0]].tolist()]

    yz,yused=zmodule(expr,yap)
    pz,pused=zmodule(expr,PURAM30)
    yap1=expr.loc["YAP1"].astype(float).values if "YAP1" in expr.index else None

    rho_obs=stats.spearmanr(w,yz)[0]
    rho_pemt=stats.spearmanr(w,pz)[0]
    rho_yap1=stats.spearmanr(w,yap1)[0] if yap1 is not None else np.nan
    partial_rho=partial_spearman(w,yz,pz)

    # Panel B: 발현 + 검출률 매칭 10,000 null (본문 §3.3과 동일: P5 방법)
    names=np.array(expr.index)
    gmean=expr.mean(axis=1).values
    det=np.asarray((expr.values>0).mean(axis=1)).ravel()
    ebin=pd.qcut(pd.Series(gmean).rank(method="first"),NBIN,labels=False).values
    dbin=pd.qcut(pd.Series(det).rank(method="first"),20,labels=False).values
    key=ebin.astype(int)*100+dbin.astype(int)
    Z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1).replace(0,np.nan),axis=0).values
    yidx=np.array([np.where(names==g)[0][0] for g in yused])
    ykey=key[yidx]
    pools={k:np.where(key==k)[0] for k in np.unique(key)}
    null=np.empty(NPERM)
    for i in range(NPERM):
        pick=[int(rng.choice(pools[k])) for k in ykey]
        null[i]=stats.spearmanr(w,np.nanmean(Z[pick],axis=0))[0]
    emp_p=(np.sum(np.abs(null)>=abs(rho_obs))+1)/(NPERM+1)
    p95=np.percentile(np.abs(null),95)
    pd.DataFrame({"null_rho":null}).to_csv(f"{a.outdir}/YAP_matched_null_distribution.csv",index=False)

    print("\n=== 재계산 확인 (P5와 대조) ===")
    print(f" YAP genes {len(yused)}/22, Puram30 {len(pused)}/30")
    print(f" rho_obs={rho_obs:.4f}  [P5 0.173]   rho_pEMT={rho_pemt:.4f}  [0.214]   rho_YAP1={rho_yap1 if yap1 is None else round(rho_yap1,4)}  [0.019]")
    print(f" matched-null mean={null.mean():.4f}  95%|rho|={p95:.4f}  perm_p={emp_p:.4f}  [P5 0.0002]")
    print(f" partial(pEMT)={partial_rho:.4f}  [P5 0.040]")

    # ---- figure ----
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(15,4.4),dpi=300); gs=fig.add_gridspec(1,3,wspace=0.32)
    axA=fig.add_subplot(gs[0]); hb=axA.hexbin(w,yz,gridsize=32,cmap="viridis",mincnt=1)
    axA.set_xlabel("WBP5  log2(TPM/10+1)"); axA.set_ylabel("YAP/TAZ target-signature score (z)")
    axA.set_title(f"Spearman $\\rho$={rho_obs:.2f}",fontsize=10)
    cb=fig.colorbar(hb,ax=axA,fraction=0.046,pad=0.04); cb.set_label("Cell count",fontsize=9)
    axB=fig.add_subplot(gs[1]); axB.hist(null,bins=40,color="#bdbdbd")
    axB.axvline(p95,color="k",ls="--",lw=1.1); axB.axvline(rho_obs,color="crimson",lw=2)
    axB.set_xlabel("Spearman $\\rho$ (WBP5 vs matched random 22-gene score)"); axB.set_ylabel("count")
    axB.set_title(f"Matched null (10,000); empirical p={emp_p:.4f}",fontsize=10)
    axC=fig.add_subplot(gs[2])
    labs=["YAP-target","YAP1 mRNA","p-EMT","YAP-target|\np-EMT(partial)"]
    vals=[rho_obs,rho_yap1,rho_pemt,partial_rho]; cols=["#2c7fb8","#cccccc","#31a354","#756bb1"]
    bars=axC.bar(labs,vals,color=cols); axC.axhline(0,color="k",lw=.6)
    axC.set_ylabel("Spearman $\\rho$ with WBP5"); axC.tick_params(axis="x",labelsize=8)
    for b,v in zip(bars,vals):
        if np.isfinite(v): axC.text(b.get_x()+b.get_width()/2,v+0.004,f"{v:.3f}",ha="center",fontsize=8)
    for ax,lab in zip([axA,axB,axC],["(A)","(B)","(C)"]):
        ax.text(-0.16,1.06,lab,transform=ax.transAxes,fontsize=15,fontweight='bold',va='top',ha='left')
    fig.savefig(f"{a.outdir}/Fig_S5_updated_300dpi.png",dpi=300,bbox_inches="tight")
    print(f"\n저장 → {os.path.abspath(a.outdir)}\\Fig_S5_updated_300dpi.png")

if __name__=="__main__":
    main()
