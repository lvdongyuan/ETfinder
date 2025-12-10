#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""


用法示例：
  # default rep 
  python main.py

  # all model
  python main.py -m all

  # -suffix
  python main.py -w ETfinder2 -j input/scored_ctx_ctx0.json -s FDDEIPF -m rep -t "Bacillus subtilis"
"""

import argparse
import json
import os
import shutil
import logging
from pathlib import Path
import csv



from paths import build_paths, ensure_dirs, PipelinePaths


from M1_ssb2ssap import (
    extract_raw_ssap, run_cdhit, parse_clstr_file, write_cluster_json,
    parse_tmp_raw_fasta, annotate_headers_with_cluster,
    reindex_from_annotated_to_minimal, build_index_from_maps
)
from M2_taxontree import (
    load_taxid_to_entries, inject_anchors_and_target,
    collect_taxids_for_tree, build_tree_from_taxids, relabel_leaves_in_text, distance_map_from_raw_nwk
)
from M3_ssaptree import (
    run_mafft_iqtree, load_index as load_index_idmap,
    relabel_tree_leaves_textual, build_balanced_groups_from_tree
)

try:
    from sp2taxid import resolve_to_taxid
except Exception:
    resolve_to_taxid = lambda s: None


def check_binaries(bins: list[str]):
    missing = [b for b in bins if shutil.which(b) is None]
    if missing:
        logging.warning("未找到外部程序：%s（相关步骤将失败）", ", ".join(missing))

def move(src: Path, dst: Path):
    src = Path(src); dst = Path(dst)
    if not src.exists():
        return
    if src.resolve() == dst.resolve():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.replace(str(src), str(dst))


def run_pipeline(
    pp: PipelinePaths,
    *,
    mode: str = "rep",                 
    target: str | None = None,         
    len_min: int = 200,
    len_max: int = 700,
    cdhit_c: float = 0.90,
    cdhit_aL: float = 0.90,
    mafft_opts: list[str] | None = None,
    iqtree_bin: str = "iqtree3",
    iqtree_opts: list[str] | None = None,
    append_refs_flag: bool = True,
    ssap_k: int = 8,
):
    """
      - mode='rep'
      - mode='all'
    """
    ensure_dirs(pp)

    if mafft_opts is None:
        mafft_opts = ["--auto", "--thread", "-1"]
    if iqtree_opts is None:
        iqtree_opts = ["-m", "JTT", "-B", "1000", "-nt", "AUTO", "-redo"]
    if mode not in ("rep", "all"):
        raise ValueError("--mode must 'rep' or 'all'")

    check_binaries(["cd-hit", "mafft", iqtree_bin])

    # ===== M1：extract → CD-HIT → JSON/FASTA/TSV =====
    tmp_raw = pp.cdhit_dir / "tmp_raw_ssap.fasta"
    cdhit_prefix = pp.cdhit_dir / f"{pp.clusters_json.stem.split('_')[0]}_cdhit"  # suffix 

    logging.info("Step1 提取 SSAP → %s", tmp_raw)
    extract_raw_ssap(str(pp.json_in), str(tmp_raw),
                     suffix=pp.clusters_json.stem.split('_')[0],  # e.g. FDDEIPF
                     min_len=len_min, max_len=len_max)

    logging.info("Step2 CD-HIT → %s.*", cdhit_prefix)
    _, clstr_file = run_cdhit(str(tmp_raw), str(cdhit_prefix), c=cdhit_c, aL=cdhit_aL)
    clusters, wp2info = parse_clstr_file(str(clstr_file))

    logging.info("Step3 cluster JSON（ctx/pfam）→ %s", pp.clusters_json)
    with open(pp.json_in, "r") as jf:
        json_data = json.load(jf)
    write_cluster_json(clusters, json_data, str(pp.clusters_json))

    logging.info("Step4 Read temporary FASTA and build sequence/host mapping")
    wp_to_seq, wp_to_hosts = parse_tmp_raw_fasta(str(tmp_raw))

    logging.info("Step5 Write FASTA with [CL*] labels (full set & representatives) → %s", pp.tmp_dir)
    annotate_headers_with_cluster(
        wp_to_seq, wp_to_hosts, wp2info,
        fa_all_annot=str(pp.fa_all_annot),
        fa_rep_annot=str(pp.fa_rep_annot)
    )

    logging.info("Step6 Minimal reindexing (>CLx-y) and append references (optional) → %s", pp.tmp_dir)
    reindex_from_annotated_to_minimal(str(pp.fa_all_annot), str(pp.fa_all_reid), wp2info, append_refs_flag=append_refs_flag)
    reindex_from_annotated_to_minimal(str(pp.fa_rep_annot), str(pp.fa_rep_reid), wp2info, append_refs_flag=append_refs_flag)

   
    move(pp.fa_all_annot, pp.fasta_with_tags)
    move(pp.fa_all_reid,  pp.fasta_reindexed)

    logging.info("Step7 index.tsv → %s", pp.index_tsv)
    build_index_from_maps(wp_to_hosts, wp2info, str(pp.index_tsv))

    # M3：SSAP tree
    logging.info("Step8 SSAP tree（mode=%s） → %s", mode, pp.ssap_tree_dir)
    fasta_for_tree = str(pp.fa_rep_reid if mode == "rep" else pp.fasta_reindexed)
    out_prefix = f"protein_tree.{ 'rep-only' if mode=='rep' else 'all' }"
    aln, tree_raw = run_mafft_iqtree(
        fasta_in=fasta_for_tree,
        out_dir=str(pp.ssap_tree_dir),
        out_prefix=out_prefix
    )
    id2rows = load_index_idmap(str(pp.index_tsv))
    ssap_tree_labeled = pp.ssap_tree_dir / (out_prefix + ".labeled.nwk")
    relabel_tree_leaves_textual(tree_raw, str(ssap_tree_labeled), id2rows)
    cl2group: dict[str, str] = {}
    if mode == "rep":
        groups_tsv = pp.ssap_tree_dir / (out_prefix + ".groups.tsv")
        cl2group = build_balanced_groups_from_tree(
            treefile=tree_raw,
            k=ssap_k,  # 你现在想要 8 个分支
            tsv_out=str(groups_tsv),
        )
        logging.info(
            "Step8b SSAP tree balanced groups: %d groups → %s",
            len(set(cl2group.values())),
            ssap_k,
            groups_tsv,
        )
    else:
        logging.info("Step8b Skip SSAP grouping for mode='all'")
    
      
    # M2：Taxon tree
    logging.info("Step9 Taxon tree（mode=%s） → %s", mode, pp.taxon_tree_dir)
    taxid_to_entries = load_taxid_to_entries(str(pp.index_tsv))
    user_tid = resolve_to_taxid(target) if target else None
    taxid_to_entries = inject_anchors_and_target(taxid_to_entries, user_tid)

    taxids_for_tree = collect_taxids_for_tree(taxid_to_entries, rep_only=(mode == "rep"))
    build_tree_from_taxids(taxids_for_tree, str(pp.taxon_raw_nwk))
    relabel_leaves_in_text(
        nwk_in=str(pp.taxon_raw_nwk),
        nwk_out=str(pp.taxon_labeled_nwk),
        taxid_to_entries=taxid_to_entries,
        max_strains=3,
        fallback_scientific=True
    )
    return cl2group

def write_index_sorted_by_distance(index_in: Path, index_out: Path, taxid_to_dist: dict[int, float]):
   
    rows = []
    with open(index_in, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = r.fieldnames or []
        out_fields = fieldnames + ([ "taxon_distance" ] if "taxon_distance" not in fieldnames else [])
        for rec in r:
            try:
                tid = int(rec.get("TaxID","") or 0)
            except ValueError:
                tid = 0
            d = taxid_to_dist.get(tid, "")
            rec["taxon_distance"] = d
            rows.append(rec)

    def sort_key(rec):
        d = rec.get("taxon_distance")
        
        return (0 if isinstance(d, (int, float)) else 1, d if d != "" else float("inf"))

    rows.sort(key=sort_key)

    with open(index_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
def write_index_group_roundrobin(index_in: Path,
                                 index_out: Path,
                                 id_to_group: dict[str, str],
                                 group_col: str = "ssap_group") -> None:
    """
    在 index_in（已按 taxon_distance 排好）基础上：
      - 为每行增加一列 group_col = G1..Gk（没有分组的留空）
      - 对有分组的行按组内 taxon_distance 先排好
      - 按“分层轮询”的方式输出：
          第 0 层：各组的第 1 条，按 distance 从小到大排
          第 1 层：各组的第 2 条，按 distance 从小到大排
          ...
      - 没有分组的行追加在最后
    """
    rows_by_group: dict[str, list[dict]] = {}
    other_rows: list[dict] = []

    with open(index_in, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = r.fieldnames or []
        for rec in r:
            ident = (rec.get("ID") or "").strip()
            gid = id_to_group.get(ident, "")
            rec[group_col] = gid
            if gid:
                rows_by_group.setdefault(gid, []).append(rec)
            else:
                other_rows.append(rec)

    # 组的顺序：按 G 后面的数字排序（G1,G2,...）；其他名字排在后面
    def _group_sort_key(g: str) -> int:
        if g.startswith("G"):
            try:
                return int(g[1:])
            except ValueError:
                return 9999
        return 9999

    ordered_groups = sorted(rows_by_group.keys(), key=_group_sort_key)

    # 每个组内部按 taxon_distance 排序（越近越靠前）
    def _dist_key(rec: dict) -> float:
        d = rec.get("taxon_distance")
        if d == "" or d is None:
            return float("inf")
        try:
            return float(d)
        except Exception:
            return float("inf")

    for g in ordered_groups:
        rows_by_group[g].sort(key=_dist_key)

    # 分层轮询：
    # 第 i 层：从每个组取第 i 条记录（若存在），收集起来按 distance 排，再输出
    out_rows: list[dict] = []
    max_len = max((len(v) for v in rows_by_group.values()), default=0)

    for layer_idx in range(max_len):
        layer_rows = []
        for g in ordered_groups:
            glist = rows_by_group[g]
            if layer_idx < len(glist):
                layer_rows.append(glist[layer_idx])
        layer_rows.sort(key=_dist_key)
        out_rows.extend(layer_rows)

    # 最后再把无分组的行加到结尾
    out_rows.extend(other_rows)

    # 确保 group_col 在字段列表里
    with open(index_in, newline="") as f:
        r0 = csv.DictReader(f, delimiter="\t")
        fieldnames = r0.fieldnames or []
    if group_col not in fieldnames:
        fieldnames = fieldnames + [group_col]

    with open(index_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

def _wrap60(seq: str, width: int = 60) -> str:
    seq = (seq or "").replace("\n", "").strip()
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width)) if seq else ""

def _first_fasta_block(ctx: dict, key: str) -> str:
    """
    取 ctx[key] 的第一条，格式化为一条标准 FASTA：
      >accession description\nSEQUENCE(60列换行)
    若无则返回空串。
    """
    if not ctx or key not in ctx or not ctx[key]:
        return ""
    e = ctx[key][0]  
    acc = (e.get("accession") or "NA").strip()
    desc = (e.get("description") or "").strip()
    seq = _wrap60(e.get("sequence") or "")
    head = f">{acc} {desc}".rstrip()
    return f"{head}\n{seq}" if seq else head

def _build_ctx_map_from_clusters(clusters_json_path: Path) -> dict[tuple[str, int], dict]:
    """
    返回 {(CL, number): ctx_dict}，包含 rep 与 members。
    """
    with open(clusters_json_path, "r") as f:
        clusters = json.load(f)
    m: dict[tuple[str,int], dict] = {}
    for cl_name, cl in clusters.items():
        rep = cl.get("rep") or {}
        if rep.get("ctx") is not None:
            # rep.number 可能不存在，稳妥起见当成 0
            rn = rep.get("number")
            try:
                rn = int(rn) if rn is not None else 0
            except Exception:
                rn = 0
            m[(cl_name, rn)] = rep["ctx"]
        for mem in cl.get("members", []) or []:
            if mem.get("ctx") is None:
                continue
            num = mem.get("number")
            try:
                num = int(num)
            except Exception:
                continue
            m[(cl_name, num)] = mem["ctx"]
    return m

def _pick_input_tsv_for_ctx(index_sorted: Path, index_raw: Path) -> Path:
    return index_sorted if index_sorted.exists() else index_raw

def write_index_sorted_with_ctx(index_in: Path, index_out: Path, ctx_map: dict[tuple[str,int], dict]):
    """
    在 index_in 基础上追加 ctx_-2 ~ ctx_+2 五列，写出 index_out。
    依赖列：优先用 (CL, Number)；若缺失则从 ID 解析 CLx-y。
    """
    rows = []
    with open(index_in, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = r.fieldnames or []
        add_cols = ["ctx_-2", "ctx_-1", "ctx_0", "ctx_+1", "ctx_+2"]
        out_fields = fields + [c for c in add_cols if c not in fields]
        for rec in r:
            cl = (rec.get("CL") or "").strip()
            num = rec.get("Number")
            # 回退：从 ID 解析
            if not cl or num in (None, ""):
                ident = (rec.get("ID") or "").strip()
                
                import re
                m = re.match(r"^(CL\d+)[\-x](\d+)$", ident)
                if m:
                    cl = m.group(1)
                    num = m.group(2)
            try:
                num_i = int(num)
            except Exception:
                num_i = None

            ctx = ctx_map.get((cl, num_i), {}) if (cl and num_i is not None) else {}

            # 只取唯一的一条 FASTA（若无则空）
            rec["ctx_-2"] = _first_fasta_block(ctx, "ctx_-2")
            rec["ctx_-1"] = _first_fasta_block(ctx, "ctx_-1")
            rec["ctx_0"]  = _first_fasta_block(ctx, "ctx_0")
            rec["ctx_+1"] = _first_fasta_block(ctx, "ctx_+1")
            rec["ctx_+2"] = _first_fasta_block(ctx, "ctx_+2")
            rows.append(rec)

    with open(index_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)   

# ---------- CLI ----------
def parse_args():
    ap = argparse.ArgumentParser(
        prog="ETfinder2",
        description="SSAP/SSB pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-w", "--work", default=".", help="工作目录")
    ap.add_argument("-j", "--json", default=None, help="输入 JSON")
    ap.add_argument("-s", "--suffix", default="FDDEIPF", help="SSB 尾巴")
    ap.add_argument("--mode", choices=["rep", "all"], default="rep",
                    help="数据模式：rep=只用代表；all=全量")
    ap.add_argument("-t", "--target", default=None, help="Taxon 目标（可选）")
    ap.add_argument("-L", "--len-min", type=int, default=200)
    ap.add_argument("-U", "--len-max", type=int, default=700)
    ap.add_argument("-c", "--cdhit-c", type=float, default=0.90)
    ap.add_argument("-a", "--cdhit-aL", type=float, default=0.90)
    ap.add_argument("-m", "--model", default="JTT",
                    help="IQ-TREE substitution model (JTT, LG, WAG, MFP, ...)")
    ap.add_argument("--iqtree-bin", default="iqtree3")
    ap.add_argument(
        "--ssap-k", type=int, default=8,
        help="Target number of RecT candidates for stratified sampling"
    )
    return ap.parse_args()

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S"
    )
    args = parse_args()                           
    os.environ["ET_IQTREE_MODEL"] = args.model    

    pp = build_paths(work_dir=Path(args.work),
                     suffix=args.suffix,
                     json_in=Path(args.json) if args.json else None)
    cl2group = run_pipeline(
        pp,
        mode=args.mode,               
        target=args.target,
        len_min=args.len_min,
        len_max=args.len_max,
        cdhit_c=args.cdhit_c,
        cdhit_aL=args.cdhit_aL,
        iqtree_bin=args.iqtree_bin,
        ssap_k=args.ssap_k,
    )

    
    user_tid = resolve_to_taxid(args.target) if args.target else None
    dist_map = distance_map_from_raw_nwk(str(pp.taxon_raw_nwk), user_tid)
    index_sorted = pp.index_tsv.with_name("index.sorted.tsv")
    write_index_sorted_by_distance(pp.index_tsv, index_sorted, dist_map)
    logging.info("  - TSV(sorted): %s", index_sorted)
    ctx_map = _build_ctx_map_from_clusters(pp.clusters_json)
    index_in_for_ctx = _pick_input_tsv_for_ctx(index_sorted, pp.index_tsv)
    index_ctx = pp.index_tsv.with_name("index.sorted.ctx.tsv")
    write_index_sorted_with_ctx(index_in_for_ctx, index_ctx, ctx_map)
    logging.info("  - TSV(sorted+ctx): %s", index_ctx)
        # NEW: 如果有 SSAP 分组信息，再做一份“组轮询排序”的 TSV
    if cl2group:
        index_group = pp.index_tsv.with_name("index.sorted.group.tsv")
        write_index_group_roundrobin(index_sorted, index_group, cl2group, group_col="ssap_group")
        logging.info("  - TSV(sorted+group_rr): %s", index_group)
    else:
        logging.info("  - No SSAP groups (mode != 'rep'); skip group-based reordering")
if __name__ == "__main__":
    main()
    