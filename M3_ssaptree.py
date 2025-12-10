#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein tree pipeline with one switch:
- USE_REP_INPUT=True  -> use rep-only reindexed FASTA (only reps on the tree)
- USE_REP_INPUT=False -> use full reindexed FASTA (reps + members on the tree)

Steps:
MAFFT --auto  ->  IQ-TREE3  ->  relabel leaves using index.tsv (text-safe, no ete3)
Output:
  <prefix>.<mode>.aln.fasta
  <prefix>.<mode>.treefile           (raw, names = Cluster{cid}x{num})
  <prefix>.<mode>.labeled.nwk        (labeled, with literal [&status=...,c=...] kept)
"""
import math
import os
import re
import subprocess
from collections import defaultdict
from typing import Dict, List, Tuple
from ete3 import Tree

# ======== 修改这里 ========
FASTA_ALL_REINDEXED = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/host_FDDEIPF_ssap.reindexed.fasta"   # 全量（rep+member）
FASTA_REP_REINDEXED = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/out_clustered/rep_reid.fasta"       # 仅 rep
INDEX_TSV           = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/index.tsv"
OUT_DIR             = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/tree"
OUT_PREFIX          = "protein_tree"
USE_REP_INPUT       = False   # <--- 核心开关：True=rep-only；False=全量
# 计算参数（按需改）
MAFFT_OPTS = ["--auto", "--thread", "-1"]
IQTREE_MODEL = os.environ.get("ET_IQTREE_MODEL", "JTT") 
IQTREE_BIN = "iqtree3"
IQTREE_OPTS  = ["-m", IQTREE_MODEL, "-B", "1000", "-nt", "AUTO", "-redo"]
# 标签展示：每个 ID 最多展示多少条宿主记录；超出折叠为 '+n'
MAX_STRAINS_PER_ID = 3
# =========================

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

# ---------- Newick 安全工具（我们自己负责“字面量”，不交给解析器） ----------
def quote_strain(s: str) -> str:
    # 用单引号包围菌株名，内部单引号替换为弯引号，避免破坏 Newick
    return "'" + str(s).replace("'", "’") + "'"

def sanitize_sim(sim: str) -> str:
    sim = (sim or "").strip()
    return re.sub(r"[^0-9.]", "", sim) if sim else ""

# ---------- 读 index.tsv ----------
def load_index(index_file: str) -> Dict[str, List[Tuple[str,str,str,str]]]:
    """
    返回：ID -> list of (Strain, Status, Role, Similarity)
    同一个 ID 可能多行（多宿主），保留全部用于标签渲染与折叠。
    """
    mapping: Dict[str, List[Tuple[str,str,str,str]]] = defaultdict(list)
    with open(index_file, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        for k in ("ID", "Strain", "Status", "Role", "Similarity"):
            if k not in cols:
                raise RuntimeError(f"index.tsv 缺列：{k}")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            id_  = parts[cols["ID"]].strip()
            if not id_:
                continue
            strain = parts[cols["Strain"]].strip()
            status = parts[cols["Status"]].strip()
            role   = parts[cols["Role"]].strip()
            sim    = parts[cols["Similarity"]].strip()
            if strain:
                tup = (strain, status, role, sim)
                if tup not in mapping[id_]:
                    mapping[id_].append(tup)
    return mapping

# ---------- 标签渲染（单叶=单个 ID） ----------
def build_label_for_id(id_: str,
                       rows: List[Tuple[str,str,str,str]],
                       max_strains: int = MAX_STRAINS_PER_ID) -> str:
    """
    目标：Cluster7x0|'StrainA'[&status=1,c=rep]_'StrainB'[&status=0,c=97.3]_'+2'
    - rep 节点：c=rep
    - member 节点：c=<Similarity>（若无则 c=member）
    """
    def score(r):
        s, st, role, sim = r
        is_rep = 1 if (role or "").lower() == "rep" else 0
        s1     = 1 if str(st) == "1" else 0
        try:
            simv = float(sanitize_sim(sim)) if sanitize_sim(sim) else -1.0
        except Exception:
            simv = -1.0
        return (is_rep, s1, simv)

    rows_sorted = sorted(rows, key=score, reverse=True)
    shown = rows_sorted[:max_strains] if (max_strains and max_strains > 0) else rows_sorted
    hidden = rows_sorted[len(shown):]

    bits = []
    for strain, st, role, sim in shown:
        attrs = []
        if st != "":
            attrs.append(f"status={st}")
        if (role or "").lower() == "rep":
            attrs.append("c=rep")
        else:
            sim_clean = sanitize_sim(sim)
            attrs.append(f"c={sim_clean}" if sim_clean else "c=member")
        attr_text = ",".join(attrs)
        bits.append(f"{quote_strain(strain)}[&{attr_text}]")

    if hidden:
        bits.append(quote_strain(f"+{len(hidden)}"))

    strain_part = "_".join(bits) if bits else "''"
    return f"{id_}|{strain_part}"

# ---------- MAFFT → IQ-TREE ----------
def run_mafft_iqtree(fasta_in: str, out_dir: str, out_prefix: str):
    ensure_dir(out_dir)
    aln = os.path.join(out_dir, out_prefix + ".aln.fasta")
    tree_prefix = os.path.join(out_dir, out_prefix)

    # MAFFT
    cmd_mafft = ["mafft"] + MAFFT_OPTS + [fasta_in]
    print("[RUN]", " ".join(cmd_mafft), ">", aln)
    with open(aln, "w") as fo:
        subprocess.run(cmd_mafft, check=True, stdout=fo)

    # IQ-TREE3
    cmd_iq = [IQTREE_BIN, "-s", aln] + IQTREE_OPTS + ["-pre", tree_prefix]
    print("[RUN]", " ".join(cmd_iq))
    subprocess.run(cmd_iq, check=True)

    treefile = tree_prefix + ".treefile"
    print(f"[OK] MSA: {aln}")
    print(f"[OK] Tree: {treefile}")
    return aln, treefile

# ---------- 逐叶“字面量”安全替换（不经解析器） ----------
def relabel_tree_leaves_textual(tree_in: str,
                                tree_out: str,
                                id2rows: Dict[str, List[Tuple[str, str, str, str]]]):
    """
    只在“叶子位点”把原始极简ID（如 CL3-2）替换为带注释的标签：
      CL3-2|'StrainA'[&status=1,c=rep]_'StrainB'[&status=0,c=97.3]_'+2'

    设计要点：
    - 纯“文本级”替换，不用解析器，避免 []=& 被库“消毒”。
    - 只匹配叶子 token：它们出现在 '(' 或 ',' 或 行首 之后，且在 ':' 或 ')' 或 ',' 之前。
      这样不会匹配到内部节点名（它们出现在 ')' 之后）。
    - 单次正则扫描，回调中按需替换，O(n)。
    - 统计命中/遗漏，便于核对。
    """
    # 读取原始 newick 文本
    with open(tree_in, "r", encoding="utf-8") as f:
        nwk_txt = f.read()

    # 统计用
    mapped = 0          # 成功替换的叶子次数（可能同一个ID出现多次）
    mapped_ids = set()  # 至少被替换过一次的ID集合
    missed_ids = set()  # 在树中出现但不在 id2rows 的ID集合

    # 预先做一个存在性集合，加速判断
    known_ids = set(id2rows.keys())

    # 只匹配叶子位置的 token：
    #  group(1) = 分隔符（行首"" 或 '(' 或 ','），保留在输出中
    #  group(2) = 叶子名（不含 : ( ) , 这些结构字符）
    leaf_token_pat = re.compile(r'(^|[\(,])([^:(),\s]+)(?=[:),])')

    def _repl(m: re.Match) -> str:
        nonlocal mapped
        prefix = m.group(1)
        name   = m.group(2)

        if name in known_ids:
            rows = id2rows[name]
            label = build_label_for_id(name, rows)
            mapped += 1
            mapped_ids.add(name)
            return prefix + label
        else:
            missed_ids.add(name)
            return prefix + name

    # 执行替换
    labeled_txt = leaf_token_pat.sub(_repl, nwk_txt)

    # 写出结果
    with open(tree_out, "w", encoding="utf-8") as f:
        f.write(labeled_txt)

    # 打印统计
    print(f"[OK] 替换完成：{tree_out}")
    print(f"     mapped_substitutions={mapped}, "
          f"unique_mapped_ids={len(mapped_ids)}, "
          f"unique_missed_ids={len(missed_ids)}")
    if missed_ids:
        # 只列一小部分，避免刷屏
        sample = ", ".join(list(missed_ids)[:10])
        more = "" if len(missed_ids) <= 10 else f" ...(+{len(missed_ids)-10})"
        print(f"     missed examples: {sample}{more}")

def _compute_leaf_distances(tree_path: str):
    """
    读取 Newick 树，返回 (leaf_names, 距离矩阵)，距离是树上的拓扑距离（patristic distance）。
    """
    t = Tree(tree_path, format=1)
    leaves = t.get_leaves()
    names = [lf.name for lf in leaves]
    n = len(leaves)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = t.get_distance(leaves[i], leaves[j])
            dist[i][j] = d
            dist[j][i] = d
    return names, dist


def _farthest_point_seeds(dist, k: int):
    """
    在距离矩阵上做 farthest-point sampling，选出 k 个“彼此尽量远”的种子索引。
    """
    n = len(dist)
    if n == 0:
        return []

    seeds = [0]  # 随便先选一个点做第一个 seed

    while len(seeds) < k and len(seeds) < n:
        best_idx = None
        best_min_d = -1.0
        for i in range(n):
            if i in seeds:
                continue
            dmin = min(dist[i][s] for s in seeds)
            if dmin > best_min_d:
                best_min_d = dmin
                best_idx = i
        if best_idx is None:
            break
        seeds.append(best_idx)

    return seeds


def _balanced_assign(dist, seeds, capacity: int):
    """
    在容量限制下，把所有点分配给最近的 seed：
    - 先把每个 seed 自己放进组里
    - 其他点按“先远后近”顺序，优先分配给最近且没满的 seed
    """
    n = len(dist)
    groups = {s: [s] for s in seeds}

    # 每个点到“最近种子”的距离，用来排序（离中心越远越优先被分配）
    min_dist_to_seed = []
    for i in range(n):
        dmin = min(dist[i][s] for s in seeds)
        min_dist_to_seed.append((dmin, i))
    min_dist_to_seed.sort(reverse=True)

    assigned = set(seeds)

    for _, i in min_dist_to_seed:
        if i in assigned:
            continue
        # 按“离这个点的距离”给种子排个序
        seed_order = sorted(seeds, key=lambda s: dist[i][s])
        placed = False
        for s in seed_order:
            if len(groups[s]) < capacity:
                groups[s].append(i)
                assigned.add(i)
                placed = True
                break
        # 理论上不会太常见：都满了时，硬塞进最近的那个
        if not placed:
            s0 = min(seeds, key=lambda s: dist[i][s])
            groups[s0].append(i)
            assigned.add(i)

    return groups


def build_balanced_groups_from_tree(treefile: str,
                                    k: int = 8,
                                    tsv_out: str | None = None) -> dict[str, str]:
    """
    基于 SSAP 树，对叶节点（CLx-y）做“距离+容量”平衡分组：
      - 在树距离空间里选 k 个种子叶子（尽量远）
      - 用容量约束的最近邻把其余叶子分到各组
    返回：{leaf_id: group_id}，group_id 为 G1..Gk
    若指定 tsv_out，则写一个两列的 TSV：leaf_id, group_id
    """
    names, dist = _compute_leaf_distances(treefile)
    n = len(names)
    if n == 0:
        return {}

    k_eff = min(k, n)
    capacity = math.ceil(n / k_eff)

    seeds = _farthest_point_seeds(dist, k_eff)
    groups_idx = _balanced_assign(dist, seeds, capacity)

    # 按组大小从大到小排序，给 G1..Gk 编号
    seed_list = list(groups_idx.keys())
    seed_list.sort(key=lambda s: len(groups_idx[s]), reverse=True)

    id_to_group: dict[str, str] = {}
    for gi, s in enumerate(seed_list, start=1):
        gid = f"G{gi}"
        for idx in groups_idx[s]:
            id_to_group[names[idx]] = gid

    if tsv_out:
        import csv
        with open(tsv_out, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["leaf_id", "group_id"])
            for leaf, gid in sorted(id_to_group.items()):
                w.writerow([leaf, gid])

    return id_to_group

# ---------- 主流程 ----------
def main():
    fasta_in = FASTA_REP_REINDEXED if USE_REP_INPUT else FASTA_ALL_REINDEXED
    mode = "rep-only" if USE_REP_INPUT else "all"
    if not os.path.exists(fasta_in):
        raise FileNotFoundError(f"{fasta_in} not found")
    ensure_dir(OUT_DIR)
    prefix = f"{OUT_PREFIX}.{mode}"

    # 1) MAFFT → IQ-TREE3
    aln, tree_raw = run_mafft_iqtree(fasta_in, OUT_DIR, prefix)

    # 2) 读 index → 文本级逐叶渲染标签
    id2rows = load_index(INDEX_TSV)
    tree_labeled = os.path.join(OUT_DIR, prefix + ".labeled.nwk")
    relabel_tree_leaves_textual(tree_raw, tree_labeled, id2rows)

    print("[DONE] 原始树:", tree_raw)
    print("[DONE] 替换树:", tree_labeled)

if __name__ == "__main__":
    main()