#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TaxID → 物种树构建与带注释叶标签渲染（兼容 CL 前缀；安全文本替换，不用解析器写回）
- 从 index.tsv 读取：ID(=CL{cid}-{num}), Strain, TaxID, Status, Role, Similarity
- 构建 NCBI 物种树（叶名=taxid），再把叶子的 taxid 文本替换为：
    CL3-1|'StrainA'[&status=1,c=rep]_'StrainB'[&status=0,c=97.3]_'+2'
- 额外的“member-only / 额外 rep 的 CL 块”折叠进外层注释属性：
    [&member_clusters_count=..., member_clusters='...']
    [&extra_rep_clusters_count=..., extra_rep_clusters='...']
  折叠文本内部的注释方括号会被转成全角［］，以便**保留 status/c 内容**且不干扰 Newick 解析。
"""

import os
import re
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple, Optional
from ete3 import NCBITaxa, Tree

# ====== 修改这里 ======
INDEX_TSV   = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/index.tsv"
OUT_TREE    = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/tree/taxid_tree.nwk"
OUT_LABELED = "/Users/lvdongyuan/Desktop/Bioinformatics/ETfinder/analyze_ssap_ssb/ETfinder2/tree/taxid_tree.labeled.nwk"

# 精简锚点（可选）
ANCHOR_TAXIDS = {
    562:    "⭐Escherichia coli",
    155892: "⭐Caulobacter crescentus",
    1718:   "⭐Corynebacterium glutamicum",
    1902:   "⭐Streptomyces coelicolor",
    817:    "⭐Bacteroides fragilis",
    1423:   "⭐Bacillus subtilis",
}

# —— 用户输入（可为空）：给“物种名/俗名/下划线形式/或 taxid 字符串”
USER_TARGET_INPUT: Optional[str] = None   # 例："Caulobacter crescentus" 或 "155892"

# 策略
ANCHOR_POLICY = 'skip_if_in_index'   # 或 'merge'
TARGET_POLICY = 'merge'              # 或 'skip_if_in_index'

# 展示控制
MAX_STRAINS_PER_CLUSTER = 3          # 每个 CL 展示的菌株条目上限
FALLBACK_SCIENTIFIC     = True       # 找不到映射时回退为学名

# 折叠策略
COLLAPSE_MEMBER_CLUSTERS_TO_ANNOT   = True   # 仅 member 的 CL 整段折叠进注释
COLLAPSE_EXTRA_REP_CLUSTERS_TO_ANNOT = True  # 额外 rep 的 CL 也折叠

# 仅用含 rep 的物种参与构树（否则 anchor/target 也会被纳入）
USE_REP_ONLY_FOR_TREE = True

def distance_map_from_raw_nwk(nwk_path: str, target_taxid: int | None) -> dict[int, float]:
    """
    读取 *原始* 物种树（叶名=TaxID 的 newick），返回 target_taxid 到所有叶的树距。
    目标不在树或为 None → 返回空 dict。
    """
    if not target_taxid:
        return {}
    t = Tree(nwk_path, format=1)
    leaves = {}
    for lf in t.iter_leaves():
        name = (lf.name or "").strip()
        if name.isdigit():
            leaves[int(name)] = lf
    if target_taxid not in leaves:
        return {}
    anchor = leaves[target_taxid]
    return {tid: t.get_distance(anchor, lf) for tid, lf in leaves.items()}

# ============ 小工具 ============
def ensure_parent_dir(p: str):
    d = os.path.dirname(os.path.abspath(p))
    if d:
        os.makedirs(d, exist_ok=True)

def quote_strain(s: str) -> str:
    """把菌株名包成 Newick 安全的 '...'（内部单引号转为弯引号）"""
    return "'" + str(s).replace("'", "’") + "'"

def sanitize_sim(sim: str) -> str:
    sim = (sim or "").strip()
    return re.sub(r"[^0-9.]", "", sim) if sim else ""

def quote_attr(s: str) -> str:
    """把注释里的属性值包成 '...'，内部单引号转为弯引号，去掉换行"""
    if s is None:
        return "''"
    s = str(s).replace("'", "’").replace("\n", " ").strip()
    return f"'{s}'"

def neutralize_ann_for_attr(s: str) -> str:
    """
    把内层 Newick 注释块 `[&key=val,...]` 转写成花括号 `{key=val,...}`，
    以便作为属性值安全嵌入，并保留全部内容。
    只做括号形状替换，不动内容与逗号/等号。
    """
    if not s:
        return s
    t = s.replace("\n", " ").strip()
    # 单层替换：把每个 `[& ... ]` 变成 `{ ... }`
    return re.sub(r"\[\&([^\[\]]*)\]", r"{\1}", t)



def escape_ann_for_attr(s: str) -> str:
    """
    折叠到属性值里的“内层文本块”安全化：
    - **保留**内部 [&status=...,c=...] 内容
    - 仅把 '[' ']' 替换为全角，避免外层注释误解析
    - 清理换行与首尾空白
    """
    if not s:
        return s
    t = s.replace("\n", " ").strip()
    return t.replace("[", "［").replace("]", "］")

# ============ 读取 index.tsv ============
def load_taxid_to_entries(index_file: str) -> Dict[str, OrderedDict]:
    """
    TaxID(str) -> OrderedDict(
        ID(str: CL{cid}-{num}) -> list[(strain, status, role, sim)]
    )
    """
    with open(index_file, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        for k in ("ID", "Strain", "TaxID", "Status"):
            if k not in cols:
                raise RuntimeError(f"index.tsv 缺少列：{k}")
        has_role = "Role" in cols
        has_sim  = "Similarity" in cols

        mapping: Dict[str, OrderedDict] = defaultdict(OrderedDict)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                taxid = str(int(parts[cols["TaxID"]]))
            except ValueError:
                continue

            cid    = parts[cols["ID"]].strip()      # e.g., CL3-1
            strain = parts[cols["Strain"]].strip()
            status = parts[cols["Status"]].strip()
            role   = parts[cols["Role"]].strip() if has_role else ""
            sim    = parts[cols["Similarity"]].strip() if has_sim else ""

            if not cid:
                continue
            if cid not in mapping[taxid]:
                mapping[taxid][cid] = []
            quad = (strain, status, role, sim)
            if strain and quad not in mapping[taxid][cid]:
                mapping[taxid][cid].append(quad)
    return mapping

# ============ 注入 anchor / target ============
def inject_anchors_and_target(mapping: Dict[str, OrderedDict],
                              user_taxid: Optional[int]) -> Dict[str, OrderedDict]:
    ncbi = NCBITaxa()
    anchor_map = dict(ANCHOR_TAXIDS)

    # 若用户目标与锚点重合：锚点让路
    if user_taxid is not None:
        anchor_map.pop(int(user_taxid), None)

    # 注入锚点（ID=物种名，role=anchor）
    for tid, label in anchor_map.items():
        s_tid = str(int(tid))
        already = s_tid in mapping and len(mapping[s_tid]) > 0
        if already and ANCHOR_POLICY == 'skip_if_in_index':
            continue
        if s_tid not in mapping:
            mapping[s_tid] = OrderedDict()
        if label not in mapping[s_tid]:
            mapping[s_tid][label] = []
        quad = ("", "", "anchor", "")
        if quad not in mapping[s_tid][label]:
            mapping[s_tid][label].append(quad)

    # 注入目标（ID=🎯Species，role=target）
    if user_taxid is not None:
        s_tid = str(int(user_taxid))
        try:
            sci = ncbi.get_taxid_translator([int(user_taxid)])[int(user_taxid)]
            label = f"🎯{sci}"
        except Exception:
            label = f"🎯taxid_{user_taxid}"

        already = s_tid in mapping and len(mapping[s_tid]) > 0
        if not (already and TARGET_POLICY == 'skip_if_in_index'):
            if s_tid not in mapping:
                mapping[s_tid] = OrderedDict()
            if label not in mapping[s_tid]:
                mapping[s_tid][label] = []
            quad = ("", "", "target", "")
            if quad not in mapping[s_tid][label]:
                mapping[s_tid][label].append(quad)

    return mapping

# ============ 选择构树 taxid ============
def collect_taxids_for_tree(mapping: Dict[str, OrderedDict],
                            rep_only: bool) -> List[int]:
    selected = set()
    for s_tid, entry_dict in mapping.items():
        roles = [ (r[2] or "").lower() for rows in entry_dict.values() for r in rows ]
        has_anchor_target = any(x in ("anchor","target") for x in roles)
        has_rep = any(x == "rep" for x in roles)
        if (not rep_only) or has_rep or has_anchor_target:
            try:
                selected.add(int(s_tid))
            except ValueError:
                pass
    return sorted(selected)

# ============ 构原始树 ============
def build_tree_from_taxids(taxids: List[int], out_file: str):
    ensure_parent_dir(out_file)
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(sorted(taxids))
    tree.write(format=1, outfile=out_file)
    print(f"[INFO] 原始物种树已保存: {out_file}（tips={len(taxids)}）")

# ============ 渲染单个 TaxID 的叶标签 ============
def build_label_for_taxid(entry_dict: OrderedDict,
                          max_strains: int = MAX_STRAINS_PER_CLUSTER) -> str:
    """
    生成单个 taxid 的叶标签（支持折叠“member-only”与“额外 rep”的 CL 块）。
    折叠文本里的内部注释从 `[&key=val,...]` 转写为 `{key=val,...}`，内容保留，避免破坏外层 Newick 注释。
    返回示例：
      CL3-1|'StrainA'[&status=1,c=rep]_'StrainB'[&status=0,c=97.3]_'+2'
        [&member_clusters_count=1,member_clusters='CL8-0|’S’{status=0,c=91.2} || ...']
        [&extra_rep_clusters_count=1,extra_rep_clusters='CL9-0|’S’{status=1,c=rep}']
    """
    # --- 把内层注释 `[&...]` 转写成 `{...}`，仅改括号形状，内容一字不丢 ---
    def neutralize_ann_for_attr(s: str) -> str:
        if not s:
            return s
        t = s.replace("\n", " ").strip()
        return re.sub(r"\[\&([^\[\]]*)\]", r"{\1}", t)

    # 1) 分离 anchor/target（总是可见）与普通 CL
    anchor_chunks = []
    normal_clusters = []
    for ident, rows in entry_dict.items():
        roles = {(r[2] or "").lower() for r in rows}
        if "anchor" in roles or "target" in roles:
            cval = "anchor" if "anchor" in roles else "target"
            anchor_chunks.append(f"{ident}[&c={cval}]")
        else:
            normal_clusters.append((ident, rows))

    # 2) CL 排序：含 rep > 最高相似度 > 条目数
    def cluster_score(rows):
        has_rep = any((r[2] or "").lower() == "rep" for r in rows)
        best_sim = max(
            (float(sanitize_sim(r[3])) if sanitize_sim(r[3]) else -1.0) for r in rows
        ) if rows else -1.0
        return (1 if has_rep else 0, best_sim, len(rows))

    normal_clusters_sorted = sorted(
        normal_clusters, key=lambda x: cluster_score(x[1]), reverse=True
    )

    # 3) 渲染单个 CL 文本块
    def render_cluster_block(cid: str, rows: list) -> tuple[str, bool]:
        """返回 (block_text, has_rep_in_this_taxid)"""
        def score_row(r):
            status, role, sim = r[1], (r[2] or "").lower(), r[3]
            is_rep = 1 if role == "rep" else 0
            s1 = 1 if str(status) == "1" else 0
            try:
                simv = float(sanitize_sim(sim)) if sanitize_sim(sim) else -1.0
            except Exception:
                simv = -1.0
            return (is_rep, s1, simv)

        rows_sorted = sorted(rows, key=score_row, reverse=True)
        shown = rows_sorted[:max_strains] if (max_strains and max_strains > 0) else rows_sorted
        hidden_rows = rows_sorted[len(shown):]

        bits = []
        has_rep_here = any((r[2] or "").lower() == "rep" for r in rows)
        for s, st, role, sim in shown:
            if not s:
                continue
            attrs = []
            if st != "":
                attrs.append(f"status={st}")
            if (role or "").lower() == "rep":
                attrs.append("c=rep")
            else:
                sim_clean = sanitize_sim(sim)
                attrs.append(f"c={sim_clean}" if sim_clean else "c=member")
            bits.append(f"{quote_strain(s)}[&{','.join(attrs)}]")

        if hidden_rows:
            hidden_names = [h[0] for h in hidden_rows if h and h[0]]
            bits.append(
                f"{quote_strain(f'+{len(hidden_rows)}')}[&fold=strain,hidden={quote_attr('|'.join(hidden_names))}]"
            )

        strain_part = "_".join(bits) if bits else "''"
        return f"{cid}|{strain_part}", has_rep_here

    # 4) 渲染所有普通 CL，区分含 rep 与 member-only
    rendered = [(cid,) + render_cluster_block(cid, rows)
                for cid, rows in normal_clusters_sorted]
    visible_rep_blocks   = [text for cid, text, has_rep in rendered if has_rep]
    hidden_member_blocks = [text for cid, text, has_rep in rendered if not has_rep] \
                           if COLLAPSE_MEMBER_CLUSTERS_TO_ANNOT else []

    # 5) 组装：anchor/target 在前，主块为第一个含 rep 的 CL，其它折叠进注释
    chunks: List[str] = []
    chunks.extend(anchor_chunks)

    if visible_rep_blocks:
        first = visible_rep_blocks[0]

        if hidden_member_blocks:
            merged_hidden = neutralize_ann_for_attr(" || ".join(hidden_member_blocks))
            first = (
                f"{first}[&member_clusters_count={len(hidden_member_blocks)},"
                f"member_clusters={quote_attr(merged_hidden)}]"
            )

        extra_reps = visible_rep_blocks[1:]
        if extra_reps:
            if COLLAPSE_EXTRA_REP_CLUSTERS_TO_ANNOT:
                merged_extra = neutralize_ann_for_attr(" || ".join(extra_reps))
                first = (
                    f"{first}[&extra_rep_clusters_count={len(extra_reps)},"
                    f"extra_rep_clusters={quote_attr(merged_extra)}]"
                )
                chunks.append(first)
            else:
                chunks.append(first)
                chunks.extend(extra_reps)
        else:
            chunks.append(first)
    else:
        # 没有 rep：兜底展示第一块
        if rendered:
            chunks.append(rendered[0][1])

    # anchor/target 和可见块之间用 " || " 分隔
    return " || ".join(chunks)

# ============ 文本级仅替换“叶子 taxid” ============
# 只匹配叶子位置的 token：
#  group(1) = 分隔符（行首"" 或 '(' 或 ','），保留在输出中
#  group(2) = 叶子名（这里是 taxid 数字）
LEAF_TOKEN_PAT = re.compile(r'(^|[\(,])(\d+)(?=[:),])')

def relabel_leaves_in_text(nwk_in: str,
                           nwk_out: str,
                           taxid_to_entries: Dict[str, OrderedDict],
                           max_strains: int = MAX_STRAINS_PER_CLUSTER,
                           fallback_scientific: bool = True):
    ensure_parent_dir(nwk_out)
    with open(nwk_in, "r", encoding="utf-8") as f:
        nwk_txt = f.read()

    # 先收集所有叶子 taxid，便于一次性翻译学名（可选）
    all_taxid_tokens = {m.group(2) for m in LEAF_TOKEN_PAT.finditer(nwk_txt)}
    sci = {}
    if fallback_scientific and all_taxid_tokens:
        ints = [int(x) for x in all_taxid_tokens if x.isdigit()]
        if ints:
            try:
                sci = NCBITaxa().get_taxid_translator(ints)
            except Exception:
                sci = {}

    known_ids = set(taxid_to_entries.keys())
    mapped = 0
    missed_ids = set()

    def _repl(m: re.Match) -> str:
        nonlocal mapped
        prefix = m.group(1)
        tid    = m.group(2)   # string
        if tid in known_ids:
            label = build_label_for_taxid(taxid_to_entries[tid], max_strains=max_strains)
            mapped += 1
            return prefix + label
        else:
            # 回退显示学名
            try:
                name = sci.get(int(tid))
            except Exception:
                name = None
            if name:
                return prefix + name
            missed_ids.add(tid)
            return prefix + tid

    labeled_txt = LEAF_TOKEN_PAT.sub(_repl, nwk_txt)
    with open(nwk_out, "w", encoding="utf-8") as f:
        f.write(labeled_txt)

    print(f"[OK] 替换完成：{nwk_out} (mapped={mapped}, missed={len(missed_ids)})")
    if missed_ids:
        example = ", ".join(list(missed_ids)[:8])
        more = f" ...(+{len(missed_ids)-8})" if len(missed_ids) > 8 else ""
        print(f"[NOTE] 这些 taxid 未出现在 index.tsv：{example}{more}")

# ============ 主程序 ============
if __name__ == "__main__":
    # 解析用户目标（可选）
    try:
        from sp2taxid import resolve_to_taxid  # 你本地的小工具
    except Exception:
        resolve_to_taxid = lambda s: None

    user_tid = resolve_to_taxid(USER_TARGET_INPUT) if USER_TARGET_INPUT else None
    if USER_TARGET_INPUT and not user_tid:
        print(f"[WARN] 无法解析用户输入：{USER_TARGET_INPUT}；将不加入目标物种。")

    # Step 1: 读取 index.tsv
    print("[STEP] 读取 index.tsv …")
    taxid_to_entries = load_taxid_to_entries(INDEX_TSV)

    # Step 2: 注入 anchor / target
    print("[STEP] 注入 anchor/target …")
    taxid_to_entries = inject_anchors_and_target(taxid_to_entries, user_tid)

    # Step 3: 选取参与构树的 taxid
    print("[STEP] 选择构树 taxid …")
    taxids_for_tree = collect_taxids_for_tree(taxid_to_entries, rep_only=USE_REP_ONLY_FOR_TREE)
    print(f"[INFO] 纳入构树的 taxid 数：{len(taxids_for_tree)} | 模式：{'rep-only' if USE_REP_ONLY_FOR_TREE else 'all'}")

    # Step 4: 构建原始 Newick（叶名=taxid）
    print("[STEP] 构建原始 Newick …")
    build_tree_from_taxids(taxids_for_tree, OUT_TREE)

    # Step 5: 文本级替换 taxid → 带注释标签
    print("[STEP] 渲染带注释叶标签 …")
    relabel_leaves_in_text(
        nwk_in=OUT_TREE,
        nwk_out=OUT_LABELED,
        taxid_to_entries=taxid_to_entries,
        max_strains=MAX_STRAINS_PER_CLUSTER,
        fallback_scientific=FALLBACK_SCIENTIFIC
    )

    print("[DONE] 原始树:", OUT_TREE)
    print("[DONE] 替换后:", OUT_LABELED)