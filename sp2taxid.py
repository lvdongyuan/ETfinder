#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sp2taxid.py
将任意“物种名/俗名/TaxID 字符串”解析为 NCBI TaxID。
- 输入可以是: "Escherichia coli", "E_coli", "562", "Bacillus_subtilis", 等
- 多义时优先选择 rank: species > strain > subspecies > genus
- 可作为库 import: resolve_to_taxid("Caulobacter crescentus") -> 155892
- 可命令行使用: python sp2taxid.py "Caulobacter crescentus"
"""

import sys
from typing import Optional, List
from ete3 import NCBITaxa

_PREF_ORDER = ["species", "strain", "subspecies", "genus"]

def _pick_best(taxids: List[int], ncbi: NCBITaxa) -> Optional[int]:
    if not taxids:
        return None
    ranks = ncbi.get_rank(taxids)
    best = None
    best_score = 10**9
    for tid in taxids:
        r = ranks.get(tid, "")
        try:
            score = _PREF_ORDER.index(r)
        except ValueError:
            score = 999
        if score < best_score:
            best_score = score
            best = tid
    return best if best is not None else taxids[0]

def resolve_to_taxid(q: Optional[str]) -> Optional[int]:
    """
    返回 int taxid 或 None
    """
    if not q:
        return None
    s = str(q).strip()
    if not s:
        return None
    # 数字直接返回
    if s.isdigit():
        return int(s)

    # 归一化：下划线→空格
    name = s.replace("_", " ").strip()

    ncbi = NCBITaxa()
    # 1) 原样
    try:
        m = ncbi.get_name_translator([name])
    except Exception:
        m = {}
    taxids = list(m.values())[0] if m else []

    # 2) 小写再试
    if not taxids:
        name2 = name.lower()
        try:
            m2 = ncbi.get_name_translator([name2])
            if m2:
                taxids = list(m2.values())[0]
        except Exception:
            pass

    # 3) 首字母大写（粗暴修正）
    if not taxids and " " in name:
        parts = name.split()
        name3 = parts[0].capitalize() + " " + " ".join(parts[1:])
        try:
            m3 = ncbi.get_name_translator([name3])
            if m3:
                taxids = list(m3.values())[0]
        except Exception:
            pass

    if not taxids:
        return None

    best = _pick_best(taxids, ncbi)
    return int(best) if best is not None else None

def main():
    if len(sys.argv) < 2:
        print("Usage: sp2taxid.py <species_name_or_taxid>")
        sys.exit(1)
    q = " ".join(sys.argv[1:])
    tid = resolve_to_taxid(q)
    if tid is None:
        print("NA")
        sys.exit(2)
    print(tid)

if __name__ == "__main__":
    main()