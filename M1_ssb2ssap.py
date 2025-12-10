#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
M1_ssb2ssap.py
—————————————
从 scored_ctx_ctx0.json 提取 SSAP → CD-HIT → clusters.json / index.tsv / FASTA
最小改动增强：
- 在 Step1（extract_raw_ssap）从每个 ori_* 抽取 GCF/GCA（装配号），写入 FASTA 头“宿主块”的第4段
- parse_tmp_raw_fasta 同时兼容 3 段（无 GCF）与 4 段（含 GCF）
- build_index_from_maps 在 index.tsv 中新增 GCF 列

供 main.py 引用的函数：
  extract_raw_ssap, run_cdhit, parse_clstr_file, write_cluster_json,
  parse_tmp_raw_fasta, annotate_headers_with_cluster,
  reindex_from_annotated_to_minimal, build_index_from_maps
"""

from __future__ import annotations
import os, re, json, subprocess
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# ========== 小工具 ==========
_ALLOWED = re.compile(r"[^A-Za-z0-9._ ]+-")

def sanitize_label(s: str) -> str:
    if not s:
        return ""
    s = s.strip()
    s = _ALLOWED.sub("_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")

def pick_gcf(ori: dict) -> str:
    """
    尽可能从 ori_* 元数据里抽出装配号，返回 'GCF_xxx(.y)' 或 'GCA_xxx(.y)'。
    优先字段：Accession / GCF / GCA / assembly_accession / asm_acc / assembly / genome；
    回退：从任意字符串里用正则再扫一遍。
    """
    keys = ["Accession", "GCF", "GCA", "assembly_accession", "asm_acc", "assembly", "genome"]
    for k in keys:
        v = ori.get(k)
        if not v: 
            continue
        m = re.search(r"(GC[FA]_\d+(?:\.\d+)?)", str(v))
        if m:
            return m.group(1)
    for v in ori.values():
        m = re.search(r"(GC[FA]_\d+(?:\.\d+)?)", str(v))
        if m:
            return m.group(1)
    return ""

def append_refs(fasta_file: str):
    """可选：在 reindexed 的输出末尾追加参考序列（EcRecT / Redβ）"""
    with open(fasta_file, "a") as out_f:
        out_f.write(
""">EcRecT
MTKQPPIAKADLQKTQGNRAPAAVKNSDVISFINQPSMKEQLAAALPRHMTAERMIRIATTEIRKVPALGNCDTMSFVSAIVQCSQLGLEPGSALGHAYLLPFGNKNEKSGKKNVQLIIGYRGMIDLARRSGQIASLSARVVREGDEFSFEFGLDEKLIHRPGENEDAPVTHVYAVARLKDGGTQFEVMTRKQIELVRSLSKAGNNGPWVTHWEEMAKKTAIRRLFKYLPVSIEIQRAVSMDEKEPLTIDPADSSVLTGEYSVIDNSEE
>Redbet
MSTALATLAGKLAERVGMDSVDPQELITTLRQTAFKGDASDAQFIALLIVANQYGLNPWTKEIYAFPDKQNGIVPVVGVDGWSRIINENQQFDGMDFEQDNESCTCRIYRKDRNHPICVTEWMDECRREPFKTREGREITGPWQSHPKRMLRHKAMIQCARLAFGFAGIYDKDEAERIVENTAYTAERQPERDITPVNDETMQEINTLLIALDKTWDDDLLPLCSQIFRRDIRASSELTQAEAVKALGFLKQKAAEQKVAA
"""
        )

# ========== Step1：提取原始 SSAP 到临时 FASTA（头部携带宿主与 GCF） ==========
def extract_raw_ssap(json_file: str, fasta_out: str,
                     suffix: str = "FDDEIPF",
                     min_len: int = 200, max_len: int = 700) -> None:
    """
    从 JSON 提取 SSAP（长度过滤 + SSB尾巴匹配），写入临时 fasta。
    每条头部形如：
      >WP_XXXX|<status,strain,taxid,gcf>|<status,strain,taxid,gcf>|...
    其中 status: 1=host, 0=host_low（或 non_host 不写入）
    """
    with open(json_file, "r") as f:
        data = json.load(f)

    kept, skipped, filtered = 0, 0, 0
    with open(fasta_out, "w") as out_f:
        for ssap_id, ssap_data in data.items():
            ctx_ssap = ssap_data.get("ctx_ssap", {})
            ctx0 = ctx_ssap.get("ctx_0", [])
            if not ctx0:
                skipped += 1
                continue
            ssap_acc = ctx0[0].get("accession") or ssap_id
            ssap_seq = ctx0[0].get("sequence", "")
            if not ssap_seq:
                skipped += 1
                continue
            if not (min_len <= len(ssap_seq) <= max_len):
                filtered += 1
                continue

            # 收集各 ori_* 的宿主命中（优先 host，其次 host_low），写入四段式（末段是 GCF，可为空）
            host_blocks: List[str] = []
            for key, ori in ssap_data.items():
                if not str(key).startswith("ori_"):
                    continue
                strain = ori.get("strain", "Unknown")
                taxid  = ori.get("TaxID", "NA")
                gcf    = pick_gcf(ori)

                best = None  # (status_str '1'/'0', block_text)
                for ssb_key, ssb in ori.items():
                    if not str(ssb_key).startswith("ssb_"):
                        continue
                    status = ssb.get("status", "")
                    ctx0_s = (ssb.get("ctx_0") or [])
                    if not ctx0_s:
                        continue
                    ssb_seq = ctx0_s[0].get("sequence", "")
                    if not ssb_seq:
                        continue
                    # 尾巴匹配
                    if ssb_seq.endswith(suffix):
                        if status == "host":
                            best = ("1", f"1,{strain},{taxid},{gcf}")
                            break
                        elif status == "host_low":
                            best = best or ("0", f"0,{strain},{taxid},{gcf}")
                if best:
                    host_blocks.append(best[1])

            if host_blocks:
                kept += 1
                header = ">" + ssap_acc + "|" + "|".join(host_blocks)
                out_f.write(f"{header}\n{ssap_seq}\n")

    print(f"[M1] 提取完成：{kept} 条 → {fasta_out}（跳过缺失 {skipped}；长度过滤 {filtered}）")

# ========== Step2：CD-HIT ==========
def run_cdhit(in_fasta: str, out_prefix: str, c: float = 0.90, aL: float = 0.90) -> Tuple[str, str]:
    """
    运行 cd-hit；返回 (主输出fasta, .clstr) 的路径。
    """
    cmd = ["cd-hit", "-i", in_fasta, "-o", out_prefix, "-c", str(c), "-aL", str(aL), "-T", "0", "-M", "0"]
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_prefix, out_prefix + ".clstr"

# ========== Step3：解析 .clstr ==========
def parse_clstr_file(clstr_file: str, strict: bool = False):
    """
    只解析 accession ∈ {WP_, YP_, XP_} 的条目。
    返回：
      clusters[cid] = {"rep": (acc, identity or None, number),
                       "members": [(acc, identity or None, number), ...]}
      wp2info[acc] = (cid, number, role, similarity)
    """
    pat_acc = re.compile(r">(?:WP|YP|XP)_\d+(?:\.\d+)?")
    clusters: Dict[int, dict] = {}
    wp2info: Dict[str, Tuple[int,int,str,Optional[str]]] = {}
    cid = None
    skipped = 0

    with open(clstr_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cid = int(line.split()[1])
                clusters[cid] = {"rep": None, "members": []}
                continue
            # 条目序号
            mnum = re.match(r"^(\d+)\s+", line)
            number = int(mnum.group(1)) if mnum else 0
            # 只抓 WP/YP/XP
            macc = pat_acc.search(line)
            if not macc:
                skipped += 1
                if strict: 
                    raise ValueError(f"非 WP/YP/XP 行：{line}")
                else:
                    continue
            acc = macc.group(0)[1:]  # 去掉 >
            # 相似度（可能没有）
            msim = re.search(r"at\s+([\d\.]+)%", line)
            sim = msim.group(1) if msim else None
            if line.endswith("*"):
                clusters[cid]["rep"] = (acc, None, number)
                wp2info[acc] = (cid, number, "rep", None)
            else:
                clusters[cid]["members"].append((acc, sim, number))
                wp2info[acc] = (cid, number, "member", sim)

    if skipped:
        print(f"[M1] parse_clstr_file: 跳过 {skipped} 条非 WP/YP/XP 记录")
    return clusters, wp2info

# ========== Step4：生成 cluster JSON ==========
def write_cluster_json(clusters: dict, json_data: dict, out_json: str) -> None:
    """
    生成带 ctx/pfam 的 cluster JSON（键名用 CL{cid}）
    通过遍历 json_data 中的 ctx_ssap.*.*.accession == WP_XXXX 来回填上下文
    """
    def extract_ctx_for_wp(wp_id: str):
        for _, content in json_data.items():
            ctx_ssap = content.get("ctx_ssap", {})
            for ctx_key, ctx_entries in ctx_ssap.items():
                if isinstance(ctx_entries, list):
                    for entry in ctx_entries:
                        if isinstance(entry, dict) and entry.get("accession") == wp_id:
                            return {"ctx": ctx_ssap, "pfams": ctx_ssap.get("context_elements", [])}
        return None

    result = {}
    for cid, cl in clusters.items():
        entry = {"rep": None, "members": []}
        if cl.get("rep"):
            rep_wp, _, rep_no = cl["rep"]
            r = extract_ctx_for_wp(rep_wp)
            if r:
                entry["rep"] = {
                    "id": rep_wp,
                    "identity": None,
                    "number": rep_no,
                    "ctx": r["ctx"],
                    "pfams": r["pfams"]
                }
        for mwp, msim, mno in cl.get("members", []):
            r = extract_ctx_for_wp(mwp)
            if r:
                entry["members"].append({
                    "id": mwp,
                    "identity": msim,
                    "number": mno,
                    "ctx": r["ctx"],
                    "pfams": r["pfams"]
                })
        result[f"CL{cid}"] = entry

    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[M1] cluster JSON → {out_json}")

# ========== Step5：解析临时 FASTA（WP → 序列 / 宿主四元组） ==========
def parse_tmp_raw_fasta(fa_path: str):
    """
    解析 Step1 生成的临时 FASTA：
      >WP_XXXX | status,strain,taxid[,gcf] | status,strain,taxid[,gcf] ...
    返回：
      wp_to_seq:   WP -> AA 序列
      wp_to_hosts: WP -> [(status, strain, taxid, gcf), ...]   # gcf 可能为空串
    """
    wp_to_seq: Dict[str, str] = {}
    wp_to_hosts: Dict[str, List[Tuple[str,str,str,str]]] = defaultdict(list)
    wp, buf = None, []
    with open(fa_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if wp and buf:
                    wp_to_seq[wp] = "".join(buf); buf = []
                header = line.strip().lstrip(">")
                wp = header.split("|", 1)[0]
                for block in header.split("|")[1:]:
                    block = block.strip()
                    if not block:
                        continue
                    parts = [p.strip() for p in block.split(",")]
                    if len(parts) >= 3:
                        status, strain, taxid = parts[:3]
                        gcf = parts[3] if len(parts) >= 4 else ""
                        wp_to_hosts[wp].append((status, strain, taxid, gcf))
            else:
                buf.append(line.strip())
        if wp and buf:
            wp_to_seq[wp] = "".join(buf)
    return wp_to_seq, wp_to_hosts

# ========== Step6：写“带 [CL*] 标签”的 FASTA ==========
def annotate_headers_with_cluster(wp_to_seq: Dict[str,str],
                                  wp_to_hosts: Dict[str, List[Tuple[str,str,str,str]]],
                                  wp2info: Dict[str, Tuple[int,int,str,Optional[str]]],
                                  fa_all_annot: str,
                                  fa_rep_annot: str) -> None:
    """
    输出两份 FASTA：
      - 全量：>WP|host... [CL{cid}_{number}_rep] / [CL{cid}_{number}_{sim或NA}]
      - 仅 reps：>WP|host... [CL{cid}_{number}_rep]
    """
    # 全量
    with open(fa_all_annot, "w") as fa:
        for wp, seq in wp_to_seq.items():
            tag = ""
            if wp in wp2info:
                cid, number, role, sim = wp2info[wp]
                if role == "rep":
                    tag = f"[CL{cid}_{number}_rep]"
                else:
                    sim_s = sim if sim is not None else "NA"
                    tag = f"[CL{cid}_{number}_{sim_s}]"
            hosts = "|".join([",".join(h) for h in wp_to_hosts.get(wp, [])])
            header = f">{wp}"
            if hosts:
                header += "|" + hosts
            if tag:
                header += tag
            fa.write(f"{header}\n{seq}\n")

    # reps-only
    rep_set = {wp for wp, (_, _, role, _) in wp2info.items() if role == "rep"}
    with open(fa_rep_annot, "w") as fr:
        for wp in sorted(rep_set):
            if wp not in wp_to_seq:
                continue
            seq = wp_to_seq[wp]
            cid, number, _, _ = wp2info[wp]
            tag = f"[CL{cid}_{number}_rep]"
            hosts = "|".join([",".join(h) for h in wp_to_hosts.get(wp, [])])
            header = f">{wp}"
            if hosts:
                header += "|" + hosts
            header += tag
            fr.write(f"{header}\n{seq}\n")

    print(f"[M1] FASTA(带标签)：all → {fa_all_annot}  |  reps → {fa_rep_annot}")

# ========== Step7：极简 reindex（>CLx-y），可选追加参考 ==========
def reindex_from_annotated_to_minimal(in_fasta: str, out_fasta: str,
                                      wp2info: Dict[str, Tuple[int,int,str,Optional[str]]],
                                      append_refs_flag: bool = False) -> None:
    """
    把“带 CL 标签”的 fasta 头改成极简统一 ID：仅 `>CL{cid}-{number}`。
    """
    with open(in_fasta, "r") as fin, open(out_fasta, "w") as fout:
        wp = None
        seq_buf: List[str] = []
        for line in fin:
            if line.startswith(">"):
                if wp is not None:
                    if wp not in wp2info:
                        raise RuntimeError(f"未找到 WP={wp} 的 cluster 信息")
                    cid, number, _, _ = wp2info[wp]
                    fout.write(f">CL{cid}-{number}\n")
                    seq = "".join(seq_buf)
                    fout.write(seq + ("" if seq.endswith("\n") else "\n"))
                header = line.strip().lstrip(">")
                wp = header.split("|", 1)[0].split("[", 1)[0]
                seq_buf = []
            else:
                seq_buf.append(line.strip())
        if wp is not None:
            if wp not in wp2info:
                raise RuntimeError(f"未找到 WP={wp} 的 cluster 信息")
            cid, number, _, _ = wp2info[wp]
            fout.write(f">CL{cid}-{number}\n")
            seq = "".join(seq_buf)
            fout.write(seq + ("" if seq.endswith("\n") else "\n"))

    if append_refs_flag:
        append_refs(out_fasta)
    print(f"[M1] reindexed → {out_fasta}")

# ========== Step8：index.tsv ==========
def build_index_from_maps(wp_to_hosts: Dict[str, List[Tuple[str,str,str,str]]],
                          wp2info: Dict[str, Tuple[int,int,str,Optional[str]]],
                          index_out: str) -> None:
    """
    生成 index.tsv：
    列：ID, WP, GCF, Strain, TaxID, Status, CL, Number, Role, Similarity
    - ID = CL{cid}-{number}
    - CL = CL{cid}
    - GCF 来自宿主块第 4 段（可能为空）
    """
    with open(index_out, "w") as idx:
        idx.write("ID\tWP\tGCF\tStrain\tTaxID\tStatus\tCL\tNumber\tRole\tSimilarity\n")
        for wp, host_list in wp_to_hosts.items():
            if wp not in wp2info:
                continue
            cid, number, role, sim = wp2info[wp]
            sim = "" if (role == "rep" or sim is None) else str(sim)
            id_text = f"CL{cid}-{number}"
            cl_text = f"CL{cid}"

            if not host_list:
                idx.write(f"{id_text}\t{wp}\t\t\t\t\t{cl_text}\t{number}\t{role}\t{sim}\n")
            else:
                for rec in host_list:
                    status, strain, taxid = rec[:3]
                    gcf = rec[3] if len(rec) >= 4 else ""
                    clean_strain = sanitize_label(strain)
                    idx.write(
                        f"{id_text}\t{wp}\t{gcf}\t{clean_strain}\t{taxid}\t{status}\t{cl_text}\t{number}\t{role}\t{sim}\n"
                    )
    print(f"[M1] index.tsv → {index_out}")

# ========== 可选：本文件做脚本执行时的自测入口 ==========
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="M1: JSON→CD-HIT→FASTA/JSON/TSV（含 GCF 抽取）")
    ap.add_argument("-j", "--json", required=True, help="scored_ctx_ctx0.json")
    ap.add_argument("-o", "--outdir", required=True, help="输出目录（会创建）")
    ap.add_argument("-s", "--suffix", default="FDDEIPF", help="SSB 尾巴匹配后缀")
    ap.add_argument("-L", "--len-min", type=int, default=200)
    ap.add_argument("-U", "--len-max", type=int, default=700)
    ap.add_argument("-c", "--cdhit-c", type=float, default=0.90)
    ap.add_argument("-a", "--cdhit-aL", type=float, default=0.90)
    ap.add_argument("--append-refs", action="store_true", help="reindexed FASTA 末尾追加参考序列")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmp_raw = os.path.join(args.outdir, "tmp_raw_ssap.fasta")
    print(">> Step1")
    extract_raw_ssap(args.json, tmp_raw, suffix=args.suffix, min_len=args.len_min, max_len=args.len_max)

    print(">> Step2")
    cdhit_prefix = os.path.join(args.outdir, f"{args.suffix}_cdhit")
    _, clstr = run_cdhit(tmp_raw, cdhit_prefix, c=args.cdhit_c, aL=args.cdhit_aL)

    print(">> Step3")
    clusters, wp2info = parse_clstr_file(clstr)

    print(">> Step4")
    with open(args.json, "r") as jf:
        json_data = json.load(jf)
    clusters_json = os.path.join(args.outdir, f"{args.suffix}_ssap_clusters.json")
    write_cluster_json(clusters, json_data, clusters_json)

    print(">> Step5")
    wp_to_seq, wp_to_hosts = parse_tmp_raw_fasta(tmp_raw)

    print(">> Step6")
    fa_all = os.path.join(args.outdir, "raw_with_clusters.fasta")
    fa_rep = os.path.join(args.outdir, "rep_with_clusters.fasta")
    annotate_headers_with_cluster(wp_to_seq, wp_to_hosts, wp2info, fa_all, fa_rep)

    print(">> Step7")
    fa_all_reid = os.path.join(args.outdir, "raw_reid.fasta")
    fa_rep_reid = os.path.join(args.outdir, "rep_reid.fasta")
    reindex_from_annotated_to_minimal(fa_all, fa_all_reid, wp2info, append_refs_flag=args.append_refs)
    reindex_from_annotated_to_minimal(fa_rep, fa_rep_reid, wp2info, append_refs_flag=args.append_refs)

    print(">> Step8")
    index_tsv = os.path.join(args.outdir, "index.tsv")
    build_index_from_maps(wp_to_hosts, wp2info, index_tsv)

    print("[DONE]")