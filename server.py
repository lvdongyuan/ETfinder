#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ETfinder2 · Flask backend（稳、干净、最小依赖）
- 调度 main.py 管线：
    python main.py -s <suffix> --mode <rep|all> -t <target?> -c <cdhit_c> -m <IQTREE_MODEL>
- /task/<id> + /stream/<id>：SSE 实时日志
- /result/<id>：进化树 + 嵌入式 SSAP/SSB Context Viewer（clusters.json → Jinja）
- /files/*：仅开放 output/ 目录（防穿越）
- /clusters/<suffix>：独立 cluster 浏览页（可选）
- 若存在 index.sorted.tsv，则按 taxon_distance 对 clusters 排序
"""

from __future__ import annotations

import csv
import json
import queue
import shlex
import subprocess
import threading
import time
import uuid
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Iterable
from collections import defaultdict

from flask import (
    Flask, request, redirect, url_for, render_template, Response,
    send_from_directory, abort, jsonify
)

# 你的路径系统
from paths import build_paths, ensure_dirs

# ===== 基础配置 =====
WORK_DIR = Path(__file__).resolve().parent
PYEXE    = sys.executable
MAIN_PY  = str(WORK_DIR / "main.py")

app = Flask(
    __name__,
    template_folder=str(WORK_DIR / "templates"),
    static_folder=str(WORK_DIR / "static"),
)

# 任务登记表：task_id → 状态/参数/结果/线程
TASKS: Dict[str, Dict[str, Any]] = {}

# ====== 小工具 ======
def _make_task_id() -> str:
    return uuid.uuid4().hex[:12]

def _safe_float(v: Optional[str], default: float) -> float:
    try:
        return float(v) if v is not None and str(v).strip() != "" else float(default)
    except Exception:
        return float(default)

def _pick_index_for_ui(pp) -> Path:
    """优先 index.sorted.group.tsv → index.sorted.tsv → index.tsv。"""
    idx = Path(pp.index_tsv)
    grp = idx.with_name("index.sorted.group.tsv")
    srt = idx.with_name("index.sorted.tsv")
    if grp.exists():
        return grp
    if srt.exists():
        return srt
    return idx

def load_cluster_order_from_tsv(pp) -> list[str]:
    """
    从 index.sorted.group.tsv / index.sorted.tsv / index.tsv 里
    按行顺序抽取 CL 的顺序（去重保留第一次出现）。
    这反映了你在 main.py 里做好的最终排序。
    """
    tsv = _pick_index_for_ui(pp)
    if not tsv.exists():
        return []
    order: list[str] = []
    seen: set[str] = set()
    with open(tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            cl = (rec.get("CL") or "").strip()
            if not cl:
                continue
            if cl not in seen:
                seen.add(cl)
                order.append(cl)
    return order

def build_cl_distance_map(pp) -> Dict[str, float]:
    """
    读取 index(.sorted).tsv，构造 CL → 最小 taxon_distance。
    未知距离 → +inf（排序靠后）
    """
    cl2dist: Dict[str, float] = {}
    tsv = _pick_index_for_ui(pp)
    if not tsv.exists():
        return cl2dist
    with open(tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            cl = (rec.get("CL") or "").strip()
            if not cl:
                continue
            try:
                d = float(rec.get("taxon_distance", ""))
            except Exception:
                d = float("inf")
            if cl not in cl2dist or d < cl2dist[cl]:
                cl2dist[cl] = d
    return cl2dist

def build_wp2gcf(pp) -> Dict[str, list[str]]:
    """
    读取 index(.sorted).tsv，构造 WP → 去重后的 GCF 列表。
    兼容列名：优先 'GCF'；若无则尝试 'Accession'。
    """
    mapping: dict[str, set[str]] = defaultdict(set)
    tsv = _pick_index_for_ui(pp)
    if not tsv.exists():
        return {}
    with open(tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = r.fieldnames or []
        col_gcf = "GCF" if "GCF" in cols else ("Accession" if "Accession" in cols else None)
        if not col_gcf:
            return {}
        for rec in r:
            wp = (rec.get("WP") or "").strip()
            g  = (rec.get(col_gcf) or "").strip()
            if wp and g:
                mapping[wp].add(g)
    # 排序稳定
    return {wp: sorted(list(gs)) for wp, gs in mapping.items()}

def sort_clusters_by_distance(clusters: Dict[str, Any],
                              cl2dist: Dict[str, float]) -> Iterable[tuple[str, Any]]:
    """按与 target 的树距升序排序；未知距离排最后。"""
    return sorted(
        clusters.items(),
        key=lambda kv: (cl2dist.get(kv[0], float("inf")), kv[0])
    )

# ====== SSE 推流 ======
def _stream_sse_lines(q: "queue.Queue[str]", get_status, get_result):
    """把任务日志逐条推送到浏览器；附带 done / error 事件。"""
    yield "retry: 1000\n\n"  # 断线重连等待(ms)
    last_heartbeat = time.time()
    while True:
        try:
            line = q.get(timeout=0.5)
            yield f"data: {line.rstrip()}\n\n"
        except queue.Empty:
            now = time.time()
            if now - last_heartbeat > 5:
                yield ": heartbeat\n\n"
                last_heartbeat = now
        st = get_status()
        if st == "done":
            payload = json.dumps(get_result() or {})
            yield f"event: done\ndata: {payload}\n\n"
            break
        if st == "error":
            payload = json.dumps(get_result() or {"error": "unknown"})
            yield f"event: error\ndata: {payload}\n\n"
            break

# ====== 后台任务 ======
def _run_task_worker(task_id: str, suffix: str, mode: str, target: Optional[str]):
    t = TASKS[task_id]
    q: "queue.Queue[str]" = t["q"]

    # 输出路径（按 suffix）
    pp = build_paths(WORK_DIR, suffix=suffix)
    ensure_dirs(pp)

    # 取表单选项
    iq_model = t["params"].get("iqtree_model", "MFP")
    cdhit_c  = float(t["params"].get("cdhit_c", 0.90))
    ssap_k   = int(t["params"].get("ssap_k", 8))

    # 结果文件的约定路径（与 main.py 统一）
    ssap_prefix    = f"protein_tree.{'rep-only' if mode=='rep' else 'all'}"
    ssap_labeled   = pp.ssap_tree_dir / f"{ssap_prefix}.labeled.nwk"
    taxon_labeled  = pp.taxon_labeled_nwk
    clusters_json  = pp.clusters_json

    # 任务日志
    log_file = (pp.logs_dir / f"task_{task_id}.log")
    t["log_file"] = log_file
    log_file.parent.mkdir(parents=True, exist_ok=True)

    # ✅ 关键：--mode 传 rep/all；-m 传 IQ-TREE 模型
    cmd = [PYEXE, MAIN_PY,
           "-s", suffix,
           "--mode", mode,
           "-c", f"{cdhit_c:.4f}",
           "-m", iq_model,
           "--ssap-k", str(ssap_k)]
    if target:
        cmd += ["-t", target]

    # 把最终命令写回 params（用于 UI 展示/复现）
    t["params"]["cmd"] = " ".join(shlex.quote(x) for x in cmd)

    q.put(f"[INFO] starting pipeline: {' '.join(cmd)}")
    try:
        with open(log_file, "w", buffering=1) as lf:
            proc = subprocess.Popen(
                cmd, cwd=str(WORK_DIR),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, universal_newlines=True
            )
            for line in iter(proc.stdout.readline, ""):
                lf.write(line)
                q.put(line.rstrip("\n"))
            code = proc.wait()
            lf.write(f"\n[EXIT] returncode={code}\n")
            q.put(f"[EXIT] returncode={code}")

        if code != 0:
            t["status"] = "error"
            t["result"] = {"returncode": code, "log": str(log_file)}
            return

        # 汇总结果路径，传给模板
        index_sorted = pp.index_tsv.with_name("index.sorted.tsv")
        index_sorted_ctx = pp.index_tsv.with_name("index.sorted.ctx.tsv")
        index_sorted_group = pp.index_tsv.with_name("index.sorted.group.tsv")

        res = {
            "suffix": suffix,
            "mode": mode,
            "target": target or "",
            "clusters_json": str(clusters_json),
            "ssap_labeled": str(ssap_labeled),
            "taxon_labeled": str(taxon_labeled),
            "fasta_with_tags": str(pp.fasta_with_tags),
            "fasta_reindexed": str(pp.fasta_reindexed),
            "log": str(log_file),
            "cli": t["params"]["cmd"],
            "cdhit_c": cdhit_c,
            "iqtree_model": iq_model,
            "index_sorted": str(index_sorted),
            "index_sorted_ctx": str(index_sorted_ctx),
            # NEW:
            "index_sorted_group": str(index_sorted_group) if index_sorted_group.exists() else "",
            "ssap_k": ssap_k,
        }
        t["result"] = res
        t["status"] = "done"

    except Exception as e:
        t["status"] = "error"
        t["result"] = {"error": str(e), "log": str(log_file)}
        q.put(f"[ERROR] {e!r}")

# ====== Jinja 过滤器（viewer 用）======
@app.template_filter("wrap60")
def wrap60(s):
    import textwrap
    if not s:
        return ""
    return "\n".join(textwrap.wrap(s, 60))

@app.template_filter("dedupe_desc")
def dedupe_desc(desc, accession):
    if not desc:
        return ""
    acc = str(accession or "").strip()
    d = str(desc)
    if acc and d.startswith(acc):
        return d[len(acc):].lstrip()
    return d

# ====== 路由 ======
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

@app.route("/run", methods=["POST"])
def run_pipeline():
    suffix = (request.form.get("suffix") or "FDDEIPF").strip()
    target = (request.form.get("target") or "").strip() or None
    mode   = (request.form.get("mode") or "rep").strip()
    if mode not in ("rep", "all"):
        mode = "rep"

    # 只暴露 -c 与 IQ-TREE 模型
    cdhit_c = _safe_float(request.form.get("cdhit_c"), 0.90)
    model   = (request.form.get("iqtree_model") or "MFP").strip()
    if model not in {"MFP","JTT","LG","WAG","Dayhoff","VT","BLOSUM62"}:
        model = "MFP"
    ssap_k_raw = request.form.get("ssap_k", "").strip()
    try:
        ssap_k = int(ssap_k_raw) if ssap_k_raw else 8
    except ValueError:
        ssap_k = 8
    if ssap_k < 1:
        ssap_k = 1

    task_id = _make_task_id()
    TASKS[task_id] = {
        "status": "running",
        "q": queue.Queue(),
        "params": {
            "suffix": suffix,
            "target": target or "",
            "mode": mode,
            "cdhit_c": cdhit_c,
            "iqtree_model": model,
            "ssap_k": ssap_k,
            "cmd": "",  # 由 worker 填充最终命令
        },
        "result": None,
        "log_file": None,
        "thread": None,
    }

    th = threading.Thread(target=_run_task_worker, args=(task_id, suffix, mode, target), daemon=True)
    TASKS[task_id]["thread"] = th
    th.start()
    return redirect(url_for("task_page", task_id=task_id))

@app.route("/task/<task_id>")
def task_page(task_id):
    if task_id not in TASKS:
        abort(404)
    return render_template("task.html", task_id=task_id, params=TASKS[task_id]["params"])

@app.route("/tasks/<task_id>.json")
def task_json(task_id):
    t = TASKS.get(task_id)
    if not t:
        abort(404)
    return jsonify({
        "task_id": task_id,
        "status": t["status"],
        "params": t["params"],
        "result": t.get("result"),
        "log_file": str(t.get("log_file") or ""),
    })

@app.route("/stream/<task_id>")
def stream(task_id):
    if task_id not in TASKS:
        abort(404)
    t = TASKS[task_id]
    return Response(
        _stream_sse_lines(t["q"], lambda: t["status"], lambda: t["result"]),
        mimetype="text/event-stream"
    )

@app.route("/result/<task_id>")
def result_page(task_id):
    t = TASKS.get(task_id)
    if not t:
        abort(404)

    status = t["status"]
    result = t["result"] or {}

    # 读取预览（前 2000 字符）
    def _read_preview(p: Optional[str]) -> str:
        try:
            pth = Path(p) if p else None
            if pth and pth.exists():
                return pth.read_text(encoding="utf-8", errors="ignore")[:2000]
        except Exception:
            pass
        return ""

    ssap_preview  = _read_preview(result.get("ssap_labeled"))
    taxon_preview = _read_preview(result.get("taxon_labeled"))

    # 把绝对路径转成 /files/<relpath>
    out_root = (WORK_DIR / "output").resolve()
    def _href(abs_path: Optional[str]):
        if not abs_path:
            return None
        try:
            rel = str(Path(abs_path).resolve().relative_to(out_root))
            return url_for("files", subpath=rel)
        except Exception:
            return None

    ssap_url        = _href(result.get("ssap_labeled"))
    taxon_url       = _href(result.get("taxon_labeled"))
    log_url         = _href(result.get("log"))
    clusters_url    = _href(result.get("clusters_json"))
    fasta_tags_url  = _href(result.get("fasta_with_tags"))
    fasta_reid_url  = _href(result.get("fasta_reindexed"))
    index_sorted_url     = _href(result.get("index_sorted"))
    index_sorted_ctx_url = _href(result.get("index_sorted_ctx"))
    index_sorted_group_url = _href(result.get("index_sorted_group"))
    # 嵌入 viewer 的数据（即使没有也保证存在）
    clusters: Dict[str, Any] = {}
    clusters_pairs: Iterable[tuple[str, Any]] = []
    try:
        cj = result.get("clusters_json")
        if cj and Path(cj).exists():
            with open(cj, "r") as f:
                clusters = json.load(f)
        if clusters:
            try:
                pp = build_paths(WORK_DIR, suffix=result.get("suffix", ""))

                # ① 优先用 TSV 里的顺序（index.sorted.group.tsv → ...）
                order = load_cluster_order_from_tsv(pp)

                if order:
                    # 用 TSV 顺序构建 (cl_name, cl_data) 列表
                    order_set = set(order)
                    clusters_pairs = [
                        (cl_name, clusters[cl_name])
                        for cl_name in order
                        if cl_name in clusters
                    ]
                    # 把 TSV 里没出现的剩余 CL（例如某些没有 TaxID 的）放到最后（按名字排一下）
                    rest = [
                        (cl_name, data)
                        for cl_name, data in clusters.items()
                        if cl_name not in order_set
                    ]
                    rest.sort(key=lambda kv: kv[0])
                    clusters_pairs.extend(rest)
                else:
                    # ② 如果没有任何顺序信息，就退回到 distance 排序
                    cl2dist = build_cl_distance_map(pp)
                    clusters_pairs = list(sort_clusters_by_distance(clusters, cl2dist))
            except Exception:
                clusters_pairs = list(clusters.items())
    except Exception:
        clusters, clusters_pairs = {}, []

    # WP → GCF 列表，用于卡片显示 origin
    try:
        pp_for_map = build_paths(WORK_DIR, suffix=result.get("suffix", ""))
        wp2gcf = build_wp2gcf(pp_for_map)
    except Exception:
        wp2gcf = {}

    return render_template(
        "result.html",
        task_id=task_id,
        status=status,
        result=result,
        # 预览
        ssap_preview=ssap_preview,
        taxon_preview=taxon_preview,
        # 树文件 URL
        ssap_url=ssap_url,
        taxon_url=taxon_url,
        # 下载区
        log_url=log_url,
        clusters_url=clusters_url,
        fasta_tags_url=fasta_tags_url,
        fasta_reid_url=fasta_reid_url,
        # viewer 数据
        clusters=clusters,
        clusters_pairs=clusters_pairs,
        wp2gcf=wp2gcf, 
        index_sorted_url=index_sorted_url,
        index_sorted_ctx_url=index_sorted_ctx_url,
        index_sorted_group_url=index_sorted_group_url, # <--- 供模板显示 origin
    )

# ====== 独立 clusters 浏览（可选）======
@app.route("/clusters/<suffix>")
def clusters_page(suffix):
    pp = build_paths(WORK_DIR, suffix=suffix)
    if not pp.clusters_json.exists():
        abort(404)
    with open(pp.clusters_json, "r") as f:
        clusters = json.load(f)

    order = load_cluster_order_from_tsv(pp)
    if order:
        order_set = set(order)
        clusters_pairs = [
            (cl_name, clusters[cl_name])
            for cl_name in order
            if cl_name in clusters
        ]
        rest = [
            (cl_name, data)
            for cl_name, data in clusters.items()
            if cl_name not in order_set
        ]
        rest.sort(key=lambda kv: kv[0])
        clusters_pairs.extend(rest)
    else:
        cl2dist = build_cl_distance_map(pp)
        clusters_pairs = list(sort_clusters_by_distance(clusters, cl2dist))

    return render_template("ssap.html", clusters=clusters, clusters_pairs=clusters_pairs)

@app.route("/data/<suffix>")
def data_json(suffix):
    pp = build_paths(WORK_DIR, suffix=suffix)
    if not pp.clusters_json.exists():
        abort(404)
    with open(pp.clusters_json, "r") as f:
        data = json.load(f)
    return app.response_class(json.dumps(data), mimetype="application/json")

# ====== 静态文件下载（仅 output/ 下的文件）======
@app.route("/files/<path:subpath>")
def files(subpath):
    out_dir = (WORK_DIR / "output").resolve()
    file_path = (out_dir / subpath).resolve()
    if not str(file_path).startswith(str(out_dir)) or not file_path.exists():
        abort(404)
    return send_from_directory(out_dir, subpath, as_attachment=False)

@app.route("/health")
def health():
    return jsonify({"ok": True})

if __name__ == "__main__":
    # 外网访问可改 host="0.0.0.0"
    app.run(debug=False, port=5050, host="0.0.0.0")