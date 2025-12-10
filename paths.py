# paths.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

@dataclass
class PipelinePaths:
    work_dir: Path
    input_dir: Path
    output_dir: Path
    tmp_dir: Path          # 放 4 个阶段性 FASTA
    fasta_dir: Path        # 放最终两份 FASTA
    tree_dir: Path
    ssap_tree_dir: Path
    taxon_tree_dir: Path
    cdhit_dir: Path
    logs_dir: Path

    # 依赖 suffix 的文件（运行时填充）
    json_in: Path
    clusters_json: Path
    index_tsv: Path
    # 中间 FASTA（tmp 目录）
    fa_all_annot: Path
    fa_rep_annot: Path
    fa_all_reid: Path
    fa_rep_reid: Path
    # 最终 FASTA（fasta 目录）
    fasta_with_tags: Path
    fasta_reindexed: Path
    # Taxon tree 输出
    taxon_raw_nwk: Path
    taxon_labeled_nwk: Path
    # 你也可以在这里追加 SSAP tree 前缀，但我们在 main 里动态生成

def ensure_dirs(pp: PipelinePaths) -> None:
    for p in [
        pp.input_dir, pp.output_dir, pp.tmp_dir, pp.fasta_dir,
        pp.tree_dir, pp.ssap_tree_dir, pp.taxon_tree_dir,
        pp.cdhit_dir, pp.logs_dir
    ]:
        p.mkdir(parents=True, exist_ok=True)

def build_paths(
    work_dir: Path | str,
    suffix: str,
    json_in: Path | str | None = None,
) -> PipelinePaths:
    wd = Path(work_dir).resolve()
    input_dir  = wd / "input"
    output_dir = wd / "output"
    tmp_dir    = output_dir / "tmp"
    fasta_dir  = output_dir / "fasta"
    tree_dir   = output_dir / "tree"
    ssap_dir   = tree_dir / "ssap_tree"
    taxon_dir  = tree_dir / "taxon_tree"
    cdhit_dir  = wd / "cdhit"
    logs_dir   = wd / "logs"

    # 依 suffix 命名的输出
    clusters_json = output_dir / f"{suffix}_ssap_clusters.json"
    index_tsv     = output_dir / "index.tsv"

    # 4 个阶段性 FASTA（都在 tmp 里）
    fa_all_annot  = tmp_dir / "raw_with_clusters.fasta"
    fa_rep_annot  = tmp_dir / "rep_with_clusters.fasta"
    fa_all_reid   = tmp_dir / "raw_reid.fasta"
    fa_rep_reid   = tmp_dir / "rep_reid.fasta"

    # 最终名（在 fasta 里）
    fasta_with_tags = fasta_dir / f"host_{suffix}_ssap.fasta"
    fasta_reindexed = fasta_dir / f"host_{suffix}_ssap.reindexed.fasta"

    # taxon tree（固定名）
    taxon_raw_nwk     = taxon_dir / "taxid_tree.nwk"
    taxon_labeled_nwk = taxon_dir / "taxid_tree.labeled.nwk"

    return PipelinePaths(
        work_dir=wd,
        input_dir=input_dir,
        output_dir=output_dir,
        tmp_dir=tmp_dir,
        fasta_dir=fasta_dir,
        tree_dir=tree_dir,
        ssap_tree_dir=ssap_dir,
        taxon_tree_dir=taxon_dir,
        cdhit_dir=cdhit_dir,
        logs_dir=logs_dir,
        json_in=Path(json_in) if json_in else input_dir / "scored_ctx_ctx0.json",
        clusters_json=clusters_json,
        index_tsv=index_tsv,
        fa_all_annot=fa_all_annot,
        fa_rep_annot=fa_rep_annot,
        fa_all_reid=fa_all_reid,
        fa_rep_reid=fa_rep_reid,
        fasta_with_tags=fasta_with_tags,
        fasta_reindexed=fasta_reindexed,
        taxon_raw_nwk=taxon_raw_nwk,
        taxon_labeled_nwk=taxon_labeled_nwk,
    )