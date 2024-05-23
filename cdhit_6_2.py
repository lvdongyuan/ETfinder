import subprocess

def run_cd_hit(input_fasta, output_fasta, c=0.95):
    """
    使用 CD-HIT 对输入的 FASTA 文件进行聚类
    :param input_fasta: 输入的 FASTA 文件路径
    :param output_fasta: 输出的 FASTA 文件路径
    :param c: 序列相似性阈值，默认是 0.9
    """
    cd_hit_cmd = [
        'cd-hit',
        '-i', input_fasta,
        '-o', output_fasta,
        '-c', str(c)  # 设置序列相似性阈值
    ]
    
    try:
        # 执行 CD-HIT 命令
        subprocess.run(cd_hit_cmd, check=True)
        print(f"CD-HIT clustering completed. Results saved to {output_fasta}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running CD-HIT. Error message: {e}")
