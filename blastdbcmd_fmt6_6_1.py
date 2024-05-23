import os
import subprocess

# 指定ID列表的文件路径
id_list_path = '/Volumes/ex/recETresults/paired_recT_ids.txt'
# 指定BLAST数据库的路径
db_path = '/Volumes/ex/recETresults/lddeipfdb/lddeipf'
# 输出FASTA文件的路径
fasta_output_path = '/Volumes/ex/recETresults/paired_recTseq.txt'

# 确保输出文件所在目录存在
os.makedirs(os.path.dirname(fasta_output_path), exist_ok=True)

# 读取ID列表
with open(id_list_path, 'r') as id_file:
    seq_ids = [line.strip() for line in id_file if line.strip()]

# 用blastdbcmd从数据库中提取序列
with open(fasta_output_path, 'w') as fasta_out:
    for seq_id in seq_ids:
        print(f"Retrieving sequence for ID: {seq_id}")  # 打印当前处理的序列ID
        blastdbcmd_cmd = [
            'blastdbcmd',
            '-db', db_path,
            '-entry', seq_id,
            '-outfmt', '%f',
            '-line_length', '9999'  # 设置大行长度避免换行
        ]
        proc = subprocess.Popen(blastdbcmd_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.communicate()
        
        if proc.returncode == 0:
            fasta_out.write(stdout.strip() + '\n')
        else:
            print(f"Error retrieving {seq_id}: {stderr}")

print(f"Sequences have been extracted to {fasta_output_path}")
