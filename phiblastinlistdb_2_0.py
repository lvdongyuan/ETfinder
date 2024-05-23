import os
import subprocess

def blast_and_extract(db_directory, query_file, blast_results_file, no_hits_file, patternfile, max_hits):
    
    if not os.path.exists(os.path.dirname(blast_results_file)):
        os.makedirs(os.path.dirname(blast_results_file))
    if not os.path.exists(os.path.dirname(no_hits_file)):
        os.makedirs(os.path.dirname(no_hits_file))

    db_files = [f for f in os.listdir(db_directory) if f.endswith('.pin')]
    db_paths = [os.path.join(db_directory, os.path.splitext(f)[0]) for f in db_files]

    total_dbs = len(db_paths)
    dbs_with_no_hits = 0

    with open(blast_results_file, 'w') as blast_out, open(no_hits_file, 'w') as no_hits_out:
        for db_path in db_paths:
            print(f"Processing database: {db_path.split('/')[-1]}")
            blast_cmd = [
                'psiblast',
                '-evalue', '0.005',
                '-query', query_file,
                '-db', db_path,
                '-outfmt', '6 sseqid sseq',
                '-max_target_seqs', str(max_hits),
                '-phi_pattern', patternfile
            ]
            result = subprocess.run(blast_cmd, text=True, capture_output=True, check=True)
            if result.stdout:
                # 过滤空行和不包含制表符的行
                lines = [line for line in result.stdout.strip().split('\n') if '\t' in line]

                blast_out.write(f"Results for database: {db_path.split('/')[-1]}\n")
                blast_out.write(result.stdout)
                blast_out.write('-' * 40 + '\n')
            else:
                # 记录没有结果的数据库
                no_hits_out.write(f"{db_path.split('/')[-1]}\n")
                dbs_with_no_hits += 1

    print(f"BLAST results have been saved to {blast_results_file}.")
    print(f"Total databases processed: {total_dbs}")
    print(f"Databases with no hits: {dbs_with_no_hits}")
    print(f"List of databases with no hits have been saved to {no_hits_file}.")


if __name__ == "__main__":
    query_file = '/Users/lvdongyuan/Downloads/alphaproteobacteria/alphaproteodb/file/ssb.faa'
    db_directory = '/Users/lvdongyuan/Desktop/alphadb'
    blast_results_file = '/Volumes/ex/recETresults/ssb_lddeipf.txt'
    no_hits_file = '/Volumes/ex/testresults/no_hits_lddeipf.txt'  # 文件记录没有结果的数据库
    patternfile = '/Users/lvdongyuan/Downloads/alphaproteobacteria/alphaproteodb/file/pattern_lddeipf.txt'
    max_hits = 1

    # Function Call
    blast_and_extract(db_directory, query_file, blast_results_file, no_hits_file, patternfile, max_hits)
