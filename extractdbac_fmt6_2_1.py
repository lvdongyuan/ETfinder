import re

def extract_database_ids(input_data_path, output_path):
    # 存储提取到的所有数据库标识
    extracted_ids = []

    # 定义正则表达式，匹配格式为"GCA_xxxxxxx.x"的访问编号
    accession_pattern = re.compile(r'GCA_\d+\.\d+')

    # 打开文件并逐行读取数据
    with open(input_data_path, 'r') as file:
        for line in file:
            # 检查行是否以"Results for database:"开始
            if line.startswith("Results for database:"):
                # 提取冒号之后的内容（除去首尾空白）
                database_id = line.strip().split("Results for database: ")[1]
                # 使用正则表达式查找匹配的Accession编号
                match = accession_pattern.search(database_id)
                if match:
                    extracted_ids.append(match.group() + '\n')

    # 将所有提取到的数据库标识写入输出文件
    with open(output_path, 'w') as output_file:
        output_file.writelines(extracted_ids)

    print(f"Extracted {len(extracted_ids)} database identifiers.")

if __name__ == " __main__":
    input_data_path = '/Volumes/ex/recETresults/ssb_lddeipf.txt'
    output_path = '/Volumes/ex/recETresults/lddeipfdbac.txt'

    # 调用函数以提取并保存数据库标识
    extract_database_ids(input_data_path, output_path)
