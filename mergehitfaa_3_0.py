import os

def merge_faa_files(file_list_path, output_file_path, directory_path):
    # 打开输出文件
    with open(output_file_path, 'w') as output_file:
        # 读取包含文件名的列表文件
        with open(file_list_path, 'r') as file_list:
            for file_name in file_list:
                file_name = file_name.strip() + '.faa'  # 确保文件名格式正确
                file_path = os.path.join(directory_path, file_name)  # 构建完整路径
                if os.path.exists(file_path):
                    with open(file_path, 'r') as file:
                        output_file.write(file.read() + '\n')  # 读取并写入内容
                        print(f"Merged: {file_name}")
                else:
                    print(f"File not found: {file_name}")

if __name__ == "__main__":
    file_list_path = '/Volumes/ex/recETresults/extractdbac.txt'
    output_file_path = '/Volumes/ex/recETresults/lddeipf.faa'
    directory_path = '/Users/lvdongyuan/Downloads/alphaproteobacteria/gzdb'

    # 调用函数合并文件
    merge_faa_files(file_list_path, output_file_path, directory_path)
