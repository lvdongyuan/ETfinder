import os

def read_ids_from_file(file_path):
    """从文件中读取ID，并返回一个列表"""
    with open(file_path, 'r') as file:
        ids = [line.strip() for line in file if line.strip()]
    return ids

def extract_number(id):
    """从ID中提取数字部分"""
    num_part = id.split('.')[0]  # 先按'.'分割，取前半部分
    numeric_part = ''.join([char for char in num_part if char.isdigit()])  # 提取数字部分
    return int(numeric_part)  # 转换为整数

def extract_prefix(id):
    """从ID中提取字母前缀部分"""
    num_part = id.split('.')[0]  # 先按'.'分割，取前半部分
    prefix_part = ''.join([char for char in num_part if not char.isdigit()])  # 提取非数字部分
    return prefix_part

def process_ids(file_path1, file_path2, output_dir):
    """从两个文件中读取ID，合并、排序并找出相邻的ID对，同时单独保存每种类型的ID"""
    # 从两个文件中读取ID
    rect_ids = read_ids_from_file(file_path1)
    rece_ids = read_ids_from_file(file_path2)
    all_ids = [(id, 'RecT') for id in rect_ids] + [(id, 'RecE') for id in rece_ids]  # 合并列表并附带隐含标签

    # 按照字母前缀和数字顺序排序ID
    all_ids.sort(key=lambda x: (extract_prefix(x[0]), extract_number(x[0])))

    # 初始化集合用于存储成对的ID
    paired_recT = set()
    paired_recE = set()

    # 打开文件用于写入结果
    consecutive_output_path = os.path.join(output_dir, 'consecutive_id_2.txt')
    with open(consecutive_output_path, 'w') as output:
        previous_id = None
        previous_number = None
        previous_prefix = None
        previous_tag = None  # Initialize previous_tag
        for current_id, current_tag in all_ids:
            current_number = extract_number(current_id)
            current_prefix = extract_prefix(current_id)
            if previous_id is not None and previous_prefix == current_prefix and abs(current_number - previous_number) <= 2 and previous_tag != current_tag:
                output.write(f"{previous_id[0]}\n{current_id}\n---------\n")
                if previous_tag == 'RecT' and current_tag == 'RecE':
                    paired_recT.add(previous_id[0])
                    paired_recE.add(current_id)
                elif previous_tag == 'RecE' and current_tag == 'RecT':
                    paired_recE.add(previous_id[0])
                    paired_recT.add(current_id)
            previous_id = (current_id, current_tag)
            previous_number = current_number
            previous_prefix = current_prefix
            previous_tag = current_tag  # Update previous_tag

    # 将成对的ID写入到指定的文件中
    paired_recT_output_path = os.path.join(output_dir, 'paired_recT_ids.txt')
    with open(paired_recT_output_path, 'w') as recT_output:
        for id in paired_recT:
            recT_output.write(f"{id}\n")

    paired_recE_output_path = os.path.join(output_dir, 'paired_recE_ids.txt')
    with open(paired_recE_output_path, 'w') as recE_output:
        for id in paired_recE:
            recE_output.write(f"{id}\n")
