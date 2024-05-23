def read_ids_from_file(file_path):
    """读取文件中的ID，并存储在集合中"""
    with open(file_path, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    return ids

def compare_ids(file_path1, file_path2, output_union_path):
    """比较两个文件中的ID，并输出它们的并集、交集和差集的数量"""
    # 读取两个文件中的ID
    ids1 = read_ids_from_file(file_path1)
    ids2 = read_ids_from_file(file_path2)

    # 计算并集、交集和差集
    union = ids1.union(ids2)
    intersection = ids1.intersection(ids2)
    difference1 = ids1.difference(ids2)
    difference2 = ids2.difference(ids1)

    # 输出结果
    print(f"Number of IDs in Union: {len(union)}")
    print(f"Number of IDs in Intersection: {len(intersection)}")
    print(f"Number of IDs in Difference (File1 - File2): {len(difference1)}")
    print(f"Number of IDs in Difference (File2 - File1): {len(difference2)}")
    with open(output_union_path, 'w') as output_file:
        for id in union:
            output_file.write(id + '\n')
