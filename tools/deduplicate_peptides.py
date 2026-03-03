import argparse
import os
import sys
import difflib

# 增加递归深度限制，以防数据量极大时递归爆栈
sys.setrecursionlimit(10000)

def get_args():
    parser = argparse.ArgumentParser(
        description="肿瘤新抗原去冗余工具 (聚类版：相似序列全员待定)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='输入文件路径 (.faa)')
    parser.add_argument('-o', '--output_prefix', help='输出前缀')
    parser.add_argument('-t', '--threshold', type=float, default=0.85, 
                        help='相似度阈值 (0.0-1.0), 默认 0.85')
    parser.add_argument('--max-x', type=float, default=0.5,
                        help='序列中允许的最大 X 比例 (0.0-1.0)，超过此比例直接丢弃。默认 0.5')
    return parser.parse_args()

def read_peptide_blocks(file_path):
    """读取FASTA Block"""
    with open(file_path, 'r') as f:
        block = []
        for line in f:
            block.append(line.strip())
            if len(block) == 4:
                yield block
                block = []

def is_similar(seq1, seq2, threshold):
    """判断相似性：包含关系 或 序列相似度"""
    # 包含关系 (Substring)
    if seq1 in seq2 or seq2 in seq1:
        return True, "Substring"
    
    # 序列相似度
    ratio = difflib.SequenceMatcher(None, seq1, seq2).ratio()
    if ratio >= threshold:
        return True, f"Score: {ratio:.2f}"
    
    return False, None

def main():
    args = get_args()
    input_file = args.input
    
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        prefix = os.path.splitext(os.path.basename(input_file))[0]

    file_clean = f"{prefix}_clean.faa"
    file_pending = f"{prefix}_pending.faa"
    file_dropped = f"{prefix}_dropped.faa" # 存放完全重复的
    file_junk = f"{prefix}_junk.faa"       # 存放X太多的

    # 1. 读取并进行“完全去重” (Exact Deduplication)
    # 使用字典: key=VarSeq, value=Block
    # 如果遇到完全相同的序列，直接视为 Dropped (保留第一份)
    unique_map = {}
    dropped_count = 0
    junk_count = 0
    
    print(f"[*] 第一步：读取文件并进行完全去重/质量过滤...")
    
    # 用于写入 Dropped 和 Junk 的临时操作
    with open(file_dropped, 'w') as f_drop, open(file_junk, 'w') as f_junk:
        for block in read_peptide_blocks(input_file):
            var_seq = block[3]
            
            # A. 检查垃圾序列 (X 含量)
            x_ratio = var_seq.count('X') / len(var_seq) if len(var_seq) > 0 else 0
            if x_ratio > args.max_x:
                f_junk.write('\n'.join(block) + '\n')
                junk_count += 1
                continue

            # B. 检查完全重复
            # 为了区分不同基因产生的相同肽段（极为罕见但理论存在），
            # 这里的Key可以只用序列，或者序列+基因名。
            # 既然是去冗余，通常认为序列一样就是重复，保留第一个Header即可。
            if var_seq in unique_map:
                f_drop.write('\n'.join(block) + '\n')
                dropped_count += 1
            else:
                unique_map[var_seq] = block

    unique_records = list(unique_map.values())
    n_unique = len(unique_records)
    print(f"    - 输入总 Block 数（去重后有效）: {n_unique}")
    print(f"    - 丢弃完全重复: {dropped_count}")
    print(f"    - 丢弃低质量(X): {junk_count}")

    # 2. 构建相似性图 (Adjacency List)
    # 复杂度 O(N^2)，对于几千条肽段是可以接受的
    print(f"[*] 第二步：构建相似性网络 (N={n_unique})...")
    adj = {i: [] for i in range(n_unique)}
    
    for i in range(n_unique):
        for j in range(i + 1, n_unique):
            seq_i = unique_records[i][3]
            seq_j = unique_records[j][3]
            
            sim, reason = is_similar(seq_i, seq_j, args.threshold)
            if sim:
                adj[i].append(j)
                adj[j].append(i) # 无向图

    # 3. 提取连通分量 (Connected Components)
    # 每一个分量就是一个“簇”
    visited = set()
    clusters = [] # List of lists of indices

    for i in range(n_unique):
        if i not in visited:
            # 开始 BFS/DFS 寻找连通分量
            component = []
            stack = [i]
            visited.add(i)
            while stack:
                node = stack.pop()
                component.append(node)
                for neighbor in adj[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            clusters.append(component)

    # 4. 写入文件
    print(f"[*] 第三步：写入结果...")
    clean_count = 0
    pending_group_count = 0
    pending_seq_count = 0

    with open(file_clean, 'w') as f_clean, open(file_pending, 'w') as f_pend:
        for cluster in clusters:
            # 如果簇里只有一个元素 -> Clean
            if len(cluster) == 1:
                idx = cluster[0]
                block = unique_records[idx]
                f_clean.write('\n'.join(block) + '\n')
                clean_count += 1
            else:
                # 如果簇里有多个元素 -> 全部 Pending
                pending_group_count += 1
                f_pend.write(f"# --- Pending Group {pending_group_count} (Size: {len(cluster)}) ---\n")
                
                for idx in cluster:
                    block = unique_records[idx]
                    # 为了方便你看，把Header稍微改一下，加上pending标记
                    # 但不要改乱了，就在Var header后面加个注释
                    # 也可以选择不改 Header，靠分隔符区分
                    f_pend.write('\n'.join(block) + '\n')
                    pending_seq_count += 1
                
                f_pend.write("# ------------------------------------------\n")

    print("-" * 50)
    print("处理完成！统计如下：")
    print(f"[Clean]   无任何相似对象 (进入 {file_clean}): {clean_count} 条")
    print(f"[Pending] 存在相似序列 (进入 {file_pending}): {pending_seq_count} 条 (共 {pending_group_count} 组)")
    print(f"[Dropped] 完全重复序列 (进入 {file_dropped}): {dropped_count} 条")
    print(f"[Junk]    低质量X序列  (进入 {file_junk})   : {junk_count} 条")
    print("-" * 50)

if __name__ == "__main__":
    main()
