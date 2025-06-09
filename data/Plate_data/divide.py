

def generate_kirchhoff_plate_input_with_load(nrows, ncols, q, dx, dy, E=2.0e11, nu=0.3, thickness=0.01):
    """
    生成带均布载荷q的Kirchhoff板单元输入文件，自动计算等效节点载荷。

    q: 均布载荷大小，单位力/面积，向下为负数
    其余参数同之前
    """
    n_nodes_x = ncols + 1
    n_nodes_y = nrows + 1
    n_nodes = n_nodes_x * n_nodes_y
    n_elements = nrows * ncols

    lines = []
    lines.append("Kirchhoff Plate Uniform Load Example")
    lines.append(f"{n_nodes} 1 1 1")

    # 节点定义
    # 节点定义
    for j in range(n_nodes_y):
        for i in range(n_nodes_x):
            node_id = j * n_nodes_x + i + 1
            x = i * dx
            y = j * dy

            # 判断是否是边界节点
            is_boundary = (i == 0 or i == n_nodes_x - 1 or j == 0 or j == n_nodes_y - 1)
            bc_flag = "1 1 1" if is_boundary else "0 0 0"

            lines.append(f"{node_id} {'1 1'} {bc_flag} {'1'} {x:.6f} {y:.6f} 0")

    # 计算等效节点载荷
    # 每个节点自由度对应一个编号：
    # 假设1为垂直力，2、3为力矩，载荷格式：节点号 自由度 载荷值
    # 首先准备一个字典，存储每个节点的载荷向量（3自由度）
    nodal_loads = {i: [0.0, 0.0, 0.0] for i in range(1, n_nodes+1)}

    A = dx * dy

    # 均布载荷q等效分配到单元4节点，每个单元按编号顺序分配载荷
    for row in range(nrows):
        for col in range(ncols):
            n1 = row * n_nodes_x + col + 1
            n2 = n1 + 1
            n4 = n1 + n_nodes_x
            n3 = n4 + 1

            # 每个节点等效载荷
            loads = [
                [ q * A / 4,  q * A * dy / 24, - q * A * dx / 24],
                [ q * A / 4, -q * A * dy / 24,  -q * A * dx / 24],
                [ q * A / 4, q * A * dy / 24, q * A * dx / 24],
                [ q * A / 4,  -q * A * dy / 24, q * A * dx / 24],
            ]
            nodes = [n1, n2, n3, n4]

            for node, load_vec in zip(nodes, loads):
                nodal_loads[node][0] += load_vec[0]
                nodal_loads[node][1] += load_vec[1]
                nodal_loads[node][2] += load_vec[2]

    # 写载荷部分
    # 统计非零载荷数
    load_lines = []
    count_loads = 0
    for node_id, load_vec in nodal_loads.items():
        for dof in range(3):
            val = load_vec[dof]
            if abs(val) > 1e-12:
                count_loads += 1
                # 节点号 自由度 载荷值
                # 自由度编号从1开始
                load_lines.append(f"{node_id} {dof+3} {val:.8f}")

    # 载荷头：载荷编号(写1) 载荷数
    lines.append(f"1 {count_loads}")
    lines.extend(load_lines)
    lines.append(f"6 {n_elements} 1")
    # 材料参数
    lines.append(f"1 {E} {nu} {thickness}")

    # 单元定义
    for row in range(nrows):
        for col in range(ncols):
            elem_id = row * ncols + col + 1
            n1 = row * n_nodes_x + col + 1
            n2 = n1 + 1
            n4 = n1 + n_nodes_x
            n3 = n4 + 1
            lines.append(f"{elem_id} {n1} {n2} {n3} {n4} 1")

    return "\n".join(lines)


if __name__ == "__main__":
    nrows = 32
    ncols = 32
    q = -1.0  # 负号表示向下均布载荷，单位N/m^2
    dx=8/ncols
    dy=8/nrows
    content = generate_kirchhoff_plate_input_with_load(nrows, ncols,q,dx,dy)
    with open("plate_example32x32.dat", "w") as f:
        f.write(content)
    print("带均布载荷的输入文件 kirchhoff_plate_with_load.dat 已生成")