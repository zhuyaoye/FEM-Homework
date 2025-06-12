import numpy as np
import math
import matplotlib.pyplot as plt
import os
import re
def NmatPlate(eta, psi, nen, ae, be):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	# parent coordinates at nodes
	eta_I = np.array([-1, 1, 1, -1])
	psi_I = np.array([-1, -1, 1, 1])
	
	N = np.zeros((1, 12))
	
	for i in range(nen):
		N[0,3*i] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(2 + eta_I[i]*eta + psi_I[i]*psi - eta**2 - psi**2)
				
		N[0,3*i+1] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(-be * psi_I[i] * (1 - psi**2))
				
		N[0,3*i+2] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(ae * eta_I[i] * (1 - eta**2))

	return N


def BmatPlate(eta, psi, nen, ae,be):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt by explicit expression

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	# parent coordinates at nodes
	eta_val = np.array([-1, 1, 1, -1])
	psi_val = np.array([-1, -1, 1, 1])
	
	#Calculate the B_M matrix
	B = np.zeros((3, 12))
	for i in range(nen):
		B[:, 3*i:3*i+3] = 1.0 / (4 * ae * be) * np.array([[ \
								-3*be/ae*eta_val[i]*eta*(1+psi_val[i]*psi), \
								0, \
								-be*eta_val[i]*(1+3*eta_val[i]*eta)*(1+psi_val[i]*psi)], \
								[-3*ae/be*psi_val[i]*psi*(1+eta_val[i]*eta), \
								ae*psi_val[i]*(1+3*psi_val[i]*psi)*(1+eta_val[i]*eta), \
								0], \
								[eta_val[i]*psi_val[i]*(4-3*eta**2-3*psi**2), \
								be*eta_val[i]*(3*psi**2+2*psi_val[i]*psi-1), \
								ae*psi_val[i]*(1-2*eta_val[i]*eta-3*eta**2)]])

	# Compute Jacobian determination
	detJ = ae * be

	return B, detJ


def compute_exact_deflection(x, y, ae, be, D=18315.01831501832, nu=0.3):
    """
    返回给定 (x, y) 点处的精确挠度值 w 及其二阶导数 w_xx, w_yy
    
    参数:
        x, y : 可为标量或 NumPy 数组
        ae, be : 半边长（单元边长的一半）
        D : 板刚度
        nu : 泊松比
    返回:
        w, w_xx, w_yy
    """
    q = -1
    a = ae * 2
    b = be * 2
    pi = math.pi
    K = -4 * q * a**2 / pi**3
    m_all = np.array([1, 3, 5, 7])
    E = np.zeros(8, float)
    E[m_all[0]] = 0.3722 * K
    E[m_all[1]] = -0.0380 * K
    E[m_all[2]] = -0.0178 * K
    E[m_all[3]] = -0.0085 * K

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    w = np.zeros_like(x)
    w_xx = np.zeros_like(x)
    w_yy = np.zeros_like(x)

    for m in m_all:
        a_m = m * pi * b / (2 * a)
        b_m = m * pi * a / (2 * b)

        A1 = 4 * q * a**4 / (pi**5 * D) * (-1)**((m-1)//2) / m**5
        B1 = (a_m * np.tanh(a_m) + 2) / (2 * np.cosh(a_m))
        C1 = 1 / (2 * np.cosh(a_m))
        D1 = m * pi / a

        A2 = -a**2 / (2 * pi**2 * D) * E[m] * (-1)**((m-1)//2) / (m**2 * np.cosh(a_m))
        B2 = a_m * np.tanh(a_m)
        D2 = m * pi / a

        A3 = -b**2 / (2 * pi**2 * D) * E[m] * (-1)**((m-1)//2) / (m**2 * np.cosh(b_m))
        B3 = b_m * np.tanh(b_m)
        D3 = m * pi / b

        cosD1x = np.cos(D1 * x)
        cosD2x = np.cos(D2 * x)
        cosD3y = np.cos(D3 * y)
        sinhD1y = np.sinh(D1 * y)
        coshD1y = np.cosh(D1 * y)
        sinhD2y = np.sinh(D2 * y)
        coshD2y = np.cosh(D2 * y)
        sinhD3x = np.sinh(D3 * x)
        coshD3x = np.cosh(D3 * x)

        # 挠度
        w1 = A1 * cosD1x * (1 - B1 * coshD1y + C1 * D1 * y * sinhD1y)
        w2 = A2 * cosD2x * (D2 * y * sinhD2y - B2 * coshD2y)
        w3 = A3 * cosD3y * (D3 * x * sinhD3x - B3 * coshD3x)
        w += w1 + w2 + w3

        # 曲率
        w1_xx = A1 * (-D1**2 * cosD1x) * (1 - B1 * coshD1y + C1 * D1 * y * sinhD1y)
        w2_xx = A2 * (-D2**2 * cosD2x) * (D2 * y * sinhD2y - B2 * coshD2y)
        w3_xx = A3 * cosD3y * ((2 - B3) * D3**2 * coshD3x + D3**3 * x * sinhD3x)
        w_xx += w1_xx + w2_xx + w3_xx

        w1_yy = A1 * cosD1x * ((2 * C1 - B1) * D1**2 * coshD1y + C1 * D1**3 * y * sinhD1y)
        w2_yy = A2 * cosD2x * ((2 - B2) * D2**2 * coshD2y + D2**3 * y * sinhD2y)
        w3_yy = A3 * (-D3**2 * cosD3y) * (D3 * x * sinhD3x - B3 * coshD3x)
        w_yy += w1_yy + w2_yy + w3_yy

    return w, w_xx, w_yy

def compute_global_L2_error(mesh_nodes, elements, displacements, ae, be):
    """
    mesh_nodes: (N_nodes, 2) 节点坐标
    elements: (N_elem, 4) 每个单元4个节点编号
    displacements: (N_nodes,) 挠度 w 值
    """
    L2_error = 0.0
    
    gauss_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]
    gauss_weights = [1.0, 1.0]
    
    nen = 4
    
    for elem in elements:
        # 节点坐标 (4, 2)
        coords = mesh_nodes[elem, :]
        # 挠度 (4,)
        w_e = displacements[elem]
        
        # 几何参数（正方形元素）
        x_min = np.min(coords[:, 0])
        x_max = np.max(coords[:, 0])
        y_min = np.min(coords[:, 1])
        y_max = np.max(coords[:, 1])
        
        ae_e = (x_max - x_min)/2
        be_e = (y_max - y_min)/2
        
        for i, eta in enumerate(gauss_pts):
            for j, psi in enumerate(gauss_pts):
                wgt = gauss_weights[i] * gauss_weights[j]
                
                # 形函数矩阵
                N = NmatPlate(eta, psi, nen, ae_e, be_e)  # shape (1, 12)
                Nw = N[0, ::3]  # 提取 w 对应的列
                
                w_h = np.dot(Nw, w_e)  # 有限元解
                
                # 对应物理坐标
                x_eta = np.array([-1, 1, 1, -1])
                y_psi = np.array([-1, -1, 1, 1])
                x_gp = np.dot(0.25 * (1 + x_eta * eta) * (1 + y_psi * psi), coords[:, 0])-4
                y_gp = np.dot(0.25 * (1 + x_eta * eta) * (1 + y_psi * psi), coords[:, 1])-4
                
                w_ex = compute_exact_deflection(x_gp, y_gp, 4, 4)[0]
                #print(w_h)
                #print(compute_exact_deflection(0.0,0.0, 4, 4)[0])
                err_sq = (w_h - w_ex)**2
                #print(err_sq)
                
                detJ = ae_e * be_e  # 正方形单元
               
                L2_error += err_sq * detJ * wgt
    
    return np.sqrt(L2_error)

def read_node_file(filepath):
    """
    读取节点位移文件，返回节点坐标数组和对应挠度（Z-DISPLACEMENT）数组。
    假设：
    - 节点编号顺序对应规则网格节点顺序
    - 坐标在mesh生成时自行定义
    - 这里只读Z-DISPLACEMENT作为挠度w
    """
    node_ids = []
    w_vals = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # 正则匹配浮点数
    float_pattern = re.compile(r"[-+]?\d*\.\d+e[-+]?\d+|[-+]?\d+\.?\d*")
    for line in lines:
        if line.strip() == '' or line.lower().startswith('node'):
            continue
        parts = float_pattern.findall(line)
        if len(parts) >= 4:
            node_id = int(parts[0])
            # X, Y, Z displacements
            # 只提取Z方向挠度（第4列）
            z_disp = float(parts[1])
            node_ids.append(node_id)
            w_vals.append(z_disp)
    
    return np.array(node_ids), np.array(w_vals)

def generate_mesh(nx, ny, length_x=8.0, length_y=8.0):
    """
    生成规则矩形网格节点坐标和单元拓扑
    节点编号从左下角开始，先从左到右，再逐行向上
    nx, ny: x 和 y 方向的节点数
    返回：
        nodes: (N_nodes, 2) 节点坐标
        elements: (N_elem, 4) 每个四边形单元的 4 个节点编号（从 0 开始）
    """
    nx+=1
    ny+=1
    x = np.linspace(0, length_x, nx)
    y = np.linspace(0, length_y, ny)
    xv, yv = np.meshgrid(x, y)
    # 先从左到右，再从下到上（按行优先）
    nodes = np.vstack([xv.ravel(), yv.ravel()]).T

    elements = []
    for j in range(ny - 1):  # 最后一行不能再组成单元
        for i in range(nx - 1):  # 最后一列不能再组成单元
            n1 = j * nx + i
            n2 = n1 + 1
            n3 = n1 + nx + 1
            n4 = n1 + nx
            elements.append([n1, n2, n3, n4])
    
    return nodes, np.array(elements)
def main():
    # 3个划分情况
    grids = {
        #'2x2': '/Users/zhuyaoye/Desktop/2x2.txt',
        '4x4': '4x4.txt',
        '8x8': '8x8.txt',
        '16x16': '16x16.txt',
        '32x32': '32x32.txt'
        
    }

    length_x = 8.0
    length_y = 8.0

    L2_errors = []
    h_values = []

    for key, filepath in grids.items():
        # 解析网格节点数
        nx, ny = map(int, key.split('x'))
        
        # 读取节点数据
        node_ids, w_disp = read_node_file(filepath)
        
        # 生成网格坐标与单元
        nodes, elements = generate_mesh(nx, ny, length_x, length_y)
        
        # 挠度w的排序要与节点坐标对应，这里假设节点编号顺序一致且从1开始
        # 文件中节点编号从1开始，代码中从0开始，直接按顺序对应即可
        w_disp_sorted = w_disp[node_ids - 1]
        
        # 计算L2误差
        ae = length_x / nx / 2  # 单元半边长
        be = length_y / ny  / 2
        
        L2_error = compute_global_L2_error(nodes, elements, w_disp_sorted, ae, be)
        L2_errors.append(L2_error)

        # 单元尺寸 h 取单元边长，横纵一致
        h = length_x / nx 
        h_values.append(h)

        #print(f"{key} mesh: h = {h:.4f}, L2 error = {L2_error:.6e}")
    # 转换为 log 空间
    log_h = np.log10(h_values)
    log_error = np.log10(L2_errors)

    # 线性拟合
    slope, intercept = np.polyfit(log_h, log_error, 1)

    print(slope)

    # 绘制双对数收敛曲线
    plt.figure()
    plt.loglog(h_values, L2_errors, '-o')
    plt.xlabel('Element size h')
    plt.ylabel('L2 norm error')
    plt.title('L2 norm error vs element size (log-log)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

if __name__ == "__main__":
    main()