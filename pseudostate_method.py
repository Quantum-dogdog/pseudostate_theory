from jplephem.spk import SPK
import numpy as np
import math
from constant import *
from erti import Twobody


np.set_printoptions(precision=5)  #小数点后两位
kernel = SPK.open('de405.bsp')



'''
position, velocity = kernel[0,3].compute_and_differentiate(jd)  #3是地球barycenter 
print('x y z', position)    #单位km

velocity_per_second = velocity / 86400.0
print('vx vy vx',velocity_per_second)  #单位km/s
'''

'''
以下代码已与GMAT交叉验证，本代码中的儒略日对应GMAT中当日零时，坐标对应EarthMJ2000Eq下的月球坐标.
'''

def gregorian_to_julian_day(year, month, day):
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + A // 4
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return JD

# 示例
year = 2020
month = 12
day = 1
jd = gregorian_to_julian_day(year, month, day)
print(f"The Julian Day for {year}-{month}-{day} is {jd}")



print()


position = kernel[3,301].compute(jd)      
position -= kernel[3,399].compute(jd) 
print('moon x y z', position)    #单位km

diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
print('diyue_juli', diyue_juli)
vector_A = position
moon_position_1 = vector_A


print()


position, velocity = kernel[3,301].compute_and_differentiate(jd+5)  #3是地球barycenter
v1 = velocity
position, velocity = kernel[3,399].compute_and_differentiate(jd+5)  #3是地球barycenter
v1 -= velocity
sudu = v1 / 86400.0
print('vx vy vx',sudu)  #单位km/s






position = kernel[3,301].compute(jd+5)   #5d后moon位置    
position -= kernel[3,399].compute(jd+5) 
print('moon x y z', position)    #单位km
diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
print('diyue_juli', diyue_juli)
zuobiao = diyue_juli
vector_B = position
moon_position_2 = vector_B



# 计算点积
dot_product = sum(a*b for a, b in zip(vector_A, vector_B))

# 计算向量的模
norm_A = math.sqrt(sum(a**2 for a in vector_A))
norm_B = math.sqrt(sum(b**2 for b in vector_B))

# 计算夹角的余弦值
cos_theta = dot_product / (norm_A * norm_B)

# 计算夹角（以弧度为单位）
theta_radians = math.acos(cos_theta)

# 将弧度转换为度
theta_degrees = math.degrees(theta_radians)

print(f"夹角的弧度: {theta_radians}")
print(f"夹角的度数: {theta_degrees}")




# 计算单位向量
unit_vector_b = vector_B / np.linalg.norm(vector_B)

r_park = r_earth + 200              #表示200km的LEO轨道上


r_e0 = - r_park * unit_vector_b                  

print('r_e,0',r_e0)
'''
r_e0_gaodu = np.sqrt(r_e0[0]**2 + r_e0[1]**2 + r_e0[2]**2)
print('r_e0_gaodu', r_e0_gaodu)

'''

# 将 vb 归一化作为 x 轴

va = vector_A / np.linalg.norm(vector_A)
vb = vector_B / np.linalg.norm(vector_B)
x_axis= vb


# 计算法向量（z轴）
z_axis = np.cross(va, vb)
z_axis = z_axis / np.linalg.norm(z_axis)  # 确保z轴也是单位向量
# 计算y轴
y_axis = np.cross(z_axis, x_axis)    #x轴向右，y轴向上
y_axis = y_axis / np.linalg.norm(y_axis)  # 确保y轴也是单位向量

# 确保坐标轴正交
assert np.allclose(np.dot(x_axis, y_axis), 0)
assert np.allclose(np.dot(y_axis, z_axis), 0)
assert np.allclose(np.dot(z_axis, x_axis), 0)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# 重新设置坐标轴范围以适应月球位置的实际范围
max_range = np.array([moon_position_1.max(), moon_position_2.max()]).max() * 1.1

# 创建 3D 图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 画出坐标轴
ax.quiver(0, 0, 0, x_axis[0], x_axis[1], x_axis[2], length=max_range, color='r', arrow_length_ratio=0.1, label='X Axis')
ax.quiver(0, 0, 0, y_axis[0], y_axis[1], y_axis[2], length=max_range, color='g', arrow_length_ratio=0.1, label='Y Axis')
ax.quiver(0, 0, 0, z_axis[0], z_axis[1], z_axis[2], length=max_range, color='b', arrow_length_ratio=0.1, label='Z Axis')

# 标出两个 moon xyz 位置
ax.scatter(moon_position_1[0], moon_position_1[1], moon_position_1[2], color='m', label='Moon Position 1')
ax.scatter(moon_position_2[0], moon_position_2[1], moon_position_2[2], color='c', label='Moon Position 2')

# 计算并绘制 re0 向量
re0_length = np.linalg.norm(r_e0)  # 计算 re0 向量的长度
re0_direction = r_e0 / re0_length  # 归一化 re0 向量以获取方向
# 绘制 re0 向量
ax.quiver(0, 0, 0, re0_direction[0], re0_direction[1], re0_direction[2], length=re0_length, color='y', arrow_length_ratio=0.3, label='RE0 Vector')
'''
from mpl_toolkits.mplot3d import art3d

# 在原点绘制一个半径为 6378 的蓝色球体，增加立体感
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6378 * np.outer(np.cos(u), np.sin(v))
y = 6378 * np.outer(np.sin(u), np.sin(v))
z = 6378 * np.outer(np.ones(np.size(u)), np.cos(v))

# 将球体添加到图形中，调整透明度和网格线密度
ax.plot_surface(x, y, z, color='blue', rstride=1, cstride=1, alpha=0.7, linewidth=0)
'''



new_origin = moon_position_2
# 画出新坐标系的坐标轴
ax.quiver(new_origin[0], new_origin[1], new_origin[2], x_axis[0], x_axis[1], x_axis[2], length=max_range, color='m', arrow_length_ratio=0.1, label='New X Axis')
ax.quiver(new_origin[0], new_origin[1], new_origin[2], y_axis[0], y_axis[1], y_axis[2], length=max_range, color='c', arrow_length_ratio=0.1, label='New Y Axis')
ax.quiver(new_origin[0], new_origin[1], new_origin[2], z_axis[0], z_axis[1], z_axis[2], length=max_range, color='y', arrow_length_ratio=0.1, label='New Z Axis')

# 标记新原点
ax.scatter(new_origin[0], new_origin[1], new_origin[2], color='k', s=5, label='New Origin')

# 设置图形属性
ax.set_xlim([-max_range, max_range])
ax.set_ylim([-max_range, max_range])
ax.set_zlim([-max_range, max_range])
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')
ax.set_title('Custom Coordinate System with Moon Positions')
ax.legend()

# 显示图形
#plt.show()

######################################################################
#                                                                    #
#分割线                                                              #
#                                                                    #
######################################################################




def conic2distance(r_i, v_i, rs):
    # 计算 r_i 的模
    ri = np.linalg.norm(r_i)
    # 计算 v_i 的模
    vi = np.linalg.norm(v_i)
    
    # 确保初始位置大于月球半径
    if ri < r_moon:
        raise ValueError("错误：初始位置小于月球半径。")
    
    # 计算逃逸速度 (km/s)
    v_esc = np.sqrt(2 * mu_moon / ri)
    
    # 确保初始速度大于逃逸速度
    if vi < v_esc:
        raise ValueError("错误：初始速度不大于逃逸速度。")
    h = ri * vi
    e =h**2/(mu_moon * ri) -1
    if e < 1:
        raise ValueError("错误：不是双曲线。")
    theta_inf = math.degrees(math.acos(-1/e))    #度
    jiao = find_jiao(h, e, theta_inf, rs)
    right_side = math.sqrt((e - 1) / (e + 1)) * math.tan(math.radians(jiao / 2))

    # 使用双曲反正切函数（atanh）计算 F
    F = 2 * math.atanh(right_side)
    M_h = e * math.sinh(F)-F
    ts = M_h * h**3 /((e**2 - 1)**1.5 * mu_moon**2)
    ri_i = r_i / np.linalg.norm(r_i)
    ri_rotated = rotate_vector(ri_i, jiao) 
    r_s = [x * rs for x in ri_rotated]
    v_chuizhi = h / rs
    v_shuiping = mu_moon * e * math.sin(math.radians(jiao))/h
    v_sp = [x * v_shuiping for x in ri_rotated]
    # 计算垂直向量
    qq = np.array([-ri_i[1], ri_i[0]])

    # 计算垂直向量的长度
    length = np.linalg.norm(qq)

    # 计算单位向量
    unit_vector = qq / length


    v_cz = v_chuizhi * unit_vector
    v_s = [-v_sp[0] - v_cz[0], -v_sp[1] - v_cz[1], 0]  #这一块正负的定义顺时针逆时针的搞迷糊了，但是这样写才能通过GMAT交叉验证
    r_s = convert_to_3d_vector(r_s)
    #print(ri_rotated)
    print('jiao',jiao)  #这个倒是可以解释成逆时针方向为负
    ts = abs(ts)
    return ts, r_s, v_s
def find_jiao(h, e, theta_inf, rs):
    jiao_found = None

    # 遍历角度jiao从1到100
    for jiao in np.arange(-abs(theta_inf), 0, step=0.1):
        # 计算r_loop
        r_loop = h**2 / (mu_moon * (1 + e * math.cos(math.radians(jiao))))
        
        # 检查条件是否满足
        if abs(r_loop) - rs < 0.1:
            jiao_found = jiao
            break  # 找到满足条件的jiao后退出循环

    # 循环结束后检查是否找到了满足条件的jiao
    if jiao_found is not None:
        return jiao_found
    else:
        # 如果没有找到满足条件的jiao，可以返回一个提示或者None
        return None
 
def rotate_vector(ri, theta):
    # 将角度转换为弧度
    theta_rad = math.radians(theta)
    
    # 旋转矩阵
    rotation_matrix = [
        [math.cos(theta_rad), math.sin(theta_rad)],
        [-math.sin(theta_rad), math.cos(theta_rad)]
    ]
    
    # 向量
    xi = ri[0]
    yi = ri[1]
    # 应用旋转矩阵
    x_rotated = rotation_matrix[0][0] * xi + rotation_matrix[0][1] * yi
    y_rotated = rotation_matrix[1][0] * xi + rotation_matrix[1][1] * yi
    
    # 返回旋转后的向量
    return [x_rotated, y_rotated]

def convert_to_3d_vector(r_s):
    # 假设 r_s 是一个二维向量 [x, y]
    x = r_s[0]
    y = r_s[1]
    # 返回一个三维向量，其中第三维为0
    return [x, y, 0]





r_park = r_moon + 300              #表示300km的轨道上
#r_i = np.array([r_park, 0, 0])
r_i = np.array([1018, 1764, 0])
print('r_i',r_i)
#print(y_axis)
#v_i = np.array([0, 3.2, 0])
v_i = np.array([2.03, 2.3, 0])
print('v_i',v_i)
rs = 24 * r_earth
ts, r_s, v_s = conic2distance(r_i, v_i, rs)
print(f"{ts}\n{r_s}\n{v_s}")
#算出来的结果
#60974.07831904006
#[np.float64(-44217.92550631152), np.float64(146549.78984306968), 0]
#[np.float64(-0.676993315848212), np.float64(2.2011419058573565), 0]
#GMAT的结果
#64082.71979405545 
#X  =  -44217.925506501 km Y  =   146066.98191424 km
#VX =  -0.7197437663753 km/sec VY =   2.2301163742187 km/sec
r_s_array = np.array(r_s)
v_s_array = np.array(v_s)                     
v_istar = v_s_array
r_istar = r_s_array + ts * v_s_array
R_m  = np.array([zuobiao, 0, 0])

# 在新坐标系下表示速度向量
sudu_new_coords = np.array([
    np.dot(sudu, x_axis),
    np.dot(sudu, y_axis),
    np.dot(sudu, z_axis)
])


V_m = sudu_new_coords
print('V_m',V_m)
R_i2star = R_m + r_istar
V_i2star = V_m + v_istar
print(f"{R_i2star}\n{V_i2star}")
print()
tb = Twobody()
h, i, e, raan, aop, ta = tb.statetoelement(R_i2star, V_i2star, mu_earth)
print(f"h{h}\ni{i}\ne{e}\nraan{raan*57.3}\naop{aop*57.3}\nta{ta*57.3}")

#以上都是对的，以下应该也是对的，只是只适合椭圆


R_EI, V_EI = tb.elementtostate(h, i, e, raan, aop, 0)
print(f"{R_EI}\n{V_EI}")
R = np.linalg.norm(R_EI)
V = np.linalg.norm(V_EI)
print()
print(f"{R}\n{V}")
