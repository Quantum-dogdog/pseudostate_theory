from jplephem.spk import SPK
import numpy as np
import math
from constant import *



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


position = kernel[3,301].compute(jd+5)   #5d后moon位置    
position -= kernel[3,399].compute(jd+5) 
print('moon x y z', position)    #单位km
diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
print('diyue_juli', diyue_juli)

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
plt.show()

