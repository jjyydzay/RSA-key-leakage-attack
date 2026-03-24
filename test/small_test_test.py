import subprocess

# 定义 Sage 文件的路径
sage_file = 'small_exponent_test.sage'

# 调用 SageMath 执行该文件
result = subprocess.run(['sage', sage_file], capture_output=True, text=True)

# 打印输出
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)