import numpy as np
import fractions

def find_variables():
    for _ in range(10000):  # 尝试10000次
        # 生成随机有理数 delta, beta, kappa 在 (0, 1) 的开区间
        delta = np.random.uniform(0, 1)
        beta = np.random.uniform(0, 1)
        kappa = np.random.uniform(0, 1)

        # 确保 beta > delta + kappa
        if beta <= delta + kappa:
            continue

        # 计算第一个 if 的条件
        first_condition = delta < (5 / 6 - (1 / 3) * np.sqrt(1 + 6 * beta))

        # 计算第二个 if 的条件
        second_condition = (
            (delta <= (3 - 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa) and 
             beta <= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa)) or
            (delta <= (1 / 3 + (1 / 3) * beta - (1 / 3) * np.sqrt(4 * beta ** 2 + 2 * beta - 2)) and 
             beta >= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa))
        )

        # 检查条件
        if not first_condition and second_condition:
            return delta, beta, kappa

    return None  # 如果没有找到符合条件的变量

def validate_variables(delta, beta, kappa):
    # 确保输入的变量在 (0, 1) 的开区间
    if not (0 < delta < 1) or not (0 < beta < 1) or not (0 < kappa < 1):
        return "变量必须在 (0, 1) 的开区间内"

    # 检查 beta 是否大于 delta + kappa
    if not (beta > delta + kappa):
        return "不符合要求：beta 应该大于 delta + kappa"

    # 计算第一个 if 的条件
    first_condition = delta < (5 / 6 - (1 / 3) * np.sqrt(1 + 6 * beta))

    # 计算第二个 if 的条件
    second_condition = (
        (delta <= (3 - 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa) and 
         beta <= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa)) or
        (delta <= (1 / 3 + (1 / 3) * beta - (1 / 3) * np.sqrt(4 * beta ** 2 + 2 * beta - 2)) and 
         beta >= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa))
    )

    # 检查条件
    if not first_condition and second_condition:
        return "验证符合要求：变量满足所有条件"
    else:
        return "验证不符合要求：变量不满足第二个 if 条件"



# result = find_variables()

# if result:
#     delta, beta, kappa = result
#     print(f"找到的变量：delta = {delta}, beta = {beta}, kappa = {kappa}")
#     print(validate_variables(delta, beta, kappa))
# else:
#     print("没有找到满足条件的变量")

