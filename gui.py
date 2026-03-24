import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import matplotlib
import matplotlib.patheffects as PathEffects

from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
from matplotlib import rcParams
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from click import progressbar
from tkinter import ttk, messagebox
from attacks.rsa.partial_key_exposure import attack
from sage.rings.integer import Integer

def option_changed(event):
    selected_option = combo_option.get()
    if selected_option == None:
        entry_traversal_bit_length.config(state='disabled')
        entry_traversal_bit_length.delete(0, tk.END)  # 清空输入框
    else:
        entry_traversal_bit_length.config(state='normal')
        entry_traversal_bit_length.delete(0, tk.END)  # 清空输入框
        entry_traversal_bit_length.insert(0, '10')  # 设置默认值为 10

        # 当选择 'LSB' 或 'MSB' 时，显示提示框
        if selected_option in ['MSB', 'LSB']:
            messagebox.showinfo("提示", "正在使用额外遍历攻击,并进行并行化加速")

matplotlib.use("TkAgg")

# 定义嵌入 Tkinter 的绘图函数
def plot_dm_dl(ax, d_bit_length, dm_bit_length, dl_bit_length):
    # 清空之前的图形
    ax.clear()

    # 设置坐标轴的范围
    ax.set_xlim(0, d_bit_length)
    ax.set_ylim(-1, 1)
    ax.tick_params(axis='x', which='both', bottom=True, top=False)  # 启用底部刻度


    # 计算各个区域的起点和终点
    dm_start = 0
    dm_end = dm_bit_length
    dl_start = d_bit_length - dl_bit_length
    dl_end = d_bit_length

    # 绘制 d 区域 (浅灰色背景)
    ax.add_patch(patches.Rectangle((0, -0.5), d_bit_length, 0.5, edgecolor='black', facecolor='lightgray', lw=2, alpha=0.5))
   
    # 绘制 d_m 区域 (红色)及其标签
    if dm_bit_length > 0:
        ax.add_patch(patches.Rectangle((dm_start, -0.5), dm_bit_length, 0.5, edgecolor='black', facecolor='darkred', lw=2))
        ax.text(dm_start + dm_bit_length / 2, -0.7, r'$d_m$', fontsize=12, ha='center')

    # 绘制 d_l 区域 (绿色)及其标签
    if dl_bit_length > 0:
        ax.add_patch(patches.Rectangle((dl_start, -0.5), dl_bit_length, 0.5, edgecolor='black', facecolor='green', lw=2))
        ax.text(dl_start + dl_bit_length / 2, -0.7, r'$d_l$', fontsize=12, ha='center')

    #  d 
    ax.text(d_bit_length / 2, 0.7, "d", fontsize=14, ha='center', color='black')
    
    # 去掉竖坐标轴的刻度和刻度线
    ax.set_yticks([])  # 去掉竖坐标轴的刻度值
    ax.spines['left'].set_visible(False)  # 去掉左边的竖坐标轴的轴线
    ax.spines['right'].set_visible(False)  

    # 更新绘图
    canvas.draw()

def reset_for_new_attack():
    # 清空所有的结果显示区域
    text_p_output.config(state='normal')
    text_p_output.delete("1.0", tk.END)
    text_p_output.config(state='disabled')

    text_q_output.config(state='normal')
    text_q_output.delete("1.0", tk.END)
    text_q_output.config(state='disabled')

    text_d_output.config(state='normal')
    text_d_output.delete("1.0", tk.END)
    text_d_output.config(state='disabled')

    # 重置状态显示
    result_label_status.config(text="")


def run_attack():
    # 调用 `run_attack()` 函数时，先调用 `reset_for_new_attack()` 以清除旧数据
    reset_for_new_attack() 

    # 第一个弹窗：提示最优多项式构造完成
    messagebox.showinfo("提示", "最优多项式构造完成\n\n\n                 ")

    # 第二个弹窗：提示攻击中
    messagebox.showinfo("提示", "攻击中，请等待...\n\n\n                 ")

    # 获取用户输入的值
    N = text_N.get("1.0", tk.END).strip()  # 使用 Text 小部件获取 N
    e = text_e.get("1.0", tk.END).strip()  # 使用 Text 小部件获取 e
    d_bit_length = int(text_d_bit_length.get("1.0", tk.END).strip())
    dl = int(text_dl.get("1.0", tk.END).strip())
    dl_bit_length = int(text_dl_bit_length.get("1.0", tk.END).strip())
    dm = int(text_dm.get("1.0", tk.END).strip())
    dm_bit_length = int(text_dm_bit_length.get("1.0", tk.END).strip())
    
    # 获取 m 和 t 的值，如果没有输入则使用默认值
    m_input = entry_m.get()
    t_input = entry_t.get()
    
    m = int(m_input) if m_input else 4  # 如果没有输入，默认为 4
    t = int(t_input) if t_input else 3  # 如果没有输入，默认为 3
    
    option = combo_option.get()
    # print(option,type(option),"option")
    # 只有当 option 不为 None 时，才需要 traversal_bit_length
    if option != "None":
        traversal_bit_length = int(entry_traversal_bit_length.get())
    else:
        option = None
        traversal_bit_length = None
    print("traversal_bit_length =", traversal_bit_length)
    
    #调用函数生成图形，显示 dm 和 dl 的比例关系
    plot_dm_dl(ax,d_bit_length, dm_bit_length, dl_bit_length)
    
    try:
        nn = int(N).bit_length()
        beta = d_bit_length / nn
        if (beta > 0.25 and beta < 0.293 and dl_bit_length < 1 and dm_bit_length < 1):
            ee = Integer(e)
        else:
            ee = int(e)
	
	# In BDF 4.2, we need factor_e = false
        n = int(N).bit_length()
        t_ = int(e).bit_length() - 1
        if 0 <= t_ <= n / 2 and dm_bit_length >= n - t_:
            factor_e = False
        else:
            factor_e = True
        start_time = time.time()  # 记录开始时间
        p_, q_, d_ = attack(int(N), ee, d_bit_length, dl, dl_bit_length, dm, dm_bit_length, factor_e=factor_e, m=m, t=t, traversal_bit_length=traversal_bit_length, option=option)
        end_time = time.time()  # 记录结束时间
        elapsed_time = end_time - start_time  # 计算消耗的时间
    
        # 输出耗时到状态标签
        result_label_status.config(text=f"攻击执行状态：耗时 {elapsed_time:.2f} 秒")
        # 检查攻击结果
        if p_ is None or q_ is None or d_ is None:
            # 如果结果无效，则抛出异常
            raise ValueError("攻击失败")  

    except Exception as ex:
        print(f"Error during attack: {ex}")
        # 显示攻击失败的提示框
        messagebox.showinfo("攻击失败", "建议您更换更大的 m 和 t ,比如在原来的基础上分别加 1.")
        return  # 结束运行
    
    # 显示攻击结果 p_, q_, d_ 在对应文本框中
    text_p_output.config(state='normal')
    text_p_output.delete("1.0", tk.END)
    text_p_output.insert(tk.END, f"{p_}")
    text_p_output.config(state='disabled')

    text_q_output.config(state='normal')
    text_q_output.delete("1.0", tk.END)
    text_q_output.insert(tk.END, f"{q_}")
    text_q_output.config(state='disabled')

    text_d_output.config(state='normal')
    text_d_output.delete("1.0", tk.END)
    text_d_output.insert(tk.END, f"{d_}")
    text_d_output.config(state='disabled')

    # 更新状态显示，提示攻击执行完毕
    messagebox.showinfo("提示", "已完成攻击\n\n\n                 ")  

# 创建主窗口
root = tk.Tk()
root.title("RSA Key Attack")
root.geometry("1550x1150")  

# 创建 Matplotlib 图形对象，并嵌入到 Tkinter
fig, ax = plt.subplots(figsize=(6,2))  # 定义绘图区域
canvas = FigureCanvasTkAgg(fig, master=root)  # 创建画布，嵌入到 Tkinter
canvas.get_tk_widget().grid(row=7, column=2, columnspan=4, padx=10, pady=10)  # 将绘图嵌入到界面中的右下角


# 处理窗口关闭
def on_closing():
    if tk.messagebox.askokcancel("提示","你确定要关闭窗口嘛？\n\n\n "):
        root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)

# 字体大小
font_size = 18

# 创建标签和输入框
labels = ['N', 'e', 'd_bit_length', 'dl', 'dl_bit_length', 'dm', 'dm_bit_length', 'm', 't']
text_widgets = []

# 增加标签和输入框之间的间距
label_padding_x = 30  # 标签的左右间距
label_padding_y = 10  # 缩小标签和输入框之间的上下间距

# 创建滚动条的辅助函数
def create_scrollable_text_widget(row, column, height=2):
    frame = tk.Frame(root)
    frame.grid(row=row, column=column, padx=label_padding_x, pady=label_padding_y)
    
    text_widget = tk.Text(frame, font=("Arial", font_size), height=height, width=20)
    text_widget.pack(side=tk.LEFT)
    
    scrollbar = tk.Scrollbar(frame, command=text_widget.yview)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    text_widget.config(yscrollcommand=scrollbar.set)
    
    return text_widget

# 比例图，内容为 "私钥部分泄露示意图"
title_label = ttk.Label(root, text="私钥部分泄露示意图", font=("Arial", 20), foreground="black")
title_label.grid(row=6, column=2, padx=label_padding_x, pady=label_padding_y)

# 系统名称
title_label = ttk.Label(root, text="RSA部分私钥泄露攻击系统", font=("Arial", 26), foreground="blue")
title_label.grid(row=0, column=2, padx=label_padding_x, pady=label_padding_y)


# 输入参数
title_label = ttk.Label(root, text="输入参数", font=("Arial", 20), foreground="black")
title_label.grid(row=1, column=1, padx=label_padding_x, pady=label_padding_y)

# N 和 e 使用带滚动条的 Text 小部件
lbl_N = tk.Label(root, text="公共模数N", font=("Arial", font_size))
lbl_N.grid(row=2, column=0, padx=label_padding_x, pady=label_padding_y)
text_N = create_scrollable_text_widget(2, 1)

lbl_e = tk.Label(root, text="公钥e", font=("Arial", font_size))
lbl_e.grid(row=3, column=0, padx=label_padding_x, pady=label_padding_y)
text_e = create_scrollable_text_widget(3, 1)

# d_bit_length 在 e 的正下方
lbl_d_bit_length = tk.Label(root, text="私钥d比特长度", font=("Arial", font_size))
lbl_d_bit_length.grid(row=4, column=0, padx=label_padding_x, pady=label_padding_y)
text_d_bit_length = create_scrollable_text_widget(4, 1)

# 创建 dl 和 dm 及其对应的滚动条 Text 输入框
lbl_dl = tk.Label(root, text="私钥d的LSB泄露信息dl", font=("Arial", font_size))
lbl_dl.grid(row=5, column=0, padx=label_padding_x, pady=label_padding_y)
text_dl = create_scrollable_text_widget(5, 1)

lbl_dl_bit_length = tk.Label(root, text="dl比特长度", font=("Arial", font_size))
lbl_dl_bit_length.grid(row=6, column=0, padx=label_padding_x, pady=label_padding_y)
text_dl_bit_length = create_scrollable_text_widget(6, 1)

lbl_dm = tk.Label(root, text="私钥d的MSB泄露信息dm", font=("Arial", font_size))
lbl_dm.grid(row=7, column=0, padx=label_padding_x, pady=label_padding_y)
text_dm = create_scrollable_text_widget(7, 1)

lbl_dm_bit_length = tk.Label(root, text="dm比特长度", font=("Arial", font_size))
lbl_dm_bit_length.grid(row=8, column=0, padx=label_padding_x, pady=label_padding_y)
text_dm_bit_length = create_scrollable_text_widget(8, 1)

lbl_m = tk.Label(root, text="m(矩阵维度信息)", font=("Arial", font_size))
lbl_m.grid(row=9, column=0, padx=label_padding_x, pady=label_padding_y)
entry_m = tk.Entry(root, font=("Arial", font_size))
entry_m.grid(row=9, column=1, padx=label_padding_x, pady=label_padding_y)
# 设置 m 的默认值为 4
entry_m.insert(0, '4')  

lbl_t = tk.Label(root, text="t(矩阵维度信息)", font=("Arial", font_size))
lbl_t.grid(row=10, column=0, padx=label_padding_x, pady=label_padding_y)
entry_t = tk.Entry(root, font=("Arial", font_size))
entry_t.grid(row=10, column=1, padx=label_padding_x, pady=label_padding_y)
# 设置 t 的默认值为 3
entry_t.insert(0, '3')  

# 选项 option 的下拉框
lbl_option = tk.Label(root, text="额外遍历选项", font=("Arial", font_size))
lbl_option.grid(row=11, column=0, padx=label_padding_x, pady=label_padding_y)
combo_option = ttk.Combobox(root, values=[None, "MSB", "LSB"], font=("Arial", font_size))
combo_option.grid(row=11, column=1, padx=label_padding_x, pady=label_padding_y)
combo_option.current(0)

# 绑定事件，当选择更改时调用 option_changed 函数
combo_option.bind("<<ComboboxSelected>>", option_changed)

# traversal_bit_length 的输入框（默认禁用）
lbl_traversal_bit_length = tk.Label(root, text="额外遍历比特数", font=("Arial", font_size))
lbl_traversal_bit_length.grid(row=12, column=0, padx=label_padding_x, pady=label_padding_y)
entry_traversal_bit_length = tk.Entry(root, font=("Arial", font_size), state='disabled')  # 初始状态禁用
entry_traversal_bit_length.grid(row=12, column=1, padx=label_padding_x, pady=label_padding_y)

# 运行攻击按钮
#btn_run = tk.Button(root, text="运行攻击", command=run_attack, font=("Arial", font_size))
btn_run = tk.Button(root, text="运行攻击", command=run_attack, font=("Arial", font_size), width=20, height=2, padx=10, pady=5)
btn_run.grid(row=10, column=2, columnspan=2, pady=20)

# 显示攻击执行状态
result_label_status = tk.Label(root, text="", font=("Arial", font_size))
result_label_status.grid(row=12, column=2, columnspan=2)

# 右侧布局的结果显示部分
lbl_result_section = tk.Label(root, text="攻击结果", font=("Arial", font_size, "bold"))
lbl_result_section.grid(row=1, column=3, padx=label_padding_x, pady=label_padding_y)

# p_ 结果显示及滚动条文本框
lbl_p = tk.Label(root, text="素因子p", font=("Arial", font_size))
lbl_p.grid(row=2, column=2, padx=label_padding_x, pady=label_padding_y)
text_p_output = create_scrollable_text_widget(2, 3)

# q_ 结果显示及滚动条文本框
lbl_q = tk.Label(root, text="素因子q", font=("Arial", font_size))
lbl_q.grid(row=3, column=2, padx=label_padding_x, pady=label_padding_y)
text_q_output = create_scrollable_text_widget(3, 3)

# d_ 结果显示及滚动条文本框
lbl_d = tk.Label(root, text="私钥d", font=("Arial", font_size))
lbl_d.grid(row=4, column=2, padx=label_padding_x, pady=label_padding_y)

text_d_output = create_scrollable_text_widget(4, 3)

# 启动主事件循环
root.mainloop()
