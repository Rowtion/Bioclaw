#!/usr/bin/env python3
"""
Bioclaw GUI 管理器
简单图形界面，适合小白用户
"""

import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import subprocess
import os
import threading
import webbrowser

class BioclawGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Bioclaw 管理器")
        self.root.geometry("600x500")
        self.root.resizable(False, False)
        
        # 设置样式
        self.style = ttk.Style()
        self.style.configure('Title.TLabel', font=('Helvetica', 20, 'bold'))
        self.style.configure('Subtitle.TLabel', font=('Helvetica', 12))
        self.style.configure('Action.TButton', font=('Helvetica', 12), padding=10)
        
        self.create_widgets()
        self.check_status()
    
    def create_widgets(self):
        # 标题
        title_frame = ttk.Frame(self.root, padding="20")
        title_frame.pack(fill='x')
        
        ttk.Label(title_frame, text="🧬 Bioclaw", style='Title.TLabel').pack()
        ttk.Label(title_frame, text="生物科研环境管理器", style='Subtitle.TLabel').pack()
        
        # 状态显示
        status_frame = ttk.LabelFrame(self.root, text="运行状态", padding="10")
        status_frame.pack(fill='x', padx=20, pady=10)
        
        self.status_var = tk.StringVar(value="正在检查...")
        ttk.Label(status_frame, textvariable=self.status_var, font=('Helvetica', 14)).pack()
        
        # 操作按钮
        button_frame = ttk.Frame(self.root, padding="20")
        button_frame.pack(fill='x')
        
        # 启动按钮
        self.start_btn = ttk.Button(
            button_frame, 
            text="🚀 启动", 
            command=self.start_bioclaw,
            style='Action.TButton'
        )
        self.start_btn.pack(fill='x', pady=5)
        
        # 停止按钮
        self.stop_btn = ttk.Button(
            button_frame, 
            text="🛑 停止", 
            command=self.stop_bioclaw,
            style='Action.TButton'
        )
        self.stop_btn.pack(fill='x', pady=5)
        
        # 打开 RStudio
        self.rstudio_btn = ttk.Button(
            button_frame, 
            text="📊 打开 RStudio", 
            command=lambda: self.open_browser("http://localhost:8787")
        )
        self.rstudio_btn.pack(fill='x', pady=5)
        
        # 打开 Jupyter
        self.jupyter_btn = ttk.Button(
            button_frame, 
            text="📝 打开 JupyterLab", 
            command=lambda: self.open_browser("http://localhost:8888")
        )
        self.jupyter_btn.pack(fill='x', pady=5)
        
        # 日志显示
        log_frame = ttk.LabelFrame(self.root, text="日志", padding="10")
        log_frame.pack(fill='both', expand=True, padx=20, pady=10)
        
        self.log_text = scrolledtext.ScrolledText(log_frame, height=8, wrap=tk.WORD)
        self.log_text.pack(fill='both', expand=True)
        self.log_text.insert('end', "欢迎使用 Bioclaw！\n")
        self.log_text.insert('end', "密码: bioclaw\n")
        self.log_text.config(state='disabled')
        
        # 底部按钮
        bottom_frame = ttk.Frame(self.root, padding="10")
        bottom_frame.pack(fill='x', side='bottom')
        
        ttk.Button(bottom_frame, text="❓ 帮助", command=self.show_help).pack(side='left', padx=5)
        ttk.Button(bottom_frame, text="🔄 检查更新", command=self.check_update).pack(side='left', padx=5)
        ttk.Button(bottom_frame, text="❌ 退出", command=self.root.quit).pack(side='right', padx=5)
    
    def check_status(self):
        """检查 Bioclaw 运行状态"""
        def check():
            try:
                result = subprocess.run(
                    ["docker-compose", "ps"],
                    cwd=os.path.expanduser("~/.bioclaw"),
                    capture_output=True,
                    text=True
                )
                if "Up" in result.stdout:
                    self.status_var.set("✅ 运行中")
                    self.log("Bioclaw 正在运行")
                else:
                    self.status_var.set("⏹️  已停止")
                    self.log("Bioclaw 未运行")
            except:
                self.status_var.set("❓ 未安装")
                self.log("Bioclaw 可能未安装")
        
        threading.Thread(target=check, daemon=True).start()
    
    def start_bioclaw(self):
        """启动 Bioclaw"""
        def start():
            self.log("正在启动 Bioclaw...")
            self.start_btn.config(state='disabled')
            
            try:
                result = subprocess.run(
                    ["docker-compose", "up", "-d"],
                    cwd=os.path.expanduser("~/.bioclaw"),
                    capture_output=True,
                    text=True
                )
                
                if result.returncode == 0:
                    self.log("✅ 启动成功！")
                    self.log("访问地址:")
                    self.log("  RStudio: http://localhost:8787")
                    self.log("  JupyterLab: http://localhost:8888")
                    self.status_var.set("✅ 运行中")
                    messagebox.showinfo("成功", "Bioclaw 已启动！\n\nRStudio: http://localhost:8787\nJupyterLab: http://localhost:8888")
                else:
                    self.log(f"❌ 启动失败: {result.stderr}")
                    messagebox.showerror("错误", f"启动失败:\n{result.stderr}")
            except Exception as e:
                self.log(f"❌ 错误: {str(e)}")
                messagebox.showerror("错误", str(e))
            finally:
                self.start_btn.config(state='normal')
        
        threading.Thread(target=start, daemon=True).start()
    
    def stop_bioclaw(self):
        """停止 Bioclaw"""
        def stop():
            self.log("正在停止 Bioclaw...")
            self.stop_btn.config(state='disabled')
            
            try:
                result = subprocess.run(
                    ["docker-compose", "down"],
                    cwd=os.path.expanduser("~/.bioclaw"),
                    capture_output=True,
                    text=True
                )
                
                if result.returncode == 0:
                    self.log("✅ 已停止")
                    self.status_var.set("⏹️  已停止")
                    messagebox.showinfo("成功", "Bioclaw 已停止")
                else:
                    self.log(f"❌ 停止失败: {result.stderr}")
            except Exception as e:
                self.log(f"❌ 错误: {str(e)}")
            finally:
                self.stop_btn.config(state='normal')
        
        threading.Thread(target=stop, daemon=True).start()
    
    def open_browser(self, url):
        """打开浏览器"""
        webbrowser.open(url)
        self.log(f"已打开浏览器: {url}")
    
    def show_help(self):
        """显示帮助"""
        help_text = """Bioclaw 使用帮助

1. 首次使用：
   - 点击"启动"按钮
   - 等待 30 秒左右
   - 点击"打开 RStudio"或"打开 JupyterLab"

2. 默认密码：
   - 用户名：任意
   - 密码：bioclaw

3. 数据保存位置：
   - 输入数据：~/.bioclaw/data/
   - 分析结果：~/.bioclaw/outputs/

4. 常见问题：
   - 查看 FAQ.md 文件
   - 访问 GitHub Issues
"""
        messagebox.showinfo("帮助", help_text)
    
    def check_update(self):
        """检查更新"""
        def update():
            self.log("正在检查更新...")
            try:
                result = subprocess.run(
                    ["git", "pull"],
                    cwd=os.path.expanduser("~/.bioclaw"),
                    capture_output=True,
                    text=True
                )
                
                if "Already up to date" in result.stdout:
                    self.log("✅ 已是最新版")
                    messagebox.showinfo("更新", "已是最新版本")
                else:
                    self.log("✅ 更新完成，请重启 Bioclaw")
                    messagebox.showinfo("更新", "更新完成，请重新启动 Bioclaw")
            except Exception as e:
                self.log(f"❌ 检查更新失败: {str(e)}")
        
        threading.Thread(target=update, daemon=True).start()
    
    def log(self, message):
        """添加日志"""
        self.log_text.config(state='normal')
        self.log_text.insert('end', message + '\n')
        self.log_text.see('end')
        self.log_text.config(state='disabled')

def main():
    root = tk.Tk()
    app = BioclawGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
