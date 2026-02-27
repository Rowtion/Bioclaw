# Bioclaw 图形界面

Bioclaw 提供了图形界面（GUI），适合不熟悉命令行的用户。

## 启动 GUI

### macOS 用户（推荐）

双击 `Bioclaw.app` 即可启动图形界面。

### 命令行启动

```bash
# 如果已安装 Bioclaw
bioclaw-gui

# 或首次运行（自动安装+启动 GUI）
python3 bioclaw-gui.py
```

## GUI 功能

图形界面提供以下功能：

- 🚀 **一键启动** - 点击按钮启动 Bioclaw
- 🛑 **一键停止** - 点击按钮停止 Bioclaw
- 📊 **状态显示** - 实时显示运行状态
- 🔗 **快速访问** - 点击按钮打开 RStudio/JupyterLab
- 📝 **日志查看** - 实时显示操作日志
- 🔄 **检查更新** - 一键更新到最新版

## 界面说明

```
┌─────────────────────────────────────┐
│  🧬 Bioclaw                         │
│  生物科研环境管理器                  │
├─────────────────────────────────────┤
│  运行状态: ✅ 运行中                 │
├─────────────────────────────────────┤
│  [🚀 启动]                          │
│  [🛑 停止]                          │
│  [📊 打开 RStudio]                  │
│  [📝 打开 JupyterLab]               │
├─────────────────────────────────────┤
│  日志:                              │
│  Bioclaw 正在运行                   │
│  RStudio: http://localhost:8787    │
├─────────────────────────────────────┤
│  [❓ 帮助] [🔄 检查更新] [❌ 退出]   │
└─────────────────────────────────────┘
```

## 密码

默认密码：**bioclaw**

在 RStudio 和 JupyterLab 登录时使用。

## 首次使用

1. 启动 Bioclaw GUI
2. 点击"启动"按钮
3. 等待状态显示"运行中"
4. 点击"打开 RStudio"或"打开 JupyterLab"
5. 输入密码 `bioclaw` 登录

## 故障排除

### GUI 无法启动

确保已安装 Python 3：
```bash
python3 --version
```

如果没有安装：
- macOS: `brew install python3`
- Ubuntu: `sudo apt-get install python3 python3-tk`

### 提示"未找到 Bioclaw"

请先运行安装脚本：
```bash
bash install.sh
```

## 与命令行对比

| 功能 | GUI | 命令行 |
|------|-----|--------|
| 启动 | 点击按钮 | `bioclaw start` |
| 停止 | 点击按钮 | `bioclaw stop` |
| 查看状态 | 实时显示 | `bioclaw status` |
| 打开浏览器 | 点击按钮 | 手动输入 URL |
| 适合人群 | 小白用户 | 高级用户 |

**推荐：** 日常使用 GUI，高级操作使用命令行。
