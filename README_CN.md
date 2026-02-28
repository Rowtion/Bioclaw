<p align="center">
  <img src="assets/logo.svg" alt="Bioclaw Logo" width="180">
</p>

<h1 align="center">Bioclaw</h1>

<p align="center">
  <strong>开源生物科研环境一键安装包</strong><br>
  <em>10分钟搭建完整生物科研环境</em>
</p>

<p align="center">
  <a href="https://github.com/Rowtion/Bioclaw/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT">
  </a>
  <a href="#安装">
    <img src="https://img.shields.io/badge/安装-一键-green.svg" alt="Install: One-Click">
  </a>
  <a href="https://opencode.dev">
    <img src="https://img.shields.io/badge/Powered%20by-Opencode-orange.svg" alt="Powered by Opencode">
  </a>
  <a href="README.md">
    <img src="https://img.shields.io/badge/Language-中文%20%7C%20English-red.svg" alt="Language">
  </a>
</p>

<p align="center">
  <a href="README.md">English</a> • 
  <strong>简体中文</strong> • 
  <a href="#快速开始">快速开始</a> • 
  <a href="#特性">特性</a>
</p>

---

## 🎯 Bioclaw 是什么？

**Bioclaw** 是一个开源集成安装包，一条命令即可搭建完整的生物科研环境。

**10 分钟内，你将获得：**
- 🧬 **RStudio Server** - 专业统计分析环境
- 🐍 **JupyterLab** - 交互式 Python 数据科学笔记本
- 🤖 **AI 助手** - 通过聊天软件进行自然语言生物分析
- 📦 **预装工具包** - 150+ 经人工筛选的生物科研技能

**适合人群：**
- 🔬 需要标准化分析环境的研究人员
- 🎓 学习生物信息学的学生
- 👥 需要共享可复现工作流的团队
- 🚀 讨厌配置地狱的任何人

---

## ✨ 特性

<table>
<tr>
<td width="50%">

### 🚀 一行命令安装
```bash
curl -fsSL https://.../install.sh | bash
```
无需 Docker 知识，无需复杂配置，开箱即用。

</td>
<td width="50%">

### 🐳 完全容器化
一切都在 Docker 容器中运行：
- 与系统隔离
- 跨机器可复现
- 易于更新/回滚

</td>
</tr>
<tr>
<td width="50%">

### 💬 AI 驱动分析
与你的数据对话：
- "分析这个基因表达数据"
- "画一个火山图"
- "搜索 CRISPR 相关文献"

</td>
<td width="50%">

### 📊 专业工具
预配置环境：
- R 4.3 + Bioconductor (DESeq2, Seurat)
- Python 3 + scanpy, biopython
- 150+ 生物科研技能

</td>
</tr>
</table>

---

## 🚀 快速开始

### 前置条件
- **macOS** 10.14+ 或 **Linux** (Ubuntu 20.04+)
- **Docker Desktop** ([安装指南](https://docs.docker.com/get-docker/))

### 安装

```bash
# 下载并安装
curl -fsSL https://raw.githubusercontent.com/Rowtion/Bioclaw/main/install.sh | bash

# 启动 Bioclaw
bioclaw start
```

**完成！**

5-10 分钟后，访问你的科研环境：

| 服务 | 地址 | 密码 |
|---------|-----|----------|
| **RStudio** | http://localhost:8787 | `bioclaw` |
| **JupyterLab** | http://localhost:8888 | `bioclaw` |

---

## 📸 截图

<p align="center">
  <img src="docs/screenshots/rstudio.png" alt="RStudio" width="45%">
  &nbsp;&nbsp;
  <img src="docs/screenshots/jupyter.png" alt="JupyterLab" width="45%">
</p>

*左：RStudio Server 运行 DESeq2 分析 | 右：JupyterLab 运行 scanpy*

---

## 🎬 演示

```bash
# 示例 1：启动并分析数据
$ bioclaw start
✅ Bioclaw 已启动！
# 在浏览器中打开 http://localhost:8787

# 示例 2：使用 AI 助手（在 Slack/飞书中）
用户："用opencode分析我的单细胞数据"
AI："正在使用scanpy进行分析..."
[生成 UMAP 图，保存到 ./outputs/]

# 示例 3：完成后停止
$ bioclaw stop
✅ Bioclaw 已停止
```

---

## 🏗️ 架构

```
┌─────────────────────────────────────────────────────────────┐
│                        你的电脑                              │
│  ┌─────────────┐    ┌──────────────┐    ┌────────────────┐ │
│  │   浏览器     │    │  Slack/飞书   │    │    终端        │ │
│  └──────┬──────┘    └───────┬──────┘    └───────┬────────┘ │
│         │                   │                    │          │
│         └───────────────────┼────────────────────┘          │
│                             │                               │
│                    ┌────────┴────────┐                      │
│                    │  OpenClaw       │                      │
│                    │  (AI 网关)      │                      │
│                    └────────┬────────┘                      │
│                             │ HTTP                          │
│                    ┌────────┴────────┐                      │
│                    │  Opencode       │                      │
│                    │  (端口 4096)    │                      │
│                    └────────┬────────┘                      │
│                             │                               │
│         ┌───────────────────┴───────────────────┐          │
│         │           Docker 环境                  │          │
│         │  ┌───────────────────────────────┐    │          │
│         │  │  RStudio (:8787)              │    │          │
│         │  │  • R + Bioconductor           │    │          │
│         │  └───────────────────────────────┘    │          │
│         │  ┌───────────────────────────────┐    │          │
│         │  │  JupyterLab (:8888)           │    │          │
│         │  │  • Python + scanpy            │    │          │
│         │  └───────────────────────────────┘    │          │
│         └───────────────────────────────────────┘          │
└─────────────────────────────────────────────────────────────┘
```

---

## 🛠️ 安装方式

### 方式 1：一行命令安装（推荐）
```bash
curl -fsSL https://raw.githubusercontent.com/Rowtion/Bioclaw/main/install.sh | bash
```

### 方式 2：手动安装
```bash
# 克隆仓库
git clone https://github.com/Rowtion/Bioclaw.git ~/.bioclaw
cd ~/.bioclaw

# 运行安装
bash install.sh
```

### 方式 3：图形界面（macOS）
```bash
# 安装
bash install.sh

# 启动 GUI
open ~/.bioclaw/Bioclaw.app
```

---

## 📚 使用指南

### 基础命令

```bash
# 启动 Bioclaw
bioclaw start

# 查看状态
bioclaw status

# 查看日志
bioclaw logs

# 停止 Bioclaw
bioclaw stop

# 更新到最新版本
bioclaw update
```

### 处理数据

**上传数据：**
```bash
# 复制文件到 data 目录
cp my_data.csv ~/.bioclaw/data/

# 在 RStudio/JupyterLab 中访问
# 路径: /home/rstudio/data/
```

**保存结果：**
```bash
# 结果自动保存到
~/.bioclaw/outputs/
```

### AI 助手使用

1. 配置 Opencode：
```bash
~/.opencode/bin/opencode auth login
```

2. 在 Slack/飞书中发送消息：
```
"帮我用DESeq2做差异表达分析"
"画一个前50个基因的heatmap"
"搜索PubMed上COVID-19疫苗相关论文"
```

---

## 🔧 高级配置

### 修改默认密码

编辑 `.env` 文件：
```bash
cd ~/.bioclaw
nano .env

# 修改以下项：
RSTUDIO_PASSWORD=your_secure_password
JUPYTER_TOKEN=your_secure_token
```

重启：
```bash
docker-compose restart
```

### 启用远程访问

参见 [docs/REMOTE_ACCESS.md](docs/REMOTE_ACCESS.md)

### 添加自定义技能

将你的技能放入：
```
~/.bioclaw/scientific-skills/
```

---

## 🐛 故障排除

### 安装问题

**Q: 找不到 Docker**
```bash
# macOS: 从 https://docs.docker.com/desktop/install/mac-install/ 安装 Docker Desktop
# Ubuntu: sudo apt-get install docker.io
```

**Q: 端口已被占用**
```bash
# 检查什么在使用端口 8787
lsof -i :8787

# 结束进程或在 docker-compose.yml 中修改端口
```

**Q: 构建失败**
```bash
# 检查 Docker 守护进程是否运行
docker info

# 重新构建
cd ~/.bioclaw && docker-compose build --no-cache
```

### 运行问题

**Q: 无法访问 localhost:8787**
- 确保 Bioclaw 正在运行：`bioclaw status`
- 检查防火墙设置
- 尝试 http://127.0.0.1:8787

**Q: 忘记密码**
- 默认：`bioclaw`
- 在 `~/.bioclaw/.env` 中修改

更多问题？查看 [FAQ.md](FAQ.md)

---

## 🤝 贡献

欢迎贡献！

### 贡献方式
- 🐛 报告 bug
- 💡 建议新功能
- 📝 改进文档
- 🔧 提交 PR

### 开发设置

```bash
# Fork 并克隆
git clone https://github.com/YOUR_USERNAME/Bioclaw.git
cd Bioclaw

# 测试更改
bash install.sh
```

参见 [CONTRIBUTING.md](CONTRIBUTING.md) 了解指南。

---

## 📄 许可证

本项目采用 [MIT 许可证](LICENSE)。

---

## 🙏 致谢

Bioclaw 离不开以下项目：

- [OpenClaw](https://github.com/openclaw/openclaw) - AI 对话网关
- [Opencode](https://opencode.dev) - 代码执行环境
- [Docker](https://docker.com) - 容器平台
- [RStudio](https://rstudio.com) - 统计计算 IDE
- [Project Jupyter](https://jupyter.org) - 交互式计算

---

<p align="center">
  <strong>用 ❤️ 为研究社区打造</strong><br>
  <a href="https://github.com/Rowtion/Bioclaw">⭐ 在 GitHub 上给我们加星</a> • 
  <a href="https://github.com/Rowtion/Bioclaw/issues">🐛 报告问题</a>
</p>
