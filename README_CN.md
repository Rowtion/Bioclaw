<p align="center">
  <img src="assets/logo.svg" alt="Bioclaw Logo" width="200">
</p>

<h1 align="center">Bioclaw</h1>

<p align="center">
  <strong>开源科学 AI 研究助手</strong><br>
  OpenClaw × Opencode × 精选生物科研技能
</p>

<p align="center">
  <a href="https://github.com/Rowtion/Bioclaw/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT">
  </a>
  <a href="https://github.com/Rowtion/Bioclaw">
    <img src="https://img.shields.io/badge/Skills-150+-green.svg" alt="Skills: 150+">
  </a>
  <a href="https://opencode.dev">
    <img src="https://img.shields.io/badge/Powered%20by-Opencode-orange.svg" alt="Powered by Opencode">
  </a>
</p>

<p align="center">
  <a href="README_EN.md">English</a> | <strong>中文</strong>
</p>

---

## 目录

- [概述](#概述)
- [特性](#特性)
- [架构](#架构)
- [快速开始](#快速开始)
- [配置](#配置)
- [可用技能](#可用技能)
- [故障排除](#故障排除)
- [贡献](#贡献)
- [许可证](#许可证)

## 概述

**Bioclaw** 是一个开源的科学 AI 研究助手，让研究人员能够通过自然语言对话执行复杂的数据分析。它集成了三个强大的组件：

- **OpenClaw**: 对话式 AI 网关，连接到你常用的消息平台（飞书、WhatsApp、Slack、Discord）
- **Opencode**: 轻量级、自托管的代码执行环境
- **精选生物科研技能库**: 经人工筛选的150+高质量生物科研技能，涵盖生物信息学、数据分析、文献检索等

### 为什么选择 Bioclaw？

| 特性 | 优势 |
|---------|---------|
| 🔬 **研究就绪** | 150+ 人工筛选的高质量生物科研技能 |
| 💬 **对话式交互** | 通过你已在使用的消息应用进行交互 |
| 🔒 **自托管** | 你的数据保留在自己的机器上 |
| 🚀 **可扩展** | 轻松添加自定义技能 |
| 🐳 **容器化** | 通过 Docker 实现可复现的分析环境 |

## 特性

### 核心能力

- **临床研究**: 试验设计、患者数据分析、诊断支持
- **生物信息学**: 基因组学、转录组学、蛋白质组学工作流
- **数据分析**: 统计分析、可视化、机器学习
- **文献综合**: 文献检索、系统综述、证据分级
- **药物研究**: 药物发现、ADMET 预测、靶点分析

### 支持的平台

- 飞书
- WhatsApp
- Slack
- Discord
- 以及更多通过 OpenClaw 插件支持的平台

## 架构

```
┌─────────────────────────────────────────────────────────────────┐
│                         用户层                                   │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────────────┐    │
│  │  飞书   │  │ WhatsApp│  │  Slack  │  │     Discord     │    │
│  └────┬────┘  └────┬────┘  └────┬────┘  └────────┬────────┘    │
└───────┼────────────┼────────────┼────────────────┼─────────────┘
        │            │            │                │
        └────────────┴────────────┴────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                      OpenClaw 网关                               │
│                      (消息路由层)                                │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                   opencode-bridge 技能                          │
│              (OpenClaw → Opencode 集成)                         │
└──────────────────────────┬──────────────────────────────────────┘
                           │ HTTP API
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Opencode 服务器                              │
│              (localhost:4096, 代码执行)                         │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              精选生物科研技能库                                  │
│                      (150+ 技能)                                │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                 Docker 分析环境                                  │
│                                                                  │
│  ┌─────────────────────┐    ┌─────────────────────┐             │
│  │   RStudio 服务器    │    │     JupyterLab      │             │
│  │      (:8787)        │    │      (:8888)        │             │
│  │                     │    │                     │             │
│  │  • R 4.3.3          │    │  • Python 3         │             │
│  │  • Bioconductor     │    │  • scanpy           │             │
│  │  • DESeq2, Seurat   │    │  • biopython        │             │
│  │  • tidyverse        │    │  • scikit-learn     │             │
│  └─────────────────────┘    └─────────────────────┘             │
└─────────────────────────────────────────────────────────────────┘
```

## 快速开始

### 前提条件

- [Docker](https://docs.docker.com/get-docker/) 和 docker-compose
- [Git](https://git-scm.com/)

### 一行命令安装

```bash
# 克隆并安装所有内容
git clone --recurse-submodules https://github.com/Rowtion/Bioclaw.git
cd Bioclaw
bash setup.sh
```

就这些！`setup.sh` 脚本会自动：
- ✅ 安装 OpenClaw 和 Opencode（如果尚未安装）
- ✅ 集成 150+ 经人工筛选的生物科研技能
- ✅ 构建并启动 Docker 服务（RStudio + JupyterLab）
- ✅ 配置 OpenClaw 使用 Bioclaw 身份
- ✅ 在 4096 端口启动 Opencode 服务器

### 配置 API 访问

安装完成后，设置你的 AI 模型提供商：

```bash
# 配置 API 密钥
~/.opencode/bin/opencode auth login

# 验证安装
openclaw gateway restart
```

**注意：** `auth login` 命令将引导你选择和配置模型提供商（Anthropic Claude、OpenAI、OpenRouter 等）。

### 开始使用

通过你连接的消息平台发送消息：

**生物信息学：**
```
使用 scanpy 分析我的单细胞 RNA-seq 数据
```

**数据分析：**
```
对我的计数矩阵进行差异表达分析
```

**文献研究：**
```
搜索 PubMed 上最近的 CRISPR 论文并创建摘要
```

**编程：**
```
帮我用 matplotlib 可视化这些数据
```

结果将保存到 `./outputs/` 目录，可以在以下位置查看：
- **RStudio 服务器**: http://localhost:8787（密码：`bioclaw`）
- **JupyterLab**: http://localhost:8888（token：`bioclaw`）

## 配置

### 1. 配置 Opencode API 访问

Opencode 通过自身的认证系统管理 API 密钥：

```bash
# 登录到你的模型提供商
~/.opencode/bin/opencode auth login

# 或者如果 opencode 在你的 PATH 中
opencode auth login

# 列出已配置的提供商
opencode auth list
```

支持的提供商：Anthropic（Claude）、OpenAI，以及通过 OpenRouter 的其他提供商。

### 2. 环境文件（可选）

`.env` 文件在安装过程中自动创建。仅在需要自定义时编辑：

```bash
nano .env
```

常见自定义：
```bash
# 更改默认密码（生产环境推荐）
RSTUDIO_PASSWORD=your_secure_password

# 为长时间分析调整 Opencode 超时
OPENCODE_TIMEOUT=600
```

## 可用技能

Bioclaw 包含 **150+ 经人工筛选的高质量生物科研技能** ：

| 类别 | 技能数 | 示例 |
|----------|--------|----------|
| 🔬 **生物信息学** | 25+ | 基因组学、转录组学、蛋白质组学、代谢组学 |
| 📊 **数据分析** | 30+ | 统计分析、可视化、机器学习 |
| 📚 **文献与检索** | 20+ | 文献检索、论文分析、引文管理 |
| 🧪 **实验室工具** | 15+ | 实验方案设计、试剂计算、实验跟踪 |
| 💻 **编程** | 20+ | Python、R、数据处理、自动化 |
| 🗄️ **数据库** | 15+ | PubMed、ChEMBL、ClinicalTrials.gov 集成 |
| 📝 **出版** | 10+ | 图表创建、手稿格式化、同行评审 |
| 🔧 **工具** | 10+ | 文件转换、数据清理、API 工具 |


## 故障排除

### Opencode 连接问题

```bash
# 检查 Opencode 是否运行
curl http://localhost:4096/status

# 查看 Opencode 日志
~/.opencode/bin/opencode logs

# 重启 Opencode
~/.opencode/bin/opencode serve --port 4096
```

### Docker 服务无法启动

```bash
# 查看容器日志
docker-compose logs

# 重新构建镜像（无缓存）
docker-compose down
docker-compose build --no-cache
docker-compose up -d

# 检查容器状态
docker-compose ps
```

### 技能未触发

1. 验证 bridge 技能是否已安装：
   ```bash
   ls ~/.openclaw/workspace/skills/opencode-bridge/
   ```

2. 检查 OpenClaw 技能识别：
   ```bash
   openclaw skills check
   ```

3. 确保消息中包含触发词：
   - `用opencode`、`opencode分析`、`医学分析`、`生物信息学`
   - `clinical analysis`、`research analysis`、`run analysis`

## 项目结构

```
Bioclaw/
├── docker/
│   ├── Dockerfile              # R + Python 分析环境
│   └── entrypoint.sh           # 服务启动脚本
│
├── skills/
│   └── opencode-bridge/        # OpenClaw → Opencode 桥接
│       ├── SKILL.md
│       └── scripts/
│           └── bridge.py       # 桥接实现
│
├── docs/
│   ├── API.md                  # API 文档
│   ├── ARCHITECTURE.md         # 系统架构详情
│   └── QUICKSTART_ZH.md        # 中文快速开始
│
├── data/                       # 数据目录（挂载到容器）
├── outputs/                    # 分析输出（挂载到容器）
│
├── docker-compose.yml          # Docker 服务配置
├── setup.sh                    # 安装脚本
├── .env.template               # 环境变量模板
├── .gitignore
├── CLAUDE.md                   # AI 代理上下文
├── CONTRIBUTING.md             # 贡献指南
├── LICENSE                     # MIT 许可证
├── README.md                   # 英文文档
└── README_CN.md                # 本文档（中文）
```

## 贡献

我们欢迎贡献！请参阅 [CONTRIBUTING.md](CONTRIBUTING.md) 了解指南。

贡献方式：
- 🐛 报告 bug
- 💡 建议新功能
- 📝 改进文档
- 🔧 提交 bug 修复
- ✨ 添加新的科学计算技能

## 许可证

本项目采用 [MIT 许可证](LICENSE)。

## 致谢

- [OpenClaw](https://github.com/openclaw/openclaw) - 对话式 AI 网关
- [Opencode](https://opencode.dev) - 代码执行环境
- [精选技能 科学技能库](https://github.com/精选技能-AI/claude-scientific-skills) - 技能库

## 链接

- 🏠 **仓库**: https://github.com/Rowtion/Bioclaw
- 📖 **文档**: 参见 `docs/` 目录
- 🐛 **问题**: https://github.com/Rowtion/Bioclaw/issues
- 💬 **讨论**: https://github.com/Rowtion/Bioclaw/discussions

---

<p align="center">
  Made with ❤️ for scientific researchers worldwide
</p>
