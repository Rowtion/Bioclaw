<p align="center">
  <img src="assets/logo.svg" alt="Bioclaw Logo" width="200">
</p>

<h1 align="center">Bioclaw</h1>

<p align="center">
  <strong>开源生物科研集成安装包</strong><br>
  一键搭建完整的生物科研 AI 环境
</p>

---

## 两步安装

### 第 1 步：安装 Docker

**macOS:**
1. 访问 https://docs.docker.com/desktop/install/mac-install/
2. 下载并安装 Docker Desktop
3. 打开 Docker Desktop 应用，等待显示 "Docker Desktop is running"

**Ubuntu:**
```bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo systemctl start docker
sudo usermod -aG docker $USER
# 注销并重新登录
```

### 第 2 步：安装 Bioclaw

```bash
curl -fsSL https://raw.githubusercontent.com/Rowtion/Bioclaw/main/install.sh | bash
```

**安装过程约 5-10 分钟**，会自动：
- 下载 Bioclaw 代码
- 构建 Docker 镜像
- 启动所有服务
- 创建 `bioclaw` 命令

---

## 使用

### 启动
```bash
bioclaw start
```

### 停止
```bash
bioclaw stop
```

### 查看状态
```bash
bioclaw status
```

### 访问科研环境

- **RStudio**（数据分析）：http://localhost:8787
- **JupyterLab**（Python 编程）：http://localhost:8888

**默认密码：** `bioclaw`

---

## 这是什么？

Bioclaw 是一个**集成安装包**，帮你快速搭建：

- **OpenClaw** - AI 对话网关（连接飞书/Slack/WhatsApp）
- **Opencode** - 代码执行环境
- **Docker 工具包** - 预装 RStudio + JupyterLab + 生物科研工具

**适合人群：**
- 想快速搭建生物科研环境的用户
- 不想折腾配置的小白用户
- 需要标准化分析环境的团队

---

## 项目结构

```
~/.bioclaw/
├── docker-compose.yml     # Docker 配置
├── setup.sh               # 安装脚本
├── data/                  # 数据目录
├── outputs/               # 分析结果
└── scientific-skills/     # 生物科研技能库
```

---

## 需要帮助？

- GitHub: https://github.com/Rowtion/Bioclaw
- Issues: https://github.com/Rowtion/Bioclaw/issues
