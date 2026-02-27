# 快速开始指南

AIpoch-claw 中文文档

## 系统要求

- Docker + docker-compose
- Git
- Node.js 22+ (用于 Opencode)
- OpenClaw (已安装)

## 安装步骤

### 1. 克隆仓库

```bash
git clone https://github.com/Rowtion/AIpoch-claw.git
cd AIpoch-claw
```

### 2. 配置环境

```bash
# 复制环境模板
cp .env.template .env

# 编辑 .env 文件，添加你的 API 密钥
nano .env
```

### 3. 运行安装脚本

```bash
bash setup.sh
```

### 4. 启动 Opencode

```bash
# 在另一个终端中运行
~/.opencode/bin/opencode serve --port 4096
```

### 5. 开始使用

在 Feishu 或其他 OpenClaw 连接的平台中发送消息：

- `用opencode分析这个临床数据`
- `opencode帮我做单细胞RNA测序分析`
- `医学分析：计算样本量和功效`

## 访问分析环境

启动后可以通过浏览器访问：

| 服务 | 地址 | 默认密码 |
|------|------|---------|
| RStudio Server | http://localhost:8787 | aipoch |
| JupyterLab | http://localhost:8888 | aipoch |

## 技能分类

AIPOCH Medical Research Skills 包含以下类别：

- 🏥 临床研究 (Clinical Research)
- 🔬 实验研究 (Experimental Research)  
- 📊 医学数据分析 (Medical Data Analysis)
- 🧬 生物信息学 (Bioinformatics)
- 📚 文献与证据综合 (Literature & Evidence)
- 💊 药物研究 (Pharmaceutical)
- 🎓 教育培训 (Education)
- 🧾 基金与策略 (Grant & Strategy)

## 故障排除

### Opencode 无法连接

```bash
# 检查服务状态
curl http://localhost:4096/status

# 如果未运行，手动启动
~/.opencode/bin/opencode serve --port 4096
```

### Docker 服务启动失败

```bash
# 查看日志
docker-compose logs

# 重建并重启
docker-compose down
docker-compose build --no-cache
docker-compose up -d
```

### 技能未触发

检查 OpenClaw 是否能识别技能：

```bash
openclaw skills check
```

确保 `opencode-bridge` 技能已正确安装到 `~/.openclaw/workspace/skills/` 目录。

## 自定义技能

你可以添加自己的技能到 `scientific-skills/` 目录：

1. 创建新目录：`scientific-skills/my-skill/`
2. 编写 `SKILL.md` 文件
3. 重启 OpenClaw 或刷新技能缓存

## 更多信息

- 项目主页：https://github.com/Rowtion/AIpoch-claw
- OpenClaw：https://github.com/openclaw/openclaw
- Opencode：https://opencode.dev
- AIPOCH：https://AIPOCH.com
