---
name: opencode-bridge
description: Bridge skill that forwards messages to local Opencode server (port 4096) for execution using AIPOCH Medical Research Skills. Use when user wants to run medical research analysis, bioinformatics workflows, or any task requiring the specialized medical research skills. Triggers words: "用opencode", "opencode分析", "医学分析", "生物信息学", "run analysis", "analyze data", "perform research".
triggers:
  - "用opencode"
  - "opencode分析"
  - "医学分析"
  - "生物信息学"
  - "clinical analysis"
  - "research analysis"
  - "run analysis"
  - "analyze data"
  - "perform research"
  - "使用opencode"
---

# Opencode Bridge Skill

将 OpenClaw 消息转发到本地 Opencode 服务执行，使用 AIPOCH Medical Research Skills (200+ 医学研究技能)。

## When to Use This Skill

当用户需要进行以下操作时触发：

- **医学数据分析**: "帮我分析这个临床数据"
- **生物信息学工作流**: "用opencode做单细胞分析"
- **文献综述**: "搜索并总结最近的CRISPR论文"
- **研究设计**: "帮我设计一个随机对照试验"
- **统计计算**: "计算样本量","生存分析"
- **药物发现**: "虚拟筛选","ADMET预测"

## How to Use

### 1. 确保 Opencode 服务运行

```bash
# 检查服务状态
curl http://localhost:4096/status

# 如果未运行，启动它
~/.opencode/bin/opencode serve --port 4096
```

### 2. 使用技能

用户发送消息后，技能会自动转发给 Opencode。

### 3. 查看结果

- 短时间任务：直接在对话中返回
- 长时间任务：返回 session ID，用户可在 http://localhost:4096 查看进度

## Available Skills (AIPOCH Medical Research)

安装后，Opencode 将拥有 200+ 医学研究技能，包括：

| 类别 | 示例技能 |
|------|---------|
| 🏥 Clinical Research | 临床试验设计、患者数据分析、诊断支持 |
| 🔬 Experimental Research | 实验方案设计、protocol优化 |
| 📊 Medical Data Analysis | 统计分析、可视化、机器学习 |
| 🧬 Bioinformatics | 基因组学、转录组学、蛋白质组学 |
| 📚 Literature & Evidence | 文献搜索、系统综述、证据合成 |
| 💊 Pharmaceutical | 药物发现、ADMET预测、靶点分析 |

## Environment Variables

```bash
OPENCODE_URL=http://localhost:4096  # Opencode 服务地址
OPENCODE_TIMEOUT=300                 # 请求超时时间(秒)
```

## Scripts

- `scripts/bridge.py` - 核心桥接脚本

## Architecture

```
Feishu/OpenClaw
    ↓
opencode-bridge skill (this)
    ↓ HTTP API
Opencode Server (:4096)
    ↓
AIPOCH Medical Research Skills (200+)
    ↓
Docker Analysis Environment
```
