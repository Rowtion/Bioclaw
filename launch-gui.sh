#!/bin/bash
# Bioclaw GUI 启动器
# 如果没有安装，先安装；如果已安装，启动 GUI

BIoclaw_DIR="$HOME/.bioclaw"

# 检查是否已安装
if [ ! -d "$BIoclaw_DIR" ]; then
    echo "🚀 Bioclaw 未安装，开始安装..."
    echo ""
    
    # 运行安装脚本
    if [ -f "install.sh" ]; then
        bash install.sh
    else
        echo "❌ 未找到 install.sh，请先下载 Bioclaw"
        exit 1
    fi
    
    echo ""
    echo "✅ 安装完成！启动图形界面..."
    sleep 2
fi

# 启动 GUI
cd "$BIoclaw_DIR"
python3 bioclaw-gui.py
