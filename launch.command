#!/bin/bash
# Bioclaw 启动器 - 双击运行

# 打开终端并运行
cd "$HOME/.bioclaw"

# 检查是否在运行
if docker-compose ps | grep -q "Up"; then
    echo "✅ Bioclaw 已经在运行！"
    echo ""
    echo "访问地址:"
    echo "   RStudio:    http://localhost:8787"
    echo "   JupyterLab: http://localhost:8888"
    echo ""
    echo "按回车键打开浏览器..."
    read
    open http://localhost:8787
else
    echo "🚀 正在启动 Bioclaw..."
    docker-compose up -d
    echo ""
    echo "✅ 启动完成！"
    echo ""
    echo "访问地址:"
    echo "   RStudio:    http://localhost:8787"
    echo "   JupyterLab: http://localhost:8888"
    echo "   密码: bioclaw"
    echo ""
    echo "按回车键打开浏览器..."
    read
    open http://localhost:8787
fi
