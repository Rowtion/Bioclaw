#!/bin/bash
# Bioclaw 卸载脚本

echo "⚠️  这将卸载 Bioclaw 并删除所有数据！"
echo ""
read -p "确定要继续吗？输入 'yes' 确认: " confirm

if [ "$confirm" != "yes" ]; then
    echo "已取消卸载"
    exit 0
fi

echo ""
echo "🗑️  正在卸载 Bioclaw..."

# 停止并删除容器
echo "  停止 Docker 容器..."
cd "$HOME/.bioclaw" 2>/dev/null && $DOCKER_COMPOSE down 2>/dev/null || true

# 删除目录
echo "  删除项目文件..."
rm -rf "$HOME/.bioclaw"

# 删除命令
echo "  删除快捷命令..."
rm -f "$HOME/.local/bin/bioclaw"

echo ""
echo "✅ Bioclaw 已卸载"
echo ""
echo "注意：Docker 镜像和卷未被删除。如需清理，请运行:"
echo "  docker system prune -a"
