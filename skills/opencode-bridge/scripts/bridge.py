#!/usr/bin/env python3
"""
OpenClaw → Opencode Bridge

将 OpenClaw 消息转发给 Opencode 执行
"""

import os
import sys
import json
import time
import requests
from typing import Optional, Dict, Any

# 配置
OPENCODE_URL = os.environ.get("OPENCODE_URL", "http://localhost:4096")
DEFAULT_TIMEOUT = int(os.environ.get("OPENCODE_TIMEOUT", "300"))  # 5分钟


class OpencodeBridge:
    """Opencode 桥接客户端"""
    
    def __init__(self, base_url: str = OPENCODE_URL):
        self.base_url = base_url.rstrip("/")
        self.session_id: Optional[str] = None
        
    def health_check(self) -> tuple[bool, str]:
        """检查 Opencode 服务是否可用
        
        Returns:
            tuple: (是否可用, 错误信息)
        """
        try:
            resp = requests.get(f"{self.base_url}/status", timeout=5)
            resp.raise_for_status()
            return True, ""
        except requests.exceptions.ConnectionError:
            return False, "无法连接到 Opencode。请运行: bioclaw start"
        except requests.exceptions.Timeout:
            return False, "连接 Opencode 超时。服务可能正在启动中，请稍后再试"
        except Exception as e:
            return False, f"Opencode 服务异常: {e}"
    
    def list_sessions(self) -> list:
        """获取所有 session 列表"""
        try:
            resp = requests.get(f"{self.base_url}/session", timeout=10)
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.ConnectionError:
            print("❌ 无法连接到 Opencode 服务", file=sys.stderr)
            return []
        except requests.exceptions.Timeout:
            print("❌ 获取 session 列表超时", file=sys.stderr)
            return []
        except Exception as e:
            print(f"⚠️  获取 session 列表失败: {e}", file=sys.stderr)
            return []
    
    def create_session(self, title: str = "OpenClaw Bridge Session") -> Optional[str]:
        """创建新 session"""
        try:
            resp = requests.post(
                f"{self.base_url}/session",
                json={"title": title},
                timeout=10
            )
            resp.raise_for_status()
            data = resp.json()
            self.session_id = data.get("id")
            return self.session_id
        except requests.exceptions.ConnectionError:
            print("❌ 创建 session 失败: 无法连接到 Opencode", file=sys.stderr)
            return None
        except requests.exceptions.Timeout:
            print("❌ 创建 session 超时", file=sys.stderr)
            return None
        except Exception as e:
            print(f"⚠️  创建 session 失败: {e}", file=sys.stderr)
            return None
    
    def get_or_create_session(self) -> Optional[str]:
        """获取现有 session 或创建新的"""
        # 检查是否有活跃 session
        sessions = self.list_sessions()
        
        # 按更新时间排序，找最近活跃的
        if sessions:
            sessions.sort(key=lambda x: x.get("time", {}).get("updated", 0), reverse=True)
            self.session_id = sessions[0]["id"]
            print(f"✅ 复用现有 session: {self.session_id}")
            return self.session_id
        
        # 没有则创建新的
        return self.create_session()
    
    def send_message(self, message: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """发送消息到 Opencode session"""
        sid = session_id or self.session_id
        if not sid:
            sid = self.get_or_create_session()

        if not sid:
            return {"error": "无法获取或创建 session"}

        try:
            resp = requests.post(
                f"{self.base_url}/session/{sid}/message",
                json={
                    "parts": [{"type": "text", "text": message}]
                },
                timeout=DEFAULT_TIMEOUT
            )
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.ConnectionError:
            return {"error": "发送消息失败: 无法连接到 Opencode 服务"}
        except requests.exceptions.Timeout:
            return {"error": f"请求超时 (>{DEFAULT_TIMEOUT}s)。任务可能仍在后台运行，可以通过 session 查看结果"}
        except Exception as e:
            return {"error": f"发送消息失败: {e}"}
    
    def get_messages(self, session_id: Optional[str] = None, limit: int = 10) -> list:
        """获取 session 消息历史"""
        sid = session_id or self.session_id
        if not sid:
            return []

        try:
            resp = requests.get(
                f"{self.base_url}/session/{sid}/message?limit={limit}",
                timeout=10
            )
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.ConnectionError:
            print("❌ 获取消息失败: 无法连接到 Opencode", file=sys.stderr)
            return []
        except requests.exceptions.Timeout:
            print("❌ 获取消息超时", file=sys.stderr)
            return []
        except Exception as e:
            print(f"⚠️  获取消息失败: {e}", file=sys.stderr)
            return []
    
    def extract_last_response(self, messages: list) -> str:
        """从消息列表中提取最后的助手回复"""
        # 倒序查找助手消息
        for msg in reversed(messages):
            if msg.get("role") == "assistant":
                content = msg.get("content", [])
                texts = []
                for part in content:
                    if part.get("type") == "text":
                        texts.append(part.get("text", ""))
                return "\n".join(texts)
        return ""


def main():
    """主函数 - 从命令行接收消息"""
    if len(sys.argv) < 2:
        print("用法: python bridge.py '你的消息'")
        sys.exit(1)

    message = sys.argv[1]

    # 创建桥接客户端
    bridge = OpencodeBridge()

    # 第 1 步：健康检查
    is_healthy, error_msg = bridge.health_check()
    if not is_healthy:
        print(f"❌ {error_msg}")
        print("💡 提示: 运行 'bioclaw start' 启动所有服务")
        sys.exit(1)

    # 第 2 步：确保有 session
    if not bridge.get_or_create_session():
        print("❌ 无法创建或获取 session")
        sys.exit(1)

    print(f"📤 发送消息: {message[:50]}...")

    # 第 3 步：发送消息
    result = bridge.send_message(message)

    if "error" in result:
        print(f"❌ {result['error']}")
        sys.exit(1)

    # 第 4 步：等待并获取回复
    time.sleep(2)
    messages = bridge.get_messages(limit=5)
    response = bridge.extract_last_response(messages)

    if response:
        print(response)
    else:
        print("⏳ 处理中，请稍后通过 session 查看结果")
        print(f"Session ID: {bridge.session_id}")


if __name__ == "__main__":
    main()
