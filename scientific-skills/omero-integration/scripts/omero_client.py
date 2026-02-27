#!/usr/bin/env python3
"""
OMERO 显微镜数据管理工具
用于访问和分析显微镜图像数据
"""

from typing import Optional, Dict, List, Tuple, Any, Union
from dataclasses import dataclass
from pathlib import Path
import numpy as np

# 尝试导入 omero
try:
    from omero.gateway import BlitzGateway
    from omero.model import RectangleI, EllipseI, PolygonI, LineI, MaskI
    from omero.rtypes import rint, rlong, rstring, rdouble
    OMERO_AVAILABLE = True
except ImportError:
    OMERO_AVAILABLE = False
    print("警告: omero-py 未安装。运行: uv pip install omero-py")


@dataclass
class OMEROCredentials:
    """OMERO 凭据"""
    username: str
    password: str
    host: str
    port: int = 4064
    
    def to_dict(self) -> Dict:
        return {
            'username': self.username,
            'host': self.host,
            'port': self.port
        }


class OMEROClient:
    """OMERO 客户端"""
    
    def __init__(self, credentials: OMEROCredentials):
        """
        初始化 OMERO 客户端
        
        Args:
            credentials: OMERO 凭据
        """
        if not OMERO_AVAILABLE:
            raise ImportError("omero-py 未安装。运行: uv pip install omero-py")
        
        self.credentials = credentials
        self.conn: Optional[BlitzGateway] = None
    
    def connect(self) -> bool:
        """
        连接到 OMERO 服务器
        
        Returns:
            连接是否成功
        """
        try:
            self.conn = BlitzGateway(
                self.credentials.username,
                self.credentials.password,
                host=self.credentials.host,
                port=self.credentials.port
            )
            return self.conn.connect()
        except Exception as e:
            print(f"连接失败: {e}")
            return False
    
    def disconnect(self):
        """断开连接"""
        if self.conn:
            self.conn.close()
            self.conn = None
    
    def __enter__(self):
        """上下文管理器入口"""
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器出口"""
        self.disconnect()
    
    def list_projects(self) -> List[Dict]:
        """
        列出所有项目
        
        Returns:
            项目列表
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        projects = []
        for project in self.conn.listProjects():
            projects.append({
                'id': project.getId(),
                'name': project.getName(),
                'description': project.getDescription()
            })
        return projects
    
    def list_datasets(self, project_id: Optional[int] = None) -> List[Dict]:
        """
        列出数据集
        
        Args:
            project_id: 项目 ID (None 则列出所有)
        
        Returns:
            数据集列表
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        datasets = []
        
        if project_id:
            project = self.conn.getObject("Project", project_id)
            if project:
                for dataset in project.listChildren():
                    datasets.append(self._dataset_to_dict(dataset))
        else:
            for dataset in self.conn.listDatasets():
                datasets.append(self._dataset_to_dict(dataset))
        
        return datasets
    
    def _dataset_to_dict(self, dataset) -> Dict:
        """将 Dataset 对象转换为字典"""
        return {
            'id': dataset.getId(),
            'name': dataset.getName(),
            'description': dataset.getDescription(),
            'image_count': dataset.countChildren()
        }
    
    def list_images(self, dataset_id: int) -> List[Dict]:
        """
        列出数据集中的图像
        
        Args:
            dataset_id: 数据集 ID
        
        Returns:
            图像列表
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        dataset = self.conn.getObject("Dataset", dataset_id)
        if not dataset:
            return []
        
        images = []
        for image in dataset.listChildren():
            images.append(self._image_to_dict(image))
        
        return images
    
    def _image_to_dict(self, image) -> Dict:
        """将 Image 对象转换为字典"""
        return {
            'id': image.getId(),
            'name': image.getName(),
            'size_x': image.getSizeX(),
            'size_y': image.getSizeY(),
            'size_z': image.getSizeZ(),
            'size_c': image.getSizeC(),
            'size_t': image.getSizeT(),
            'pixel_size_x': image.getPixelSizeX(),
            'pixel_size_y': image.getPixelSizeY()
        }
    
    def get_image(self, image_id: int):
        """
        获取图像对象
        
        Args:
            image_id: 图像 ID
        
        Returns:
            图像对象
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        return self.conn.getObject("Image", image_id)
    
    def get_pixel_data(self, image_id: int, z: int = 0, t: int = 0, 
                       channel: Optional[int] = None) -> np.ndarray:
        """
        获取图像像素数据
        
        Args:
            image_id: 图像 ID
            z: Z 平面索引
            t: 时间帧索引
            channel: 通道索引 (None 则获取所有通道)
        
        Returns:
            像素数据数组
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        image = self.conn.getObject("Image", image_id)
        if not image:
            raise ValueError(f"图像 {image_id} 不存在")
        
        # 获取像素数据
        pixels = image.getPrimaryPixels()
        
        if channel is not None:
            plane = pixels.getPlane(z, channel, t)
        else:
            # 获取所有通道
            planes = []
            for c in range(image.getSizeC()):
                planes.append(pixels.getPlane(z, c, t))
            plane = np.stack(planes, axis=-1)
        
        return np.array(plane)
    
    def get_thumbnail(self, image_id: int, size: Tuple[int, int] = (256, 256)) -> np.ndarray:
        """
        获取图像缩略图
        
        Args:
            image_id: 图像 ID
            size: 缩略图大小
        
        Returns:
            缩略图数组
        """
        if not self.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        image = self.conn.getObject("Image", image_id)
        if not image:
            raise ValueError(f"图像 {image_id} 不存在")
        
        thumbnail = image.getThumbnail(size)
        return np.array(thumbnail)


class OMEROImageAnalyzer:
    """OMERO 图像分析器"""
    
    def __init__(self, client: OMEROClient):
        """
        初始化图像分析器
        
        Args:
            client: OMERO 客户端
        """
        self.client = client
    
    def analyze_intensity(self, image_id: int, z: int = 0, t: int = 0,
                         channel: int = 0) -> Dict:
        """
        分析图像强度
        
        Args:
            image_id: 图像 ID
            z: Z 平面
            t: 时间帧
            channel: 通道
        
        Returns:
            强度统计
        """
        data = self.client.get_pixel_data(image_id, z, t, channel)
        
        return {
            'mean': float(np.mean(data)),
            'std': float(np.std(data)),
            'min': float(np.min(data)),
            'max': float(np.max(data)),
            'median': float(np.median(data))
        }
    
    def create_mip(self, image_id: int, channel: int = 0, t: int = 0) -> np.ndarray:
        """
        创建最大强度投影 (MIP)
        
        Args:
            image_id: 图像 ID
            channel: 通道
            t: 时间帧
        
        Returns:
            MIP 图像
        """
        image = self.client.get_image(image_id)
        size_z = image.getSizeZ()
        
        # 获取所有 Z 平面
        planes = []
        for z in range(size_z):
            plane = self.client.get_pixel_data(image_id, z, t, channel)
            planes.append(plane)
        
        # 最大强度投影
        stack = np.stack(planes, axis=0)
        mip = np.max(stack, axis=0)
        
        return mip


class OMEROMetadataManager:
    """OMERO 元数据管理器"""
    
    def __init__(self, client: OMEROClient):
        """
        初始化元数据管理器
        
        Args:
            client: OMERO 客户端
        """
        self.client = client
    
    def add_key_value(self, object_type: str, object_id: int, 
                     key_values: Dict[str, str]):
        """
        添加键值对注释
        
        Args:
            object_type: 对象类型 ('Image', 'Dataset', 'Project')
            object_id: 对象 ID
            key_values: 键值对字典
        """
        if not self.client.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        obj = self.client.conn.getObject(object_type, object_id)
        if not obj:
            raise ValueError(f"对象 {object_type}:{object_id} 不存在")
        
        # 创建 MapAnnotation
        from omero.model import MapAnnotationI
        from omero.rtypes import rstring
        
        map_ann = MapAnnotationI()
        map_ann.setNs(rstring('openmicroscopy.org/omero/client/mapAnnotation'))
        
        # 添加键值对
        from omero.model import NamedValue
        nv_list = []
        for key, value in key_values.items():
            nv = NamedValue(key, value)
            nv_list.append(nv)
        
        map_ann.setMapValue(nv_list)
        
        # 保存
        update_service = self.client.conn.getUpdateService()
        map_ann = update_service.saveAndReturnObject(map_ann)
        
        # 链接到对象
        link_class = getattr(omero.model, f'{object_type}AnnotationLinkI')
        link = link_class()
        link.setParent(obj._obj)
        link.setChild(map_ann)
        update_service.saveObject(link)
    
    def add_tag(self, object_type: str, object_id: int, tag_name: str):
        """
        添加标签
        
        Args:
            object_type: 对象类型
            object_id: 对象 ID
            tag_name: 标签名称
        """
        if not self.client.conn:
            raise RuntimeError("未连接到 OMERO 服务器")
        
        from omero.model import TagAnnotationI
        
        # 创建标签
        tag = TagAnnotationI()
        tag.setTextValue(rstring(tag_name))
        
        update_service = self.client.conn.getUpdateService()
        tag = update_service.saveAndReturnObject(tag)
        
        # 链接到对象
        obj = self.client.conn.getObject(object_type, object_id)
        link_class = getattr(omero.model, f'{object_type}AnnotationLinkI')
        link = link_class()
        link.setParent(obj._obj)
        link.setChild(tag)
        update_service.saveObject(link)


def demo():
    """演示功能"""
    print("=" * 60)
    print("OMERO 显微镜数据管理工具演示")
    print("=" * 60)
    
    if not OMERO_AVAILABLE:
        print("\n错误: omero-py 未安装")
        print("请运行: uv pip install omero-py")
        return
    
    print("\nOMERO 工具可用功能:")
    print()
    print("1. OMEROClient - 连接和数据访问:")
    print("   • connect() - 连接到 OMERO 服务器")
    print("   • list_projects() - 列出项目")
    print("   • list_datasets() - 列出数据集")
    print("   • list_images() - 列出图像")
    print("   • get_pixel_data() - 获取像素数据")
    print("   • get_thumbnail() - 获取缩略图")
    print()
    print("2. OMEROImageAnalyzer - 图像分析:")
    print("   • analyze_intensity() - 强度分析")
    print("   • create_mip() - 最大强度投影")
    print()
    print("3. OMEROMetadataManager - 元数据管理:")
    print("   • add_key_value() - 添加键值对")
    print("   • add_tag() - 添加标签")
    print()
    print("=" * 60)


if __name__ == "__main__":
    demo()
