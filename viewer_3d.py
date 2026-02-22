"""
3D Viewer Module — Wbudowany viewer PyQt6 + pyqtgraph OpenGL
===============================================================
Zastępuje trimesh.show() (wymaga pyglet) wbudowanym viewerem.

Features:
  - Obracanie (LPM), Zoom (scroll), Przesuwanie (PPM)
  - Reset View
  - Decimation for large meshes (smooth interaction)

Author: Claude (Anthropic)
Date: 2026-02-20
"""

import numpy as np
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QComboBox
)
from PyQt6.QtCore import Qt
import pyqtgraph.opengl as gl
import trimesh
import logging

logger = logging.getLogger(__name__)


class MeshViewer3D(QWidget):
    """Built-in 3D mesh viewer using pyqtgraph OpenGL."""

    MAX_FACES_DISPLAY = 200_000  # Decimate if more for smooth interaction

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mesh = None
        self._mesh_item = None
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Controls bar
        controls_widget = QWidget()
        controls_widget.setMaximumHeight(50)
        controls = QHBoxLayout(controls_widget)

        self.btn_reset = QPushButton("Reset View")
        self.btn_reset.clicked.connect(self._reset_view)
        controls.addWidget(self.btn_reset)

        self.btn_top = QPushButton("Top")
        self.btn_top.clicked.connect(lambda: self._set_view('top'))
        controls.addWidget(self.btn_top)

        self.btn_front = QPushButton("Front")
        self.btn_front.clicked.connect(lambda: self._set_view('front'))
        controls.addWidget(self.btn_front)

        self.btn_side = QPushButton("Side")
        self.btn_side.clicked.connect(lambda: self._set_view('side'))
        controls.addWidget(self.btn_side)

        controls.addStretch()

        self.info_label = QLabel("No mesh loaded")
        controls.addWidget(self.info_label)

        layout.addWidget(controls_widget)

        # 3D View
        self.gl_widget = gl.GLViewWidget()
        self.gl_widget.setBackgroundColor('w')  # White background
        layout.addWidget(self.gl_widget, stretch=1)

        # Add grid
        grid = gl.GLGridItem()
        grid.setSize(10, 10, 1)
        grid.setSpacing(1, 1, 1)
        grid.setColor((200, 200, 200, 100))
        self.gl_widget.addItem(grid)

    def set_mesh(self, mesh: trimesh.Trimesh):
        """Load a mesh into the viewer."""
        self._mesh = mesh

        # Remove old mesh
        if self._mesh_item is not None:
            self.gl_widget.removeItem(self._mesh_item)
            self._mesh_item = None

        if mesh is None:
            self.info_label.setText("No mesh loaded")
            return

        # Decimate if too large
        display_mesh = mesh
        decimated = False
        if len(mesh.faces) > self.MAX_FACES_DISPLAY:
            try:
                display_mesh = mesh.simplify_quadric_decimation(
                    face_count=self.MAX_FACES_DISPLAY
                )
                decimated = True
            except Exception:
                # Fallback: subsample faces
                indices = np.random.choice(
                    len(mesh.faces), self.MAX_FACES_DISPLAY, replace=False
                )
                display_mesh = mesh.submesh([indices], append=True)
                decimated = True

        # Get mesh data
        vertices = np.array(display_mesh.vertices, dtype=np.float32)
        faces = np.array(display_mesh.faces, dtype=np.uint32)

        # Center the mesh
        center = vertices.mean(axis=0)
        vertices -= center

        # Compute face normals for lighting
        if display_mesh.face_normals is not None:
            normals = np.array(display_mesh.face_normals, dtype=np.float32)
        else:
            normals = None

        # Color: light gray with slight shading
        colors = np.ones((len(faces), 4), dtype=np.float32)
        colors[:, 0:3] = 0.75  # Light gray

        # Simple directional shading
        if normals is not None:
            light_dir = np.array([0.3, 0.3, 1.0])
            light_dir /= np.linalg.norm(light_dir)
            shade = np.abs(normals @ light_dir)
            shade = 0.4 + 0.6 * shade  # Range 0.4-1.0
            colors[:, 0] = shade
            colors[:, 1] = shade
            colors[:, 2] = shade

        # Create mesh item
        mesh_item = gl.GLMeshItem(
            vertexes=vertices,
            faces=faces,
            faceColors=colors,
            smooth=True,
            drawEdges=False,
        )
        self.gl_widget.addItem(mesh_item)
        self._mesh_item = mesh_item

        # Set camera
        bbox = vertices.max(axis=0) - vertices.min(axis=0)
        max_dim = max(bbox)
        self.gl_widget.setCameraPosition(distance=max_dim * 2.5)

        # Update info
        info = f"Faces: {len(mesh.faces):,}"
        if decimated:
            info += f" (display: {len(display_mesh.faces):,})"
        info += f" | Vertices: {len(mesh.vertices):,}"
        info += f" | Watertight: {'Yes' if mesh.is_watertight else 'No'}"
        self.info_label.setText(info)

    def _reset_view(self):
        """Reset camera to default position."""
        if self._mesh is not None:
            bbox = self._mesh.bounds
            size = bbox[1] - bbox[0]
            max_dim = max(size)
            self.gl_widget.setCameraPosition(
                distance=max_dim * 2.5,
                elevation=30,
                azimuth=45
            )

    def _set_view(self, view: str):
        """Set predefined camera angle."""
        if self._mesh is None:
            return
        bbox = self._mesh.bounds
        size = bbox[1] - bbox[0]
        dist = max(size) * 2.5

        if view == 'top':
            self.gl_widget.setCameraPosition(distance=dist, elevation=90, azimuth=0)
        elif view == 'front':
            self.gl_widget.setCameraPosition(distance=dist, elevation=0, azimuth=0)
        elif view == 'side':
            self.gl_widget.setCameraPosition(distance=dist, elevation=0, azimuth=90)
