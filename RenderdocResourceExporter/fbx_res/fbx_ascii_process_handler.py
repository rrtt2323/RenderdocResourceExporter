'''
FBX ASCII 格式处理器
'''

import inspect


class FbxAsciiProcessHandler:

    def __init__(self, config):
        self.__dict__.update(config) # 更新到参数字典里

    # fbx ASCII 格式的模板
    FBX_ASCII_TEMPLETE = """
        ; FBX 7.3.0 project file
        ; ----------------------------------------------------

        ; Object definitions
        ;------------------------------------------------------------------

        Definitions:  {

            ObjectType: "Geometry" {
                Count: 1
                PropertyTemplate: "FbxMesh" {
                    Properties70:  {
                        P: "Primary Visibility", "bool", "", "",1
                    }
                }
            }

            ObjectType: "Model" {
                Count: 1
                PropertyTemplate: "FbxNode" {
                    Properties70:  {
                        P: "Visibility", "Visibility", "", "A",1
                    }
                }
            }
        }

        ; Object properties
        ;------------------------------------------------------------------

        Objects:  {
            Geometry: 2035541511296, "Geometry::", "Mesh" {
                Vertices: *%(vertices_num)s {
                    a: %(vertices)s
                } 
                PolygonVertexIndex: *%(polygons_num)s {
                    a: %(polygons)s
                } 
                GeometryVersion: 124
                %(LayerElementNormal)s
                %(LayerElementBiNormal)s
                %(LayerElementTangent)s
                %(LayerElementColor)s
                %(LayerElementUV)s
                %(LayerElementUV2)s
                Layer: 0 {
                    Version: 100
                    %(LayerElementNormalInsert)s
                    %(LayerElementBiNormalInsert)s
                    %(LayerElementTangentInsert)s
                    %(LayerElementColorInsert)s
                    %(LayerElementUVInsert)s
                    
                }
                Layer: 1 {
                    Version: 100
                    %(LayerElementUV2Insert)s
                }
            }
            Model: 2035615390896, "Model::%(model_name)s", "Mesh" {
                Properties70:  {
                    P: "DefaultAttributeIndex", "int", "Integer", "",0
                }
            }
        }

        ; Object connections
        ;------------------------------------------------------------------

        Connections:  {
            
            ;Model::pCube1, Model::RootNode
            C: "OO",2035615390896,0
            
            ;Geometry::, Model::pCube1
            C: "OO",2035541511296,2035615390896

        }

        """

    def run(self):
        # 把所有 run_ 前缀的方法都执行一遍
        for name, func in inspect.getmembers(self, inspect.isroutine):
            if name.startswith("run_"):
                func()
        
        fbx = self.FBX_ASCII_TEMPLETE % self.ARGS
        return fbx


    # run_ 方法 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def run_vertices(self):
        vertices = [
            str(v) for values in self.vertex_data[self.POSITION].values() for v in values[:3]
        ]
        self.ARGS["vertices"] = ",".join(vertices)
        self.ARGS["vertices_num"] = len(vertices)


    def run_polygons(self):
        if self.ENGINE == "unreal":
            polygons = []
            indices = []
            for i, idx in enumerate(self.idx_dict):
                if i % 3 == 0:
                    indices.append(idx - self.min_poly)
                elif i % 3 == 1:
                    indices.append(idx - self.min_poly + 1)
                elif i % 3 == 2:
                    indices.append(idx - self.min_poly)
                    polygons.append(str(indices[0]))
                    polygons.append(str(indices[2]))
                    polygons.append(str(-indices[1]))
                    indices = []
        else:
            polygons = [
                str(idx - self.min_poly)
                if i % 3
                else str(-(idx - self.min_poly + 1))
                for i, idx in enumerate(self.idx_dict, 1)
            ]

        self.ARGS["polygons"] = ",".join(polygons)
        self.ARGS["polygons_num"] = len(polygons)


    def run_normals(self):
        if not self.vertex_data.get(self.NORMAL):
            return

        if self.ENGINE == "unreal":
            normals = []
            indices = []
            for i, idx in enumerate(self.idx_dict):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    _normals = self.vertex_data[self.NORMAL]
                    normals.extend([str(v) for v in _normals[indices[0]][:3]])
                    normals.extend([str(v) for v in _normals[indices[2]][:3]])
                    normals.extend([str(v) for v in _normals[indices[1]][:3]])
                    indices = []
        else:
            normals = [
                str(v) for values in self.value_dict[self.NORMAL] for v in values[:3]
            ]

        self.ARGS["LayerElementNormal"] = """
            LayerElementNormal: 0 {
                Version: 101
                Name: ""
                MappingInformationType: "ByPolygonVertex"
                ReferenceInformationType: "Direct"
                Normals: *%(normals_num)s {
                    a: %(normals)s
                } 
            }
        """ % {
            "normals": ",".join(normals),
            "normals_num": len(normals),
        }

        self.ARGS["LayerElementNormalInsert"] = """
            LayerElement:  {
                    Type: "LayerElementNormal"
                TypedIndex: 0
            }
        """


    def run_binormals(self):
        if not self.vertex_data.get(self.BINORMAL):
            return

        if self.ENGINE == "unreal":
            binormals = []
            indices = []
            for i, idx in enumerate(self.idx_dict):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    _binormals = self.vertex_data[self.BINORMAL]
                    binormals.extend([str(v) for v in _binormals[indices[0]][:3]])
                    binormals.extend([str(v) for v in _binormals[indices[2]][:3]])
                    binormals.extend([str(v) for v in _binormals[indices[1]][:3]])
                    indices = []
        else:
            binormals = [
                str(-v) for values in self.value_dict[self.BINORMAL] for v in values[:3]
            ]

        self.ARGS[
            "LayerElementBiNormal"
        ] = """
            LayerElementBinormal: 0 {
                Version: 101
                Name: "map1"
                MappingInformationType: "ByVertice"
                ReferenceInformationType: "Direct"
                Binormals: *%(binormals_num)s {
                    a: %(binormals)s
                } 
                BinormalsW: *%(binormalsW_num)s {
                    a: %(binormalsW)s
                } 
            }
        """ % {
            "binormals": ",".join(binormals),
            "binormals_num": len(binormals),
            "binormalsW": ",".join(["1" for i in range(self.idx_len)]),
            "binormalsW_num": self.idx_len,
        }
        self.ARGS[
            "LayerElementBiNormalInsert"
        ] = """
            LayerElement:  {
                    Type: "LayerElementBinormal"
                TypedIndex: 0
            }
        """


    def run_tangents(self):
        if not self.vertex_data.get(self.TANGENT):
            return

        if self.ENGINE == "unreal":
            tangents = []
            indices = []
            for i, idx in enumerate(self.idx_dict):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    _tangents = self.vertex_data[self.TANGENT]
                    tangents.extend([str(v) for v in _tangents[indices[0]][:3]])
                    tangents.extend([str(v) for v in _tangents[indices[2]][:3]])
                    tangents.extend([str(v) for v in _tangents[indices[1]][:3]])
                    indices = []
        else:
            tangents = [
                str(v) for values in self.value_dict[self.TANGENT] for v in values[:3]
            ]

        self.ARGS[
            "LayerElementTangent"
        ] = """
            LayerElementTangent: 0 {
                Version: 101
                Name: "map1"
                MappingInformationType: "ByPolygonVertex"
                ReferenceInformationType: "Direct"
                Tangents: *%(tangents_num)s {
                    a: %(tangents)s
                } 
            }
        """ % {
            "tangents": ",".join(tangents),
            "tangents_num": len(tangents),
        }

        self.ARGS[
            "LayerElementTangentInsert"
        ] = """
                LayerElement:  {
                    Type: "LayerElementTangent"
                    TypedIndex: 0
                }
        """


    def run_color(self):
        if not self.vertex_data.get(self.COLOR):
            return

        if self.ENGINE == "unreal":
            colors = []
            indices = []
            for i, idx in enumerate(self.idx_dict):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    _colors = self.vertex_data[self.COLOR]
                    colors.extend([str(v) for v in _colors[indices[0]]])
                    colors.extend([str(v) for v in _colors[indices[2]]])
                    colors.extend([str(v) for v in _colors[indices[1]]])
                    indices = []
        else:
            colors = [
                # str(v) if i % 4 else "1"
                str(v)
                for values in self.value_dict[self.COLOR]
                for i, v in enumerate(values, 1)
            ]

        self.ARGS[
            "LayerElementColor"
        ] = """
            LayerElementColor: 0 {
                Version: 101
                Name: "colorSet1"
                MappingInformationType: "ByPolygonVertex"
                ReferenceInformationType: "IndexToDirect"
                Colors: *%(colors_num)s {
                    a: %(colors)s
                } 
                ColorIndex: *%(colors_indices_num)s {
                    a: %(colors_indices)s
                } 
            }
        """ % {
            "colors": ",".join(colors),
            "colors_num": len(colors),
            "colors_indices": ",".join([str(i) for i in range(self.idx_len)]),
            "colors_indices_num": self.idx_len,
        }
        self.ARGS[
            "LayerElementColorInsert"
        ] = """
            LayerElement:  {
                Type: "LayerElementColor"
                TypedIndex: 0
            }
        """


    def run_uv(self):
        if not self.vertex_data.get(self.UV):
            return

        if self.ENGINE == "unreal":
            uvs_indices = []
            indices = []
            for i, idx in enumerate(self.idx_list):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    uvs_indices.append(str(indices[0]))
                    uvs_indices.append(str(indices[2]))
                    uvs_indices.append(str(indices[1]))
                    indices = []
            uvs_indices = ",".join(uvs_indices)
        else:
            uvs_indices = self.idx_data

        uvs = [
            # NOTE flip y axis
            str(1 - v if i else v)
            for values in self.vertex_data[self.UV].values()
            for i, v in enumerate(values)
        ]

        self.ARGS[
            "LayerElementUV"
        ] = """
            LayerElementUV: 0 {
                Version: 101
                Name: "map1"
                MappingInformationType: "ByPolygonVertex"
                ReferenceInformationType: "IndexToDirect"
                UV: *%(uvs_num)s {
                    a: %(uvs)s
                } 
                UVIndex: *%(uvs_indices_num)s {
                    a: %(uvs_indices)s
                } 
            }
        """ % {
            "uvs": ",".join(uvs),
            "uvs_num": len(uvs),
            "uvs_indices": uvs_indices,
            "uvs_indices_num": self.idx_len,
        }

        self.ARGS[
            "LayerElementUVInsert"
        ] = """
            LayerElement:  {
                Type: "LayerElementUV"
                TypedIndex: 0
            }
        """


    def run_uv2(self):
        if not self.vertex_data.get(self.UV2):
            return

        uvs_indices = ""
        if self.ENGINE == "unreal":
            uvs_indices = []
            indices = []
            for i, idx in enumerate(self.idx_list):
                if i % 3 == 0:
                    indices.append(idx)
                elif i % 3 == 1:
                    indices.append(idx)
                elif i % 3 == 2:
                    indices.append(idx)
                    uvs_indices.append(str(indices[0]))
                    uvs_indices.append(str(indices[2]))
                    uvs_indices.append(str(indices[1]))
                    indices = []
            uvs_indices = ",".join(uvs_indices)
        else:
            uvs_indices = self.idx_data

        uvs = [
            # NOTE flip y axis
            str(1 - v if i else v)
            for values in self.vertex_data[self.UV2].values()
            for i, v in enumerate(values)
        ]

        self.ARGS[
            "LayerElementUV2"
        ] = """
            LayerElementUV: 1 {
                Version: 101
                Name: "map2"
                MappingInformationType: "ByPolygonVertex"
                ReferenceInformationType: "IndexToDirect"
                UV: *%(uvs_num)s {
                    a: %(uvs)s
                } 
                UVIndex: *%(uvs_indices_num)s {
                    a: %(uvs_indices)s
                } 
            }
        """ % {
            "uvs": ",".join(uvs),
            "uvs_num": len(uvs),
            "uvs_indices": uvs_indices,
            "uvs_indices_num": self.idx_len,
        }

        self.ARGS[
            "LayerElementUV2Insert"
        ] = """
            LayerElement:  {
                Type: "LayerElementUV"
                TypedIndex: 1
            }
        """
